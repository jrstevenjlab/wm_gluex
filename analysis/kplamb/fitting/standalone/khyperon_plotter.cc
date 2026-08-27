#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/PlotGenerator.h"

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "FSRootDataReader.h"
#include "KHyperon.h"

using namespace std;

namespace
{

struct KHyperonObservables
{
  double mass;
  double cosThetax;
  double cosThetay;
  double cosThetaz;
  double cosThetaHyp;
  double phiHyp;
  double PhiHyp;
};

enum KHyperonProjectionIndex
{
  kHyperonMass = 0,
  kCosThetax,
  kCosThetay,
  kCosThetaz,
  kCosThetaHyp,
  kphiHyp,
  kPhiHyp,
  kNumKHyperonPlots
};

struct PlotSet
{
  TH1F* mass;
  TH1F* cosThetaHyp;
  TH1F* phiHyp;
  TH1F* cosThetaX;
  TH1F* cosThetaY;
  TH1F* cosThetaZ;
  TH1F* Phi;
  double sumWeights;

  PlotSet() :
    mass( NULL ),
    cosThetaHyp( NULL ),
    phiHyp( NULL ),
    cosThetaX( NULL ),
    cosThetaY( NULL ),
    cosThetaZ( NULL ),
    Phi( NULL ),
    sumWeights( 0.0 )
  {}

  void scale( double factor )
  {
    if( factor == 1.0 ) return;

    mass->Scale( factor );
    cosThetaHyp->Scale( factor );
    phiHyp->Scale( factor );
    cosThetaX->Scale( factor );
    cosThetaY->Scale( factor );
    cosThetaZ->Scale( factor );
    Phi->Scale( factor );
  }
};

KHyperonObservables computeKHyperonObservables( Kinematics* kin );

class KHyperonPlotGenerator : public PlotGenerator
{
public:

  KHyperonPlotGenerator( const FitResults& results ) :
    PlotGenerator( results )
  {
    createHistograms();
  }

protected:

  void projectEvent( Kinematics* kin ) override
  {
    projectEvent( kin, "" );
  }

  void projectEvent( Kinematics* kin, const string& reactionName ) override
  {
    (void)reactionName;

    const KHyperonObservables obs = computeKHyperonObservables( kin );
    fillHistogram( kHyperonMass, obs.mass );
    fillHistogram( kCosThetax, obs.cosThetax );
    fillHistogram( kCosThetay, obs.cosThetay );
    fillHistogram( kCosThetaz, obs.cosThetaz );
    fillHistogram( kCosThetaHyp, obs.cosThetaHyp );
    fillHistogram( kphiHyp, obs.phiHyp );
    fillHistogram( kPhiHyp, obs.PhiHyp );
  }

private:

  void createHistograms()
  {
    bookHistogram( kHyperonMass,
                   new Histogram1D( 500, 0.5, 2.0, "MHyp", "Invariant Mass of Hyperon" ) );
    bookHistogram( kCosThetax,
                   new Histogram1D( 200, -1.0, 1.0, "cosThetax", "cos( #theta_{x} ) of Hyperon decay" ) );
    bookHistogram( kCosThetay,
                   new Histogram1D( 200, -1.0, 1.0, "cosThetay", "cos( #theta_{y} ) of Hyperon decay" ) );
    bookHistogram( kCosThetaz,
                   new Histogram1D( 200, -1.0, 1.0, "cosThetaz", "cos( #theta_{z} ) of Hyperon decay" ) );
    bookHistogram( kCosThetaHyp,
                   new Histogram1D( 200, -1.0, 1.0, "cosThetaHyp", "cos( #theta ) of Hyperon decay" ) );
    bookHistogram( kphiHyp,
                   new Histogram1D( 180, -TMath::Pi(), TMath::Pi(), "phiHyp", "#phi of Hyperon decay" ) );
    bookHistogram( kPhiHyp,
                   new Histogram1D( 180, -TMath::Pi(), TMath::Pi(), "PhiHyp", "#Phi of Hyperon decay" ) );
  }
};

string
defaultOutputFileName( const string& fitFile )
{
  const string suffix = ".fit";

  if( fitFile.size() >= suffix.size() &&
      fitFile.compare( fitFile.size() - suffix.size(), suffix.size(), suffix ) == 0 ){
    return fitFile.substr( 0, fitFile.size() - suffix.size() ) + ".plots.root";
  }

  return fitFile + ".plots.root";
}

string
defaultParameterCsvFileName( const string& fitFile )
{
  const string suffix = ".fit";

  if( fitFile.size() >= suffix.size() &&
      fitFile.compare( fitFile.size() - suffix.size(), suffix.size(), suffix ) == 0 ){
    return fitFile.substr( 0, fitFile.size() - suffix.size() ) + "_parameters.csv";
  }

  return fitFile + "_parameters.csv";
}

string
csvNumber( double value )
{
  if( !std::isfinite( value ) ){
    return "nan";
  }

  ostringstream stream;
  stream << setprecision( 15 ) << value;
  return stream.str();
}

void
writeSelectedParameterCsv( const FitResults& results, const string& csvFile )
{
  static const array< string, 5 > selectedParameters = {
    "Sigma", "Ox", "P", "T", "Oz"
  };

  const vector< string >& names = results.parNameList();
  const vector< double >& values = results.parValueList();
  const vector< vector< double > >& covariance = results.errorMatrix();

  ofstream csv( csvFile.c_str() );
  if( !csv ){
    cerr << "Could not create parameter CSV file " << csvFile << endl;
    return;
  }

  csv << "parameter,value,error" << endl;

  for( array< string, 5 >::const_iterator par = selectedParameters.begin();
       par != selectedParameters.end();
       ++par ){

    vector< string >::const_iterator found = find( names.begin(), names.end(), *par );
    if( found == names.end() ){
      cerr << "WARNING: parameter " << *par
           << " was not found in the fit results; writing nan to "
           << csvFile << endl;
      csv << *par << ",nan,nan" << endl;
      continue;
    }

    const size_t index = static_cast< size_t >( found - names.begin() );
    const double value = ( index < values.size() ?
                           values[index] :
                           numeric_limits< double >::quiet_NaN() );

    double error = numeric_limits< double >::quiet_NaN();
    if( index < covariance.size() &&
        index < covariance[index].size() &&
        covariance[index][index] >= 0.0 ){
      error = sqrt( covariance[index][index] );
    }

    csv << *par << ","
        << csvNumber( value ) << ","
        << csvNumber( error ) << endl;
  }

  cout << endl << "Wrote selected parameter CSV to " << csvFile << endl;
}

void
printParameterSummary( const FitResults& results )
{
  const vector< string >& names = results.parNameList();
  const vector< double >& values = results.parValueList();

  cout << "Parameters:" << endl;
  for( size_t i = 0; i < names.size(); ++i ){
    cout << "  " << setw( 32 ) << left << names[i]
         << " = " << setw( 14 ) << right << values[i];

    try{
      cout << " +/- " << results.parError( names[i] );
    }
    catch( ... ){
      cout << " +/- n/a";
    }

    cout << endl;
  }
}

void
labelMatrixAxes( TH2D& hist, const vector< string >& labels )
{
  for( size_t i = 0; i < labels.size(); ++i ){
    hist.GetXaxis()->SetBinLabel( static_cast< int >( i + 1 ), labels[i].c_str() );
    hist.GetYaxis()->SetBinLabel( static_cast< int >( i + 1 ), labels[i].c_str() );
  }
}

void
fillCovarianceHistogram( TH2D& hist,
                         const vector< string >& labels,
                         const vector< vector< double > >& covariance )
{
  const size_t nRows = min( labels.size(), covariance.size() );

  for( size_t i = 0; i < nRows; ++i ){
    const size_t nCols = min( labels.size(), covariance[i].size() );
    for( size_t j = 0; j < nCols; ++j ){
      hist.SetBinContent( static_cast< int >( i + 1 ),
                          static_cast< int >( j + 1 ),
                          covariance[i][j] );
    }
  }
}

void
fillCorrelationHistogram( TH2D& hist,
                          const vector< string >& labels,
                          const vector< vector< double > >& covariance )
{
  const size_t nRows = min( labels.size(), covariance.size() );

  for( size_t i = 0; i < nRows; ++i ){
    const double sigmaI = ( i < covariance.size() && i < covariance[i].size() ?
                            sqrt( fabs( covariance[i][i] ) ) : 0.0 );

    const size_t nCols = min( labels.size(), covariance[i].size() );
    for( size_t j = 0; j < nCols; ++j ){
      const double sigmaJ = ( j < covariance.size() && j < covariance[j].size() ?
                              sqrt( fabs( covariance[j][j] ) ) : 0.0 );

      double corr = 0.0;
      if( sigmaI > 0.0 && sigmaJ > 0.0 ){
        corr = covariance[i][j] / ( sigmaI * sigmaJ );
      }
      if( i == j && sigmaI > 0.0 ){
        corr = 1.0;
      }

      hist.SetBinContent( static_cast< int >( i + 1 ),
                          static_cast< int >( j + 1 ),
                          corr );
    }
  }
}

TH1F*
bookComparisonHistogram( vector< unique_ptr< TH1F > >& ownedHistograms,
                         const string& name,
                         const string& title,
                         int nBins,
                         double low,
                         double high )
{
  ownedHistograms.push_back( unique_ptr< TH1F >( new TH1F( name.c_str(), title.c_str(), nBins, low, high ) ) );
  ownedHistograms.back()->Sumw2();
  return ownedHistograms.back().get();
}

PlotSet
bookComparisonPlotSet( vector< unique_ptr< TH1F > >& ownedHistograms,
                       const string& prefix,
                       const string& category,
                       const string& categoryTitle )
{
  PlotSet plots;

  plots.mass = bookComparisonHistogram( ownedHistograms,
                                        prefix + "MHyperon" + category,
                                        "MHyperon " + categoryTitle + ";Invariant Mass of Hyperon;Events",
                                        500,
                                        0.5,
                                        2.0 );
  plots.cosThetaHyp = bookComparisonHistogram( ownedHistograms,
                                               prefix + "cosThetaHyp" + category,
                                               "cosThetaHyp " + categoryTitle + ";cos(#theta);Events",
                                               200,
                                               -1.0,
                                               1.0 );
  plots.phiHyp = bookComparisonHistogram( ownedHistograms,
                                           prefix + "phiHyp" + category,
                                           "phiHyp " + categoryTitle + ";#phi;Events",
                                           180,
                                           -TMath::Pi(),
                                           TMath::Pi() );
  plots.cosThetaX = bookComparisonHistogram( ownedHistograms,
                                             prefix + "cosThetaX" + category,
                                             "cosThetaX " + categoryTitle + ";cos(#theta_{x});Events",
                                             200,
                                             -1.0,
                                             1.0 );
  plots.cosThetaY = bookComparisonHistogram( ownedHistograms,
                                             prefix + "cosThetaY" + category,
                                             "cosThetaY " + categoryTitle + ";cos(#theta_{y});Events",
                                             200,
                                             -1.0,
                                             1.0 );
  plots.cosThetaZ = bookComparisonHistogram( ownedHistograms,
                                             prefix + "cosThetaZ" + category,
                                             "cosThetaZ " + categoryTitle + ";cos(#theta_{z});Events",
                                             200,
                                             -1.0,
                                             1.0 );
  plots.Phi = bookComparisonHistogram( ownedHistograms,
                                       prefix + "Phi" + category,
                                       "Phi " + categoryTitle + ";#Phi;Events",
                                       180,
                                       -TMath::Pi(),
                                       TMath::Pi() );

  return plots;
}

std::array< TH1*, kNumKHyperonPlots >
makeComparisonTargets( PlotSet& plots )
{
  return { plots.mass,
           plots.cosThetaX,
           plots.cosThetaY,
           plots.cosThetaZ,
           plots.cosThetaHyp,
           plots.phiHyp,
           plots.Phi };
}

std::array< TH1*, kNumKHyperonPlots >
makeSummaryTargets( TH1D& massHist,
                    TH1D& cosThetaxHist,
                    TH1D& cosThetayHist,
                    TH1D& cosThetazHist,
                    TH1D& cosThetaHypHist,
                    TH1D& phiHypHist,
                    TH1D& PhiHypHist )
{
  return { &massHist,
           &cosThetaxHist,
           &cosThetayHist,
           &cosThetazHist,
           &cosThetaHypHist,
           &phiHypHist,
           &PhiHypHist };
}

KHyperonObservables
computeKHyperonObservables( Kinematics* kin )
{
  TLorentzVector beam = kin->particle( 0 );
  TLorentzVector k = kin->particle( 1 );
  TLorentzVector y1 = kin->particle( 2 );
  TLorentzVector y2 = kin->particle( 3 );

  TLorentzVector hyperon = y1 + y2;
  TLorentzRotation hyperonBoost( -hyperon.BoostVector() );

  TLorentzVector k_hyperon = hyperonBoost * k;
  TLorentzVector y1_hyperon = hyperonBoost * y1;

  // normal to the production plane (formed by beam and kaon)
  TVector3 y = ( beam.Vect().Unit().Cross( -k.Vect().Unit() ) ).Unit();

  // choose helicity frame: z-axis opposite kaon in hyperon rest frame
  TVector3 z = -1.0 * k_hyperon.Vect().Unit();
  TVector3 x = y.Cross( z ).Unit();
  TVector3 angles( ( y1_hyperon.Vect() ).Dot( x ),
                   ( y1_hyperon.Vect() ).Dot( y ),
                   ( y1_hyperon.Vect() ).Dot( z ) );

  TVector3 eps( 1.0, 0.0, 0.0 ); // reference beam polarization vector at 0 degrees

  KHyperonObservables obs;
  obs.mass = hyperon.M();
  obs.cosThetax = cos( y1_hyperon.Vect().Angle( x ) );
  obs.cosThetay = cos( y1_hyperon.Vect().Angle( y ) );
  obs.cosThetaz = cos( y1_hyperon.Vect().Angle( z ) );
  obs.cosThetaHyp = angles.CosTheta();
  obs.phiHyp = angles.Phi();
  obs.PhiHyp = atan2( y.Dot( eps ), beam.Vect().Unit().Dot( eps.Cross( y ) ) );
  return obs;
}

void
addProjectedHistograms( PlotGenerator& plotGen,
                        const string& reactionName,
                        unsigned int type,
                        const std::array< TH1*, kNumKHyperonPlots >& targets,
                        double* sumWeights = NULL )
{
  for( size_t i = 0; i < targets.size(); ++i ){
    Histogram* proj = plotGen.projection( static_cast< unsigned int >( i ), reactionName, type );
    if( !proj || !targets[i] ) continue;

    if( i == 0 && sumWeights != NULL ){
      *sumWeights += proj->entries();
    }

    unique_ptr< TH1 > rootHist( proj->toRoot() );
    targets[i]->Add( rootHist.get() );
  }
}

void
copyPlotSetToSummaryHistograms( const PlotSet& plots,
                                TH1D& massHist,
                                TH1D& cosThetaxHist,
                                TH1D& cosThetayHist,
                                TH1D& cosThetazHist,
                                TH1D& cosThetaHypHist,
                                TH1D& phiHypHist,
                                TH1D& PhiHypHist )
{
  massHist.Add( plots.mass );
  cosThetaxHist.Add( plots.cosThetaX );
  cosThetayHist.Add( plots.cosThetaY );
  cosThetazHist.Add( plots.cosThetaZ );
  cosThetaHypHist.Add( plots.cosThetaHyp );
  phiHypHist.Add( plots.phiHyp );
  PhiHypHist.Add( plots.Phi );
}

void
scaleToMatchData( PlotSet& dataPlots, PlotSet& accPlots, PlotSet& genPlots )
{
  if( dataPlots.sumWeights <= 0.0 || accPlots.sumWeights <= 0.0 ){
    return;
  }

  const double scale = dataPlots.sumWeights / accPlots.sumWeights;
  accPlots.scale( scale );
  genPlots.scale( scale );
}

} // namespace

int
main( int argc, char* argv[] )
{
  if( argc < 2 ){
    cout << "Usage: khyperon_plotter <fit file> [output root file] [parameter csv file]" << endl;
    return 0;
  }

  const string fitFile = argv[1];
  const string outputFile = ( argc > 2 ? argv[2] : defaultOutputFileName( fitFile ) );
  const string parameterCsvFile = ( argc > 3 ? argv[3] : defaultParameterCsvFileName( fitFile ) );

  FitResults results( fitFile );

  if( !results.valid() ){
    cerr << "Could not load fit results from " << fitFile << endl;
    return 1;
  }

  cout << endl << " *** KHyperon Fit Summary *** " << endl << endl;
  cout << "Fit file: " << fitFile << endl;
  cout << "Likelihood: " << results.likelihood() << endl;
  cout << "MINUIT command status: " << results.lastMinuitCommandStatus() << endl;
  cout << "MINUIT error matrix status: " << results.eMatrixStatus() << endl;
  cout << "MINUIT best minimum: " << results.bestMinimum() << endl;
  cout << endl;

  printParameterSummary( results );
  writeSelectedParameterCsv( results, parameterCsvFile );

  cout << endl;
  cout << "Amplitudes:" << endl;
  const vector< string > amps = results.ampList();
  for( vector< string >::const_iterator amp = amps.begin(); amp != amps.end(); ++amp ){
    complex< double > prod = results.productionParameter( *amp );
    cout << "  " << *amp << "  prod = (" << prod.real() << ", " << prod.imag() << ")";
    cout << "  scale = " << results.ampScale( *amp ) << endl;
  }

  const vector< string >& parNames = results.parNameList();
  const vector< vector< double > >& covariance = results.errorMatrix();

  if( !parNames.empty() && !covariance.empty() ){
    TH2D covarianceHist( "covariance_matrix",
                         "Covariance matrix;parameter;parameter",
                         static_cast< int >( parNames.size() ),
                         0,
                         static_cast< int >( parNames.size() ),
                         static_cast< int >( parNames.size() ),
                         0,
                         static_cast< int >( parNames.size() ) );
    TH2D correlationHist( "correlation_matrix",
                          "Correlation matrix;parameter;parameter",
                          static_cast< int >( parNames.size() ),
                          0,
                          static_cast< int >( parNames.size() ),
                          static_cast< int >( parNames.size() ),
                          0,
                          static_cast< int >( parNames.size() ) );

    covarianceHist.SetStats( false );
    correlationHist.SetStats( false );
    labelMatrixAxes( covarianceHist, parNames );
    labelMatrixAxes( correlationHist, parNames );
    fillCovarianceHistogram( covarianceHist, parNames, covariance );
    fillCorrelationHistogram( correlationHist, parNames, covariance );

    AmpToolsInterface::registerAmplitude( KHyperon() );
    AmpToolsInterface::registerDataReader( FSRootDataReader() );

    TH1::AddDirectory( kFALSE );
    TH1D massHist( "MHyp", "Invariant Mass of Hyperon", 500, 0.5, 2.0 );
    TH1D cosThetaxHist( "cosThetax", "cos( #theta_{x} ) of Hyperon decay", 200, -1.0, 1.0 );
    TH1D cosThetayHist( "cosThetay", "cos( #theta_{y} ) of Hyperon decay", 200, -1.0, 1.0 );
    TH1D cosThetazHist( "cosThetaz", "cos( #theta_{z} ) of Hyperon decay", 200, -1.0, 1.0 );
    TH1D cosThetaHypHist( "cosThetaHyp", "cos( #theta ) of Hyperon decay", 200, -1.0, 1.0 );
    TH1D phiHypHist( "phiHyp", "#phi of Hyperon decay", 180, -TMath::Pi(), TMath::Pi() );
    TH1D PhiHypHist( "PhiHyp", "#Phi of Hyperon decay", 180, -TMath::Pi(), TMath::Pi() );

    massHist.Sumw2();
    cosThetaxHist.Sumw2();
    cosThetayHist.Sumw2();
    cosThetazHist.Sumw2();
    cosThetaHypHist.Sumw2();
    phiHypHist.Sumw2();
    PhiHypHist.Sumw2();

    vector< unique_ptr< TH1F > > ownedComparisonHistograms;
    PlotSet dataPlots = bookComparisonPlotSet( ownedComparisonHistograms, "", "dat", "data" );
    PlotSet bkgndPlots = bookComparisonPlotSet( ownedComparisonHistograms, "", "bkgnd", "background" );
    PlotSet accPlots = bookComparisonPlotSet( ownedComparisonHistograms, "", "acc", "accepted MC" );
    PlotSet genPlots = bookComparisonPlotSet( ownedComparisonHistograms, "", "gen", "generated MC" );

    KHyperonPlotGenerator plotGen( results );
    const vector< string > reactions = plotGen.reactions();

    const array< TH1*, kNumKHyperonPlots > dataTargets = makeComparisonTargets( dataPlots );
    const array< TH1*, kNumKHyperonPlots > bkgndTargets = makeComparisonTargets( bkgndPlots );
    const array< TH1*, kNumKHyperonPlots > accTargets = makeComparisonTargets( accPlots );
    const array< TH1*, kNumKHyperonPlots > genTargets = makeComparisonTargets( genPlots );

    cout << endl << "Projecting KHyperon histograms with AmpTools PlotGenerator:" << endl;
    for( vector< string >::const_iterator reaction = reactions.begin();
         reaction != reactions.end();
         ++reaction ){

      cout << "  " << *reaction << endl;
      addProjectedHistograms( plotGen, *reaction, PlotGenerator::kData, dataTargets, &dataPlots.sumWeights );
      addProjectedHistograms( plotGen, *reaction, PlotGenerator::kBkgnd, bkgndTargets, &bkgndPlots.sumWeights );
      addProjectedHistograms( plotGen, *reaction, PlotGenerator::kAccMC, accTargets, &accPlots.sumWeights );
      addProjectedHistograms( plotGen, *reaction, PlotGenerator::kGenMC, genTargets, &genPlots.sumWeights );
    }

    copyPlotSetToSummaryHistograms( dataPlots,
                                    massHist,
                                    cosThetaxHist,
                                    cosThetayHist,
                                    cosThetazHist,
                                    cosThetaHypHist,
                                    phiHypHist,
                                    PhiHypHist );

    scaleToMatchData( dataPlots, accPlots, genPlots );

    TFile outFile( outputFile.c_str(), "RECREATE" );
    if( outFile.IsZombie() ){
      cerr << "Could not create output ROOT file " << outputFile << endl;
      return 1;
    }

    outFile.cd();
    covarianceHist.Write();
    correlationHist.Write();
    massHist.Write();
    cosThetaxHist.Write();
    cosThetayHist.Write();
    cosThetazHist.Write();
    cosThetaHypHist.Write();
    phiHypHist.Write();
    PhiHypHist.Write();

    for( vector< unique_ptr< TH1F > >::iterator hist = ownedComparisonHistograms.begin();
         hist != ownedComparisonHistograms.end();
         ++hist ){
      (*hist)->Write();
    }

    outFile.Write();
    outFile.Close();

    cout << endl;
    cout << "Wrote summary histograms to " << outputFile << endl;
    cout << "  covariance_matrix and correlation_matrix are labeled by parNameList() order" << endl;
    cout << "  angle histograms: MHyp, cosThetax, cosThetay, cosThetaz, cosThetaHyp, phiHyp, PhiHyp" << endl;
    cout << "  comparison histograms use MHyperon{dat,bkgnd,acc,gen}, cosThetaHyp{dat,bkgnd,acc,gen}, ..." << endl;
  }
  else{
    cerr << "No floating parameters were found in the fit results; skipping ROOT output." << endl;
  }

  return 0;
}
