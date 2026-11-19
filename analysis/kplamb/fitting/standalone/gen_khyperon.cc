// Standalone generator for gamma p -> K+ Lambda, Lambda -> p pi-
//
// Modeled on halld_sim/src/programs/Simulation/gen_amp, but with no GlueX
// software dependence beyond ROOT and AmpTools:
//
//   1. throw the beam energy (flat, 1/E bremsstrahlung, or from a histogram)
//   2. throw -t from exp(-b_t|t|) + uFrac exp(-b_u|u|) (forward t-channel and
//      backward u-channel peaks) and build the two-body K+ Lambda final state,
//      rejecting kinematically forbidden -t
//   3. decay Lambda -> p pi- isotropically in the Lambda rest frame
//   4. accept/reject using the KHyperon amplitude evaluated with AmpTools
//      (skipped with -f to produce flat phase space for genmc/accmc)
//
// The output tree uses the flat FSRoot layout with MC-prefixed branches
// (MCEnPB, MCEnP2 = K+, MCEnP1a = p, MCEnP1b = pi-, MCEnP1 = Lambda), so the
// files can be read directly by FSRootDataReader as data, genmc or accmc.

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/Kinematics.h"

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector3.h"

#include "KHyperon.h"

using namespace std;

namespace {

// PDG masses (GeV)
const double kMassProton = 0.938272;
const double kMassKPlus  = 0.493677;
const double kMassLambda = 1.115683;
const double kMassPiMinus = 0.139570;

enum BeamShape { kBeamFlat, kBeamBrem, kBeamHist };

// Production model for the unpolarized cross section:
//   dsigma/dt ~ exp( tSlope * t ) + uFrac * exp( uSlope * u )
// i.e. a forward (t-channel) peak and a backward (u-channel) peak.  t and u
// are negative except very near the kinematic limits, so the slopes are
// positive numbers in GeV^-2.  uFrac is the height of the backward peak at
// u = 0 relative to the forward peak at t = 0.
struct ProductionModel
{
  double tSlope = 0.0;
  double uSlope = 0.0;
  double uFrac = 0.0;

  double weight( double t, double u ) const
  {
    return exp( tSlope * t ) + uFrac * exp( uSlope * u );
  }
};

const double kSumMass2 = kMassProton*kMassProton + kMassKPlus*kMassKPlus + kMassLambda*kMassLambda;

void
usage( const char* prog )
{
  cout << endl << " Usage: " << prog << " -o <output.root> [options]" << endl << endl
       << "   -c <file>      AmpTools config with a single reaction using the KHyperon" << endl
       << "                  amplitude (not needed with -f)" << endl
       << "   -o <file>      output ROOT file" << endl
       << "   -n <N>         number of events [10000]" << endl
       << "   -s <seed>      random seed (0 = unique seed from TRandom3) [0]" << endl
       << "   -a <E>         minimum photon energy in GeV [8.2]" << endl
       << "   -b <E>         maximum photon energy in GeV [8.8]" << endl
       << "   -beam <shape>  photon spectrum: flat or brem (1/E) [flat]" << endl
       << "   -flux <file> <hist>  sample the photon spectrum from a TH1 (restricted to -a/-b)" << endl
       << "   -tmin <val>    minimum -t in GeV^2 [0.0]" << endl
       << "   -tmax <val>    maximum -t in GeV^2 [3.0]" << endl
       << "   -t <slope>     t-channel slope b_t in exp(-b_t|t|); 0 gives flat -t [0.0]" << endl
       << "   -u <slope>     u-channel slope b_u in exp(-b_u|u|) [0.0]" << endl
       << "   -uFrac <val>   u-channel peak height relative to t-channel at t = u = 0;" << endl
       << "                  0 disables the u-channel term [0.0]" << endl
       << "                  (use e.g. -tmax 20 to cover the backward region)" << endl
       << "   -tree <name>   output tree name [ntFSGlueX_100000000_100000]" << endl
       << "   -batch <N>     events per AmpTools batch [100000]" << endl
       << "   -f             generate flat phase space (no amplitude weighting)" << endl
       << endl;
}

double
breakupMomentum( double M, double m1, double m2 )
{
  double a = ( M*M - ( m1 + m2 )*( m1 + m2 ) );
  double b = ( M*M - ( m1 - m2 )*( m1 - m2 ) );
  if( a < 0.0 || b < 0.0 ) return -1.0;
  return sqrt( a * b ) / ( 2.0 * M );
}

// Throws gamma p -> K+ Lambda, Lambda -> p pi- events.
class KLambdaGenerator
{
public:

  KLambdaGenerator( TRandom3& rng, double eLow, double eHigh, BeamShape shape,
                    TH1* fluxHist, double tMin, double tMax,
                    const ProductionModel& model ) :
    m_rng( rng ), m_eLow( eLow ), m_eHigh( eHigh ), m_shape( shape ),
    m_flux( fluxHist ), m_tMin( tMin ), m_tMax( tMax ),
    m_model( model ), m_useU( model.uFrac > 0.0 ), m_wMax( 0.0 ),
    m_target( 0, 0, 0, kMassProton )
  {
    if( m_useU ) m_wMax = findWeightMaximum();
  }

  // particle order: beam, K+, p, pi- (matches KHyperon and FSRootDataReader)
  Kinematics* generate()
  {
    const long kMaxTries = 100000000;
    for( long iTry = 0; iTry < kMaxTries; ++iTry ){

      const double eBeam = throwBeamEnergy();
      // t-channel only: throw exp(-b_t|t|) directly.  With a u-channel term,
      // throw -t flat and get the shape from accept/reject below.
      const double absT = ( m_useU ? m_rng.Uniform( m_tMin, m_tMax ) : throwAbsT() );

      TLorentzVector beam( 0, 0, eBeam, eBeam );
      TLorentzVector cm = beam + m_target;
      const double W = cm.M();
      if( W < kMassKPlus + kMassLambda ) continue;

      // CM-frame quantities for t = (beam - K+)^2
      const CMKinematics cmk( W );
      const double cosTheta = cmk.cosTheta( -absT );
      if( cosTheta < -1.0 || cosTheta > 1.0 ) continue; // kinematically forbidden -t

      if( m_useU ){
        // u = (gamma - Lambda)^2 = sum(m^2) - s - t
        const double w = m_model.weight( -absT, kSumMass2 - W*W + absT );
        if( w > m_wMax ){
          cout << "  WARNING: production weight " << w << " exceeds ceiling " << m_wMax
               << "; raising it" << endl;
          m_wMax = 1.1 * w;
        }
        if( w < m_rng.Uniform() * m_wMax ) continue;
      }

      const double pK = cmk.pK;
      const double eK = cmk.eK;

      TVector3 kMom;
      kMom.SetMagThetaPhi( pK, acos( cosTheta ), m_rng.Uniform( -M_PI, M_PI ) );
      TLorentzVector kaon( kMom, eK );
      TLorentzVector lambda( -kMom, sqrt( pK*pK + kMassLambda*kMassLambda ) );

      // Lambda -> p pi- isotropic in the Lambda rest frame
      const double q = breakupMomentum( kMassLambda, kMassProton, kMassPiMinus );
      TVector3 pMom;
      pMom.SetMagThetaPhi( q, acos( m_rng.Uniform( -1.0, 1.0 ) ), m_rng.Uniform( -M_PI, M_PI ) );
      TLorentzVector proton( pMom, sqrt( q*q + kMassProton*kMassProton ) );
      TLorentzVector pion( -pMom, sqrt( q*q + kMassPiMinus*kMassPiMinus ) );
      proton.Boost( lambda.BoostVector() );
      pion.Boost( lambda.BoostVector() );

      // CM -> lab
      const TVector3 cmBoost = cm.BoostVector();
      kaon.Boost( cmBoost );
      proton.Boost( cmBoost );
      pion.Boost( cmBoost );

      vector< TLorentzVector > parts;
      parts.push_back( beam );
      parts.push_back( kaon );
      parts.push_back( proton );
      parts.push_back( pion );
      return new Kinematics( parts, 1.0 );
    }

    cerr << "ERROR: could not generate a kinematically allowed event after "
         << kMaxTries << " tries; check -a/-b and -tmin/-tmax" << endl;
    exit( 1 );
    return NULL;
  }

private:

  // two-body gamma p -> K+ Lambda kinematics in the CM frame
  struct CMKinematics
  {
    double pGam, pK, eK;

    explicit CMKinematics( double W ) :
      pGam( breakupMomentum( W, 0.0, kMassProton ) ),
      pK( breakupMomentum( W, kMassKPlus, kMassLambda ) ),
      eK( sqrt( pK*pK + kMassKPlus*kMassKPlus ) )
    {}

    // t = mK^2 - 2 (Egam EK - pGam pK cos)
    double t( double cosTheta ) const
    { return kMassKPlus*kMassKPlus - 2.0*( pGam*eK - pGam*pK*cosTheta ); }

    double cosTheta( double t ) const
    { return ( t - kMassKPlus*kMassKPlus + 2.0*pGam*eK ) / ( 2.0*pGam*pK ); }
  };

  // Upper bound of the production weight in the beam-energy and -t windows.
  // Each term is monotonic, so bound it by its value at the edge where it
  // peaks: the t term at the largest t, the u term at the largest u (the
  // most backward allowed point at the lowest beam energy in the scan).
  double findWeightMaximum() const
  {
    const int nE = 21;
    double tHiMax = -1.0e9, uHiMax = -1.0e9;
    for( int iE = 0; iE < nE; ++iE ){
      const double eBeam = m_eLow + ( m_eHigh - m_eLow ) * iE / ( nE - 1 );
      const double s = kMassProton*kMassProton + 2.0*kMassProton*eBeam;
      const double W = sqrt( s );
      if( W < kMassKPlus + kMassLambda ) continue;

      const CMKinematics cmk( W );
      const double tLo = max( cmk.t( -1.0 ), -m_tMax );
      const double tHi = min( cmk.t( 1.0 ), -m_tMin );
      if( tLo > tHi ) continue;

      tHiMax = max( tHiMax, tHi );
      uHiMax = max( uHiMax, kSumMass2 - s - tLo );
    }

    if( tHiMax < -1.0e8 ){
      cerr << "ERROR: no kinematically allowed region for -a/-b and -tmin/-tmax" << endl;
      exit( 1 );
    }
    return 1.05 * ( exp( m_model.tSlope * tHiMax ) +
                    m_model.uFrac * exp( m_model.uSlope * uHiMax ) );
  }

  double throwBeamEnergy()
  {
    if( m_shape == kBeamHist ){
      double e;
      do { e = m_flux->GetRandom( &m_rng ); } while( e < m_eLow || e > m_eHigh );
      return e;
    }
    if( m_shape == kBeamBrem ){
      // dN/dE ~ 1/E  ->  E = Elow * (Ehigh/Elow)^u
      return m_eLow * pow( m_eHigh / m_eLow, m_rng.Uniform() );
    }
    return m_rng.Uniform( m_eLow, m_eHigh );
  }

  double throwAbsT()
  {
    const double slope = m_model.tSlope;
    if( slope <= 0.0 ) return m_rng.Uniform( m_tMin, m_tMax );
    // inverse CDF of exp(-b t) truncated to [tMin, tMax]
    const double a = exp( -slope * m_tMin );
    const double b = exp( -slope * m_tMax );
    return -log( a - m_rng.Uniform() * ( a - b ) ) / slope;
  }

  TRandom3& m_rng;
  double m_eLow, m_eHigh;
  BeamShape m_shape;
  TH1* m_flux;
  double m_tMin, m_tMax;
  ProductionModel m_model;
  bool m_useU;
  double m_wMax;
  TLorentzVector m_target;
};

// Flat FSRoot-style output with MC-prefixed branches plus diagnostics.
class FSRootMCWriter
{
public:

  FSRootMCWriter( const string& fileName, const string& treeName ) :
    m_file( TFile::Open( fileName.c_str(), "RECREATE" ) ), m_count( 0 )
  {
    if( !m_file || m_file->IsZombie() ){
      cerr << "ERROR: could not create " << fileName << endl;
      exit( 1 );
    }
    m_file->cd();
    m_tree = new TTree( treeName.c_str(), treeName.c_str() );

    // order in the arrays: B, 2 (K+), 1a (p), 1b (pi-), 1 (Lambda)
    const char* labels[kN] = { "B", "2", "1a", "1b", "1" };
    for( int i = 0; i < kN; ++i ){
      const string l = labels[i];
      m_tree->Branch( ( "MCEnP" + l ).c_str(), &m_en[i], ( "MCEnP" + l + "/D" ).c_str() );
      m_tree->Branch( ( "MCPxP" + l ).c_str(), &m_px[i], ( "MCPxP" + l + "/D" ).c_str() );
      m_tree->Branch( ( "MCPyP" + l ).c_str(), &m_py[i], ( "MCPyP" + l + "/D" ).c_str() );
      m_tree->Branch( ( "MCPzP" + l ).c_str(), &m_pz[i], ( "MCPzP" + l + "/D" ).c_str() );
    }

    m_hE      = new TH1D( "gen_EBeam", "Beam energy;E_{#gamma} (GeV);Events", 200, 0.0, 12.0 );
    // -t and -u cover the full angular range at GlueX energies, 0.01 GeV^2 bins
    m_hT      = new TH1D( "gen_t", "-t;-t (GeV^{2});Events", 1600, 0.0, 16.0 );
    m_hU      = new TH1D( "gen_u", "-u;-u (GeV^{2});Events", 1700, -1.0, 16.0 );
    m_hCosCM  = new TH1D( "gen_cosThetaCM", "K^{+} cos#theta_{CM};cos#theta_{CM}^{K^{+}};Events",
                          200, -1.0, 1.0 );
    m_hMLam   = new TH1D( "gen_MLambda", "M(p#pi^{-});M (GeV);Events", 200, 1.10, 1.13 );
    m_hCosX   = new TH1D( "gen_cosThetaX", "cos#theta_{x};cos#theta_{x};Events", 100, -1.0, 1.0 );
    m_hCosY   = new TH1D( "gen_cosThetaY", "cos#theta_{y};cos#theta_{y};Events", 100, -1.0, 1.0 );
    m_hCosZ   = new TH1D( "gen_cosThetaZ", "cos#theta_{z};cos#theta_{z};Events", 100, -1.0, 1.0 );
    m_hPhi    = new TH1D( "gen_Phi", "#Phi;#Phi (rad);Events", 100, -M_PI, M_PI );
    m_hCosYPhi = new TH2D( "gen_cosThetaY_vs_Phi", "cos#theta_{y} vs #Phi;#Phi (rad);cos#theta_{y}",
                           50, -M_PI, M_PI, 50, -1.0, 1.0 );
  }

  ~FSRootMCWriter()
  {
    // histograms are detached (TH1::AddDirectory(kFALSE)), so write them explicitly
    m_file->cd();
    m_tree->Write();
    TH1* hists[] = { m_hE, m_hT, m_hU, m_hCosCM, m_hMLam,
                     m_hCosX, m_hCosY, m_hCosZ, m_hPhi, m_hCosYPhi };
    for( TH1* h : hists ){
      h->Write();
      delete h;
    }
    m_file->Close();
    delete m_file;
  }

  void write( const Kinematics& kin )
  {
    const TLorentzVector beam = kin.particle( 0 );
    const TLorentzVector k  = kin.particle( 1 );
    const TLorentzVector p  = kin.particle( 2 );
    const TLorentzVector pi = kin.particle( 3 );
    const TLorentzVector lam = p + pi;

    const TLorentzVector vecs[kN] = { beam, k, p, pi, lam };
    for( int i = 0; i < kN; ++i ){
      m_en[i] = vecs[i].E();
      m_px[i] = vecs[i].Px();
      m_py[i] = vecs[i].Py();
      m_pz[i] = vecs[i].Pz();
    }
    m_tree->Fill();
    ++m_count;

    // diagnostics, same frame definitions as KHyperon::calcUserVars
    TLorentzRotation lamBoost( -lam.BoostVector() );
    TLorentzVector kLam = lamBoost * k;
    TLorentzVector pLam = lamBoost * p;
    TVector3 y = ( beam.Vect().Unit().Cross( -k.Vect().Unit() ) ).Unit();
    TVector3 z = -1.0 * kLam.Vect().Unit();
    TVector3 x = y.Cross( z ).Unit();
    TVector3 eps( 1.0, 0.0, 0.0 );
    const double Phi = atan2( y.Dot( eps ), beam.Vect().Unit().Dot( eps.Cross( y ) ) );
    const double cosY = cos( pLam.Vect().Angle( y ) );

    // production variables: t = (p - Lambda)^2, u = (gamma - Lambda)^2, and
    // the K+ polar angle relative to the beam in the gamma p CM frame
    const TLorentzVector target( 0, 0, 0, kMassProton );
    TLorentzVector kCM = k;
    kCM.Boost( -( beam + target ).BoostVector() );

    m_hE->Fill( beam.E() );
    m_hT->Fill( -( target - lam ).M2() );
    m_hU->Fill( -( beam - lam ).M2() );
    m_hCosCM->Fill( kCM.Vect().Unit().Dot( beam.Vect().Unit() ) );
    m_hMLam->Fill( lam.M() );
    m_hCosX->Fill( cos( pLam.Vect().Angle( x ) ) );
    m_hCosY->Fill( cosY );
    m_hCosZ->Fill( cos( pLam.Vect().Angle( z ) ) );
    m_hPhi->Fill( Phi );
    m_hCosYPhi->Fill( Phi, cosY );
  }

  long count() const { return m_count; }

private:

  static const int kN = 5;

  TFile* m_file;
  TTree* m_tree;
  long m_count;

  double m_en[kN], m_px[kN], m_py[kN], m_pz[kN];

  TH1D *m_hE, *m_hT, *m_hU, *m_hCosCM, *m_hMLam, *m_hCosX, *m_hCosY, *m_hCosZ, *m_hPhi;
  TH2D* m_hCosYPhi;
};

} // namespace

int
main( int argc, char* argv[] )
{
  string cfgFile, outName;
  string fluxFile, fluxHistName;
  string treeName = "ntFSGlueX_100000000_100000";
  long nEvents = 10000;
  int batchSize = 100000;
  unsigned int seed = 0;
  double eLow = 8.2, eHigh = 8.8;
  double tMin = 0.0, tMax = 3.0;
  BeamShape shape = kBeamFlat;
  ProductionModel model;
  bool genFlat = false;

  for( int i = 1; i < argc; ++i ){
    const string arg = argv[i];
    const bool hasNext = ( i + 1 < argc );

    if( arg == "-h" || arg == "--help" ){ usage( argv[0] ); return 0; }
    else if( arg == "-c" && hasNext ) cfgFile = argv[++i];
    else if( arg == "-o" && hasNext ) outName = argv[++i];
    else if( arg == "-n" && hasNext ) nEvents = atol( argv[++i] );
    else if( arg == "-s" && hasNext ) seed = static_cast< unsigned int >( atol( argv[++i] ) );
    else if( arg == "-a" && hasNext ) eLow = atof( argv[++i] );
    else if( arg == "-b" && hasNext ) eHigh = atof( argv[++i] );
    else if( arg == "-tmin" && hasNext ) tMin = atof( argv[++i] );
    else if( arg == "-tmax" && hasNext ) tMax = atof( argv[++i] );
    else if( arg == "-t" && hasNext ) model.tSlope = atof( argv[++i] );
    else if( arg == "-u" && hasNext ) model.uSlope = atof( argv[++i] );
    else if( arg == "-uFrac" && hasNext ) model.uFrac = atof( argv[++i] );
    else if( arg == "-tree" && hasNext ) treeName = argv[++i];
    else if( arg == "-batch" && hasNext ) batchSize = atoi( argv[++i] );
    else if( arg == "-f" ) genFlat = true;
    else if( arg == "-beam" && hasNext ){
      const string s = argv[++i];
      if( s == "flat" ) shape = kBeamFlat;
      else if( s == "brem" ) shape = kBeamBrem;
      else { cerr << "Unknown beam shape " << s << endl; usage( argv[0] ); return 1; }
    }
    else if( arg == "-flux" && i + 2 < argc ){
      fluxFile = argv[++i];
      fluxHistName = argv[++i];
      shape = kBeamHist;
    }
    else{
      cerr << "Unrecognized or incomplete argument: " << arg << endl;
      usage( argv[0] );
      return 1;
    }
  }

  if( outName.empty() || ( cfgFile.empty() && !genFlat ) ){
    cerr << "An output file (-o) and, unless -f is given, a config file (-c) are required." << endl;
    usage( argv[0] );
    return 1;
  }
  if( eLow <= 0.0 || eHigh <= eLow || tMax <= tMin || tMin < 0.0 || nEvents <= 0 || batchSize <= 0 ){
    cerr << "Invalid beam energy, t range, event count or batch size." << endl;
    return 1;
  }
  if( model.tSlope < 0.0 || model.uSlope < 0.0 || model.uFrac < 0.0 ){
    cerr << "-t, -u and -uFrac must be non-negative." << endl;
    return 1;
  }
  if( model.uFrac > 0.0 && model.uSlope <= 0.0 ){
    cerr << "-uFrac > 0 needs a positive u-channel slope (-u <slope>)." << endl;
    return 1;
  }

  TH1::AddDirectory( kFALSE );
  TRandom3 rng( seed );
  cout << "TRandom3 seed: " << rng.GetSeed() << endl;

  unique_ptr< TH1 > fluxHist;
  if( shape == kBeamHist ){
    TFile* f = TFile::Open( fluxFile.c_str() );
    TH1* h = ( f ? dynamic_cast< TH1* >( f->Get( fluxHistName.c_str() ) ) : NULL );
    if( !h ){
      cerr << "ERROR: could not read histogram " << fluxHistName << " from " << fluxFile << endl;
      return 1;
    }
    fluxHist.reset( static_cast< TH1* >( h->Clone( "gen_flux" ) ) );
    f->Close();
    delete f;
  }

  cout << "Production model: exp(-" << model.tSlope << " |t|)";
  if( model.uFrac > 0.0 ) cout << " + " << model.uFrac << " exp(-" << model.uSlope << " |u|)";
  cout << ",  " << tMin << " < -t < " << tMax << " GeV^2" << endl;

  KLambdaGenerator generator( rng, eLow, eHigh, shape, fluxHist.get(), tMin, tMax, model );

  // --- flat phase space: no AmpTools needed ---
  if( genFlat ){
    cout << "Generating " << nEvents << " flat phase-space events into " << outName << endl;
    FSRootMCWriter writer( outName, treeName );
    for( long i = 0; i < nEvents; ++i ){
      Kinematics* kin = generator.generate();
      writer.write( *kin );
      delete kin;
    }
    cout << "Wrote " << writer.count() << " events" << endl;
    return 0;
  }

  // --- amplitude-weighted generation ---
  AmpToolsInterface::registerAmplitude( KHyperon() );

  ConfigFileParser parser( cfgFile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  AmpToolsInterface ati( cfgInfo, AmpToolsInterface::kMCGeneration );

  if( cfgInfo->reactionList().size() != 1 ){
    cerr << "ERROR: " << cfgFile << " must define exactly one reaction" << endl;
    return 1;
  }
  const string reaction = cfgInfo->reactionList()[0]->reactionName();

  cout << endl << "Generating " << nEvents << " events for reaction " << reaction
       << " into " << outName << endl;

  FSRootMCWriter writer( outName, treeName );
  double maxInten = -1.0;
  long nThrown = 0;

  while( writer.count() < nEvents ){

    ati.clearEvents();
    for( int i = 0; i < batchSize; ++i ){
      Kinematics* kin = generator.generate();
      ati.loadEvent( kin, i, batchSize );
      delete kin;
    }

    const double batchMax = ati.processEvents( reaction );
    if( maxInten < 0.0 ){
      // fix the accept/reject ceiling from the first batch (with a safety
      // margin) so all batches are sampled from the same distribution
      maxInten = 1.2 * batchMax;
      cout << "  maximum intensity for accept/reject: " << maxInten << endl;
    }
    else if( batchMax > maxInten ){
      cout << "  WARNING: batch intensity " << batchMax << " exceeds ceiling "
           << maxInten << "; raising it (consider a larger -batch)" << endl;
      maxInten = 1.2 * batchMax;
    }

    for( int i = 0; i < batchSize && writer.count() < nEvents; ++i ){
      ++nThrown;
      if( ati.intensity( i ) > rng.Uniform() * maxInten ){
        Kinematics* kin = ati.kinematics( i );
        writer.write( *kin );
        delete kin;
      }
    }

    cout << "  " << writer.count() << " / " << nEvents << " events accepted" << endl;
  }

  cout << "  accept/reject efficiency: "
       << static_cast< double >( writer.count() ) / nThrown << endl;

  return 0;
}
