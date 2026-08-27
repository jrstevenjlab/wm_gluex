#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <set>
#include <vector>

#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/AmpToolsInterface.h"

#include "MinuitInterface/MinuitParameter.h"
#include "MinuitInterface/MinuitParameterManager.h"

#include "FSRootDataReader.h"
#include "KHyperon.h"

using namespace std;

namespace {

struct TrialResult
{
  double likelihood;
  int status;
  int eMatrixStatus;
  bool good;
  map< string, double > startValues;
  map< string, double > finalValues;

  TrialResult() :
    likelihood( numeric_limits< double >::infinity() ),
    status( -1 ),
    eMatrixStatus( -1 ),
    good( false ),
    startValues(),
    finalValues()
  {}
};

void
usage()
{
  cout << "Usage:\n"
       << "  fit -c <config file> [-r <trials>] [-m <max iterations>] [-s <seed file>]\n"
       << endl;
}

map< string, pair< double, double > >
parRangeMap( const ConfigurationInfo* cfgInfo )
{
  map< string, pair< double, double > > ranges;
  vector< vector< string > > args = cfgInfo->userKeywordArguments( "parRange" );

  for( vector< vector< string > >::const_iterator it = args.begin();
       it != args.end();
       ++it ){

    if( it->size() != 3 ) continue;

    const string& parName = (*it)[0];
    double low = atof( (*it)[1].c_str() );
    double high = atof( (*it)[2].c_str() );
    ranges[parName] = make_pair( low, high );
  }

  return ranges;
}

map< string, double >
snapshotParameterValues( const MinuitParameterManager& parManager )
{
  map< string, double > values;
  for( const MinuitParameter* par : parManager ){
    values[ par->name() ] = par->value();
  }

  return values;
}

set< string >
floatingParameterNames( const MinuitParameterManager& parManager )
{
  set< string > names;
  for( const MinuitParameter* par : parManager ){
    if( par->floating() ){
      names.insert( par->name() );
    }
  }

  return names;
}

void
applyParameterValues( MinuitParameterManager& parManager,
                      const map< string, double >& values )
{
  if( values.empty() ) return;

  for( MinuitParameter* par : parManager ){
    map< string, double >::const_iterator it = values.find( par->name() );
    if( it != values.end() ){
      // Bulk parameter seeding must stay silent here.  Notifying observers
      // while the fit state is being staged can force AmpTools to rebuild the
      // covariance matrix before MINUIT has synchronized its floating-parameter
      // bookkeeping, which triggers an assertion in updateParCovariance().
      *par->valuePtr() = it->second;
    }
  }
}

map< string, double >
randomizeParameters( const map< string, double >& baselineValues,
                     const map< string, pair< double, double > >& ranges,
                     const set< string >& floatingNames,
                     mt19937& rng,
                     bool randomize )
{
  map< string, double > appliedValues = baselineValues;
  if( !randomize ) return appliedValues;

  uniform_real_distribution< double > uniform01( 0.0, 1.0 );
  uniform_real_distribution< double > productionUniform( -500.0, 500.0 );
  for( map< string, pair< double, double > >::const_iterator range = ranges.begin();
       range != ranges.end();
       ++range ){

    map< string, double >::iterator value = appliedValues.find( range->first );
    if( value == appliedValues.end() ) continue;
    if( floatingNames.find( value->first ) == floatingNames.end() ) continue;

    double low = range->second.first;
    double high = range->second.second;
    value->second = ( low == high ? low : low + ( high - low ) * uniform01( rng ) );
  }

  for( map< string, double >::iterator value = appliedValues.begin();
       value != appliedValues.end();
       ++value ){

    if( floatingNames.find( value->first ) == floatingNames.end() ) continue;
    if( ranges.find( value->first ) != ranges.end() ) continue;

    if( value->first.find( "::" ) != string::npos ){
      value->second = productionUniform( rng );
    }
  }

  return appliedValues;
}

void
applyParameterValues( ConfigurationInfo* cfgInfo,
                      const map< string, double >& values )
{
  if( values.empty() ) return;

  vector< ParameterInfo* > pars = cfgInfo->parameterList();
  for( vector< ParameterInfo* >::iterator par = pars.begin();
       par != pars.end();
       ++par ){

    map< string, double >::const_iterator it = values.find( (*par)->parName() );
    if( it != values.end() ){
      (*par)->setValue( it->second );
    }
  }
}

TrialResult
runFitTrial( AmpToolsInterface& ati,
             const map< string, double >& startValues,
             int maxCalls )
{
  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  applyParameterValues( fitManager->parameterManager(), startValues );

  fitManager->migradMinimization();

  TrialResult result;
  result.likelihood = ati.likelihood();
  result.status = fitManager->status();
  result.eMatrixStatus = fitManager->eMatrixStatus();
  result.good = ( result.status == 0 || result.eMatrixStatus == 3 );
  result.startValues = startValues;
  result.finalValues = snapshotParameterValues( fitManager->parameterManager() );

  return result;
}

void
syncConfigurationInfo( AmpToolsInterface& ati )
{
  ConfigurationInfo* cfgInfo = ati.configurationInfo();
  MinuitParameterManager& parManager = ati.minuitMinimizationManager()->parameterManager();
  map< string, double > values = snapshotParameterValues( parManager );

  applyParameterValues( cfgInfo, values );

  vector< ReactionInfo* > reactions = cfgInfo->reactionList();
  for( vector< ReactionInfo* >::iterator reaction = reactions.begin();
       reaction != reactions.end();
       ++reaction ){

    const string& reactionName = (*reaction)->reactionName();
    IntensityManager* intenMan = ati.intensityManager( reactionName );
    if( !intenMan ) continue;

    vector< AmplitudeInfo* > amps = cfgInfo->amplitudeList( reactionName );
    for( vector< AmplitudeInfo* >::iterator amp = amps.begin();
         amp != amps.end();
         ++amp ){

      const string& ampName = (*amp)->fullName();
      complex< double > prodFactor = intenMan->productionFactor( ampName );
      double scale = static_cast< double >( intenMan->getScale( ampName ) );
      (*amp)->setValue( ( scale == 0.0 ? prodFactor : prodFactor / scale ) );
    }

    vector< PDFInfo* > pdfs = cfgInfo->pdfList( reactionName );
    for( vector< PDFInfo* >::iterator pdf = pdfs.begin();
         pdf != pdfs.end();
         ++pdf ){

      const string& pdfName = (*pdf)->fullName();
      complex< double > prodFactor = intenMan->productionFactor( pdfName );
      double scale = static_cast< double >( intenMan->getScale( pdfName ) );
      (*pdf)->setValue( ( scale == 0.0 ? prodFactor.real() : prodFactor.real() / scale ) );
    }
  }
}

} // namespace

int
main( int argc, char* argv[] )
{
  string cfgFile;
  string seedFile;
  int numTrials = 1;
  int maxCalls = 5000;

  for( int i = 1; i < argc; ++i ){
    string arg = argv[i];

    if( arg == "-h" || arg == "--help" ){
      usage();
      return 0;
    }
    else if( arg == "-c" && i + 1 < argc ){
      cfgFile = argv[++i];
    }
    else if( arg == "-r" && i + 1 < argc ){
      numTrials = atoi( argv[++i] );
    }
    else if( arg == "-m" && i + 1 < argc ){
      maxCalls = atoi( argv[++i] );
    }
    else if( arg == "-s" && i + 1 < argc ){
      seedFile = argv[++i];
    }
    else if( !arg.empty() && arg[0] != '-' && cfgFile.empty() ){
      cfgFile = arg;
    }
    else{
      cerr << "Unrecognized argument: " << arg << endl;
      usage();
      return 1;
    }
  }

  if( cfgFile.empty() ){
    usage();
    return 1;
  }

  cout << endl << " *** Performing the KHyperon Fit *** " << endl << endl;
  cout << "Config file: " << cfgFile << endl;
  cout << "Random trials: " << numTrials << endl;
  cout << "Max iterations: " << maxCalls << endl;
  if( !seedFile.empty() ){
    cout << "Seed file: " << seedFile << endl;
  }
  cout << endl;

  AmpToolsInterface::registerAmplitude( KHyperon() );
  AmpToolsInterface::registerDataReader( FSRootDataReader() );

  ConfigFileParser parser( cfgFile );
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  AmpToolsInterface ATI( cfgInfo );
  MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
  fitManager->setPrecision( 1E-13 );
  fitManager->setStrategy( 1 );
  if( maxCalls > 0 ){
    fitManager->setMaxIterations( maxCalls );
  }

  map< string, double > baselineValues =
    snapshotParameterValues( fitManager->parameterManager() );
  map< string, pair< double, double > > ranges = parRangeMap( cfgInfo );
  if( numTrials > 1 && ranges.empty() ){
    cout << "WARNING: no parRange keywords were found; AmpPars will stay at the baseline values, but production amplitudes will still be randomized." << endl;
  }

  mt19937 rng( random_device{}() );
  int trialCount = max( 1, numTrials );
  set< string > floatingNames = floatingParameterNames( fitManager->parameterManager() );

  auto makeStartValues = [&]() {
    return randomizeParameters( baselineValues, ranges, floatingNames, rng, true );
  };

  if( trialCount == 1 ){
    map< string, double > startValues = makeStartValues();
    TrialResult result = runFitTrial( ATI, startValues, maxCalls );
    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << result.likelihood << endl;
    cout << "MINUIT STATUS: " << result.status << "  E-MATRIX STATUS: " << result.eMatrixStatus << endl;
    syncConfigurationInfo( ATI );
    ATI.finalizeFit();
    if( !seedFile.empty() ){
      ATI.fitResults()->writeSeed( seedFile );
    }
    return 0;
  }

  TrialResult best;
  bool haveBest = false;

  for( int trial = 0; trial < trialCount; ++trial ){
    map< string, double > startValues = makeStartValues();

    TrialResult result = runFitTrial( ATI, startValues, maxCalls );

    cout << "Trial " << ( trial + 1 ) << "/" << trialCount
         << " likelihood = " << result.likelihood
         << " status = " << result.status
         << " eMatrix = " << result.eMatrixStatus << endl;

    if( !haveBest ||
        ( result.good && !best.good ) ||
        ( result.good == best.good && result.likelihood < best.likelihood ) ){
      best = result;
      haveBest = true;
    }
  }

  cout << endl;
  cout << "Best trial likelihood = " << best.likelihood << endl;
  applyParameterValues( fitManager->parameterManager(), best.finalValues );
  fitManager->migradMinimization();
  syncConfigurationInfo( ATI );
  cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ATI.likelihood() << endl;
  cout << "MINUIT STATUS: " << fitManager->status()
       << "  E-MATRIX STATUS: " << fitManager->eMatrixStatus() << endl;

  ATI.finalizeFit();
  if( !seedFile.empty() ){
    ATI.fitResults()->writeSeed( seedFile );
  }

  return 0;
}
