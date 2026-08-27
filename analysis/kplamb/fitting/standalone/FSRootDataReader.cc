#include "FSRootDataReader.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

namespace {

bool
treeHasBranches( TTree* tree, const string& prefix, const vector< string >& labels )
{
  if( !tree ) return false;

  const string beamBranches[4] = {
    prefix + "EnPB",
    prefix + "PxPB",
    prefix + "PyPB",
    prefix + "PzPB"
  };

  for( unsigned int i = 0; i < 4; ++i ){
    if( tree->GetBranch( beamBranches[i].c_str() ) == NULL ) return false;
  }

  for( vector< string >::const_iterator lab = labels.begin();
       lab != labels.end();
       ++lab ){

    if( tree->GetBranch( ( prefix + "EnP" + *lab ).c_str() ) == NULL ) return false;
    if( tree->GetBranch( ( prefix + "PxP" + *lab ).c_str() ) == NULL ) return false;
    if( tree->GetBranch( ( prefix + "PyP" + *lab ).c_str() ) == NULL ) return false;
    if( tree->GetBranch( ( prefix + "PzP" + *lab ).c_str() ) == NULL ) return false;
  }

  return true;
}

} // namespace

FSRootDataReader::FSRootDataReader() :
  UserDataReader< FSRootDataReader >(),
  m_inFile( NULL ),
  m_inTree( NULL ),
  m_weightFile( NULL ),
  m_weightTree( NULL ),
  m_eventCounter( 0 ),
  m_expectedFinalState( 0 ),
  m_useMC( false ),
  m_useInlineWeight( false ),
  m_useExternalWeight( false ),
  m_branchPrefix( "" ),
  m_weightBranchName( "Weight" ),
  m_eBeam( 0.0 ),
  m_pxBeam( 0.0 ),
  m_pyBeam( 0.0 ),
  m_pzBeam( 0.0 ),
  m_weight( 1.0 ),
  m_externalWeight( 1.0 )
{
}

FSRootDataReader::FSRootDataReader( const vector< string >& args ) :
  UserDataReader< FSRootDataReader >( args ),
  m_inFile( NULL ),
  m_inTree( NULL ),
  m_weightFile( NULL ),
  m_weightTree( NULL ),
  m_eventCounter( 0 ),
  m_expectedFinalState( 0 ),
  m_useMC( false ),
  m_useInlineWeight( false ),
  m_useExternalWeight( false ),
  m_branchPrefix( "" ),
  m_weightBranchName( "Weight" ),
  m_eBeam( 0.0 ),
  m_pxBeam( 0.0 ),
  m_pyBeam( 0.0 ),
  m_pzBeam( 0.0 ),
  m_weight( 1.0 ),
  m_externalWeight( 1.0 )
{
  assert( args.size() >= 3 );

  TH1::AddDirectory( kFALSE );

  m_expectedFinalState = static_cast< unsigned int >( atoi( args[2].c_str() ) );
  assert( m_expectedFinalState > 0 );
  assert( m_expectedFinalState < Kinematics::kMaxParticles );

  vector< string > trailing;
  for( unsigned int i = 3; i < args.size(); ++i ){
    if( args[i] == "MC" || args[i] == "mc" ){
      m_useMC = true;
    }
    else{
      trailing.push_back( args[i] );
    }
  }

  m_inFile = TFile::Open( args[0].c_str() );
  assert( m_inFile != NULL );

  string treeName = ( args[1].size() ? args[1] : string( "kin" ) );
  m_inTree = dynamic_cast< TTree* >( m_inFile->Get( treeName.c_str() ) );
  assert( m_inTree != NULL );

  if( trailing.size() == 3 ){
    m_useExternalWeight = true;
    m_weightFile = TFile::Open( trailing[0].c_str() );
    assert( m_weightFile != NULL );
    m_weightTree = dynamic_cast< TTree* >( m_weightFile->Get( trailing[1].c_str() ) );
    assert( m_weightTree != NULL );
    m_weightBranchName = trailing[2];
  }
  else{
    assert( trailing.empty() );
  }

  assert( configureBranchLayout() );
  bindMainTreeBranches();

  if( !m_useExternalWeight ){
    if( hasBranch( m_branchPrefix + "Weight" ) ){
      m_useInlineWeight = true;
      m_weightBranchName = m_branchPrefix + "Weight";
      m_inTree->SetBranchAddress( m_weightBranchName.c_str(), &m_weight );
    }
    else if( hasBranch( "Weight" ) ){
      m_useInlineWeight = true;
      m_weightBranchName = "Weight";
      m_inTree->SetBranchAddress( m_weightBranchName.c_str(), &m_weight );
    }
  }
  else{
    assert( m_weightTree->GetBranch( m_weightBranchName.c_str() ) != NULL );
    m_weightTree->SetBranchAddress( m_weightBranchName.c_str(), &m_externalWeight );
  }
}

FSRootDataReader::~FSRootDataReader()
{
  if( m_inFile != NULL ) m_inFile->Close();
  if( m_weightFile != NULL ) m_weightFile->Close();
}

bool
FSRootDataReader::hasBranch( const string& name ) const
{
  return ( m_inTree != NULL && m_inTree->GetBranch( name.c_str() ) != NULL );
}

vector< vector< string > >
FSRootDataReader::candidateLabelOrders( unsigned int nFinalState ) const
{
  vector< vector< string > > candidates;

  if( nFinalState == 3 ){
    candidates.push_back( vector< string >{ "2", "1a", "1b" } );
    candidates.push_back( vector< string >{ "1", "2", "3" } );
    candidates.push_back( vector< string >{ "1", "2a", "2b" } );
    candidates.push_back( vector< string >{ "1", "1a", "1b" } );
  }
  else if( nFinalState == 2 ){
    candidates.push_back( vector< string >{ "1", "2" } );
    candidates.push_back( vector< string >{ "1a", "1b" } );
    candidates.push_back( vector< string >{ "2a", "2b" } );
  }
  else{
    vector< string > labels;
    for( unsigned int i = 1; i <= nFinalState; ++i ){
      labels.push_back( to_string( i ) );
    }
    candidates.push_back( labels );
  }

  return candidates;
}

string
FSRootDataReader::branchName( const string& stem, const string& label ) const
{
  return m_branchPrefix + stem + label;
}

bool
FSRootDataReader::configureBranchLayout()
{
  vector< string > prefixes;
  if( m_useMC ){
    prefixes.push_back( "MC" );
    prefixes.push_back( "" );
  }
  else{
    prefixes.push_back( "" );
    prefixes.push_back( "MC" );
  }

  vector< vector< string > > candidates = candidateLabelOrders( m_expectedFinalState );
  for( vector< string >::const_iterator prefix = prefixes.begin();
       prefix != prefixes.end();
       ++prefix ){

    for( vector< vector< string > >::const_iterator cand = candidates.begin();
         cand != candidates.end();
         ++cand ){

      if( treeHasBranches( m_inTree, *prefix, *cand ) ){
        m_branchPrefix = *prefix;
        m_finalStateLabels = *cand;
        return true;
      }
    }
  }

  cerr << "FSRootDataReader ERROR: could not identify the branch layout for "
       << m_inTree->GetName() << " in " << m_inFile->GetName() << endl;
  cerr << "  expected final-state count = " << m_expectedFinalState << endl;
  cerr << "  looked for either reconstructed or MC-prefixed beam/final-state branches" << endl;
  return false;
}

void
FSRootDataReader::bindMainTreeBranches()
{
  assert( m_finalStateLabels.size() == m_expectedFinalState );

  m_inTree->SetBranchAddress( ( m_branchPrefix + "EnPB" ).c_str(), &m_eBeam );
  m_inTree->SetBranchAddress( ( m_branchPrefix + "PxPB" ).c_str(), &m_pxBeam );
  m_inTree->SetBranchAddress( ( m_branchPrefix + "PyPB" ).c_str(), &m_pyBeam );
  m_inTree->SetBranchAddress( ( m_branchPrefix + "PzPB" ).c_str(), &m_pzBeam );

  for( unsigned int i = 0; i < m_finalStateLabels.size(); ++i ){
    const string& label = m_finalStateLabels[i];

    m_inTree->SetBranchAddress( branchName( "EnP", label ).c_str(), &m_e[i+1] );
    m_inTree->SetBranchAddress( branchName( "PxP", label ).c_str(), &m_px[i+1] );
    m_inTree->SetBranchAddress( branchName( "PyP", label ).c_str(), &m_py[i+1] );
    m_inTree->SetBranchAddress( branchName( "PzP", label ).c_str(), &m_pz[i+1] );
  }
}

void
FSRootDataReader::resetSource()
{
  cout << "Resetting source " << m_inTree->GetName()
       << " in " << m_inFile->GetName() << endl;
  m_eventCounter = 0;
}

Kinematics*
FSRootDataReader::getEvent()
{
  if( m_eventCounter >= static_cast< unsigned int >( m_inTree->GetEntries() ) ){
    return NULL;
  }

  const unsigned int entry = m_eventCounter++;
  m_inTree->GetEntry( entry );

  if( m_useExternalWeight ){
    if( entry < static_cast< unsigned int >( m_weightTree->GetEntries() ) ){
      m_weightTree->GetEntry( entry );
    }
    else{
      m_externalWeight = 1.0;
    }
  }

  vector< TLorentzVector > particleList;
  particleList.push_back( TLorentzVector( m_pxBeam, m_pyBeam, m_pzBeam, m_eBeam ) );

  for( unsigned int i = 0; i < m_expectedFinalState; ++i ){
    particleList.push_back( TLorentzVector( m_px[i+1], m_py[i+1], m_pz[i+1], m_e[i+1] ) );
  }

  double weight = 1.0;
  if( m_useExternalWeight ) weight = m_externalWeight;
  else if( m_useInlineWeight ) weight = m_weight;

  return new Kinematics( particleList, weight );
}

unsigned int
FSRootDataReader::numEvents() const
{
  return ( m_inTree ? static_cast< unsigned int >( m_inTree->GetEntries() ) : 0 );
}
