#if !defined(FSROOTDATAREADER)
#define FSROOTDATAREADER

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/UserDataReader.h"

#include <string>
#include <vector>

class TFile;
class TTree;

using std::string;
using std::vector;

class FSRootDataReader : public UserDataReader< FSRootDataReader >
{

public:

  FSRootDataReader();
  FSRootDataReader( const vector< string >& args );
  ~FSRootDataReader();

  string name() const { return "FSRootDataReader"; }

  virtual Kinematics* getEvent();
  virtual void resetSource();
  virtual unsigned int numEvents() const;

  virtual bool hasWeight() { return m_useInlineWeight || m_useExternalWeight; }

private:

  bool hasBranch( const string& name ) const;
  vector< vector< string > > candidateLabelOrders( unsigned int nFinalState ) const;
  bool configureBranchLayout();
  void bindMainTreeBranches();

  string branchName( const string& stem, const string& label ) const;

  TFile* m_inFile;
  TTree* m_inTree;
  TFile* m_weightFile;
  TTree* m_weightTree;

  unsigned int m_eventCounter;
  unsigned int m_expectedFinalState;

  bool m_useMC;
  bool m_useInlineWeight;
  bool m_useExternalWeight;

  string m_branchPrefix;
  vector< string > m_finalStateLabels;
  string m_weightBranchName;

  double m_eBeam;
  double m_pxBeam;
  double m_pyBeam;
  double m_pzBeam;

  double m_e[ Kinematics::kMaxParticles ];
  double m_px[ Kinematics::kMaxParticles ];
  double m_py[ Kinematics::kMaxParticles ];
  double m_pz[ Kinematics::kMaxParticles ];

  double m_weight;
  double m_externalWeight;
};

#endif
