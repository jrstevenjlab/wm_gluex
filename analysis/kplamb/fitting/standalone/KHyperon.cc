#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"

#include "IUAmpTools/Kinematics.h"
#include "KHyperon.h"

KHyperon::KHyperon( const vector< string >& args ) :
  UserAmplitude< KHyperon >( args )
{
  assert( args.size() == 8 || args.size() == 9 );

  alpha  = AmpParameter( args[0] );
  Sigma  = AmpParameter( args[1] );
  Ox     = AmpParameter( args[2] );
  P      = AmpParameter( args[3] );
  T      = AmpParameter( args[4] );
  Oz     = AmpParameter( args[5] );

  polAngle = AmpParameter( args[6] );

  // need to register any free parameters so the framework knows about them
  registerParameter( alpha );
  registerParameter( Sigma );
  registerParameter( Ox );
  registerParameter( P );
  registerParameter( T );
  registerParameter( Oz );

  registerParameter( polAngle );

  // Two possibilities to initialize this amplitude:
  // 1: 8 arguments, fixed polarization
  //    Usage: amplitude <reaction>::<sum>::<ampName> KHyperon alpha Sigma Ox P T Oz <polAngle> <polFraction>
  if(args.size() == 8) {
    polFraction = atof(args[7].c_str());
    cout << "Fitting with constant polarization " << polFraction << endl;
  }
  // 2: 8 arguments, read polarization from histogram <hist> in file <rootFile>
  //    Usage: amplitude <reaction>::<sum>::<ampName> KHyperon alpha Sigma Ox P T Oz <polAngle> <rootFile> <hist>
  else if(args.size() == 9) {
    polFraction = 0.;
    TFile* f = new TFile( args[7].c_str() );
    polFrac_vs_E = (TH1D*)f->Get( args[8].c_str() );
    assert( polFrac_vs_E != NULL );
    cout << "Fitting with polarization from " << polFrac_vs_E->GetName() << endl;
  }
}


complex< GDouble >
KHyperon::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {

  GDouble cosThetaX = userVars[kCosThetaX];
  GDouble cosThetaY = userVars[kCosThetaY];
  GDouble cosThetaZ = userVars[kCosThetaZ];
  GDouble phi = polAngle*0.017453293 + userVars[kPhi]; // rotate Phi (in rad)
  GDouble Pgamma = userVars[kPgamma];

  // CLAS paper intensity formulation (DOI 10.1103/physrevc.93.065201)
  GDouble I = 1.0 + alpha*cosThetaY*P;
  I -= Pgamma*cos(2.0*phi) * (Sigma + alpha*cosThetaY*T);
  I += Pgamma*sin(2.0*phi) * alpha * (cosThetaX*Ox + cosThetaZ*Oz);

  return complex< GDouble > ( sqrt(fabs(I)) );
}

void
KHyperon::calcUserVars( GDouble** pKin, GDouble* userVars ) const {

  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector k  ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] );
  TLorentzVector y1 ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] );
  TLorentzVector y2 ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] );

  TLorentzVector hyperon = y1 + y2;
  TLorentzRotation hyperonBoost( -hyperon.BoostVector() );

  TLorentzVector beam_hyperon = hyperonBoost * beam; // beam photon in hyperon rest frame
  TLorentzVector k_hyperon = hyperonBoost * k;       // kaon in hyperon rest frame
  TLorentzVector y1_hyperon = hyperonBoost * y1;     // proton in hyperon rest frame

  // normal to the production plane (formed by beam and kaon)
  TVector3 y = (beam.Vect().Unit().Cross(-k.Vect().Unit())).Unit();

  // choose helicity frame: z-axis opposite kaon in hyperon rest frame
  TVector3 z = -1. * k_hyperon.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles( (y1_hyperon.Vect()).Dot(x),
                   (y1_hyperon.Vect()).Dot(y),
                   (y1_hyperon.Vect()).Dot(z) );

  userVars[kCosThetaX] = cos(y1_hyperon.Vect().Angle(x));
  userVars[kCosThetaY] = cos(y1_hyperon.Vect().Angle(y));
  userVars[kCosThetaZ] = cos(y1_hyperon.Vect().Angle(z));

  TVector3 eps(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees
  userVars[kPhi] = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble Pgamma;
  if(polFraction > 0.) { // for fitting with constant polarization
    Pgamma = polFraction;
  }
  else{
    int bin = polFrac_vs_E->GetXaxis()->FindBin(pKin[0][0]);
    if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
      Pgamma = 0.;
    }
    else Pgamma = polFrac_vs_E->GetBinContent(bin);
  }
  userVars[kPgamma] = Pgamma;
}

#ifdef GPU_ACCELERATION
void
KHyperon::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUKHyperon_exec( dimGrid, dimBlock, GPU_AMP_ARGS, Sigma, Ox, P, T, Oz,
                           polAngle );
}

#endif //GPU_ACCELERATION
