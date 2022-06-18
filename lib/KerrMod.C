#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrMod.h"
#include "GyotoWorldline.h"
#include "GyotoError.h"
#include "GyotoProperty.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std ; 
using namespace Gyoto ;
using namespace Gyoto::Metric ;

GYOTO_PROPERTY_START(KerrMod,
		     "Metric around a rotating black-hole, in spherical Boyer-Lindquist coordinates.")
GYOTO_PROPERTY_DOUBLE(KerrMod, Spin, spin,
		      "Spin parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(KerrMod, Alpha, alpha,
		      "Alpha parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(KerrMod, Sigma, sigma,
		      "Beta parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(KerrMod, Kappa, kappa,
		      "Kappa parameter (adimensioned, 0).")
GYOTO_PROPERTY_DOUBLE(KerrMod, HorizonSecurity, horizonSecurity,
		      "Thickness of sink layer around horizon (geometrical units, 0.01).")
GYOTO_PROPERTY_BOOL(KerrMod, GenericIntegrator, SpecificIntegrator,
		    genericIntegrator,
		    "Which version of the Legacy integrator should be used (specific).")
GYOTO_PROPERTY_DOUBLE(KerrMod, DiffTol, difftol,
		      "Tuning parameter for the specific Legacy integrator (0.01).")
GYOTO_PROPERTY_END(KerrMod, Generic::properties)



#define GYOTO_MIN_THETA 1e-10

					       
KerrMod::KerrMod() :
  Generic(GYOTO_COORDKIND_SPHERICAL, "KerrMod"),
  spin_(0.), a2_(0.), a3_(0.), a4_(0.),
  alpha_(0.), al2_(0.), al3_(0.), al4_(0.),
  sigma_(0.), sg2_(0.), sg3_(0.), sg3_(0.),
  kappa_(0.), kp2_(0.), kp3_(0.), kp4_(0.),
  difftol_(GYOTO_KERRMOD_DEFAULT_DIFFTOL),
  rsink_(2.+GYOTO_KERR_HORIZON_SECURITY),
  drhor_(GYOTO_KERR_HORIZON_SECURITY), generic_integrator_(false)
{}

// default copy constructor should be fine 
KerrMod * KerrMod::clone () const { return new KerrMod(*this); }

// Mutators
void KerrMod::spin(const double a) {
  spin_=a;
  a2_=spin_*spin_;
  a3_=a2_*spin_;
  a4_=a2_*a2_;
  rsink_=1.+sqrt(1.-a2_)+drhor_; //CAMBIAR EL HORIZONTE
  tellListeners();
}

void KerrMod::alpha(const double al) {
  alpha_=al;
  al2_=alpha_*alpha_;
  al3_=al2_*alpha_;
  al4_=al2_*al2_;
  rsink_=1.+sqrt(1.-a2_)+drhor_; //CAMBIAR EL HORIZONTE
  tellListeners();
}

void KerrMod::sigma(const double sg) {
  sigma_=sg;
  sg2_=sigma_*sigma_;
  sg3_=sg2_*sigma_;
  sg4_=sg2_*sg2_;
  rsink_=1.+sqrt(1.-a2_)+drhor_; //CAMBIAR EL HORIZONTE
  tellListeners();
}

void KerrMod::kappa(const double kp) {
  kappa_=kp;
  kp2_=kappa_*kappa_;
  kp3_=kp2_*kappa_;
  kp4_=kp2_*kp2_;
  rsink_=1.+sqrt(1.-a2_)+drhor_; //CAMBIAR EL HORIZONTE
  tellListeners();
}

// Accessors
double KerrMod::spin() const { return spin_ ; }
double KerrMod::alpha() const { return alpha_ ; }
double KerrMod::sigma() const { return sigma_ ; }
double KerrMod::kappa() const { return kappa_ ; }

double KerrMod::difftol() const { return difftol_;}
void KerrMod::difftol(double t) {difftol_=t;}

void KerrMod::horizonSecurity(const double drhor) {
  drhor_=drhor;
  rsink_=1.+sqrt(1.-a2_)+drhor_; //CAMBIAR HORIZONTE
  tellListeners();
}
double KerrMod::horizonSecurity() const {return drhor_; }

void KerrMod::genericIntegrator(bool t)
{
  generic_integrator_=t;
  tellListeners();
}
bool KerrMod::genericIntegrator() const {return generic_integrator_;}

int KerrMod::isStopCondition(double const * const coord) const {
  return coord[1] < rsink_ ;
}

//Prograde marginally stable orbit, SE DEBE MODIFICAR
double KerrMod::getRms() const {
  double aa=spin_;
  double  z1 = 1. + pow((1. - a2_),1./3.)*(pow((1. + aa),1./3.) + pow((1. - aa),1./3.)); 
  double  z2 = pow(3.*a2_ + z1*z1,1./2.);

  return (3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),1./2.));
}

//Prograde marginally bound orbit SE DEBE CAMBIAR
double KerrMod::getRmb() const {
  return 2.-spin_+2.*sqrt(1.-spin_);
}

//SE DEBE CAMBIAR
double KerrMod::getSpecificAngularMomentum(double rr) const {
  // this is l = -u_phi/u_t for a circular equatorial 4-velocity
  double aa=spin_, sqrtr=sqrt(rr);
  return (rr*rr-2.*aa*sqrtr+aa*aa)/(pow(rr,1.5)-2.*sqrtr+aa);
}

double KerrMod::getPotential(double const pos[4], double l_cst) const {
  // this is W = -ln(|u_t|) for a circular equatorial 4-velocity
  double  gtt = gmunu(pos,0,0);
  double  gtp = gmunu(pos,0,3);
  double  gpp = gmunu(pos,3,3);
  double  Omega = -(gtp + l_cst * gtt)/(gpp + l_cst * gtp) ;
  
  double  W = 0.5 * log(abs(gtt + 2. * Omega * gtp + Omega*Omega * gpp)) 
    - log(abs(gtt + Omega * gtp)) ;
  return  W ;
}

void MyKerr::gmunu(double g[4][4], const double * pos) const {
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double r2=r*r;
  double sigma=r2+a2_*cth2;
  double delta=r2-2.*r+a2_;
  double M=al2_*(1.+ sigma_/r + kappa_/r2);
  
  int mu, nu;
  for (mu=0; mu<4; ++mu)
    for (nu=0; nu<4; ++nu)
      g[mu][nu]=0.;

  g[0][0] = -1.+2.*r/sigma;
  g[1][1] = sigma/delta;
  g[2][2] = sigma;
  g[3][3] = (r2+a2_+2.*r*a2_*sth2/sigma)*sth2;
  g[0][3] = g[3][0] = (M-2*r/sigma)*spin_*sth2;
}

double MyKerr::gmunu(const double * pos, int mu, int nu) const {
  double r = pos[1];
  double sth2, cth2;
  sincos(pos[2], &sth2, &cth2);
  sth2*=sth2; cth2*=cth2;
  double r2=r*r;
  double sigma=r2+a2_*cth2;
  double delta=r2-2.*r+a2_;
  double M=al2_*(1.+ sigma_/r + kappa_/r2);

  if ((mu==0) && (nu==0)) return -(1.-2.*r/sigma);
  if ((mu==1) && (nu==1)) return sigma/delta;
  if ((mu==2) && (nu==2)) return sigma;
  if ((mu==3) && (nu==3))
    return (r2+a2_+2.*r*a2_*sth2/sigma)*sth2;
  if (((mu==0) && (nu==3)) || ((mu==3) && (nu==0))){
    return (M-2*r/sigma)*spin_*sth2;
  }

  return 0.;
} 



