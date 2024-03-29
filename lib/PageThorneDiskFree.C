

#include "GyotoPhoton.h"
#include "GyotoPageThorneDiskFree.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoMod.h"
#include "GyotoKerrKS.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>
#include <cstring>
#include <time.h> 

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(PageThorneDiskFree)
GYOTO_PROPERTY_DOUBLE(PageThorneDiskFree, Mdot, mdot)
GYOTO_PROPERTY_BOOL(PageThorneDiskFree, UniFlux, NonUniFlux, uniFlux)
GYOTO_PROPERTY_END(PageThorneDiskFree, ThinDisk::properties)

void PageThorneDiskFree::mdot(double v) {
  mdot_=v;
}
double PageThorneDiskFree::mdot() const { return mdot_; }

void PageThorneDiskFree::uniFlux(bool t) {uniflux_=t;}
bool PageThorneDiskFree::uniFlux() const {return uniflux_;}

PageThorneDiskFree::PageThorneDiskFree() :
  ThinDisk("PageThorneDiskFree"), aa_(0.), aa2_(0.),
  x0_(0.), x1_(0.), x2_(0.), x3_(0.), mdot_(1.),
  uniflux_(0), spectrumBB_(NULL)
{
  if (debug()) cerr << "DEBUG: PageThorneDiskFree Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

PageThorneDiskFree::PageThorneDiskFree(const PageThorneDiskFree& o) :
  ThinDisk(o), aa_(o.aa_), aa2_(o.aa2_),
  x0_(o.x0_), x1_(o.x1_), x2_(o.x2_), x3_(o.x3_),
  mdot_(o.mdot_), uniflux_(o.uniflux_), 
  spectrumBB_(NULL)
{
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
  if (gg_) gg_->hook(this);
}
PageThorneDiskFree* PageThorneDiskFree::clone() const
{ return new PageThorneDiskFree(*this); }

bool PageThorneDiskFree::isThreadSafe() const {
  // spectrumBB_ is not a Property
  return ThinDisk::isThreadSafe()
    && (!spectrumBB_ || spectrumBB_->isThreadSafe());
}

PageThorneDiskFree::~PageThorneDiskFree() {
  GYOTO_DEBUG<<endl;
  if (gg_) gg_->unhook(this);
}

void PageThorneDiskFree::updateSpin() {
  if (!gg_) return;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    aa_ = static_cast<SmartPointer<Metric::Mod> >(gg_) -> spin();
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    aa_ = static_cast<SmartPointer<Metric::KerrKS> >(gg_) -> spin();
    break;
  default:
    GYOTO_ERROR("PageThorneDiskFree::getSpin(): unknown COORDKIND");
  }
  aa2_=aa_*aa_;
  double z1 =1.+pow((1.-aa2_),1./3.)*(pow((1.+ aa_),1./3.)+pow((1.-aa_),1./3.));
  double z2 = pow(3.*aa2_ + z1*z1,0.5);
  double acosaao3= acos(aa_)/3.;

  x0_ = sqrt((3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),0.5)));
  x1_ = 2.*cos(acosaao3 - M_PI/3.);
  x2_ = 2.*cos(acosaao3 + M_PI/3.); 
  x3_ = -2.*cos(acosaao3);

  if (rin_==0.) rin_=(3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
}

void PageThorneDiskFree::metric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kin = gg->kind();
  if (kin != "Mod" && kin != "Mod")
    GYOTO_ERROR
      ("PageThorneDiskFree::metric(): metric must be KerrBL or KerrKS");
  ThinDisk::metric(gg);
  updateSpin();
  gg->hook(this);
}

double PageThorneDiskFree::emission(double nu_em, double dsem,
				    state_t const &,
				    double const coord_obj[8]) const{
  double Ibolo=bolometricEmission(nu_em,dsem,coord_obj);
  /*
    From Ibolo, find T, then Bnu(T)
   */
  /*double mass=gg_->mass()*1e3; // in cgs
  double c6=GYOTO_C_CGS*GYOTO_C_CGS*GYOTO_C_CGS
    *GYOTO_C_CGS*GYOTO_C_CGS*GYOTO_C_CGS;
  double g2m2=GYOTO_G_CGS*GYOTO_G_CGS*mass*mass;
  Ibolo*=mdot_*c6/g2m2; // Ibolo in cgs*/ // --> this is now done in boloEm 
  //F = sigma * T^4 (and F=pi*I)
  double TT=pow(Ibolo*M_PI/GYOTO_STEFANBOLTZMANN_CGS,0.25);
  spectrumBB_->temperature(TT);
  double Iem=(*spectrumBB_)(nu_em); // in SI
  //cout << "r T nu Iem = " << coord_obj[1] << " " << TT << " " << nu_em << " " << Iem << endl;
  if (Iem < 0.) GYOTO_ERROR("In PageThorneDiskFree::emission"
			   " blackbody emission is negative!");
  return Iem;
}

Quantity_t PageThorneDiskFree::getDefaultQuantities() {
  return GYOTO_QUANTITY_USER4;
}

double PageThorneDiskFree::bolometricEmission(double /* nuem */, double dsem,
				    double const coord_obj[8]) const{
  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  // Important remark: this emision function gives I(r),
  // not I_nu(r). And I(r)/nu^4 is conserved.
  //  cout << "r hit= " << coord_obj[1] << endl;
  if (uniflux_) return 1;
  double xx;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    xx=sqrt(coord_obj[1]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xx=pow(coord_obj[1]*coord_obj[1]+coord_obj[2]*coord_obj[2]-aa2_, 0.25);
    break;
  default:
    GYOTO_ERROR("Unknown coordinate system kind");
    xx=0;
  }
  
  // the above formula assume M=1 (x=sqrt(r/M)=sqrt(r))
  double x2=xx*xx;
  
  double ff=
    3./(2.)*1./(xx*xx*(xx*xx*xx-3.*xx+2.*aa_))
    *( 
      xx-x0_-3./2.*aa_*log(xx/x0_)
      -3.*(x1_-aa_)*(x1_-aa_)/(x1_*(x1_-x2_)*(x1_-x3_))*log((xx-x1_)
							    /(x0_-x1_)) 
      -3.*(x2_-aa_)*(x2_-aa_)/(x2_*(x2_-x1_)*(x2_-x3_))*log((xx-x2_)
							    /(x0_-x2_)) 
      -3.*(x3_-aa_)*(x3_-aa_)/(x3_*(x3_-x1_)*(x3_-x2_))*log((xx-x3_)
							    /(x0_-x3_))
       );
  // f of Page&Thorne 1974 eq 15n, in units of 1/M

  double Iem=ff/(4.*M_PI*M_PI*x2); // natural-units value
  /*
    Assuming isotropic emission: 
    flux at r = (I at r)* \int cos\theta dOmega = pi*I
    thus intensity is only: 1/pi * flux;
    the flux F(r) is given by eq 11b in Page-Thorne,
    so it is F(r) = Mdot/(4pi) * 1/r * f 
    (see eq 15d for why the exp term is 1/r).
    So finally, indeed, Iem, the bolometric intensity, is given
    by ff/(4.*M_PI*M_PI*x2)
  */
  if (gg_->mass()!=1. and mdot_!=1.){ // non-default values for M and Mdot
    // the cgs value of I is c^6/G^2*Mdot/M^2 * Iem
    double mass=gg_->mass()*1e3; // in cgs
    double c6=GYOTO_C_CGS*GYOTO_C_CGS*GYOTO_C_CGS
      *GYOTO_C_CGS*GYOTO_C_CGS*GYOTO_C_CGS;
    double g2m2=GYOTO_G_CGS*GYOTO_G_CGS*mass*mass;
    Iem*=mdot_*c6/g2m2; // this incorporates in particular
    // the 1/M factor in Page-Thorne eq 15n which is not used
    // in the definition of ff above
  }

  if (flag_radtransf_) Iem *= dsem;
  GYOTO_DEBUG_EXPR(Iem);
  
  return Iem*GYOTO_INU_CGS_TO_SI; // in SI
  
}

void PageThorneDiskFree::processHitQuantities(Photon* ph, state_t const &coord_ph_hit,
				     double const *coord_obj_hit, double dt,
				     Properties* data) const {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  /*
      NB: freqObs is the observer's frequency chosen in
      Screen::getRayCoord for the actual computation of the geodesic ;
      the physical value of nuobs will be used in spectrum
      computations by resorting to the xml specifications of the user
      (see below) ; this freqObs is used to transform the null
      worldline parameter dlambda (see below)
  */
  double freqObs=ph->freqObs(); // this is a useless quantity, always 1
  SmartPointer<Spectrometer::Generic> spr = ph -> spectrometer();
  size_t nbnuobs = spr() ? spr -> nSamples() : 0 ;
  double const * const nuobs = nbnuobs ? spr -> getMidpoints() : NULL;
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(&coord_ph_hit[0],coord_obj_hit+4,
				    &coord_ph_hit[4]);// / 1.; 
                                       //this is nu_em/nu_obs
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  if (uniflux_) ggred=1.;
  double dsem = dlambda*ggredm1; // *1.
  double inc =0.;
  if (data) {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "data requested. " 
	      << ", ggredm1=" << ggredm1
	      << ", ggred=" << ggred
	      << endl;
#endif

    if (data->redshift) {
      *data->redshift=ggred;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->redshift);
#endif
    }
    if (data->time) {
      *data->time=coord_ph_hit[0];
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->time);
#endif
    }
    if (data->impactcoords) {
      if (coord_ph_hit.size() > 8) GYOTO_ERROR("ImpactCoords is incompatible with parallel transport");
      memcpy(data->impactcoords, coord_obj_hit, 8 * sizeof(double));
      memcpy(data->impactcoords+8, &coord_ph_hit[0], 8 * sizeof(double));
    }
#if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "dlambda = (dt="<< dt << ")/(tdot="<< coord_ph_hit[4]
		<< ") = " << dlambda << ", dsem=" << dsem << endl;
#endif
    if (data->intensity) GYOTO_ERROR("unimplemented");
    if (data->user4) {
      inc = (bolometricEmission(freqObs*ggredm1, dsem, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* ggred*ggred*ggred*ggred; // I/nu^4 invariant
      *data->user4 += inc;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->user4);
#endif

    }
    if (data->binspectrum) GYOTO_ERROR("unimplemented");
    if (data->spectrum)  {
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	double nuem=nuobs[ii]*ggredm1;
	inc = (emission(nuem, dsem, coord_ph_hit, coord_obj_hit))
	  * (ph -> getTransmission(size_t(-1)))
	  * ggred*ggred*ggred; // Inu/nu^3 invariant
	data->spectrum[ii*data->offset] += inc;
	//cout << "in spec stored= " << ggred << " " << inc << endl;
      }
    }
    /* update photon's transmission */
    ph -> transmit(size_t(-1),
		   transmission(freqObs*ggredm1, dsem,coord_ph_hit, coord_obj_hit));
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "NO data requested!" << endl;
#   endif
  }
}

void PageThorneDiskFree::tell(Hook::Teller* msg) {
  if (msg==gg_) updateSpin();
}
