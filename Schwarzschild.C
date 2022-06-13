#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
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

// Property list:

GYOTO_PROPERTY_START(Schwarzschild,
		     "Schwarzschild space-time")
GYOTO_PROPERTY_BOOL(Schwarzschild, GenericIntegrator, SpecificIntegrator,
		    genericIntegrator,
		    "Which version of the Legacy integrator should be used (specific).")
GYOTO_PROPERTY_DOUBLE(Schwarzchild, DiffTol, difftol,
		      "Tuning parameter for the specific Legacy integrator (0.01).")
GYOTO_PROPERTY_END(Schwarchild, Generic::properties)

// Constructor

Schwarzschild::Schwarzschild():
    Generic(GYOTO_COORDKIND_SPHERICAL, "Schwarzschild"),
    generic_integrator_(false)
{}

// Copy constructor
Schwarzschild * Schwarzschild::clone () const { return new Schwarzschild(*this); }

// Metric in covariant form
void Schwarzschild::gmunu(double g[4][4], const double * pos) const
{
    double r = pos[1];
    double sth = sin(pos[2]);
    double sth2 = sth*sth;
    double r2 = r*r;
    
    int mu, nu;
    for (mu=0; mu<4; ++mu)
        for  (nu=0; nu<4; ++nu)
            g[mu][nu]=g[nu][mu]=0;

    g[0][0] = -(1.-2./r);
    g[1][1] = 1./(1.-2./r); 
    g[2][2] = r2;
    g[3][3] = r2*sth2;  
}

double Schwarzschild::gmunu(const double * pos, int mu, int nu) const
{
    double r = pos[1];
    double sth = sin(pos[2]);
    double sth2 = sth*sth;
    double r2 = r*r;
    
    if ((mu==0) && (nu==0)) return -(1.-2./r);
    if ((mu==1) && (nu==1)) return 1./(1.-2./r);
    if ((mu==2) && (nu==2)) return r2;
    if ((mu==3) && (nu==3)) return r2*sth2;
    
    return 0.;
}

//Metric in contravariant form
void Schwarzschild::gmunu_up(double g[4][4], const double * pos) const
{
    double r = pos[1];
    double sth = sin(pos[2]);
    double sth2 = sth*sth;
    double r2 = r*r;
    
    int mu, nu;
    for (mu=0; mu<4; ++mu)
        for  (nu=0; nu<4; ++nu)
            g[mu][nu]=g[nu][mu]=0;

    g[0][0] = -1./(1.-2./r);
    g[1][1] = (1.-2./r); 
    g[2][2] = 1./r2;
    g[3][3] = 1./(r2*sth2);  
}

double Schwarzschild::gmunu_up(const double * pos, int mu, int nu) const
{
    double r = pos[1];
    double sth = sin(pos[2]);
    double sth2 = sth*sth;
    double r2 = r*r;
    
    if ((mu==0) && (nu==0)) return -1./(1.-2./r);
    if ((mu==1) && (nu==1)) return (1.-2./r);
    if ((mu==2) && (nu==2)) return 1./r2;
    if ((mu==3) && (nu==3)) return 1./(r2*sth2);
    
    return 0.;
}


int Schwarzschild::christoffel(double dst[4][4][4], double const pos[4]) const
{
  int a, mu, nu;
  for (a=0; a<4; ++a)
    for(mu=0; mu<4; ++mu)
      for(nu=0; nu<4; ++nu)
	dst[a][mu][nu]=0.;

  double r = pos[1];
  double sth, cth;
  sincos(pos[2], &sth, &cth);
  double 
    sth2=sth*sth, cth2 = cth*cth,
    ctgth=cth/sth;
  double r2 = r*r;
  double r3 = r2*2;

  dst[0][0][1] = dst[0][1][0] = 1./(r*(r-2.));
  dst[1][0][0] = (r-2.)/r3;
  dst[1][2][2] = -(r-2.);
  dst[1][1][1] = - dst[0][0][1];  
  dst[1][3][3] = -(r-2.)*sth2;
  dst[2][3][3] = -sth*cth;
  dst[2][1][2] = dst[2][2][1] = 1./r;
  dst[3][1][3] = dst[3][3][1] = dst[2][1][2]
  dst[3][2][3] = dst[3][3][2] = ctgth;

  return 0.;
}