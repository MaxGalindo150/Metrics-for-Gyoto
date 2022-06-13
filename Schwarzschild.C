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

// Metric in covarint form
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

//Metric in contravarint form
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

    g[0][0] = -1./(1.-2./r);
    g[1][1] = (1.-2./r); 
    g[2][2] = 1./r2;
    g[3][3] = 1./(r2*sth2);  
}

double Schwarzschild::gmunu(const double * pos, int mu, int nu) const
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
