/**
 * \file GyotoPageThorneDiskFree.h
 * \brief A geometrically thin, optically thick disk
 *
 *  This class describes a disk contained in the z=0 (equatorial)
 *  plane, extending from r=r_ISCO to r=infinity.  The flux emitted
 *  at radius r is given by Page & Thorne (1974, ApJ 191:499,
 *  Eqs. 11b, 14, 15).
 *
 *  Only bolometric flux is implemented (as quantity User4), no
 *  spectral resolution.
 */



#ifndef __GyotoPageThorneDiskFree_H_ 
#define __GyotoPageThorneDiskFree_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class PageThorneDiskFree; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>
#include <GyotoBlackBodySpectrum.h>

/**
 * \class Gyoto::Astrobj::PageThorneDiskFree
 * \brief Geometrically thin disk in Kerr metric
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane, extending from r=r_ISCO to r=infinity.  The flux emitted
 *   at radius r is given by Page & Thorne (1974, ApJ 191:499,
 *   Eqs. 11b, 14, 15).
 *
 *   The metric must be either KerrBL or KerrKS. Emission, Spectrum
 *   and BinSpectrum are <STRONG>not</STRONG> provide, the only
 *   intensity provided is provided, as quantity User4 and it is the
 *   default quantity returned if nothing is requested. The other
 *   quantities implemented in ThinDisk are also provided.
 *
 */
class Gyoto::Astrobj::PageThorneDiskFree
: public Astrobj::ThinDisk,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PageThorneDiskFree>;
 private:
  double aa_; ///< Generic::gg_ spin parameter, monitored by tell()
  double aa2_; ///< aa_<SUP>2</SUP>
  double x0_; ///< Value cached for bolometricEmission()
  double x1_; ///< Value cached for bolometricEmission()
  double x2_; ///< Value cached for bolometricEmission()
  double x3_; ///< Value cached for bolometricEmission()
  double mdot_; ///< accretion rate (for BB spectrum computation)
  bool uniflux_; ///< Flag for uniform flux = 1
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< disk black body
  ///< emission law
  
  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

  PageThorneDiskFree(); ///< Standard constructor
  
  PageThorneDiskFree(const PageThorneDiskFree& ) ;///< Copy constructor
  virtual PageThorneDiskFree* clone () const; ///< Cloner
  
  virtual ~PageThorneDiskFree() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  using ThinDisk::metric;
  virtual void metric(SmartPointer<Metric::Generic>);
  ///< Set metric, checking that it is either KerrBL or KerrKS

  /// Set #mdot_ to v
  void mdot(double v);
  double mdot() const;
  void uniFlux(bool t) ;
  bool uniFlux() const ;

 private:
  virtual void updateSpin() ;
  ///< Get spin from metric, which must be KerrBL or KerrKS

 public:
  using ThinDisk::emission;
  /**
   * \brief Not implemented
   * Throws a Gyoto::Error
   */
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph, double const c_obj[8]=NULL) const;

  /**
   * \brief Bolometric emission
   *
   * Similar to Generic::emission(), but bolometric.
   */
  virtual double bolometricEmission(double nuem, double dsem,
				    double const c_obj[8]) const;

  /**
   * \brief 
   * processHitQuantities fills the requested data in Impact. For
   * PageThorneDiskFree, only fill User4, which corresponds to bolometric
   * intensity.
   */
  virtual void processHitQuantities(Photon* ph, state_t const &coord_ph_hit,
                                   double const *coord_obj_hit, double dt,
                                   Astrobj::Properties* data) const;

  Gyoto::Quantity_t getDefaultQuantities();

  // Hook::Listener API //
 public:
  /**
   * \brief Update PageThorneDiskFree::aa_
   *
   * Calls updateSpin().
   *
   * See Hook::Listener::tell()
   */
  virtual void tell(Gyoto::Hook::Teller *msg);

};

#endif
