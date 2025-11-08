//
//  nonbond_functions.c
//
//  Library of nonbond functions
//  used in MDG4A chip
//
//  Copyright (c) 2016 by Makoto Taiji/RIKEN QBiC
//  Copyright (c) 2019 by Makoto Taiji/RIKEN BDR
//  taiji@riken.jp
//
//  Force returns absolute magnitude of force.
//  Function evaluator should return force value/r.
//
package rial.table

import scala.math._
import java.lang.Math
import org.apache.commons.math3.special.Erf._

private[rial] object NonbondFunctions {

  ////////////////////////////////////////////////////////////////////////
  // Ewald normal

  def ewaldPotentialCore (alpha : Double) (r: Double) = {
    erfc(alpha*r)
  }

  def ewaldPotential (alpha : Double) (r: Double) = {
    erfc(alpha*r)/r
  }

  def ewaldForceCore (alpha : Double) (r: Double) = {
    val x=alpha*r
    erfc(x) + 2.0*x*exp(-x*x)/sqrt(Pi)
  }

  def ewaldForce (alpha : Double) (r: Double) = {
    ewaldForceCore(alpha)(r)/(r*r)
  }

  def ewaldPotentialExcl (alpha : Double) (r: Double) = {
    -erf(alpha*r)/r
  }

  def ewaldForceExcl (alpha : Double) (r: Double) = {
    val x=alpha*r
    -erf(x)/(r*r) + 2.0*alpha*exp(-x*x)/sqrt(Pi)/r
  }

  def ewaldPotentialExcl14 (alpha : Double) (r: Double) = {
    ewaldPotential(alpha)(r) - 0.5/r
  }

  def ewaldForceExcl14 (alpha : Double) (r: Double) = {
   ewaldForce(alpha)(r) - 0.5/(r*r)
  }

  ////////////////////////////////////////////////////////////////////////
  // Zero dipole

  private def zdUtil (r : Double, rc: Double, alpha : Double) = {
    ewaldPotential(alpha)(r)+0.5*ewaldForce(alpha)(rc)*r*r
  }

  def zeroDipolePotential(rc : Double) (alpha : Double) (r : Double) = {
    if (r>=rc) 0.0
    else zdUtil(r, rc, alpha) - zdUtil(rc, rc, alpha)
  }

  def zeroDipoleForce (rc : Double) (alpha : Double) (r : Double) = {
    if (r>=rc) 0.0
    ewaldForce(alpha)(r) - ewaldForce(alpha)(rc)
  }

  ////////////////////////////////////////////////////////////////////////
  // Soft core (linear exploration from 1/r)

  def coulombSoftcorePotential (r0 : Double) (r: Double) = {
    if (r>=r0) 1.0/r
    else {
      val grad = 1.0/(r0*r0)
      val dfdr = -2.0*pow(r0,-3.0)
      1.0/r0 - grad*(r-r0) - 0.5*dfdr*(r-r0)*(r-r0)
    }
  }

  def coulombSoftcoreForce (r0 : Double) (r: Double) = {
    if (r>=r0) pow(r, -2.0)
    else {
      val dfdr = -2.0*pow(r0,-3.0)
      1.0/(r0*r0) + dfdr * (r-r0)
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Soft core (linear exploration from Ewald)

  def ewaldSoftcorePotential (alpha : Double) (r0 : Double) (r: Double) = {
    val y = ewaldPotential(alpha)(r)
    if (r>=r0) y
    else y - 1.0/r + coulombSoftcorePotential(r0)(r)
  }

  def ewaldSoftcoreForce (alpha : Double) (r0 : Double) (r: Double) = {
    val y = ewaldForce(alpha)(r)
    if (r>=r0) y;
    else y - pow(r,-2.0) + coulombSoftcoreForce(r0)(r)
  }

  def ewaldSoftcorePotentialExcl (alpha : Double) (r0 : Double) (r: Double) = {
    ewaldPotentialExcl(alpha)(r)
  }

  def ewaldSoftcoreForceExcl (alpha : Double) (r0 : Double) (r: Double) = {
    ewaldForceExcl(alpha)(r)
  }

  def ewaldSoftcorePotentialExcl14 (alpha : Double) (r0 : Double) (r: Double) = {
    if (r>=r0) ewaldPotentialExcl14(alpha)(r)
    else ewaldPotentialExcl(alpha)(r) + 0.5*coulombSoftcorePotential(r0)(r)
  }

  def ewaldSoftcoreForceExcl14 (alpha : Double) (r0 : Double) (r: Double) = {
    if (r>=r0) ewaldForceExcl14(alpha)(r)
    else ewaldForceExcl(alpha)(r) + 0.5*coulombSoftcoreForce(r0)(r)
  }

  ////////////////////////////////////////////////////////////////////////
  // Soft core (potential linear exploration from 1/r)

  def coulombCfSoftcorePotential (r0 : Double) (r: Double) = {
    if (r>=r0) 1.0/r
    else 1.0/r0 - (r-r0)/(r0*r0)
  }

  def coulombCfSoftcoreForce (r0 : Double) (r: Double) = {
    if (r>=r0) 1.0/(r*r)
    else       1.0/(r0*r0)
  }

  ////////////////////////////////////////////////////////////////////////
  // Soft core (force propotional to r below r0)

  private def lfSoftcorePotential (r0: Double) (f0r0:Double) (phi0:Double) (r:Double) = {
    phi0 + 0.5*f0r0*(r0*r0 - r*r)
  }

  private def lfSoftcoreForce (f0r0:Double) (r:Double) = {
    f0r0*r
  }

  def coulombLfSoftcorePotential (r0 : Double) (r: Double) = {
    if (r>=r0) 1.0/r
    else {
      val pot = 1.0/r0
      val grad = pot*pot*pot
      lfSoftcorePotential(r0)(grad)(pot)(r)
    }
  }

  def coulombLfSoftcoreForce (r0 : Double) (r: Double) = {
    if (r>=r0) 1.0/(r*r)
    else lfSoftcoreForce(1.0/(r0*r0*r0))(r)
  }

  ////////////////////////////////////////////////////////////////////////
  // Soft core (Gromacs)

  // Actually, rA = (alpha sigma^6 lambda^p + r^6)^(1/6)
  // Here, alpha=lambda=1, so give sigma to include alpha/lambda
  // Gromacs suggests alpha=0.7 for p=1 or alpha=0.5 for p=2.
  // In the latter case, alpha * (lambda)^p = 1/8 for lambda=0.5.
  // Thus, set sigma <- sigma / sqrt(2) in the case.

  private def gromacsSoftcoreR (r : Double, sigma : Double) = {
    pow(pow(sigma,6.0)+pow(r,6.0), 1.0/6.0)
  }

  def coulombGromacsSoftcorePotential (sigma: Double) (r: Double) = {
    val ra = gromacsSoftcoreR(r, sigma)
    1.0/ra
  }

  def coulombGromacsSoftcoreForce (sigma: Double) (r: Double) = {
    val ra = gromacsSoftcoreR(r, sigma)
    1.0/(ra*ra) * pow(r/ra, 5.0)
  }


  ////////////////////////////////////////////////////////////////////////
  // Cutoff (Gromacs)

  ////////////////////////////////////////////////////////////////////////
  // van der Waals

  def ljPotential(r: Double) = {
    pow(r, -12.0) - pow(r, -6.0)
  }

  def ljForce (r: Double) = {
    12.0*pow(r, -13.0) - 6.0*pow(r, -7.0)
  }

}
