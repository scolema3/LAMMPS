/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(diffraction,FixDIFF)

#else

#ifndef LMP_FIX_DIFF_H
#define LMP_FIX_DIFF_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixDIFF : public Fix {
 public:
  FixDIFF(class LAMMPS *, int, char **);
  ~FixDIFF();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_vector(int);
  void reset_timestep(bigint);


 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int which;
  char *ids;
  
  int nrows;

  int ave,nwindow,nsum,startstep;
  int overwrite;
  long filepos;

  int norm,iwindow,window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;

  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, char **);

  bigint nextvalid();
  
  FILE    *fp_xyz;
  FILE    *fp_vtk;   
  char    *filename;         // user-specified file
  bool    steptime;           

  int     nOutput; 
  bool    vtkflag;
  bool    xyzflag;  

  double  Threshold;          // Intensity threshold value

  double  prd_inv[3];        // Inverse spacing of unit cell

  // XRD
  class ComputeXRD *compute_xrd;
  bool    xrdflag;
  double  lambda_xrd;        // Radiation wavelength (distance units)
  bool    manual_xrd;        // Turn on manual reciprocal map
  double  Kmax_xrd;          //
  double  dK_xrd[3];         // spacing of reciprocal points in each dimension  
  int     Knmax_xrd[3];      // maximum integer value for K points in each dimension
  int     Knmin_xrd[3];      // minimum integer value for K points in each dimension
  int     Dim_xrd[3];  
  double  Max2Theta;
  double  Min2Theta;
  double  ang;

  // SAED
  class ComputeSAED *compute_saed;
  bool    saedflag;
  double  lambda_saed;       // Radiation wavelength (distance units)
  bool    manual_saed;       // Turn on manual reciprocal map
  double  dK_saed[3];        // spacing of reciprocal points in each dimension
  double  Kmax_saed;         // Maximum reciprocal distance to exploren
  int     Knmax_saed[3];     // maximum integer value for K points in each dimension
  int     Knmin_saed[3];     // minimum integer value for K points in each dimension
  int     Dim_saed[3];       // integer dimension of vtk 
  double  Zone[3];           // Zone axis to view SAED
  double  R_Ewald;           // Radius of Ewald sphere (distance units)
  double  dR_Ewald;          // Thickness of Ewald sphere slice
  int     KnSlice[6];        // min 0-2 max 2-5 hkl index using zone
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix ave/time does not exist

Self-explanatory.

E: Fix ID for fix ave/time does not exist

Self-explanatory.

E: Invalid fix ave/time off column

Self-explantory.

E: Fix ave/time compute does not calculate a scalar

Self-explantory.

E: Fix ave/time compute does not calculate a vector

Self-explantory.

E: Fix ave/time compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ave/time compute does not calculate an array

Self-explanatory.

E: Fix ave/time compute array is accessed out-of-range

An index for the array is out of bounds.

E: Fix ave/time fix does not calculate a scalar

Self-explanatory.

E: Fix ave/time fix does not calculate a vector

Self-explanatory.

E: Fix ave/time fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/time not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/time
is requesting a value on a non-allowed timestep.

E: Fix ave/time fix does not calculate an array

Self-explanatory.

E: Fix ave/time fix array is accessed out-of-range

An index for the array is out of bounds.

E: Variable name for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time variable is not equal-style variable

Self-explanatory.

E: Fix ave/time cannot use variable with vector mode

Variables produce scalar values.

E: Fix ave/time columns are inconsistent lengths

Self-explanatory.

E: Fix ave/time cannot set output array intensive/extensive from these inputs

One of more of the vector inputs has individual elements which are
flagged as intensive or extensive.  Such an input cannot be flagged as
all intensive/extensive when turned into an array by fix ave/time.

E: Cannot open fix ave/time file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Fix ave/time missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
