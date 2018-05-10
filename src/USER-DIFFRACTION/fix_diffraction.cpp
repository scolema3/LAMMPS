/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Shawn Coleman (ARL)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fix_diffraction.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "compute_saed.h"
#include "compute_xrd.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{COMPUTE};
enum{ONE,RUNNING,WINDOW};
enum{FIRST,MULTI};

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

/* ---------------------------------------------------------------------- */

FixDIFF::FixDIFF(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix diffraction command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);

  global_freq = nfreq;

  nvalues = 0;
  int iarg = 6;

  // identify the diffraction computes to print
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0){
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues != 1) error->all(FLERR,"Illegal fix diffraction command");

  // parse optional flags
  options(narg,arg);

  nOutput = 0;
  ids = NULL;
  nvalues = 0;
  iarg = 6;

  double  c_xrd[3];          // Parameters controlling resolution of reciprocal space explored
  double  c_saed[3];         // Parameters controlling resolution of reciprocal space explored

  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ) {

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      char *ptr = strchr(suffix,'[');
      if (ptr) error->all(FLERR,"Illegal fix diffraction command");
      n = strlen(suffix) + 1;
      ids = new char[n];
      strcpy(ids,suffix);
      delete [] suffix;

      int icompute = modify->find_compute(ids);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix diffraction does not exist");

      Compute *compute = modify->compute[icompute];

      // SAED
      if (strcmp(compute->style,"saed") == 0) {
        saedflag = true;

        // gather varialbes from specified compute_saed
        compute_saed = (ComputeSAED*) modify->compute[icompute];
        double *saed_var = compute_saed->saed_var;
        lambda_saed   = saed_var[0];
        Kmax_saed     = saed_var[1];
        Zone[0]       = saed_var[2];
        Zone[1]       = saed_var[3];
        Zone[2]       = saed_var[4];
        c_saed[0]     = saed_var[5];
        c_saed[1]     = saed_var[6];
        c_saed[2]     = saed_var[7];
        dR_Ewald      = saed_var[8];
        double manual_double = saed_var[9];
        manual_saed = false;
        if (manual_double == 1) manual_saed = true;

        // Standard error check for compute saed
        if (compute->vector_flag == 0)
            error->all(FLERR,"Fix diffraction compute does not calculate a vector");
        if (compute->extvector != 0)
          error->all(FLERR,"Illegal fix diffraction command");
        nrows = compute->size_vector;
      }

      // XRD
      else if (strcmp(compute->style,"xrd") == 0) {
        xrdflag = true;

        // Gather varialbes from specified compute_xrd
        compute_xrd = (ComputeXRD*) modify->compute[icompute];
        double *xrd_var = compute_xrd->xrd_var;
        lambda_xrd     = xrd_var[0];
        Max2Theta      = xrd_var[1];
        Min2Theta      = xrd_var[2];
        c_xrd[0]       = xrd_var[3];
        c_xrd[1]       = xrd_var[4];
        c_xrd[2]       = xrd_var[5];
        double manual_double = xrd_var[6];
        manual_xrd = false;
        if (manual_double == 1) manual_xrd = true;

        // Standard error check for compute xrd
        if (compute->array_flag == 0)
          error->all(FLERR,"Fix diffraction compute does not calculate a array");
        if (2 != compute->size_array_cols)
          error->all(FLERR,"Fix diffraction compute array is accessed out-of-range");
        if (compute->extarray != 0)
          error->all(FLERR,"Illegal fix diffraction command");
        nrows = compute->size_array_rows;

      } else error->all(FLERR,"Illegal fix diffraction command");

      nvalues++;
      iarg++;
    } else break;
  }

  // Setup Histogram
  if (histflag)
    bin = new double[nbins];


  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix diffraction command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix diffraction command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix diffraction command");


  // allocate memory for averaging
  vector = vector_total = NULL;
  vector_list = NULL;

  if (ave == WINDOW)
    memory->create(vector_list,nwindow,nvalues,"diffraction:vector_list");
  memory->create(vector,nrows,"diffraction:vector");
  memory->create(vector_total,nrows,"diffraction:vector_total");
  extlist = NULL;

  vector_flag = 1;
  size_vector = nrows;


  // Common parameters for xrd and saed computes
  int *periodicity = domain->periodicity;
  double *prd;
  double ave_inv = 0.0;
  prd = domain->prd;

  if (periodicity[0]){
    prd_inv[0] = 1 / prd[0];
    ave_inv += prd_inv[0];
  }
  if (periodicity[1]){
    prd_inv[1] = 1 / prd[1];
    ave_inv += prd_inv[1];
  }
  if (periodicity[2]){
    prd_inv[2] = 1 / prd[2];
    ave_inv += prd_inv[2];
  }

  ave_inv = ave_inv / (periodicity[0] + periodicity[1] + periodicity[2]);
  if (!periodicity[0]) prd_inv[0] = ave_inv;
  if (!periodicity[1]) prd_inv[1] = ave_inv;
  if (!periodicity[2]) prd_inv[2] = ave_inv;



  // SAED specific maps
  if (saedflag == true) {

    // Zone flag to capture entire recrocal space volume
    if (  (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    } else {
        R_Ewald = (1 / lambda_saed);
        double Rnorm = R_Ewald/ sqrt(Zone[0] * Zone[0] +
                       Zone[1] * Zone[1] +  Zone[2]* Zone[2]);
        Zone[0] = Zone[0] * Rnorm;
        Zone[1] = Zone[1] * Rnorm;
        Zone[2] = Zone[2] * Rnorm;
    }

    // Use manual_saed mapping of reciprocal lattice
    if (manual_saed) {
      for (int i=0; i<3; i++) prd_inv[i] = 1.0;
    }

    // Find integer dimensions of the reciprocal lattice box bounds
    if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
      for (int i=0; i<3; i++) {
        dK_saed[i] = prd_inv[i]*c_saed[i];
        Knmax_saed[i] = ceil(Kmax_saed / dK_saed[i]);
        Knmin_saed[i] = -Knmax_saed[i];
      }
    } else {

      for (int i=0; i<3; i++) {
        Knmax_saed[i] = -10000;
        Knmin_saed[i] =  10000;
      }
      double dinv2 = 0.0;
      double r = 0.0;
      double K[3];
      int Ksearch[3];
      for (int i=0; i<3; i++) {
        dK_saed[i] = prd_inv[i]*c_saed[i];
        Ksearch[i] = ceil(Kmax_saed / dK_saed[i]);
      }

      for (int k = -Ksearch[2]; k <= Ksearch[2]; k++) {
        for (int j = -Ksearch[1]; j <= Ksearch[1]; j++) {
          for (int i = -Ksearch[0]; i <= Ksearch[0]; i++) {
            K[0] = i * dK_saed[0];
            K[1] = j * dK_saed[1];
            K[2] = k * dK_saed[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (dinv2 < Kmax_saed * Kmax_saed) {
              r = 0.0;
              for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
              r = sqrt(r);
              if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
                if ( i < Knmin_saed[0] ) Knmin_saed[0] = i;
                if ( j < Knmin_saed[1] ) Knmin_saed[1] = j;
                if ( k < Knmin_saed[2] ) Knmin_saed[2] = k;
                if ( i > Knmax_saed[0] ) Knmax_saed[0] = i;
                if ( j > Knmax_saed[1] ) Knmax_saed[1] = j;
                if ( k > Knmax_saed[2] ) Knmax_saed[2] = k;
              }
            }
          }
        }
      }
    }

    if ( vtkflag == true ) {
      // Finding dimensions for vtk files
      for (int i=0; i<3; i++) {
        if ( ( (Knmin_saed[i] > 0) && (Knmax_saed[i] > 0) ) || ( (Knmin_saed[i] < 0) && (Knmax_saed[i] < 0) ) ){
          Dim_saed[i] = abs( (int) Knmin_saed[i] ) + abs( (int) Knmax_saed[i] );
        } else Dim_saed[i] = abs( (int) Knmin_saed[i] ) + abs( (int) Knmax_saed[i] ) + 1;
      }
    }
  }

  // XRD specific maps
  if (xrdflag == true) {

    Kmax_xrd = 2 * sin(Max2Theta) / lambda_xrd;

    // Use manual mapping of reciprocal lattice
    if (manual_xrd) {
      for (int i=0; i<3; i++)
        prd_inv[i] = 1.0;
    }

    // Find reciprocal spacing and integer dimensions
    for (int i=0; i<3; i++) {
      dK_xrd[i] = prd_inv[i]*c_xrd[i];
      Knmax_xrd[i] = ceil(Kmax_xrd / dK_xrd[i]);
      Knmin_xrd[i] = -Knmax_xrd[i];
    }


   // Finding dimensions for vtk files
    for (int i=0; i<3; i++) {
      if ( ( (Knmin_xrd[i] > 0) && (Knmax_xrd[i] > 0) ) || ( (Knmin_xrd[i] < 0) && (Knmax_xrd[i] < 0) ) ){
        Dim_xrd[i] = abs( (int) Knmin_xrd[i] ) + abs( (int) Knmax_xrd[i] );
      } else Dim_xrd[i] = abs( (int) Knmin_xrd[i] ) + abs( (int) Knmax_xrd[i] ) + 1;
    }
  }

  // initialization

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  for (int i = 0; i < nrows; i++)
     vector_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

}

/* ---------------------------------------------------------------------- */

FixDIFF::~FixDIFF()
{
  delete [] extlist;
  memory->destroy(vector);
  memory->destroy(vector_list);
  memory->destroy(vector_total);
  if (fp_xyz && me == 0) fclose(fp_xyz);
  if (fp_vtk && me == 0) fclose(fp_vtk);
  fp_xyz = NULL;
  fp_vtk = NULL;
}

/* ---------------------------------------------------------------------- */

int FixDIFF::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDIFF::init()
{
  // set current indices for all computes,fixes,variables


  int icompute = modify->find_compute(ids);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix diffraction does not exist");

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixDIFF::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixDIFF::end_of_step()
{
  // skip if not step which requires doing something
  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixDIFF::invoke_vector(bigint ntimestep)
{
  // zero if first step
  int icompute = modify->find_compute(ids);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix diffraction does not exist");

  if (irepeat == 0) {
    for (int i = 0; i < nrows; i++)
       vector[i] = 0.0;
    if (histflag) 
      for (int i = 0; i < nbins; i++) bin[i] = 0.0;
  }


  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add
  modify->clearstep_compute();

  // invoke compute if not previously invoked
  Compute *compute = modify->compute[icompute];

  if (saedflag == true) {

    if (!(compute->invoked_flag & INVOKED_VECTOR)) {
      compute->compute_vector();
      compute->invoked_flag |= INVOKED_VECTOR;
    }

    double *cvector = compute->vector;
    for (int i = 0; i < nrows; i++)
      vector[i] = cvector[i];
  }

  if (xrdflag == true) {
    if (!(compute->invoked_flag & INVOKED_ARRAY)) {
      compute->compute_array();
      compute->invoked_flag |= INVOKED_ARRAY;
    }
    double **carray = compute->array;
    int icol = 1;
    for (int i = 0; i < nrows; i++)
      vector[i] = carray[i][icol];
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for ( int i = 0; i < nrows; i++)
    vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (int i = 0; i < nrows; i++) vector_total[i] = vector[i];
    norm = 1;
  } else if (ave == RUNNING) {
    for (int i = 0; i < nrows; i++) vector_total[i] += vector[i];
    norm++;
  } else if (ave == WINDOW) {
    for (int i = 0; i < nrows; i++) {
      vector_total[i] += vector[i];
      if (window_limit) vector_total[i] -= vector_list[iwindow][i];
        vector_list[iwindow][i] = vector[i];
    }
    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }


  // Writing to file
  // Index files using timestep or by output count
  int indexP;
  if (steptime == true) indexP = ntimestep;
  else indexP = nOutput;

  // XYZ format initializations
  if (xyzflag && me == 0) {
    // Iterate filenames for multiple snapshots
    if ( nOutput > 0 ) fclose(fp_xyz);
    fp_xyz = NULL;
    char nName [128];
    sprintf(nName,"%s.%d.xyz",filename,indexP);
    fp_xyz = fopen(nName,"w");
    if (fp_xyz == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix diffraction-xyz file %s",nName);
      error->one(FLERR,str);
    }
  }

  // VTK format initializations
  if (vtkflag && me == 0) {
    // Iterate filenames for multiple snapshots
    if ( nOutput > 0 ) fclose(fp_vtk);
    fp_vtk = NULL;
    char nName [128];
    sprintf(nName,"%s.%d.vtk",filename,indexP);
    fp_vtk = fopen(nName,"w");
    if (fp_vtk == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix diffraction-vtk file %s",nName);
      error->one(FLERR,str);
    }
  }


  // HISTOGRAM
  if (histflag && me == 0) {
    // Iterate filenames for multiple snapshots
    if ( nOutput > 0 ) fclose(fp_hist);
    fp_hist = NULL;
    char nName [128];
    sprintf(nName,"%s.%d.hist",filename,indexP);
    fp_hist = fopen(nName,"w");
    if (fp_hist == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix diffraction-hist file %s",nName);
      error->one(FLERR,str);
    }
  }


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


  // HIST format for SAED
  if (fp_hist && saedflag && me == 0) {

    double coord = 0.0;
    double value = 0.0;
    double lo = 0.0;
    double hi = Kmax_saed;
    double binsize = (hi-lo)/nbins;
    double bininv = 1.0/binsize;


    // Finding the intersection of the reciprocal space and Ewald sphere
    int NROW1 = 0;
    double dinv2 = 0.0;
    double r = 0.0;
    double K[3];

    // Zone flag to capture entire reciprocal space volume
    if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
      for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
        for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
          for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
            K[0] = i * dK_saed[0];
            K[1] = j * dK_saed[1];
            K[2] = k * dK_saed[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (dinv2 < Kmax_saed * Kmax_saed) {
              value = sqrt(dinv2);
              int ibin = static_cast<int> ((value-lo)*bininv);
              ibin = MIN(ibin,nbins-1);
              bin[ibin] += vector_total[NROW1]/norm;
              NROW1++;
            }
          }
        }
      }
    } else {
      for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
        for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
          for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
            K[0] = i * dK_saed[0];
            K[1] = j * dK_saed[1];
            K[2] = k * dK_saed[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (dinv2 < hi) {
              r = 0.0;
              for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
              r = sqrt(r);
              if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
                value = sqrt(dinv2);
                int ibin = static_cast<int> ((value-lo)*bininv);
                ibin = MIN(ibin,nbins-1);
                bin[ibin] += vector_total[NROW1]/norm;
                NROW1++;
              }
            }
          }
        }
      }
    }

    // Print header information
    fprintf(fp_hist,"#LAMBDA %g\n", lambda_saed);
    fprintf(fp_hist,"#ASPECT_RATIO %g %g %g\n", dK_saed[0], dK_saed[1], dK_saed[2]);
    fprintf(fp_hist,"#NBINS %g\n", nbins);

    filepos = ftell(fp_hist);
    if (overwrite) fseek(fp_hist,filepos,SEEK_SET);

    for (int i = 0; i < nbins; i++) {
        coord = lo + (i+0.5)*binsize;
        fprintf(fp_hist,"%d %g %g\n",i+1,coord,bin[i]);
    }

 } //END HIST SAED


  // HIST format for XRD
  if (fp_hist && saedflag && me == 0) {

    double coord = 0.0;
    double value = 0.0;
    double lo = Min2Theta;
    double hi = Max2Theta;
    double binsize = (hi-lo)/nbins;
    double bininv = 1.0/binsize;


    // Finding the intersection of the reciprocal space and Ewald sphere
    int NROW1 = 0;
    double dinv2 = 0.0;
    double r = 0.0;
    double K[3];

    for (int i = -Knmax_xrd[0]; i <= Knmax_xrd[0]; i++) {
      for (int j = -Knmax_xrd[1]; j <= Knmax_xrd[1]; j++) {
        for (int k = -Knmax_xrd[2]; k <= Knmax_xrd[2]; k++) {
          K[0] = i * dK_xrd[0];
          K[1] = j * dK_xrd[1];
          K[2] = k * dK_xrd[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if  (4 >= dinv2 * lambda_xrd * lambda_xrd ) {
	    ang = asin(lambda_xrd * sqrt(dinv2) / 2);
            if ( (ang <= Max2Theta) & (ang >= Min2Theta) ) {
                value = ang;
                int ibin = static_cast<int> ((value-lo)*bininv);
                ibin = MIN(ibin,nbins-1);
                bin[ibin] += vector_total[NROW1]/norm;
                NROW1++;
            }
          }
        }
      }
    }



    // Print header information
    fprintf(fp_hist,"#LAMBDA %g\n", lambda_saed);
    fprintf(fp_hist,"#ASPECT_RATIO %g %g %g\n", dK_saed[0], dK_saed[1], dK_saed[2]);
    fprintf(fp_hist,"#NBINS %g\n", nbins);

    filepos = ftell(fp_hist);
    if (overwrite) fseek(fp_hist,filepos,SEEK_SET);


    for (int i = 0; i < nbins; i++) {
        coord = lo + (i+0.5)*binsize;
        fprintf(fp_hist,"%d %g %g\n",i+1,coord,bin[i]);
    }

 } //END HIST XRD


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


 
  // XYZ format for SAED
  if (fp_xyz && saedflag && me == 0) {

    // Print header information
    fprintf(fp_xyz,"#LAMBDA %g\n", lambda_saed);
    fprintf(fp_xyz,"#ASPECT_RATIO %g %g %g\n", dK_saed[0], dK_saed[1], dK_saed[2]);

    filepos = ftell(fp_xyz);
    if (overwrite) fseek(fp_xyz,filepos,SEEK_SET);

    // Finding the intersection of the reciprocal space and Ewald sphere
    int NROW1 = 0;
    double dinv2 = 0.0;
    double r = 0.0;
    double K[3];

    if ( (Threshold > 0)) {
      fprintf(fp_xyz,"#THRESHOLD %g\n", Threshold);
      // Zone flag to capture entire reciprocal space volume
      if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
        for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
          for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
            for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
              K[0] = i * dK_saed[0];
              K[1] = j * dK_saed[1];
              K[2] = k * dK_saed[2];
              dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
              if (dinv2 < Kmax_saed * Kmax_saed) {
                if ( (vector_total[NROW1]/norm >= Threshold ))
                  fprintf(fp_xyz,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1]/norm);
                NROW1++;
              }
            }
          }
        }
      } else {
        for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
          for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
            for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
              K[0] = i * dK_saed[0];
              K[1] = j * dK_saed[1];
              K[2] = k * dK_saed[2];
              dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
              if (dinv2 < Kmax_saed * Kmax_saed) {
                r = 0.0;
                for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
                r = sqrt(r);
                if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
                  NROW1++;
                  if ( (vector_total[NROW1]/norm >= Threshold ))
                    fprintf(fp_xyz,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1]/norm);
                }
              }
            }
          }
        }
      }
    } else { // No Threshold
      // Zone flag to capture entire reciprocal space volume
      if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
        for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
          for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
            for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
              K[0] = i * dK_saed[0];
              K[1] = j * dK_saed[1];
              K[2] = k * dK_saed[2];
              dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
              if (dinv2 < Kmax_saed * Kmax_saed) {
                fprintf(fp_xyz,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1]/norm);
                NROW1++;
              }
            }
          }
        }
      } else {
        for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
          for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
            for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
              K[0] = i * dK_saed[0];
              K[1] = j * dK_saed[1];
              K[2] = k * dK_saed[2];
              dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
              if (dinv2 < Kmax_saed * Kmax_saed) {
                r = 0.0;
                for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
                r = sqrt(r);
                if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
                  fprintf(fp_xyz,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1]/norm);
                  NROW1++;
                }
              }
            }
          }
        }
      }
    }
  } // END - XYZ format for SAED


  // XYZ format for XRD
  if (fp_xyz && xrdflag && me == 0) {

    // Header information with relp spacing
    fprintf(fp_xyz,"#LAMBDA %g\n", lambda_xrd);
    fprintf(fp_xyz,"#ASPECT_RATIO %g %g %g\n", dK_xrd[0], dK_xrd[1], dK_xrd[2]);

    filepos = ftell(fp_xyz);
    if (overwrite) fseek(fp_xyz,filepos,SEEK_SET);

    int NROW1 = 0;
    double dinv2 = 0.0;
    double K[3];
    if ( (Threshold > 0)) {
      fprintf(fp_xyz,"#THRESHOLD %g\n", Threshold);
      for (int i = -Knmax_xrd[0]; i <= Knmax_xrd[0]; i++) {
        for (int j = -Knmax_xrd[1]; j <= Knmax_xrd[1]; j++) {
          for (int k = -Knmax_xrd[2]; k <= Knmax_xrd[2]; k++) {
            K[0] = i * dK_xrd[0];
            K[1] = j * dK_xrd[1];
            K[2] = k * dK_xrd[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (4 >= dinv2 * lambda_xrd * lambda_xrd ) {
       	      ang = asin(lambda_xrd * sqrt(dinv2) / 2);
              if ( (ang <= Max2Theta) & (ang >= Min2Theta) ) {
                if ( (vector_total[NROW1]/norm >= Threshold ))
                  fprintf(fp_xyz,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1]/norm);
                NROW1++;
	          }
            }
          }
        }
      }
    }
    else { // No Threshold
      for (int i = -Knmax_xrd[0]; i <= Knmax_xrd[0]; i++) {
        for (int j = -Knmax_xrd[1]; j <= Knmax_xrd[1]; j++) {
          for (int k = -Knmax_xrd[2]; k <= Knmax_xrd[2]; k++) {
            K[0] = i * dK_xrd[0];
            K[1] = j * dK_xrd[1];
            K[2] = k * dK_xrd[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if  (4 >= dinv2 * lambda_xrd * lambda_xrd ) {
       	      ang = asin(lambda_xrd * sqrt(dinv2) / 2);
              if ( (ang <= Max2Theta) & (ang >= Min2Theta) ) {
                fprintf(fp_xyz,"%g %g %g %g\n",K[0],K[1],K[2],vector_total[NROW1]/norm);
                NROW1++;
              }
            }
          }
        }
      }
    }
  } // END - XYZ format for XRD



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


  // VTK format for SAED
  if (fp_vtk && saedflag && me == 0) {
    fprintf(fp_vtk,"# vtk DataFile Version 3.0 c_%s\n",ids);
    fprintf(fp_vtk,"Image data set\n");
    fprintf(fp_vtk,"ASCII\n");
    fprintf(fp_vtk,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp_vtk,"DIMENSIONS %d %d %d\n", Dim_saed[0],  Dim_saed[1], Dim_saed[2]);
    fprintf(fp_vtk,"ASPECT_RATIO %g %g %g\n", dK_saed[0], dK_saed[1], dK_saed[2]);
    fprintf(fp_vtk,"ORIGIN %g %g %g\n", Knmin_saed[0] * dK_saed[0],  Knmin_saed[1] * dK_saed[1], Knmin_saed[2] * dK_saed[2]);
    fprintf(fp_vtk,"POINT_DATA %d\n",  Dim_saed[0] *  Dim_saed[1] * Dim_saed[2] );
    fprintf(fp_vtk,"SCALARS intensity float\n");
    fprintf(fp_vtk,"LOOKUP_TABLE default\n");

    filepos = ftell(fp_vtk);
    if (overwrite) fseek(fp_vtk,filepos,SEEK_SET);

   // Finding the intersection of the reciprocal space and Ewald sphere
    int NROW1 = 0;
    double dinv2 = 0.0;
    double r = 0.0;
    double K[3];

    // Zone flag to capture entire reciprocal space volume
    if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
      for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
        for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
          for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
            K[0] = i * dK_saed[0];
            K[1] = j * dK_saed[1];
            K[2] = k * dK_saed[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (dinv2 < Kmax_saed * Kmax_saed) {
               fprintf(fp_vtk,"%g\n",vector_total[NROW1]/norm);
               NROW1++;
            } else
              fprintf(fp_vtk,"%d\n",-1);
          }
        }
      }
    } else {
      for (int k = Knmin_saed[2]; k <= Knmax_saed[2]; k++) {
        for (int j = Knmin_saed[1]; j <= Knmax_saed[1]; j++) {
          for (int i = Knmin_saed[0]; i <= Knmax_saed[0]; i++) {
            K[0] = i * dK_saed[0];
            K[1] = j * dK_saed[1];
            K[2] = k * dK_saed[2];
            dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
            if (dinv2 < Kmax_saed * Kmax_saed) {
              r=0.0;
              for (int m=0; m<3; m++) r += pow(K[m] - Zone[m],2.0);
              r = sqrt(r);
              if  ( (r >  (R_Ewald - dR_Ewald) ) && (r < (R_Ewald + dR_Ewald) ) ){
               fprintf(fp_vtk,"%g\n",vector_total[NROW1]/norm);
               NROW1++;
              } else
                fprintf(fp_vtk,"%d\n",-1);
            } else
              fprintf(fp_vtk,"%d\n",-1);
          }
        }
      }
    }
  } // END - VTK format for SAED

  // VTK format for XRD
  if (fp_vtk && xrdflag && me == 0) {
    fprintf(fp_vtk,"# vtk DataFile Version 3.0 c_%s\n",ids);
    fprintf(fp_vtk,"Image data set\n");
    fprintf(fp_vtk,"ASCII\n");
    fprintf(fp_vtk,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp_vtk,"DIMENSIONS %d %d %d\n", Dim_xrd[0],  Dim_xrd[1], Dim_xrd[2]);
    fprintf(fp_vtk,"ASPECT_RATIO %g %g %g\n", dK_xrd[0], dK_xrd[1], dK_xrd[2]);
    fprintf(fp_vtk,"ORIGIN %g %g %g\n", Knmin_xrd[0] * dK_xrd[0],  Knmin_xrd[1] * dK_xrd[1], Knmin_xrd[2] * dK_xrd[2]);
    fprintf(fp_vtk,"POINT_DATA %d\n",  Dim_xrd[0] *  Dim_xrd[1] * Dim_xrd[2] );
    fprintf(fp_vtk,"SCALARS intensity float\n");
    fprintf(fp_vtk,"LOOKUP_TABLE default\n");

    filepos = ftell(fp_vtk);

    if (overwrite) fseek(fp_vtk,filepos,SEEK_SET);

   // Finding the intersection of the reciprocal space and Ewald sphere
    int NROW1 = 0;
    double dinv2 = 0.0;
    double K[3];
    for (int i = -Knmax_xrd[0]; i <= Knmax_xrd[0]; i++) {
      for (int j = -Knmax_xrd[1]; j <= Knmax_xrd[1]; j++) {
        for (int k = -Knmax_xrd[2]; k <= Knmax_xrd[2]; k++) {
          K[0] = i * dK_xrd[0];
          K[1] = j * dK_xrd[1];
          K[2] = k * dK_xrd[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if  (4 >= dinv2 * lambda_xrd * lambda_xrd ) {
         	  ang = asin(lambda_xrd * sqrt(dinv2) / 2);
            if ( (ang <= Max2Theta) & (ang >= Min2Theta) ) {
              NROW1++;
              fprintf(fp_vtk,"%g\n",vector_total[NROW1-1]/norm);
  	        } else
              fprintf(fp_vtk,"%d\n",-1);
          } else
            fprintf(fp_vtk,"%d\n",-1);
        }
      }
    }
  } // END - VTK format for XRD

  nOutput++;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixDIFF::compute_vector(int i)
{
  if (norm) {
    return vector_total[i]/norm;
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixDIFF::options(int narg, char **arg)
{
  // option defaults

  fp_vtk = NULL;
  fp_xyz = NULL;
  fp_hist = NULL;

  ave = ONE;
  startstep = 0;
  overwrite = 0;
  Threshold = -1;
  nbins = 100; 

  vtkflag = false;
  xyzflag = false;
  histflag = false;

  steptime = false;

  bool formatflag = false;

  // optional args
  int iarg = 6 + nvalues;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix diffraction command");
      int n = strlen(arg[iarg+1]) + 1;
      filename = new char[n];
      strcpy(filename,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      formatflag = true;
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix diffraction command");
      if (strcmp(arg[iarg+1],"vtk") == 0) vtkflag=true;
      else if (strcmp(arg[iarg+1],"xyz") == 0) xyzflag=true;
      else if (strcmp(arg[iarg+1],"hist") == 0) histflag=true;
      else  error->all(FLERR,"Illegal fix diffraction command");
      iarg += 2;
   } else if (strcmp(arg[iarg],"nbins") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix diffraction command");
      nbins = atof(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix diffraction command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix diffraction command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix diffraction command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix diffraction command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix diffraction command");
      startstep = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"timestep") == 0) {
      steptime = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"threshold") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix diffraction command");
      Threshold = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix diffraction command");
  }

  if ( formatflag == false)
    error->all(FLERR,"Illegal fix diffraction command - define format");

}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixDIFF::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixDIFF::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix diffraction missed timestep");
}
