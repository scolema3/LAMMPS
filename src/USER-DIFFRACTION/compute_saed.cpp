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
   Contributing authors: Shawn Coleman & Douglas Spearot (Arkansas)
------------------------------------------------------------------------- */

// Attempting compatibility with USER-INTEL
#ifdef LMP_INTEL_OFFLOAD
#define _LMP_INTEL_OFFLOAD
#include "offload.h"
#endif

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "math_const.h"
#include "compute_saed.h"
#include "compute_saed_consts.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "stdio.h"
#include "string.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_compute_saed_c[] =
  "compute_saed command:\n\n"
  "@Article{Coleman13,\n"
  " author = {S. P. Coleman, D. E. Spearot, L. Capolungo},\n"
  " title = {Virtual diffraction analysis of Ni [010] symmetric tilt grain boundaries},\n"
  " journal = {Modelling and Simulation in Materials Science and Engineering},\n"
  " year =    2013,\n"
  " volume =  21,\n"
  " pages =   {055020}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

ComputeSAED::ComputeSAED(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_saed_c);

  int ntypes = atom->ntypes;
  int natoms = group->count(igroup);
  int dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int triclinic = domain->triclinic;
  me = comm->me;
  
  // Checking errors specific to the compute
  if (dimension == 2) 
    error->all(FLERR,"Compute SAED does not work with 2d structures");
  if (narg < 4+ntypes) 
    error->all(FLERR,"Illegal Compute SAED Command");

  vector_flag = 1;
  extvector = 0;
 
  // Store radiation wavelength
  lambda = atof(arg[3]);
  if (lambda < 0)
    error->all(FLERR,"Compute SAED: Wavelength must be greater than zero");

  // Define atom types for atomic scattering factor coefficients
  int iarg = 4;
  ztype = new int[ntypes];
  for (int i = 0; i < ntypes; i++){
    ztype[i] = SAEDmaxType + 1;
  }
  for (int i=0; i<ntypes; i++) {
    for(int j = 0; j < SAEDmaxType; j++){
      if (strcasecmp(arg[iarg],SAEDtypeList[j]) == 0) {
        ztype[i] = j;
       }
     }
    if ( ztype[i] == SAEDmaxType + 1 )
        error->all(FLERR,"Compute SAED: Invalid ASF atom type");
    iarg++;
  }

  // Set defaults for optional args
  Kmax = 1.70;
  Zone[0] = 1; Zone[1] = 0; Zone[2] = 0;
  c[0] = 1; c[1] = 1; c[2] = 1;
  dR_Ewald = 0.01 / 2;
  manual = false;
  double manual_double = 0; 
  complex = false;
  double complex_double = 0;
  echo = false;
  ratio = 0.5;

  // Process optional args
  while (iarg < narg) {
    
    if (strcmp(arg[iarg],"Kmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      Kmax = atof(arg[iarg+1]);
      if (Kmax / 2 < 0 || Kmax / 2 > 6) 
        error->all(FLERR,"Compute SAED: |K|max/2 must be between 0 and 6 ");
      iarg += 2;
        
    } else if (strcmp(arg[iarg],"Zone") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      Zone[0] = atof(arg[iarg+1]);
      Zone[1] = atof(arg[iarg+2]);
      Zone[2] = atof(arg[iarg+3]);
      iarg += 4;
      
    } else if (strcmp(arg[iarg],"c") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      c[0] = atof(arg[iarg+1]);
      c[1] = atof(arg[iarg+2]);
      c[2] = atof(arg[iarg+3]);
      if (c[0] < 0 || c[1] < 0 || c[2] < 0) 
        error->all(FLERR,"Compute SAED: dKs must be greater than 0");  
      iarg += 4;
      
    } else if (strcmp(arg[iarg],"dR_Ewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute SAED Command");
      dR_Ewald = atof(arg[iarg+1]);
      if (dR_Ewald < 0) 
        error->all(FLERR,"Compute SAED: dR_Ewald slice must be greater than 0"); 
      iarg += 2;

    } else if (strcmp(arg[iarg],"ratio") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal Compute XRD Command");
      ratio = atof(arg[iarg+1]);
      if (ratio < 0) 
        error->all(FLERR,"Ratio value must be greater than or equal to zero");
      if (ratio > 1) 
        error->all(FLERR,"Ratio value must be less than or equal to one");
      iarg += 2;

    } else if (strcmp(arg[iarg],"echo") == 0) {
      echo = true;
      iarg += 1;

    } else if (strcmp(arg[iarg],"manual") == 0) {
      manual = true;
      manual_double = 1; 
      iarg += 1;        
                
    } else error->all(FLERR,"Illegal Compute SAED Command");
  }

  // Zone flag to capture entire reciprocal space volume
  if (  (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
  } else {
      R_Ewald = (1 / lambda); 
      double Rnorm = R_Ewald/ sqrt(Zone[0] * Zone[0] + 
                     Zone[1] * Zone[1] +  Zone[2]* Zone[2]);
      Zone[0] = Zone[0] * Rnorm; 
      Zone[1] = Zone[1] * Rnorm; 
      Zone[2] = Zone[2] * Rnorm; 
  }

  // Procedure to determine how many rows are needed given the constraints on 2theta
  // Calculating spacing between reciprocal lattice points 
  // Using distance based on periodic repeating distance
  if (!manual) {  
    if (!periodicity[0] && !periodicity[1] && !periodicity[2])
      error->all(FLERR,"Compute SAED must have at least one periodic boundary unless manual spacing specified");
    if (triclinic == 1) 
      error->all(FLERR,"Compute SAED does not work with triclinic structures -- turn on manual spacing"); 

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

    // Using the average inverse dimensions for non-periodic direction
    ave_inv = ave_inv / (periodicity[0] + periodicity[1] + periodicity[2]);
    if (!periodicity[0]){
      prd_inv[0] = ave_inv; 
    } 
    if (!periodicity[1]){
      prd_inv[1] = ave_inv;
    } 
    if (!periodicity[2]){
      prd_inv[2] = ave_inv;
    }
  }

  // Use manual mapping of reciprocal lattice 
  if (manual) {
    for (int i=0; i<3; i++) {
      prd_inv[i] = 1.0;
    }
  } 
  
  // Find reciprocal spacing and integer dimensions
  for (int i=0; i<3; i++) {
    dK[i] = prd_inv[i]*c[i];
    Knmax[i] = ceil(Kmax / dK[i]);
  }
 
  // Finding the intersection of the reciprocal space and Ewald sphere
  int n = 0;
  double dinv2, r2, EmdR2, EpdR2;
  double K[3];
  
  // Zone flag to capture entire reciprocal space volume
  if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) n++;
        } 
      }
    } 
  } else {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            r2 = 0.0;
            for (int m=0; m<3; m++)
              r2 += pow(K[m] - Zone[m],2.0);
            EmdR2 = pow(R_Ewald - dR_Ewald,2);
            EpdR2 = pow(R_Ewald + dR_Ewald,2);          
            if  ( (r2 >  EmdR2 ) && (r2 < EpdR2 ) ) {
              n++;
            }
          }
        } 
      }
    } 
  }

  if (me == 0) {
    if (screen &&  echo)
      fprintf(screen,"-----\nCompute SAED id:%s, # of atoms:%d, # of relp:%d\n",id,natoms,n);
      fprintf(screen,"Reciprocal point spacing in k1,k2,k3 = %g %g %g\n-----\n",
              dK[0], dK[1], dK[2]);
  }  

  nRows = n;
  size_vector = n;
  memory->create(vector,size_vector,"saed:vector");
  memory->create(store_tmp,3*size_vector,"saed:store_tmp");

  // Create vector of variables to be passed to fix saed/vtk
  saed_var[0] = lambda;
  saed_var[1] = Kmax;
  saed_var[2] = Zone[0];
  saed_var[3] = Zone[1];
  saed_var[4] = Zone[2];
  saed_var[5] = c[0];
  saed_var[6] = c[1];
  saed_var[7] = c[2];
  saed_var[8] = dR_Ewald;
  saed_var[9] = manual_double;
}

/* ---------------------------------------------------------------------- */

ComputeSAED::~ComputeSAED()
{
  memory->destroy(vector);
  memory->destroy(store_tmp);
  delete ztype;
}

/* ---------------------------------------------------------------------- */

void ComputeSAED::init()
{
  double dinv2, r2, EmdR2, EpdR2;
  double K[3];
  int n = 0;
  
  // Zone 0 0 0 flag to capture entire reciprocal space volume
  if ( (Zone[0] == 0) && (Zone[1] == 0) && (Zone[2] == 0) ){
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            store_tmp[3*n] = i;
            store_tmp[3*n+1] = j;
            store_tmp[3*n+2] = k;
            n++;
          }
        } 
      } 
    }
  } else {
    for (int k = -Knmax[2]; k <= Knmax[2]; k++) {
      for (int j = -Knmax[1]; j <= Knmax[1]; j++) {
        for (int i = -Knmax[0]; i <= Knmax[0]; i++) {
          K[0] = i * dK[0];
          K[1] = j * dK[1];
          K[2] = k * dK[2];
          dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
          if (dinv2 < Kmax * Kmax) {
            r2=0.0;
            for (int m=0; m<3; m++)
              r2 += pow(K[m] - Zone[m],2.0);
            EmdR2 = pow(R_Ewald - dR_Ewald,2);
            EpdR2 = pow(R_Ewald + dR_Ewald,2);          
            if  ( (r2 >  EmdR2 ) && (r2 < EpdR2 ) ) {
              store_tmp[3*n] = i;
              store_tmp[3*n+1] = j;
              store_tmp[3*n+2] = k;
              n++;
            }
          }
        }
      } 
    } 
  }
  if (n != nRows)  error->all(FLERR,"Compute SAED Nrows inconsistent");
}

/* ---------------------------------------------------------------------- */

void ComputeSAED::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen,"-----\nComputing SAED intensities for %s\n", id);
  }
  
  
  double t0 = MPI_Wtime();
  double *Fvec = new double[2*nRows]; // Strct factor  (real & imaginary)
 
  ntypes = atom->ntypes;
  int nlocal = atom->nlocal;
  int *type  = atom->type;
  int natoms = group->count(igroup);
  int *mask = atom->mask;

// ==========================================================
// Begin Setups (5 total)
// ==========================================================

  // 1- Setup for vectorizing compute (removing mask/group dependency on xlocal)
  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     nlocalgroup++;
    }
  }

  double *xlocal = new double [3*nlocalgroup];
  int *typelocal = new int [nlocalgroup];

  nlocalgroup = 0;
  for (int ii = 0; ii < nlocal; ii++) {
    if (mask[ii] & groupbit) {
     xlocal[3*nlocalgroup] = atom->x[ii][0];
     xlocal[3*nlocalgroup+1] = atom->x[ii][1];
     xlocal[3*nlocalgroup+2] = atom->x[ii][2];
     typelocal[nlocalgroup]=type[ii];
     nlocalgroup++;
    }
  } 
 
  // 2- Setup for blocking compute to remain on L1 & L2 Cache
  int wblockCPU = 8192;
  int nblockCPU = floor(nlocalgroup/wblockCPU);
  if ( nblockCPU == 0 ) wblockCPU=nlocalgroup;

#ifdef _LMP_INTEL_OFFLOAD   
  int wblockMIC = 4096;
  int nblockMIC = floor(nlocalgroup/wblockMIC);
  if ( nblockMIC == 0 ) wblockMIC=nlocalgroup;
#endif  
  int iistart = 0;
  int iistop = 0;

  if (me == 0 && echo) {
      if (screen)
        fprintf(screen," - %d block sections on CPU (Block size = %d) \n", nblockCPU, wblockCPU);
#ifdef _LMP_INTEL_OFFLOAD           
        fprintf(screen," - %d block sections on MIC (Block size = %d) \n", nblockMIC, wblockMIC);
#endif
  }  
  
  // 3- Setup for using OpenMP   
  int *ztype_omp = ztype;
#if defined(_OPENMP)
  int *store_omp = store_tmp;
  int nthreadCPU = comm->nthreads;
  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," - each MPI is using %d OMP threads on CPU",nthreadCPU);
  }
#endif  

  // 4- Setup for using MIC offloading    
  int nRowsMIC = 0;
  int nRowsCPU = nRows;

#ifdef _LMP_INTEL_OFFLOAD
  int nthreadMIC = 1;
  int nprocMIC = 0 ;   
  char sig0;
  float MIC2CPU_ratio = ratio;
  nRowsMIC = MIC2CPU_ratio*nRows;
  nRowsCPU = nRows-nRowsMIC;
  if (getenv("MIC_OMP_NUM_THREADS") == NULL) {
     nthreadMIC = 240; 
    if (me == 0)
     error->warning(FLERR,"MIC_OMP_NUM_THREADS environment is not set (returning to default 240).");
  }
   
#pragma offload target(mic:0) \
  nocopy(ztype_omp : length(ntypes) alloc_if(1) free_if(0)) \
  nocopy(xlocal : length(3*nlocalgroup) alloc_if(1) free_if(0)) \
  nocopy(typelocal : length(nlocalgroup) alloc_if(1) free_if(0)) \
  nocopy(store_omp : length(3*nRowsMIC) alloc_if(1) free_if(0)) \
  nocopy(Fvec : length(2*nRowsMIC) alloc_if(1) free_if(0)) signal(&sig0)
{  
  nthreadMIC = omp_get_max_threads(); 
  nprocMIC = omp_get_num_procs();
} 

  if (me == 0 && echo) {
    printf(" - nRowsMIC = %d  (Ratio=%0.2f) \n",nRowsMIC,ratio);
    printf(" - nRowsCPU = %d  (Ratio=%0.2f) \n",nRowsCPU,1-ratio);
  }
#endif

  // 5- Setup determining parameter set to use from maximum S=sin(theta)/lambda
  double Smax = Kmax / 2;
  int offset = 0;             // offset the ASFSAED matrix for appropriate value
  if (Smax > 2)  offset = 10;  
  
// ==========================================================
// End Setups
// ==========================================================

// ==========================================================
// Begin MIC Offload Region
// ==========================================================
#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_wait target(mic:0) wait(&sig0)
char signal_var;
#pragma offload target(mic:0) \
   in(ztype_omp : length(ntypes) alloc_if(0) free_if(1)) \
   in(xlocal : length(3*nlocalgroup) alloc_if(0) free_if(1)) \
   in(typelocal : length(nlocalgroup) alloc_if(0) free_if(1)) \
   in(store_omp : length(3*nRowsMIC) alloc_if(0) free_if(1)) \
   out(Fvec : length(2*nRowsMIC) alloc_if(0) free_if(1)) \
   signal(&signal_var)
{  // Start of the MIC region
#pragma omp parallel num_threads(nthreadMIC)
  {
    double *f = new double[ntypes];    // atomic structure factor by type
    int typei = 0;                     // atom type
    double Fatom1 = 0.0;               // structure factor per atom
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)

    double MY_PI  = 3.14159265358979323846; // pi
    double K[3];                       // Position in reciprocal space
    double dinv2 = 0.0;                // inverse spacing squared
    double dinv  = 0.0;                // inverse spacing
    double SinTheta_lambda  = 0.0;     // sin(theta)/lambda

    double *inners = new double[nlocalgroup];  // (2pi* dot(K,xlocal)

    if (me == 0 && echo && omp_get_thread_num() == 0) {
     printf("Inside the parallel region on MIC, there are %d openMP thread(s) available for SAED\n",omp_get_max_threads());
     printf(" nthreadsCPU = %d\n",nthreadCPU);
     printf(" nthreadsMIC = %d\n",nthreadMIC);
     printf(" wblockMIC = %d, nblockMIC = %d", wblockMIC,nblockMIC);     
    }

#pragma omp for
      for (int n = 0; n < nRowsMIC; n++) {
        int i = store_omp[3*n];
        int j = store_omp[3*n+1];
        int k = store_omp[3*n+2];
        K[0] = i * dK[0];
        K[1] = j * dK[1];
        K[2] = k * dK[2];

        dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
        dinv = sqrt(dinv2);
        SinTheta_lambda = 0.5*dinv;

        Fatom1 = 0.0;
        Fatom2 = 0.0;

        // Calculate the atomic structure factor by type
      for (int ii = 0; ii < ntypes; ii++){
        f[ii] = 0;
        for (int C = 0; C < 5; C++){
          int D = C + offset;
          f[ii] += ASFSAED[ztype_omp[ii]][D] * exp(-1 * ASFSAED[ztype_omp[ii]][5+D] * SinTheta_lambda * SinTheta_lambda);
        }
      }
      
        // Evaluate the structure factor equation -- looping over all atoms
        // Blocking the compute to remain in L1 & L2 cache
        for (int block = 0; block < nblockMIC; block++){
          iistart = block*wblockMIC;
          iistop =  iistart+wblockMIC;
          for (int ii = iistart; ii < iistop; ii++){
               inners[ii-iistart] = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                         K[2] * xlocal[3*ii+2]);
          }
          for (int ii = iistart; ii < iistop; ii++){
               typei=typelocal[ii]-1;
               Fatom1 += f[typei] * cos(inners[ii-iistart]);
               Fatom2 += f[typei] * sin(inners[ii-iistart]);
          }
        }
        // Finishing atoms not in block
        iistart = nblockMIC*wblockMIC;
        iistop =  iistart+wblockMIC;
        for (int ii = nblockMIC*wblockMIC; ii < nlocalgroup; ii++){
             inners[ii-iistart] = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                       K[2] * xlocal[3*ii+2]);
        }
        for (int ii = nblockMIC*wblockMIC; ii < nlocalgroup; ii++){
            typei=typelocal[ii]-1;
            Fatom1 += f[typei] * cos(inners[ii-iistart]);
            Fatom2 += f[typei] * sin(inners[ii-iistart]);
        }
        
        Fvec[2*n] = Fatom1;
        Fvec[2*n+1] = Fatom2;
      } // End of pragma omp for region
    delete [] f;
  } // End of pragma omp parallel region 
}  // End of MIC region
#endif 

// ==========================================================
// END MIC Offload Region
// ==========================================================


// ==========================================================
// Begin CPU Concurrent Activity Region
// ==========================================================  
  int m = 0;
  double frac = 0.1;
  double tCPU0 = MPI_Wtime();
#if defined(_OPENMP)  
#pragma omp parallel num_threads(nthreadCPU)
#endif
  {
    double *f = new double[ntypes];    // atomic structure factor by type
    int typei = 0;
    double Fatom1 = 0.0;               // structure factor per atom
    double Fatom2 = 0.0;               // structure factor per atom (imaginary)
    double K[3];
    double dinv2 = 0.0;
    double dinv  = 0.0;
    double SinTheta_lambda = 0.0;
     
    double *inners = new double[nlocalgroup];

#if defined(_OPENMP)  
    if (me == 0 && echo && omp_get_thread_num() == 0) {
       printf("Inside the CPU parallel region, there are %d openMP thread(s) available\n",omp_get_max_threads());
    }
#pragma omp for
#endif
    for (int n = nRowsMIC; n < nRows; n++) {
      int i = store_tmp[3*n];
      int j = store_tmp[3*n+1];
      int k = store_tmp[3*n+2];
      K[0] = i * dK[0];
      K[1] = j * dK[1];
      K[2] = k * dK[2];
      
      dinv2 = (K[0] * K[0] + K[1] * K[1] + K[2] * K[2]);
      dinv = sqrt(dinv2);
      SinTheta_lambda = 0.5*dinv;
              
      Fatom1 = 0.0;
      Fatom2 = 0.0;

      // Calculate the atomic structure factor by type
      // determining parameter set to use based on sin(theta)/lambda <> 2
      for (int ii = 0; ii < ntypes; ii++){
        f[ii] = 0;
        for (int C = 0; C < 5; C++){
          int D = C + offset;
          f[ii] += ASFSAED[ztype_omp[ii]][D] * exp(-1*ASFSAED[ztype_omp[ii]][5+D] * SinTheta_lambda * SinTheta_lambda);
        }
      }

        // Evaluate the structure factor equation -- looping over all atoms
        // Blocking the compute to remain in L1 & L2 cache
        for (int block = 0; block < nblockCPU; block++){
          iistart = block*wblockCPU;
          iistop =  iistart+wblockCPU;
          for (int ii = iistart; ii < iistop; ii++){
               inners[ii-iistart] = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                         K[2] * xlocal[3*ii+2]);
          }
          for (int ii = iistart; ii < iistop; ii++){
               typei=typelocal[ii]-1;
               Fatom1 += f[typei] * cos(inners[ii-iistart]);
               Fatom2 += f[typei] * sin(inners[ii-iistart]);
          }
        }
        // Finishing atoms not in block
        iistart = nblockCPU*wblockCPU;
        iistop =  iistart+wblockCPU;
        for (int ii = nblockCPU*wblockCPU; ii < nlocalgroup; ii++){
             inners[ii-iistart] = 2 * MY_PI * (K[0] * xlocal[3*ii] + K[1] * xlocal[3*ii+1] +
                       K[2] * xlocal[3*ii+2]);
        }
        for (int ii = nblockCPU*wblockCPU; ii < nlocalgroup; ii++){
            typei=typelocal[ii]-1;
            Fatom1 += f[typei] * cos(inners[ii-iistart]);
            Fatom2 += f[typei] * sin(inners[ii-iistart]);
        }

      Fvec[2*n] = Fatom1;
      Fvec[2*n+1] = Fatom2;

      // reporting progress of calculation
      if ( echo ) {
#if defined(_OPENMP)       
        #pragma omp critical
#endif        
        {
          if ( m == round(frac * nRowsCPU) ) {
            if (me == 0 && screen) fprintf(screen," %0.0f%% -",frac*100);
            frac += 0.1;
          }
          m++;
        }
      }
    } // End of pragma omp for region
    delete [] f;
    delete [] inners;
  }  // End of pragma omp parallel region on CPU
  double tCPU1 = MPI_Wtime();
  
// ==========================================================
// End CPU Concurrent Activity Region
// ==========================================================



// ==========================================================
// Gather computed data from all sources
// ==========================================================
#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_wait target(mic:0) wait(&signal_var)
#endif
  double *scratch = new double[2*nRows];
  MPI_Allreduce(Fvec,scratch,2*nRows,MPI_DOUBLE,MPI_SUM,world);

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int i = 0; i < nRows; i++) {
    vector[i] = (scratch[2*i] * scratch[2*i] + scratch[2*i+1] * scratch[2*i+1]) / natoms;
  }

/*  TESTING: this is where changes to add a "complex" flag will need to go

!!! In ADDITION !!! nRows will need to be modified to account for the extra data being saved.

  if (complex = 0) {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < nRows; i++) {
      vector[i] = (scratch[2*i] * scratch[2*i] + scratch[2*i+1] * scratch[2*i+1]) / natoms;
    }
  }
  if complex {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < nRows; i++) {
      vector[2*i] = scratch[2*i];
      vector[2*i+1] = scratch[2*i+1];
    }
  }  
   
    TESTING */
    
  delete [] xlocal;
  delete [] typelocal;  
  delete [] scratch;
  delete [] Fvec;

  // compute memory usage per processor
  double bytes = nRows * sizeof(double); //vector
  bytes +=  4.0 * nRows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * nlocalgroup * sizeof(double); // xlocal
  bytes += nlocalgroup * sizeof(int); // typelocal
  bytes += 3.0 * nRows * sizeof(int); // store_temp

  double t2 = MPI_Wtime();  
  if (me == 0 && echo) {
    if (screen)
      fprintf(screen," 100%%\nTime ellapsed during compute_xrd = %0.2f sec using %0.2f Mbytes/processor", t2-t0, bytes/1024.0/1024.0);
      fprintf(screen," \n -time ellapsed within CPU loop = %0.2f sec", tCPU1-tCPU0);
#ifdef _LMP_INTEL_OFFLOAD      
      fprintf(screen," \n -time waiting for MIC to finish = %0.2f sec\n-----\n", (t2-t0)-(tCPU1-tCPU0));
#endif      
  }
}

/* ----------------------------------------------------------------------
 memory usage of arrays
 ------------------------------------------------------------------------- */

double ComputeSAED::memory_usage()
{
  double bytes = nRows * sizeof(double); //vector
  bytes +=  4.0 * nRows * sizeof(double); //Fvec1 & 2, scratch1 & 2
  bytes += ntypes * sizeof(double); // f
  bytes += 3.0 * nlocalgroup * sizeof(double); // xlocal
  bytes += nlocalgroup * sizeof(int); // typelocal
  bytes += 3.0 * nRows * sizeof(int); // store_temp
  
  return bytes;
}

