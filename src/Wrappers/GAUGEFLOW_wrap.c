/*
    Copyright 2016 Evan Solomon Weinberg

    This file (GAUGEFLOW_wrap.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file GAUGEFLOW_wrap.c
   @brief this wraps routines which measure observables at specific times
          along the non-adaptive Wilson flow
 */

#include "Mainfile.h"
#include "clover.h"
#include "GAUGEFLOW_wrap.h" // to declare the gaugeflow wrapper
#include "geometry.h" // for the function get_vec_from_origin
#include "GLU_timer.h" // for the timer
#include "plan_ffts.h" // config space correlator is convolution
#include "POLY.h" // for the static potential
#include "Qsusc.h" // for the topological susceptibility correlator
#include "SM_wrap.h"      // for the smearing wrapper

// what're we writing?
static GLU_bool
gaugeflow( )
{
  fprintf( stdout , "[GAUGEFLOW] " ) ;
  return GLU_TRUE ;
}

// Spit out symmE, Q to file.
void
gaugeflow_output_gauge( struct site *__restrict lat , 
                        double flow_time , 
                        const struct gaugeflow_info GAUGEFLOWINFO ,
                        const struct cut_info CUTINFO )
{
#ifndef HAVE_FFTW3_H
  fprintf( stderr , "[FFT] Need to have FFTW to get correlator.\n");
  return ;
#endif
  
  size_t i, j;
  int site [ ND ];
  double complex polyloop [ ND ]; // holds the polyakov loop. 
  char fout [ 256 ];

  fftw_plan forward , backward ; 

  // Set up outputs. 
  double qtop_sum = 0.;
  double symmE_sum = 0.;

  double *qtop_real = malloc( LVOLUME * sizeof ( double ) ) ;
  double *symmE_real = malloc( LVOLUME * sizeof ( double ) ) ;

  GLU_complex *glue = malloc ( LVOLUME * sizeof ( GLU_complex ) ) ;

  GLU_complex *glue_p = malloc ( LVOLUME * sizeof ( GLU_complex ) ) ;

  // Get gluonic quantities.
  compute_Gmunu_array_both( &qtop_sum, &symmE_sum, qtop_real, symmE_real, lat );

  fprintf( stdout , "[GLUONIC] t %.15e SymmE %.15e qtop %.15e \n", flow_time, symmE_sum, qtop_sum ) ; 

  // Init parallel threads.
  if ( parallel_ffts() == GLU_FAILURE ) {
    fprintf( stderr , "[PAR] Problem with initialising the OPENMP FFTW routines \n" ) ;
    goto memfree;
  }

  small_create_plans_DFT( &forward, &backward, glue, glue_p, ND);

  // Start doing FFTs!
  if ( GAUGEFLOWINFO.type == GFLOW_ALL || GAUGEFLOWINFO.type == GFLOW_SYMME ) {
    #pragma omp parallel for private(i)
    PFOR( i = 0; i < LVOLUME ; i++ ) {
      glue[i] = symmE_real[i];
    }

    fftw_execute( forward ) ; 
    
    // do the convolution.
    #pragma omp parallel for private(i)
    PFOR( i = 0; i < LVOLUME ; i++ ) {
      glue_p[i] *= conj( glue_p[i] ) / (double)LVOLUME ; 
    }

    fftw_execute( backward );  

    // print things to file
    sprintf( fout , "%ssymmE_wflow_%.5f.dat", CUTINFO.where, flow_time );
    FILE *outfile = fopen( fout, "w" ) ; 

    for(i=0; i<LVOLUME; i++) {
      get_vec_from_origin(site, i, ND);
      for(j=0; j<ND; j++) {
        fprintf( outfile, "%d ", site[j] );
      }
      fprintf( outfile, "%.15e\n", creal(glue[i]) );
    }
    fclose( outfile );
  }


  // Start doing FFTs!
  if ( GAUGEFLOWINFO.type == GFLOW_ALL || GAUGEFLOWINFO.type == GFLOW_TOPO ) {
    #pragma omp parallel for private(i)
    PFOR( i = 0; i < LVOLUME ; i++ ) {
      glue[i] = qtop_real[i];
    }

    fftw_execute( forward ) ; 
    
    // do the convolution.
    #pragma omp parallel for private(i)
    PFOR( i = 0; i < LVOLUME ; i++ ) {
      glue_p[i] *= conj( glue_p[i] ) / (double)LVOLUME ; 
    }

    fftw_execute( backward );  

    // print things to file
    sprintf( fout , "%sqtop_wflow_%.5f.dat", CUTINFO.where, flow_time );
    FILE *outfile = fopen( fout, "w" ) ; 

    for(i=0; i<LVOLUME; i++) {
      get_vec_from_origin(site, i, ND);
      for(j=0; j<ND; j++) {
        fprintf( outfile, "%d ", site[j] );
      }
      fprintf( outfile, "%.15e\n", creal(glue[i]) );
    }
    fclose( outfile );
  }
    
    if ( GAUGEFLOWINFO.type == GFLOW_ALL || GAUGEFLOWINFO.type == GFLOW_POLYAKOV) {
        fprintf( stdout, "[GLUONIC] t %.15e", flow_time);
        for ( i = 0; i < ND; i++ )
        {
            // as defined in src/Field/POLY.c
            polyloop[i] = poly_all( lat, i ) ; 
            fprintf( stdout, " poly%zur %.15e poly%zui %.15e", i, creal(polyloop[i])/(NC*LVOLUME), i, cimag(polyloop[i])/(NC*LVOLUME));
        }
        
        
        fprintf( stdout, "\n"); 
    }

memfree:
  free(qtop_real);
  free(symmE_real);
  free(glue);
  free(glue_p);

}


// wrapper for the various measurements we want done along the flow
void
GAUGEFLOW_wrap_struct( struct site *__restrict lat ,
                       const struct gaugeflow_info GAUGEFLOWINFO ,
                       const struct cut_info CUTINFO , 
                       const struct sm_info SMINFO )
{
  size_t i, j;
  size_t counter = 0;
  size_t total_iters;
  struct sm_info SMINFO_mod;
  size_t sort_steps[ 256 ];
  start_timer( ) ;
  

  // Sort the times we want to measure at.
  for (i = 0; i < GAUGEFLOWINFO.nmeas; i++)
  {
    sort_steps[i] = GAUGEFLOWINFO.meassteps[i];
  }

  // Just do a lazy bubble sort.
  {
    for (i = 0; i < GAUGEFLOWINFO.nmeas; i++)
    {
      for (j = 0; j < GAUGEFLOWINFO.nmeas-1; j++)
      {
        if (sort_steps[j] > sort_steps[j+1])
        {
          size_t tmp = sort_steps[j];
          sort_steps[j] = sort_steps[j+1];
          sort_steps[j+1] = tmp;
        }
      }
    }
  }

  // Copy the smearing info to something mutable.
  SMINFO_mod.dir = SMINFO.dir;
  SMINFO_mod.smiters = SMINFO.smiters;
  SMINFO_mod.type = SMINFO.type;


  // What if the first time is zero?
  if ( sort_steps[counter] == 0)
  {
    gaugeflow();
    printf("Doing a measurement at iter %zu, step %zu.\n", counter, sort_steps[counter]);

    gaugeflow_output_gauge( lat , 0.0 , GAUGEFLOWINFO, CUTINFO );

    counter++;
  }

  // What's the total number of iterations?
  total_iters = SMINFO.smiters; 

  while (counter < GAUGEFLOWINFO.nmeas) {
    // How far should we flow?
    if (counter == 0) { // then the first time wasn't 0.
      SMINFO_mod.smiters = sort_steps[counter];
    }
    else
    {
      SMINFO_mod.smiters = sort_steps[counter] - sort_steps[counter-1];
    }
    
    // Flow the amount we need to!
    SM_wrap_struct( lat , SMINFO_mod ) ;
    
    gaugeflow();
    printf("Doing a measurement at iter %zu, step %zu.\n", counter, sort_steps[counter]);

    gaugeflow_output_gauge( lat , Latt.sm_alpha[0]*sort_steps[counter] , GAUGEFLOWINFO, CUTINFO );

    counter++;
  }

  // Flow the rest of the way.
  SMINFO_mod.smiters = total_iters - sort_steps[counter-1];
  SM_wrap_struct( lat , SMINFO_mod ) ;



  ///////////// CREATE CUTS ////////////
  /*
  switch( CUTINFO.dir ) {
  case INSTANTANEOUS_GLUONS :
    cuts_spatial( lat , CUTINFO ) ;
    break ;
  case CONFIGSPACE_GLUONS :
    cuts_struct_configspace( lat , CUTINFO , SMINFO ) ;
    break ;
  case SMEARED_GLUONS :
    cuts_struct_smeared( lat , CUTINFO , SMINFO ) ;
    break ;
  case STATIC_POTENTIAL :
    Coul_staticpot( lat , CUTINFO , SMINFO ) ;
    break ;
  case TOPOLOGICAL_SUSCEPTIBILITY :
    compute_Qsusc( lat , CUTINFO , SMINFO ) ;
    break ;
  case GLUON_PROPS :
  case EXCEPTIONAL :
  case NONEXCEPTIONAL :
  case FIELDS :
    // all of the other momentum space routines are in here
    cuts_struct( lat , CUTINFO ) ;
    break ;
  default :
    return ;
  }
  */

  print_time( ) ;

  return ;
}
