/*
    Copyright 2013 Renwick James Hudspith

    This file (adaptive_flow.c) is part of GLU.

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
   @file adaptive_flow.c
   @brief the (two step) adaptive rk4 wilson flow routine

   Slows down, performing fine measurements at ~t_0 and ~w_0
   WFLOW_STOP and TMEAS_STOP are defined in wflowfuncs.h
 */

#include "Mainfile.h"

#include "geometry.h"     // init_navig is called for the temporary
#include "plaqs_links.h"  // clover and plaquette measurements
#include "wflowfuncs.h"   // wilson flow general routines

// enable this if we are doing a straight shot for a specific time
// to avoid measuring the topological charge and clover and all that,
// be careful when T > 10, will need to increase ADAPTIVE_EPS
//#define TIME_ONLY

static const double mnineOseventeen = -0.52941176470588235294 ; // -9.0/17.0

/**
   @fn static inline double adaptfmax( const double a , const double b )
   @brief the maximum of two numbers
   Is probably unsafe because Infs and Nans are not considered, although
   what would you do with them?
 */
static inline double
adaptfmax( a , b )
     const double a ;
     const double b ;
{
  return ( b < a ? a : b ) ;
}

/**
   @fn static inline double adaptfmin( const double a , const double b )
   @brief the minimum of two numbers
   Again probably unsafe
 */
static inline double
adaptfmin( a , b )
     const double a ;
     const double b ;
{
  return ( a < b ? a : b ) ;
}

/**
   @enum adaptive_control
   @brief when to break our adaptive algorithm if we have done this many halvings and still have no result
 */
enum adaptive_control{ ADAPTIVE_BIG_NUMBER = 20 } ;

// Adaptive stepsize version 
int 
flow4d_adaptive_RK( struct site *__restrict lat , 
		    const int smiters ,
		    const int DIR ,
		    const int SIGN ,
		    const int SM_TYPE )
{  
  ////// USUAL STARTUP INFORMATION /////////
  print_GG_info( SM_TYPE , RK4_ADAPTIVE ) ;

  // Martin Luescher's choice for the variable t_0
  static const double T0_STOP = 0.3 ;

  // the error between the two plaquettes
  const double ADAPTIVE_EPS = 2.5E-6 ;
  // Standard shrink and factor from NRC
  const double ADAPTIVE_SHRINK = -0.25 ;
  // Standard growth and factor from NRC
  const double ADAPTIVE_GROWTH = -0.20 ;
  // define adaptive safe
  const double ADAPTIVE_SAFE = 0.9 ;
  // adaptive error conserving
  const double ADAPTIVE_ERRCON = powl( 5./ADAPTIVE_SAFE , 1./ADAPTIVE_GROWTH ) ;

  printf( "[WFLOW] Adaptive Error :: %e \n" , ADAPTIVE_EPS ) ;
  printf( "[WFLOW] Adaptive ErrCon :: %f \n" , ADAPTIVE_ERRCON ) ;
  printf( "[WFLOW] Adaptive Safety Factor :: %g \n" , ADAPTIVE_SAFE ) ;
  printf( "[WFLOW] Adaptive growth factor :: %g \n" , ADAPTIVE_GROWTH ) ;
  printf( "[WFLOW] Adaptive shrink factor :: %g \n\n" , ADAPTIVE_SHRINK ) ; 

  //////////////////////////////////////////

  struct spt_site_herm *Z = malloc( LVOLUME * sizeof ( struct spt_site_herm ) ) ;
  struct spt_site *lat2 = malloc( LCU * sizeof ( struct spt_site ) ) ;
#ifdef IMPROVED_SMEARING
  struct spt_site *lat3 = malloc( 2*LCU * sizeof ( struct spt_site ) ) ;
  struct spt_site *lat4 = malloc( 2*LCU * sizeof ( struct spt_site ) ) ;
#else
  struct spt_site *lat3 = malloc( LCU * sizeof ( struct spt_site ) ) ;
  struct spt_site *lat4 = malloc( LCU * sizeof ( struct spt_site ) ) ;
#endif

  // Set up a temporary lattice of the semi-steps
  struct site *lat_two = malloc( LVOLUME * sizeof ( struct site ) ) ;
  init_navig( lat_two ) ;

  // set up the step sizes ...
  double delta_t = SIGN * Latt.sm_alpha[0] , yscal = 0. , qtop_out = 0. , t = 0.0 ;

  lattice_gmunu( lat , &qtop_out , &yscal ) ;
  printf("[WFLOW] {err} %1.3e {t} %f {dt} %e {p} %1.10f {q} %1.10f {ttGG} 0. {W} 0. \n" ,
	 0.0 , t , delta_t , yscal , qtop_out ) ;

  // counters for the derivative ...
  double flow = 0. , flow_next = 0. ;
  int count = 0 , OK_STEPS = 0 , NOTOK_STEPS = 0 ; 

  // stop us from stopping the flow too early ...
#ifndef verbose
  int first_deriv_flag = GLU_FALSE ;
#endif

  // set up a control for the measurement to stop us from running over the desired value
  // We then step cautiously through to W0, ensuring that at least 8 "fine" measurements
  // around WFLOW_STOP have been made so that we can fit linearly to get W0
  const double WFLOW_MEAS_CONTROL = WFLOW_STOP * 0.1 ;
  int T0_MEAS = GLU_FALSE ;

  for( count = 1 ; count <= smiters ; count++ ) { 
    int counter = 0 ;
    double errmax = 10. ;
    double new_plaq = 0. ;
    while( ( errmax > 1.0 ) && ( counter < ADAPTIVE_BIG_NUMBER ) ) {
      int i ;
      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	memcpy( &lat_two[i] , &lat[i] , sizeof( struct site ) ) ; 
      }

      // Step forward in two halves ...
      const double rk1 = mnineOseventeen * delta_t ;
      const double rk2 = delta_t ;
      const double rk3 = ( -delta_t ) ;  

      // step forward once and write into lat_two
      step_distance_memcheap( lat_two , lat2 , lat3 , lat4 , Z , 
			      rk1 , rk2 , rk3 , SM_TYPE ) ;
      
      // compute the one-step first comparison
      const double old_plaq = av_plaquette( lat_two ) ;

      #pragma omp parallel for private(i)
      PFOR( i = 0 ; i < LVOLUME ; i++ ) {
	memcpy( &lat_two[i] , &lat[i] , sizeof( struct site ) ) ; 
      } 

      // and step forward twice and write into lat_two
      const double half_rk1 = 0.5 * rk1 ;
      const double half_rk2 = 0.5 * rk2 ;
      const double half_rk3 = 0.5 * rk3 ;

      step_distance_memcheap( lat_two , lat2 , lat3 , lat4 , Z , 
			      half_rk1 , half_rk2 , half_rk3 , 
			      SM_TYPE ) ;

      step_distance_memcheap( lat_two , lat2 , lat3 , lat4 , Z , 
			      half_rk1 , half_rk2 , half_rk3 , 
			      SM_TYPE ) ;
	  
      // compute the error I will use the average plaquette ...
      new_plaq = (double)av_plaquette( lat_two ) ;

      errmax = fabsl( ( new_plaq - old_plaq ) / ( yscal ) ) ;
      errmax /= ADAPTIVE_EPS ;

      // Break the while loop if conditions are satisfied
      if( errmax < 1.0 ) {
	if( counter < 1 ) { OK_STEPS ++ ; } 
	break ;
      }
      // Increment the counter for not adequate steps
      NOTOK_STEPS ++ ;

      // shorten the delta_t ...
      const double del_temp =  ADAPTIVE_SAFE * delta_t * pow( errmax , ADAPTIVE_SHRINK ) ; 
      // set up a tolerance s.t del_temp doesn't go too crazy, only really used when starting guess is bad
      const double tol = 0.1 * delta_t ; 
      delta_t = ( 0. < del_temp ? adaptfmax( del_temp , tol ) : adaptfmin( del_temp , tol ) ) ;

      // Print if we are in trouble
      if( counter == ADAPTIVE_BIG_NUMBER - 1 ) {
	printf( "[WFLOW] Not stepping to required accuracy after < %d > attempts! \n" , counter ) ;
      }

      // Increment our counter ...
      counter ++ ;
    }

    // set up a scaling parameter to control the adaptation uses a first order finite difference def ...
    const double yscal_new = new_plaq ;
    yscal = 2.0 * yscal_new - yscal ;

    // If we get stuck in a rut of updating by zero we leave
    if( fabs( delta_t ) < DBL_MIN ) {
      printf( "[WFLOW] No update made delta_t :: %1.5e \nLeaving ... \n" , delta_t ) ;  
      return GLU_FAILURE ;
    }

    // rewrite lat .. 
    int i ;
    #pragma omp parallel for private(i)
    PFOR( i = 0 ; i < LVOLUME ; i++ ) {
      memcpy( &lat[i] , &lat_two[i] , sizeof( struct site ) ) ; 
    }

    t += delta_t ; // add one time step to the overall time 

#ifndef TIME_ONLY
    // derivative bit, should be used as a guide for when to stop
    // measureing for W_0 but not used to measure it, 
    // use ttGG and take the derivative!
    double deriv = 0. ;
    #ifndef verbose
    if( t > MEAS_START ) {
    #endif
      printf( "[WFLOW] {err} %1.3e " , errmax * ADAPTIVE_EPS ) ;
      deriv = deriv_euler( lat , &flow , &flow_next , t , delta_t ) ;
      #ifndef verbose
      if( first_deriv_flag == GLU_FALSE )  {
	first_deriv_flag = GLU_TRUE ;
	deriv = 0. ;
      }
      #endif
    #ifndef verbose
    }
    #endif

    // If we are close to our goal we reduce the step size
    if( ( WFLOW_STOP - deriv ) <= WFLOW_MEAS_CONTROL ) {
      delta_t = 0.01 ; // Set this to something nice and low
      // put this little routine here so that we do not do the half step again
      const double rk1 = mnineOseventeen * delta_t ;
      const double rk2 = delta_t ;
      const double rk3 = ( -delta_t ) ;  
      int fine_meas = 0 ;
      while( ( deriv < WFLOW_STOP ) || ( fine_meas < 8 ) )  {
	t += delta_t ; // add one time step to the overall time 
	// regurgitate the previous error, will be much smaller!
	printf( "[WFLOW W0] {err} %1.3e " , errmax * ADAPTIVE_EPS ) ;
	step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
				Z , rk1 , rk2 , rk3 , SM_TYPE ) ;
	deriv = deriv_euler( lat , &flow , &flow_next , t , delta_t ) ;
	flow = flow_next ;
	fine_meas++ ;
      }
      // leave the big loop over smearing iterations
      break ;
      // we have one for t_0 too
    } else if ( ( ( T0_STOP - flow_next ) <= T0_STOP * 0.1 ) && 
		T0_MEAS == GLU_FALSE ) {
      const double prev_dt = delta_t ;
      delta_t = 0.01 ;
      const double rk1 = mnineOseventeen * delta_t ;
      const double rk2 = delta_t ;
      const double rk3 = ( -delta_t ) ;  
      int fine_meas = 0 ;
      while( ( flow_next < T0_STOP ) || ( fine_meas < 8 ) )  {
	t += delta_t ;
	printf( "[WFLOW T0] {err} %1.3e " , errmax * ADAPTIVE_EPS ) ;
	step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
				Z , rk1 , rk2 , rk3 , SM_TYPE ) ;
	deriv = deriv_euler( lat , &flow , &flow_next , t , delta_t ) ;
	flow = flow_next ;
	fine_meas++ ;
      }
      T0_MEAS = GLU_TRUE ;
      delta_t = prev_dt ;
      // do not break from here ...
    } else {
      // Increase the step size ...
      if( errmax > ADAPTIVE_ERRCON ) {
	delta_t = ADAPTIVE_SAFE * delta_t * pow( errmax , ADAPTIVE_GROWTH ) ;
      } else {
	delta_t = ADAPTIVE_SAFE * 5.0 * delta_t ;
      }
    }
#else
    printf( "[WFLOW] {err} %1.3e {t} %f {dt} %g\n" , errmax * ADAPTIVE_EPS , t , delta_t ) ;
    if( errmax > ADAPTIVE_ERRCON ) {
      delta_t = ADAPTIVE_SAFE * delta_t * pow( errmax , ADAPTIVE_GROWTH ) ;
      // be too large for the exponentiation to handle ...
      // delta_t = ( delta_t > 0.4 ) ? 0.4 : delta_t ;
    } else {
      delta_t = ADAPTIVE_SAFE * 5.0 * delta_t ;
    }
#endif

    // If we are wilson-flowing to a specific time, we compute a negative flow-time correction
    if( t > TMEAS_STOP ) {
      const double delta_tcorr = TMEAS_STOP - t ; 
      t = TMEAS_STOP ;
      const double rk1 = mnineOseventeen * delta_tcorr ;
      const double rk2 = delta_tcorr ;
      const double rk3 = ( -delta_tcorr ) ;  
      step_distance_memcheap( lat , lat2 , lat3 , lat4 , 
			      Z , rk1 , rk2 , rk3 , SM_TYPE ) ;
      // regurgitate the previous error
      printf( "[WFLOW] {err} %1.3e " , errmax * ADAPTIVE_EPS ) ;
      deriv_euler( lat , &flow , &flow_next , t , delta_tcorr ) ;
      break ;
    }
   // end of RK step
  }

  // Print out the stepping information
  printf( "\n[WFLOW] Inadequate steps :: %d \n" , NOTOK_STEPS ) ;
  printf( "[WFLOW] Adequate steps :: %d \n" , OK_STEPS ) ;

  // free our fields
  free( Z ) ;
  free( lat2 ) ;
  free( lat3 ) ;
  free( lat4 ) ;
  free( lat_two ) ;

  return GLU_SUCCESS ;
}

// if we have set the code in "TIME_ONLY" mode we make sure we clean it up
#ifdef TIME_ONLY
  #undef TIME_ONLY
#endif
