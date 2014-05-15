/*
    Copyright 2013 Renwick James Hudspith

    This file (cut_output.c) is part of GLU.

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
   @file cut_output.c
   @brief output functions for the cutting procedures
 */

#include "Mainfile.h"

#include "geometry.h"  // geometry for computing momenta
#include "GLU_bswap.h" // byte swapping

// big endian file writers, set up to mimic the c-standard fwrite but cannot
// use a multiple of the array length for the size value
static void
FWRITE( void *arr ,
	const int size , // MUST BE sizeof(int) or sizeof(double)
	const int stride ,
	FILE *__restrict file )
{
  // I enforce the re byte-swap in case we want to look at the data again
  if( size == sizeof(int) ) {
    if( !WORDS_BIGENDIAN ) { bswap_32( stride , arr ) ; }
    fwrite( arr , size , stride , file ) ;
    if( !WORDS_BIGENDIAN ) { bswap_32( stride , arr ) ; }
  } else if ( size == sizeof( double ) ) {
    if( !WORDS_BIGENDIAN ) { bswap_64( stride , arr ) ; }
    fwrite( arr , size , stride , file ) ;
    if( !WORDS_BIGENDIAN ) { bswap_64( stride , arr ) ; }
  } else {
    printf( "[CUT] I do not understand your byte-swapping size %d \n" , size ) ;
  }
  return ;
}

// gives us some nominal information on the cutting procedure
int
check_psq( const struct cut_info CUTINFO )
{
  int small = Latt.dims[0] ;
  int i ;
  for( i = 1 ; i < ND ; i++ ) {
    if( Latt.dims[i] < small ) {
      small = Latt.dims[i] ;
    }
  }

  // Type is the type of cut we are using
  // 0 is normal hypercubic
  // 1 is normal psq
  // 2 is a cylinder cut
  // 3 is with an angle

  printf( "\n[CUTS] " ) ;
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    printf( "PSQ Cutting Procedure \n" ) ;
  } else {
    switch( CUTINFO.type ) 
      {
      case HYPERCUBIC_CUT :
	printf( "HYPERCUBIC Cutting Procedure \n" ) ;
	break ;
      case PSQ_CUT :
	printf( "PSQ Cutting Procedure \n" ) ;
	break ;
      case CYLINDER_CUT :
	printf( "Cylinder Cutting Procedure \n" ) ;
	printf( "[CUTS] Cylinder Width in Lattice units :: %g \n" , CUTINFO.cyl_width ) ; 
	printf( "[CUTS] Physical Momentum Cap :: %g \n" , TWOPI*CUTINFO.cyl_width/small  ) ; 
	break ;
      case CYLINDER_AND_CONICAL_CUT :
	printf( "Cylinder and Conical Cutting Procedure \n" ) ;
	printf( "[CUTS] Cylinder Radius in Lattice units :: %g \n" , CUTINFO.cyl_width ) ; 
	printf( "[CUTS] Physical Momentum Cap :: %g \n" , TWOPI*CUTINFO.cyl_width/small ) ; 
	printf( "[CUTS] Cone Angle :: %g \n" , (double)CUTINFO.angle ) ; 
	break ;
      default:
	printf( "I don't recognise the type.. Leaving \n" ) ;
	return GLU_FAILURE ;
      }
  }
  printf( "[CUTS] Maximum p^2 :: %d \n\n" , CUTINFO.max_mom ) ;
  if( CUTINFO.definition == LOG_DEF ) {
    printf( "[CUTS] Using the LOGARITHMIC field definition \n" ) ;
  } else {
    printf( "[CUTS] Using the LINEAR field definition \n" ) ;
  }

  return GLU_SUCCESS ;
}

// Decides the shape of our output string
char *
output_str_struct( const struct cut_info CUTINFO )
{
  char *str = malloc( 1024 * sizeof ( char ) ) ;
  
  char tmp[ 256 ] ;
  sprintf( tmp , "%s" , CUTINFO.where ) ;
  #ifndef CONDOR_MODE
  printf( "[CUTS] Our file is located @ %s \n", CUTINFO.where ) ;
  #else
  sprintf( tmp , "./" ) ;
  #endif

  // each mode has a slightly different start so that you know what you are getting into 
  switch( CUTINFO.dir )
    {
    case INSTANTANEOUS_GLUONS :
      sprintf( str , "%sG2spG2t" , tmp ) ;
      break ;
    case SMEARED_GLUONS :
      sprintf( str , "%sG2_G2SM" , tmp ) ;
      break ;
    case GLUON_PROPS :
      sprintf( str , "%sG2F2" , tmp ) ; 
      break ;
    case CONFIGSPACE_GLUONS :
      sprintf( str , "%sCSPACE" , tmp ) ; 
      break ;
    case FIELDS :
      sprintf( str , "%sMATS" , tmp) ; 
      break ;
    case EXCEPTIONAL :
      sprintf( str , "%sMOMgg" , tmp ) ; 
      break ;
    case NONEXCEPTIONAL :
      #ifdef PROJ_GRACEY
      sprintf( str , "%spgracey", tmp ) ;
      #else
      sprintf( str , "%spboucaud", tmp ) ;
      #endif
      break ;
    case STATIC_POTENTIAL :
      sprintf( str , "%sSTATPOT" , tmp ) ;
      break ;
    case TOPOLOGICAL_SUSCEPTIBILITY :
      sprintf( str , "%sQSUSC" , tmp ) ;
      break ;
    default : // default behaviour is the 3 gluon vertex
      sprintf( str , "%sMOMgg" , tmp ) ; 
      break ;
    }
  // color info
  sprintf( str , "%s_SU%d_" , str , NC ) ;

  // print out the dimensions of the lattice used for the cut
  int mu ;
  for( mu = 0 ; mu < ND - 1 ; mu++ ) {
    sprintf( str , "%s%dx" , str , Latt.dims[mu] ) ;
  }
  sprintf( str , "%s%d_" , str , Latt.dims[ ND - 1 ] ) ;

  if( CUTINFO.dir == CONFIGSPACE_GLUONS ) {
    if( CUTINFO.definition == LOG_DEF ) {
      sprintf( str , "%sLOG.%d.bin" , str , Latt.flow ) ;
    } else {
      sprintf( str , "%sLIN.%d.bin" , str , Latt.flow ) ;
    }
    return str ;
  }

  /////// MOMENTUM SPACE STUFF FROM HERE ///////
  sprintf( str , "%s%d_" , str , CUTINFO.max_mom ) ;

  // static potential measurement has T-detail here
  if( CUTINFO.dir == STATIC_POTENTIAL ) { 
    sprintf( str , "%sT%d_" , str , CUTINFO.max_t ) ;
  }

  // switch on the correct cut, only PSQ is allowed for the non-exceptional
  momentum_cut_def type = CUTINFO.type ;
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    type = PSQ_CUT ;
  }

  // print out the cutting type performed
  switch( type )
    {
    case HYPERCUBIC_CUT :
      sprintf( str , "%scutHYP" , str ) ;
      break;
    case PSQ_CUT :
      sprintf( str , "%scutPSQ" , str ) ;
      break;
    case CYLINDER_CUT :
      sprintf( str , "%scutCYL_w%1.2f" , str , CUTINFO.cyl_width ) ;
      break ;
    case CYLINDER_AND_CONICAL_CUT :
      sprintf( str , "%scutCON_%d" , str , CUTINFO.angle ) ;
      break ; 
    default:
      break ; 
    }

  // print out the momentum def used for the projectors, only used for
  // the three point functions
  if( CUTINFO.dir != STATIC_POTENTIAL && 
      CUTINFO.dir != TOPOLOGICAL_SUSCEPTIBILITY ) {
    #ifdef SIN_MOM
    sprintf( str , "%s_SINMOM" , str ) ;
    #else
    sprintf( str , "%s_PSQMOM" , str ) ;
    #endif
  }

  // have PROJ_GRACEY == 0 for the 6 combined
  #ifdef PROJ_GRACEY
  if( CUTINFO.dir == NONEXCEPTIONAL ) {
    #if PROJ_GRACEY > -1 && PROJ_GRACEY < 15
    sprintf( str , "%s_PROJ%d" , str , PROJ_GRACEY ) ;
    #endif
  }
  #endif
  
  // field definition for momentum space glue
  if( CUTINFO.dir != STATIC_POTENTIAL && 
      CUTINFO.dir != TOPOLOGICAL_SUSCEPTIBILITY ) {
    if( CUTINFO.definition == LOG_DEF ) {
      sprintf( str , "%s_LOG" , str ) ;
    } else {
      sprintf( str , "%s_LIN" , str ) ;
    }
  }

  // and finish up with the configuration number
  sprintf( str , "%s.%d.bin" , str , Latt.flow ) ;

  return str ;
}

// Write the gluonic two point function to the file Ap
void
write_g2_to_list( FILE *__restrict Ap , 
		  double *__restrict g2 , 
		  int num_mom[ 1 ] ) 
{
#ifdef ASCII_CHECK
  int i ;
  printf( "%d\n" , num_mom[0] ) ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { printf("%e \n", g2[i] ) ; }
  printf("%d \n", num_mom[ 0 ] ) ;
#else
  // write NUM_MOM[0] , G2
  FWRITE( num_mom , sizeof(int) , 1 , Ap ) ;
  FWRITE( g2 , sizeof(double) , num_mom[0] , Ap ) ;
#endif
  return ;
}

// Write the gluonic two and three point functions to the file Ap
void
write_g2g3_to_list( FILE *__restrict Ap , 
		    double *__restrict g2 , 
		    double *__restrict g3 , 
		    int num_mom[ 1 ] ) 
{
#ifdef ASCII_CHECK
  printf("%d \n", num_mom[ 0 ] ) ;
  int i ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { printf("%e \n", g2[i] ) ; }
  printf("%d \n", num_mom[ 0 ] ) ;
  for( i = 0 ; i < num_mom[0] ; i++ ) { printf("%e \n", g3[i] ) ; }
#else
  // write NUM_MOM[0] , G2 , NUM_MOM[0] , G3
  FWRITE( num_mom , sizeof(int) , 1 , Ap ) ;
  FWRITE( g2 , sizeof(double) , num_mom[0] , Ap ) ;
  FWRITE( num_mom , sizeof(int) , 1 , Ap ) ;
  FWRITE( g3 , sizeof(double) , num_mom[0] , Ap ) ;
#endif
  return ;
}

void
write_lattice_fields( FILE *__restrict Ap , 
		      const struct site *__restrict A , 
		      const struct veclist *__restrict list , 
		      int num_mom[1] )
{
  // write the length of the array
  FWRITE( num_mom , sizeof(int) , 1 , Ap ) ;

  //write lattice fields
  const int stride = ND * NCNC * 2 ; // the 2 is for the real and imaginary bits
  double *aout = malloc( stride * sizeof( double ) ) ;

  int i ;
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    int mu ;
    //int b = 0 ;
    #pragma omp parallel for private(mu)
    PFOR( mu = 0 ; mu < ND ; mu++ ) {
      int a ;
      for( a = 0 ; a < NCNC ; a++ ) {
	const int b = 2 * ( a + mu * NCNC ) ;
	aout[ b ] = (double)creal( A[ list[i].idx ].O[mu][a] ) ;
	aout[ b + 1 ] = (double)cimag( A[ list[i].idx ].O[mu][a] ) ;
      }
    }
    FWRITE( aout , sizeof(double) , stride , Ap ) ;
  }
  free( aout ) ;
  return ;
}

// Support for writing our momentum list
void
write_mom_veclist( FILE *__restrict Ap , 
		   int *__restrict num_mom , 
		   const struct veclist *__restrict list ,
		   const int DIR )
{
  const int stride = ( DIR + 1 ) * num_mom[ 0 ] ;
  int *kall = (int*)malloc( stride * sizeof(int) ) ;

  FWRITE( num_mom , sizeof(int) , 1 , Ap ) ;

  int i ;
  //write momenta
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < num_mom[0] ; i++ ) {
    int n = i * ( DIR + 1 ) ;
    kall[ n ] = DIR ; 
    n++ ; 
    int a ;
    for( a = 0  ;  a < DIR  ;  a++ ) {
      kall[ n ] = list[i].MOM[a] ; 
      n++ ; 
    }
  }
  FWRITE( kall , sizeof(int) , stride , Ap ) ;
  free( kall ) ;
  return ;
}

// write out rsq's for the correlator
void
write_rr_values( FILE *__restrict Ap ,
		 int size[1] ,
		 const int *__restrict rsq ,
		 const int max_r2 ,
		 const int ARR_SIZE )
{
  int i ;
#ifdef ASCII_CHECK
  printf( "%d \n" , size[0] ) ;
  // loop the possible rsq's
  for( i = 0 ; i < ARR_SIZE ; i++ ) {
    if( rsq[i] < max_r2 ) {
      int rr[1] = { rsq[i] } ;
      printf( "%d\n" , rr[0] ) ;
    }
  }
#endif
  // output to a file ... for reading in with UKhadron or something
  // write out the length of the array ...
  FWRITE( size , sizeof(int) , 1 , Ap ) ;
  // loop the possible rsq's
  for( i = 0 ; i < ARR_SIZE ; i++ ) {
    if( rsq[i] < max_r2 ) {
      int rr[1] = { rsq[i] } ;
      FWRITE( rr , sizeof(int) , 1 , Ap ) ;
    }
  }
}

// Support for writing our triplets momentum
void
write_triplet_mom_list( FILE *__restrict Ap , 
			int *__restrict num_mom , 
			int *__restrict *__restrict momentum ,
			int *__restrict *__restrict triplet ,
			const int DIR )
{
  const int stride = ( ND + 1 ) * num_mom[0] ;
  int *kall = malloc( stride * sizeof(int) ) ;

#ifdef ASCII_CHECK
  printf("%d\n", num_mom[0] ) ;
#endif

  FWRITE( num_mom , sizeof(int) , 1 , Ap ) ;

  // start writing the momenta to the list, we just write the first momentum of the triplet...
  int i ;
#pragma omp parallel for private(i)
  PFOR( i = 0 ; i < num_mom[0] ; i++ ) {
    //const int trip = i % 3 ;
    int n = i * ( ND + 1 ) ;
    kall[ n ] = ND ; 
    n++ ; 
    int mu ;
    for( mu = 0 ; mu < ND ; mu++ ) {
      kall[n] = momentum[ triplet[ i ][ 0 ] ][ mu ] ;
      n++ ;
    }
  }

#ifdef ASCII_CHECK
  int mu ;
  int b = 0 ;
  for( i = 0 ; i < num_mom[0] ; i++ ) {
    for( mu = 0 ; mu < ND + 1 ; mu++ ) {
      printf("%d ", kall[ b ] ) ; 
      b++ ;
    }
    printf( "\n" ) ;
  }
  printf("%d \n", kall[b] ) ;
  printf("%d \n", num_mom[0] ) ;
#endif
  FWRITE( kall , sizeof(int) , stride , Ap ) ;
  free( kall ) ;
  return ;
}

// write the timeslices, LT[0] is the size of the array
void
write_tslice_list( FILE *__restrict Ap , 
		   int *__restrict LT )
{
  int T[ LT[0] ] , t ;
#ifdef ASCII_CHECK
  printf( "%d \n" , LT[0] ) ;
  for( t = 0 ; t < LT[0] ; t++ ) { printf( "%d \n" , t ) ; }
  printf( "%d \n" , LT[0] ) ;
#endif
  FWRITE( LT , sizeof(int) , 1 , Ap ) ;
  for( t = 0 ; t < LT[0] ; t++ ) { T[ t ] = t ; }
  FWRITE( T , sizeof(int) , LT[0] , Ap ) ;
  return ;
}
