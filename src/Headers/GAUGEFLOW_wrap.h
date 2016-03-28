/*
    Copyright 2016 Evan Solomon Weinberg

    This file (GAUGEFLOW_wrap.h) is part of GLU.

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
   @file GAUGEFLOW_wrap.h
   @brief Wrapper for the performing measurements along the flow.
 */

#ifndef GLU_GAUGEFLOW_WRAP
#define GLU_GAUGEFLOW_WRAP

/**
   @fn void GAUGEFLOW_wrap_struct( struct site *__restrict lat , const struct gaugeflow_info GAUGEFLOWINFO , const struct sm_info SMINFO )
   @brief Wrapper function

   @param lat :: Gauge fields
   @param GAUGEFLOWINFO :: gauge flow information
   @param SMINFO :: smearing information
 */
void
GAUGEFLOW_wrap_struct( struct site *__restrict lat , 
		       const struct gaugeflow_info GAUGEFLOWINFO ,
		       const struct sm_info SMINFO ) ;

#endif
