 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: run blatn and blastx, and obtain local alignment for reduced space searching
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/blastn.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _blastn_h
#define _blastn_h

#include "sequence.h"
#include "match.h"
#include "parameters.h"
#include "tube.h"
#include "model_parameters.h"

// note : all function do not report an error if fed with or returning zero matches

Tube<int> convert_matches_to_tube(// input
				  const int    length_x,
				  const int    length_y,
				  const int    n_of_matches,
				  const Match* matches, 			
				  const int    min_radius   // minimum value of tube radius
    );

int get_blastn_results(// input
		       Sequence* sequence_x,
		       Sequence* sequence_y,
		       // output
		       int*    n_of_matches,
		       Match** matches,
		       model_parameters* const MP);

int get_tblastx_results(// input
    Sequence* sequence_x,
    Sequence* sequence_y,
    // output
    int*    n_of_matches,
    Match** matches,
    model_parameters* const MP);

int find_maximal_subset_of_matches(//input
				   const int    n_of_matches,
				   const Match* matches,
				   // output
				   int*    new_n_of_matches,
				   Match** new_matches);

int discard_blastn_matches_below_score_threshold(// input
						 const Score blastn_score_threshold,
						 // input and output
						 int* n_of_matches,
						 Score** scores_of_matches,
						 int**   x_start_positions,
						 int**   x_end_positions,
						 int**   y_start_positions,
						 int**   y_end_positions);

int get_blastn_results(// input                                   
		       Sequence* sequence_x,
		       Sequence* sequence_y,
		       // output
		       int*    n_of_matches,
		       Score** scores_of_matches,
		       int**   x_start_positions,
		       int**   x_end_positions,
		       int**   y_start_positions,
		       int**   y_end_positions,
		       model_parameters* const MP);

int get_tblastx_results(// input  
    Sequence* sequence_x,
    Sequence* sequence_y,
    // output
    int*    n_of_matches,
    Score** scores_of_matches,
    int**   x_start_positions,
    int**   x_end_positions,
    int**   y_start_positions,
    int**   y_end_positions,
    model_parameters* const MP);


int select_blastn_matches(// input 
			  const int sequences_in_same_orientation, // 1 = yes, 0 == no
			  // input & output                 
			  int* n_of_matches,
			  Score** scores_of_matches,
			  int**   x_start_positions,
			  int**   x_end_positions,
			  int**   y_start_positions,
			  int**   y_end_positions);

#endif
