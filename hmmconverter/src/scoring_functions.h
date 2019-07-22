 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare scoring functions with which sequences can be assinged scores
 */

#ifndef _scoring_functions_h
#define _scoring_functions_h

#include <stdio.h> 
#include "multi_array_2.h" 
#include "define.h"
#include "sequence.h"
#include "hmm_state.h"
#include "hmm.h"

// for special emission
int set_scores_for_annotated_sequence(// input and output
				      Sequence* const sequence,
				      // input
				      const int            x_or_y, // 0 = x, 1 = y
				      Hmm* const pairhmm);

// for special emission
int set_scores_for_annotated_sequence_for_given_state(// input
						      Sequence*      const sequence,
						      const int            x_or_y, // 0 = x, 1 = y
						      const Hmm_State* const state);

#endif
