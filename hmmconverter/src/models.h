/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define different pairhmm models
 */

#ifndef _models_h
#define _models_h

#include "multi_array_2.h"
#include "define.h"
#include "model_parameters.h"
#include "hmm.h"
#include "transitionprob.h"
#include "emissionprob.h"

int get_previous_states( const int state_number,
			 int**  const connection_number,
			 int*** const connection_map,
			 int* number_of_previous_states,
			 array<int>* previous_states,
			 model_parameters* const MP );

int get_next_states( const int state_number,
		     int**  const connection_number,
		     int*** const connection_map,
		     int* number_of_next_states,
		     array<int>* next_states,
		     model_parameters* const MP );


int get_connection_map(const char* filename, int** const array_connection_number,int*** const array_connection_map, model_parameters* const MP);

int set_up_topology_of_model(const char* filename,
			     Hmm* const real,
			     model_parameters* const MP);

double evaluate_transition_expression(const char* post, TransitionProb* const TP, int& check);

int set_up_transition_probs_of_model(Hmm* const real,
				     const char* filename,
				     model_parameters* const MP);


Prob sum_over(Hmm* const real,
	      EmissionProb* const EP,
	      const int From,
	      const int CurState,
	      const int index,
	      const bool state,
	      int& check);

int set_up_emission_probs_of_model(Hmm* const real,
				   const char* filename,
				   model_parameters* const MP,
				   EmissionProb* const EP);

int derive_emission_probs(Hmm* const real, 
			  model_parameters* const MP,
			  EmissionProb* const EP);

int derive_transition_probs_from_transition_parameters(Hmm* const real, 
						       model_parameters* const MP, 
						       TransitionProb* const TP);

int get_hmm(const char* filename,
    Hmm* const real,
    model_parameters* const MP,
    TransitionProb* const TP,
    EmissionProb* const EP);

char* postfix(const char* exp, int& check);

#endif

