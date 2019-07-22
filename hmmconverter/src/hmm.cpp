/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: header-file for pairhmm class

   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/hmm.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
*/
#include <fstream>
#include <iostream>
#include <string.h>

#ifndef _INTEL
#include <vector.h>
#include <vector>
#include <map.h>
#endif

#ifdef _SGI
#include <strstream>
#endif
#ifdef _ACARI
#include <strstream>
#endif
#ifndef _SGI

#ifndef _ACARI
#include <backward/strstream>
#endif

#endif
#include <math.h>
#include "tube.h"
#include "hmm.h"
#include "sequence.h"
#include "hmm_state.h"

Hmm::Hmm()
{
    // this constructur shall not be used
    model=NULL;
    mirrored=0;
    alphabet=0;
    number_of_states=0;
    steps=0;
    score=0;
    sequence_of_scores=NULL;
    sequence_of_states=NULL;
    sequence_of_xsteps=NULL;
    sequence_of_ysteps=NULL;

    condensed_steps=0;
    condensed_sequence_of_scores=NULL;
    condensed_sequence_of_states=NULL;    
    condensed_sequence_of_xsteps=NULL;
    condensed_sequence_of_ysteps=NULL;
}

Hmm::Hmm(model_parameters* const MP)
{

    mirrored=0;
    alphabet= MP->get_Alphabet_size();
    number_of_states= MP->get_Number_of_States();
    
    model = new Hmm_State[number_of_states];
    
//  for (int j=0; j<i; j++)
//    {model[j]=NULL;}

    steps=0;
    score=0;
    sequence_of_scores=NULL;
    sequence_of_states=NULL;
    sequence_of_xsteps=NULL;
    sequence_of_ysteps=NULL;
    condensed_steps=0;
    condensed_sequence_of_scores=NULL;
    condensed_sequence_of_states=NULL;    
    condensed_sequence_of_xsteps=NULL;
    condensed_sequence_of_ysteps=NULL;
}

Hmm::Hmm(const Hmm &p)
{
    mirrored=0;
    alphabet=0;
    number_of_states=0;
    model=NULL;
    
    *this=p;
}

Hmm::Hmm(const Hmm &p, int mirror, int change_strand)
{
    // NOTE: mirror has to be 0 (do not mirror) or -1 (mirror)
    //       change_strand has to be 0 (do not change strand) or 1 (yes, do change strand)
    
    mirrored=0;
    alphabet=0;
    number_of_states=0;
    model=NULL;
    
    if ((mirror==-1)  &&    // if mirror of p is to be created
	(this != &p))       // if p is not equal to what we try to create
    {
	mirrored=(p.mirrored+1)%2; // make sure that mirrored is 0 or 1
	alphabet=p.alphabet;
	number_of_states=p.number_of_states;
	
	if (model) delete [] model;
	model = NULL;
	model = new Hmm_State[number_of_states];
	
	int x_or_y_indices_first = 1;

	if (p.mirrored == 0) {
	    x_or_y_indices_first = -1;
	}
      
	{      
	    for (int i=0; i<number_of_states; i++)
	    {
		// if (change_strand == 1):
		// permute emission probs such that (effect of a + b + c) :
		//
		// prob of old emission indices: (x_1,x_2,x_3,y_1,y_2,y_3) = 
		// prob of new emission indices: (x_3',x_2',x_1',y_3',y_2',y_1')
		
		// a) mirror state
		//    emission indices: (x_1,x_2,x_3,y_1,y_2,y_3) => (y_3,y_2,y_1,x_3,x_2,x_1)
		
		model[i].get_mirrored_copy_of(p.model[number_of_states-1-i]);
		
		if (change_strand == 1) {
		    
		    Hmm_State copy(model[i]);
		    
		    // b) permute block of x-indices <-> block of y-indices 
		    //    emission indices: (y_3,y_2,y_1,x_3,x_1,x_1) => (x_3,x_2,x_1,y_3,y_2,y_1)
		    
		    copy.permute_x_and_y_indices_of_emission_probs(x_or_y_indices_first); 
		    
		    // c) complement state
		    //    emission indices: (x_3,x_2,x_1,y_3,y_2,y_1) => (x_3',x_2',x_1',y_3',y_2',y_1')
		    
		    model[i].get_complemented_copy_of(copy);
		}
		// reset all prob and score arrays
		// except the emission prob array which gets set in function get_mirrored_copy_of(state)
		
		model[i].reset_transition_probs();
		model[i].reset_transition_scores();
		model[i].reset_emission_scores();
	    }
	}

	// don't copy results of any alignment
	
	steps=0;
	score=0;
	sequence_of_scores=NULL;
	sequence_of_states=NULL;
	sequence_of_xsteps=NULL;      
	sequence_of_ysteps=NULL;
	condensed_steps=0;
	condensed_sequence_of_scores=NULL;
	condensed_sequence_of_states=NULL;    
	condensed_sequence_of_xsteps=NULL;
	condensed_sequence_of_ysteps=NULL;      
	
	// set transition probs
	
	Prob transition_prob=0;
	
	for (int i=0; i<number_of_states; i++)
	{	  

	    for (int j=0; j<number_of_states; j++)
	    {
		// transition prob j -> i of model becomes transition prob number_of_states-1-i ->
		// number_of_states-1-j of mirrored model
		
		transition_prob=0;

		if (model[i].is_state_next_state(j) == 1)
		{
		    transition_prob=p.model[number_of_states-1-j].get_transition_prob(number_of_states-1-i);
		    model[i].set_transition_prob(j, transition_prob);

		}
	    }	      
	}
    }
    else if (mirror == 0) // if a unmirrored copy of p is requested
    {
	*this=p;
	
	if (change_strand == 1) {
	    
	    for (int i=0; i<number_of_states; i++) {
		
		// complement state
		// emission indices: (x_1,x_2,x_3,y_1,y_2,y_3) => (x_1',x_2',x_3',y_1',y_2',y_3')
		
		model[i].get_complemented_copy_of(p.model[i]);
		
		// reset all prob and score arrays
		// except the emission prob array which gets set in function get_mirrored_copy_of(state)
		// and the transition prob array
		
		model[i].reset_transition_scores();
		model[i].reset_emission_scores();
	    }
	    
	    // don't copy results of any alignment
	    
	    steps=0;
	    score=0;
	    sequence_of_scores=NULL;
	    sequence_of_states=NULL;
	    sequence_of_xsteps=NULL;      
	    sequence_of_ysteps=NULL;
	    condensed_steps=0;
	    condensed_sequence_of_scores=NULL;
	    condensed_sequence_of_states=NULL;    
	    condensed_sequence_of_xsteps=NULL;
	    condensed_sequence_of_ysteps=NULL;      
	    
	} // if (change_strand == 1)
    }
    
}

Hmm::Hmm(const Hmm &p_1, const int number_of_merge_state_1,
	 const Hmm &p_2, const int number_of_merge_state_2,
	 const Prob merge_prob) 
{
  // initialise this Hmm

    model=NULL;
    mirrored=0;
    alphabet=0;
    number_of_states=0;
    steps=0;
    score=0;
    sequence_of_scores=NULL;
    sequence_of_states=NULL;
    sequence_of_xsteps=NULL;
    sequence_of_ysteps=NULL;
    
    condensed_steps=0;
    condensed_sequence_of_scores=NULL;
    condensed_sequence_of_states=NULL;    
    condensed_sequence_of_xsteps=NULL;
    condensed_sequence_of_ysteps=NULL;
    
    // make checks

    int check = 0;
    
    if (p_1.get_alphabet() != p_2.get_alphabet()) {
	cout << "ERROR: Hmm::Hmm: alphabet of Hmm p_1 (" << p_1.get_alphabet()
	     << ") != alphabet of Hmm p_2 (" << p_2.get_alphabet() << ").\n" << flush;
	check++;
    }
    if ((number_of_merge_state_1 < 0) || (number_of_merge_state_1 > (p_1.get_number_of_states()-1))) {
	cout << "ERROR: Hmm::Hmm: number_of_merge_state_1 (" << number_of_merge_state_1
	     << ") out of range [0, " << (p_1.get_number_of_states()-1) << "].\n" << flush;
	check++;
  }
    if ((number_of_merge_state_2 < 0) || (number_of_merge_state_2 > (p_2.get_number_of_states()-1))) {
	cout << "ERROR: Hmm::Hmm: number_of_merge_state_2 (" << number_of_merge_state_2
	 << ") out of range [0, " << (p_2.get_number_of_states()-1) << "].\n" << flush;
	check++;
    }
    if ((merge_prob < 0.0) || (merge_prob > 1.0)) {
	cout << "ERROR: Hmm::Hmm: merge_prob (" << merge_prob << ") out of range [0,1].\n";
	check++;
    }


    if (check == 0) {
	
	int i = 0;
	int k = 0;
	
	const int  jump = 0; 
	
	// NOTE: - jump = 0 means that the two merge states will be connected, but none of them discarded
	//                  from the list of states
	//       - jump = 1 means that the two merge states will be merged into one state (number_of_merge_state_1)
	//                  and that the other merge state (number_of_merge_state_2) will be discarded
	//                  This option has not been completely implemented into this function. What is still
	//                  missing is the updating of the number_of_merge_state_1 state (see ***** below).
	//       - make jump an input parameter once also the jump = 1 option has been implemented

	const int n_of_states_1 = p_1.get_number_of_states();
	const int n_of_states_2 = p_2.get_number_of_states();    
	
	int var_n_of_states = n_of_states_1 + n_of_states_2 - 2;  // do not double count start and end state 
	int var_shift       = var_n_of_states - n_of_states_2;
    
	if (jump == 1) {
	    var_n_of_states -= 1; // state number_of_merge_state_2 will be removed
	    var_shift       += 1;
	}
	
	const int n_of_states   = var_n_of_states;
	const int shift         = var_shift;
	
	// ======================================================================
	// start : initialse Hmm private variables 
	// ======================================================================
	
	model                        = NULL;
	mirrored                     = max(p_1.mirrored, p_2.mirrored);
	alphabet                     = p_1.alphabet;
	number_of_states             = n_of_states;
	steps                        = 0;
	score                        = 0;
	sequence_of_scores           = NULL;
	sequence_of_states           = NULL;
	sequence_of_xsteps           = NULL;
	sequence_of_ysteps           = NULL;    
	condensed_steps              = 0;
	condensed_sequence_of_scores = NULL;
	condensed_sequence_of_states = NULL;    
	condensed_sequence_of_xsteps = NULL;
	condensed_sequence_of_ysteps = NULL;
    
	model = new Hmm_State[n_of_states];

	for (i=0; i<n_of_states_1; i++) { // loop over all p_1 states
	    
	    if (check == 0) {
		
		Hmm_State new_state = p_1.model[i];
		check += new_state.modify_state(n_of_states,            // new_n_of_states
						0,                      // shift
						number_of_merge_state_1,
						number_of_merge_state_1);
	    
		if (check == 0) {
		    if (i == (n_of_states_1 - 1)) {
			
			model[n_of_states - 1] = new_state;  // end state of p_1 becomes end state of new pairhmm
		    
		    }
		    else {
			
			model[i] = new_state;
		   
		    }
		}
		else {
		    cout << "ERROR: Hmm::Hmm: in function Hmm_State::modify_state for p_1 state " 
			 << i << ".\n" << flush;
		    check++;
		    break;
		}
	    }
	}
	
	if (check == 0) {
	    for (i=1; i<n_of_states_2-1; i++) { // loop over all p_2 states except start and end state and merge_state_2
		
		if (check == 0) {
		    
		    Hmm_State new_state = p_2.model[i];
		    
		    if ((i != number_of_merge_state_2) && (jump == 1)) {
			
			check += new_state.modify_state(n_of_states,            // new_n_of_states
							shift, 
							number_of_merge_state_2,
							number_of_merge_state_1);
		    }
		    else if (jump == 0) {
			
			check += new_state.modify_state(n_of_states,            // new_n_of_states
							shift, 
							number_of_merge_state_2,
							number_of_merge_state_2);
		    }
		    
		    if (check == 0) {
			if ((i != number_of_merge_state_2) && (jump == 1)) {
			    if (i < number_of_merge_state_2)  {
				model[i + shift] = new_state;	      
			    }
			    else {
				model[i + shift - 1] = new_state;	      

			    }
			}
			else if (jump == 0) {
			    model[i + shift] = new_state;	      

			}
		    }
		    else {
			cout << "ERROR: Hmm::Hmm: in function Hmm_State::modify_state for p_2 state " 
			     << i << ".\n" << flush;
			check++;
			break;
		    }
		}
	    }
	}
	
	// ======================================================================
	// end : initialse Hmm private variables 
	// ======================================================================

	// ======================================================================
	// start : modify start state
	// ======================================================================
	//
	// - start state never has previous states
	// - start state never has child_state_x and child_state_y as it is a silent state
	// - start state never has child states for transition probs
	// - start state has only next transitions and some of them may be special, so 
	//   we have to update the following :
	//
	//     int number_of_next_states;   
	//     array<int>  numbers_of_next_states; 
	//     array<int>  special_flags_of_transitions_to_next_states; 
	//     array<Prob> transition_probs;      
	
	if (check == 0) {
	    


	    Prob sum_new_probs = 0.0;  // sum of all transition probs which get added from p_2 (may be < 1)
	    Prob sum_all_probs = 0.0;  // sum of all transition probs in model[0] state
	    
	    int  n_of_new_transitions         = 0;
	    
	    // loop over all states that have been added from p_2 and get sum of new transition probs
	    // that shall be added to model[0] state
	    
	    
	    for (i=n_of_states_1 - 1; i<n_of_states - 1; i++) {
		
		n_of_new_transitions += model[i].is_state_previous_state(0);
		
		// i - number of state in new pairhmm
		// k - number of state in p_2 pairhmm
		
		if ( ((i <= (n_of_states_1 - 1 + number_of_merge_state_2 - 2)) && (jump == 1)) ||
		     (jump == 0) ) {
		    
		    k = i - n_of_states_2 + 2;
		}
		else {
		    k = i - n_of_states_2 + 3;
		}
		
		sum_new_probs += p_2.model[0].transition_probs.GetElement(k);
		
	    }
	    
	    
	    if ((n_of_new_transitions > 0) && (sum_new_probs == 0.0)) {
		cout << "ERROR: Hmm::Hmm:: sum_new_probs = " << sum_new_probs 
		     << " == 0, but n_of_new_transitions (" << n_of_new_transitions << ") > 0.\n" << flush;
		check++;
	    }
	    
	    if ((check == 0) && (n_of_new_transitions > 0)) {
		
		// * update mirrored variable
		// --------------------------------------------------
		
		model[0].mirrored = max(p_1.model[0].mirrored, p_2.model[0].mirrored);
	  
		// * update number_of_next_states 
		// --------------------------------------------------
		
		model[0].number_of_next_states                       += n_of_new_transitions;
		//model[0].number_of_special_transitions_to_next_states = 0;
		
		// update :
		// --------------------------------------------------
		// 
		//  int number_of_next_states;   
		//  array<int>  numbers_of_next_states; 
		//  array<int>  special_flags_of_transitions_to_next_states; 
		//  array<Prob> transition_probs;      
		//
		// 1.) if p_1 and p_2 are unmirrored :
		//   
		//  normalise trans_probs such that their sum = 1 using sum_all_probs
		//
		// 2.) else
		//
		//  do not normalise trans_probs i.e. use sum_all_probs := 1
		
		Prob check_sum_probs = 0.0;
	
		array<int>  old_numbers_of_next_states    = model[0].numbers_of_next_states;
		array<Prob> old_transition_probs          = model[0].transition_probs;
		
		array<int>  new_numbers_of_next_states(1);
		array<Prob> new_transition_probs(1);
		
		new_numbers_of_next_states.SetDimension(0, model[0].number_of_next_states);
		new_transition_probs.SetDimension(0, n_of_states);
		
		sum_all_probs = sum_new_probs; // overall normalisation factor
		
		const int old_number_of_next_states = old_numbers_of_next_states.GetDimension(0);
	  
		for (i=0; i<old_number_of_next_states; i++) {
		    
		    sum_all_probs += old_transition_probs.GetElement(i);
		}
	  
		if (! ((p_1.mirrored == 0) && (p_2.mirrored == 0))) {
		    
		    sum_all_probs = 1.0;
		    
		    // do not normalise transition probs if one of the pairhmms or both 
		    // are mirrored
		}
	       
		
		// copy old next state numbers and special flags
		
		for (i=0; i<old_number_of_next_states; i++) {
		    new_numbers_of_next_states.SetElement(i, old_numbers_of_next_states.GetElement(i));
		}
		
		// copy old transition probs
		
		for (i=0; i<n_of_states_1; i++) {
		    
		    // map transition prob of start to p_1.end_state to that of start to new end_state
		    if   (i == (n_of_states_1 - 1)) {k = n_of_states - 1;}
		    else                            {k = i;}
		    new_transition_probs.SetElement(k, old_transition_probs.GetElement(i)/ sum_all_probs);
		    check_sum_probs += new_transition_probs.GetElement(k);
		}
		
		// copy new state numbers, special flags and transition probs
		
		int count = 0;
		
		for (i=n_of_states_1 - 1; i<n_of_states - 1; i++) {
		    if (check == 0) {
			if (model[i].is_state_previous_state(0) == 1) {
			    
			    // i - number of state in new pairhmm
			    // k - number of state in p_2 pairhmm
			    
			    if ( ((i <= (n_of_states_1 - 1 + number_of_merge_state_2 - 2)) && (jump == 1)) ||
				 (jump == 0) ) {
				
				k = i - n_of_states_2 + 2;
			    }
			    else {
				k = i - n_of_states_2 + 3;
			    }
			    
			    new_numbers_of_next_states.SetElement(old_number_of_next_states + count, i);
			    new_transition_probs.SetElement(i, p_2.model[0].transition_probs.GetElement(k)/sum_all_probs);

			    check_sum_probs += new_transition_probs.GetElement(i);
			    count++;
			}
		    }
		}
		
		if ((p_1.mirrored == 0) && (p_2.mirrored == 0) && 
		    (abs(1.0 - check_sum_probs) > Max_deviation)) {
		    cout << "ERROR: Hmm::Hmm: abs(check_sum_probs (" << check_sum_probs 
			 << ") - 1.0) > Max_deviation = " << Max_deviation << "\n" << flush;
		    check++;
		}

		// * update arrays 
		// --------------------------------------------------
		
		model[0].numbers_of_next_states                        = new_numbers_of_next_states;
		model[0].transition_probs                              = new_transition_probs;
	    }
	}

	// ======================================================================
	// end : start state
	// ======================================================================
	
	// ======================================================================
	// Start state : modify end state
	// ======================================================================
	//
	// - end state never has next states
	// - end state never has child_state_x and child_state_y as it is a silent state
	// - end state never has child states for transition probs
	// - end state has array of transition_probs, but they are all 0
	// - end state has only previous transitions and some of them may be special, so 
	//   we have to update the following :
	//
	//     int number_of_previous_states;   
	//     array<int>  numbers_of_previous_states; 
	//     array<int>  special_flags_of_transitions_to_previous_states; 
	
	if (check == 0) {
	    

	    int  n_of_new_transitions         = 0;
	    
	    // loop over all states that have been added from p_2 and get sum of new transition probs
	    // that shall be added to model[n_of_states-1] state
	   	    
	    for (i=n_of_states_1 - 1; i<n_of_states - 1; i++) {
		
		// i - number of state in new pairhmm
		// k - number of state in p_2 pairhmm
		
		if ( ((i <= (n_of_states_1 - 1 + number_of_merge_state_2 - 2)) && (jump == 1)) ||
		     (jump == 0) ) {
		    
		    k = i - n_of_states_2 + 2;
		}
		else {
		    k = i - n_of_states_2 + 3;
		}
	  
		n_of_new_transitions += p_2.model[n_of_states_2-1].is_state_previous_state(k);
		
	    }
	   
	    
	    if ((check == 0) && (n_of_new_transitions > 0)) {
		
		// * update mirrored variable
		// --------------------------------------------------
		
		model[n_of_states-1].mirrored = max(p_1.model[n_of_states_1-1].mirrored, 
						    p_2.model[n_of_states_2-1].mirrored);
		
		// * update number_of_previous_states 
		// --------------------------------------------------
		
		model[n_of_states-1].number_of_previous_states                       += n_of_new_transitions;
		
		// update :
		// --------------------------------------------------
		// 
		//  int number_of_previous_states;   
		//  array<int>  numbers_of_previous_states; 
		//  array<int>  special_flags_of_transitions_to_previous_states; 
		//

		array<int>  old_numbers_of_previous_states = model[n_of_states-1].numbers_of_previous_states;
		
		array<int>  new_numbers_of_previous_states(1);
		
		new_numbers_of_previous_states.SetDimension(0, model[n_of_states-1].number_of_previous_states);
		
		const int old_number_of_previous_states = old_numbers_of_previous_states.GetDimension(0);
		
		// copy old previous state numbers and special flags

		for (i=0; i<old_number_of_previous_states; i++) {
		    new_numbers_of_previous_states.SetElement(i, old_numbers_of_previous_states.GetElement(i));
		    
		}
		
		// copy new state numbers, special flags and transition probs
	
		int count = 0;
		

		for (i=n_of_states_1 - 1; i<n_of_states - 1; i++) {
		    if (check == 0) {
			
			// i - number of state in new pairhmm
			// k - number of state in p_2 pairhmm
			
			if ( ((i <= (n_of_states_1 - 1 + number_of_merge_state_2 - 2)) && (jump == 1)) ||
			     (jump == 0) ) {
			    
			    k = i - n_of_states_2 + 2;
			}
			else {
			    k = i - n_of_states_2 + 3;
			}
			
			if (p_2.model[n_of_states_2-1].is_state_previous_state(k) == 1) {
			    
			    new_numbers_of_previous_states.SetElement(old_number_of_previous_states + count, i);
			    
			    count++;
			}
		    }
		}
		
		// * update arrays 
		// --------------------------------------------------
		
		model[n_of_states-1].numbers_of_previous_states                        = new_numbers_of_previous_states;
		
	    }
	}
	
	// get_number_of_special_transitions_to_previous_states(void) const;
	
	// *****
	
	if ((check == 0) && (jump == 0)) {

	    // ======================================================================
	    // Start : link merge states
	    // ======================================================================
	    //
	    

	    // ----------------------------------------------------------------------
	    // loop over two merge states: 
	    // ----------------------------------------------------------------------

	    // rescale some of the transition probs:
	    //
	    // - keep transition probs to start state and end state fixed
	    // - rescale only non-special transition probs
	    // - rescale the outgoing (if mirrored(number_of_merge_state_1) == 0) or ingoing 
	    //   (if mirrored(number_of_merge_state_1) == 1) non-special transition
	    //   probs such that their sum including the new transition prob to the other merge state 
	    //   remains the same as before
	    // 
	    //   NOTE: if there are no non-special transition probs, report error

	    Prob rescale_factor_1 = 0.0;
	    Prob rescale_factor_2 = 0.0;

	    for (k = 0; k<2; k++) {
		
		int this_merge_state  = 0;
		int other_merge_state = 0;

		if (k==0) {

		    this_merge_state  = number_of_merge_state_1;        
		    other_merge_state = number_of_merge_state_2 + shift;
		}
		else if (k == 1) {
		    
		    this_merge_state  = number_of_merge_state_2 + shift;        
		    other_merge_state = number_of_merge_state_1;
		}

		// implement transition to other_merge_state into this_merge_state
		
		const int mirrored_this_merge_state = model[this_merge_state].get_mirrored();
		
		Prob sum_non_special_probs = 0.0;
		
		if (mirrored_this_merge_state == 0) {
		    
		    // loop over all outgoing non-special transitions except those to the start state, the End state
		    // and the other_merge_state
		    
		    for (i=1; i<n_of_states-1; i++) {
	      
			if ((model[this_merge_state].is_state_next_state(i) == 1)                 &&
			    (i != other_merge_state)) {
			    
			    sum_non_special_probs += model[this_merge_state].get_transition_prob(i);
			    
			}
		    }

		}
		else if (mirrored_this_merge_state == 1) {
		    
		    // loop over all ingoing non-special transitions except those from the Start state, the End state
		    // and the other_merge_state
		    
		    for (i=1; i<n_of_states-1; i++) { 
			
			if ((model[i].is_state_next_state(this_merge_state) == 1) &&  
			    (i != other_merge_state)) {
			    
			    sum_non_special_probs += model[i].get_transition_prob(this_merge_state);

			}
		    }
		}
		
		// ERROR if sum_non_special_probs == 0
		
		if (sum_non_special_probs == 0.0) {
		    if (mirrored_this_merge_state == 0) {
			cout << "ERROR: Hmm::Hmm: sum_non_special_probs of outgoing transition probs of "
			     << "state this_merge_state [" << this_merge_state << "] == 0.\n" << flush;
		    }
		    else if (mirrored_this_merge_state == 1) {
			cout << "ERROR: Hmm::Hmm: sum_non_special_probs of ingoing transition probs of "
			     << "state this_merge_state [" << this_merge_state << "] == 0.\n" << flush;
		    }
		    check++;
		}
		
		if (check == 0) {

		    const Prob rescale_factor = sum_non_special_probs / (sum_non_special_probs + merge_prob);
		    
		    if (k == 0) {
			rescale_factor_1 = rescale_factor;
		    }
		    else if (k == 1) {
			rescale_factor_2 = rescale_factor;
		    }
		    
		    // add transition to state number_of_merge_state_2 
		    
		    if (mirrored_this_merge_state == 0) {
 	       
			// this outgoing transition is included in the following rescaling loop over ingoing
			// transitions and therefore is not rescaled here
			
			model[this_merge_state].set_transition_prob(other_merge_state, 
								    merge_prob);
		    }	
		    else if (mirrored_this_merge_state == 1) {
			
			// this outgoing transition is not included in the following rescaling loop over ingoing 
			// transitions and has to be rescaled separately here
			
			model[this_merge_state].set_transition_prob(other_merge_state, 
								    merge_prob * rescale_factor);
		    }	
		    
		    // rescale the non-special transition probs
		    
		    Prob  check_sum_non_special_probs = 0.0;
		    
		    if (mirrored_this_merge_state == 0) {
			
			// loop over all outgoing non-special transitions except those to the start and the End state	  
			
			for (i=1; i<n_of_states-1; i++) { 
			    
			    if (model[this_merge_state].is_state_next_state(i) == 1) {
				
	
				model[this_merge_state].
				    set_transition_prob(i, 
							model[this_merge_state].get_transition_prob(i) * rescale_factor);
			
				check_sum_non_special_probs += model[this_merge_state].get_transition_prob(i);

			    }
			}

		    } // if mirrored_this_merge_state == 0
		    
		    else if (mirrored_this_merge_state == 1) {
	    
			// loop over all ingoing non-special transitions except those to the Start state and the End state	  
			
			for (i=1; i<n_of_states-1; i++) { 
			    
			    if ((model[i].is_state_next_state(this_merge_state) == 1) &&   
				(i != other_merge_state)) {

				model[i].
				    set_transition_prob(this_merge_state, 
							model[i].get_transition_prob(this_merge_state) * rescale_factor);
				
				check_sum_non_special_probs += model[i].get_transition_prob(this_merge_state);
		    }
			}
			
			check_sum_non_special_probs += model[this_merge_state].get_transition_prob(other_merge_state);
			

		    } // if mirrored_this_merge_state == 1
		    
		    if (abs(sum_non_special_probs - check_sum_non_special_probs) > Max_deviation) {
			
			cout << "ERROR: Hmm::Hmm: | sum_non_special_probs (" << sum_non_special_probs 
			     << ") - check_sum_non_special_probs (" << check_sum_non_special_probs
			     << ") | > Max_deviation = " << Max_deviation << ".\n" << flush;
			check++;
		    }
		    	    
		} // if check == 0
		
	    } // loop over two merge states
	    

	    if (abs(rescale_factor_1 - rescale_factor_2) > Max_deviation) {
		cout << "ERROR: Hmm::Hmm: | rescale_factor_1 (" << rescale_factor_1 
		     << ") - rescale_factor_2 (" << rescale_factor_2 
		     << ") | = " << abs(rescale_factor_1 - rescale_factor_2) 
		     << " > Max_deviation = " << Max_deviation << "\n";
		check++;
	    }
	    
	    // ======================================================================
	    // end : link merge states
	    // ======================================================================
	    //
	    
	} // if check == 0 and jump == 0
	
    } // if check == 0
    
    if (check != 0) {
	
	this->~Hmm();
    }
    
}

Hmm::~Hmm()
{
    mirrored=0;
    alphabet=0;
    number_of_states=0;
    
    steps=0;
    score=0;

    if (sequence_of_scores) delete [] sequence_of_scores;
    sequence_of_scores=NULL;

    if (sequence_of_states) delete [] sequence_of_states;
    sequence_of_states=NULL;
    
    if (sequence_of_xsteps) delete [] sequence_of_xsteps;
    sequence_of_xsteps=NULL;
    
    if (sequence_of_ysteps) delete [] sequence_of_ysteps;
    sequence_of_ysteps=NULL;
    
    condensed_steps=0;
    
    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
    condensed_sequence_of_scores=NULL;
    
    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
    condensed_sequence_of_states=NULL;
    
    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
    condensed_sequence_of_xsteps=NULL;
    
    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
    condensed_sequence_of_ysteps=NULL;
    
}

int Hmm::copy_rectangle_from_strip_to_next_strip(Score**** const strip,
						     const int x_start, const int x_end,
						     const int y_start, const int y_end,
						     const int x_start_copy, const int x_end_copy,
						     const int y_start_copy, const int y_end_copy,
						     Score**** const next_strip,
						     const int x_start_next, const int x_end_next,
						     const int y_start_next, const int y_end_next,
						     const int x_start_copy_next, const int x_end_copy_next,
						     const int y_start_copy_next, const int y_end_copy_next,
						     const int direction)
{
    int check=0;

    // check existence of strip and next_strip
    
    if ((*strip)==NULL)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: strip is NULL.\n" << flush;
	check+=1;
    }

    if ((*next_strip)==NULL)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: next_strip is NULL.\n" << flush;
	check+=1;
    }

    // check correct order of all (Start,End) pairs, i.e. that start <= end
    
    if (x_start>x_end)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_start (" << x_start << ") > x_end ("
	     << x_end << ").\n" << flush;
	check+=1;
    }

    if (y_start>y_end)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_start (" << y_start << ") > y_end ("
	     << y_end << ").\n" << flush;
	check+=1;
    }

    if (x_start_next>x_end_next)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_start_next (" << x_start_next 
	     << ") > x_end_next (" << x_end_next << ").\n" << flush;
	check+=1;
    }

    if (y_start_next>y_end_next)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_start_next (" << y_start_next
	     << ") > y_end_next (" << y_end_next << ").\n" << flush;
	check+=1;
    }

    if (x_start_copy>x_end_copy)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_start_copy (" << x_start_copy 
	  << ") > x_end_copy (" << x_end_copy << ").\n" << flush;
	check+=1;
    }

    if (y_start_copy>y_end_copy)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_start_copy (" << y_start_copy 
	     << ") > y_end_copy (" << y_end_copy << ").\n" << flush;
	check+=1;
    }

    if (x_start_copy_next>x_end_copy_next)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_start_copy_next (" << x_start_copy_next 
	     << ") > x_end_copy_next (" << x_end_copy_next << ").\n" << flush;
	check+=1;
    }

    if (y_start_copy_next>y_end_copy_next)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_start_copy_next (" << y_start_copy_next
	     << ") > y_end_copy_next (" << y_end_copy_next << ").\n" << flush;
	check+=1;
    }
  
    // check that each small rectangle lies inside corresponding larger rectangle
  
    if ((x_start_copy<x_start) || (x_start_copy>x_end))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_start_copy (" << x_start_copy 
	     << ") out of range [x_start (" << x_start << "), x_end (" << x_end << ")].\n" << flush;
	check+=1;
    }

    if ((x_end_copy<x_start) || (x_end_copy>x_end))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_end_copy (" << x_end_copy 
	     << ") out of range [x_start (" << x_start << "), x_end (" << x_end << ")].\n" << flush;
	check+=1;
    }

    if ((y_start_copy<y_start) || (y_start_copy>y_end))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_start_copy (" << y_start_copy 
	     << ") out of range [y_start (" << y_start << "), y_end (" << y_end << ")].\n" << flush;
	check+=1;
    }

    if ((y_end_copy<y_start) || (y_end_copy>y_end))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_end_copy (" << y_end_copy 
	     << ") out of range [y_start (" << y_start << "), y_end (" << y_end << ")].\n" << flush;
	check+=1;
    }

    if ((x_start_copy_next<x_start_next) || (x_start_copy_next>x_end_next))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_start_copy_next (" << x_start_copy_next 
	     << ") out of range [x_start_next (" << x_start_next << "), x_end_next (" << x_end_next << ")].\n" << flush;
	check+=1;
    }
    
    if ((x_end_copy_next<x_start_next) || (x_end_copy_next>x_end_next))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: x_end_copy_next (" << x_end_copy_next 
	     << ") out of range [x_start_next (" << x_start_next << "), x_end_next (" << x_end_next << ")].\n" << flush;
	check+=1;
    }

    if ((y_start_copy_next<y_start_next) || (y_start_copy_next>y_end_next))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_start_copy_next (" << y_start_copy_next 
	     << ") out of range [y_start_next (" << y_start_next << "), y_end_next (" << y_end_next << ")].\n" << flush;
	check+=1;
    }

    if ((y_end_copy_next<y_start_next) || (y_end_copy_next>y_end_next))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: y_end_copy_next (" << y_end_copy_next 
	     << ") out of range [y_start_next (" << y_start_next << "), y_end_next (" << y_end_next << ")].\n" << flush;
	check+=1;
    }

    // check that two small rectangles that are to be copied onto eachother have the same dimensions

    if ((x_end_copy-x_start_copy) != (x_end_copy_next-x_start_copy_next))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: (x_end_copy-x_start_copy) ("
	     << (x_end_copy-x_start_copy) << ") != (x_end_copy_next-x_start_copy_next) ("
	     << (x_end_copy_next-x_start_copy_next) << ").\n" << flush;
	check+=1;
    }
    
    if ((y_end_copy-y_start_copy) != (y_end_copy_next-y_start_copy_next))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: (y_end_copy-y_start_copy) ("
	     << (y_end_copy-y_start_copy) << ") != (y_end_copy_next-y_start_copy_next) ("
	     << (y_end_copy_next-y_start_copy_next) << ").\n" << flush;
	check+=1;
    }
    
    // check values of direction

    if ((direction!=1) && (direction!=-1))
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: direction (" << direction 
	     << ") must be either 1 or -1.\n" << flush;
	check+=1;
    }

    // check this Hmm

    if (number_of_states<1)
    {
	cout << "ERROR Hmm::copy_rectangle_from_strip_to_next_strip: number_of_states (" 
	     << number_of_states << ") <1.\n" << flush;
	check+=1;}
    
    if (check==0)
    {
	// determine coordinates of both rectangles within their respective strips
	
	int x_start_copy_strip = 0;
	int x_end_copy_strip   = 0;
	
	int y_start_copy_strip = y_start_copy-y_start; 
	int y_end_copy_strip   = y_end_copy-y_start;   
	
	int x_start_copy_next_strip = 0;
	int x_end_copy_next_strip   = 0;
	
	int y_start_copy_next_strip = y_start_copy_next-y_start_next; 

	if (direction==1)
	{
	    x_start_copy_strip      = x_start_copy-x_start;
	    x_end_copy_strip        = x_end_copy-x_start;  
	    x_start_copy_next_strip = x_start_copy_next-x_start_next;  
	    x_end_copy_next_strip   = x_end_copy_next-x_start_next;    

	}
	else if (direction==-1)
	{
	    x_start_copy_strip      = x_end-x_start_copy;
	    x_end_copy_strip        = x_end-x_end_copy;  
	    x_start_copy_next_strip = x_end_next-x_start_copy_next;  
	    x_end_copy_next_strip   = x_end_next-x_end_copy_next;    

	}
	
	// copy things over from rectangle in strip to different rectangle in next_strip

	int x_offset=x_start_copy_next_strip-x_start_copy_strip;
	int y_offset=y_start_copy_next_strip-y_start_copy_strip;

	for (int state=0; state<number_of_states; state++)
	{
	    for (int x=x_start_copy_strip; x<x_end_copy_strip+1; x++)
	    {

		for (int y=y_start_copy_strip; y<y_end_copy_strip+1; y++)
		{
		    (*next_strip)[state][x+x_offset][y+y_offset]=(*strip)[state][x][y];

		}
	    }
	}
	
	if (direction==-1)
	{
	    // transfer next_strip if direction==-1
	    //
	    //       state -> number_of_states-1-state
	    //       x     -> (length_of_next_strip_in_x_direction-1)-x

	  int x_length_next=x_end_next-x_start_next+1;
	  int y_length_next=y_end_next-y_start_next+1;
	  
	    Score*** copy_strip= new Score**[number_of_states];
      
	    for (int state=0; state<number_of_states; state++)
	    {
		copy_strip[state]= new Score*[x_length_next];
		
		for (int i=0; i<x_length_next; i++)
		{
		    copy_strip[state][i]= new Score[y_length_next];
		    
		    for (int j=0; j<y_length_next; j++)
		    {	     
			copy_strip[state][i][j]= (*next_strip)[state][i][j];
		    }
		}	 
	    }

	    {	  
		for (int state=0; state<number_of_states; state++)
		{
		    for (int i=0; i<x_length_next; i++)
		    {
			for (int j=0; j<y_length_next; j++)
			{	     
			    (*next_strip)[number_of_states-1-state][x_length_next-1-i][j]=copy_strip[state][i][j];
			}
		    }	 
		}
	    }

	    if (copy_strip)
	    {
		for (int state=0; state<number_of_states; state++)
		{
		    for (int i=0; i<x_length_next; i++)
		    {
			if (copy_strip[state][i]) delete [] copy_strip[state][i];
			copy_strip[state][i] = NULL;
		    }
		    if (copy_strip[state]) delete [] copy_strip[state];
		    copy_strip[state] = NULL;
		}
		if (copy_strip) delete [] copy_strip;      
		copy_strip=NULL;
	    }
	}
    }
    return(check);
}

int Hmm::allocate_memory_for_strip(const Hmm *s, 
				       const int x_margin,
				       const int y_start, const int y_end,
				       Score**** const strip)
{
    int check=0;
    
    // check that strip is NULL
    
    if ((*strip) != NULL) {
	cout << "ERROR Hmm::allocate_memory_for_strip: strip has to be NULL.\n" << flush;
	check+=1;
    }
    
    // check x_margin value
    
    if (x_margin<1)
    {
	cout << "ERROR Hmm::allocate_memory_for_strip: x_margin (" << x_margin 
	     << ") out of range. Must be > 0.\n" << flush;
	check+=1;
    }
    
    // check y_start, y_end values
    
    if (y_start > y_end)
    {
	cout << "ERROR Hmm::allocate_memory_for_strip: y_start (" << y_start
	     << ") > y_end (" << y_end << ")\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::allocate_memory_for_strip: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }

    // check existence of Hmm

    if (s==NULL)
    {
	cout << "ERROR Hmm::allocate_memory_for_strip: Hmm s (" << s
	     << ") NULL\n" << flush;
	check+=1;
    }

    if (check==0)
    {
	const int n_of_states=s->get_number_of_states();
	
	(*strip)= new Score**[n_of_states]; 
	
	for (int state=0; state<n_of_states; state++)
	{
	    (*strip)[state]= new Score*[x_margin];
	    
	    for (int i=0; i<x_margin; i++)
  	    {
		(*strip)[state][i]= new Score[y_end-y_start+1]; 
	      
		for (int j=0; j<y_end-y_start+1; j++)
		{	     
		    (*strip)[state][i][j]= Logzero;
		}
	    }	 
	}
    }
    return(check);
}

int Hmm::check_that_there_are_valid_values_in_strip(Score**** const strip,
							const Hmm *s, 
							const int x_margin,
							const int y_start, const int y_end)
							
{
    int check=0;
    
    // check existence of Hmm
    
    if (s==NULL)
    {
	cout << "ERROR Hmm::check_that_there_are_valid_values_in_strip: Hmm s is NULL.\n" << flush;
	check+=1;
    }
    
    // check existence of strip
    
    if ((*strip) == NULL)
    {
	cout << "ERROR Hmm::check_that_there_are_valid_values_in_strip: strip is NULL.\n" << flush;
	check+=1;
    }
    
    // check x_margin value

    if (x_margin<1)
    {
	cout << "ERROR Hmm::check_that_there_are_valid_values_in_strip: x_margin (" << x_margin 
	     << ") out of range. Must be > 0.\n" << flush;
	check+=1;
    }
  
    // check y_start, y_end values

    if (y_start > y_end)
    {
	cout << "ERROR Hmm::check_that_there_are_valid_values_in_strip: y_start (" << y_start
	     << ") > y_end (" << y_end << ")\n" << flush;
	check+=1;
    }

    if (y_start<0)
    {
	cout << "ERROR Hmm::check_that_there_are_valid_values_in_strip: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }


    if (check==0)
    {
	int n_of_states=0;
	n_of_states=s->get_number_of_states();
      
	check=1;

	for (int state=0; state<n_of_states; state++)
	{
	    for (int i=0; i<x_margin; i++)
	    {
		for (int j=0; j<y_end-y_start+1; j++)
		{	     
		    if ((*strip)[state][i][j]>Logzero)  // if there is one valid element exit and declare strip to be o.k.
		    {
			check=0;
			break;
		    }
		}
		if (check==0) {break;}
	    }	 
	    if (check==0) {break;}
	}

	if (check!=0)
	{
	    cout << "ERROR Hmm::check_that_there_are_valid_values_in_strip: all elements in strip are Logzero.\n" << flush;
	}
    }
    return(check);
}

int Hmm::delete_memory_for_strip(const Hmm *s, 
				     const int x_margin,
				     Score**** const strip)
{
    int check=0;
    
    if ((*strip) != NULL) {
	
	// check x_margin value
	
	if (x_margin<1) {
	    cout << "ERROR Hmm::delete_memory_for_strip: x_margin (" << x_margin 
		 << ") out of range. Must be > 0.\n" << flush;
	    check+=1;
	}
  
	// check existence of Hmm

	if (s==NULL) {
	    cout << "ERROR Hmm::delete_memory_for_strip: Hmm s (" << s
		 << ") NULL\n" << flush;
	    check+=1;
	}
    
	if (check==0) {

	    const int n_of_states=s->get_number_of_states();

	    for (int state=0; state<n_of_states; state++) {
		for (int i=0; i<x_margin; i++) {
		    if ((*strip)[state][i]) delete [] (*strip)[state][i];
		    (*strip)[state][i] = NULL;
		}	 
		if ((*strip)[state]) delete [] (*strip)[state];
		(*strip)[state] = NULL;
	    }
	    if (*strip) delete [] (*strip);
	    (*strip)=NULL;
	}
    }
    return(check);
}

int Hmm::allocate_memory_for_viterbi_rectangle(const Hmm *s, 
						   const int x_start, const int x_end,
						   const int y_start, const int y_end,
						   const int direction,
						   Score**** const viterbi_rectangle)
{
    int check=0;

    if ((*viterbi_rectangle) != NULL) {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: viterbi_rectangle has to be NULL.\n" << flush;
	check+=1;
    }

    // check x_start, x_end and y_start, y_end values

    if (x_start > x_end)
    {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: x_start (" << x_start 
	     << ") > x_end (" << x_end << ")\n" << flush;
	check+=1;
    }
    if (x_start<0)
    {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: x_start (" << x_start 
	     << ") <0\n" << flush;
	check+=1;
    }
    if (y_start > y_end)
    {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: y_start (" << y_start
	     << ") > y_end (" << y_end << ")\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    // check existence of Hmm

    if (s==NULL)
    {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: Hmm s (" << s
	     << ") NULL\n" << flush;
	check+=1;
    }

    // check direction
    
    if ( ! ((direction==1) || (direction==-1)))
    {
	cout << "ERROR Hmm::allocate_memory_for_viterbi_rectangle: direction (" << direction
	     << ") != 1 or -1.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {
	const int n_of_states=s->get_number_of_states();
	int x_length=x_end-x_start+2;
	int y_length=y_end-y_start+2;
	
	if (direction==-1) 
	{
	    if (x_start != 0)
	    {
		x_length=x_end-x_start+3;
	    }
	    if (y_start != 0)
	    {
		y_length=y_end-y_start+3;
	    }
	}

	// allocate memory for the viterbi_rectangles of all states      
	
	(*viterbi_rectangle)= new Score**[n_of_states];
	
	for (int state=0; state<n_of_states; state++)
	{
	    (*viterbi_rectangle)[state]= new Score*[x_length]; 
	    
	    for (int i=0; i<x_length; i++)
	    {
		(*viterbi_rectangle)[state][i]= new Score[y_length];
		
		for (int j=0; j<y_length; j++)
		{	     
		    (*viterbi_rectangle)[state][i][j]= Logzero;
		}
	    }	
	}

    }
    return(check);
}

int Hmm::allocate_memory_for_state_path(const int x_start, const int x_end,
					    const int y_start, const int y_end,
					    Score** const local_sequence_of_scores,
					    int** const local_sequence_of_states,
					    int** const local_sequence_of_xsteps,
					    int** const local_sequence_of_ysteps)
{
    int check=0;
    
    if ((*local_sequence_of_scores) != NULL) {
	cout << "ERROR Hmm::allocate_memory_for_state_path: local_sequence_of_scores has to be NULL.\n" << flush;
	check+=1;
    }
    if ((*local_sequence_of_states) != NULL) {
	cout << "ERROR Hmm::allocate_memory_for_state_path: local_sequence_of_states has to be NULL.\n" << flush;
	check+=1;
    }
    if ((*local_sequence_of_xsteps) != NULL) {
	cout << "ERROR Hmm::allocate_memory_for_state_path: local_sequence_of_xsteps has to be NULL.\n" << flush;
	check+=1;
    }
    if ((*local_sequence_of_ysteps) != NULL) {
	cout << "ERROR Hmm::allocate_memory_for_state_path: local_sequence_of_ysteps has to be NULL.\n" << flush;
	check+=1;
    }
  
    // check x_start, x_end and y_start, y_end values
    
    if (x_start > x_end)
    {
	cout << "ERROR Hmm::allocate_memory_for_state_path: x_start (" << x_start 
	     << ") > x_end (" << x_end << ")\n" << flush;
	check+=1;
    }
    if (x_start<0)
    {
	cout << "ERROR Hmm::allocate_memory_for_state_path: x_start (" << x_start 
	     << ") <0\n" << flush;
	check+=1;
    }
    if (y_start > y_end)
    {
	cout << "ERROR Hmm::allocate_memory_for_state_path: y_start (" << y_start
	     << ") > y_end (" << y_end << ")\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::allocate_memory_for_state_path: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

	// allocate memory for state path
	
	(*local_sequence_of_scores)= new Score[(x_end-x_start)+(y_end-y_start)];
	(*local_sequence_of_states)= new int[(x_end-x_start)+(y_end-y_start)];
	(*local_sequence_of_xsteps)= new int[(x_end-x_start)+(y_end-y_start)];
	(*local_sequence_of_ysteps)= new int[(x_end-x_start)+(y_end-y_start)];  
	
	for (int i=0; i < ((x_end-x_start)+(y_end-y_start)); i++)
	{
	    (*local_sequence_of_scores)[i]=Logzero;
	    (*local_sequence_of_states)[i]=0;
	    (*local_sequence_of_xsteps)[i]=0;
	    (*local_sequence_of_ysteps)[i]=0;
	}
    }
    return(check);
}

int Hmm::allocate_memory_for_state_path(const int max_length_of_state_path,
					    const int default_steps,
					    const Score default_overall_score,
					    const int default_int,
					    const Score default_score) { 

    int check = 0;
    
    if (max_length_of_state_path < 1) {
	
	cout << "ERROR Hmm::allocate_memory_for_state_path: max_length_of_state_path ("
	     << max_length_of_state_path << ") < 1. Don't do anything.\n";
	check+=1;
    }
    if (check == 0) {
	
	// first delete state path, if it exists
	
	this->reset_variables_for_state_path();
	
	// now allocate new memory and assign desired values
	
	score = default_overall_score;
	steps = default_steps;
	
	sequence_of_scores = new Score[max_length_of_state_path];
	sequence_of_states = new int[max_length_of_state_path];
	sequence_of_xsteps = new int[max_length_of_state_path];
	sequence_of_ysteps = new int[max_length_of_state_path];
	
	int i;
	
	for (i=0; i< max_length_of_state_path; i++) {
	    
	    sequence_of_states[i] = default_int;   
	    sequence_of_scores[i] = default_score;    
	    sequence_of_xsteps[i] = default_int;   
	    sequence_of_ysteps[i] = default_int;   
	}
    }
    return(check);
}

int Hmm::delete_memory_for_state_path(Score** const local_sequence_of_scores,
					  int** const local_sequence_of_states,
					  int** const local_sequence_of_xsteps,
					  int** const local_sequence_of_ysteps)
{


    // deleting local_sequence_of_scores and local_sequence_of_states
    
    if ((*local_sequence_of_scores)!=NULL) delete [] (*local_sequence_of_scores);
    (*local_sequence_of_scores)=NULL;
    
    if ((*local_sequence_of_states)!=NULL) delete [] (*local_sequence_of_states);
    (*local_sequence_of_states)=NULL;
    
    if ((*local_sequence_of_xsteps)!=NULL) delete [] (*local_sequence_of_xsteps);
    (*local_sequence_of_xsteps)=NULL;
    
    if ((*local_sequence_of_ysteps)!=NULL) delete [] (*local_sequence_of_ysteps);
    (*local_sequence_of_ysteps)=NULL;
    
    return(0);
}

int Hmm::delete_memory_for_viterbi_rectangle(const Hmm *s, 
						 const int x_start, const int x_end,
						 const int direction,
						 Score**** const viterbi_rectangle)
{
    int check=0;
    
    if ((*viterbi_rectangle)!=NULL) {
	
	const int n_of_states=s->get_number_of_states();
	
	// check x_start, x_end values
	
	if (x_start > x_end) {
	    cout << "ERROR Hmm::delete_memory_for_viterbi_rectangle: x_start (" << x_start 
		 << ") > x_end (" << x_end << ")\n" << flush;
	    check+=1;
	}
	if (x_start<0) {
	    cout << "ERROR Hmm::delete_memory_for_viterbi_rectangle: x_start (" << x_start 
		 << ") <0\n" << flush;
	    check+=1;
	}
	// check existence of Hmm
	
	if (s==NULL) {
	    cout << "ERROR Hmm::delete_memory_for_viterbi_rectangle: Hmm s (" << s
		 << ") NULL\n" << flush;
	    check+=1;
	}
    
	// check direction
	
	if ( ! ((direction==1) || (direction==-1))) {
	    cout << "ERROR Hmm::delete_memory_for_viterbi_rectangle: direction (" << direction
		 << ") != 1 or -1.\n" << flush;
	    check+=1;
	}
	
	if (check==0) {

	    int x_length=x_end-x_start+2;
	    
	    if (direction==-1) {
		if (x_start != 0) {
		    x_length=x_end-x_start+3;
		}
	    }

	    for (int state=0; state<n_of_states; state++) {
		for (int i=0; i<x_length; i++) {
		    if ((*viterbi_rectangle)[state][i]) delete [] (*viterbi_rectangle)[state][i];
		    (*viterbi_rectangle)[state][i] = NULL;
		}	 
		if ((*viterbi_rectangle)[state]) delete [] (*viterbi_rectangle)[state];
		(*viterbi_rectangle)[state] = NULL;
	    }
	    if ((*viterbi_rectangle)) delete [] (*viterbi_rectangle);
	    (*viterbi_rectangle)=NULL;
	}
    }
    return(check);
}

int Hmm::calculate_viterbi_rectangle(Score*** const viterbi_rectangle,
					 const Hmm *s, 
					 const Hmm *s_unmirrored,
					 const Sequence *x, const int x_start, const int x_end,
					 const Sequence *y, const int y_start, const int y_end,
					 const int* start_states, 
					 const int number_of_start_states, 
					 int x_margin, 
					 int y_margin, 
					 const int direction)
{
    int check=0;
    int n_of_states=0;
    n_of_states=s->get_number_of_states();

    // check viterbi_rectangle
    
    if (viterbi_rectangle==NULL)
    {
	cout << "ERROR Hmm::calculate_viterbi_rectangle: viterbi_rectangle (" << viterbi_rectangle
	     << ") = NULL\n" << flush;
	check+=1;
    }
    
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::calculate_viterbi_rectangle : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL \n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::calculate_viterbi_rectangle : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL \n" << flush;
	check+=1;
    }
  
    // Check Hmm
    
    if (s->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : number of states = " 
	     << s->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (s->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : state 0 of type != Start \n" << flush;
	cout << "s->model[0].get_letters_to_read() : "<<s->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (s->model[n_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : last state of type != End \n" << flush;
	cout << "s->model["<<n_of_states-1<<"].get_letters_to_read() : "<<s->model[n_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (int i=0; i<n_of_states; i++)
    {
	if (s->model[i].get_alphabet()!=alphabet) 
	{
	    cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << s->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_mirrored()!=s->get_mirrored()) 
	{
	    cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << s->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << s->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_number_of_states() != n_of_states)
	{
	    cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << s->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << n_of_states << "\n" << flush;
	    check+=1;
	}

	if ((i!=0) && (i!=(n_of_states-1)) &&(s->model[i].get_letters_to_read()==0 ))
	{
	    cout << "ERROR: Hmm::calculate_viterbi_rectangle :  Hmm : state i = " 
		 << i << " should emit! \n "<<flush;
	    check+=1;
	}
	
    }
    
    // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::calculate_viterbi_rectangle: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::calculate_viterbi_rectangle: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::calculate_viterbi_rectangle: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::calculate_viterbi_rectangle: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }

    // either x_margin and y_margin are used and viterbi_rectangle contains an initialised
    // rectangle of values or number_of_start_states is used to initialise the calculation
    // of the viterbi_rectangle, i.e. the two options are mutually exclusive
    
    if ( ((x_margin>0) || (y_margin>0)) && (number_of_start_states!=0))
    {
	cout << "ERROR: Hmm::calculate_viterbi_rectangle : x_margin (" << x_margin 
	     << ") > 0 or y_margin (" << y_margin << ") > 0, but number_of_start_states ("
	     << number_of_start_states << ") != 0. Set x_margin=y_margin=0 if you want"
	     << " to use Start states OR set number_of_start_states=0 if you want to Start "
	     << " the calculation of the viterbi_rectangle with an initialised sub-rectangle.\n" << flush;
	check+=1;
    }
    else if ( ((x_margin>0) || (y_margin>0)) && (number_of_start_states==0)) 
	// if sub-rectangle is used for initialisation
    {
	// x_margin has to be 0 for x_start=0
	// y_margin has to be 0 for y_start=0

	if ((x_start==0) && (x_margin!=0))
	{
	    cout << "ERROR: Hmm::calculate_viterbi_rectangle : x_margin (" << x_margin 
		 << ") has to be zero for x_start (" 
		 << x_start << ") = 0.\n" << flush;
	    check+=1;
	}
	
	if ((y_start==0) && (y_margin!=0))
	{
	    cout << "ERROR: Hmm::calculate_viterbi_rectangle : y_margin (" << y_margin 
		 << ") has to be zero for y_start (" 
		 << y_start << ") = 0.\n" << flush;
	    check+=1;
	}
	
	// check that values of x_margin and y_margin are sensible
      
	if (x_margin < 0)
	{
	    cout << "ERROR Hmm::calculate_viterbi_rectangle: x_margin (" << x_margin
		 << ") < 0\n" << flush;
	    check+=1;
	}
	
	if (y_margin < 0)
	{
	    cout << "ERROR Hmm::calculate_viterbi_rectangle: y_margin (" << y_margin
		 << ") < 0\n" << flush;
	    check+=1;
	}
      
	// check that x_margin and y_margin are compatible with values of x_start, x_end and y_start and y_end
      
	if ((x_end-x_start+1) < x_margin)
	{
	    cout << "ERROR Hmm::calculate_viterbi_rectangle: x_margin (" << x_margin
		 << ") > (x_end (" << x_end << ") - x_start (" << x_start << ") + 1) = " 
		 << x_end-x_start+1 << "\n" << flush;
	    check+=1;
	}
      
	if ((y_end-y_start+1) < y_margin)
	{
	    cout << "ERROR Hmm::calculate_viterbi_rectangle: y_margin (" << y_margin
		 << ") > (y_end (" << y_end << ") - y_start (" << y_start << ") + 1) = " 
		 << y_end-y_start+1 << "\n" << flush;
	    check+=1;
	}
    }
    else if ((x_margin==0) && (y_margin==0)) // if start_states shall be used for initialisation
    {
	// check start states
	
	if (number_of_start_states<1)
	{
	    cout << "ERROR Hmm::calculate_viterbi_rectangle: number_of_start_sites (" 
		 << number_of_start_states << ") <1.\n" << flush;
	    check+=1;
	}
	else
	{
	    if ((x_start==0) && (y_start==0) && (direction==1))
	    {
		if ((number_of_start_states!=1) || (start_states[0]!=0))
		{
		    cout << "ERROR Hmm::calculate_viterbi_rectangle: direction==1: for "
			 << "x_start=y_start=0 only "
			 << "one start_state is possible, Start state (i.e. the number_of_start_states ("
			 << number_of_start_states << ") has to be 1 and start_states[0] ("
			 << start_states[0] << ") has to be 0).\n" << flush;
		    check+=1;
		}
	    }
	    if ((x_end==x->length()) && (y_end==y->length()) && (direction==-1))
	    {
		if ((number_of_start_states!=1) || (start_states[0]!=n_of_states-1))
		{
		    cout << "ERROR Hmm::calculate_viterbi_rectangle: direction==-1: for x_end=x->length()"
			 << " and y_end=y->length() only one end_state is possible, the End state "
			 << "(i.e. the number_of_start_states (" << number_of_start_states 
			 << ") has to be 1 and start_states[0] (" << start_states[0] 
			 << ") has to be number_of_states-1 (" << n_of_states-1 << "))\n" << flush;
		    check+=1;
		}
	    }
	    if (!  (((x_start==0) && (y_start==0)) || ((x_end==x->length()) && (y_end==y->length()))))
		// if we don't start at (0,0)
	    {
		// check that each start state lies in the allowed range (1...n_of_states-2)
		// i.e. Start and End states are not allowed 
	      
		for (int i=0; i<number_of_start_states; i++)
		{
		    if ((start_states[i]<1) || (start_states[i]>n_of_states-2))
		    {
			cout << "ERROR Hmm::calculate_viterbi_rectangle: start_states[" << i 
			     << "] (" << start_states[i] << ") out of range [1, " << n_of_states-2
			     << "].\n" << flush;
			check+=1;
		    }
		}
	    }
	}
    }
    
    // check direction
    
    if ((direction!=1) && (direction!=-1))
    {
	cout << "ERROR Hmm::calculate_viterbi_rectangle: direction (" << direction 
	     << ") has to be 1 or -1.\n" << flush;
	check+=1;
    }
    
    // enter program if all check were successfully passed
    
    if (check==0)
    {
	int i = 0;

      // variables which will only be needed to calculate alignment

	int shift_index = bitshift(alphabet);

	int j;
	int xsteps, ysteps;
	int deltax, deltay;
	
	int max_state, max_deltax, max_deltay;
	Score max_score;
	
	int test_state, test_deltax, test_deltay;
	Score test_score;
	Score pre_test_score;
	
	int states=n_of_states;
	
	int start_x=0;
	int end_x=0;
	int start_y=0;
	int end_y=0;
	int offset=0;
	
	if (direction==1)
	{
	    offset=0;
	    
	    start_x=x_start; 
	    end_x=x_end+1;
	    start_y=y_start;
	    end_y=y_end+1; 
	}
	else if (direction==-1) // permute start and end points if running backwards
	{
	    offset=1;

	    start_x=x_end-1;
	    start_y=y_end-1;
	    
	    end_x=x_start-2;
	    end_y=y_start-2;
	}
	
	if ((number_of_start_states>0) && (x_margin==0) && (y_margin==0))
	{
	    // implement start constraints into viterbi_rectangle

	    
	    for (i=0; i<number_of_start_states; i++)
	    {

		if (direction==1)
		{	
		    viterbi_rectangle[start_states[i]][start_x-x_start][start_y-y_start]=0.;

		}
		else if (direction==-1)
		{
		    viterbi_rectangle[n_of_states-1-start_states[i]]
			[start_x+offset-x_start]
			[start_y+offset-y_start]=0.;

		}
	    }
	}    
	for (xsteps=start_x; xsteps != end_x; xsteps+=direction) 
	{	 
	    int max_count = 0;

	  // start of ysteps loop can be a function of xsteps value

	  int start_of_ysteps_loop=0;

	  if (abs(xsteps-start_x) < x_margin)
	  {
	      start_of_ysteps_loop=start_y+direction*y_margin;
	  }
	  else
	  {
	      start_of_ysteps_loop=start_y;
	  }
	  
	  for ( ysteps=start_of_ysteps_loop; ysteps != end_y; ysteps+=direction) 
	  {

	      if ((direction==1) && (xsteps==0) && (ysteps==0))
	      {
		  viterbi_rectangle[0][xsteps-x_start][ysteps-y_start]=0;

	      }
	      else if ((direction==-1) && ((x->length()-1-xsteps)==0) && ((y->length()-1-ysteps)==0))	       
	      {
		  viterbi_rectangle[0][xsteps+offset-x_start][ysteps+offset-y_start]=0.;

	      }
	      else
	      {

		  for (int new_state=1; new_state<states-1; new_state++)
		  {		

		      max_score=Logzero;
		      max_state=0;
		      max_deltax=0;
		      max_deltay=0;
		
		      deltax=0;
		      deltay=0;
		      
		      deltax=s->model[new_state].get_letters_to_read_x();
		      deltay=s->model[new_state].get_letters_to_read_y();

		      if ( ((direction==1) && 
			    (xsteps-start_x>=deltax) && (ysteps-start_y>=deltay)) ||
			   ((direction==-1) && 
			    ((start_x-xsteps)>=deltax) && ((start_y-ysteps)>=deltay)) )
		      {

			  int linear_index_emission = 0;
			  int* indices     = NULL;
			  indices = new int[s->model[new_state].get_number_of_dimensions_of_emission_scores()+1];
			  
			  int xcount=0;
			  int ycount=0;
			
			  if (direction==1)
			  {
			      xcount=deltax;
			      ycount=0;
			  }
			  else if (direction==-1)
			  {
			      xcount=0;
			      ycount=deltay;
			  }
			  for (i=ycount; i<deltax+ycount; i++)
			  {
			      indices[i+1] = x->letter(xsteps - direction*(deltax+ycount-i));

			  }
			  for (j=xcount; j<deltay+xcount; j++)
			  {
			      indices[j+1] = y->letter(ysteps - direction*(deltay+xcount-j)); 

			  }

			  int number_of_old_states = s->model[new_state].get_number_of_previous_states();
			  
			  for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
			  {
			      int old_state = s->model[new_state].get_number_of_previous_state(old_state_number);
			      
			      pre_test_score=Logzero;
			      test_score=Logzero;
			      test_state=0;
			      test_deltax=0;
			      test_deltay=0;
			      
			      int x_position = 0;
			      int y_position = 0;
			      
			      if (direction == 1)
			      {
				  x_position = xsteps - deltax;
				  y_position = ysteps - deltay;
			      }
			      else if (direction == -1)
			      {
				  x_position = xsteps + 1;
				  y_position = ysteps + 1;
			      }

			      indices[0] = old_state;
			      
			      linear_index_emission = indices[1];

			      for (int p=2; p<s->model[new_state].get_number_of_dimensions_of_emission_scores()+1; p++)
			      {
				  linear_index_emission = (linear_index_emission << shift_index) | indices[p];

			      }
			      

			      pre_test_score=s->new_get_transition_score(s_unmirrored, 
									 old_state, new_state, 
									 x, x_position, y, y_position) +
				  s->model[new_state].get_emission_score(x, x_position, y, y_position,
									 linear_index_emission);

			      if (pre_test_score > Logzero)
			      {
				  test_state=old_state;
				  test_score=pre_test_score+
				      viterbi_rectangle[old_state][xsteps-direction*deltax+offset-x_start][ysteps-direction*deltay+offset-y_start];
				  
				  test_deltax=deltax;
				  test_deltay=deltay;

				  if (test_score>max_score)
				  {

				      max_score=test_score;
				      max_state=test_state;
				      max_deltax=test_deltax;
				      max_deltay=test_deltay;

				  }
			      }
			  } // for loop old_state
			  if ( max_score>Logzero )
			  {
			      if (viterbi_rectangle[new_state][xsteps+offset-x_start][ysteps+offset-y_start] !=
				  Logzero) {
				  
				  cout << "ERROR : viterbi_rectangle[" << new_state 
				       << "][" << xsteps+offset-x_start
				       << "][" << ysteps+offset-y_start
				       << "] = "
				       << viterbi_rectangle[new_state][xsteps+offset-x_start][ysteps+offset-y_start]
				       << " != Logzero\n"
				       << flush;
			      }
			      viterbi_rectangle[new_state]
				  [xsteps+offset-x_start]
				  [ysteps+offset-y_start]=max_score;

			      max_count++;

			  }
			  
			  if (indices) delete [] indices;
			  indices = NULL;


		      } // if ((xsteps>=deltax) && (ysteps>=deltay))
		      else
		      {

		      }
		  } // for loop new_state
	      } // if not ( (xsteps==0) && (ysteps==0) )
	  } // for loop ysteps

	} // for loop xsteps
	
	if (direction==1)
	{
	    xsteps-=1;
	    ysteps-=1;
	}
	else if (direction==-1)
	{
	    xsteps+=1;
	    ysteps+=1;
	}
 
	if (  ((direction==1)          &&
	       (xsteps==x->length()) &&   
	       (ysteps==y->length())) 
	      ||
	      ((direction==-1)         &&
	       (xsteps==-1) &&  
	       (ysteps==-1))  )
	{	

	    max_score=Logzero;
	    max_state=0;
	    max_deltax=0;
	    max_deltay=0;
	    
	    deltax=0;
	    deltay=0;
	    
	    for (int state=0; state<states; state++)  // loop over all possible previous states
	    {
		pre_test_score=Logzero;
		test_score=Logzero;
		test_state=0;
		test_deltax=0;
		test_deltay=0;
		
		int x_position = 0;
		int y_position = 0;
		
		if (direction == 1)
		{
		    x_position = xsteps - deltax;
		    y_position = ysteps - deltay;
		}
		else if (direction == -1)
		{
		    x_position = xsteps + 1;
		    y_position = ysteps + 1;
		}
		
		pre_test_score = s->new_get_transition_score(s_unmirrored, 
							     state, states-1, 
							     x, x_position, y, y_position);

		if (pre_test_score > Logzero)
		{
		    test_state=state;
		    test_score=pre_test_score+
			viterbi_rectangle[state][xsteps-direction*deltax+offset-x_start][ysteps-direction*deltay+offset-y_start];
		    test_deltax=0;
		    test_deltay=0;

		    if (test_score>max_score)
		    {

			max_score=test_score;
			max_state=test_state;
			max_deltax=test_deltax;
			max_deltay=test_deltay;

		    }
		}
	    }
	    if ( max_score>Logzero )
	    {
		viterbi_rectangle[states-1][xsteps+offset-x_start][ysteps+offset-y_start]=max_score;
	    }
	}      
    }
    return(check);
}

int Hmm::retrieve_state_path_from_viterbi_rectangle(Score*** const viterbi_rectangle,
							const Hmm *s, 
							const Hmm *s_unmirrored,
							const Sequence *x, const int x_start, const int x_end,
							const Sequence *y, const int y_start, const int y_end,
							const int end_state, 
							const int direction,
							int& local_steps,
							Score& local_score,
							Score* const local_sequence_of_scores,
							int* const local_sequence_of_states,
							int* const local_sequence_of_xsteps,
							int* const local_sequence_of_ysteps)
{							 
    // note: this function may only used on un-mirrored pairhmm (this is checked) and on a viterbi_rectangle
    //       which has been filled using this un-mirrored pairhmm (not checked)
    
    int check=0;
    int n_of_states=0;
    n_of_states=s->get_number_of_states();
    
    // check viterbi_rectangle
    
    if (viterbi_rectangle==NULL)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: viterbi_rectangle (" 
	    << viterbi_rectangle << ") = NULL\n" << flush;
	check+=1;
    }
 
    // check local_sequences that will be filled with the state path
    
    if (local_sequence_of_scores==NULL)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: local_sequence_of_scores (" 
	     << local_sequence_of_scores << ") = NULL\n" << flush;
	check+=1;
    }
    
    if (local_sequence_of_states==NULL)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: local_sequence_of_states (" 
	     << local_sequence_of_states << ") = NULL\n" << flush;
	check+=1;
    }

    if (local_sequence_of_xsteps==NULL)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: local_sequence_of_xsteps (" 
	     << local_sequence_of_xsteps << ") = NULL\n" << flush;
	check+=1;
    }

    if (local_sequence_of_ysteps==NULL)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: local_sequence_of_ysteps (" 
	     << local_sequence_of_ysteps << ") = NULL\n" << flush;
	check+=1;
    }
    
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle : length of sequence x (" 
	     << x->length() << ") < 1 or x (" << x << ") NULL \n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle : length of sequence y (" 
	     << y->length() << ") < 1 or y (" << y << ") NULL \n" << flush;
	check+=1;
    }
    
    
  // check Hmm
    
    if (s->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : number of states = " 
	     << s->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (s->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : state 0 of type != Start \n" << flush;
	cout << "s->model[0].get_letters_to_read() : "<<s->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (s->model[n_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : last state of type != End \n" << flush;
	cout << "s->model["<<n_of_states-1<<"].get_letters_to_read() : "<<s->model[n_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (int i=0; i<n_of_states; i++)
    {
	if (s->model[i].get_alphabet() != alphabet) 
	{
	    cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << s->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_mirrored() != s->get_mirrored()) 
	{
	    cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << s->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << s->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << s->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << n_of_states << "\n" << flush;
	    check+=1;
	}
	
	if ((i!=0) && (i!=(n_of_states-1)) && (s->model[i].get_letters_to_read()==0 ))
	{
	    cout << "ERROR: Hmm::retrieve_state_path_from_viterbi_rectangle :  Hmm : state i = " 
		 << i << " should emit ! \n" << flush;
	    check+=1;
	}

    }
    
    // check x_start, x_end and y_start, y_end values

    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    // check end state
    
    if ((x_end==x->length()) && (y_end==y->length()))
    {
	if ((direction==1) && (end_state!=n_of_states-1))
	{
	    cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: direction==1: end_state (" 
		 << end_state 
		 << ") for x_end (" << x_end << ") = x->length() (" << x->length() << ") and y_end ("
		 << y_end << ") = y->length() (" << y->length() << ") has to be End state ("
		 << n_of_states-1 << ").\n" << flush;
	    check+=1;
	}
	
    }
    else // if we don't end at x_end=x->length() and y_end=y->length()
    {
	// check that each end state lies in the allowed range (1...n_of_states-2)
	// i.e. Start and End states are not allowed (if direction==1)
	
	if ((direction==1) && ((end_state<1) || (end_state>n_of_states-2)))
	{
	    cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: end_state (" 
		<< end_state << ") out of range [1, " << n_of_states-2 << "].\n" << flush;
	    check+=1;
	}
    }
    
    if ((x_start==0) && (y_start==0))
    {
	if ((direction==-1) && (end_state!=0))
	{
	    cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: direction==-1: end_state (" 
		 << end_state 
		 << ") for x_start (" << x_start << ") = 0 and y_start ("
		 << y_start << ") = 0 has to be Start state, 0.\n" << flush;
	    check+=1;
	}
    }
    
    // check direction
    
    if ((direction!=1) && (direction!=-1))
    {
	cout << "ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: direction (" << direction 
	     << ") has to be 1 or -1.\n" << flush;
	check+=1;
    }
    
    // enter program if all check were successfully passed
    
    if (check==0)
    {

	int shift_index = bitshift(alphabet);
	
	local_score=Logzero;
	local_steps=0;

	int test_state=0;
	Score test_score=Logzero;
	Score pre_test_score=Logzero;
	Score test_score_2=Logzero;
	int new_state=0;
	
	Score max_score=Logzero;
	Score max_score_2=Logzero;
	int max_state=0;
	int deltax=0;
	int deltay=0;

	int start_x=0;
	int end_x=0;
	int start_y=0;
	int end_y=0;
	int offset=0;
	
	if (direction==1)
	{
	    offset=0;
	    
	    start_x=x_start; 
	    end_x=x_end+1;
	    start_y=y_start;
	    end_y=y_end+1; 
	}
	else if (direction==-1) // permute start and end points if running backwards
	{
	    offset=1;
	    
	    start_x=x_end-1;
	    start_y=y_end-1;
	    
	    end_x=x_start-2;
	    end_y=y_start-2;
	}
	
	// calculate actual xsteps and ysteps values
	
	int xsteps=end_x;
	int ysteps=end_y;

	xsteps+=(-1)*direction;
	ysteps+=(-1)*direction;

	// implement stop constraint
      
	if (direction==1)
	{
	    new_state=end_state;
	}
	else if (direction==-1)
	{
	    new_state=n_of_states-1-end_state;
	} 
	
	local_score=viterbi_rectangle[new_state][xsteps+offset-x_start][ysteps+offset-y_start];

	if (local_score>Logzero) // if this state path is possible
	{      
	    Score remaining_score=local_score;

	    local_sequence_of_states[0]=new_state;
	    local_sequence_of_scores[0]=local_score;
	    
	    deltax=s->model[new_state].get_letters_to_read_x();
	    deltay=s->model[new_state].get_letters_to_read_y();
  
	  int boundaries_reached=0;

	  while (( ! (((xsteps == start_x) && (ysteps == start_y)) || (boundaries_reached==1))) &&
		 (check == 0))
	    // as long as we haven't reached the start or the boundaries of rectangle
	    {
		max_score=Logzero;
		max_score_2=Logzero;
		max_state=0;  
		if (((xsteps-direction*deltax+offset-x_start) > -1) && 
		    ((ysteps-direction*deltay+offset-y_start) > -1) && 
		    ((xsteps-direction*deltax+offset-x_start) < (x_end-x_start+1)) && 
		    ((ysteps-direction*deltay+offset-y_start) < (y_end-y_start+1)))
		    // 
		    // if old_state -> new_state transition is possible within dimensions of this rectangle     
		{
		    int number_of_old_states = s->model[new_state].get_number_of_previous_states();
		    
		    for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
		    {
			int old_state = s->model[new_state].get_number_of_previous_state(old_state_number);
			
			test_score=Logzero;
			test_score_2=Logzero;
			test_state=0;

			int linear_index_emission = 0;
			int* indices     = NULL;
			indices = new int[model[new_state].get_number_of_dimensions_of_emission_scores()+1];
			
			int xcount=0;
			int ycount=0;
			
			if (direction==1) {
			    xcount=deltax;
			    ycount=0;
			}
			else if (direction==-1) {
			    xcount=0;
			    ycount=deltay;
			}
			for (int i=ycount; i<deltax+ycount; i++) {
			    indices[i+1] = x->letter(xsteps - direction*(deltax+ycount-i));

			}
			for (int j=xcount; j<deltay+xcount; j++) {
			    indices[j+1] = y->letter(ysteps - direction*(deltay+xcount-j)); 

			}
			indices[0] = old_state;
			
			int x_position = 0;
			int y_position = 0;
			
			if (direction == 1) {
			    x_position = xsteps - deltax;
			    y_position = ysteps - deltay;
			}
			else if (direction == -1) {
			    x_position = xsteps + 1;
			    y_position = ysteps + 1;
			}

			if (deltax+deltay > 0)
			{
			    linear_index_emission = indices[1];
			    for (int p=2; p<s->model[new_state].get_number_of_dimensions_of_emission_scores()+1; p++)
			    {
				linear_index_emission = (linear_index_emission << shift_index) | indices[p];
			    }
			    
			    pre_test_score=s->new_get_transition_score(s_unmirrored, 
								       old_state, new_state, 
								       x, x_position, y, y_position) +
				s->model[new_state].get_emission_score(x, x_position, y, y_position,
								       linear_index_emission);
			}
			else
			{
			    pre_test_score=s->new_get_transition_score(s_unmirrored, 
								       old_state, new_state, 
								       x, x_position, y, y_position);
			}
			
			test_state=old_state;
			test_score_2=viterbi_rectangle[old_state][xsteps-direction*deltax+offset-x_start][ysteps-direction*deltay+offset-y_start];
			test_score=pre_test_score + test_score_2;

			if (test_score>max_score)
			{

			    max_score=test_score;
			    max_score_2=test_score_2;
			    max_state=test_state;

			}

			if (indices) delete [] indices;
			indices = NULL;
			
		    } // for loop over old_state
		    if ( max_score>Logzero )
		    {
			local_sequence_of_scores[local_steps]=remaining_score-max_score_2;
			local_sequence_of_states[local_steps+1]=max_state;
			local_sequence_of_scores[local_steps+1]=max_score_2;		  

			new_state=max_state;
			xsteps-=direction*deltax;
			ysteps-=direction*deltay;
			deltax=s->model[new_state].get_letters_to_read_x();
			deltay=s->model[new_state].get_letters_to_read_y();
			remaining_score=max_score_2;
			local_steps++;

		    }
		    else // if no old_state had score > Logzero
		    {
			cout << "   ERROR Hmm::retrieve_state_path_from_viterbi_rectangle: traceback: all previous states have score Logzero.\n" << flush;
			check+=1;
		    }

		} // if transition old_state -> new_state is possible within the dimensions of this rectangle 
		else
		{

		    boundaries_reached=1;
		}
	    }
	  
	  if (check == 0)
	  {
	      local_sequence_of_states[local_steps]=new_state;
	      local_sequence_of_scores[local_steps]=0; 

	      Score copy_score=0;
	      int copy_state=0;
	      
	      if (direction==1)
	      {
		  int end_of_loop=static_cast<int>(ceil(static_cast<float>(local_steps)/2.));

		  copy_score=local_sequence_of_scores[local_steps];
		  copy_state=local_sequence_of_states[local_steps];
		  
		  local_sequence_of_states[local_steps]=local_sequence_of_states[0];
		  local_sequence_of_scores[local_steps]=local_sequence_of_scores[0];
	      
		  for (int i=1; i<end_of_loop; i++) // put sequence_of_scores/states in right order
		  {
		      local_sequence_of_states[0]=local_sequence_of_states[i];
		      local_sequence_of_scores[0]=local_sequence_of_scores[i];
		      
		      local_sequence_of_states[i]=local_sequence_of_states[local_steps-i];
		      local_sequence_of_scores[i]=local_sequence_of_scores[local_steps-i];
		      
		      local_sequence_of_states[local_steps-i]=local_sequence_of_states[0];
		      local_sequence_of_scores[local_steps-i]=local_sequence_of_scores[0];
		  }
		  
		  local_sequence_of_states[0]=copy_state;
		  local_sequence_of_scores[0]=copy_score;  
	      }  
	      
	      Score score_sum=0;
	      int steps_x=xsteps; 
	      int steps_y=ysteps; 
	      
	      int backward_state=0;

	      for (int i=0; i<local_steps+1; i++)               
	      {
		  score_sum+=local_sequence_of_scores[i];
		  
		  if (direction==1)
		  {
		      backward_state=local_sequence_of_states[i];
		  }
		  else if (direction==-1)
		  {
		      local_sequence_of_states[i]=n_of_states-1-local_sequence_of_states[i];
		      backward_state=n_of_states-1-local_sequence_of_states[i];
		  }	      
		  
		  local_sequence_of_xsteps[i]=steps_x;
		  local_sequence_of_ysteps[i]=steps_y;
		  
		  if (i < local_steps)
		  {
		      steps_x+=s->model[local_sequence_of_states[i+1]].get_letters_to_read_x();
		      steps_y+=s->model[local_sequence_of_states[i+1]].get_letters_to_read_y();
		  } 
	      }
	  } // if check == 0
	} // if state path is possible, i.e. if local_score>Logzero
	else
	{
	    cout << "NOTE: Hmm::retrieve_state_path: no state path can be found for end state " 
		 << end_state << "\n" << flush;
	    check+=1;
	}
    }
    return(check);
}

int Hmm::get_strips_and_merge(// input values
				  const Hmm *mirror,         
				  const Sequence *x, const int x_start, const int x_end,
				  const Sequence *y, const int y_start, const int y_end,
				  const int start_state, const int end_state,
				  const int x_midpoint, const int min_strip_width,
				  // output values
				  int* max_i, int* max_j, int* max_state_forward, int* max_state_backward)

{
    int check=0;
    
    // check Hmm
  
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : Hmm : mirrored (" << mirrored 
	     << ") must be 0.\n" << flush;
	check+=1;
    }
    
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : state 0 of type != Start \n" << flush;
	cout<<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : last state of type != End \n" << flush;
	cout<<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : alphabet of state i = " << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : mirrored of state i = " << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	   
	    if ((i!=0) && (i!=(number_of_states-1)) && (model[i].get_letters_to_read()==0))
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : state i = " << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : mirrored (" 
	     << mirror->get_mirrored() << ") must be 1.\n" << flush;
	check+=1;
    }
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : state 0 of typeStart \n" << flush;
	cout<<"mirror->model[0].get_letters_to_read() : "<<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : last state of type != End \n" << flush;
	cout<<"mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }
    
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!=mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	  
	    if ((i!=0) && (i!=(number_of_states-1)) && (mirror->model[i].get_letters_to_read()==0))
	    {
		cout << "ERROR: Hmm::get_strips_and_merge : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    // check that this Hmm and the mirror Hmm belong together
    
    if (number_of_states != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::get_strips_and_merge : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (int i=1; i<number_of_states-1; i++)
	{
	    if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read()){
		cout << "ERROR: Hmm::get_strips_and_merge : letters_to_read of  model[" << i << "] (" 
		     << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		     << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::get_strips_and_merge: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::get_strips_and_merge: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::get_strips_and_merge: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::get_strips_and_merge: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    // check start and end states
    
    if ((x_start==0) && (y_start==0) && (start_state!=0))
    {
	cout << "ERROR Hmm::get_strips_and_merge: start_state (" << start_state 
	     << ") for x_start=y_start=0 should be 0.\n" << flush;
	check+=1;
    }
    if ((x_start==x->length()) && (y_start==y->length()) && (end_state!=number_of_states-1))
    {
	cout << "ERROR Hmm::get_strips_and_merge: end_state (" << end_state 
	     << ") for x_start = x->length() (" << x->length() << ") and y_start = y->length() (" 
	     << y->length() << ") should be 0.\n" << flush;
	check+=1;
    }
    if ((start_state<0) || (start_state>number_of_states-1))
    {
	cout << "ERROR Hmm::get_strips_and_merge: start_state (" << start_state << ") out of range [0, " 
	     << number_of_states-1 << "]\n" << flush;
	check+=1;
    }
    if ((end_state<0) || (end_state>number_of_states-1))
    {
	cout << "ERROR Hmm::get_strips_and_merge: end_state (" << end_state << ") out of range [0, " 
	     << number_of_states-1 << "]\n" << flush;
	check+=1;
    }
    
    // check min_strip_width
    
    if (min_strip_width<1)
    {
	cout << "ERROR Hmm::get_strips_and_merge: min_strip_width (" << min_strip_width 
	     << ") has to be >= 1.\n" << flush;
	check+=1;
    }
    
  // check x_midpoint
    
    if ((x_midpoint-(min_strip_width+1))<0)
    {
	cout << "ERROR Hmm::get_strips_and_merge: x_midpoint (" << x_midpoint << ") - (min_strip_width+1) ("
	     << min_strip_width+1 << ") = " << x_midpoint-(min_strip_width+1) << " has to be >=0.\n" << flush;
	check+=1;
    }
    if ((x_midpoint<x_start) || (x_midpoint>x_end))
    {
	cout << "ERROR Hmm::get_strips_and_merge: x_midpoint (" << x_midpoint 
	     << ") out of range [x_start (" << x_start << "), x_end (" << x_end << ")].\n" << flush;
	check+=1;
    }

    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::get_strips_and_merge : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::get_strips_and_merge : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }

    // declare and initialise return values
    
    *max_i=0;
    *max_j=0;
    *max_state_forward=0;
    *max_state_backward=0;
    
    if (check==0)
    {

	// for each state initialize its viterbi_strip and set strip's dimensions

      
	Score*** forward_strip= new Score**[number_of_states];
	Score*** backward_strip= new Score**[number_of_states];
	
	for (int state=0; state<number_of_states; state++)
	{
	    forward_strip[state]= new Score*[min_strip_width+1];
	    backward_strip[state]= new Score*[min_strip_width+1];
	    
	    for (int i=0; i<min_strip_width+1; i++)
	    {
		forward_strip[state][i]= new Score[y_end-y_start+1];
		backward_strip[state][i]= new Score[y_end-y_start+1];
		
		for (int j=0; j<y_end-y_start+1; j++)
		{	     
		    forward_strip[state][i][j]= Logzero;
		    backward_strip[state][i][j]= Logzero;
		}
	    }	 
	}
	
	// implement constraints into forward and into backward strip
	
	forward_strip[start_state][min_strip_width][0]=0;  
	backward_strip[number_of_states-1-end_state][min_strip_width][y_end-y_start]=0; 

	int x_margin_forward = 0;
	if (x_start != 0) {x_margin_forward = min_strip_width+1;}
	
	int x_margin_backward = 0;
	if (x_midpoint-min_strip_width-1 != 0) {x_margin_backward = min_strip_width+1;}
	
	if (check == 0) {
	    check += get_strip(&forward_strip, 
			       this,
			       NULL, // no unmirrored Hmm needed
			       x, x_start, x_midpoint-1, 
			       y, y_start, y_end, 
			       min_strip_width+1,
			       x_margin_forward, 0,
			       1);
	    if (check != 0) {
		cout << "ERROR Hmm::get_strips_and_merge : error occurred when calculating the forward strip "
		     << "(x_start (" << x_start << "), x_midpoint-1 (" << x_midpoint-1 
		     << ")) (y_start (" << y_start << "), y_end (" << y_end << "))\n" << flush;	  
	    }
	}
	
	
	if (check == 0) {
	    check += get_strip(&backward_strip, 
			       mirror,
			       this, // unmirrored Hmm needed
			       x, x_midpoint-min_strip_width-1, x_end, 
			       y, y_start, y_end, 
			       min_strip_width+1,
			       x_margin_backward, 0,
			       -1);
	    if (check != 0) {
		cout << "ERROR Hmm::get_strips_and_merge : error occurred when calculating the backward strip "
		     << "(x_midpoint-min_strip_width-1 (" << x_midpoint-min_strip_width-1 
		     << "), x_end (" << x_end
		     << ")) (y_start (" << y_start << "), y_end (" << y_end << "))\n" << flush;	  
	    }
	}
	
	if (check==0)
	{
	    // check results
	 
	    Score max_score=Logzero;
	    Score test_score=Logzero;
	    Score transition_score=Logzero;
	    Score forward_score=Logzero;
	    Score backward_score=Logzero;
	    
	    *max_i=0;
	    *max_j=0;
	    *max_state_forward=0;
	    *max_state_backward=0;
	    
	    Score max_transition_score=0;
	    Score max_forward_score=Logzero;
	    Score max_backward_score=Logzero;
	    
	    for (int i=0; i<min_strip_width+1; i++)
	    {
		for (int j=0; j<y_end-y_start+1; j++) 
		{

		    for (int new_state=0; new_state<number_of_states; new_state++)
		    {

			int number_of_old_states = model[new_state].get_number_of_previous_states();
			  
			for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
			{
			    int old_state = model[new_state].get_number_of_previous_state(old_state_number);
			    
			    transition_score = Logzero;
			    
			    int x_position = (i + x_midpoint - min_strip_width) - this->model[new_state].get_letters_to_read_x();
			    int y_position = (j + y_start + 1)                  - this->model[new_state].get_letters_to_read_y();


			    // note: this must be un-mirrored pairhmm
	    
			    transition_score = this->new_get_transition_score(NULL, 
									      old_state, new_state,
									      x, x_position,
									      y, y_position);

			    if (transition_score > Logzero)
			    {
				test_score=Logzero;
				
				test_score=forward_strip[old_state][i][j]+
				    backward_strip[new_state][i][j]+
				    transition_score;
				
				forward_score=forward_strip[old_state][i][j];
				backward_score=backward_strip[new_state][i][j];

				if (test_score > max_score)
				{
				    max_score=test_score;
				    *max_i=i;
				    *max_j=j;
				    *max_state_forward=old_state;
				    *max_state_backward=new_state;
				    
				    max_forward_score=forward_score;
				    max_backward_score=backward_score;
				    max_transition_score=transition_score;

				}
			    }
			} // for old_state
		    } // for new_state
		} // for j (y values)
	    } // for i (x value in strip)

	}
	
	// for each state delete its viterbi_strip

	if (forward_strip && backward_strip)
	{      
	    for (int state=0; state<number_of_states; state++)
	    {
		for (int i=0; i<min_strip_width+1; i++)
		{
		    if (forward_strip[state][i])  delete [] forward_strip[state][i];
		    if (backward_strip[state][i]) delete [] backward_strip[state][i];
		    forward_strip[state][i]  = NULL;
		    backward_strip[state][i] = NULL;
		}
		if (forward_strip[state])  delete [] forward_strip[state];
		if (backward_strip[state]) delete [] backward_strip[state];
		forward_strip[state]  = NULL;
		backward_strip[state] = NULL;
	    }
	    if (forward_strip)  delete [] forward_strip;
	    if (backward_strip) delete [] backward_strip;
	    forward_strip=NULL;
	    backward_strip=NULL;
	}
    }
    return(check);
}

int Hmm::new_get_strips_and_merge(// input values
				      Score***         forward_viterbi_strip,
				      Score***         backward_viterbi_strip, 
				      const Hmm*   mirror,
				      const Sequence*  x, const int x_start, const int x_end,
				      const Sequence*  y, const int y_start, const int y_end,
				      const int*       start_states, const int number_of_start_states,
				      const int*       end_states,   const int number_of_end_states,
				      const int        x_midpoint,				      
				      const int        x_margin,
				      const int        y_margin,
				      // in- and output values
				      int* min_strip_width,
				      // output values
				      int* max_i, int* max_j, int* max_state_forward, int* max_state_backward)
{
    // function returns 0 if o.k., else not o.k.
    // note : 
    //
    //
    // note: - if both forward_ and backward_viterbi_strip are used make sure that they
    //         have the same strip width (this cannot be checked !)
    // 
    //       - this function can either be used 
    // (
    //       1.) with start_states (in which case number_of_start_states > 0). Memory
    //       for forward_viterbi_rectangle_array is allocated and deleted within function.
    // 
    //  or
    //       
    //       2.) with a pre-initialised forward_viterbi_rectangle_array (in which case
    //       number_of start_states == 0).
    //       Memory for forward_viterbi_rectangle_array is allocated and deleted outside function.
    //
    // )
    // and
    // (
    //       3.) with end_states (in which case number_of_end_states > 0). Memory
    //       for backward_viterbi_rectangle_array is allocated and deleted within function.
    // 
    //  or
    //       
    //       4.) with a pre-initialised backward_viterbi_rectangle_array (in which case 
    //       number_of end_states == 0).
    //       Memory for backward_viterbi_rectangle_array is allocated and deleted outside function.
    // )
    
    int check=0;
    int i, j = 0;
    
    // check whether pre-initialised forward_viterbi_rectangle or start_states shall be used
    
    if ((forward_viterbi_strip == NULL)    &&
	(number_of_start_states == 0))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: forward_viterbi_strip (" << forward_viterbi_strip
	     << ") = NULL and number_of_start_states (" << number_of_start_states 
	     << ") == 0. Use this function either with a pre-initialised forward_viterbi_strip != NULL "
	     << "or with nonempty sets of start_states.\n" << flush;
	check+=1;
    }
    else if ((forward_viterbi_strip != NULL)    &&
	   (number_of_start_states != 0))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: forward_viterbi_strip (" << forward_viterbi_strip
	     << ") != NULL and number_of_start_states (" << number_of_start_states 
	     << ") != 0. Use this function either with a pre-initialised forward_viterbi_strip != NULL "
	     << "or with nonempty sets of start_states.\n" << flush;
	check+=1;
    }
    else // if input values are not completely inconsistent
    {
	if ((forward_viterbi_strip      == NULL)   &&
	    (number_of_start_states != 0))
	{
	    // if function shall be used with start_states
	    // check start states
	    
	    if (number_of_start_states<1)
	    {
		cout << "ERROR Hmm::new_get_strips_and_merge: number_of_start_sites (" 
		     << number_of_start_states << ") <1.\n" << flush;
		check+=1;
	    }
	    else
	    {
		if ((x_start==0) && (y_start==0))
		{
		    if ((number_of_start_states!=1) || (start_states[0]!=0))
		    {
			cout << "ERROR Hmm::new_get_strips_and_merge: or "
			     << "x_start=y_start=0 only "
			     << "one start_state is possible, the Start state (i.e. the number_of_start_states ("
			     << number_of_start_states << ") has to be 1 and start_states[0] ("
			     << start_states[0] << ") has to be 0).\n" << flush;
			check+=1;
		    }
		}
		else // if we don't Start at (0,0)
		{
		    // check that each start state lies in the allowed range (1...n_of_states-2)
		    // i.e. Start and End states are not allowed 
		    
		    for (i=0; i<number_of_start_states; i++)
		    {
			if ((start_states[i]<1) || (start_states[i]>this->get_number_of_states()-2))
			{
			    cout << "ERROR Hmm::new_get_strips_and_merge: start_states[" << i 
				 << "] (" << start_states[i] << ") out of range [1, " 
				 << this->get_number_of_states()-2
				 << "].\n" << flush;
			    check+=1;
			}
		    }
		}
	    }
	}
	
	// check x_margin and y_margin values
	
	if (! (((number_of_start_states > 0) && (number_of_end_states > 0)))) 
	{
	    if (x_margin < 0)
	    {
		cout << "ERROR Hmm::new_get_strips_and_merge: x_margin (" << x_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	    
	    if (y_margin < 0)
	    {
		cout << "ERROR Hmm::new_get_strips_and_merge: y_margin (" << y_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	    
	    // check that x_margin and y_margin are compatible with values of x_start, x_end and y_start and y_end
	    
	    if ((x_end-x_start+1) < x_margin)
	    {
		cout << "ERROR Hmm::new_get_strips_and_merge: x_margin (" << x_margin
		     << ") > (x_end (" << x_end << ") - x_start (" << x_start << ") + 1) = " 
		     << x_end-x_start+1 << "\n" << flush;
		check+=1;
	    }
	
	    if ((y_end-y_start+1) < y_margin)
	    {
		cout << "ERROR Hmm::new_get_strips_and_merge: y_margin (" << y_margin
		     << ") > (y_end (" << y_end << ") - y_start (" << y_start << ") + 1) = " 
		     << y_end-y_start+1 << "\n" << flush;
		check+=1;
	    }
	}
    }
    
    // check whether pre-initialised backward_viterbi_rectangle or end_states shall be used
    
    if ((backward_viterbi_strip == NULL)    &&
	(number_of_end_states == 0))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: backward_viterbi_strip (" << backward_viterbi_strip
	     << ") = NULL and number_of_end_states (" << number_of_end_states 
	     << ") == 0. Use this function either with a pre-initialised backward_viterbi_strip != NULL "
	     << "or with nonempty sets of end_states.\n" << flush;
	check+=1;
    }
    else if ((backward_viterbi_strip != NULL)    &&
	     (number_of_end_states != 0))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: backward_viterbi_strip (" << backward_viterbi_strip
	     << ") != NULL and number_of_end_states (" << number_of_end_states 
	     << ") != 0. Use this function either with a pre-initialised backward_viterbi_strip != NULL "
	     << "or with nonempty sets of end_states.\n" << flush;
	check+=1;
    }
    else // if input values are not completely inconsistent
    {
	if ((backward_viterbi_strip      == NULL)   &&
	    (number_of_end_states != 0))
	{
	  // if function shall be used with end_states
	  // check start  states
	    
	    if (number_of_end_states<1)
	    {
		cout << "ERROR Hmm::new_get_strips_and_merge: number_of_end_sites (" 
		     << number_of_end_states << ") <1.\n" << flush;
		check+=1;
	    }
	    else
	    {
		// check end_states
		
		if ((x_end==x->length()) && (y_end==y->length()))
		{
		    if (! ((number_of_end_states == 1) && (end_states[0] == this->get_number_of_states()-1)))
		    {
			cout << "ERROR Hmm::new_get_strips_and_merge: number_of_end_states (" 
			     << number_of_end_states 
			     << ") for x_end (" << x_end << ") = x->length() (" << x->length() 
			     << ") and y_end (" << y_end << ") = y->length() (" << y->length() 
			     << ") has to be 1 and end_states[0] ("
			     << end_states[0] << ") has to be End state (" 
			     << this->get_number_of_states()-1 << ").\n" << flush;
			check+=1;
		    }
		}
		else // if we don't stop at ends of both sequences
		{
		    // check that end_states lie in the allowed range (1...n_of_states-2)
		    // i.e. Start  and End states are not allowed
		    
		    for (i=0; i<number_of_end_states; i++) 
		    {
			if ((end_states[i]<1) || (end_states[i]>this->get_number_of_states()-2))
			{
			    cout << "ERROR Hmm::new_get_strips_and_merge: end_states[" << i << "] "
				 << end_states[i] << ") out of range [1, " 
				 << this->get_number_of_states()-2 << "].\n" << flush;
			    check+=1;
			}
		    }
		}
	    }
	}
    }
    
    
    // check Hmm
    
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : Hmm : mirrored (" << mirrored 
	     << ") must be 0.\n" << flush;
	check+=1;
    }
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : state 0 != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : last state != End \n" << flush;
	cout <<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : alphabet of state i = " << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored of state i = " << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1))&&(model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : state i = " << i << " should emit! \n";
		check+=1;
	    }
	   
	}
    }
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : mirrored (" 
	     << mirror->get_mirrored() << ") must be 1.\n" << flush;
	check+=1;
    }
    
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : state 0 of type  != Start \n" << flush;
	cout<<"mirror->model[0].get_letters_to_read() : "
	    <<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : last state != End \n" << flush;
	cout<<"mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!=mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1))&& (mirror->model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	   
	}
    }
    // check that this Hmm and the mirror Hmm belong together
    
    if (number_of_states != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (i=1; i<number_of_states-1; i++)
	{
	    if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::new_get_strips_and_merge : letters_to_read of model[" << i << "] (" 
		     << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		     << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    // check min_strip_width
    
    if ((number_of_start_states > 0) && 
	(number_of_end_states > 0)   &&
	((*min_strip_width)<1))
    {
	cout << "ERROR Hmm::new_get_strips_and_merge: min_strip_width (" << (*min_strip_width) 
	     << ") has to be > 0.\n" << flush;
	check+=1;
    }
    
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::new_get_strips_and_merge : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

	// initialise return values
	
	*max_i=0;
	*max_j=0;
	*max_state_forward=0;
	*max_state_backward=0;
	
	// internal variables
	
	int internal_strip_width       = 0;
	int internal_x_margin_forward  = 0;
	int internal_x_margin_backward = 0;
	int internal_y_margin_forward  = 0;
	int internal_y_margin_backward = 0;
	
	internal_strip_width = max((*min_strip_width)+1, x_margin);
	
	if (number_of_start_states == 0) { // if forward_rectangle is used with initialised strip
	    
	    if (x_start != 0) {
		internal_x_margin_forward = x_margin; 
	    }
	    if (y_start != 0) {
		internal_y_margin_forward = y_margin;
	    }
	}
	if (number_of_end_states == 0) { // if backward_rectangle is used with initialised strip
	    
	    if (x_end != x->length()) {
		internal_x_margin_backward = x_margin;
	    }
	    if (y_end != y->length()) {
		internal_y_margin_backward = y_margin;
	    }
	}
	
	// assign value to output variable
	
	(*min_strip_width) = internal_strip_width - 1;

	Score*** forward_strip  = NULL;
	Score*** backward_strip = NULL;
	
	int state = 0;
	
	// create forward_strip and implement constraints into forward_strip
	
	if (number_of_start_states == 0) {
	    
	    forward_strip  = forward_viterbi_strip;
	}
	else {
	    
	    // for each state initialize its viterbi_strip and set strip's dimensions

	    forward_strip= new Score**[number_of_states];
	    
	    for (state=0; state<number_of_states; state++) {
		
		forward_strip[state]= new Score*[internal_strip_width];
		
		for (i=0; i<internal_strip_width; i++) {
		    
		    forward_strip[state][i]= new Score[y_end-y_start+1];
		    
		    for (j=0; j<y_end-y_start+1; j++) {
			
			forward_strip[state][i][j]= Logzero;
			
		    }
		}	 
	    }
	    
	    // implement constraints into forward_strip

	    for (i=0; i<number_of_start_states; i++) {
		
		forward_strip[start_states[i]][internal_strip_width-1][0]=0; 

	    }
	} // else: number_of_start_states != 0
	
	// create backward_strip and implement constraints into backward_strip
	
	if (number_of_end_states == 0) {
	    
	    backward_strip  = backward_viterbi_strip;
	}
	else {
	    
	    // for each state initialize its viterbi_strip and set strip's dimensions

	    backward_strip= new Score**[number_of_states];
	    
	    for (state=0; state<number_of_states; state++) {
		
		backward_strip[state]= new Score*[internal_strip_width];
		for (i=0; i<internal_strip_width; i++) {
		    
		    backward_strip[state][i]= new Score[y_end-y_start+1];
		    for (j=0; j<y_end-y_start+1; j++) {
			
			backward_strip[state][i][j]= Logzero;
		    }
		}	 
	    }
	    
	    // implement constraints into backward_strip

	    for (i=0; i<number_of_end_states; i++) {
		
		backward_strip[number_of_states-1-end_states[i]][internal_strip_width-1][y_end-y_start]=0;

	    }
	} // else : number_of_end_states != 0 
	
	// get strips for two rectangles

	if (check == 0) {

	    check+=get_strip(&forward_strip, 
			     this,
			     NULL, // no unmirrored Hmm needed
			     x, x_start, x_midpoint-1, 
			     y, y_start, y_end, 
			     internal_strip_width,
			     internal_x_margin_forward, 
			     internal_y_margin_forward, 
			     1);

	    if (check != 0) {
		cout << "ERROR Hmm::new_get_strips_and_merge : error occurred when calculating the forward strip "
		     << "(x_start (" << x_start << "), x_midpoint-1 (" << x_midpoint-1 
		     << ")) (y_start (" << y_start << "), y_end (" << y_end << "))\n" << flush;	  
	    }
	}
	
	if (check == 0) {
	    
	    check+=get_strip(&backward_strip, 
			     mirror,
			     this, // unmirrored Hmm needed
			     x, x_midpoint-internal_strip_width, x_end, 
			     y, y_start, y_end, 
			     internal_strip_width,
			     internal_x_margin_backward, 
			     internal_y_margin_backward, 
			     -1);
	    
	    if (check != 0) {
		cout << "ERROR Hmm::new_get_strips_and_merge : error occurred when calculating the backward strip "
		     << "(x_midpoint-min_strip_width-1 (" << x_midpoint-internal_strip_width
		     << "), x_end (" << x_end
		     << ")) (y_start (" << y_start << "), y_end (" << y_end << "))\n" << flush;	  
	    }
	}
	
	if (check==0)
	{
	    
	    Score max_score=Logzero;
	    Score test_score=Logzero;
	    Score transition_score=Logzero;
	    Score forward_score=Logzero;
	    Score backward_score=Logzero;
	    
	    *max_i=0;
	    *max_j=0;
	    *max_state_forward=0;
	    *max_state_backward=0;
	    
	    Score max_transition_score=0;
	    Score max_forward_score=Logzero;
	    Score max_backward_score=Logzero;
	    
	    for (i=0; i<internal_strip_width; i++)
	    {
		for (j=0; j<y_end-y_start+1; j++) 
		{

		    for (int new_state=0; new_state<number_of_states; new_state++)
		    {

			int number_of_old_states = model[new_state].get_number_of_previous_states();
			
			for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
			{
			    int old_state = model[new_state].get_number_of_previous_state(old_state_number);
			    
			    transition_score = Logzero;
			    
			    int x_position = (i + x_midpoint - internal_strip_width-1) - this->model[new_state].get_letters_to_read_x();
			    int y_position = (j + y_start + 1) - this->model[new_state].get_letters_to_read_y();


			    // note: this must be un-mirrored pairhmm
			    
			    transition_score = this->new_get_transition_score(NULL, 
									      old_state, new_state,
									      x, x_position,
									      y, y_position);

			    if (transition_score > Logzero)
			    {
				test_score=Logzero;
				
				test_score=forward_strip[old_state][i][j]+
				    backward_strip[new_state][i][j]+
				    transition_score;
				
				forward_score=forward_strip[old_state][i][j];
				backward_score=backward_strip[new_state][i][j];

				if (test_score > max_score)
				{
				    max_score=test_score;
				    *max_i=i;
				    *max_j=j;
				    *max_state_forward=old_state;
				    *max_state_backward=new_state;
				    
				    max_forward_score=forward_score;
				    max_backward_score=backward_score;
				    max_transition_score=transition_score;

				}
			    }
			} // for old_state
		    } // for new_state
		} // for j (y values)
	    } // for i (x value in strip)

	}
	
	// delete forward_strip if it was created within function
	
	if (number_of_start_states != 0) {

	    for (state=0; state<number_of_states; state++) {
		for (i=0; i<internal_strip_width; i++) {
		    if (forward_strip[state][i])  delete [] forward_strip[state][i];
		    forward_strip[state][i]  = NULL;
		}
		if (forward_strip[state])  delete [] forward_strip[state];
		forward_strip[state]  = NULL;
	    }
	    if (forward_strip)  delete [] forward_strip;
	    forward_strip=NULL;
	}
	
	// delete backward_strip if it was created within function
	
	if (number_of_end_states != 0) {

	    for (state=0; state<number_of_states; state++) {
		for (i=0; i<internal_strip_width; i++) {
		    if (backward_strip[state][i])  delete [] backward_strip[state][i];
		    backward_strip[state][i]  = NULL;
		}
		if (backward_strip[state])  delete [] backward_strip[state];
		backward_strip[state]  = NULL;
	    }
	    if (backward_strip)  delete [] backward_strip;
	    backward_strip=NULL;
	}
    }
    if (check != 0) {
	
	// initialise return values
	
	*max_i=0;
	*max_j=0;
	*max_state_forward=0;
	*max_state_backward=0;
    }
    return(check);
}

int Hmm::get_strip(Score**** const viterbi_strip, 
		       const Hmm *s, 
		       const Hmm *unmirrored_s,
		       const Sequence *x, const int x_start, const int x_end,
		       const Sequence *y, const int y_start, const int y_end, 
		       const int strip_width, // strip_width-1 = min_strip_width
		       const int x_margin,
		       const int y_margin,
		       const int direction)
{
    int check=0;

    // check strip_width
    
    if (strip_width<2)
    {
	cout << "ERROR Hmm::get_strip: strip_width (" << strip_width << ") has to be >= 2.\n" << flush;
	check+=1;
    }
    
    // x_margin has to be 0 for x_start=0
    // y_margin has to be 0 for y_start=0
    
    if (direction == 1) {
	if ((x_start==0) && (x_margin!=0))
	{
	    cout << "ERROR: Hmm::get_strip : direction == 1 : x_margin (" << x_margin 
		 << ") has to be zero for x_start (" << x_start << ") = 0.\n" << flush;
	    check+=1;
	}
    
	if ((y_start==0) && (y_margin!=0))
	{
	    cout << "ERROR: Hmm::get_strip : direction == 1 : y_margin (" << y_margin 
		 << ") has to be zero for y_start (" << y_start << ") = 0.\n" << flush;
	    check+=1;
	}
    }
    else if (direction == -1) {
	if ((x_end==x->length()) && (x_margin!=0))
	{
	    cout << "ERROR: Hmm::get_strip : direction == -1 : x_margin (" << x_margin 
		 << ") has to be zero for x_end (" << x_end << ") = x->length() (" << x->length() << ").\n" << flush;
	    check+=1;
	}
	
	if ((y_end==y->length()) && (y_margin!=0))
	{
	    cout << "ERROR: Hmm::get_strip : direction == -1 : y_margin (" << y_margin 
		 << ") has to be zero for y_end (" << y_end << ") = y->length() (" << y->length() << ").\n" << flush;
	    check+=1;
	}
    }
    
    // check that x_margin and y_margin have sensible values

    if (x_margin<0)
    {
	cout << "ERROR: Hmm::get_strip : x_margin (" << x_margin << ") < 0.\n" << flush;
	check+=1;
    }
    
    if (y_margin<0)
    {
	cout << "ERROR: Hmm::get_strip : y_margin (" << y_margin << ") < 0.\n" << flush;
	check+=1;
    }

    // check x_margin and y_margin and their compatibility with x_start, x_end
    // y_start and y_end
    
    if ((x_end-x_start+1) < x_margin)
    {
	cout << "ERROR: Hmm::get_strip : x_margin (" << x_margin << ") > (x_end (" << x_end
	     << ") - x_start (" << x_start << ") + 1) = " << x_end-x_start+1 << "\n" << flush;
	check+=1;
    }
    
    if ((y_end-y_start+1) < y_margin)
    {
	cout << "ERROR: Hmm::get_strip : y_margin (" << y_margin << ") > (y_end (" << y_end
	     << ") - y_start (" << y_start << ") + 1) = " << y_end-y_start+1 << "\n" << flush;
	check+=1;
    }

    // check viterbi_strip

    if (viterbi_strip == NULL)
    {
	cout << "ERROR: Hmm::get_strip : viterbi_strip is NULL.\n" << flush;
	check+=1;
    }

    // check Hmm
  
    if (s->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::get_strip :  Hmm : number of states = " 
	     << s->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (s->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::get_strip :  Hmm : state 0 of != Start \n" << flush;
	cout << "s->model[0].get_letters_to_read() : "<<s->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (s->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::get_strip :  Hmm : last state of type != End \n" << flush;
	cout<<"s->model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<s->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (int i=0; i<number_of_states; i++)
    {
	if (s->model[i].get_alphabet()!=alphabet) 
	{
	    cout << "ERROR: Hmm::get_strip :  Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << s->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_mirrored()!=s->get_mirrored()) 
	{
	    cout << "ERROR: Hmm::get_strip :  Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << s->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << s->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::get_strip :  Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << s->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}

	if ((i!=0) && (i!=(number_of_states-1))&& (s->model[i].get_letters_to_read()==0))
	{
	    cout << "ERROR: Hmm::get_strip :  Hmm : state i = " 
		 << i << " should emit! \n" << flush;
	    check+=1;
	}
	
    }

    // check x_start, x_end and y_start, y_end values

    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::get_strip: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::get_strip: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::get_strip: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::get_strip: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    // check value of direction
    
    if ((direction!=1) && (direction!=-1))
    {
	cout << "ERROR Hmm::get_strip: direction (" << direction << ") has to be 1 or -1.\n" << flush;
	check+=1;
    }
    
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::get_strip : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::get_strip : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

	const int n_of_states=s->get_number_of_states();

	// variables which will only be needed to calculate alignment
	
	int shift_index = bitshift(alphabet);
	int i,j;
	int xsteps, ysteps;
	int deltax, deltay;
	int x_strip = 0;
	int max_state, max_deltax, max_deltay;
	Score max_score;
	
	int test_state, test_deltax, test_deltay;
	Score test_score;
	Score pre_test_score;

	int states=n_of_states;

	int start_x=0;
	int end_x=0;
	int start_y=0;
	int end_y=0;
	int offset=0;
	
	int y_length=y_end-y_start+1;

	if (direction==1)
	{
	    offset=0;
	    
	    start_x=x_start; 
	    end_x=x_end+1; 
	    start_y=y_start;
	    end_y=y_end+1; 
	}
	else if (direction==-1) // permute start and end points if running backwards
	{
	    offset=1;
	    
	    start_x=x_end-1;
	    end_x=x_start-2;
	    start_y=y_end-1;
	    end_y=y_start-2;
	}
	for ( xsteps=start_x; xsteps != end_x; xsteps+=direction) 
	{	 
	    if ((abs(xsteps - start_x) < strip_width) && (x_margin != 0)) {
		x_strip = abs(xsteps - start_x);
	    }
	    else {
		x_strip = strip_width-1;
	    }

	    if (((x_margin == 0) && (xsteps != start_x)) ||
		((x_margin != 0) && (abs(xsteps - start_x) >= strip_width)))
	    {

		for (int state=0; state<n_of_states-1; state++)
		{
		    Score* pointer=(*viterbi_strip)[state][0];

		    for (int i=0; i<strip_width-1; i++)
		    {
			(*viterbi_strip)[state][i]=(*viterbi_strip)[state][i+1];

		    }
		    (*viterbi_strip)[state][strip_width-1]=pointer;

		    for (int j=0; j<y_length; j++)
		    {	     
			(*viterbi_strip)[state][strip_width-1][j]= Logzero;
		    }

		}

	    }

	    // start of loop over ysteps is dependent on value of xsteps

	    int start_of_ysteps_loop=0;
	    
	    if (abs(xsteps-start_x) < x_margin)
	    {
		start_of_ysteps_loop=start_y+direction*y_margin;

	    }
	    else
	    {
		start_of_ysteps_loop=start_y;

	    }

	    for ( ysteps=start_of_ysteps_loop; ysteps != end_y; ysteps+=direction) 
	    {

		if ((direction==1) && (xsteps==0) && (ysteps==0))
		{
		    (*viterbi_strip)[0][strip_width-1][0-y_start]=0.;

		}
		else if ((direction==-1) && ((x->length()-1-xsteps)==0) && ((y->length()-1-ysteps)==0))	       
		{
		    (*viterbi_strip)[0][strip_width-1][ysteps+offset-y_start]=0.;

		}
		else
		{

#ifdef _CONSTRAINT_CHECK			    
		
		    int   max_old_state      = 0;
		    Score max_emit_score     = 0;
		    Score max_trans_score    = 0;
		    Score max_this_score     = 0;
		    Score max_this_previous_score = 0;
		    Score max_pre_test_score = 0;
		
#endif // #ifdef _CONSTRAINT_CHECK			    

		    for (int new_state=1; new_state<states-1; new_state++)
		    {		     
			max_score=Logzero;
			max_state=0;
			max_deltax=0;
			max_deltay=0;
			
			deltax=0;
			deltay=0;
			
			deltax=s->model[new_state].get_letters_to_read_x();
			deltay=s->model[new_state].get_letters_to_read_y();

			if ( ((direction==1) && 
			      (xsteps-start_x>=deltax) && (ysteps-start_y>=deltay)) ||
			     ((direction==-1) && 
			      ((start_x-xsteps)>=deltax) && ((start_y-ysteps)>=deltay)) )
			{
			    int linear_index_emission = 0;
			    int* indices     = NULL;
			    indices = new int[s->model[new_state].get_number_of_dimensions_of_emission_scores()+1];
			    int* positions   = NULL;
			    positions = new int[s->model[new_state].get_number_of_dimensions_of_emission_scores()+1];
			    
			    int xcount=0;
			    int ycount=0;

			    if (direction==1)
			    {
				xcount=deltax;
				ycount=0;
			    }
			    else if (direction==-1)
			    {
				xcount=0;
				ycount=deltay;
			    }
			    for (i=ycount; i<deltax+ycount; i++)
			    {
				indices[i+1]   = x->letter(xsteps - direction*(deltax+ycount-i));
				positions[i+1] = xsteps - direction*(deltax+ycount-i);

			    }
			    for (j=xcount; j<deltay+xcount; j++)
			    {
				indices[j+1]   = y->letter(ysteps - direction*(deltay+xcount-j));
				positions[j+1] = ysteps - direction*(deltay+xcount-j);

			    }
			    int number_of_old_states = s->model[new_state].get_number_of_previous_states();
			    
			    for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
			    {
				int old_state = s->model[new_state].get_number_of_previous_state(old_state_number);
				pre_test_score=Logzero;
				test_score=Logzero;
				test_state=0;
				test_deltax=0;
				test_deltay=0;

				indices[0] = old_state;
				linear_index_emission = indices[1];
				for (int p=2; p<s->model[new_state].get_number_of_dimensions_of_emission_scores()+1; p++)
				{
				    linear_index_emission = (linear_index_emission << shift_index) | indices[p];
				}

				int x_position = 0;
				int y_position = 0;
				
				if (direction == 1)
				{
				    x_position = xsteps - deltax;
				    y_position = ysteps - deltay;
				}
				else if (direction == -1)
				{
				    x_position = xsteps + 1;
				    y_position = ysteps + 1;
				}

				pre_test_score = s->new_get_transition_score(this, 
									     old_state, new_state, 
									     x, x_position, y, y_position) +
				    s->model[new_state].get_emission_score(x, x_position, y, y_position,
									   linear_index_emission);

#ifdef _CONSTRAINT_CHECK      		    

				Score this_trans = 
				    s->new_get_transition_score(this, 
								old_state, new_state, 
								x, x_position, y, y_position);
				
				Score this_emit = s->model[new_state].get_emission_score(x, x_position, y, y_position,
											 linear_index_emission);
				Score this_score          = pre_test_score;
				Score this_previous_score = (*viterbi_strip)[old_state]
				    [x_strip-deltax]
				    [ysteps-direction*deltay+offset-y_start];


				if ((this_emit == Logzero) && (this_trans == Logzero)) {

				    if (pre_test_score > Logzero) {
					cout << "get_strip : (xsteps, ysteps, old_state, new_state) = ("
					     << xsteps << ", " << ysteps << ", " << old_state << ", " << new_state 
					     << ") : pre_test_score = this_trans (" << this_trans 
					     << ") + this_emit (" << this_emit << ") = " << pre_test_score << "\n";
					cout << "ERROR : get_strip : pre_test_score = " << pre_test_score 
					     << " > Logzero = "
					     << Logzero << " but this_trans (" << this_trans 
					     << ") == Logzero and this_emit (" << this_emit << ") == Logzero.\n";
					exit;
				    }
				}
				
				if (this_score > Logzero) {
				    
				    this_score += this_previous_score;
				}
				
#endif // #ifdef _CONSTRAINT_CHECK			    

				if (pre_test_score > Logzero)
				{
				    test_state=old_state;
				    test_score=pre_test_score+ 
					(*viterbi_strip)[old_state][x_strip-deltax][ysteps-direction*deltay+offset-y_start];
				    test_deltax=deltax;
				    test_deltay=deltay;

				    if (test_score>max_score)
				    {

					max_score=test_score;
					max_state=test_state;
					max_deltax=test_deltax;
					max_deltay=test_deltay;

#ifdef _CONSTRAINT_CHECK			    

					max_old_state   = old_state;
					max_emit_score  = this_emit; 
					max_trans_score = this_trans;
					max_this_score  = this_score;
					max_this_previous_score = this_previous_score;
					max_pre_test_score = pre_test_score;

#endif // #ifdef _CONSTRAINT_CHECK			    
				    }
				}
			    } // for loop old_state
			    
			    if (indices) delete [] indices;
			    indices = NULL;
			    if (positions) delete [] positions;
			    positions = NULL;

			    if ( max_score>Logzero )
			    {
				(*viterbi_strip)[new_state][x_strip][ysteps+offset-y_start]=max_score;
				
#ifdef _CONSTRAINT_CHECK			    
				
				Score reference_score = x->get_sp_score(xsteps, 
									ysteps,
									new_state);
				
				if (this->mirrored == 0) {
				    if ((reference_score != Logzero) && (reference_score != max_score)) {
					
					cout << "ERROR: get_strip : (xsteps, ysteps, max_old_state, new_state) = ("
					     << xsteps << ", " << ysteps << ", " << max_old_state << ", " << new_state 
					     << ") : max_score = " << max_score 
					     << " = max_pre_test_score (" << max_pre_test_score 
					     << ") + max_this_previous_score(" 
					     << max_this_previous_score << ") != reference_score = "
					     << reference_score << "\n";
					exit;
				    }
				    else if ((reference_score != Logzero) && (reference_score == max_score)) {
					
					cout << "get_strip : (xsteps, ysteps, max_old_state, new_state) = ("
					     << xsteps << ", " << ysteps << ", " << max_old_state << ", " << new_state 
					     << ") : max_score = " << max_score 
					     << " = " << max_pre_test_score << " + " 
					     << max_this_previous_score << " == reference_score = "
					     << reference_score << "\n";
				    }
				}
				
#endif // #ifdef _CONSTRAINT_CHECK			    

			    }
			} // if ((xsteps>=deltax) && (ysteps>=deltay))
			else
			{

			}
		    } // for loop new_state
		} // if not ( (xsteps==0) && (ysteps==0) )
	  } // for loop ysteps
	} // for loop xsteps      

	
	if (  ((direction==1)          &&
	       (xsteps==x->length()+1) &&   
	       (ysteps==y->length()+1)) 
	      ||
	      ((direction==-1)         &&
	       (xsteps==-2) &&  
	       (ysteps==-2))  )
	{	

	    if (direction==1)
	    {
		xsteps=x->length();
		ysteps=y->length();
	    }
	    else if (direction==-1)
	    {
		xsteps=-1;
		ysteps=-1;
	    }
	    
	    max_score=Logzero;
	    max_state=0;
	    max_deltax=0;
	    max_deltay=0;
	    
	    deltax=0;
	    deltay=0;
	
	    for (int state=0; state<states-1; state++)  
		// loop over all possible previous states except the End state
	    {
		pre_test_score=Logzero;
		test_score=Logzero;
		test_state=0;
		test_deltax=0;
		test_deltay=0;

		int x_position = 0;
		int y_position = 0;

		if (direction == 1)
		{
		    x_position = xsteps - deltax;
		    y_position = ysteps - deltay;
		}
		else if (direction == -1)
		{
		    x_position = xsteps + 1;
		    y_position = ysteps + 1;
		}
		
		pre_test_score=s->new_get_transition_score(this, 
							   state, states-1,
							   x, x_position, y, y_position);

		if (pre_test_score > Logzero)
		{
		    test_state=state;
		    test_score=pre_test_score+
			(*viterbi_strip)[state][strip_width-1-deltax][ysteps-direction*deltay+offset-y_start];
		    test_deltax=0;
		    test_deltay=0;

		    if (test_score>max_score)
		    {

			max_score=test_score;
			max_state=test_state;
			max_deltax=test_deltax;
			max_deltay=test_deltay;
		    }
		}
	    }
	    if ( max_score>Logzero )
	    {
		(*viterbi_strip)[states-1][strip_width-1][ysteps+offset-y_start]=max_score;

	    }
	}
	
	if (direction==-1)
	{

	    // reverse viterbi_strip : 
	    //
	    //        state -> number_of_states-1-state
	    //        x     -> strip_width-1-x
	    
	    Score*** copy_strip= new Score**[number_of_states];
	    
	    for (int state=0; state<number_of_states; state++)
	    {
		copy_strip[state]= new Score*[strip_width];
		
		for (int i=0; i<strip_width; i++)
		{
		    copy_strip[state][i]= new Score[y_length];
		    
		    for (int j=0; j<y_length; j++)
		    {	     
			copy_strip[state][i][j]= (*viterbi_strip)[state][i][j];
		    }
		}	 
	    }

	    {
		for (int state=0; state<number_of_states; state++)
		{
		    for (int i=0; i<strip_width; i++)
		    {
			for (int j=0; j<y_length; j++)
			{	     
			    (*viterbi_strip)[number_of_states-1-state][strip_width-1-i][j]=copy_strip[state][i][j];
			}
		    }	 
		}
	    }

	    // delete copy_strip

	    if (copy_strip)
	    {
		for (int state=0; state<number_of_states; state++)
		{
		    for (int i=0; i<strip_width; i++)
		    {
			if (copy_strip[state][i]) delete [] copy_strip[state][i];
			copy_strip[state][i] = NULL;
		    }
		    if (copy_strip[state]) delete [] copy_strip[state];
		    copy_strip[state] = NULL;
		}
		if (copy_strip) delete [] copy_strip;      
		copy_strip=NULL;
	    }
	}
    }
    return(check);
}

int Hmm::memory_viterbi(const Hmm* mirror,		    
			    const Sequence *x, const int x_start, const int x_end,
			    const Sequence *y, const int y_start, const int y_end,
			    const int start_state, const int end_state,
			    const int min_strip_width, const int max_area)
{
    int check=0;
    // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::memory_viterbi: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::memory_viterbi: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::memory_viterbi: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::memory_viterbi: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    // check start_state and end_state
    
    if ((x_start==0) && (y_start==0) && (start_state!=0))
    {
	cout << "ERROR Hmm::memory_viterbi: start_state (" << start_state 
	     << ") for x_start=y_start=0 should be 0.\n" << flush;
	check+=1;
    }
    if ((x_start==x->length()) && (y_start==y->length()) && (end_state!=number_of_states-1))
    {
	cout << "ERROR Hmm::memory_viterbi: end_state (" << end_state 
	     << ") for x_start = x->length() (" << x->length() << ") and y_start = y->length() (" 
	     << y->length() << ") should be 0.\n" << flush;
	check+=1;
    }
    if ((start_state<0) || (start_state>number_of_states-1))
    {
	cout << "ERROR Hmm::memory_viterbi: start_state (" << start_state << ") out of range [0, " 
	     << number_of_states-1 << "]\n" << flush;
	check+=1;
    }
    if ((end_state<0) || (end_state>number_of_states-1))
    {
	cout << "ERROR Hmm::memory_viterbi: end_state (" << end_state << ") out of range [0, " 
	     << number_of_states-1 << "]\n" << flush;
	check+=1;
    }
    
 
    // check Hmm
    
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::memory_viterbi : mirrored (" << mirrored 
	     << ") of Hmm must be 0.\n" << flush;
	check+=1;
    }  
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::memory_viterbi : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::memory_viterbi : state 0 != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::memory_viterbi : last state != End \n" << flush;
	cout <<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::memory_viterbi : alphabet of state i = " << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::memory_viterbi : mirrored of state i = " << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::memory_viterbi : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1))&&(model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::memory_viterbi : state i = " << i << " should emit! \n" << flush;
		check+=1;
	    }
	   
	}
    }
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : mirrored (" 
	     << mirror->get_mirrored() << ") must be 1.\n" << flush;
	check+=1;
    }
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : state != Start \n" << flush;
	cout <<"mirror->model[0].get_letters_to_read() : "
	     <<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : last state != End \n" << flush;
	cout <<"mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!=mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }	   
	    if ((i!=0) && (i!=(number_of_states-1)) && (mirror->model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::memory_viterbi : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    // check that this Hmm and the mirror Hmm belong together

    if (number_of_states != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::memory_viterbi : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (int i=1; i<number_of_states-1; i++)
	{
	    if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::memory_viterbi : letters_to_read of model[" << i << "] (" 
		  << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		<< number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    
    // check min_strip_width 
    
    if (min_strip_width<1)
    {
	cout << "ERROR Hmm::memory_viterbi: min_strip_width (" << max_area << ") has to be >= 1.\n" << flush;
      check+=1;
    }
    
    // check max_area
    
    if (max_area<1)
    {
	cout << "ERROR Hmm::memory_viterbi: max_area (" << max_area << ") has to be >= 1.\n" << flush;
	check+=1;
    }

    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::memory_viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }

    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::memory_viterbi : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {


	if ((x_end-x_start+1)*(y_end-y_start+1) < max_area)
	{
	    // use standard viterbi 	  
	    
	    int start_states[1]={start_state};
	    const int number_of_start_states=1;
	    int x_start_traceback     = 0;
	    int y_start_traceback     = 0; 
	    int start_state_traceback = 0;

	    check+=viterbi_rectangle(NULL, // do not use external viterbi_rectangle
				     this,
				     x, x_start, x_end,
				     y, y_start, y_end,
				     start_states, end_state,
				     number_of_start_states,
				     0,0, // x_margin, y_margin
				     0,   // offset for add_local_to_global_solution
				     // output values (won't be used any further in this function)
				     &x_start_traceback, &y_start_traceback,
				     &start_state_traceback);
	}
	else
	{
	    // use memory_viterbi


	    // split area up in two rectangles 
	    
	    int x_midpoint=x_start+static_cast<int>(floor(static_cast<double>(x_end-x_start+min_strip_width+1)/2.));

	    // check whether x_midpoint value can be used
	    
	    int max_i=0;
	    int max_j=0;
	    int max_state_forward=0;
	    int max_state_backward=0;
	    
	    if ((x_midpoint-(min_strip_width+1))>=0)
	    {
		check+=get_strips_and_merge(mirror, 
					    x, x_start, x_end,
					    y, y_start, y_end,
					    start_state, end_state,
					    x_midpoint, min_strip_width,
					    &max_i, &max_j, &max_state_forward, &max_state_backward);
		
		if (check==0)
		{
		    // start new iteration
		    
		    int x_value=x_midpoint-1-(min_strip_width-max_i);
		    int y_value=y_start+max_j;

	    
		    // lower left rectangle
		    
		    check+=memory_viterbi(mirror,		    
					  x, x_start, x_value,
					  y, y_start, y_value,
					  start_state, max_state_forward,
					  min_strip_width, max_area);
		    
		    // upper right rectangle
		    
		    check+=memory_viterbi(mirror,		   
					  x, x_value+model[max_state_backward].get_letters_to_read_x(), x_end,
					  y, y_value+model[max_state_backward].get_letters_to_read_y(), y_end,
					  max_state_backward, end_state,
					  min_strip_width, max_area);
		}
		else
		{
		    cout << "ERROR: Hmm::memory_viterbi : error occured in get_strips_and_merge function.\n" << flush;
		    check+=1;
		}
	    }
	    else
	    {
		cout << "ERROR: Hmm::memory_viterbi : x_midpoint (" << x_midpoint << ") - (min_strip_width("
		     << min_strip_width << ") + 1 ) = " << x_midpoint-(min_strip_width+1) 
		     << " < 0. Cannot continue calculation.\n" << flush;
		check+=1;
	    }
	}
    }
    return(check);
}

int Hmm::new_memory_viterbi(// input values
				Score*** viterbi_strip, 
				const Hmm* mirror,
				const Sequence *x, const int x_start, const int x_end,
				const Sequence *y, const int y_start, const int y_end, 
				const int* start_states, const int end_state,
				const int number_of_start_states,
				const int x_margin, const int y_margin,
				const int offset,
				const int min_strip_width,
				const int max_area,
				// output values
				int* x_start_traceback,
				int* y_start_traceback,
				int* start_state_traceback)
{
    // function returns 0 if o.k., else not o.k.
    //
    // note: this function can either be used 
    //
    //       1.) with start_states and an end_state (in which case
    //       viterbi_strip == NULL and x_margin and y_margin == 0). Memory
    //       for viterbi_strip is allocated and deleted within function.
    //       
    //       2.) with a pre-initialised viterbi_strip (in which case start_states
    //        == NULL and end_state == 0 and number_of start_states and number_of_end_states == 0).
    //       Memory for viterbi_strip is allocated and deleted outside function.
    
    int check=0;
    
    const int n_of_states= this->get_number_of_states();
    
    // check whether pre-initialised viterbi_rectangle or start_ and end_states shall be used
    
    if ((viterbi_strip == NULL)    &&
	(number_of_start_states == 0))
    {
	cout << "ERROR Hmm::new_memory_viterbi: viterbi_strip (" << viterbi_strip
	     << ") = NULL and number_of_start_states (" << number_of_start_states 
	     << ") == 0. Use this function either with a pre-initialised viterbi_strip != NULL "
	     << "or with nonempty sets of start_ and end_states.\n" << flush;
	check+=1;
    }
    else if ((viterbi_strip != NULL)    &&
	     (number_of_start_states != 0))
    {
	cout << "ERROR Hmm::new_memory_viterbi: viterbi_strip (" << viterbi_strip
	     << ") != NULL and number_of_start_states (" << number_of_start_states 
	     << ") != 0. Use this function either with a pre-initialised viterbi_strip != NULL "
	     << "or with nonempty sets of start_ and end_states.\n" << flush;
	check+=1;
    }
    else // if input values are not completely inconsistent
    {
	if ((viterbi_strip      == NULL)   &&
	    (number_of_start_states != 0))
	{
	    // if function shall be used with start_ and end_states

	    // check start states
	    
	    if (number_of_start_states<1)
	    {
		cout << "ERROR Hmm::new_memory_viterbi: number_of_start_sites (" 
		     << number_of_start_states << ") <1.\n" << flush;
		check+=1;
	    }
	    else
	    {
		if ((x_start==0) && (y_start==0))
		{
		    if ((number_of_start_states!=1) || (start_states[0]!=0))
		    {
			cout << "ERROR Hmm::new_memory_viterbi: or "
			     << "x_start=y_start=0 only "
			     << "one start_state is possible, Start state (i.e. the number_of_start_states ("
			     << number_of_start_states << ") has to be 1 and start_states[0] ("
			     << start_states[0] << ") has to be 0).\n" << flush;
			check+=1;
		    }
		}
		else // if we don't start at (0,0)
		{
		    // check that each start state lies in the allowed range (1...n_of_states-2)
		    // i.e.  Start and end states are not allowed 
		    
		    for (int i=0; i<number_of_start_states; i++)
		    {
			if ((start_states[i]<1) || (start_states[i]>n_of_states-2))
			{
			    cout << "ERROR Hmm::new_memory_viterbi: start_states[" << i 
				 << "] (" << start_states[i] << ") out of range [1, " << n_of_states-2
				 << "].\n" << flush;
			    check+=1;
			}
		    }
		}
	    }
	    
	    // check end state
	    
	    if ((x_end==x->length()) && (y_end==y->length()))
	    {
		if (end_state != n_of_states-1)
		{
		    cout << "ERROR Hmm::new_memory_viterbi: end_state (" 
			 << end_state << ") for x_end (" << x_end << ") = x->length() (" << x->length() 
			 << ") and y_end (" << y_end << ") = y->length() (" << y->length() 
			 << ") has to be End state (" << n_of_states-1 << ").\n" << flush;
		    check+=1;
		}
	    }
	    else // if we don't stop at ends of both sequences
	    {
		// check that end state lies in the allowed range (1...n_of_states-2)
		// i.e. start and End states are not allowed
		
		if ((end_state<1) || (end_state>n_of_states-2))
		{
		    cout << "ERROR Hmm::new_memory_viterbi: end_state ("
			 << end_state << ") out of range [1, " << n_of_states-2 << "].\n" << flush;
		    check+=1;
		}
	    }
	}
	
	if ((viterbi_strip      != NULL)   &&
	    (number_of_start_states == 0))
	{
	    // if function shall be used with pre-initialised viterbi_strip
	    
	    // check values of x_margin and y_margin

	    if ((x_start==0) && (x_margin!=0))
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : x_margin (" << x_margin 
		     << ") has to be zero for x_start (" 
		     << x_start << ") = 0.\n" << flush;
		check+=1;
	    }
	    
	    if ((y_start==0) && (y_margin!=0))
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : y_margin (" << y_margin 
		     << ") has to be zero for y_start (" 
		     << y_start << ") = 0.\n" << flush;
		check+=1;
	    }
	    
	    if (x_margin < 0)
	    {
		cout << "ERROR Hmm::new_memory_viterbi: x_margin (" << x_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	    
	    if (y_margin < 0)
	    {
		cout << "ERROR Hmm::new_memory_viterbi: y_margin (" << y_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	    
	    // check that x_margin and y_margin are compatible with values of x_start, x_end and y_start and y_end
	    
	    if ((x_end-x_start+1) < x_margin)
	    {
		cout << "ERROR Hmm::new_memory_viterbi: x_margin (" << x_margin
		     << ") > (x_end (" << x_end << ") - x_start (" << x_start << ") + 1) = " 
		     << x_end-x_start+1 << "\n" << flush;
		check+=1;
	    }
	    
	    if ((y_end-y_start+1) < y_margin)
	    {
		cout << "ERROR Hmm::new_memory_viterbi: y_margin (" << y_margin
		     << ") > (y_end (" << y_end << ") - y_start (" << y_start << ") + 1) = " 
		     << y_end-y_start+1 << "\n" << flush;
		check+=1;
	    }
	}
    }

    // check this Hmm
    
    if (this->get_mirrored() != 0)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : Hmm : mirrored (" 
	     << this->get_mirrored() << ") must be 0.\n" << flush;
	check+=1;}
    if (this->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : Hmm : number of states = " 
	     << this->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (this->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : Hmm : state 0 of type != Start \n" << flush;
	cout <<"this->model[0].get_letters_to_read() : "<<this->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (this->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : Hmm : last state of type != End \n" << flush;
	check+=1;
    }  
    for (int i=0; i<number_of_states; i++)
    {
	if (this->model[i].get_alphabet()!=alphabet) 
	{
	    cout << "ERROR: Hmm::new_memory_viterbi : Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << this->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (this->model[i].get_mirrored()!=this->get_mirrored()) 
	{
	    cout << "ERROR: Hmm::new_memory_viterbi : Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << this->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << this->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (this->model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::new_memory_viterbi : Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << this->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}
	if ((i!=0) && (i!=(number_of_states-1)) && (this->model[i].get_letters_to_read()==0))
	{
	    cout << "ERROR: Hmm::new_memory_viterbi : Hmm : state i = " 
		 << i << " should emit! \n" << flush;
	    check+=1;
	}
	
    }

    // check that this Hmm and the mirror Hmm belong together
    
    if (this->get_number_of_states() != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::new_memory_viterbi : number_of_states of this Hmm (" 
	     << this->get_number_of_states()
	     << ") != number_of_states of mirrored Hmm (" 
	     << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (int i=1; i<number_of_states-1; i++)
	{
	    if (this->model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : letters_to_read of model[" << i 
		     << "] in this Hmm (" 
		  << this->model[i].get_letters_to_read() << ") != state type of mirrored Hmm model["
		     << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : mirrored (" << mirror->get_mirrored()
	     << ") of mirrored Hmm must be 1.\n" << flush;
	check+=1;
    }
    
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : state 0 != Start \n" << flush;
	cout <<"mirror->model[0].get_letters_to_read() : "<<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : last state of type != End \n" << flush;
	cout<< "mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=this->model[i].get_alphabet()) 
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		  << mirror->model[i].get_alphabet() << " != alphabet of this Hmm = " 
		     << this->model[i].get_alphabet() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!= mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != this->model[i].get_number_of_states())
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of this Hmm = " 
		     << this->model[i].get_number_of_states() << "\n" << flush;
		check+=1;
	    }
	   
	    if ((i!=0) && (i!=(number_of_states-1))&& (mirror->model[i].get_letters_to_read()))
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    
    // check type of Scoremodel requested
    
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::new_memory_viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::new_memory_viterbi : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }
    
  // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::new_memory_viterbi: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::new_memory_viterbi: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::new_memory_viterbi: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::new_memory_viterbi: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }

    if ((offset != 0) && (offset != -1))
    {
	cout << "ERROR Hmm::new_memory_viterbi: offset (" << offset << ") must be either 0 or -1.\n" << flush;
	check+=1;}
    
    // check min_strip_width 
    
    if (min_strip_width<1)
    {
	cout << "ERROR Hmm::new_memory_viterbi: min_strip_width (" << max_area << ") has to be >= 1.\n" << flush;
	check+=1;
    }
    
    // check max_area
    
    if (max_area<1)
    {
	cout << "ERROR Hmm::new_memory_viterbi: max_area (" << max_area << ") has to be >= 1.\n" << flush;
	check+=1;
    }

    if (check==0)
    {

	// initialise values of output variables
	// in the call of this function the output variables below are either filled by a call to viterbi_rectangle
	// (if this is the final iteration) or by a call to new_memory_viterbi (if there is a next iteration)
	
	(*x_start_traceback)     = 0;
	(*y_start_traceback)     = 0; 
	(*start_state_traceback) = 0;      
	
	int unused_x_start_traceback     = 0;      
	int unused_y_start_traceback     = 0;      
	int unused_start_state_traceback = 0;      
	
	int internal_x_margin = 0;
	int internal_y_margin = 0;
	int internal_max_area = max_area;
	
    TEST_RECTANGLE_AGAIN:
	
	if ((x_end-x_start+1)*(y_end-y_start+1) < internal_max_area)
	{
	    // use standard viterbi 	  

	    // if this rectangle is to be initialised with a strip 
	    
	    if (number_of_start_states == 0) {
		
		internal_x_margin = x_margin;
		internal_y_margin = y_margin;
	    }
	    else {
		
		internal_x_margin = 0;
		internal_y_margin = 0;
	    }
	    
	    Score*** internal_viterbi_rectangle_array = NULL;
	    
	    if (number_of_start_states == 0) {	  

		check += allocate_memory_for_viterbi_rectangle(this,
							       x_start, x_end,
							       y_start, y_end,
							       1, // = direction 
							       &internal_viterbi_rectangle_array);
		
		if ((check != 0) || (internal_viterbi_rectangle_array == NULL)) {
		    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
			 << "Hmm::allocate_memory_for_viterbi_rectangle.\n";
		}
		if (check == 0) {	  

		    check += copy_rectangle_from_strip_to_next_strip(// source strip:
			&viterbi_strip,
			// coordinates of strip
			x_start, x_start+x_margin-1,
			y_start, y_end,
			// coordinates to be copied
			x_start, x_start+x_margin-1,
			y_start, y_start+y_margin-1,
			// target strip:
			&internal_viterbi_rectangle_array,
			// coordinates of strip
			x_start, x_end,
			y_start, y_end,
			// coordinates to be copied to
			x_start, x_start+x_margin-1,
			y_start, y_start+y_margin-1,
			// direction of two strips
			1);

		    if (check != 0) {
			cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
			     << "Hmm::copy_rectangle_from_strip_to_next_strip.\n" << flush;
		    }
		} // if check == 0
	    } // if number_of_start_states == 0
	    
	    if (check == 0) {

		check+=viterbi_rectangle(internal_viterbi_rectangle_array,  // may be NULL
					 this,
					 x, x_start, x_end,
					 y, y_start, y_end,
					 start_states, end_state,   
					 number_of_start_states,   // number_of_start_states may be 0
					 internal_x_margin, 
					 internal_y_margin, 
					 offset,                   // offset for add_local_to_global_solution
				     // output values
					 x_start_traceback, 
					 y_start_traceback,
					 start_state_traceback);
		if (check != 0) {
		    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
			 << "Hmm::viterbi_rectangle.\n" << flush;
		}
	    }
	    if (internal_viterbi_rectangle_array != NULL) {	  

		check += delete_memory_for_viterbi_rectangle(this,
							     x_start, x_end,
							     1,
							     &internal_viterbi_rectangle_array);
		
		if ((check != 0) || (internal_viterbi_rectangle_array != NULL)) {
		    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
			 << "Hmm::delete_memory_for_viterbi_rectangle.\n" << flush;
		}
	    }
	}
	else
	{
	    // use memory_viterbi

	    // split area up in two rectangles 

	    int x_midpoint=x_start+static_cast<int>(floor(static_cast<double>(x_end-x_start+min_strip_width+1)/2.));

	    // check if x_midpoint value can be used
	    
	    int internal_min_strip_width = min_strip_width;
	    int max_i=0;
	    int max_j=0;
	    int max_state_forward=0;
	    int max_state_backward=0;

	    const int end_states[1]        = {end_state};
	    const int number_of_end_states = 1;
	    
	    if ((x_midpoint-(min_strip_width+1))>=0)
	    {
		// prepare calculation of lower left rectangle by
		// copying viterbi_strip into viterbi_strip_copy
		// before is used in function get_strips_and_merge
	      
		Score*** viterbi_strip_copy     = NULL;
		Score*** internal_viterbi_strip = NULL;
		
		if (number_of_start_states == 0) {	  

		    cout << "3 allocate_memory_for_strip:\n" << flush;
		    cout << " x " << x_margin << " y_1 " << y_start << " y_2 " << y_end << "strip " 
			 << viterbi_strip_copy << "\n" << flush;
		    
		    check += allocate_memory_for_strip(this,
						       x_margin,
						   y_start, y_end,
						       &viterbi_strip_copy);
		    
		    if ((check != 0) || (viterbi_strip_copy == NULL)) {
			cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
			     << "Hmm::allocate_memory_for_strip.\n" << flush;
		    }
		    
		    if (check == 0) {

			check += copy_rectangle_from_strip_to_next_strip(// source strip:
			    &viterbi_strip,
			    // coordinates of strip
			    x_start, x_start+x_margin-1,
			    y_start, y_end,
			    // coordinates to be copied
			    x_start, x_start+x_margin-1,
			    y_start, y_end,
			    // target strip:
			    &viterbi_strip_copy,
			    // coordinates of strip
			    x_start, x_start+x_margin-1,
			    y_start, y_end,
			    // coordinates to be copied to
			    x_start, x_start+x_margin-1,
			    y_start, y_end,
			    // direction of two strips
			    1);

			if (check != 0) {
			    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
				 << "Hmm::copy_rectangle_from_strip_to_next_strip.\n" << flush;
			}
		    }
		} // if number_of_start_states == 0
		
		if (check == 0) {
		    
		    check+=new_get_strips_and_merge(// input
			viterbi_strip_copy, 
			// = forward_viterbi_strip (may be NULL)
			NULL,                    
			// = backward_viterbi_strip
			mirror, 
			x, x_start, x_end,
			y, y_start, y_end,
						start_states, number_of_start_states,
			end_states, number_of_end_states,
			x_midpoint, 					      
			x_margin, y_margin,
			// in- and output
			&internal_min_strip_width,
			// output
			&max_i, &max_j, &max_state_forward, &max_state_backward);
		    
		    if (check != 0) {
			cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
			     << "Hmm::new_get_strips_and_merge.\n" << flush;
		    }
		} // if check == 0
		
		if (viterbi_strip_copy != NULL) {

		    check += delete_memory_for_strip(this,
						     x_margin,
						     &viterbi_strip_copy);
		    
		    if ((check != 0) || (viterbi_strip_copy != NULL)) {
			cout << "ERROR: Hmm::new_memory_viterbi: message 1 : error occurred in function "
			     << "Hmm::delete_memory_for_strip.\n" << flush;
		    }
	      }
		
		if (check==0)
		{
		    // start new iteration
		    
		    int x_value=x_midpoint-1-(internal_min_strip_width-max_i);
		    int y_value=y_start+max_j;

		    // check if two new rectangles have some minimum size, if not
		    // merge them into one and adjust new value of internal_max_area
		    // accordingly

		    if (((x_value - x_start + 1) < max(min_strip_width, x_margin)) ||
			((y_value - y_start + 1) < max(min_strip_width, y_margin)) ||
			((x_end - x_value+model[max_state_backward].get_letters_to_read_x() + 1) < 
			 max(min_strip_width, x_margin)) ||
			((y_end - y_value+model[max_state_backward].get_letters_to_read_y() + 1) < 
			 max(min_strip_width, y_margin))) {

			internal_max_area = (x_end - x_start + 1) * (y_end - y_start + 1) + 1;
			goto TEST_RECTANGLE_AGAIN;
		    }
		    
		    if (number_of_start_states == 0) {

			cout << "4 allocate_memory_for_strip:\n" << flush;
			cout << " x " << x_margin << " y_1 " << y_start << " y_2 " << y_value << "strip " 
			     << internal_viterbi_strip << "\n" << flush;
			
			
			check += allocate_memory_for_strip(this,
							   x_margin,
							   y_start, y_value,
							   &internal_viterbi_strip);
		  
			cout << "after 4 allocate_memory_for_strip:\n" << flush;
			
			if ((check != 0) || (internal_viterbi_strip == NULL)) {
			    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
				 << "Hmm::allocate_memory_for_strip.\n" << flush;
			}
			
			if (check == 0) {

			    check += copy_rectangle_from_strip_to_next_strip(// source strip:
				&viterbi_strip,
				// coordinates of strip
				x_start, x_start+x_margin-1,
				y_start, y_end,
				// coordinates to be copied
				x_start, x_start+x_margin-1,
				y_start, y_start+y_margin-1,
				// target strip:
				&internal_viterbi_strip,
				// coordinates of strip
				x_start, x_value,
				y_start, y_value,
				// coordinates to be copied to
				x_start, x_start+x_margin-1,
				y_start, y_start+y_margin-1,
				// direction of two strips
				1);

			    if (check != 0) {
				cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
				     << "Hmm::copy_rectangle_from_strip_to_next_strip.\n" << flush;
			    }
			} // if check == 0
		    } // if number_of_start_states == 0
		    
		    if (check == 0) {
			
			// lower left rectangle

			check += new_memory_viterbi(internal_viterbi_strip, // may be NULL 
						    mirror,
						    x, x_start, x_value,
						    y, y_start, y_value,
						    start_states, max_state_forward,
						    number_of_start_states,
						    x_margin, y_margin,
						    0, // = offset (no overlap with upper right rectangle)
						    min_strip_width,
						    internal_max_area,
						    // output values
						    x_start_traceback, 
						    y_start_traceback,
						    start_state_traceback);


			if (check != 0) {
			    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
				 << "Hmm::new_memory_viterbi for lower left rectangle.\n" << flush;
			}
		    } // if check == 0
		    
		    if (check == 0) {

		    // upper right rectangle

			const int backward_start_states[1] = {max_state_backward};
			
			check += new_memory_viterbi(NULL, // = viterbi_strip
						    mirror,
						    x, x_value+model[max_state_backward].get_letters_to_read_x(), x_end,
						    y, y_value+model[max_state_backward].get_letters_to_read_y(), y_end,
						    backward_start_states, end_state,
						    1, // = number_of_start_states
						    x_margin, y_margin,
						    offset,
						    min_strip_width,
						    internal_max_area,
						    // output values
						    &unused_x_start_traceback, 
						    &unused_y_start_traceback,
						    &unused_start_state_traceback);


			if (check != 0) {
			    cout << "ERROR: Hmm::new_memory_viterbi: error occurred in function "
				 << "Hmm::new_memory_viterbi for upper right rectangle.\n" << flush;
			}
		    } // if check == 0
		}
		if (internal_viterbi_strip != NULL) {

		    check += delete_memory_for_strip(this,
						     x_margin,
						     &internal_viterbi_strip);
		    
		    if ((check != 0) || (internal_viterbi_strip != NULL)) {
			cout << "ERROR: Hmm::new_memory_viterbi: message 2 error occurred in function "
			     << "Hmm::delete_memory_for_strip.\n" << flush;
		    }
		}
	    }
	    else
	    {
		cout << "ERROR: Hmm::new_memory_viterbi : x_midpoint (" << x_midpoint << ") - (min_strip_width("
		     << min_strip_width << ") + 1 ) = " << x_midpoint-(min_strip_width+1) 
		     << " < 0. Cannot continue calculation.\n" << flush;
		check+=1;
	    }
	}
    }
    if (check != 0) {
	
	// set output variables to default values
	
	(*x_start_traceback)     = 0;
	(*y_start_traceback)     = 0; 
	(*start_state_traceback) = 0;      
    }
    return(check);
}

int Hmm::add_local_solution(const Sequence *x, 
				const Sequence *y, 
				const int& local_steps,
				const Score& local_score,
				const Score* const local_sequence_of_scores,
				const int* const local_sequence_of_states,			 
				const int* const local_sequence_of_xsteps,
				const int* const local_sequence_of_ysteps,
				const int direction)
{
    // NOTE: - local solution may not have an overlap with existing solution
    //       - local_score is not used, the elements of local_sequence_of_scores 
    //         are used to calculate the new score
    //       - make sure that this function 
    //         uses only local_sequence information which was derived from a matrix
    //         calculated by an un-mirrored pairhmm
    
    int check=0;
    int max_index=x->length()+y->length()-1;
    
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::add_local_solution : mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }
    

    // check sequences

    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::add_local_solution : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::add_local_solution : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }

    // check x_start, x_end and y_start, y_end values
    
    if (sequence_of_xsteps == NULL)
    {
	cout << "ERROR Hmm::add_local_solution: sequence_of_xsteps is NULL. First allocate memory.\n" << flush;
	check+=1;
    }
    else
    {
	if ((local_sequence_of_xsteps[local_steps] > x->length()) || 
	    (local_sequence_of_xsteps[local_steps] < local_sequence_of_xsteps[0]))
	{
	    cout << "ERROR Hmm::add_local_solution: local_sequence_of_xsteps[" << local_steps << "] (" 
		 << local_sequence_of_xsteps[local_steps] 
		 << ") out of range, may take values in [local_sequence_of_xsteps[0] (" 
		 << local_sequence_of_xsteps[0]
		 << "), x->length() (" << x->length() << ")].\n" << flush;
	    check+=1;
	}
	if (local_sequence_of_xsteps[0] <0)
	{
	    cout << "ERROR Hmm::add_local_solution: local_sequence_of_xsteps[0] (" << local_sequence_of_xsteps[0] 
		 << ") out of range (must be >0).\n" << flush;
	    check+=1;
	}
    }
    
    if (sequence_of_ysteps == NULL)
    {
	cout << "ERROR Hmm::add_local_solution: sequence_of_ysteps is NULL. First allocate memory.\n" << flush;
	check+=1;
    }
    else
    {
	if ((local_sequence_of_ysteps[local_steps] > y->length()) || 
	    (local_sequence_of_ysteps[local_steps] < local_sequence_of_ysteps[0]))
	{
	    cout << "ERROR Hmm::add_local_solution: local_sequence_of_ysteps[" << local_steps << "] (" 
		 << local_sequence_of_ysteps[local_steps] 
		 << ") out of range, may take values in [local_sequence_of_ysteps[0] (" 
		 << local_sequence_of_ysteps[0]
		 << "), y->length() (" << y->length() << ")].\n" << flush;
	    check+=1;
	}
	if (local_sequence_of_ysteps[0]<0)
	{
	    cout << "ERROR Hmm::add_local_solution: local_sequence_of_ysteps[0] (" << local_sequence_of_ysteps[0] 
	       << ") out of range (must be >0).\n" << flush;
	    check+=1;
	}
    }
    
    // check value of direction
    
    if (direction!=1)
    {
	cout << "ERROR Hmm::add_local_solution: direction (" << direction << ") must be 1.\n" << flush;
	check+=1;
    }
    
    // check value of local_steps

    if (local_steps<0)
    {
	cout << "ERROR Hmm::add_local_solution: local_steps (" << local_steps << ") must be >=0.\n" << flush;
	check+=1;
    }
    
    // check value of local_score
    
    if (local_score==Logzero)
    {
	cout << "ERROR Hmm::add_local_solution: local_score (" << local_score << ") is Logzero.\n" << flush;
	check+=1;
    }
    
    // check local sequences
    
    if (local_sequence_of_scores == NULL)
    {
	cout << "ERROR Hmm::add_local_solution: local_sequence_of_scores is NULL\n" << flush;
	check+=1;
    }
    if (local_sequence_of_states == NULL)
    {
	cout << "ERROR Hmm::add_local_solution: local_sequence_of_states is NULL\n" << flush;
	check+=1;
    }
    if (sequence_of_scores == NULL)
    {
	cout << "ERROR Hmm::add_local_solution: sequence_of_scores is NULL. First allocate memory.\n" << flush;
	check+=1;
    }
    if (sequence_of_states == NULL)
    {
	cout << "ERROR Hmm::add_local_solution: sequence_of_states is NULL. First allocate memory.\n" << flush;
	check+=1;
    }
    
    // check that new sequence has no overlap with existing sequence
    
    for (int i=0; i<max_index-1; i++)
    {
	if ((local_sequence_of_xsteps[0]-sequence_of_xsteps[i] == 0) &&
	    (local_sequence_of_ysteps[0]-sequence_of_ysteps[i] == 0) && 
	    (local_sequence_of_states[0] == sequence_of_states[i]))
	{
	    cout << "ERROR Hmm::add_local_solution: cannot implement local_sequence if it has overlap"
		 << " with existing sequence:\n" << flush;
	    cout << "local_sequence_of_xsteps[0] (" << local_sequence_of_xsteps[0] <<") == sequence_of_xsteps[" 
		 << i << "] (" << sequence_of_xsteps[i] 
		 << ") && \n" << flush;
	    cout << "local_sequence_of_ysteps[0] (" << local_sequence_of_ysteps[0] <<") == sequence_of_ysteps[" 
		 << i << "] (" << sequence_of_ysteps[i] 
		 << ") && \n" << flush;
	    cout << "local_sequence_of_states[0] (" << local_sequence_of_states[0] << ") == sequence_of_states[" 
		 << i << "] (" << sequence_of_states[i] << ")\n" << flush;
	    check+=1;
	}
	else if ((local_sequence_of_xsteps[local_steps]-sequence_of_xsteps[i] == 0) &&
		 (local_sequence_of_ysteps[local_steps]-sequence_of_ysteps[i] == 0) && 
		 (local_sequence_of_states[local_steps] == sequence_of_states[i]))
	{
	    cout << "ERROR Hmm::add_local_solution: cannot implement local_sequence if it has overlap"
		 << " with existing sequence:\n" << flush;
	    cout << "local_sequence_of_xsteps[" << local_steps << "] (" << local_sequence_of_xsteps[local_steps] 
		 <<") == sequence_of_xsteps[" << i << "] (" << sequence_of_xsteps[i] 
		 << ") && \n" << flush;
	    cout << "local_sequence_of_ysteps[" << local_steps << "] (" << local_sequence_of_ysteps[local_steps]
		 <<") == sequence_of_ysteps[" << i << "] (" << sequence_of_ysteps[i] 
		 << ") && \n" << flush;
	    cout << "local_sequence_of_states[" << local_steps << "] (" 
		 << local_sequence_of_states[local_steps] << ") == sequence_of_states[" 
		 << i << "] (" << sequence_of_states[i] << ")\n" << flush;
	    check+=1;
	}
    }
    
    if (check==0)
    {

	int shift_index = bitshift(alphabet);

	// look for place where new local solution has to be implemented into the global solution

	int insert_i=0;
	int move=0;
	int last_entry=0;
	
	int prefix=0;
	int append_or_insert=0;

      // if nothing has been added so far, add new local solution starting at [0]

	if ((sequence_of_xsteps[0]==-1) && (sequence_of_ysteps[0]==-1)) 
	{
	    insert_i=0;
	    move=0;
	    last_entry=0;
	}
	else 
	{
	    for (int i=0; i<max_index-1; i++)
	    {
		// if local solution fits in at end of position [i] of existing solution, insert local solution at [i+1] 
		
		if ((local_sequence_of_xsteps[0]-sequence_of_xsteps[i] >= model[local_sequence_of_states[0]].get_letters_to_read_x()) &&
		    (local_sequence_of_ysteps[0]-sequence_of_ysteps[i] >= model[local_sequence_of_states[0]].get_letters_to_read_y()))
		{
		    // if existing solution ends at position [i]
		    
		    if ((sequence_of_xsteps[i]!=-1) && (sequence_of_xsteps[i+1]==-1) &&
			(sequence_of_ysteps[i]!=-1) && (sequence_of_ysteps[i+1]==-1))
		    {
			insert_i=i+1;
			append_or_insert=1;

		    }
		    // if exising solution does not end at position [i], but local solution fits in between
		    // existing solution at position [i] and [i+1], move existing solution from position [i+1]
		    // on and then insert local solution at position [i+1]
		    
		    else if ( ((sequence_of_xsteps[i]!=-1) && (sequence_of_xsteps[i+1]!=-1) &&
			       (sequence_of_ysteps[i]!=-1) && (sequence_of_ysteps[i+1]!=-1))
			      &&			  
			      (((sequence_of_xsteps[i+1]-local_sequence_of_xsteps[local_steps]) 
				>= model[sequence_of_states[i+1]].get_letters_to_read_x()) &&
			       ((sequence_of_ysteps[i+1]-local_sequence_of_ysteps[local_steps]) 
				>= model[sequence_of_states[i+1]].get_letters_to_read_y())) )
		    {
			insert_i=i+1;
			move=1;
			append_or_insert=1;

		    }
		}	      
		
		// determine position of last entry in existing solution
		
		if ((sequence_of_xsteps[i]!=-1)   && (sequence_of_ysteps[i]!=-1) &&
		    (sequence_of_xsteps[i+1]==-1) && (sequence_of_ysteps[i+1]==-1))
		{
		    last_entry=i;

		    break;
		}	      
	    } // for loop over i 
	    
	    // if local solution fits in at start of existing solution
	    
	    if ((sequence_of_xsteps[0]-local_sequence_of_xsteps[local_steps] >= model[sequence_of_states[0]].get_letters_to_read_x()) &&
		(sequence_of_ysteps[0]-local_sequence_of_ysteps[local_steps] >= model[sequence_of_states[0]].get_letters_to_read_y()))
	    {
		insert_i=0;
		move=1;		 
		prefix=1;

	    }
	}
	
	// determine final number of steps

	if (last_entry==0)
	{
	    steps=last_entry+local_steps;

	}
	else if (last_entry>0)
	{
	    steps=last_entry+local_steps+1;

	}

	// if there is space to implement the local solution
	
	if (((steps-1) < max_index) && (check==0))
	{
	    // move entries where new local solution has to be inserted, if necessary
	    
	    if (move==1)
	    {

		if ( (last_entry+local_steps) < max_index ) // if there is space to implement the local solution
		{
		    // shift all entries with indices from insert_i to last_entry incl. local_steps+1 forward

		    for (int i=last_entry; i>insert_i-1; i--)
		    {
			sequence_of_scores[i+local_steps+1]=sequence_of_scores[i];
			sequence_of_states[i+local_steps+1]=sequence_of_states[i];
			sequence_of_xsteps[i+local_steps+1]=sequence_of_xsteps[i];
			sequence_of_ysteps[i+local_steps+1]=sequence_of_ysteps[i];

		    }
		}
		else
		{
		    cout << "ERROR Hmm::add_local_solution: no space left to implement local solution."
			 << "(last_entry (" << last_entry << ") + local_steps (" << local_steps 
			 << ") >= max_index (" << max_index << ")).\n" << flush;
		    check+=1;
		}
	    }
	    
	    if (check==0)
	    {
		// put local solution in entries with indices from insert_i to insert_i+local_steps+1 inclusive

		for (int i=0; i<local_steps+1; i++)
		{
		    score+=local_sequence_of_scores[i];
		    
		    sequence_of_scores[insert_i+i]=local_sequence_of_scores[i];
		    sequence_of_states[insert_i+i]=local_sequence_of_states[i];
		    sequence_of_xsteps[insert_i+i]=local_sequence_of_xsteps[i]; 
		    sequence_of_ysteps[insert_i+i]=local_sequence_of_ysteps[i]; 

		}
	  
	      // close gaps to the left and right, if possible
	  
		Score gap_score=Logzero;
		
		if ( (insert_i>0)
		     &&
		     (sequence_of_xsteps[insert_i-1]!=-1) 
		     &&
		     ((sequence_of_xsteps[insert_i-1]+model[sequence_of_states[insert_i]].get_letters_to_read_x())
		      == sequence_of_xsteps[insert_i]) 
		     &&
		     (sequence_of_ysteps[insert_i-1]!=-1) 
		     &&
		     ((sequence_of_ysteps[insert_i-1]+model[sequence_of_states[insert_i]].get_letters_to_read_y())
		      == sequence_of_ysteps[insert_i]) )
		{

		  // read letters of sequence into index

		    gap_score=Logzero;
		    
		    int* indices     = NULL;
		    indices = new int[model[sequence_of_states[insert_i]].get_number_of_dimensions_of_emission_scores()+1];
		    
		    int deltax=model[sequence_of_states[insert_i]].get_letters_to_read_x();
		    int deltay=model[sequence_of_states[insert_i]].get_letters_to_read_y();
		    
		    int xsteps=sequence_of_xsteps[insert_i];
		    int ysteps=sequence_of_ysteps[insert_i];

		    for (int i=0; i<deltax; i++)
		    {
			indices[i+1]=x->letter(xsteps - (deltax-i));

		    }
		    for (int j=deltax; j<deltay+deltax; j++)
		    {
			indices[j+1]=y->letter(ysteps - (deltay+deltax-j));

		    }
		    
		    indices[0]=sequence_of_states[insert_i-1];
		    
		    int new_state = sequence_of_states[insert_i];
		    int old_state = sequence_of_states[insert_i-1];
		    
		    int x_position = xsteps - deltax;
		    int y_position = ysteps - deltay;

		    if ((deltax+deltay) > 0)
		    {
			int linear_index_emission = indices[1];
			for (int p=2; p<model[sequence_of_states[insert_i]].get_number_of_dimensions_of_emission_scores()+1; p++)
			{
			    linear_index_emission = (linear_index_emission << shift_index) | indices[p];
			}
			
			gap_score=this->new_get_transition_score(NULL,
								 old_state, new_state, 
								 x, x_position, y, y_position) +
			    this->model[new_state].get_emission_score(x, x_position, y, y_position,
								      linear_index_emission);
		    }
		    else
		    {
			gap_score=this->new_get_transition_score(NULL, 
								 old_state, new_state, 
								 x, x_position, y, y_position);
		    }
		    
		    if (indices) delete [] indices;
		    indices = NULL;

		    sequence_of_scores[insert_i]=gap_score; 

		    score+=gap_score;	      
		}
	  
	      // if there was something inserted and there is something to the right of the inserted local solution
	      // close gap to the right

		if ((insert_i>-1) 
		    &&
		    (sequence_of_xsteps[insert_i+local_steps+1-1]!=-1) 
		    &&
		    ((sequence_of_xsteps[insert_i+local_steps+1-1]
		      +model[sequence_of_states[insert_i+local_steps+1]].get_letters_to_read_x())
		     == sequence_of_xsteps[insert_i+local_steps+1]) 
		    &&
		    (sequence_of_ysteps[insert_i+local_steps+1-1]!=-1) 
		    &&
		    ((sequence_of_ysteps[insert_i+local_steps+1-1]
		      +model[sequence_of_states[insert_i+local_steps+1]].get_letters_to_read_y())
		     == sequence_of_ysteps[insert_i+local_steps+1]) )
		{
	    
		  // read letters of sequence into index

		    gap_score=Logzero;
		    
		    int* indices     = NULL;
		    indices = new int[model[sequence_of_states[insert_i+local_steps+1]].get_number_of_dimensions_of_emission_scores()+1];
		    int deltax=model[sequence_of_states[insert_i+local_steps+1]].get_letters_to_read_x();
		    int deltay=model[sequence_of_states[insert_i+local_steps+1]].get_letters_to_read_y();
		    
		    int xsteps=sequence_of_xsteps[insert_i+local_steps+1];
		    int ysteps=sequence_of_ysteps[insert_i+local_steps+1];
		    for (int i=0; i<deltax; i++)
		    {
			indices[i+1]=x->letter(xsteps - (deltax-i));

		    }
		    for (int j=deltax; j<deltay+deltax; j++)
		    {
			indices[j+1]=y->letter(ysteps - (deltay+deltax-j));

		    }
		  
		    indices[0]=sequence_of_states[insert_i+local_steps+1-1];
		    
		    int new_state = sequence_of_states[insert_i + local_steps + 1];
		    int old_state = sequence_of_states[insert_i + local_steps + 1 - 1];
		    
		    int x_position = xsteps - deltax;
		    int y_position = ysteps - deltay;

		    if ((deltax+deltay) > 0)
		    {
			int linear_index_emission = indices[1];
			for (int p=2; p<model[sequence_of_states[insert_i+local_steps+1]].get_number_of_dimensions_of_emission_scores()+1; p++)
			{
			    linear_index_emission = (linear_index_emission << shift_index) | indices[p];
			}
			
			gap_score=this->new_get_transition_score(NULL, 
								 old_state, new_state, 
								 x, x_position, y, y_position) +
			    this->model[new_state].get_emission_score(x, x_position, y, y_position,
								      linear_index_emission);
		    }
		    else
		    {
			gap_score=this->new_get_transition_score(NULL, 
								 old_state, new_state, 
								 x, x_position, y, y_position);
		    }
		    
		    if (indices) delete [] indices;
		    indices = NULL;

		    sequence_of_scores[insert_i+local_steps+1]=gap_score; 

		    score+=gap_score;
		}
	    }
	}
	else
	{
	    if ((steps-1) < max_index)
	    {
		cout << "ERROR Hmm::add_local_solution: no space left to implement local solution."
		     << "steps (" << steps << ") - 1 = " << steps-1 << " > max_index (" << max_index << ")).\n" << flush;
		check+=1;
	    }
	}

    }
    return(check);
}

int Hmm::add_local_state_path(const Sequence *x, 
				  const Sequence *y, 
				  const StatePath *local_state_path) { 

    // NOTE: - local_state_path is added to the existing state_path
  
    int check=0;
    int max_index = x->length() + y->length() + 2;
    
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::add_local_state_path : mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }
    
    // check sequences
    
    if (x->length() < 1) 
    {
	cout << "ERROR: Hmm::add_local_state_path : length of sequence x (" << x->length() << ") < 1.\n" << flush;
	check+=1;
    }
    // check local_state_path:
    
    // check that arrays are not NULL
    
    if (local_state_path->l_states == NULL) {
	
	cout << "ERROR Hmm::add_local_state_path: l_states is NULL.\n";
	check+=1;
    }
    if (local_state_path->l_xsteps == NULL) {
	
	cout << "ERROR Hmm::add_local_state_path: l_xsteps is NULL.\n";
	check+=1;
    }
    if (local_state_path->l_ysteps == NULL) {
	
	cout << "ERROR Hmm::add_local_state_path: l_ysteps is NULL.\n";
	check+=1;
    }
    if (local_state_path->l_scores == NULL) {
	
	cout << "ERROR Hmm::add_local_state_path: l_scores is NULL.\n";
	check+=1;
    }
    
    // check value of l_steps

    if (local_state_path->l_steps < 0) {
	
	cout << "ERROR Hmm::add_local_state_path: local_steps (" << local_state_path->l_steps << ") must be >=0.\n" << flush;
	check+=1;
    }
    
    // check value of l_score

    if (local_state_path->l_score == Logzero) {
	
	cout << "ERROR Hmm::add_local_state_path: local_score (" << local_state_path->l_score << ") is Logzero.\n" << flush;
	check+=1;
    }
    
    if (check == 0) {
	
	int i, j, p;
	int shift_index = bitshift(alphabet);
	int insert_i    = 0;
	
	// check where new solution has to be added
	
	if (! ((sequence_of_xsteps[0] == -1) && (sequence_of_ysteps[0] == -1))) { // if solution is not empty
	    for (i=0; i<max_index-1; i++) {
		
		// if existing solution ends at position [i]
		
		if ((sequence_of_xsteps[i]!=-1) && (sequence_of_xsteps[i+1]==-1) &&
		    (sequence_of_ysteps[i]!=-1) && (sequence_of_ysteps[i+1]==-1)) {
		    
		    insert_i = i+1;

		    break;
		}
	    }
	}
	
	// determine final number of steps
	
	steps = insert_i + local_state_path->l_steps;

	// add new solution at insert_i, if there is enough space left in arrays
	
	if ((steps-1) < max_index) {
	    
	    // put local solution in entries with indices from insert_i to insert_i+local_state_path->l_steps+1 inclusive

	    for (i=0; i<local_state_path->l_steps+1; i++) {
		
		//cout<<"local_state_path->l_xsteps["<<i<<"] : "<<local_state_path->l_xsteps[i]<<endl;
		
		score += local_state_path->l_scores[i];
		
		sequence_of_scores[insert_i+i] = local_state_path->l_scores[i];
		sequence_of_states[insert_i+i] = local_state_path->l_states[i];
		sequence_of_xsteps[insert_i+i] = local_state_path->l_xsteps[i]; 
		sequence_of_ysteps[insert_i+i] = local_state_path->l_ysteps[i]; 
	    }

	    // if necessary, close gaps to the left by calculating the transition score
	    
	    Score gap_score = Logzero;
	    
	    if ((insert_i>0) &&
		(sequence_of_xsteps[insert_i-1]!=-1) &&
		((sequence_of_xsteps[insert_i-1]+model[sequence_of_states[insert_i]].get_letters_to_read_x()) == sequence_of_xsteps[insert_i]) &&
		(sequence_of_ysteps[insert_i-1]!=-1) &&
		((sequence_of_ysteps[insert_i-1]+model[sequence_of_states[insert_i]].get_letters_to_read_y()) == sequence_of_ysteps[insert_i])) {
		
		int new_state = sequence_of_states[insert_i];
		int old_state = sequence_of_states[insert_i-1];
		int x_position = sequence_of_xsteps[insert_i] - model[sequence_of_states[insert_i]].get_letters_to_read_x();
		int y_position = sequence_of_ysteps[insert_i] - model[sequence_of_states[insert_i]].get_letters_to_read_y();

		gap_score                     = this->new_get_transition_score(NULL, old_state, new_state, x, x_position, y, y_position);
		sequence_of_scores[insert_i] += gap_score; 
		score                        += gap_score;	      

	    }
	}
	else { // (steps-1) < max_index
	    
	    cout << "ERROR Hmm::add_local_state_path: no space left to implement local solution."
		 << "steps (" << steps << ") - 1 = " << steps-1 << " > max_index (" << max_index << ")).\n" << flush;
	    check+=1;
	}

    }
    
    return(check);
}

int Hmm::check_consistency_of_solution(const Sequence *x, const Sequence *y)
{
    int check=0;
    
    if (mirrored != 0)
    {
	cout << "ERROR Hmm::check_consistency_of_solution: mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }

    if (sequence_of_states==NULL)
    {
      cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_states is NULL.\n" << flush;
      check+=1;
    }
    
    if (sequence_of_scores==NULL)
    {
	cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_scores is NULL.\n" << flush;
	check+=1;
    }

    if (sequence_of_xsteps==NULL)
    {
	cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_xsteps is NULL.\n" << flush;
	check+=1;
    }
    
    if (sequence_of_ysteps==NULL)
    {
	cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_ysteps is NULL.\n" << flush;
	check+=1;
    }

    if ((score==0) || (score==Logzero))
    {
	cout << "ERROR Hmm::check_consistency_of_solution: score (" << score << ") == 0 or Logzero.\n" << flush;
	check+=1;
    }

    if (steps<1)
    {
	cout << "ERROR Hmm::check_consistency_of_solution: steps (" << steps << ") < 1.\n" << flush;
	check+=1;
    }

    if (check==0)
    {
	if (sequence_of_states[0]!=0)
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_states[0] (" 
		 << sequence_of_states[0] << ") != 0\n" << flush;
	    check+=1;	  
	}
	
	if (sequence_of_scores[0]!=0)
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_scores[0] (" 
		 << sequence_of_scores[0] << ") != 0\n" << flush;
	    check+=1;	  
	}
	
	if (sequence_of_xsteps[0]!=0)
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_xsteps[0] (" 
		 << sequence_of_xsteps[0] << ") != 0\n" << flush;
	    check+=1;	  
	}
	
	if (sequence_of_ysteps[0]!=0)
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_ysteps[0] (" 
		 << sequence_of_ysteps[0] << ") != 0\n" << flush;
	    check+=1;	  
	}
      
	if (sequence_of_states[steps]!=number_of_states-1)
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_states[" << steps << "] (" 
		 << sequence_of_states[steps] << ") != number_of_states-1 (" << number_of_states-1 << ")\n" << flush;
	    check+=1;	  
	}
	
	if (sequence_of_xsteps[steps]!=x->length())
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_xsteps[" << steps<< "] (" 
		 << sequence_of_xsteps[steps] << ") != x->length() (" << x->length() << ")\n" << flush;
	    check+=1;	  
	}
	
	if (sequence_of_ysteps[steps]!=y->length())
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: sequence_of_ysteps[" << steps<< "] (" 
		 << sequence_of_ysteps[steps] << ") != y->length() (" << y->length() << ")\n" << flush;
	    check+=1;	  
	}

	int shift_index = bitshift(alphabet);
	
	Score score_sum=0;
	Score test_score=0;
	
	int old_state=0;
	int new_state=0;
	int xsteps=0;
	int ysteps=0;
	int deltax=0;
	int deltay=0;
	
	int i;

	for (i=0; i<steps; i++) 
	{
	    score_sum+=sequence_of_scores[i+1];
	    old_state=sequence_of_states[i];
	    new_state=sequence_of_states[i+1];
	    deltax=model[sequence_of_states[i+1]].get_letters_to_read_x();
	    deltay=model[sequence_of_states[i+1]].get_letters_to_read_y();

	    // check that there is a transition between old_state and new_state
	    
	    if (model[new_state].get_transition_score(old_state)>Logzero) // if old_state -> new_state
	    {                                                            // transition is allowed	       

		// check difference in xsteps and ysteps 
		
		if (sequence_of_xsteps[i]+deltax != sequence_of_xsteps[i+1])
		{
		    cout << "ERROR Hmm::check_consistency_of_solution:"
			 << " sequence_of_xsteps[" << i << "] (" << sequence_of_xsteps[i] 
			 << ") + deltax (" << deltax 
			 << ") != sequence_of_xsteps[" << i+1 << "] (" << sequence_of_xsteps[i+1] << ")\n" << flush; 
		    check+=1;	  
		}
		if (sequence_of_ysteps[i]+deltay != sequence_of_ysteps[i+1])
		{
		    cout << "ERROR Hmm::check_consistency_of_solution:"
			 << " sequence_of_ysteps[" << i << "] (" << sequence_of_ysteps[i] 
			 << ") + deltay (" << deltay 
			 << ") != sequence_of_xsteps[" << i+1 << "] (" << sequence_of_xsteps[i+1] << ")\n" << flush; 
		    check+=1;	  
		}
		
		// check score
		
		test_score=Logzero;
		xsteps=sequence_of_xsteps[i+1];
		ysteps=sequence_of_ysteps[i+1];
		
		int* indices = NULL;
		indices = new int[model[new_state].get_number_of_dimensions_of_emission_scores()+1];

		for (int k=0; k<deltax; k++)
		{
		    indices[k+1] = x->letter(xsteps - (deltax-k));

		}
		for (int j=deltax; j<deltay+deltax; j++)
		{
		    indices[j+1] = y->letter(ysteps - (deltay+deltax-j));

		}
		
		indices[0] = old_state;
		
		int x_position = xsteps - deltax;
		int y_position = ysteps - deltay;

		if ((deltax+deltay) > 0)
		{
		    int linear_index_emission = indices[1];
		    for (int p=2; p<model[new_state].get_number_of_dimensions_of_emission_scores()+1; p++)
		    {
			linear_index_emission = (linear_index_emission << shift_index) | indices[p];
		    }

		    test_score=this->new_get_transition_score(NULL, 
							      old_state, new_state, 
							      x, x_position, y, y_position) +
			this->model[new_state].get_emission_score(x, x_position, y, y_position,
								  linear_index_emission);
		}
		else
		{
		    test_score=this->new_get_transition_score(NULL, 
							      old_state, new_state, 
							      x, x_position, y, y_position);
		}
		
		if (indices) delete [] indices;
		indices = NULL;
	    
		if ( abs(sequence_of_scores[i+1]-test_score) > Max_deviation)
		{
		    cout << "ERROR Hmm::check_consistency_of_solution:"
			 << " sequence_of_scores[" << i+1 << "] (" << sequence_of_scores[i+1] 
			 << ") != test_score (" << test_score << ")\n" << flush;
		    check+=1;	  
		}	      
	    }
	    else
	    {
		cout << "ERROR Hmm::check_consistency_of_solution: a transition between sequence_of_states[" 
		     << i << "] (" << sequence_of_states[i] << ") -> sequence_of_states[" << i+1 
		     << "] (" << sequence_of_states[i+1] << ") is not allowed.\n" << flush;
		check+=1;	  
	    }
	}

	if (abs(score_sum-score)>Max_deviation)
	{
	    cout << "ERROR Hmm::check_consistency_of_solution: | score_sum (" << score_sum 
		 << ") - score (" << score << ") | = " << abs(score_sum-score) << " too large.\n" << flush;
	    check+=1;	  
	}
    }
    return(check);
}


int Hmm::viterbi_rectangle(Score*** viterbi_rectangle_array,
			       const Hmm *s, 
			       const Sequence *x, const int x_start, const int x_end,
			       const Sequence *y, const int y_start, const int y_end, 
			       const int* start_states, const int end_state,
			       const int number_of_start_states,
			       const int x_margin, const int y_margin,
			       const int offset,
			       // output values
			       int* x_start_traceback, 
			       int* y_start_traceback, 
			       int* start_state_traceback)
{
    // function returns 0 if o.k., else not o.k.
    //
    // note: this function can either be used 
    //
    //       1.) with start_states and an end_state (in which case
    //       viterbi_rectangle_array == NULL and x_margin and y_margin == 0). Memory
    //       for viterbi_rectangle_array is allocated and deleted within function.
    //       
    //       2.) with a pre-initialised viterbi_rectangle_array (in which case start_states
    //        == NULL and end_state == 0 and number_of start_states and number_of_end_states == 0).
    //       Memory for viterbi_rectangle_array is allocated and deleted outside function.
    
    int check=0;
    int n_of_states=0;
    n_of_states=s->get_number_of_states();
    
    // check whether pre-initialised viterbi_rectangle or start_ and end_states shall be used
    
    if ((viterbi_rectangle_array == NULL)    &&
	(number_of_start_states == 0))
    {
	cout << "ERROR Hmm::viterbi_rectangle: viterbi_rectangle_array (" << viterbi_rectangle_array
	     << ") = NULL and number_of_start_states (" << number_of_start_states 
	     << ") == 0. Use this function either with a pre-initialised viterbi_rectangle_array != NULL "
	     << "or with nonempty sets of start_ and end_states.\n" << flush;
	check+=1;
    }
    else if ((viterbi_rectangle_array != NULL)    &&
	     (number_of_start_states != 0))
    {
	cout << "ERROR Hmm::viterbi_rectangle: viterbi_rectangle_array (" << viterbi_rectangle_array
	     << ") != NULL and number_of_start_states (" << number_of_start_states 
	     << ") != 0. Use this function either with a pre-initialised viterbi_rectangle_array != NULL "
	     << "or with nonempty sets of start_ and end_states.\n" << flush;
	check+=1;
    }
    else // if input values are not completely inconsistent
    {
	if ((viterbi_rectangle_array      == NULL)   &&
	    (number_of_start_states != 0))
	{
	    // if function shall be used with start_ and end_states
	    
	    // check start states
	    
	    if (number_of_start_states<1)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: number_of_start_sites (" 
		     << number_of_start_states << ") <1.\n" << flush;
		check+=1;
	    }
	    else
	    {
		if ((x_start==0) && (y_start==0))
		{
		    if ((number_of_start_states!=1) || (start_states[0]!=0))
		    {
			cout << "ERROR Hmm::viterbi_rectangle: or "
			     << "x_start=y_start=0 only "
			     << "one start_state is possible, the Startstate (i.e. the number_of_start_states ("
			     << number_of_start_states << ") has to be 1 and start_states[0] ("
			     << start_states[0] << ") has to be 0).\n" << flush;
			check+=1;
		    }
		}
		else // if we don't start at (0,0)
		{
		    // check that each start state lies in the allowed range (1...n_of_states-2)
		    // i.e. Start and End states are not allowed 
		    
		    for (int i=0; i<number_of_start_states; i++)
		    {
			if ((start_states[i]<1) || (start_states[i]>n_of_states-2))
			{
			    cout << "ERROR Hmm::viterbi_rectangle: start_states[" << i 
				 << "] (" << start_states[i] << ") out of range [1, " << n_of_states-2
				 << "].\n" << flush;
			    check+=1;
			}
		    }
		}
	    }
	    
	  // check end state

	    if ((x_end==x->length()) && (y_end==y->length()))
	    {
		if (end_state != n_of_states-1)
		{
		    cout << "ERROR Hmm::viterbi_rectangle: end_state (" 
			 << end_state << ") for x_end (" << x_end << ") = x->length() (" << x->length() 
			 << ") and y_end (" << y_end << ") = y->length() (" << y->length() 
			 << ") has to be End state (" << n_of_states-1 << ").\n" << flush;
		    check+=1;
		}
	    }
	    else // if we don't stop at ends of both sequences
	    {
		// check that end state lies in the allowed range (1...n_of_states-2)
		// i.e.Start state and End states are not allowed
		
		if ((end_state<1) || (end_state>n_of_states-2))
		{
		    cout << "ERROR Hmm::viterbi_rectangle: end_state "
			 << end_state << ") out of range [1, " << n_of_states-2 << "].\n" << flush;
		    check+=1;
		}
	    }
	}

	if ((viterbi_rectangle_array      != NULL)   &&
	    (number_of_start_states == 0))
	{
	    // if function shall be used with pre-initialised viterbi_rectangle_array

	  // check values of x_margin and y_margin

	    if ((x_start==0) && (x_margin!=0))
	    {
		cout << "ERROR: Hmm::viterbi_rectangle : x_margin (" << x_margin 
		     << ") has to be zero for x_start (" 
		     << x_start << ") = 0.\n" << flush;
		check+=1;
	    }
	  
	    if ((y_start==0) && (y_margin!=0))
	    {
		cout << "ERROR: Hmm::viterbi_rectangle : y_margin (" << y_margin 
		     << ") has to be zero for y_start (" 
		     << y_start << ") = 0.\n" << flush;
		check+=1;
	    }
	  
	    if (x_margin < 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: x_margin (" << x_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	  
	    if (y_margin < 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: y_margin (" << y_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	  
	    // check that x_margin and y_margin are compatible with values of x_start, x_end and y_start and y_end
	    
	    if ((x_end-x_start+1) < x_margin)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: x_margin (" << x_margin
		     << ") > (x_end (" << x_end << ") - x_start (" << x_start << ") + 1) = " 
		     << x_end-x_start+1 << "\n" << flush;
		check+=1;
	    }
	  
	    if ((y_end-y_start+1) < y_margin)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: y_margin (" << y_margin
		     << ") > (y_end (" << y_end << ") - y_start (" << y_start << ") + 1) = " 
		     << y_end-y_start+1 << "\n" << flush;
		check+=1;
	    }
	}
    }

    // check Hmm s
    
    if (s->get_mirrored() != 0)
    {
	cout << "ERROR: Hmm::viterbi_rectangle : Hmm : mirrored (" 
	     << s->get_mirrored() << ") must be 0.\n" << flush;
	check+=1;
    }
    if (s->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::viterbi_rectangle : Hmm : number of states = " 
	     << s->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (s->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::viterbi_rectangle : Hmm : state 0 != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (s->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::viterbi_rectangle : Hmm : last state of type != End \n" << flush;
	cout << "model["<<number_of_states-1<<"].get_letters_to_read() : "<<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (int i=0; i<number_of_states; i++)
    {
	if (s->model[i].get_alphabet()!=alphabet) 
	{
	    cout << "ERROR: Hmm::viterbi_rectangle : Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << s->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_mirrored()!=s->get_mirrored()) 
	{
	    cout << "ERROR: Hmm::viterbi_rectangle : Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << s->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << s->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (s->model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::viterbi_rectangle : Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << s->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}
 
	if ((i!=0) && (i!=(number_of_states-1))  && (model[i].get_letters_to_read() ==0 ))
	{
	    cout << "ERROR: Hmm::viterbi_rectangle : Hmm : state i = " 
		 << i << " should emit ! \n"<< flush;
	    check+=1;
	}

    }
  
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::viterbi_rectangle : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::viterbi_rectangle : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }
    
    // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::viterbi_rectangle: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::viterbi_rectangle: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::viterbi_rectangle: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::viterbi_rectangle: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    if ((offset != 0) && (offset != -1))
    {
	cout << "ERROR Hmm::viterbi_rectangle: offset (" << offset << ") must be either 0 or -1.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {
	int use_internal_viterbi_rectangle = 0; // 1 = yes (use interal), 0 = no (use external)
	
	// initialise values of output variables

	(*x_start_traceback)     = 0;
	(*y_start_traceback)     = 0; 
	(*start_state_traceback) = 0;
	
	// variables for local solution
	
	int local_steps=0;
	Score local_score=0;
	Score* local_sequence_of_scores=NULL;
	int* local_sequence_of_states=NULL;
	int* local_sequence_of_xsteps=NULL;
	int* local_sequence_of_ysteps=NULL;
	
	if (viterbi_rectangle_array == NULL) // allocate memory for viterbi_rectangle_array if pre_initialised
	    // viterbi_rectangle_array is not used
	{
	    use_internal_viterbi_rectangle = 1; // use internal viterbi_rectangle_array
	    check+= allocate_memory_for_viterbi_rectangle(s,
							  x_start, x_end,
							  y_start, y_end,
							  1, // direction == 1
							  &viterbi_rectangle_array);
	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when allocating memory "
		     << "for viterbi_rectangle_array.\n" << flush;
	    }
	}
	
	if (check==0)
	{

	    check += calculate_viterbi_rectangle(viterbi_rectangle_array, 
						 s,
						 NULL, // unmirrored Hmm not needed
						 x, x_start, x_end,
						 y, y_start, y_end, 
						 start_states,  
						 number_of_start_states,
						 x_margin, y_margin, 
						 1); // direction == 1
	    
	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
		     << "calculate_viterbi_rectangle.\n" << flush;
	    }
	}
	
	if (check==0)
	{

	    check+=delete_memory_for_state_path(&local_sequence_of_scores,
						&local_sequence_of_states,
						&local_sequence_of_xsteps,
						&local_sequence_of_ysteps);
	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
		     << "delete_memory_for_state_path.\n" << flush;
	    }
	}
	
	if (check==0)
	{

	    check+=allocate_memory_for_state_path(x_start, x_end,
						  y_start, y_end,
						  &local_sequence_of_scores,
						  &local_sequence_of_states,
						  &local_sequence_of_xsteps,
						  &local_sequence_of_ysteps);
	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
		     << "allocate_memory_for_state_path.\n" << flush;
	    }
	}
	
	if (check==0)
	{

	    check+=retrieve_state_path_from_viterbi_rectangle(viterbi_rectangle_array,
							      s, 
							      NULL, // unmirrored Hmm not needed
							      x, x_start, x_end,
							      y, y_start, y_end, 
							      end_state, 
							      1, // direction == 1
							      local_steps,
							      local_score,
							      local_sequence_of_scores,
							      local_sequence_of_states,
							      local_sequence_of_xsteps,
							      local_sequence_of_ysteps);

	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
		     << "retrieve_state_path_from_viterbi_rectangle.\n" << flush;
	    }
	    else
	    {
		// fill output values
		
		(*x_start_traceback)     = local_sequence_of_xsteps[0];
		(*y_start_traceback)     = local_sequence_of_ysteps[0];
		(*start_state_traceback) = local_sequence_of_states[0];
	    }

	}
	
	if (check==0)
	{
	    check+=add_local_solution(x,
				      y,
				      local_steps+offset,
				      local_score,
				      local_sequence_of_scores,
				      local_sequence_of_states,			 
				      local_sequence_of_xsteps,
				      local_sequence_of_ysteps,
				      1); // direction == 1
	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
		     << "add_local_solution.\n" << flush;
	    }
	}
	
	if (check==0)
	{

	    check+=delete_memory_for_state_path(&local_sequence_of_scores,
						&local_sequence_of_states,
						&local_sequence_of_xsteps,
						&local_sequence_of_ysteps);
	    if (check != 0)
	    {
		cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
		     << "delete_memory_for_state_path.\n" << flush;
	    }
	}
	
	if (check==0)
	{
	    if (use_internal_viterbi_rectangle == 1) // if viterbi_rectangle_array was allocated internally
	    {

		check+=delete_memory_for_viterbi_rectangle(s, 
							   x_start, x_end,
							   1, // direction == 1
							   &viterbi_rectangle_array);
		if (check != 0)
		{
		    cout << "ERROR Hmm::viterbi_rectangle: error occured when using "
			 << "delete_memory_for_viterbi_rectangle.\n" << flush;
		}
	    }
	}
    }
    return(check);
}

int Hmm::get_sequence_labels_for_state_path(int*** const x_labels, 
						int*** const y_labels,
						const Sequence &x, 
						const Sequence &y,
						model_parameters* const MP) 
    
// requires non-condensed state path solution to exist
{
    int check = 0;
    
    if ((score > Logzero)            && 
	(steps > -1)                 &&
	(sequence_of_scores != NULL) &&     
	(sequence_of_states != NULL) &&
	(sequence_of_xsteps != NULL) &&
	(sequence_of_ysteps != NULL)) {

	int i, j, k, x_count, y_count, delta_x, delta_y, state, number_of_state_labels;
	
	number_of_state_labels = MP->get_Total_Number_of_Annotation_Labels();

	if ((*x_labels) != NULL) delete [] (*x_labels); 
	(*x_labels) = NULL;
	(*x_labels) = new int*[number_of_state_labels];
	for( i = 0; i< number_of_state_labels; i++){
	    (*x_labels)[i] = new int[x.length()];
	}
	
	if ((*y_labels) != NULL) delete [] (*y_labels); 
	(*y_labels) = NULL;
	(*y_labels) = new int*[number_of_state_labels];
	for(i=0;i<number_of_state_labels;i++){
	    (*y_labels)[i]=new int[y.length()];
	}
	
	x_count = 0;
	y_count = 0;
	
	for (i = 0; i <= steps; i++) {
	    
	    state   = sequence_of_states[i];
	    delta_x = this->model[state].get_letters_to_read_x();
	    delta_y = this->model[state].get_letters_to_read_y();
	    
	    
	    for( k = 0; k < number_of_state_labels; k++){

		for (j = 0; j < delta_x; j++) {
		    
		    (*x_labels)[k][x_count] = this->model[state].get_state_labels(k,j);
		    x_count++;
		}
		for (j = 0; j < delta_y; j++) {
		    
		    (*y_labels)[k][y_count] = this->model[state].get_state_labels(k,j);
		    y_count++;
		}
	    }
	}
	
    } 
    else { // if state path does not exists
	
	cout << "ERROR: Hmm::get_sequence_labels_for_state_path: state path does not exist.\n" << flush;
	
	int number_of_state_labels = MP->get_Total_Number_of_Annotation_Labels();
	
	if(*x_labels){
	    for(int i=0;i<number_of_state_labels;i++){    
		if ((*x_labels)[i]) delete [] (*x_labels)[i]; 
		(*x_labels)[i] = NULL;
	    }
	    if(*x_labels) delete[] (*x_labels);
	    (*x_labels) = NULL;
	}
	
	if(*y_labels){
	    for(int i=0;i<number_of_state_labels;i++){
		if ((*y_labels)[i]) delete [] (*y_labels)[i]; 
		(*y_labels)[i] = NULL;
	    }
	    if(*y_labels) delete[] (*y_labels);
	    (*y_labels) = NULL;
	}
	    
	check++;
    }
    
    return(check);
}

void Hmm::readsuperfastexptable()
{
    superfastexptable = new double[10001];
    FILE* expfile = fopen("fastexptable.txt","r");
    const int max_line_length = 100;
    char line[max_line_length];
    long i = 0;
    while(!feof(expfile))
    {
	fgets(line,max_line_length-1,expfile);
	superfastexptable[i] = atof(line);
	i++;
    }
    fclose(expfile);
    return;
}

double Hmm::superfastinterpexp(const double x)
{
    if(x<=0)
    {
	if(x<-10)
	{
	    return superfastexptable[10000];
	}else{
	    return superfastexptable[static_cast<int>(x*1000)];
	}
    }else{
	if(x>10)
	{
	    return 1/superfastexptable[10000];
	}else{
	    return 1/superfastexptable[static_cast<int>(x*1000)];
	}
    }
}

// public functions except constructors

int Hmm::get_alphabet(void) const
{
    return(alphabet);
}

int Hmm::get_mirrored(void) const
{
    return(mirrored);
}

int Hmm::get_number_of_states(void) const
{
    return(number_of_states);
}

bool Hmm::get_pair(void) const
{
    return(pair);
}

int Hmm::set_alphabet(const int a)
{
    int check = 0;
    if(a<0)
    {
	cout<<"Error: Hmm:: set_alphabet: "
	    <<"input a("<<a<<") out of range."<<endl;
	check++;
    }
    if(check)
    {
	alphabet=0;
    }else{
	alphabet=a;
    }
    return check;
}
    
int Hmm::set_mirrored(const int m)
{
    int check = 0;
    if((m!=0)&&(m!=1))
    {
	cout<<"Error: Hmm:: set_mirrored: "
	    <<"input m("<<m<<") out of range [0,1]."<<endl;
	check++;
    }
    if(check)
    {
	mirrored = 0;
    }else{
	mirrored=m;
    }
    return check; 
}
int Hmm::set_number_of_states(const int n_s)
{
    int check = 0;
    if(n_s<0)
    {
	cout<<"Error: Hmm:: set_number_of_states: "
	    <<"input n_s("<<n_s<<") out of range."<<endl;
	check++;
    }
    if(check)
    {
	number_of_states=0;
    }else{
	number_of_states=n_s;
    }
    return check;
}

int Hmm::set_pair(const bool p)
{
    pair = p;
    return 0;
}

Score Hmm::get_transition_score(
				    const Hmm *unmirrored_s,
				    const int             n_of_from_state,
				    const int             n_of_to_state,
				    const Sequence* const x,
				    const int             x_position,
				    const Sequence* const y,
				    const int             y_position) const
{
    int   check        = 0;
    Score return_score = Logzero;
    
    if ((n_of_from_state < 0) || (n_of_from_state > (number_of_states-1)))
    {
	cout << "ERROR class Hmm::get_transition_score : number of from state ("
	     << n_of_from_state << ") out of range [0, " 
	     << number_of_states-1 << "].\n" << flush;
	check++;
    }
    if ((n_of_to_state < 0) || (n_of_to_state > (number_of_states-1)))
    {
	cout << "ERROR class Hmm::get_transition_score : number of to state ("
	     << n_of_to_state << ") out of range [0, " 
	     << number_of_states-1 << "].\n" << flush;
	check++;
    }
    if (x == NULL)
    {
	cout << "ERROR class Hmm::get_transition_score : Sequence x is NULL.\n" << flush;
	check++;
    }
    if (y == NULL)
    {
	cout << "ERROR class Hmm::get_transition_score : Sequence y is NULL.\n" << flush;
	check++;
    }
  
    if (check == 0)
    {
#ifdef _DEL
	cout << "class Hmm::get_transition_score:\n" << flush;
	cout << "from " << n_of_from_state << " -> to " << n_of_to_state 
	     << " x_pos = " << x_position << " y_pos = " << y_position << "\n" << flush;
#endif

	int mirror_value = model[n_of_to_state].get_mirrored();

	
	return_score = model[n_of_to_state].get_transition_score(n_of_from_state);

	double prob_corresponding_to_return_score = 
	    calculate_number_from_score(return_score);
	
	if ((prob_corresponding_to_return_score < 0) || (prob_corresponding_to_return_score > 1))
	{
	    cout << "ERROR: class Hmm::get_transition_score: return_score ("
		 << return_score << ") for transition " << n_of_from_state << " -> "
		 << n_of_to_state << " at position (x,y) = (" << x_position << ", " << y_position 
		 << ") corresponds to prob ("
		 << prob_corresponding_to_return_score << ") which is not within [0,1].\n" << flush;
	    return_score = Logzero;
	}
    } // if check == 0

    return(return_score);
}

Score Hmm::new_get_transition_score(
					const Hmm *unmirrored_s,
					const int             n_of_from_state,
					const int             n_of_to_state,
					const Sequence* const x,
					const int             x_position,
					const Sequence* const y,
					const int             y_position) const
{

    int   check        = 0;
    Score return_score = Logzero;
    
    if ((n_of_from_state < 0) || (n_of_from_state > (number_of_states-1)))
    {
	cout << "ERROR class Hmm::new_get_transition_score : number of from state ("
	     << n_of_from_state << ") out of range [0, " 
	     << number_of_states-1 << "].\n" << flush;
	check++;
    }
    if ((n_of_to_state < 0) || (n_of_to_state > (number_of_states-1)))
    {
	cout << "ERROR class Hmm::new_get_transition_score : number of to state ("
	     << n_of_to_state << ") out of range [0, " 
	     << number_of_states-1 << "].\n" << flush;
	check++;
    }
    if (x == NULL)
    {
	cout << "ERROR class Hmm::new_get_transition_score : Sequence x is NULL.\n" << flush;
	check++;
    }
    if (y == NULL)
    {
	cout << "ERROR class Hmm::new_get_transition_score : Sequence y is NULL.\n" << flush;
	check++;
    }
	 
    if (check == 0)
    {
	
#ifdef _DEL
	cout << "class Hmm::new_get_transition_score:\n" << flush;
	cout << "from " << n_of_from_state << " -> to " << n_of_to_state 
	     << " x_pos = " << x_position << " y_pos = " << y_position << "\n" << flush;
#endif

	int mirror_value = model[n_of_to_state].get_mirrored();

	
	return_score = model[n_of_to_state].get_transition_score(n_of_from_state);

	double prob_corresponding_to_return_score = 
	    calculate_number_from_score(return_score);
	
	if ((prob_corresponding_to_return_score < 0) || (prob_corresponding_to_return_score > 1))
	{
	    cout << "ERROR: class Hmm::new_get_transition_score: return_score ("
		 << return_score << ") for transition " << n_of_from_state << " -> "
		 << n_of_to_state << " at position (x,y) = (" << x_position << ", " << y_position 
		 << ") corresponds to prob ("
		 << prob_corresponding_to_return_score << ") which is not within [0,1].\n" << flush;
	    return_score = Logzero;
	}
    } // if check == 0
    
    return(return_score);
}

int Hmm::implement_state(int i, Hmm_State* state)
{
    int check=0;

    if (state->get_mirrored() != mirrored)
    {
	cout << "ERROR class Hmm::constructor : mirrored flag of state a = " << state->get_mirrored()
	     << " does not match mirrored flag of Hmm a = " << mirrored << "\n" << flush;
	check+=1;
    }
    if (state->get_alphabet() != alphabet)
    {
	cout << "ERROR class Hmm::constructor : alphabet of state a = " << state->get_alphabet()
	     << " does not match alphabet of Hmm a = " << alphabet << "\n" << flush;
	check+=1;
    }
    if (state->get_number_of_states() != number_of_states)
    {
	cout << "ERROR class Hmm::constructor : number of states of state n = " << state->get_number_of_states()
	     << " does not match number of states of Hmm n = " << number_of_states << "\n" << flush;
	check+=1;
    }
    if ( (i<0) || (i>(number_of_states-1)) )
    {
	cout << "ERROR class Hmm::implement_state : state number out of range i = " << i 
	     << " Has to be in [0, " << number_of_states-1 << "].\n" << flush;
	check+=1;
    }
    
    
    if ((i==0) && (state->get_letters_to_read()!=0))
    {
	cout << "ERROR class Hmm::constructor : state 0 has to be of type Start\n" << flush;
	cout <<"state->get_letters_to_read() : "<<state->get_letters_to_read()<<endl;
	check+=1;
    }

    if ((i==number_of_states-1) && (state->get_letters_to_read()!=0))
    {
       cout << "ERROR class Hmm::constructor : last state has to be of type End\n" << flush;
       cout<<"state->get_letters_to_read() : "<<state->get_letters_to_read()<<endl;
       check+=1;
   }
    

    if ((i!=number_of_states-1) && (i!=0) && (state->get_letters_to_read()==0))
    {
	cout << "ERROR class Hmm::constructor : intermediate state number i = " << i 
	   << " cannot be of type Start or End.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {
	model[i]=*state;
    }
    return(check);
}

int Hmm::build_connected_pairhmm() const
{
  // for each state of pairhmm set info on next states

    int check = 0;
    int i;
    int j;
    
    for (i=0; i<number_of_states; i++)
    {
	// get set of next states for state i
	// for every state i see which state j has state i has previous state
	// => these states j are then the next states of state i
	
	int number_of_next_states = 0;
	
	for (j=0; j<number_of_states; j++)
	{
	    if (model[j].is_state_previous_state(i) == 1)
	    {
		number_of_next_states++;
	    }
	}

	array<int> numbers_of_next_states(1);
	numbers_of_next_states.SetDimension(0, number_of_next_states);
	array<int> special_flags_to_next_states(1);
	special_flags_to_next_states.SetDimension(0, number_of_next_states);
      
	int count = 0;
	
	for (j=0; j<number_of_states; j++)
	{
	    if (model[j].is_state_previous_state(i) == 1)
	    {
		numbers_of_next_states.SetElement(count, j);
		if (model[j].is_transition_to_previous_state_special(i) == 1)
		{
		    special_flags_to_next_states.SetElement(count, 1);
		}
		else
		{
		    special_flags_to_next_states.SetElement(count, 0);
		}
		count++;
	    }
	}
	
	check += model[i].set_info_on_transitions_to_next_states(numbers_of_next_states,
								 special_flags_to_next_states);
	
    }
    return(check);
}

int Hmm::check_consistency() const
{
    int check=0;
 
    if (number_of_states<2) {
	check+=1;
    }
 
    for (int j=0; j<number_of_states; j++)
    {
	if (model[j].get_alphabet() != alphabet)
	{
	    check+=1;
	}
	if (model[j].get_mirrored() != mirrored)
	{
	    check+=1;
	}
	if (model[j].get_number_of_states() != number_of_states)
	{
	    check+=1;
	}
	
	if ((j==0) && (model[j].get_letters_to_read()!=0))
	{
	    check+=1;
	}
	if ((j==number_of_states-1) && (model[j].get_letters_to_read()!=0))
	{
	    check+=1;
	}
	if ((j!=number_of_states-1) && (j!=0) && (model[j].get_letters_to_read() == 0))
	{
	    check+=1;
	}
    }
    return(check);
}


int Hmm::check_consistency_of_probs() const
{

    int probs_exist=0;
    int check = 0;
 
    for (int j=0; j<number_of_states-1; j++) // check all states except End state
    {
	if (model[j].get_number_of_dimensions_of_transition_probs()==1) // check if transition probs exist
	{
	    if (model[j].get_letters_to_read()!=0)
	    {
		if (model[j].get_number_of_dimensions_of_emission_probs()>0)	      
		{
		    probs_exist++;
		}
		else
		{
		    cout << "ERROR: Hmm::check_consistency_of_probs: state " << j 
			 << " has emission prob array not filled.\n" << flush;		  
		    check++;
		}
	    }
	    else
	    {
		probs_exist++;
	    }
	}
	else
	{
	    cout << "ERROR: Hmm::check_consistency_of_probs: state " << j 
		 << " has transition prob array not filled.\n" << flush;		  
	    check++;
	}
    }
    if (probs_exist!=number_of_states-1)
    {
	cout << "ERROR: Hmm::check_consistency_of_probs: not all states have prob arrays filled or there are silent states.\n" << flush;
	check++;
    }
    
    Prob sum_of_probs=0;
 
    if (number_of_states<2) 
    {
	cout<<"Error: Hmm: check_consistency_of_probs: number_of_states : "<<number_of_states<<endl;
	check++;
    }
    
    for (int j=0; j<number_of_states; j++)
    {
	if (model[j].get_alphabet() != alphabet)
	{
	    cout<<"Error: Hmm: check_consistency_of_probs: not consistent in number of alphabet "<<endl;
	    cout<<"model["<<j<<"].get_alphabet() : "<<model[j].get_alphabet()<<" alphabet : "<<alphabet<<endl;
	    check++;
	}
	if (model[j].get_mirrored() != mirrored)
	{
	    cout<<"Error: Hmm: check_consistency_of_probs: not consistent in mirror "<<endl;
	    cout<<"model["<<j<<"].get_mirrored() : "<<model[j].get_mirrored()<<" mirrored : "<<mirrored<<endl;
	    check++;
	}
	if (model[j].get_number_of_states() != number_of_states)
	{
	    cout<<"Error: Hmm: check_consistency_of_probs: not consistent in number_of_states "<<endl;
	    cout<<"model["<<j<<"].get_number_of_states() : "<<model[j].get_number_of_states()<<" number_of_states : "<<number_of_states<<endl;
	    check++;
	}
	
	if ((j==0) && (model[j].get_letters_to_read()!=0))
	{
	    check+=1;
	}
	if ((j==number_of_states-1) && (model[j].get_letters_to_read()!=0))
	{
	    check+=1;
	}
	
	if ((j!=number_of_states-1) && (j!=0) && (model[j].get_letters_to_read()==0))
	{
	    check+=1;
	}
    }
 
    if (!check)
    {
	for (int j=0; j<number_of_states; j++)
	{
	    sum_of_probs=0;
	    
	    // check transition_probs, when no special transition
	    
	    if (j!=number_of_states-1) // for all states except End state
	    {
		for (int k=0; k<number_of_states; k++)  
		{
		    sum_of_probs+=model[j].get_transition_prob(k);
		}
		if (abs(sum_of_probs-1.)>Max_deviation) 
		{
		    check++;
		    cout << "ERROR: Hmm::check_consistency_of_probs : transition probs of state " 
			 << j << " != 1 (sum_of_probs - 1 = " << sum_of_probs-1 << ")\n" << flush;
		}
	    }	    

	    // check emission probs for all states except silent states 
	      
	    int max=0;
	    int d=0;
	    
	    if ((j!=0) && (j!=number_of_states-1))
	    {
		sum_of_probs=0;
		max=0;
		d=0;
		
		d=model[j].get_letters_to_read();
		if (d>0) 
		{
		    max= static_cast<int>(pow( static_cast<float>(alphabet), static_cast<float>(d)));
		}
		else
		{
		    max=0;
		}
 
		if (d>0) // if this is not a silent state
		{
		    for (int i=0; i<max; i++)
		    {
			sum_of_probs+=model[j].get_emission_prob(i); 		  	      
		    }	     
		    if (abs(sum_of_probs-1.)>Max_deviation) 
		    {
			check++;
			cout << "ERROR: Hmm::check_consistency_of_probs : emission probs of state " 
			     << j << " != 1 (sum_of_probs - 1 = " << sum_of_probs-1 << ")\n" << flush;
		    }
		}
	    }
	}
    }
    return(check);
}


void Hmm::print(std::ostream &o ) const
{
    o << "number_of_states = " << number_of_states << "\n" << flush;
    o << "alphabet         = " << alphabet << "\n" << flush;
    o << "mirrored         = " << mirrored << "\n" << flush;
    o << "pair HMM         = ";
    if(pair)
    {
	o<<"yes\n"<<flush;
    }else{
	o<<"no\n"<<flush;
    }
    for (int j=0; j<number_of_states; j++)
    {
	o << "state            = " << j << "\n" << flush;
	o << "=============================\n" << flush;
	model[j].print(o);
    }
    return;
}

int Hmm::get_steps(void) const
{
    return(steps);
}
 

Score Hmm::get_score(void) const
{
    return(score);
}
 

int Hmm::get_state_in_alignment(int i) const
{
#ifdef _DEBUG
    if ((i<0) || (i>steps))
    {
	cout << "ERROR class Hmm::get_state_in_alignment : index i out of range [0, "
	     << steps << "].\n" << flush;
    }
#endif
    return(sequence_of_states[i]);
}

Score Hmm::get_local_score_in_alignment(int i) const
{
#ifdef _DEBUG
    if ((i<0) || (i>steps))
    {
	cout << "ERROR class Hmm::get_local_score_in_alignment : index i out of range [0, "
	     << steps << "].\n" << flush;
}
#endif
    return(sequence_of_scores[i]);
}

void Hmm::print_results(std::ostream &o) const
{
    Score score_sum=0;
    int steps_x=0;
    int steps_y=0;

    o << "steps          = " << steps << "\n" << flush;
    o << "score          = " << score << "\n" << flush;
    //o << "type of scores = " << MP->Scoremodel[model[0].get_score_type()] << "\n" << flush;
    o << "\n" << flush;
    for (int i=0; i<steps+1; i++)               
    {
	score_sum+=sequence_of_scores[i];
	steps_x+=model[sequence_of_states[i]].get_letters_to_read_x();
	steps_y+=model[sequence_of_states[i]].get_letters_to_read_y();
	
	o << "   sequence_of_states[" << i << "] = " << sequence_of_states[i] 
	  << "\t | xsteps : " << steps_x 
	  << "\t | ysteps : " << steps_y
	  << "\t | delta_x : " << model[sequence_of_states[i]].get_letters_to_read_x()
	  << "\t | delta_y : " << model[sequence_of_states[i]].get_letters_to_read_y()
	  << "\t | score : " << sequence_of_scores[i] 
	  << "\t | cumul_score : " << score_sum << "\n" << flush;
    }
    o << "\n" << flush;
    return;
}

void Hmm:: output1(std::ostream &o, 
		   model_parameters* const MP,
		   const char* seq_name_x,
		   const char* seq_name_y) const
{
    // input check
    if(!seq_name_x)
    {
	cout<<"Error: Hmm class:: output1: the input seq_name_x is NULL."<<endl;
	return;
    }
    if((pair)&&(!seq_name_y))
    {
	cout<<"Error: Hmm class:: output1: the input seq_name_y is NULL."<<endl;
	return;
    }

    int label_index = -1;
    int state_index = -1;
    
    o<<seq_name_x;
    if(pair)
    {
	o<<"&"<<seq_name_y;
    } 
    o<<":\n";
    for(int i = 0; i<steps+1; i++)    
    {
	o<<i<<"\t"<<sequence_of_states[i]<<"\t";
	state_index = model[sequence_of_states[i]].get_number_of_state();
	o<<MP->get_State_Name(state_index)<<"\t";
	o<<"xdim="<<model[sequence_of_states[i]].get_letters_to_read_x()<<"\t";
	if(pair)
	{
	    o<<"ydim="<<model[sequence_of_states[i]].get_letters_to_read_y()<<"\t";
	}
	for(int j=0; j<model[sequence_of_states[i]].get_number_of_state_labels(); j++)
	{
	    for(int k=0; k<model[sequence_of_states[i]].get_state_labels_dim(j); k++)
	    {
		o<<model[sequence_of_states[i]].get_state_labels_type(j);
		if(k<model[sequence_of_states[i]].get_letters_to_read_x())
		{
		    label_index = model[sequence_of_states[i]].get_state_labels(j,k);
		    o<<"_x"<<k<<"="<<MP->get_Annotation_Label_name(j,label_index)<<"\t";
		}else{
		    label_index = model[sequence_of_states[i]].get_state_labels(j,k);
		    o<<"_y"<<k-model[sequence_of_states[i]].get_letters_to_read_x()
		     <<"="<<MP->get_Annotation_Label_name(j,label_index)<<"\t";
		}		
	    }
	}
	if(i!=steps+1)
	{
	    o<<"\n";
	}
    }  
    o<<"--------------------------------------------------------------\n";
    return;
}

void Hmm:: output2(std::ostream &o, 
		   model_parameters* const MP, 
		   const char* seq_name_x,
		   const int seq_length_x,
		   const char* seq_name_y,
		   const int seq_length_y) const
{
    // input check
    if(!seq_name_x)
    {
	cout<<"Error: Hmm class:: output1: the input seq_name_x is NULL."<<endl;
	return;
    }
    if(seq_length_x<=0)
    {
	cout<<"Error: Hmm class:: output1: the input seq_legnth_x("
	    <<seq_length_x<<") out of range."<<endl;
	return;
    }
    if(pair)	
    {
	if(!seq_name_y)
	{
	    cout<<"Error: Hmm class:: output1: the input seq_name_y is NULL."<<endl;
	    return;
	}
	if(seq_length_y<=0)
	{
	    cout<<"Error: Hmm class:: output1: the input seq_legnth_y("
		<<seq_length_y<<") out of range."<<endl;
	    return;
	}
    }
    
    //int check = 0;
    int label_index1 = -1;
    int label_index2 = -1;
    int number_of_state_label1 = -1;
    int number_of_state_label2 = -1;

    int state_path_index = 0;
    int last_state_path_index = 0;
    int xdim = 0;
    int xdim_pos = 0;
    int last_xdim_pos = 0;
    int seq_start_pos = 0;
    int seq_end_pos = 0;
    bool change = false;

    // print labels for seq_x
    while((seq_end_pos<seq_length_x-1)&&(state_path_index<steps+1))
    {
	xdim = model[sequence_of_states[state_path_index]].get_letters_to_read_x();
	xdim_pos = 0;
	
	while((xdim-xdim_pos)>0) // check state labels in a state
	{
	    if(xdim-xdim_pos>1)
	    {
		last_state_path_index = state_path_index;
		last_xdim_pos = xdim_pos;
		for(int i=0; i<model[sequence_of_states[state_path_index]].get_number_of_state_labels(); i++)
		{
		    label_index1 = model[sequence_of_states[state_path_index]].get_state_labels(i,xdim_pos);
		    label_index2 = model[sequence_of_states[state_path_index]].get_state_labels(i,xdim_pos+1);
		    if(label_index1!=label_index2)
		    {
			change = true;
			break;
		    }
		}
	    }else{
		// The position is on the last dimension of a state
		last_state_path_index = state_path_index;
		last_xdim_pos = xdim_pos;
		if(last_state_path_index>=steps)
		{
		    break;
		}
		while((model[sequence_of_states[state_path_index+1]].get_letters_to_read_x()==0)
		      &&(state_path_index<steps))
		{		    
		    state_path_index++;					    
		}
		if(state_path_index>=steps)
		{
		    break;
		}
		number_of_state_label1 = model[sequence_of_states[last_state_path_index]].get_number_of_state_labels();
		number_of_state_label2 = model[sequence_of_states[state_path_index+1]].get_number_of_state_labels();		
		if(number_of_state_label1!=number_of_state_label2)
		{
		    change = true;
		}else{
		    for(int i=0; i<number_of_state_label1; i++)
		    {
			label_index1 = model[sequence_of_states[last_state_path_index]].get_state_labels(i,xdim_pos);
			label_index2 = model[sequence_of_states[state_path_index+1]].get_state_labels(i,0);
			if(label_index1!=label_index2)
			{
			    change = true;
			    break;
			}
		    }
		}		
	    }
	    if(change)
	    {
		o<<seq_name_x<<"\t"<<MP->get_Model_Name()<<"\t"<<seq_start_pos+1<<"\t"<<seq_end_pos+1<<"\t";		
		for(int i =0; i<model[sequence_of_states[last_state_path_index]].get_number_of_state_labels();i++)
		{
		    o<<model[sequence_of_states[last_state_path_index]].get_state_labels_type(i);
		    label_index1 = model[sequence_of_states[last_state_path_index]].get_state_labels(i,last_xdim_pos);
		    o<<"="<<MP->get_Annotation_Label_name(i,label_index1)<<"\t";		
		}
		o<<"\n";
		seq_end_pos++;
		seq_start_pos = seq_end_pos;
		change = false;
	    }else{
		seq_end_pos++;
	    }
	    xdim_pos++;
	}	
	state_path_index++;      
    }    
    // print the last label set for sequence x   
    o<<seq_name_x<<"\t"<<MP->get_Model_Name()<<"\t"<<seq_start_pos+1<<"\t"<<seq_end_pos+1<<"\t";		
    for(int i =0; i<model[sequence_of_states[last_state_path_index]].get_number_of_state_labels();i++)
    {
	o<<model[sequence_of_states[last_state_path_index]].get_state_labels_type(i);
	label_index1 = model[sequence_of_states[last_state_path_index]].get_state_labels(i,last_xdim_pos);
	o<<"="<<MP->get_Annotation_Label_name(i,label_index1)<<"\t";		
    }
    o<<"\n";
       
    // print labels for seq y
    if(pair)
    {
	int ydim = 0;
	int ydim_pos = 0;
	int last_ydim_pos = 0;
	last_state_path_index = 0;
	state_path_index = 0;
	seq_start_pos = 0;
	seq_end_pos = 0;
	
	while((seq_end_pos<seq_length_y-1)&&(state_path_index<steps+1))
	{
	    ydim = model[sequence_of_states[state_path_index]].get_letters_to_read_y();
	    ydim_pos = 0;
	
	    while((ydim-ydim_pos)>0) // check state labels in a state
	    {
		//change = false;
		if(ydim-ydim_pos>1)
		{
		    last_state_path_index = state_path_index;
		    last_ydim_pos = ydim_pos;
		    xdim = model[sequence_of_states[state_path_index]].get_letters_to_read_x();
		    for(int i=0; i<model[sequence_of_states[state_path_index]].get_number_of_state_labels(); i++)
		    {
			label_index1 = model[sequence_of_states[state_path_index]].get_state_labels(i,ydim_pos+xdim);
			label_index2 = model[sequence_of_states[state_path_index]].get_state_labels(i,ydim_pos+xdim+1);
			if(label_index1!=label_index2)
			{
			    change = true;
			    break;
			}
		    }
		}else{
		    // The position is on the last dimension of a state
		    last_state_path_index = state_path_index;
		    last_ydim_pos = ydim_pos;
		    if(last_state_path_index>=steps)
		    {
			break;
		    }
		    while((model[sequence_of_states[state_path_index+1]].get_letters_to_read_y()==0)
			  &&(state_path_index<steps))
		    {		    
			state_path_index++;					    
		    }
		    if(state_path_index>=steps)
		    {
			break;
		    }
		    number_of_state_label1 = model[sequence_of_states[last_state_path_index]].get_number_of_state_labels();
		    number_of_state_label2 = model[sequence_of_states[state_path_index+1]].get_number_of_state_labels();		
		    if(number_of_state_label1!=number_of_state_label2)
		    {
			change = true;
		    }else{
			
			for(int i=0; i<number_of_state_label1; i++)
			{
			    xdim = model[sequence_of_states[last_state_path_index]].get_letters_to_read_x();

			    label_index1 = model[sequence_of_states[last_state_path_index]].get_state_labels(i,ydim_pos+xdim);

			    xdim = model[sequence_of_states[state_path_index+1]].get_letters_to_read_x();
			    label_index2 = model[sequence_of_states[state_path_index+1]].get_state_labels(i,xdim);
			    if(label_index1!=label_index2)
			    {
				change = true;
				break;
			    }
			}
		    }				    		    
		}
		if(change)
		{
		    xdim = model[sequence_of_states[last_state_path_index]].get_letters_to_read_x();
		    o<<seq_name_y<<"\t"<<MP->get_Model_Name()<<"\t"<<seq_start_pos+1<<"\t"<<seq_end_pos+1<<"\t";		
		    for(int i =0; i<model[sequence_of_states[last_state_path_index]].get_number_of_state_labels();i++)
		    {
			o<<model[sequence_of_states[last_state_path_index]].get_state_labels_type(i);
			    label_index1 = model[sequence_of_states[last_state_path_index]].get_state_labels(i,last_ydim_pos+xdim);
			    o<<"="<<MP->get_Annotation_Label_name(i,label_index1)<<"\t";		
		    }
		    o<<"\n";
		    seq_end_pos++;
		    seq_start_pos = seq_end_pos;
		    change = false;
		}else{
		    seq_end_pos++;
		}
		ydim_pos++;
	    }	
	    state_path_index++;      
	}    
	// print the last label set for sequence y

	xdim = model[sequence_of_states[last_state_path_index]].get_letters_to_read_x();
	
	o<<seq_name_y<<"\t"<<MP->get_Model_Name()<<"\t"<<seq_start_pos+1<<"\t"<<seq_end_pos+1<<"\t";		
	for(int i =0; i<model[sequence_of_states[last_state_path_index]].get_number_of_state_labels();i++)
	{
	    o<<model[sequence_of_states[last_state_path_index]].get_state_labels_type(i);
	    label_index1 = model[sequence_of_states[last_state_path_index]].get_state_labels(i,last_ydim_pos+xdim);
	    o<<"="<<MP->get_Annotation_Label_name(i,label_index1)<<"\t";		
	}    
	o<<"\n";
	
    }        
    
    return;
}

void Hmm:: transition_prob_output(std::ostream &o, 
				  TransitionProb* const TP) const
{

    int FTPsize = TP->get_FTPsize();

    int i = 0;

    for(i=0; i<FTPsize; i++)
    {
	o<<TP->get_FTP_tname()<<"."<<i<<"\t"
	 <<TP->get_FTP_name(i)<<"\t"
	 <<TP->get_FTP_prob(i);
	// add pseudo-count for FTP later
	o<<endl;
    }
    
    return;
} 

void Hmm:: XML_output(const char* xmlfile,
		      std::ostream &o,
		      model_parameters* const MP,
		      TransitionProb* const TP) const
{
    if(!xmlfile)
    {
	cout<<"ERROR: hmm class:: XML_output : input xmlfile is NULL."<<endl;
	return;
    }
    FILE* xml_file = fopen(xmlfile,"rt");
    
    if(!xml_file)
    {
	cout<<"ERROR: hmm class:: XML_output : can not open file : "
	    <<xmlfile<<" ."<<endl;
	return;
    }

    char* line = new char[Max_line_length];
    
    int NoOfItems = 0;
    char** Items = new char*[Max_number_of_items];
    for(int i = 0; i<Max_number_of_items; i++)
    {
	Items[i] = new char[Max_word_length];
	strcpy(Items[i]," ");
    }

    int check = 0;

    int fromindex = -1;
    int toindex = -1;

    while((!feof(xml_file))&&(!check))
    {
	fgets(line,Max_line_length-1,xml_file);

	if(get_tag(line,"Transitions"))
	{
	    // put the transitions back
	    o<<line;
	    while(1)
	    {
		fgets(line,Max_line_length-1,xml_file);

		if(get_tag(line,"/Transitions"))
		{
		    o<<line;
		    break;
		}
		else if(get_tag(line,"from"))
		{
		    // get from index
		    check+= splitstring(line,NoOfItems,&Items,'"');
		    if(check)
		    {
			cout<<"Error : hmm class :: XML_output :in splitstring : "<<line<<endl;
			break;
		    }
		    fromindex = convert_State_id_to_int(MP,Items[1]);
		    if(fromindex<0)
		    {
			cout<<"Error: hmm class :: XML_output error format for from state : "
			    <<Items[1]<<endl;
			break;
		    }		
    		    o<<line;
		}else if(get_tag(line,"to"))
		{
		    // get to index
		    check+= splitstring(line,NoOfItems,&Items,'"');
		    if(check)
		    {
			cout<<"Error : hmm class :: XML_output :in splitstring : "<<line<<endl;
			break;
		    }
		    if(!strcmp(Items[1],"All"))
		    {
			bool train = false;
			for(int i =0; i<TP->get_TTPsize(); i++)
			{
			    if(TP->get_TTP_from(i)==fromindex)			     			       
			    {
				train = true;
				if(TP->get_TTP_trained(i))
				{
				    // output
				    o<<Items[0]<<"\"S."<<TP->get_TTP_to(i)<<"\""
				     <<Items[2]<<"\""<<TP->get_TTP_score(i)
				     <<"\""<<Items[4];
				    if(NoOfItems>5) // specified pseudoprob
				    {
					o<<"\""<<Items[5]<<"\""<<Items[6]<<endl;		   
				    }
				    o<<endl;
				}else{
				    o<<Items[0]<<"\""<<Items[1]<<"\""
				     <<Items[2]<<"\""<<TP->get_TTP_score(i)
				     <<"\""<<Items[4]<<endl;
				    break;
				}
			    }			    
			}
			if(!train)
			{
			    o<<line;
			}
		    }else{
			toindex = convert_State_id_to_int(MP,Items[1]);
			if(toindex<0)
			{
			    cout<<"Error: hmm class :: XML_output error format for to state : "
				<<Items[1]<<endl;
			    break;
			}	
			// put the value back 
			for(int i =0; i<TP->get_TTPsize(); i++)
			{
			    if((TP->get_TTP_from(i)==fromindex)
			       &&(TP->get_TTP_to(i)==toindex))
			    {
				// output
				o<<Items[0]<<"\""<<Items[1]<<"\""
				 <<Items[2]<<"\""<<TP->get_TTP_score(i)
				 <<"\""<<Items[4];
				if(NoOfItems>5) // pseudoprob specified
				{
				    o<<"\""<<Items[5]<<"\""<<Items[6]<<endl;		   
				}
				o<<endl;
				break;
			    }			    
			}
		    }
		}else if(get_tag(line,"/from")){
		    o<<line;
		    continue;
		}		
	    }
	}else if(get_tag(line,"sequence_analysis"))
	{
	    break;
	}
	else{
	    o<<line;
	}	
    } 

    // remove memory
    if(Items)
    {
	for(int i =0; i<Max_number_of_items; i++)
	{
	    if(Items[i]) delete [] Items[i];
	    Items[i] = NULL;	    	    
	}
	delete [] Items;
    }
    Items = NULL;
    return;
}

void Hmm:: emission_prob_output(std::ostream &o, 
				model_parameters* const MP,
				EmissionProb* const EP) const
{

    int FEPsize = EP->get_FEPsize();

    int i = 0;
    long j = 0;   
    long number_of_emission = 0;
    const int max_dim = 100;
    int* emission = new int[max_dim];

    for(i=0; i<FEPsize; i++)
    {
	o<<EP->get_FEP_tname()<<"."<<i<<"\t";
	if(EP->get_FEP_name(i))
	{
	    o<<EP->get_FEP_name(i)<<"\t";
	}
	o<<EP->get_FEP_dim(i)<<endl;
	
	number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_dim(i)));
	int tmp_number_of_emission = 0;
	for(j=0; j<number_of_emission; j++)
	{
	    tmp_number_of_emission = j;
	    if(EP->get_FEP_prob(i,j)!=0)
	    {
		for(int k =0; k<EP->get_FEP_dim(i); k++)
		{
		    emission[k] = 0;
		}
		int count = 0;
		while(tmp_number_of_emission>0)
		{
		    emission[count]=(tmp_number_of_emission%alphabet);
		    tmp_number_of_emission = tmp_number_of_emission/alphabet;
		    count++;
		}
		for(int k = EP->get_FEP_dim(i)-1; k>=0; k--)
		{
		    o<<MP->get_Alphabet_name(emission[k]);
		}

		o<<"\t"<<EP->get_FEP_prob(i,j);				
		if(EP->get_FEP_pseudoprob(i,j)>0)
		{
		    o<<"\t"<<EP->get_FEP_pseudoprob(i,j);
		}
		if(j==number_of_emission-1)
		{
		    if(i<FEPsize-1)
		    {
			o<<endl;
		    }
		}else{
		    o<<endl;
		}

	    }
	}
	if(i!=FEPsize-1)
	{
	    o<<endl;
	}
    }    
    if(emission) delete [] emission;
    emission = NULL;
    return;
} 

int Hmm::condense_solution(int discard_long_solution)
{
    // discard_long_solution==1 => long solution will be deleted

    int check=0;
    
    if (sequence_of_states==NULL)
    {
	cout << "ERROR Hmm::condense_solution : sequence_of_states NULL.\n" << flush;
	check+=1;
    }
    
    if (sequence_of_scores==NULL)
    {
	cout << "ERROR Hmm::condense_solution : sequence_of_scores NULL.\n" << flush;
	check+=1;
    }
    
    if (sequence_of_xsteps==NULL)
    {
	cout << "ERROR Hmm::condense_solution : sequence_of_xsteps NULL.\n" << flush;
	check+=1;
    }
    
    if (sequence_of_ysteps==NULL)
    {
	cout << "ERROR Hmm::condense_solution : sequence_of_ysteps NULL.\n" << flush;
	check+=1;
    }
  
    if (steps<1)
    {
	cout << "ERROR Hmm::condense_solution : steps (" << steps << ") < 1\n" << flush;
	check+=1;
    }
    
    if ((score==0) || (score==Logzero))
    {
	cout << "ERROR Hmm::condense_solution : score (" << score << ") either 0 or Logzero.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

	int length_of_condensed_solution=0;
	
	{      
	    for (int i=1; i<steps+1; i++)
	    {
		if (sequence_of_states[i]!=sequence_of_states[i-1])
		{
		    length_of_condensed_solution+=1;

		}
	    }
	}

	condensed_steps=0;
	if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	condensed_sequence_of_scores=NULL;
	
	if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	condensed_sequence_of_states=NULL;
	
	if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	condensed_sequence_of_xsteps=NULL;
	
	if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	condensed_sequence_of_ysteps=NULL;      
	
	condensed_sequence_of_scores=new Score[length_of_condensed_solution+1];
	condensed_sequence_of_states=new int[length_of_condensed_solution+1];
	condensed_sequence_of_xsteps=new int[length_of_condensed_solution+1];
	condensed_sequence_of_ysteps=new int[length_of_condensed_solution+1];
	
	{
	    for (int i=0; i<length_of_condensed_solution+1; i++)
	    {
		condensed_sequence_of_scores[i]=0;
		condensed_sequence_of_states[i]=0;
		condensed_sequence_of_xsteps[i]=0;
		condensed_sequence_of_ysteps[i]=0;
	    }
	}
	
	int index=0;
	int xsteps=0;
	int ysteps=0;
	Score local_score=0;

	{      
	    for (int i=1; i<steps+1; i++)
	    {
		if (sequence_of_states[i]!=sequence_of_states[i-1])
		{
		    local_score+=sequence_of_scores[i-1];
		    
		    condensed_sequence_of_states[index]=sequence_of_states[i-1];
		    condensed_sequence_of_scores[index]=local_score;
		    condensed_sequence_of_xsteps[index]=xsteps;
		    condensed_sequence_of_ysteps[index]=ysteps;	    

		    local_score=0;;
		    xsteps=sequence_of_xsteps[i];
		    ysteps=sequence_of_ysteps[i];
		    index+=1;
		}
		else
		{
		    local_score+=sequence_of_scores[i-1];	      

		}
		
		if (i==steps)
		{
		    local_score+=sequence_of_scores[i];
		    
		    condensed_sequence_of_states[index]=sequence_of_states[i];
		    condensed_sequence_of_scores[index]=local_score;
		    condensed_sequence_of_xsteps[index]=xsteps;
		    condensed_sequence_of_ysteps[index]=ysteps;	    

		}
	    }
	}    
	
	condensed_steps=length_of_condensed_solution;

	Score sum_score=0;
	
	{
	    for (int i=0; i<condensed_steps+1; i++)
	    {

		sum_score+=condensed_sequence_of_scores[i];
	    }
	}
  
	if (abs(sum_score-score)>Max_deviation)
	{
	    cout << "ERROR Hmm::condense_solution: |sum_score (" << sum_score << ") - score (" << score 
		 << ")| = " << abs(sum_score-score) << " too large.\n" << flush;
	    
	    condensed_steps=0;
	    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	    condensed_sequence_of_scores=NULL;
	    
	    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	    condensed_sequence_of_states=NULL;
	    
	    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	    condensed_sequence_of_xsteps=NULL;
	    
	    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	    condensed_sequence_of_ysteps=NULL;      
	    
	    check+=1;
	}
	
	if (discard_long_solution==1)
	{
	    // note score is used both for the solution and the condensed solution 
	    // and it is therefore not initialised here
	    
	    steps=0;
	    
	    if (sequence_of_scores) delete [] sequence_of_scores;
	    sequence_of_scores=NULL;
	    
	    if (sequence_of_states) delete [] sequence_of_states;
	    sequence_of_states=NULL;
	    
	    if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	    sequence_of_xsteps=NULL;
	    
	    if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	    sequence_of_ysteps=NULL;
	}
    }
    return(check);
}

int Hmm::viterbi(//const int s, 
    const Sequence *x, const Sequence *y)
{
    // traceback is calculated, no pointers are involved
    
    int check=0;
        
    // check pairhmm
    
#ifndef _DOUBLE_STRANDED
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::viterbi : mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }
#endif
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::viterbi : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::viterbi : state 0 != Start \n"<<endl;
	cout << "model[0].get_letters_to_read : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::viterbi : last state of != End \n" <<endl;
	cout<<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }
    
    for (int i=0; i<number_of_states; i++)
    {
	if (model[i].get_alphabet()!= alphabet) 
	{
	    cout << "ERROR: Hmm::viterbi : alphabet of state i = " << i << ", alphabet = " 
		 << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet<< "\n" << flush;
	    check+=1;
	}
	
#ifndef _DOUBLE_STRANDED
	if (model[i].get_mirrored()!=mirrored) 
	{
	    cout << "ERROR: Hmm::viterbi : mirrored of state i = " << i << ", mirrored = " 
		 << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << mirrored << "\n" << flush;
	    check+=1;
	}
#endif

	if (model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::viterbi : number of states in state i = " << i << ", number of states = " 
		 << model[i].get_number_of_states() << " != number of states of Hmm = " 
		 <<number_of_states << "\n" << flush;
	    check+=1;
	
	}
	
	if ((i!=0) && (i!=(number_of_states-1)) && (model[i].get_letters_to_read()==0))
	{
	    cout << "ERROR: Hmm::viterbi : state i = " << i << " should emit\n"<<flush;
	    check+=1;
	}
	    
    }

  // check sequences
  
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL \n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::viterbi : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL \n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

	int shift_index = bitshift(alphabet);
 
	// a) for each state

	{
	    for (int i=0; i<number_of_states; i++)
	    {
		model[i].reset_viterbi_scores();
		
		model[i].set_dimensions_viterbi_scores(x->length()+1, y->length()+1);

	    }
	}

	// b) for Hmm
 
	steps=0;
	score=0;
 
	if (sequence_of_scores) delete [] sequence_of_scores;
	sequence_of_scores=NULL;
 
	if (sequence_of_states) delete [] sequence_of_states;
	sequence_of_states=NULL;
 
	if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	sequence_of_xsteps=NULL;
	
	if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	sequence_of_ysteps=NULL;

	condensed_steps=0;
	if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	condensed_sequence_of_scores=NULL;
	
	if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	condensed_sequence_of_states=NULL;
	
	if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	condensed_sequence_of_xsteps=NULL;
	
	if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	condensed_sequence_of_ysteps=NULL;      

	
	sequence_of_scores=new Score[x->length()+y->length()];
	sequence_of_states=new int[x->length()+y->length()];
	sequence_of_xsteps=new int[x->length()+y->length()];
	sequence_of_ysteps=new int[x->length()+y->length()];

	{	
	    for (int i=0; i< (x->length()+y->length()); i++)               
	    {
		sequence_of_scores[i]=Logzero;    
		sequence_of_states[i]=0;   
		sequence_of_xsteps[i]=0;   
		sequence_of_ysteps[i]=0;    
	    }
	}
	
	// variables which will only be needed to calculate alignment
 
	int i,j;
	int xsteps, ysteps;
	int deltax, deltay;
	
	int max_state, max_deltax, max_deltay;
	Score max_score;
 
	int test_state, test_deltax, test_deltay;
	Score test_score, pre_test_score;
 
	int states=number_of_states;
 
	for ( xsteps=0; xsteps<x->length()+1; xsteps++)
	{

	    for ( ysteps=0; ysteps<y->length()+1; ysteps++)
	    {

		if ( (xsteps==0) && (ysteps==0) )
		{
		    model[0].set_viterbi_score(0,0,0);

		}
		else
		{

		    for (int new_state=1; new_state<states-1; new_state++)
		    {		     
			max_score=Logzero;
			max_state=0;
			max_deltax=0;
			max_deltay=0;
 
			deltax=0;
			deltay=0;
 
			deltax=model[new_state].get_letters_to_read_x();
			deltay=model[new_state].get_letters_to_read_y();

			if ((xsteps>=deltax) && (ysteps>=deltay))
			{

			    int  linear_index_emission = 0;
			    int* indices = NULL;
			    indices = new int[model[new_state].get_number_of_dimensions_of_emission_scores()+1];
			    int count=0;
			    for (i=0; i<deltax; i++)
			    {
				indices[i+1]=x->letter(xsteps-deltax+i);

				count+=1;
			    }
			    for (j=count; j<deltay+count; j++)
			    {
				indices[j+1]=y->letter(ysteps-(deltay+count)+j);

			    }

			    int number_of_old_states = model[new_state].get_number_of_previous_states();

			    for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
			    {
				int old_state = model[new_state].get_number_of_previous_state(old_state_number);
				pre_test_score=Logzero;
				test_score=Logzero;
				test_state=0;
				test_deltax=0;
				test_deltay=0;

				indices[0] = old_state;

				linear_index_emission = indices[1];
				for (int p=2; p<model[new_state].get_number_of_dimensions_of_emission_scores()+1; p++)
				{
				    linear_index_emission = (linear_index_emission << shift_index) | indices[p];
				}
				
				int x_position = xsteps - deltax;
				int y_position = ysteps - deltay;

				pre_test_score=this->new_get_transition_score(NULL, 
									      old_state, new_state, 
									      x, x_position, y, y_position) +
				    this->model[new_state].get_emission_score(x, x_position, y, y_position,
									      linear_index_emission);
				
				if (pre_test_score > Logzero)
				{
				    test_state=old_state;
				    test_score=pre_test_score+
					model[old_state].get_viterbi_score(xsteps-deltax, ysteps-deltay);
				    test_deltax=deltax;
				    test_deltay=deltay;

				    if (test_score>max_score)
				    {

					max_score=test_score;
					max_state=test_state;
					max_deltax=test_deltax;
					max_deltay=test_deltay;

				    }
				}
			    } // for loop old_state
			    
			    if (indices) delete [] indices;
			    indices = NULL;

			    if ( max_score>Logzero )
			    {
				model[new_state].set_viterbi_score(xsteps, ysteps, max_score);			      

			    }
			} // if ((xsteps>=deltax) && (ysteps>=deltay))
			else
			{

			}
		    } // for loop new_state
		} // if not ( (xsteps==0) && (ysteps==0) )
	    } // for loop ysteps
	} // for loop xsteps

	xsteps=x->length();
	ysteps=y->length();
	
	max_score=Logzero;
	max_state=0;
	max_deltax=0;
	max_deltay=0;
 
	deltax=0;
	deltay=0;
	
	for (int state=0; state<states-1; state++)  // loop over all possible previous states except the End state
	{
	    pre_test_score=Logzero;
	    test_score=Logzero;
	    test_state=0;
	    test_deltax=0;
	    test_deltay=0;
	    
	    int x_position = xsteps - deltax;
	    int y_position = ysteps - deltay;

	    pre_test_score=this->new_get_transition_score(NULL, 
							  state, states-1, 
							  x, x_position, y, y_position);

	    if (pre_test_score > Logzero)
	    {
		test_state=state;
		test_score=pre_test_score+
		    model[state].get_viterbi_score(xsteps-deltax, ysteps-deltay);
		test_deltax=0;
		test_deltay=0;

		if (test_score>max_score)
		{

		    max_score=test_score;
		    max_state=test_state;
		    max_deltax=test_deltax;
		    max_deltay=test_deltay;

		}
	    }
	}
	if ( max_score>Logzero )
	{
	    model[states-1].set_viterbi_score(xsteps, ysteps, max_score);

	}

	score=Logzero;
	steps=0;
	test_state=0;
	test_score=Logzero;
	Score test_score_2=Logzero;

	max_score=Logzero;
	Score max_score_2=Logzero;
	max_state=0;
	deltax=0;
	deltay=0;
	
	int new_state=number_of_states-1;


	if (model[states-1].get_viterbi_score(xsteps,ysteps) > Logzero) // if a state path is possible
	{
	    score=model[states-1].get_viterbi_score(xsteps,ysteps);
	    
	    Score remaining_score=score;
	    
	    sequence_of_states[0]=new_state;
	    sequence_of_scores[0]=score;
	    
	    deltax=model[new_state].get_letters_to_read_x();
	    deltay=model[new_state].get_letters_to_read_y();

 
	    int new_state=states-1;
 
	    while ((new_state != 0) && (check == 0)) // as long as we haven't reached the start state
	    {
		max_score=Logzero;
		max_score_2=Logzero;
		max_state=0;

		int number_of_old_states = model[new_state].get_number_of_previous_states();

		for (int old_state_number=0; old_state_number<number_of_old_states; old_state_number++)
		{
		    int old_state = model[new_state].get_number_of_previous_state(old_state_number);
		    
		    if (model[new_state].get_transition_score(old_state)>Logzero) 
			// if the transition old_state -> new_state is allowed
		    {                                                     
			pre_test_score=Logzero;
			test_score=Logzero;
			test_score_2=Logzero;
			test_state=0;

			int* indices = NULL;
			indices = new int[model[new_state].get_number_of_dimensions_of_emission_scores()+1];
			for (int k=0; k<deltax; k++)
			{
			    indices[k+1]=x->letter(xsteps - (deltax-k));

			}
			for (int j=deltax; j<deltay+deltax; j++)
			{
			    indices[j+1]=y->letter(ysteps - (deltay+deltax-j));

			}
			indices[0]=old_state;
			test_state=old_state;
			test_score_2=model[old_state].get_viterbi_score(xsteps-deltax,ysteps-deltay);
			
			int x_position = xsteps - deltax;
			int y_position = ysteps - deltay;
			
			if ((deltax+deltay) > 0)
			{
			    int linear_index_emission = indices[1];
			    for (int p=2; p<model[new_state].get_number_of_dimensions_of_emission_scores()+1; p++)
			    {
				linear_index_emission = (linear_index_emission << shift_index) | indices[p];
			    }
			    
			    pre_test_score=this->new_get_transition_score(NULL, 
									  old_state, new_state,
									  x, x_position, 
									  y, y_position) +
				model[new_state].get_emission_score(x, x_position, y, y_position,
								    linear_index_emission);
			}
			else
			{
			    pre_test_score=this->new_get_transition_score(NULL, 
									  old_state, new_state,
									  x, x_position, y, y_position);
			}
			
			test_score = pre_test_score + test_score_2;

			if (test_score>max_score)
			{

			    max_score=test_score;
			    max_score_2=test_score_2;
			    max_state=test_state;

			}
			
			if (indices) delete [] indices;
			indices = NULL;
			
		    } // if transition old_state -> new_state is allowed
		} // for loop over old_state
		
		if ( max_score>Logzero )
		{
		    sequence_of_scores[steps]=remaining_score-max_score_2;
		    sequence_of_states[steps+1]=max_state;
		    sequence_of_scores[steps+1]=max_score_2;

		    new_state=max_state;
		    xsteps-=deltax;
		    ysteps-=deltay;
		    deltax=model[new_state].get_letters_to_read_x();
		    deltay=model[new_state].get_letters_to_read_y();
		    remaining_score=max_score_2;
		    steps++;

		}
		else // if no old_state had score > Logzero
		{
		    cout << "   ERROR Hmm::viterbi: traceback: all previous states have score Logzero.\n" << flush;
		    check+=1;
		}
	    } // while (new_state != 0)


	    if (check == 0)
	    { 
		sequence_of_states[steps]=new_state;
		sequence_of_scores[steps]=remaining_score;
  
		Score copy_score=0;
		int copy_state=0;
		
		int end_of_loop=static_cast<int>(ceil(static_cast<float>(steps)/2.));

		copy_score=sequence_of_scores[steps];
		copy_state=sequence_of_states[steps];
	      
		sequence_of_states[steps]=sequence_of_states[0];
		sequence_of_scores[steps]=sequence_of_scores[0];

		{	      
		    for (int i=1; i<end_of_loop; i++) // put sequence_of_scores/states in right order
		    {
			sequence_of_states[0]=sequence_of_states[i];
			sequence_of_scores[0]=sequence_of_scores[i];
			
			sequence_of_states[i]=sequence_of_states[steps-i];
			sequence_of_scores[i]=sequence_of_scores[steps-i];
			
			sequence_of_states[steps-i]=sequence_of_states[0];
			sequence_of_scores[steps-i]=sequence_of_scores[0];
		    }	    
		}
		sequence_of_states[0]=copy_state;
		sequence_of_scores[0]=copy_score;  
		
		Score score_sum=0;
		int steps_x=0;
		int steps_y=0;

		{	
		    for (int i=0; i<steps+1; i++)               
		    {
			score_sum+=sequence_of_scores[i];
			steps_x+=model[sequence_of_states[i]].get_letters_to_read_x();
			steps_y+=model[sequence_of_states[i]].get_letters_to_read_y();
		  
			sequence_of_xsteps[i]=steps_x;
			sequence_of_ysteps[i]=steps_y;

		    }
		}

		// check consistency of solution

		check+=check_consistency_of_solution(x, y);

	    } // if check == 0
	} // if state path is possible, i.e. if score>Logzero
	else
	{
	    cout << "NOTE: Hmm::viterbi: no state path can be found.\n" << flush;
	    check+=1;
	}
	
	{
	    for (int i=0; i<number_of_states; i++)
	    {
		model[i].reset_viterbi_scores();
	    }
	}
    } // if alignment can be started

    return(check);
}

int Hmm::new_viterbi(//const int s, 
    const Sequence *x, const Sequence *y)
{
    // traceback is calculated, no pointers are involved
    
    int check=0;

    // check pairhmm

    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::new_viterbi : mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::new_viterbi : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::new_viterbi : state 0 read != Start \n" << flush;
	cout << "model[0].get_letters_to_read() = "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout <<"ERROR: Hmm::new_viterbi : last state of type != End \n" << flush;
	cout <<"model["<<number_of_states-1<<"].get_letters_to_read() : "<<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (int i=0; i<number_of_states; i++)
    {
	if (model[i].get_alphabet()!=alphabet) 
	{
	    cout << "ERROR: Hmm::new_viterbi : alphabet of state i = " << i << ", alphabet = " 
		 << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (model[i].get_mirrored()!=mirrored) 
	{
	    cout << "ERROR: Hmm::new_viterbi : mirrored of state i = " << i << ", mirrored = " 
		 << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << mirrored << "\n" << flush;
	    check+=1;
	}
	if (model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::new_viterbi : number of states in state i = " << i << ", number of states = " 
		 << model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}

	if ((i!=0) && (i!=(number_of_states-1))  && (model[i].get_letters_to_read()==0 ))
	{
	    cout << "ERROR: Hmm::new_viterbi : state i = " << i << " should emit! \n" << flush;
	    check+=1;
	}
	
    }
    
    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::new_viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL \n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::new_viterbi : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL \n" << flush;
	check+=1;
    }

    if (check==0)
    {

	steps=0; 
	score=0;
	
	if (sequence_of_scores) delete [] sequence_of_scores;
	sequence_of_scores=NULL;
      
	if (sequence_of_states) delete [] sequence_of_states;
	sequence_of_states=NULL;
	
	if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	sequence_of_xsteps=NULL;
	
	if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	sequence_of_ysteps=NULL;
	
	condensed_steps=0;
	if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	condensed_sequence_of_scores=NULL;
	
	if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	condensed_sequence_of_states=NULL;
      
	if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	condensed_sequence_of_xsteps=NULL;
	
	if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	condensed_sequence_of_ysteps=NULL;      
	
	sequence_of_scores=new Score[x->length()+y->length()];
	sequence_of_states=new int[x->length()+y->length()];
	sequence_of_xsteps=new int[x->length()+y->length()];
	sequence_of_ysteps=new int[x->length()+y->length()];
    
	for (int i=0; i< (x->length()+y->length()); i++)               
	{
	    sequence_of_states[i]=0;   
	    sequence_of_scores[i]=Logzero;    
	    sequence_of_xsteps[i]=-1;   
	    sequence_of_ysteps[i]=-1;    
	}

	// use viterbi_rectangle function
      
	int x_start=0;
	int x_end=x->length();
	int y_start=0;
	int y_end=y->length();
	
	int start_states[1]={0};
	const int number_of_start_states=1;
	int end_state = number_of_states-1;
	
	int x_start_traceback     = 0;
	int y_start_traceback     = 0;
	int start_state_traceback = 0;
	
	check+=viterbi_rectangle(NULL, // do not use external rectangle
				 this,
				 x, x_start, x_end,
				 y, y_start, y_end,
				 start_states, end_state,
				 number_of_start_states, 
				 0,0, // x_margin, y_margin 
				 0,   // offset for add_local_to_global_solution
			       // output values (won't be used any further in this function)
				 &x_start_traceback, 
				 &y_start_traceback, 
				 &start_state_traceback);
	
	check+=check_consistency_of_solution(x, y);

    } // if alignment can be started

    return(check);
}

int Hmm::Hirschberg_viterbi(//const int score_type, 
				const Hmm* mirror,			  
				const Sequence *x, const Sequence *y,
				const int max_area)
{
    int check=0;
  
    // check Hmm
  
#ifndef _DOUBLE_STRANDED
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored (" << mirrored << ") of Hmm must be 0.\n" 
	     << flush;
	check+=1;
    }
#endif
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : state 0 of type != Start\n" << flush;
	cout << "model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : last state of type != End \n" << flush;
	cout << "model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : alphabet of state i = " << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
#ifndef _DOUBLE_STRANDED
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored of state i = " << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
#endif
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1)) && (model[i].get_letters_to_read()==0))
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : state i = " << i << " should emit !\n" << flush;
		check+=1;
	    }
	   
	}
    }
    
    // check mirrored Hmm
    
#ifndef _DOUBLE_STRANDED
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored (" << mirror->get_mirrored()
	     << ") of mirrored Hmm must be 1.\n" << flush;
	check+=1;
    }
#endif
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : state 0 of type != Start \n" << flush;
	cout << "mirror->model[0].get_letters_to_read() : "<<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : last state of type != End \n" << flush;
	cout << "mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
#ifndef _DOUBLE_STRANDED
	    if (mirror->model[i].get_mirrored()!= mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
#endif
	    if (mirror->model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	   
	    if ((i!=0) && (i!=(number_of_states-1)) && (mirror->model[i].get_letters_to_read()==0))
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    // check that this Hmm and the mirror Hmm belong together
    
    if (number_of_states != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (int i=1; i<number_of_states-1; i++)
	{
	    if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : letters_to_read of model[" << i << "] (" 
		     << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		     << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }

  // check that max_area value is sensible
    
  // determine min_strip_width
    
    int min_strip_width=0;
    
    for (int i=1; i<number_of_states-1; i++) // loop over all states except Start and End state
    {
	if (max(model[i].get_letters_to_read_x(), model[i].get_letters_to_read_y())>min_strip_width)
	{
	    min_strip_width=max(model[i].get_letters_to_read_x(), model[i].get_letters_to_read_y());

	}
    }
  
  // check that max_area exceeds some minimal possible value

    if (max_area < ((y->length()+1) * 2 * (min_strip_width+1)))
    {
	cout << "ERROR: Hmm:: Hirschberg_viterbi : max_area (" << max_area
	     << ") < (y->length()+1) (" << y->length()+1 << ") * 2 * (min_strip_width+1) ("
	     << min_strip_width+1 << ") = " << (y->length()+1) * 2 * (min_strip_width+1) 
	     << ". Increase max_area value to Start the calculation.\n" << flush;      
	check+=1;
    }
    
    if (max_area<1)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : max_area (" << max_area
	     << ") < 1. Choose a larger value.\n" << flush;
	check+=1;
    }

    // check sequences

    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

      // note: in order to use the function add_local_solution, 
      //       steps has to be initialised with -1 as well as 
      //
      //	    sequence_of_states[i]=-1;   
      //	    sequence_of_scores[i]=Logzero;    
      //	    sequence_of_xsteps[i]=-1;   
      //	    sequence_of_ysteps[i]=-1;    

	steps=-1;
	score=0;
	
	if (sequence_of_scores) delete [] sequence_of_scores;
	sequence_of_scores=NULL;
	
	if (sequence_of_states) delete [] sequence_of_states;
	sequence_of_states=NULL;
	
	if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	sequence_of_xsteps=NULL;
	
	if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	sequence_of_ysteps=NULL;
	
	condensed_steps=0;
	if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	condensed_sequence_of_scores=NULL;
	
	if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	condensed_sequence_of_states=NULL;
	
	if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	condensed_sequence_of_xsteps=NULL;
	
	if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	condensed_sequence_of_ysteps=NULL;      
	
	sequence_of_scores=new Score[x->length()+y->length()];
	sequence_of_states=new int[x->length()+y->length()];
	sequence_of_xsteps=new int[x->length()+y->length()];
	sequence_of_ysteps=new int[x->length()+y->length()];

	for (int i=0; i< (x->length()+y->length()); i++)               
	{
	    sequence_of_states[i]=-1;   
	    sequence_of_scores[i]=Logzero;    
	    sequence_of_xsteps[i]=-1;   
	    sequence_of_ysteps[i]=-1;    
	}
	
	// deciding whether to use viterbi_rectangle or memory_viterbi

	if ((x->length()+1)*(y->length()+1) < max_area) // use viterbi if area < max_area
	{

	    check+=viterbi(x, y); 
	}
	else // use Hirschberg algorithm if area > max_area
	{

	    // check value of max_area
	    
	    if ( (2*(min_strip_width+1)*(y->length()+1)) < max_area)
	    {
		// use memory_viterbi
		
		const int start_states[1] = {0};
		const int number_of_start_states = 1;
		const int end_state = number_of_states-1;
		
		int x_start_traceback     = 0;
		int y_start_traceback     = 0;
		int start_state_traceback = 0;
		
		check += new_memory_viterbi(// input variables
		    NULL, // = viterbi_strip
		    mirror,
		    x, 0, x->length(),
		    y, 0, y->length(),
		    start_states, end_state,
		    number_of_start_states,
		    0, // = x_margin
		    0, // = y_margin
		    0, // = offset
		    min_strip_width,
		    max_area,
		    // output variables (here not used)
		    &x_start_traceback,    
		    &y_start_traceback,    
		    &start_state_traceback);
		
		// concatenate local state and score paths to global path

		
		if (check != 0) {
		    cout << "ERROR : Hmm::Hirschberg_viterbi : error occurred in function "
			 << "Hmm::new_memory_viterbi.\n" << flush;
		}
	    }
	    else 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi : value of max_area (" << max_area 
		     << ") =< 2 * min_strip_width+1 (" << min_strip_width+1 
		     << ") * length of sequence y +1 (" << y->length()+1 
		     << ") = " << 2*(min_strip_width+1)*(y->length()+1) 
		     << " too small to Start Hirschberg Viterbi. Increase the value of max_area.\n" << flush;
		check+=1;
	    }
	}      

	check+=check_consistency_of_solution(x, y); 
    } 
    return(check);
}

int Hmm::Hirschberg_viterbi_rectangle(Score*** viterbi_strip,
					  const Hmm *mirror,
					  const Sequence *x, const int x_start, const int x_end,
					  const Sequence *y, const int y_start, const int y_end, 
					  const int* start_states, const int end_state,
					  const int number_of_start_states,
					  const int x_margin, const int y_margin,
					  const int offset,
					  const int max_area,
					  // output values
					  int* x_start_traceback, 
					  int* y_start_traceback, 
					  int* start_state_traceback)
{
    // function returns 0 if o.k., else not o.k.
    //
    // note: this function can either be used 
    //
    //       1.) with start_states and an end_state (in which case
    //       viterbi_strip == NULL and x_margin and y_margin == 0). Memory
    //       for viterbi_strip is allocated and deleted within function.
    //       
    //       2.) with a pre-initialised viterbi_strip (in which case start_states
    //        == NULL and end_state == 0 and number_of start_states and number_of_end_states == 0).
    //       Memory for viterbi_strip is allocated and deleted outside function.

    int check=0;
    int i = 0;
    
    // check whether pre-initialised viterbi_rectangle or start_ and end_states shall be used
    
    if ((viterbi_strip == NULL)    &&
	(number_of_start_states == 0))
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: viterbi_strip (" << viterbi_strip
	     << ") = NULL and number_of_start_states (" << number_of_start_states 
	     << ") == 0. Use this function either with a pre-initialised viterbi_strip != NULL "
	   << "or with nonempty sets of start_ and end_states.\n" << flush;
	check+=1;
    }
    else if ((viterbi_strip != NULL)    &&
	     (number_of_start_states != 0))
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: viterbi_strip (" << viterbi_strip
	     << ") != NULL and number_of_start_states (" << number_of_start_states 
	     << ") != 0. Use this function either with a pre-initialised viterbi_strip != NULL "
	     << "or with nonempty sets of start_ and end_states.\n" << flush;
	check+=1;
    }
    else // if input values are not completely inconsistent
    {
	if ((viterbi_strip      == NULL)   &&
	    (number_of_start_states != 0))
	{
	    // if function shall be used with start_ and end_states
	    
	    // check Start states
	    
	    if (number_of_start_states<1)
	    {
		cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: number_of_start_sites (" 
		     << number_of_start_states << ") <1.\n" << flush;
		check+=1;
	    }
	    else
	    {
		if ((x_start==0) && (y_start==0))
		{
		    if ((number_of_start_states!=1) || (start_states[0]!=0))
		    {
			cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: or "
			     << "x_start=y_start=0 only "
			     << "one start_state is possible, the Start state (i.e. the number_of_start_states ("
			     << number_of_start_states << ") has to be 1 and start_states[0] ("
			     << start_states[0] << ") has to be 0).\n" << flush;
			check+=1;
		    }
		}
		else // if we don't start at (0,0)
		{
		    // check that each start state lies in the allowed range (1...n_of_states-2)
		    // i.e. Start and End states are not allowed 
	      
		    for (i=0; i<number_of_start_states; i++)
		    {
			if ((start_states[i]<1) || (start_states[i]>this->get_number_of_states()-2))
			{
			    cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: start_states[" << i 
				 << "] (" << start_states[i] << ") out of range [1, " 
				 << this->get_number_of_states()-2
				 << "].\n" << flush;
			    check+=1;
			}
		    }
		}
	    }
	    
	    // check end states
	    
	    if ((x_end==x->length()) && (y_end==y->length()))
	    {
		if (end_state != this->get_number_of_states()-1)
		{
		    cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: end_state (" 
			 << end_state << ") for x_end (" << x_end << ") = x->length() (" << x->length() 
			 << ") and y_end (" << y_end << ") = y->length() (" << y->length() 
			 << ") has to be End state (" << this->get_number_of_states()-1 << ").\n" << flush;
		    check+=1;
		}
	    }
	    else // if we don't stop at ends of both sequences
	    {
		// check that end state lies in the allowed range (1...n_of_states-2)
		// i.e. Start and End states are not allowed
      
		if ((end_state<1) || (end_state>this->get_number_of_states()-2))
		{
		    cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: end_state "
			 << end_state << ") out of range [1, " << this->get_number_of_states()-2 << "].\n" << flush;
		    check+=1;
		}
	    }
	}

	if ((viterbi_strip      != NULL)   &&
	    (number_of_start_states == 0))
	{
	    // if function shall be used with pre-initialised viterbi_strip
	    
	    // check values of x_margin and y_margin
	    
	    if ((x_start==0) && (x_margin!=0))
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : x_margin (" << x_margin 
		     << ") has to be zero for x_start (" 
		     << x_start << ") = 0.\n" << flush;
		check+=1;
	    }
	  
	    if ((y_start==0) && (y_margin!=0))
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : y_margin (" << y_margin 
		     << ") has to be zero for y_start (" 
		     << y_start << ") = 0.\n" << flush;
		check+=1;
	    }
	    
	    if (x_margin < 0)
	    {
		cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: x_margin (" << x_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	  
	    if (y_margin < 0)
	    {
		cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: y_margin (" << y_margin
		     << ") < 0\n" << flush;
		check+=1;
	    }
	  
	    // check that x_margin and y_margin are compatible with values of x_start, x_end and y_start and y_end
	    
	    if ((x_end-x_start+1) < x_margin)
	    {
		cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: x_margin (" << x_margin
		     << ") > (x_end (" << x_end << ") - x_start (" << x_start << ") + 1) = " 
		     << x_end-x_start+1 << "\n" << flush;
		check+=1;
	    }
	    
	    if ((y_end-y_start+1) < y_margin)
	    {
		cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: y_margin (" << y_margin
		     << ") > (y_end (" << y_end << ") - y_start (" << y_start << ") + 1) = " 
		     << y_end-y_start+1 << "\n" << flush;
		check+=1;
	    }
	}
    }
    
    // check this Hmm
    
    if (this->get_mirrored() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : mirrored (" 
	     << this->get_mirrored() << ") must be 0.\n" << flush;
	check+=1;
    }
    if (this->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : number of states = " 
	  << this->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (this->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : state 0  != Start \n" << flush;
	cout <<"this->model[0].get_letters_to_read() : "<<this->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (this->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : last state  != End \n" << flush;
	cout <<"this->model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<this->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (i=0; i<number_of_states; i++)
    {
	if (this->model[i].get_alphabet()!=alphabet) 
	{
	    cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << this->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (this->model[i].get_mirrored()!=this->get_mirrored()) 
	{
	    cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << this->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << this->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (this->model[i].get_number_of_states() != number_of_states)
	{
	    cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << this->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}

	if ((i!=0) && (i!=(number_of_states-1))  && (this->model[i].get_letters_to_read()==0))
	{
	    cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : Hmm : state i = " 
		 << i << " should emit!\n" << flush;
	    check+=1;
	}

    }
    
    // check that this Hmm and the mirror Hmm belong together
    
    if (this->get_number_of_states() != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : number_of_states of this Hmm (" 
	     << this->get_number_of_states()
	     << ") != number_of_states of mirrored Hmm (" 
	     << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (i=1; i<number_of_states-1; i++)
	{
	    if (this->model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : letters_to_read of model[" << i 
		     << "] in this Hmm (" 
		     << this->model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		  << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored (" << mirror->get_mirrored()
	     << ") of mirrored Hmm must be 1.\n" << flush;
	check+=1;
    }
    
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : state 0  != Start \n" << flush;
	cout <<"mirror->model[0].get_letters_to_read() : "<<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : last state  != End \n" << flush;
	cout << "mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	     << mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=this->model[i].get_alphabet()) 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of this Hmm = " 
		     << this->model[i].get_alphabet() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!= mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != this->model[i].get_number_of_states())
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of this Hmm = " 
		     << this->model[i].get_number_of_states() << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1))&&(mirror->model[i].get_letters_to_read() ==0))
	    {
		cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }

    // check that max_area value is sensible
    
    // determine min_strip_width
    
    int min_strip_width=0;
    
    for (i=1; i<this->get_number_of_states()-1; i++) // loop over all states except Start and End state
    {
	if (max(this->model[i].get_letters_to_read_x(), this->model[i].get_letters_to_read_y())>min_strip_width)
	{
	    min_strip_width=max(this->model[i].get_letters_to_read_x(), this->model[i].get_letters_to_read_y());

	}
    }

    // check that max_area exceeds some minimal possible value

    if (max_area < ((abs(y_end - y_start) +1 +1) * 2 * (min_strip_width+1)))
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : max_area (" << max_area
	     << ") < (abs(y_end - y_start) +1 +1) (" 
	     << abs(y_end - y_start) +1 +1 << ") * 2 * (min_strip_width+1) ("
	     << min_strip_width+1 << ") = " << (abs(y_end - y_start) +1 +1) * 2 * (min_strip_width+1) 
	     << ". Increase max_area value to Start the calculation.\n" << flush;      
	check+=1;
    }
    
    if (max_area<1)
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : max_area (" << max_area
	     << ") < 1. Choose a larger value.\n" << flush;
	check+=1;
    }

    // check sequences
    
    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL.\n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::Hirschberg_viterbi_rectangle : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL.\n" << flush;
	check+=1;
    }
    
    // check x_start, x_end and y_start, y_end values
    
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    if ((offset != 0) && (offset != -1))
    {
	cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: offset (" << offset << ") must be either 0 or -1.\n" << flush;
	check+=1;
    }
    
    if (check==0)
    {

	// initialise output variables
	
	(*x_start_traceback)     = 0;
	(*y_start_traceback)     = 0;
	(*start_state_traceback) = 0;
	
	check += new_memory_viterbi(// input variables
	    viterbi_strip,
				  mirror,
	    x, x_start, x_end,
	    y, y_start, y_end,
	    start_states, end_state,
	    number_of_start_states,
	    x_margin, y_margin,
	    offset,
	    min_strip_width,
	    max_area,
	    // output variables
	    x_start_traceback,    
	    y_start_traceback,    
	    start_state_traceback);

	if (check != 0) {
	    cout << "ERROR Hmm::Hirschberg_viterbi_rectangle: error occurred in function "
		 << "Hmm::new_memory_viterbi.\n" << flush;
	}
    }
    if (check != 0) {

	// initialise output variables
	
	(*x_start_traceback)     = 0;
	(*y_start_traceback)     = 0;
	(*start_state_traceback) = 0;
    }
    return(check);
}

int Hmm::heuristic_viterbi(const Hmm* mirror, // needed for internal calls to Hirschberg
			       const Sequence *x, const int* x_coordinates,
			       const Sequence *y, const int* y_coordinates,
			       const int number_of_xy_pairs,
			       const int x_margin,
			       const int y_margin,
			       const int max_area)
{
    // note : if, during the calculation, it turns out that there
    //        is no single state path, the calculation is aborted and the return value check is != 0
    //
    //     check that function also works on one rectangle (i.e. number_of_xy_pairs == 0)
    //     check that at least one of the values in the copied rectangle is not Logzero
    //     otherwise stop and declare that no solution could be found

    int check=0;
    int min_strip_width=0;
    
    {  
	for (int i=1; i<number_of_states-1; i++) // loop over all states except Start and End state
	{
	    if (max(model[i].get_letters_to_read_x(), model[i].get_letters_to_read_y())>min_strip_width)
	    {
		min_strip_width=max(model[i].get_letters_to_read_x(), model[i].get_letters_to_read_y());

	    }
	}
    }

    // check this Hmm
    
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : mirrored (" << mirrored 
	     << ") of Hmm must be 0.\n" << flush;
	check+=1;
    }  
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : state 0 != Start \n" << flush;
	cout << "model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : last state of type != End \n" << flush;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : alphabet of state i = " << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : mirrored of state i = " << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1))&& (model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : state i = " << i << " should emit! \n" ;
		    
		check+=1;
	    }
	   
	}
    }
    // check mirrored Hmm

    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : mirrored (" 
	     << mirror->get_mirrored() << ") must be 1.\n" << flush;
	check+=1;}
    
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : state 0 != Start \n" << flush;
	cout <<"mirror->model[0].get_letters_to_read() : "<<mirror->model[0].get_letters_to_read()<<endl;	
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : last state != End \n" << flush;
	cout <<"mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!=mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	  
	    if ((i!=0) && (i!=(number_of_states-1))&& (mirror->model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : mirrored Hmm : state i = " 
		<< i << " should emit!\n" << flush;
		check+=1;
	    }
	   
	}
    }
    // check that this Hmm and the mirror Hmm belong together
    
    if (number_of_states != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::heuristic_viterbi : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (int i=1; i<number_of_states-1; i++)
	{
	    if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : letters_to_read of  model[" << i << "] (" 
		     << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		     << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    
    // check that values of x_margin and y_margin exceed some minimal value
    
    if (x_margin < (min_strip_width+1))
    {
	cout << "ERROR: Hmm::heuristic_viterbi : x_margin (" << x_margin 
	     << ") < (min_strip_width+1) (" << min_strip_width+1 << ")\n" << flush;
	check+=1;
    }
    
    if (y_margin < (min_strip_width+1))
    {
	cout << "ERROR: Hmm::heuristic_viterbi : y_margin (" << y_margin 
	   << ") < (min_strip_width+1) (" << min_strip_width+1 << ")\n" << flush;
	check+=1;
    }

    // check that x_margin and y_margin are larger than some desired value 
    // note : the desired x_margin and y_margin values are large enough to accomodate
    //        at least three of the states with maximum number of reads per sequence (this number
    //        is equal to min_strip_width) next to eachother 
    
    if ( abs(floor(static_cast<float>(x_margin)/2.) - ceil(static_cast<float>(x_margin)/2.)) == 0) // if x_margin is even
    {
	if (x_margin < (3*min_strip_width+2))
	{
	    cout << "ERROR: Hmm::heuristic_viterbi : x_margin (" << x_margin 
		 << ") is even. x_margin < (3 * min_strip_width (" << min_strip_width << ") + 1) (" 
		 << 3*min_strip_width+1 
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    else // if x_margin is odd
    {
	if (x_margin < (3*min_strip_width))
	{
	    cout << "ERROR: Hmm::heuristic_viterbi : x_margin (" << x_margin 
		 << ") is odd. x_margin < (3 * min_strip_width (" << min_strip_width << ")) (" << 3*min_strip_width
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    
    if ( abs(floor(static_cast<float>(y_margin)/2.) - ceil(static_cast<float>(y_margin)/2.)) == 0) // if y_margin is even
    {
	if (y_margin < (3*min_strip_width+2))
	{
	    cout << "ERROR: Hmm::heuristic_viterbi : y_margin (" << y_margin 
		 << ") is even. y_margin < (3 * min_strip_width (" << min_strip_width << ") + 1) (" 
		 << 3*min_strip_width+1 
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    else // if y_margin is odd
    {
	if (y_margin < (3*min_strip_width))
	{
	    cout << "ERROR: Hmm::heuristic_viterbi : y_margin (" << y_margin 
		 << ") is odd. y_margin < (3 * min_strip_width (" << min_strip_width << ")) (" << 3*min_strip_width
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    
  
    // check Hmm
  
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : number of states = " << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : state 0 != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : last state != End \n" << flush;
	cout<<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : alphabet of state i = " << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : mirrored of state i = " << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		<< number_of_states << "\n" << flush;
		check+=1;
	    }
	   
	    if ((i!=0) && (i!=(number_of_states-1)) && (model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : state i = " << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    // check that max_area value is sensible
    
    if (max_area<1)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : max_area (" << max_area
	     << ") < 1.\n" << flush;
	check+=1;
    }
  
    // check the number of (x,y) constraints, their order and the size of rectangles they span

    if (number_of_xy_pairs<1)
    {
	cout << "ERROR: Hmm::heuristic_viterbi : number_of_xy_pairs (" << number_of_xy_pairs
	     << ") < 1.\n" << flush;
	check+=1;
    }
    else
    {
	for (int i=0; i<number_of_xy_pairs; i++)
	{
	    // check x_margin and y_margin and their compatibility with x_start, x_end
	    // y_start and y_end
	    
	    if (i==0) // check rectangle from (0,0) to first (x,y) pair
	    {
		if (((x_coordinates[i]-0+1) < x_margin+1) ||  ((y_coordinates[i]-0+1) < y_margin+1) )
		{
		    cout << "ERROR: Hmm::heuristic_viterbi : area of rectangle between (0,0) and (x_coordinates[" << i 
			 << "] (" << x_coordinates[i] << "), y_coordinates[" << i << "] (" << y_coordinates[i] 
			 << ")) cannot accumodate area of (x_margin (" << x_margin << ") + 1) * (y_margin (" << y_margin << ") + 1).\n" << flush;
		    check+=1;
		}
	    }
	    else if ((i>0) && (i<number_of_xy_pairs-1)) // check intermediate rectangle from (x,y) pair i to (x,y) pair i+1
	    {
		if ( ((x_coordinates[i+1]-x_coordinates[i]+1) < x_margin+1) ||  ((y_coordinates[i+1]-y_coordinates[i]+1) < y_margin+1) )
		{
		    cout << "ERROR: Hmm::heuristic_viterbi : area of rectangle between (x_coordinates[" << i << "] ("
			 << x_coordinates[i] << "), y_coordinates[" << i << "] (" << y_coordinates[i] << ")) and (x_coordinates[" << i+1
			 << "] (" << x_coordinates[i+1] << "), y_coordinates[" << i+1 << "] (" << y_coordinates[i+1] 
			 << ")) cannot accumodate area of (x_margin (" << x_margin << ") + 1) * (y_margin (" << y_margin << ") + 1).\n" << flush;
		    check+=1;
		}
	    }
	    else // check rectangle from last (x,y) pair to (x->length(), y->length())
	    {
		if ( ((x->length()-x_coordinates[i]+1) < x_margin+1) ||  ((y->length()-y_coordinates[i]+1) < y_margin+1) )
		{
		    cout << "ERROR: Hmm::heuristic_viterbi : area of rectangle between (x_coordinates[" << i << "] ("
			 << x_coordinates[i] << "), y_coordinates[" << i << "] (" << y_coordinates[i] << ")) and (x->length() (" << x->length()
			 << "), y->length() (" << y->length()
			 << ")) cannot accumodate area of (x_margin (" << x_margin << ") + 1) * (y_margin (" << y_margin << ") + 1).\n" << flush;
		    check+=1;
		}
	    }
	    
	    // check that no constaint is equal to 0 or x->length() or y->length(), respectively
	    
	    if (x_coordinates[i]==0 || y_coordinates[i]==0 ||
		x_coordinates[i]==x->length() || y_coordinates[i]==y->length())
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : x_coordinates[" << i << "] (" << x_coordinates[i]
		     << ") = 0 or y_coordinates[" << i << "] (" << y_coordinates[i] 
		     << ") = 0 or x_coordinates[" << i << "] (" << x_coordinates[i] << ") == x->length() ("
		     << x->length() << ") or y_coordinates[" << i << "] (" << y_coordinates[i] 
		     << ") == y->length() (" << y->length() << ")\n" << flush;
		check+=1;
	    }
	    
	    // check that 0 < x_i < x->length() and 0 < y_i < y->length()
	    
	    if ((x_coordinates[i] <= 0  ||  x_coordinates[i] >= x->length()) ||
		(y_coordinates[i] <= 0  ||  y_coordinates[i] >= y->length()))
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : x_coordinates[" << i << "] (" << x_coordinates[i]
		     << ") out of range [1, " << x->length()-1 
		     << "] or y_coordinates[" << i << "] (" << y_coordinates[i]
		     << ") out of range [1, " << y->length()-1 << "].\n" << flush;
		check+=1;
	    }
	    
	    if (i<(number_of_xy_pairs-1))
	    {
		// check that x_i < x_i+1 and y_i < y_i+1
		
		if ((x_coordinates[i]>=x_coordinates[i+1]) ||
		    (y_coordinates[i]>=y_coordinates[i+1]))
		{
		    cout << "ERROR: Hmm::heuristic_viterbi : x_coordinates[" << i << "] (" 
			 << x_coordinates[i]
			 << ") >= x_coordinates[" << i+1 << "] (" << x_coordinates[i+1] 
			 << ") or y_coordinates[" << i << "] (" << y_coordinates[i]
			 << ") >= y_coordinates[" << i+1 << "] (" << y_coordinates[i+1] << ").\n" << flush;
		    check+=1;
		}
		
		// check that each rectangle can be calculated using the Hirschberg_viterbi
		// note: the rectangle that will later to calculated using the Hirschberg_viterbi
		// are not exactly the same as those, but smaller
	      
		if (max_area < ((y_coordinates[i+1]-y_coordinates[i]+1) * 2 * (min_strip_width+1)))
		{
		    cout << "ERROR: Hmm:: heuristic_viterbi : max_area (" << max_area
			 << ") < (y_coordinates[" << i+1 << "]-y_coordinates[" << i << "]+1) (" 
			 << (y_coordinates[i+1]-y_coordinates[i]+1) << ") * 2 * (min_strip_width+1) ("
			 << min_strip_width+1 << ") = " 
			 << (y_coordinates[i+1]-y_coordinates[i]+1) * 2 * (min_strip_width+1) 
			 << ". This rectangle cannot be calculated using the Hirschberg_viterbi."
			 << " Increase value of max_area to Start the calculation.\n" << flush;      
		    check+=1;
		}
		
	    }	  
	}
    }
    
    // check sequences

    if ((x->length()<1) || (x==NULL))
    {
	cout << "ERROR: Hmm::heuristic_viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL \n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL))
    {
	cout << "ERROR: Hmm::heuristic_viterbi : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL \n" << flush;
	check+=1;
    }
    
    if (check==0)
    {
	const int number_of_rectangles=number_of_xy_pairs+1;      

	// create new list of (x,y) coordinates including (0,0) and (x->length(),y->length())
	
	int* x_values= new int[number_of_rectangles+1];
	int* y_values= new int[number_of_rectangles+1];
	
	x_values[0]=0;
	x_values[number_of_rectangles]=x->length();
	y_values[0]=0;
	y_values[number_of_rectangles]=y->length();
	
	{
	    for (int i=1; i<number_of_rectangles; i++)
	    {
		x_values[i]=x_coordinates[i-1];
		y_values[i]=y_coordinates[i-1];
	    }
	}

      // determine offsets in x and y direction

	int x_offset=static_cast<int>(floor(static_cast<float>(x_margin)/2.));
	int y_offset=static_cast<int>(floor(static_cast<float>(y_margin)/2.));

	// delete old solutions

	// note: in order to use the function add_local_solution, 
	//       steps has to be initialised with -1 as well as 
	//
	//	    sequence_of_states[i]=-1;   
	//	    sequence_of_scores[i]=Logzero;    
	//	    sequence_of_xsteps[i]=-1;   
	//	    sequence_of_ysteps[i]=-1;    
	
	steps=-1;
	score=0;
	
	if (sequence_of_scores) delete [] sequence_of_scores;
	sequence_of_scores=NULL;
	
	if (sequence_of_states) delete [] sequence_of_states;
	sequence_of_states=NULL;
	
	if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	sequence_of_xsteps=NULL;
	
	if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	sequence_of_ysteps=NULL;
	
	condensed_steps=0;
	if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	condensed_sequence_of_scores=NULL;
	
	if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	condensed_sequence_of_states=NULL;
	
	if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	condensed_sequence_of_xsteps=NULL;
	
	if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	condensed_sequence_of_ysteps=NULL;      
	
	sequence_of_scores=new Score[x->length()+y->length()];
	sequence_of_states=new int[x->length()+y->length()];
	sequence_of_xsteps=new int[x->length()+y->length()];
	sequence_of_ysteps=new int[x->length()+y->length()];

	{
	    for (int i=0; i< (x->length()+y->length()); i++)               
	    {
		sequence_of_states[i]=-1;   
		sequence_of_scores[i]=Logzero;    
		sequence_of_xsteps[i]=-1;   
		sequence_of_ysteps[i]=-1;    
	    }
	}

	int rectangle=0;

	int x_start=0;
	int x_end=0;
	int y_start=0;
	int y_end=0;
	
	int next_x_start=0;
	int next_x_end=0;
	int next_y_start=0;
	int next_y_end=0;
	
	Score**** set_of_strips= new Score***[number_of_rectangles-1];
	
	Score*** strip=NULL;
	Score*** next_strip = NULL;
	
	x_start=x_values[rectangle];
	y_start=y_values[rectangle];	      

	if (number_of_rectangles==1) // if this is the first and the last rectangle
	{
	    x_end=x_values[rectangle+1];
	    y_end=y_values[rectangle+1];

	}
	else
	{
	    x_end=x_values[rectangle+1]+x_offset;
	    y_end=y_values[rectangle+1]+y_offset;

	}
	
	check+=allocate_memory_for_strip(this,
					 x_margin,
					 y_start, y_end,
					 &strip);
	
	if ((check==0) && (number_of_rectangles>1)) 
	    //
	// if there is more than 1 rectangle and if checks so far were o.k.
	{
	    for (rectangle=0; rectangle<number_of_rectangles-1; rectangle++)
		//
		// loop over all but the last rectangle
	    {
		if (check==0)
		{

		    if (rectangle==0)
		    {
			check+=get_strip(&strip, 
					 this,
					 NULL, // no unmirrored Hmm needed
					 x, x_start, x_end,
					 y, y_start, y_end, 
					 x_margin,
					 0,
					 0,
					 1);
		    }
		    else
		    {
			check+=get_strip(&strip, 
					 this,
					 NULL, // no unmirrored Hmm needed
					 x, x_start, x_end,
					 y, y_start, y_end, 
					 x_margin,
					 x_margin,
					 y_margin,
					 1);
		    }

		    check+=check_that_there_are_valid_values_in_strip(&strip,
								      this,
								      x_margin,
								      y_start, y_end);
		    
		    if (check==0) // if strip could be calculated and contains at least one valid element 
		    {
			// check results

			set_of_strips[rectangle]=strip;

		      
			if (rectangle < number_of_rectangles-2) // if next rectangle is not the last rectangle (i.e. if there remains a rectangle
			    // to be calculated by strip-method
			{
			    // get coordinates for next rectangle

			    next_x_start=0;
			    next_x_end=0;
			    next_y_start=0;
			    next_y_end=0;
			  
			    if (rectangle+1==0)  // if next rectangle is the first rectangle (never happens)
			    {
			    }
			    else if ((rectangle+1>0) && (rectangle+1<number_of_rectangles-1)) 
				// if this is an intermeditate rectangle
			    {
				next_x_start=x_values[rectangle+1]-x_offset;
				next_x_end=x_values[rectangle+2]+x_offset;
				next_y_start=y_values[rectangle+1]-y_offset;
				next_y_end=y_values[rectangle+2]+y_offset;	      

			    }
			    else // if this is the final rectangle (and not the first one)
			    {
				next_x_start=x_values[rectangle+1]-x_offset;
				next_x_end=x_values[rectangle+2];
				next_y_start=y_values[rectangle+1]-y_offset;
				next_y_end=y_values[rectangle+2];

			    }
			    
			    // allocate memory for next_strip

	  
			    check+=allocate_memory_for_strip(this,
							     x_margin,
							     next_y_start, next_y_end,
							     &next_strip);
			    
			    if (check==0) // if next_strip could be allocated
			    {

				check+=copy_rectangle_from_strip_to_next_strip(&strip,
									       // coordinates of strip
									       x_end-x_margin+1, x_end,
									       y_start, y_end,
									       // coordinates that shall be copied
									       x_end-x_margin+1, x_end,
									       y_end-y_margin+1, y_end,
									       &next_strip,
									       // coordinates next_strip (start)
									       x_end-x_margin+1, x_end,
									       next_y_start, next_y_end,
									       // coordinates for copied rectangle
									       x_end-x_margin+1, x_end,
									       y_end-y_margin+1, y_end,
									       // direction of the two strips
									       1);

				if (check==0) // if copying was successful
				{
				    // next_strip becomes strip
	  
				    strip      = next_strip;
				    next_strip = NULL;
				    
				    x_start=next_x_start;
				    x_end=next_x_end;
				    y_start=next_y_start;
				    y_end=next_y_end;
				}
				else
				{cout << "ERROR: Hmm::heuristic_viterbi : problems copying strip of rectangle " << rectangle 
				      << " into next_strip for next rectangle.\n" << flush;}
			    }
			    else
			    {
				cout << "ERROR: Hmm::heuristic_viterbi : problems allocating memory for next_strip for rectangle "
				     << rectangle+1 << "\n" << flush;
			    }
			} // if the next rectangle is not the last rectangle
		    }
		    else
		    {
			cout << "ERROR: Hmm::heuristic_viterbi : problems calculating rectangle " << rectangle 
			     << " with strip method. Either the calculation was not successfull or the strip contains only Logzero elements.\n" << flush;
		    }
		}
		else
		{
		    cout << "ERROR: Hmm::heuristic_viterbi : problems calculating rectangle " << rectangle-1 << " using the strip method.\n" << flush;
		}
	    } // loop over rectangles
	}
	else
	{ 
	    if (check!=0) 
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : problems allocating memory for strip for rectangle " << rectangle << "\n" << flush;
	    }
	}
      
	if (check==0)
	{
	    // calculate last remaining rectangle using calculate_viterbi_rectangle function
	    // which enables direct traceback using retrieve_state_path function

	    int next_end_state=0;
	    int next_x_margin=0;
	    int next_y_margin=0;
	    
	    int* start_states=NULL;
	    int number_of_start_states=0;
	 
	    if (number_of_rectangles>1)
	    {
		next_x_start=x_values[number_of_rectangles-1]-x_offset;
		next_x_end=x_values[number_of_rectangles];
		next_y_start=y_values[number_of_rectangles-1]-y_offset;
		next_y_end=y_values[number_of_rectangles];
		
		next_end_state=number_of_states-1;
		next_x_margin=x_margin;
		next_y_margin=y_margin;
	    }
	    else if (number_of_rectangles==1)
	    {
		next_x_start=x_values[number_of_rectangles-1];
		next_x_end=x_values[number_of_rectangles];
		next_y_start=y_values[number_of_rectangles-1];
		next_y_end=y_values[number_of_rectangles];
		
		next_end_state=number_of_states-1;
		next_x_margin=0;
		next_y_margin=0;
		
		start_states= new int[1];
		start_states[0]=0;
		number_of_start_states=1;
	    }

	    
	    Score*** viterbi_rectangle_array=NULL;
	    
	    check+=allocate_memory_for_viterbi_rectangle(this, 
							 next_x_start, next_x_end,
							 next_y_start, next_y_end,
							 1,
							 &viterbi_rectangle_array);

	    if (check!=0)
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : problems with function "
		     << "allocate_memory_for_viterbi_rectangle for last rectangle.\n" << flush;}
	    
	    if (check==0)
	    {
		if (number_of_rectangles>1)
		{

		    check+=copy_rectangle_from_strip_to_next_strip(&strip,
								   // coordinates of strip
								   x_end-x_margin+1, x_end,
								   y_start, y_end,
								   // coordinates that shall be copied
								   x_end-x_margin+1, x_end,
								   y_end-y_margin+1, y_end,
								   &viterbi_rectangle_array,
								   // coordinates next_strip
								   next_x_start, next_x_end,
								   next_y_start, next_y_end,
								   // coordinates for copied rectangle
								   next_x_start, next_x_start+x_margin-1,
								   next_y_start, next_y_start+y_margin-1,
								   // direction of both strips
								   1);
		}
		
		if (check!=0)
		{
		    cout << "ERROR: Hmm::heuristic_viterbi : problems with function "
			 << "copy_rectangle_from_strip_to_next_strip to copy last strip into last rectangle.\n" << flush;
		}
	    }
	    
	    if (check==0) // if rectangle = number_of_rectangles-1 o.k.
	    {
		for (rectangle=number_of_rectangles-1; rectangle>-1; rectangle--)
		{	  
		    // ----------------------------------------------------------------------
		    if (check==0) // if all previous rectangles o.k.
		    {

			// make sure that this local solution has no overlap with previous ones
			// remove end of state path of local solution if this rectangle in not the
			// last or the only one
			
			int offset=0;
			
			if ((rectangle < number_of_rectangles-1) && (number_of_rectangles>1))
			{

			    offset=-1;
			}
		      
			// initialise output values
			
			int x_start_traceback     = 0;
			int y_start_traceback     = 0;
			int start_state_traceback = 0;
			
			// calculate rectangle

			check += viterbi_rectangle(viterbi_rectangle_array,
						   this,
						   x, next_x_start, next_x_end, 
						   y, next_y_start, next_y_end, 
						   start_states, next_end_state,
						   number_of_start_states,
						   next_x_margin, next_y_margin,
						   offset, // offset for add local to global solution
						   // output values
						   &x_start_traceback, 
						   &y_start_traceback, 
						   &start_state_traceback);

			if (check!=0)
			{
			    cout << "ERROR: Hmm::heuristic_viterbi : problems with function "
				 << "viterbi_rectangle for rectangle " << rectangle << "\n" << flush;
			}
			
			if ((check==0) && (rectangle > 0)) // if this is not the last rectangle
			{
			    check+=delete_memory_for_viterbi_rectangle(this, 
								       next_x_start, next_x_end,
								       1,
								       &viterbi_rectangle_array);
			    if (check!=0)
			    {
				cout << "ERROR: Hmm::heuristic_viterbi : problems with function "
				     << "delete_memory_for_viterbi_rectangle for rectangle " << rectangle << "\n" << flush;
			    }
			}
			
			if (check==0) 
			{
			    if (rectangle-1 > -1) // if there is a rectangle left to be calculated
			    {

				if (rectangle-1 > 0) // if next rectangle is not the last one
				{
				    next_x_start   = x_values[rectangle-1]-x_offset;
				    next_x_end     = x_start_traceback;
				    next_y_start   = y_values[rectangle-1]-y_offset;
				    next_y_end     = y_start_traceback;
				    next_end_state = start_state_traceback;

				    viterbi_rectangle_array=NULL;
			      
				    check+=allocate_memory_for_viterbi_rectangle(this, 
										 next_x_start, next_x_end,
										 next_y_start, next_y_end,
										 1,
										 &viterbi_rectangle_array);
				    if (check!=0)
				    {
					cout << "ERROR: Hmm::heuristic_viterbi : problems with function "
					     << "allocate_memory_for_viterbi_rectangle for rectangle " 
					     << rectangle << "\n" << flush;}
				    
				    if (check==0)
				    {
					if (rectangle-1 > 1) // if there are more than 2 rectangles left
					{

					    check+=copy_rectangle_from_strip_to_next_strip(&set_of_strips[rectangle-2],
											 // coordinates of strip
											   x_values[rectangle-1]+x_offset-x_margin+1, 
											   x_values[rectangle-1]+x_offset,
											   y_values[rectangle-2]-y_offset,
											   y_values[rectangle-1]+y_offset,
											   // coordinates that shall be copied
											   x_values[rectangle-1]+x_offset-x_margin+1, 
											   x_values[rectangle-1]+x_offset,
											   y_values[rectangle-1]+y_offset-y_margin+1,
											   y_values[rectangle-1]+y_offset,
											   &viterbi_rectangle_array,
											   // coordinates next_strip (convert_labels_to_int("state_type","Start",MP->get_State_Type(),MP->get_Number_of_State_Types()))
											   next_x_start, next_x_end,
											   next_y_start, next_y_end,
											   // coordinates for copied rectangle
											   next_x_start, next_x_start+x_margin-1,
											   next_y_start, next_y_start+y_margin-1,
											   // direction of both strips
											   1);


					}			 
					else if (rectangle-1 == 1) // if there are exactly 2 rectangles left
					{

				    check+=copy_rectangle_from_strip_to_next_strip(&set_of_strips[rectangle-2],
											   // coordinates of strip
											   x_values[rectangle-1]+x_offset-x_margin+1, 
											   x_values[rectangle-1]+x_offset,
											   y_values[rectangle-2],
											   y_values[rectangle-1]+y_offset,
											   // coordinates that shall be copied
											   x_values[rectangle-1]+x_offset-x_margin+1, 
											   x_values[rectangle-1]+x_offset,
											   y_values[rectangle-1]+y_offset-y_margin+1,
											   y_values[rectangle-1]+y_offset,
											   &viterbi_rectangle_array,
											   // coordinates next Start
											   next_x_start, next_x_end,
											   next_y_start, next_y_end,
											   // coordinates for copied rectangle
											   next_x_start, next_x_start+x_margin-1,
											   next_y_start, next_y_start+y_margin-1,
											   // direction of both strips
											   1);

					}			 
					if ((check!=0) && (number_of_rectangles>1))
					{
					    cout << "ERROR: Hmm::heuristic_viterbi : problems copying "
						 << "strip of rectangle " << rectangle-2 << " into viterbi_rectangle for "
						 << "rectangle " << rectangle-1 << "\n" << flush;
					}
				    }
				}
				else if (rectangle-1 == 0) // if next rectangle is the last one
				{
				    next_x_start   = x_values[rectangle-1];
				    next_x_end     = x_start_traceback;
				    next_y_start   = y_values[rectangle-1];
				    next_y_end     = y_start_traceback;
				    next_end_state = start_state_traceback;
				    
				    viterbi_rectangle_array=NULL;
				    
				    // do not allocate memory for viterbi_rectangle as this will be done
				    // inside function viterbi_rectangle

				    start_states= new int[1];
				    start_states[0]=0;
				    number_of_start_states=1;
				    
				    next_x_margin=0;
				    next_y_margin=0;
				}

			    }
			}
		    } // if all previous rectangles o.k.
		} // loop over rectangles 
	    }
	    else
	    {
		cout << "ERROR: Hmm::heuristic_viterbi : problems with rectangle " << rectangle << "\n" << flush;}
	    
	    // delete memory 
	    
	    if (start_states) delete [] start_states;
	    start_states = NULL;
	}
	
	{
	    for (int i=0; i<number_of_rectangles-1; i++)
	    {	

		delete_memory_for_strip(this,					     
					x_margin,
					&set_of_strips[i]);
	    }
	}
	
	if (set_of_strips) delete [] set_of_strips;
	set_of_strips = NULL;
	
	if (x_values) delete [] x_values;
	x_values=NULL;
	
	if (y_values) delete [] y_values;
	y_values=NULL;
	
	
	if (check==0)
	{

	    check+=check_consistency_of_solution(x,y);
	}
      
	if (check != 0)
	{
	    cout << "ERROR: Hmm::heuristic_viterbi : solution did not pass consistency check, delete it.\n" << flush;
	    
	    steps=0;
	    score=0;
	    
	    if (sequence_of_scores) delete [] sequence_of_scores;
	    sequence_of_scores=NULL;
	    
	    if (sequence_of_states) delete [] sequence_of_states;
	    sequence_of_states=NULL;
	  
	    if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	    sequence_of_xsteps=NULL;
	  
	    if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	    sequence_of_ysteps=NULL;
	  
	    condensed_steps=0;
	    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	    condensed_sequence_of_scores=NULL;
	    
	    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	    condensed_sequence_of_states=NULL;
	  
	    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	    condensed_sequence_of_xsteps=NULL;
	  
	    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	    condensed_sequence_of_ysteps=NULL;      	  
	}

    } // if heuristic_viterbi can be started
    return(check);
}


int Hmm::new_heuristic_viterbi(const Hmm* mirror, // needed for internal calls to Hirschberg
				   const Sequence *x, const int* x_coordinates,
				   const Sequence *y, const int* y_coordinates,
				   const int number_of_xy_pairs,
				   const int x_margin,
				   const int y_margin,
				   const int max_area)
{
    // note : - if, during the calculation, it turns out that there
    //          is no single state path, the calculation is aborted and the return value check is != 0
    //        - if number_of_xy_pairs == 0 the Hirschberg_viterbi is called on the only (big) rectangle
    //
    //     check that at least one of the values in the copied rectangle is not Logzero
    //     otherwise stop and declare that no solution could be found
    
    int check=0;
    int min_strip_width=0;
    int i, j = 0;
    
    for (i=1; i<number_of_states-1; i++) { // loop over all states except Start and End state
	if (max(model[i].get_letters_to_read_x(), model[i].get_letters_to_read_y())>min_strip_width) 
	{
	    min_strip_width=max(model[i].get_letters_to_read_x(), model[i].get_letters_to_read_y());

	}
    }

    // check this Hmm
    
    if (mirrored != 0) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored (" << mirrored 
	     << ") of Hmm must be 0.\n" << flush;
	check+=1;
    }  
    if (number_of_states<3) 
    {
	cout << "ERROR: Hmm::new_heuristic_viterbi : number of states = " << number_of_states 
	     << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0) 
    {
	cout << "ERROR: Hmm::new_heuristic_viterbi : state 0 != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0) 
    {
	cout << "ERROR: Hmm::new_heuristic_viterbi : last state != End \n" << flush;
	cout <<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (i=0; i<number_of_states; i++) {
	if (model[i].get_alphabet()!=alphabet) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : alphabet of state i = " << i << ", alphabet = " 
		 << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (model[i].get_mirrored()!=mirrored) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored of state i = " << i << ", mirrored = " 
		 << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << mirrored << "\n" << flush;
	    check+=1;
	}
	if (model[i].get_number_of_states() != number_of_states) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : number of states in state i = " << i 
		 << ", number of states = " 
		 << model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}

	if ((i!=0) && (i!=(number_of_states-1))&&(model[i].get_letters_to_read() == 0)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : state i = " << i << " should emit! \n" << flush;
	    check+=1;
	}

    }
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : mirrored (" 
	     << mirror->get_mirrored() << ") must be 1.\n" << flush;
	check+=1;
    }
    if (mirror->get_number_of_states()<3) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : state 0 of type != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : last state != End \n" << flush;
	cout <<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	     <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (i=0; i<number_of_states; i++) {
	if (mirror->model[i].get_alphabet()!=alphabet) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : alphabet of state i = " 
		 << i << ", alphabet = " 
		 << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (mirror->model[i].get_mirrored()!=mirror->get_mirrored()) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : mirrored of state i = " 
		 << i << ", mirrored = " 
		 << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << mirror->get_mirrored() << "\n" << flush;
	    check+=1;
	}
	if (mirror->model[i].get_number_of_states() != number_of_states) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : number of states in state i = " 
		 << i << ", number of states = " 
		 << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}
	
	if ((i!=0) && (i!=(number_of_states-1))&&(mirror->model[i].get_letters_to_read() ==0)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored Hmm : state i = " 
		 << i << " should emit! \n" << flush;
	    check+=1;
	}
	
    }
    // check that this Hmm and the mirror Hmm belong together
    
    if (number_of_states != mirror->get_number_of_states()) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    for (i=1; i<number_of_states-1; i++) {
	if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read()) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : letters_to_read of model[" << i << "] (" 
		 << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		 << number_of_states-1-i << "] (" << mirror->model[number_of_states-1-i].get_letters_to_read() 
		 << ")\n" << flush;
	    check+=1;
	}
    }  
    // check that values of x_margin and y_margin exceed some minimal value
    
    if (x_margin < (min_strip_width+1)) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : x_margin (" << x_margin 
	     << ") < (min_strip_width+1) (" << min_strip_width+1 << ")\n" << flush;
	check+=1;
    }
    if (y_margin < (min_strip_width+1)) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : y_margin (" << y_margin 
	     << ") < (min_strip_width+1) (" << min_strip_width+1 << ")\n" << flush;
	check+=1;
    }
    // check that x_margin and y_margin are larger than some desired value 
    // note : the desired x_margin and y_margin values are large enough to accomodate
    //        at least three of the states with maximum number of reads per sequence (this number
    //        is equal to min_strip_width) next to eachother 

    if ( abs(floor(static_cast<float>(x_margin)/2.) - ceil(static_cast<float>(x_margin)/2.)) == 0) {
	// if x_margin is even
	if (x_margin < (3*min_strip_width+2)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : x_margin (" << x_margin 
		 << ") is even. x_margin < (3 * min_strip_width (" << min_strip_width << ") + 1) (" 
		 << 3*min_strip_width+1 
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    else { 
	// if x_margin is odd
	if (x_margin < (3*min_strip_width)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : x_margin (" << x_margin 
		 << ") is odd. x_margin < (3 * min_strip_width (" << min_strip_width << ")) (" << 3*min_strip_width
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    if ( abs(floor(static_cast<float>(y_margin)/2.) - ceil(static_cast<float>(y_margin)/2.)) == 0) {
	// if y_margin is even
	if (y_margin < (3*min_strip_width+2)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : y_margin (" << y_margin 
		 << ") is even. y_margin < (3 * min_strip_width (" << min_strip_width << ") + 1) (" 
		 << 3*min_strip_width+1 
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    else {
	// if y_margin is odd
	if (y_margin < (3*min_strip_width)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : y_margin (" << y_margin 
		 << ") is odd. y_margin < (3 * min_strip_width (" << min_strip_width << ")) (" << 3*min_strip_width
		 << "). Choose a larger value.\n" << flush;
	    check+=1;
	}
    }
    
    // check Hmm
  
    if (mirrored != 0) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored (" << mirrored << ") must be 0.\n" << flush;
	check+=1;
    }
    if (number_of_states<3) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : number of states = " << number_of_states 
	     << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : state 0 != Start \n" << flush;
	cout <<"model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : last state != End \n" << flush;
	cout <<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	     << model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    for (i=0; i<number_of_states; i++) {
	if (model[i].get_alphabet()!=alphabet) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : alphabet of state i = " << i << ", alphabet = " 
		 << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		 << alphabet << "\n" << flush;
	    check+=1;
	}
	if (model[i].get_mirrored()!=mirrored) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : mirrored of state i = " << i << ", mirrored = " 
		 << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		 << mirrored << "\n" << flush;
	    check+=1;
	}
	if (model[i].get_number_of_states() != number_of_states) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : number of states in state i = " << i 
		 << ", number of states = " 
		 << model[i].get_number_of_states() << " != number of states of Hmm = " 
		 << number_of_states << "\n" << flush;
	    check+=1;
	}

	if ((i!=0) && (i!=(number_of_states-1)) &&(model[i].get_letters_to_read() == 0)) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : state i = " << i << " should emit! \n" << flush;
	    check+=1;
	}
	
    }
    // check that max_area value is sensible
    
    if (max_area<1) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : max_area (" << max_area
	     << ") < 1.\n" << flush;
	check+=1;
    }
    // check the number of (x,y) constraints, their order and the size of rectangles they span
    
    if (number_of_xy_pairs>0) {
	for (i=0; i<number_of_xy_pairs; i++) {
	    // check x_margin and y_margin and their compatibility with x_start, x_end
	    // y_start and y_end
	    if (i==0) { 
		// check rectangle from (0,0) to first (x,y) pair
		if (((x_coordinates[i]-0+1) < x_margin+1) ||  ((y_coordinates[i]-0+1) < y_margin+1)) { 
		    cout << "ERROR: Hmm::new_heuristic_viterbi : area of rectangle between (0,0) and (x_coordinates[" 
			 << i << "] (" << x_coordinates[i] << "), y_coordinates[" << i << "] (" << y_coordinates[i] 
			 << ")) cannot accumodate area of (x_margin (" << x_margin << ") + 1) * (y_margin (" 
			 << y_margin << ") + 1).\n" << flush;
		    check+=1;
		}
	    }
	    else if ((i>0) && (i<number_of_xy_pairs-1)) {
		// check intermediate rectangle from (x,y) pair i to (x,y) pair i+1
		if ( ((x_coordinates[i+1]-x_coordinates[i]+1) < x_margin+1) || 
		     ((y_coordinates[i+1]-y_coordinates[i]+1) < y_margin+1)) {
		    cout << "ERROR: Hmm::new_heuristic_viterbi : area of rectangle between (x_coordinates[" 
			 << i << "] (" << x_coordinates[i] << "), y_coordinates[" << i << "] (" << y_coordinates[i] 
			 << ")) and (x_coordinates[" << i+1
			 << "] (" << x_coordinates[i+1] << "), y_coordinates[" << i+1 << "] (" << y_coordinates[i+1] 
			 << ")) cannot accumodate area of (x_margin (" << x_margin << ") + 1) * (y_margin (" 
			 << y_margin << ") + 1).\n" << flush;
		    check+=1;
		}
	    }
	    else { 
		// check rectangle from last (x,y) pair to (x->length(), y->length())
		if ( ((x->length()-x_coordinates[i]+1) < x_margin+1) ||  
		     ((y->length()-y_coordinates[i]+1) < y_margin+1)) 
		{
		    cout << "ERROR: Hmm::new_heuristic_viterbi : area of rectangle between (x_coordinates[" 
			 << i << "] (" << x_coordinates[i] << "), y_coordinates[" << i << "] (" 
			 << y_coordinates[i] << ")) and (x->length() (" << x->length()
			 << "), y->length() (" << y->length()
			 << ")) cannot accumodate area of (x_margin (" << x_margin << ") + 1) * (y_margin (" 
			 << y_margin << ") + 1).\n" << flush;
		    check+=1;
		}
	    }
	    // check that no constaint is equal to 0 or x->length() or y->length(), respectively
	    
	    if (x_coordinates[i]==0 || y_coordinates[i]==0 ||
		x_coordinates[i]==x->length() || y_coordinates[i]==y->length()) {
		cout << "ERROR: Hmm::new_heuristic_viterbi : x_coordinates[" << i << "] (" << x_coordinates[i]
		     << ") = 0 or y_coordinates[" << i << "] (" << y_coordinates[i] 
		     << ") = 0 or x_coordinates[" << i << "] (" << x_coordinates[i] << ") == x->length() ("
		     << x->length() << ") or y_coordinates[" << i << "] (" << y_coordinates[i] 
		     << ") == y->length() (" << y->length() << ")\n" << flush;
		check+=1;
	    }
	    // check that 0 < x_i < x->length() and 0 < y_i < y->length()
	    
	    if ((x_coordinates[i] <= 0  ||  x_coordinates[i] >= x->length()) ||
		(y_coordinates[i] <= 0  ||  y_coordinates[i] >= y->length())) {
		cout << "ERROR: Hmm::new_heuristic_viterbi : x_coordinates[" << i << "] (" << x_coordinates[i]
		     << ") out of range [1, " << x->length()-1 
		     << "] or y_coordinates[" << i << "] (" << y_coordinates[i]
		     << ") out of range [1, " << y->length()-1 << "].\n" << flush;
		check+=1;
	    }
	    if (i<(number_of_xy_pairs-1)) {
		// check that x_i < x_i+1 and y_i < y_i+1
		if ((x_coordinates[i]>=x_coordinates[i+1]) ||
		    (y_coordinates[i]>=y_coordinates[i+1])) {
		    cout << "ERROR: Hmm::new_heuristic_viterbi : x_coordinates[" << i << "] (" 
			 << x_coordinates[i]
			 << ") >= x_coordinates[" << i+1 << "] (" << x_coordinates[i+1] 
			 << ") or y_coordinates[" << i << "] (" << y_coordinates[i]
			 << ") >= y_coordinates[" << i+1 << "] (" << y_coordinates[i+1] << ").\n" << flush;
		    check+=1;
		}
		// check that each rectangle can be calculated using the Hirschberg_viterbi
		// note: the rectangle that will later to calculated using the Hirschberg_viterbi
		// are not exactly the same as those, but smaller
		
		if (max_area < ((y_coordinates[i+1]-y_coordinates[i]+1) * 2 * (min_strip_width+1))) {
		    cout << "ERROR: Hmm:: new_heuristic_viterbi : max_area (" << max_area
			 << ") < (y_coordinates[" << i+1 << "]-y_coordinates[" << i << "]+1) (" 
			 << (y_coordinates[i+1]-y_coordinates[i]+1) << ") * 2 * (min_strip_width+1) ("
			 << min_strip_width+1 << ") = " 
			 << (y_coordinates[i+1]-y_coordinates[i]+1) * 2 * (min_strip_width+1) 
			 << ". This rectangle cannot be calculated using the Hirschberg_viterbi."
			 << " Increase value of max_area to Start the calculation.\n" << flush;      
		    check+=1;
		}
		

	    }	  
	}
    }
    // check sequences
    
    if ((x->length()<1) || (x==NULL)) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : length of sequence x (" << x->length()
	     << ") < 1 or x (" << x << ") NULL \n" << flush;
	check+=1;
    }
    if ((y->length()<1) || (y==NULL)) {
	cout << "ERROR: Hmm::new_heuristic_viterbi : length of sequence y (" << y->length()
	     << ") < 1 or y (" << y << ") NULL \n" << flush;
	check+=1;
    }
    
    if (check==0) {
	
//#ifdef _PRINT
	cout << "Hmm::New_Heuristic_Viterbi\n" << flush;
	cout << "----------------------------------------------------------------------\n" << flush;
	
	cout << "x_margin           = " << x_margin << "\n" << flush;
	cout << "y_margin           = " << y_margin << "\n" << flush;
	cout << "min_strip_width    = " << min_strip_width << "\n" << flush;
	cout << "number_of_xy_pairs = " << number_of_xy_pairs << "\n" << flush;
	for (i=0; i<number_of_xy_pairs; i++) {
	    cout << "x_coordinates[" << i << "] = " << x_coordinates[i] << "\n" << flush;
	    cout << "y_coordinates[" << i << "] = " << y_coordinates[i] << "\n" << flush;
	}
	cout << "max_area           = " << max_area << "\n" << flush;
//#endif  
	
	const int number_of_rectangles=number_of_xy_pairs+1;      

	if (number_of_rectangles == 1) {
	    
	    // use Hirschberg_viterbi if there is just one (big) rectangle
	    
	    check += Hirschberg_viterbi(mirror,			  
					x, y,
					max_area);
	    
	    if (check != 0) {
		cout << "ERROR: Hmm::new_heuristic_viterbi : error occurred in function Hmm:"
		     << "Hirschberg_viterbi when dealing with the only (big) rectangle.\n";
	    }
	}
	else { // if there is more than one rectangle
	    
	    // create new list of (x,y) coordinates including (0,0) and (x->length(),y->length())
	    
	    int* x_values= new int[number_of_rectangles+1];
	    int* y_values= new int[number_of_rectangles+1];
	    
	    x_values[0]=0;
	    x_values[number_of_rectangles]=x->length();
	    y_values[0]=0;
	    y_values[number_of_rectangles]=y->length();
    
	    for (i=1; i<number_of_rectangles; i++) {
		x_values[i]=x_coordinates[i-1];
		y_values[i]=y_coordinates[i-1];
	    }

	    // determine offsets in x and y direction

	    int x_offset=static_cast<int>(floor(static_cast<float>(x_margin)/2.));
	    int y_offset=static_cast<int>(floor(static_cast<float>(y_margin)/2.));

	    // delete old solutions

	    // note: in order to use the function add_local_solution, 
	    //       steps has to be initialised with -1 as well as 
	    //
	    //	    sequence_of_states[i]=-1;   
	    //	    sequence_of_scores[i]=Logzero;    
	    //	    sequence_of_xsteps[i]=-1;   
	    //	    sequence_of_ysteps[i]=-1;    
	    
	    steps=-1;
	    score=0;
	    
	    if (sequence_of_scores) delete [] sequence_of_scores;
	    sequence_of_scores=NULL;
	    
	    if (sequence_of_states) delete [] sequence_of_states;
	    sequence_of_states=NULL;
	    
	    if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	    sequence_of_xsteps=NULL;
	    
	    if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	    sequence_of_ysteps=NULL;
    
	    condensed_steps=0;
	    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	    condensed_sequence_of_scores=NULL;
      
	    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	    condensed_sequence_of_states=NULL;
	    
	    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	    condensed_sequence_of_xsteps=NULL;
	    
	    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	    condensed_sequence_of_ysteps=NULL;      

	    sequence_of_scores=new Score[x->length()+y->length()];
	    sequence_of_states=new int[x->length()+y->length()];
	    sequence_of_xsteps=new int[x->length()+y->length()];
	    sequence_of_ysteps=new int[x->length()+y->length()];

	    for (i=0; i< (x->length()+y->length()); i++) {
		sequence_of_states[i]=-1;   
		sequence_of_scores[i]=Logzero;    
		sequence_of_xsteps[i]=-1;   
		sequence_of_ysteps[i]=-1;    
	    }
	    
	    int rectangle=0;
	    
	    int x_start=0;
	    int x_end=0;
	    int y_start=0;
	    int y_end=0;

	    int next_x_start=0;
	    int next_x_end=0;
	    int next_y_start=0;
	    int next_y_end=0;
	    
	    int next_end_state = 0;
	    int next_x_margin  = 0;
	    int next_y_margin  = 0;

	    int* start_states = NULL;
	    int  number_of_start_states = 0;

	    Score**** set_of_strips  = new Score***[number_of_rectangles-1];
	    int*      x_start_values = new int[number_of_rectangles-1];
	    int*      x_end_values   = new int[number_of_rectangles-1];
	    int*      y_start_values = new int[number_of_rectangles-1];
	    int*      y_end_values   = new int[number_of_rectangles-1];
	    
	    Score*** strip=NULL;
	    Score*** next_strip=NULL;
	    
	    x_start=x_values[rectangle];
	    y_start=y_values[rectangle];	      

	    if (number_of_rectangles==1) { 
		// if this is the first and the last rectangle
		x_end=x_values[rectangle+1];
		y_end=y_values[rectangle+1];

	    }
	    else {
		
		x_end=x_values[rectangle+1]+x_offset;
		y_end=y_values[rectangle+1]+y_offset;

	    }

	    cout << "1 allocate_memory_for_strip\n" << flush;
	    cout << " x " << x_margin << " y_1 " << y_start << " y_2 " << y_end << " strip " << strip << "\n" << flush;
	    
	    check+=allocate_memory_for_strip(this,
					     x_margin,
					     y_start, y_end,
					     &strip);
	    
	    cout << "after 1 allocate_memory_for_strip\n" << flush;
	    
	    if (check != 0) {
		cout << "ERROR: Hmm::new_heuristic_viterbi : problems allocating memory for "
		     << "strip for rectangle 0.\n" << flush;
	    }
	    if ((check==0) && (number_of_rectangles>1)) {
		//
		// if there is more than 1 rectangle and if checks so far were o.k.
		for (rectangle=0; rectangle<number_of_rectangles-1; rectangle++) {
		    //
		    // loop over all but the last rectangle
		    if (check==0) {

			if (rectangle==0) {
			    
			    check+=get_strip(&strip, 
					     this,
					     NULL, // no unmirrored Hmm needed
					     x, x_start, x_end,
					     y, y_start, y_end, 
					     x_margin,
					     0,
					     0,
					     1);
			}
			else {
			    
			    check+=get_strip(&strip, 
					     this,
					     NULL, // no unmirrored Hmm needed
					     x, x_start, x_end,
					     y, y_start, y_end, 
					     x_margin,
					     x_margin,
					     y_margin,
					     1);
			}
			if (check != 0) {
			    cout << "ERROR: Hmm::new_heuristic_viterbi : problems in function "
				 << "Hmm::get_strip for rectangle " << rectangle << "\n" << flush;
			}
			
			if (check == 0) {
			    
			    check+=check_that_there_are_valid_values_in_strip(&strip,
									      this,
									      x_margin,
									      y_start, y_end);
			    if (check != 0) {
				cout << "ERROR: Hmm::new_heuristic_viterbi : problems in function "
				     << "check_that_there_are_valid_values_in_strip for rectangle " << rectangle << "\n" << flush;
			    }
			}	    
			
			if (check==0) {
			    
			    // if strip could be calculated and contains at least one valid element 			    
			    set_of_strips[rectangle]  = strip;
			    x_start_values[rectangle] = x_start;
			    x_end_values[rectangle]   = x_end;  
			    y_start_values[rectangle] = y_start;
			    y_end_values[rectangle]   = y_end;  

			    next_x_start=0;
			    next_x_end=0;
			    next_y_start=0;
			    next_y_end=0;
			    
			    if (rectangle+1==0) {} // if next rectangle is the first rectangle (never happens)
			    else if ((rectangle+1>0) && (rectangle+1<number_of_rectangles-1)) {
				// if this is an intermeditate rectangle
				
				next_x_start=x_values[rectangle+1]-x_offset;
				next_x_end=x_values[rectangle+2]+x_offset;
				next_y_start=y_values[rectangle+1]-y_offset;
				next_y_end=y_values[rectangle+2]+y_offset;	      

				// allocate memory for next_strip
//#ifdef _PRINT

				cout << "2 allocate_memory_for_strip\n" << flush;
				cout << " x " << x_margin << " y_1 " << next_y_start << " y_2 " << next_y_end << " strip " << next_strip 
				     << "\n" << flush;
				
//#endif
				check+=allocate_memory_for_strip(this,
								 x_margin,
								 next_y_start, next_y_end,
								 &next_strip);
				
				cout << "after 2 allocate_memory_for_strip\n" << flush;
				
				if (check != 0) {
				    cout << "ERROR: Hmm::new_heuristic_viterbi : problems allocating memory for "
					 << "next_strip for next rectangle " << rectangle+1 << "\n" << flush;
				}
				
				if (check == 0) {

				    check+=copy_rectangle_from_strip_to_next_strip(&strip,
										   // coordinates of strip
										   x_end-x_margin+1, x_end,
										   y_start, y_end,
										   // coordinates that shall be copied
										   x_end-x_margin+1, x_end,
										   y_end-y_margin+1, y_end,
										   &next_strip,
								                   // coordinates next_strip Start
										   x_end-x_margin+1, x_end,
										   next_y_start, next_y_end,
										   // coordinates for copied rectangle
										   x_end-x_margin+1, x_end,
										   y_end-y_margin+1, y_end,
										   // direction of the two strips
										   1);

				    if (check != 0) {
					cout << "ERROR: Hmm::new_heuristic_viterbi : problems copying strip of rectangle " 
					     << rectangle << " into next_strip for next rectangle.\n" << flush;
				    }
				}

				// next_strip becomes strip
		
				strip=next_strip;
				next_strip = NULL;
				
				x_start=next_x_start;
				x_end=next_x_end;
				y_start=next_y_start;
				y_end=next_y_end;

			    }
			    else {
				// if this is the final rectangle (and not the first one)

				next_x_start=x_values[rectangle+1]-x_offset;
				next_x_end=x_values[rectangle+2];
				next_y_start=y_values[rectangle+1]-y_offset;
				next_y_end=y_values[rectangle+2];
				
				next_end_state = number_of_states-1;
				next_x_margin  = x_margin;
				next_y_margin  = y_margin;

			    }
			}
			else {
			    cout << "ERROR: Hmm::new_heuristic_viterbi : problems calculating rectangle " << rectangle 
				 << " with strip method. Either the calculation was not successfull or the strip "
				 << "contains only Logzero elements.\n" << flush;
			}
		    }
		    else {
			cout << "ERROR: Hmm::new_heuristic_viterbi : problems calculating rectangle " 
			     << rectangle-1 << " using the strip method.\n" << flush;
		    }
		} // loop over rectangles
	    }
	    else {
		if (check!=0) {
		    cout << "ERROR: Hmm::new_heuristic_viterbi : problems allocating memory for strip for rectangle " 
			 << rectangle << "\n" << flush;
		}
		else {
		    
		    // if there is just one rectangle
		    
		    next_x_start=x_values[number_of_rectangles-1];
		    next_x_end=x_values[number_of_rectangles];
		    next_y_start=y_values[number_of_rectangles-1];
		    next_y_end=y_values[number_of_rectangles];
		  
		    next_end_state = number_of_states-1;
		    next_x_margin  = 0;
		    next_y_margin  = 0;
		    
		    start_states = new int[1];
		    start_states[1] = 0;
		    number_of_start_states = 0;

		}
	    }
	    
	    if (check==0) {
		
		// calculate last rectangle using calculate_viterbi_rectangle function
		// which enables direct traceback using retrieve_state_path function
		
		Score*** old_viterbi_strip = NULL;
		Score*** new_viterbi_strip = NULL;
		
		int offset                = 0;
		int x_start_traceback     = 0;
		int y_start_traceback     = 0;
		int start_state_traceback = 0;
		
		int old_x_start = 0;
		int old_x_end   = 0;
		int old_y_start = 0;
		int old_y_end   = 0;
		
		for (rectangle=number_of_rectangles-1; rectangle>-1; rectangle--) {
		    // ----------------------------------------------------------------------
		    if (check==0) { 
			// if all previous rectangles o.k.

			old_viterbi_strip = NULL;
			new_viterbi_strip = NULL;
			
			if ((number_of_rectangles > 0) && (rectangle > 0)) { 

			    old_viterbi_strip = set_of_strips[rectangle-1];
			    
			    old_x_start = x_start_values[rectangle-1];
			    old_x_end   = x_end_values[rectangle-1];  
			    old_y_start = y_start_values[rectangle-1];
			    old_y_end   = y_end_values[rectangle-1];  

			    cout << "3 allocate_memory_for_strip\n" << flush;
			    cout << " x " << x_margin << " y_1 " << next_y_start << " y_2 " << next_y_end << " strip " 
				 << new_viterbi_strip << "\n" << flush;
			    
			    check += allocate_memory_for_strip(this,
							       x_margin,
							       next_y_start, next_y_end,
							       &new_viterbi_strip);

			    cout << "after 3 allocate_memory_for_strip\n" << flush;
			    
			    if ((check != 0) || (new_viterbi_strip == NULL)) {
				cout << "ERROR: Hmm::new_heuristic_viterbi: error occurred in function "
				     << "Hmm::allocate_memory_for_strip for rectangle " << rectangle
				     << " (x_start, x_end) = (" << next_x_start << ", " << next_x_end 
				     << ") (y_start, y_end) = (" << next_y_start << ", " << next_y_end 
				     << ").\n" << flush;
			    }
			    if (check == 0) {

				check += copy_rectangle_from_strip_to_next_strip(// source strip:
				    &old_viterbi_strip,
				    // coordinates of strip
				    old_x_end-x_margin+1, old_x_end,
				    old_y_start, old_y_end,
				    // coordinates to be copied
				    old_x_end-x_margin+1, old_x_end,
				    old_y_end-y_margin+1, old_y_end,
				    // target strip:
				    &new_viterbi_strip,
				    // coordinates of strip
				    next_x_start, next_x_start+x_margin-1,
				    next_y_start, next_y_end,
				    // coordinates to be copied to
				    next_x_start, next_x_start+x_margin-1,
				    next_y_start, next_y_start+y_margin-1,
				    // direction of two strips
				    1);

				if (check != 0) {
				    cout << "ERROR: Hmm::new_heuristic_viterbi: error occurred in function "
					 << "Hmm::copy_rectangle_from_strip_to_next_strip for rectangle " << rectangle
					 << " (x_start, x_end) = (" << next_x_start << ", " << next_x_end 
					 << ") (y_start, y_end) = (" << next_y_start << ", " << next_y_end 
					 << ").\n" << flush;
				}
			    } // if check == 0
			    
			    // delete memory

			    if (old_viterbi_strip != NULL) {

				check += delete_memory_for_strip(this,
								 x_margin,
								 &old_viterbi_strip);
				
				if ((check != 0) || (old_viterbi_strip != NULL)) {
				    cout << "ERROR: Hmm::new_heuristic_viterbi: error occurred in function "
					 << "Hmm::delete_memory_for_strip for old_viterbi_strip for rectangle " << rectangle
					 << " (x_start, x_end) = (" << next_x_start << ", " << next_x_end 
					 << ") (y_start, y_end) = (" << next_y_start << ", " << next_y_end 
					 << ").\n" << flush;
				}
			    } // if old_viterbi_strip != NULL
			} // if ((number_of_rectangles > 0) && (rectangle > 0)) { 

			// make sure that this local solution has no overlap with previous ones
			// remove end of state path of local solution if this rectangle in not the
			// last or the only one
			
			offset=0;
			
			if ((rectangle < number_of_rectangles-1) && (number_of_rectangles>1)) {

			    offset=-1;
			}

			// initialise output values
			
			x_start_traceback     = 0;
			y_start_traceback     = 0;
			start_state_traceback = 0;
			
			if (check == 0) {
			    
			    // calculate rectangle
			    
			    check+=Hirschberg_viterbi_rectangle(new_viterbi_strip,
								mirror,
								x, next_x_start, next_x_end, 
								y, next_y_start, next_y_end, 
								start_states, next_end_state,
								number_of_start_states,
								next_x_margin, next_y_margin,
								offset, // offset for add local to global solution
								max_area, 
								// output values
								&x_start_traceback, 
								&y_start_traceback, 
								&start_state_traceback);

			    if (check!=0) {
				cout << "ERROR: Hmm::new_heuristic_viterbi : problems with function "
				     << "Hirschberg_viterbi_rectangle for rectangle " << rectangle
				     << " (x_start, x_end) = (" << next_x_start << ", " << next_x_end 
				     << ") (y_start, y_end) = (" << next_y_start << ", " << next_y_end 
				     << ").\n" << flush;
			    }
			} // if check == 0
			
			// delete memory

			if (new_viterbi_strip != NULL) {

			    check += delete_memory_for_strip(this,
							     x_margin,
							     &new_viterbi_strip);
			    
			    if ((check != 0) || (new_viterbi_strip != NULL)) {
				cout << "ERROR: Hmm::new_heuristic_viterbi: error occurred in function "
				     << "Hmm::delete_memory_for_strip for new_viterbi_strip for rectangle " << rectangle
		     << " (x_start, x_end) = (" << next_x_start << ", " << next_x_end 
				     << ") (y_start, y_end) = (" << next_y_start << ", " << next_y_end 
				     << ").\n" << flush;
			    }
			}

			if (check == 0) {
			    if (rectangle-1 > -1) { 
				// if there is a rectangle left to be calculated

				if (rectangle-1 > 0) { 
				    // if next rectangle is not the last one
				    
				    next_x_start   = x_values[rectangle-1]-x_offset;
				    next_x_end     = x_start_traceback;
				    next_y_start   = y_values[rectangle-1]-y_offset;
				    next_y_end     = y_start_traceback;
				    next_end_state = start_state_traceback;
				}
				else if (rectangle-1 == 0) { 
				    // if next rectangle is the last one
				    
				    next_x_start   = x_values[rectangle-1];
				    next_x_end     = x_start_traceback;
				    next_y_start   = y_values[rectangle-1];
				    next_y_end     = y_start_traceback;
				    next_end_state = start_state_traceback;

				    if (start_states) delete [] start_states;
				    start_states = NULL;
				    
				    start_states = new int[1];
				    start_states[0]=0;
				    number_of_start_states=1;
		  
				    next_x_margin=0;
				    next_y_margin=0;
				}
			    } // if there is a rectangle left to be calculated
			} // if check == 0	    
		    } // if all previous rectangles ok
		    else { 
			cout << "ERROR: Hmm::new_heuristic_viterbi : problems with rectangle " 
			     << rectangle+1 << "\n" << flush;
		    }
		} // loop over rectangles 
	    } // if check == 0
	    
	    // delete memory 
	    
	    if (start_states) delete [] start_states;
	    start_states = NULL;
	    
	    if (set_of_strips) delete [] set_of_strips;
	    set_of_strips = NULL;
      
	    if (x_start_values) delete [] x_start_values;
	    x_start_values = NULL;

	    if (x_end_values) delete [] x_end_values;
	    x_end_values = NULL;

	    if (y_start_values) delete [] y_start_values;
	    y_start_values = NULL;

	    if (y_end_values) delete [] y_end_values;
	    y_end_values = NULL;

	    if (x_values) delete [] x_values;
	    x_values=NULL;
      
	    if (y_values) delete [] y_values;
	    y_values=NULL;

	} // if there is more than one rectangle

	if (check==0) {

	    check+=check_consistency_of_solution(x,y);
	}

	if (check != 0) {
	    cout << "ERROR: Hmm::new_heuristic_viterbi : solution did not pass consistency check, "
	     << "delete it.\n" << flush;
	
	    steps=0;
	    score=0;
	    
	    if (sequence_of_scores) delete [] sequence_of_scores;
	    sequence_of_scores=NULL;
	
	    if (sequence_of_states) delete [] sequence_of_states;
	    sequence_of_states=NULL;
	    
	    if (sequence_of_xsteps) delete [] sequence_of_xsteps;
	    sequence_of_xsteps=NULL;
	    
	    if (sequence_of_ysteps) delete [] sequence_of_ysteps;
	    sequence_of_ysteps=NULL;
	    
	    condensed_steps=0;
	    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
	    condensed_sequence_of_scores=NULL;
	    
	    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
	    condensed_sequence_of_states=NULL;
	    
	    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
	    condensed_sequence_of_xsteps=NULL;
	
	    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
	    condensed_sequence_of_ysteps=NULL;      	  
	}
    } // if new_heuristic_viterbi can be started
    return(check);
}

int Hmm::calculate_scores_from_probs()
{
    // note: - calculates the transition scores and scores
    //       - does not calculate emission scores and resets the emission scores table 
    //       - only Onlyrealscore implemented so far
    
    int check = 0;

    if (!check)
    {
	
	// first reset all scores with Logzeros
	
	for (int state=0;state<number_of_states;state++)
	{
	    model[state].reset_transition_scores();
	    model[state].reset_emission_scores();
	} 
	
	int d, max;
	int* dim=NULL;
	array<int> indices(1);
	array<int> large_indices(1);
	
	Score score_value;
	Prob real_prob;
	
	// calculate transition scores
	// ======================================================================
	
	for (int to=0;to<number_of_states;to++)
	{
	    for (int from=0;from<number_of_states;from++)
	    {      
		real_prob=0;
		score_value=Logzero;
		real_prob=model[from].get_transition_prob(to);
		
		if (real_prob!=0)
		{		       
		    //score_value=log(real_prob)/Logbase;
		    score_value = log(real_prob);
		}			
		if (model[to].is_state_previous_state(from) == 1)
		{
		    check+=model[to].set_transition_score(from, score_value);
		}
		if(check){
		    cout<<"ERROR: calculate_score_from_probs, set_transition_score"<<endl;
		    break;
		}
		//model[to].set_score_type(s);		
	    } // loop over from states      
	} // loop over to states
    
	// calculate emission scores
	// ======================================================================
	
	for (int to=0;to<number_of_states;to++)
	{
	    d=model[to].get_letters_to_read();
	    if (d!=0) 
	    {
		max= static_cast<int>(pow( static_cast<float>(alphabet), static_cast<float>(d)));
	    }
	    else
	    {
		max=0;
	    }
	  
	    if (d>0) 
	    {
		dim = new int[d];
		for (int i=0; i<d; i++)
		{
		    dim[i]=alphabet;
		}
		indices.SetDimension(0,d);
	    }
	    else
	    {
		dim=NULL;
	    }
	  
	    if (max>0) // if emitting state
	    {
		for (int i=0; i<max; i++)
		{		    
		    real_prob=0;
		    score_value=Logzero;
		  
		    real_prob=model[to].get_emission_prob(i);
		    
		    if (real_prob!=0)
		    {		       
			//score_value=log(real_prob)/Logbase;
			score_value = log(real_prob);
		    }
		    
		    indices=get_indices(i, dim, d,check);		    
		    
		    if(check){
			cout<<"ERROR: calculate_score_from_probs, get_indices, indices : ";
			indices.Print(cout);
			cout<<endl;
		    }
		    
		    model[to].set_emission_score(indices, score_value);
		    //model[to].set_score_type(s);
		    
		} // loop over int i (i.e. all emission_probs of to state)
	    } // if max > 0 i.e. if emitting state
	    
	    if (dim) delete [] dim;
	    dim = NULL;
	} // loop over to states
    } // if check == 0
    return(check);
}

Hmm & Hmm::operator = (const Hmm &p)
{
    if ( this != &p)
    {
	// note (see = operator for Hmm_State) forward and 
	// backward scores of each state will not be copied,
	// Sequences x and y will not be copied

	mirrored=p.mirrored;
	alphabet=p.alphabet;
	number_of_states=p.number_of_states;
	pair = p.pair;

	if (model) delete [] model;
	model = NULL;
	model = new Hmm_State[number_of_states];
	
	for (int j=0; j<number_of_states; j++)
	{	

	    model[j]=p.model[j]; 
	}
	
	// results of any alignment will not be copied
	
	steps=0;
	score=0;
	sequence_of_scores=NULL;
	sequence_of_states=NULL;
	sequence_of_xsteps=NULL;
	sequence_of_ysteps=NULL;

	condensed_steps=0;
	condensed_sequence_of_scores=NULL;
	condensed_sequence_of_states=NULL;    
	condensed_sequence_of_xsteps=NULL;
	condensed_sequence_of_ysteps=NULL;
    }
    return(*this);
}

Hmm_State* Hmm::operator [] (int i)
{
    if ((i<0) || (i>(number_of_states-1)))
    {cout << "ERROR: operator [] : i " << i << " out of range [0, " << number_of_states << "].\n" << flush;}

    return (&model [i]);
}

void Hmm::set_information_of_fake_alignment(const Score  c_score,
						const int    c_steps,
						const int*   c_sequence_of_states,
						const int*   c_sequence_of_xsteps,
						const int*   c_sequence_of_ysteps,
						const Score* c_sequence_of_scores)
{
    score           = c_score;
    condensed_steps = c_steps;
    
    condensed_sequence_of_xsteps = new int[condensed_steps+1];
    condensed_sequence_of_ysteps = new int[condensed_steps+1];
    condensed_sequence_of_states = new int[condensed_steps+1];
    condensed_sequence_of_scores = new Score[condensed_steps+1];
    
    for (int i=0; i<condensed_steps+1; i++)
    {
	condensed_sequence_of_xsteps[i]=c_sequence_of_xsteps[i];
	condensed_sequence_of_ysteps[i]=c_sequence_of_ysteps[i];
	condensed_sequence_of_states[i]=c_sequence_of_states[i];
	condensed_sequence_of_scores[i]=c_sequence_of_scores[i];
    }
    return;
}

void Hmm::delete_information_of_fake_alignment(void)
{
    score           = 0.0;
    condensed_steps = 0;
    
    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
    
    condensed_sequence_of_xsteps = NULL;
    condensed_sequence_of_ysteps = NULL;
    condensed_sequence_of_states = NULL;
    condensed_sequence_of_scores = NULL;
    
    return;
}

Score Hmm::get_log_odds_emission_score(const int reference_state,
					   const int state,
					   const int linear_index) const
{
    int check    = 0;
    
    if ((reference_state < 0) || (reference_state > (number_of_states-1))) {
	cout << "ERROR : Hmm::get_log_odds_emission_score : reference_state (" 
	     << reference_state << ") is out of range [0," << number_of_states-1 << "].\n" << flush;
	check++;
    }
    else {
	// reference state must read one letter at a time
	
	if ((model[reference_state].get_letters_to_read_x() + model[reference_state].get_letters_to_read_y()) != 1) {
	    cout << "ERROR : Hmm::get_log_odds_emission_score : reference_state (" 
		 << reference_state << ") does not read 1 letter at a time, but "
		 << (model[reference_state].get_letters_to_read_x() + model[reference_state].get_letters_to_read_y()) << " letters.\n" << flush;
	    check++;
	}
	
#ifndef _EXPERIMENT 
    // reference state must have correct emission prob matrix
    
	if (model[reference_state].check_matrix_of_emission_probs() != 0) {
	    cout << "ERROR : Hmm::get_log_odds_emission_score : reference_state (" 
		 << reference_state << ") emission probs do not sum up to one.\n" << flush;
	    check++;
	}
#endif 
    }  

    if ((state < 0) || (state > (number_of_states-1))) {
	cout << "ERROR : Hmm::get_log_odds_emission_score : state (" 
	     << state << ") is out of range [0," << number_of_states-1 << "].\n" << flush;
	check++;
    }
    else {

	// state must not be silent
	
	if ((model[state].get_letters_to_read_x() + model[state].get_letters_to_read_y()) == 0) {
	    cout << "ERROR : Hmm::get_log_odds_emission_score : state (" 
		 << state << ") is silent.\n" << flush;
	    check++;
	}
	
	// state must have correct emission prob matrix
#ifndef _EXPERIMENT
	if (model[state].check_matrix_of_emission_probs() != 0) {
	    cout << "ERROR : Hmm::get_log_odds_emission_score : state (" 
		 << state << ") emission probs do not sum up to one.\n" << flush;
	    check++;
	}
#endif 
    }  
    
    if (check == 0) {
	int number_of_letters = model[state].get_letters_to_read_x() + model[state].get_letters_to_read_y();
	int max_linear_index  = static_cast<int>(pow(alphabet, number_of_letters)) - 1;
	
	if ((linear_index < 0) || (linear_index > max_linear_index)) {
	    cout << "ERROR : Hmm::get_log_odds_emission_score : linear_index (" << linear_index
		 << ") for emission_prob matrix of state (" << state << ") out of range [0,"
		 << max_linear_index << "].\n" << flush;
	    check++;
	}
    }
    
    Score log_odds_emission_score = Logzero;      
    
    if (check == 0)
    {
	int i = 0;
	
	Prob  random_emission_prob    = static_cast<Prob>(0.0);
	Prob  model_emission_prob     = model[state].get_emission_prob(linear_index);
	
	// get corresponding index for l_index 
	
	const int number_of_dimensions = model[state].get_letters_to_read_x() + model[state].get_letters_to_read_y();
	array<int> index(1);
	index.SetDimension(0, number_of_dimensions);
	
	convert_to_base(linear_index, alphabet, &index);
	
	// calculate random_emission_prob
	
	for (i=0; i<number_of_dimensions; i++) {
	    if (i == 0) {
		random_emission_prob  = model[reference_state].get_emission_prob(index.GetElement(i));}
	    else {
		random_emission_prob *= model[reference_state].get_emission_prob(index.GetElement(i));}
	}
	
	// calculate log_odds_emission_score
	
	if (random_emission_prob > 0.0)
	{
	    if (model_emission_prob > 0.0) {
		log_odds_emission_score = log(model_emission_prob / random_emission_prob);
	    }
	    else if (model_emission_prob == 0.0) {
		log_odds_emission_score = Logzero;}
	}
	else {
	    cout << "ERROR : Hmm::get_log_odds_emission_score : random_emission_probs is 0.\n" << flush;
	    log_odds_emission_score = Logzero;
	} 
    }
    return(log_odds_emission_score);
}

int Hmm::score_user_defined_state_path(const Sequence* x,
					   const Sequence* y,
					   const int  abs_start_x,
					   const int  abs_start_y,
					   const int  n_of_states,
					   const int* states,
					   // output
					   Sequence* sub_x,
					   Sequence* sub_y) 
{
    int check = 0;
    
    if (x == NULL) {
	cout << "ERROR: class Hmm::score_user_defined_state_path: Sequence x is NULL.\n" << flush;
	check++;
    }
    else {
	if (x->length() == 0) {
	    cout << "ERROR: class Hmm::score_user_defined_state_path: length of Sequence x is 0.\n" << flush;
	    check++;
	}

	if ((abs_start_x < x->get_start_position()) || (abs_start_x > x->get_end_position())) {
	    cout << "ERROR: class Hmm::score_user_defined_state_path: abs_start_x (" << abs_start_x
		 << ") out of range [" << x->get_start_position() << ", " << x->get_end_position()
		 << "].\n" << flush;
	    check++;
	}
    }
    if (y == NULL) {
	cout << "ERROR: class Hmm::score_user_defined_state_path: Sequence y is NULL.\n" << flush;
	check++;
    }
    else {
	if (y->length() == 0) {
	    cout << "ERROR: class Hmm::score_user_defined_state_path: length of Sequence y is 0.\n" << flush;
	    check++;
	}
	
	if ((abs_start_y < y->get_start_position()) || (abs_start_y > y->get_end_position())) {
	    cout << "ERROR: class Hmm::score_user_defined_state_path: abs_start_y (" << abs_start_y
		 << ") out of range [" << y->get_start_position() << ", " << y->get_end_position()
		 << "].\n" << flush;
	    check++;
	}
    }
    if ((states == NULL) || (n_of_states == 0)) {
	cout << "ERROR: class Hmm::score_user_defined_state_path: states array is NULL or n_of_states ("
	     << n_of_states << ") is 0.\n" << flush;
	check++;
    }
    else {
	int i = 0;
	for (i=0; i<n_of_states; i++) {
	    if ((states[i] < 0) || (states[i] > (number_of_states-1))) {
		cout << "ERROR: class Hmm::score_user_defined_state_path: states[" << i 
		     << "] is out of range [0, " << number_of_states-1 << "].\n" << flush;
		check++;
		break;
	    }
	}
    }
    
    // initialise output variables
    
    if (check == 0) {

	const int shift_index = bitshift(alphabet);
	
	int i = 0;
	int j = 0;
	int k = 0;
	int p = 0;
	
	// get subsequences for indicated state path
	// note: this are the absolute positions !
	
	int subseq_start_x = abs_start_x;
	int subseq_end_x   = subseq_start_x - 1; // removed x->get_orientation
	int subseq_start_y = abs_start_y;
	int subseq_end_y   = subseq_start_y - 1; // removed y->get_orientation
	
	int last_state_index = n_of_states-1;
	
	for (i=0; i<n_of_states; i++) {
	    
	    subseq_end_x += this->model[states[i]].get_letters_to_read_x();
	    subseq_end_y += this->model[states[i]].get_letters_to_read_y();
	    

	    if ((subseq_end_x < x->get_start_position()) ||
		(subseq_end_x > x->get_end_position())   ||
		(subseq_end_y < y->get_start_position()) ||
		(subseq_end_y > y->get_end_position()))     {

		last_state_index = i-1;
		subseq_end_x -= x->get_start_position() * this->model[states[i]].get_letters_to_read_x();
		subseq_end_y -= y->get_start_position() * this->model[states[i]].get_letters_to_read_y();

		break;
	    }
	}
	
	Sequence subseq_x;
	
	if (check == 0) {
	    check += x->get_subsequence(subseq_start_x, subseq_end_x, &subseq_x);
	    if (check != 0) {
		cout << "ERROR : class Hmm::score_user_defined_state_path: "
		     << "occurred in function Sequence::get_subsequence for sequence x.\n" << flush;
		check++;
	    }
	}
	
	Sequence subseq_y;
	
	if (check == 0) {
	    check += y->get_subsequence(subseq_start_y, subseq_end_y, &subseq_y);
	    if (check != 0) {
		cout << "ERROR : class Hmm::score_user_defined_state_path: "
		     << "occurred in function Sequence::get_subsequence for sequence y.\n" << flush;
		check++;
	    }
	}
	
	// score user-defined state path 
	
	if (check == 0) {
	    // - delete any previously existing solution
	    
	    this->reset_variables_for_state_path();
	    
	    // - initialise arrays for solution
	    
	    const int max_length_state_path = subseq_x.length() + subseq_y.length() + 2;
	    
	    sequence_of_scores = new Score[max_length_state_path];
	    sequence_of_states = new int[max_length_state_path];
	    sequence_of_xsteps = new int[max_length_state_path];
	    sequence_of_ysteps = new int[max_length_state_path];
	    
	    for (i=0; i<max_length_state_path; i++) {
		
		sequence_of_scores[i] = 0;
		sequence_of_states[i] = 0;
		sequence_of_xsteps[i] = 0;
		sequence_of_ysteps[i] = 0;
	    }
	    
	    int  last_state    	   = 0;
	    int  next_state    	   = 0;
	    int  next_delta_x  	   = 0;
	    int  next_delta_y  	   = 0;
	    int  rel_pos_x     	   = 0;
	    int  rel_pos_y     	   = 0;
	    int  last_rel_pos_x    = 0;
	    int  last_rel_pos_y    = 0;
	    int* indices           = NULL;
	    
	    int      linear_index_emission = 0;
	    Score    transition_score = Logzero;
	    Score    emission_score   = Logzero;

	    int label_x = 0;
	    int label_y = 0;

	    for (i=0; i<last_state_index+1; i++) {
		
		// ********************************************************************
		
		next_state       = states[i];
		transition_score = Logzero;
		emission_score   = Logzero;
		
		next_delta_x = this->model[next_state].get_letters_to_read_x();
		next_delta_y = this->model[next_state].get_letters_to_read_y();

		rel_pos_x += next_delta_x;
		rel_pos_y += next_delta_y;
		
		indices = new int[next_delta_x+next_delta_y+1];
		
		indices[0] = last_state;
		for (k=0; k<next_delta_x; k++) {
		    indices[k+1] = subseq_x.letter(rel_pos_x - (next_delta_x-k));

		}
		for (j=next_delta_x; j<next_delta_x+next_delta_y; j++) {
		    indices[j+1] = subseq_y.letter(rel_pos_y - (next_delta_x+next_delta_y-j));

		}
		
		linear_index_emission = indices[1];
		for (p=2; p<this->model[next_state].get_number_of_dimensions_of_emission_scores()+1; p++) {
		    linear_index_emission = (linear_index_emission << shift_index) | indices[p];
		}
		
		last_rel_pos_x = rel_pos_x - next_delta_x;
		last_rel_pos_y = rel_pos_y - next_delta_y;
		
		transition_score = this->new_get_transition_score(NULL, 
							    last_state, next_state, 
							    &subseq_x, last_rel_pos_x,
							    &subseq_y, last_rel_pos_y);
 
		if (this->model[next_state].get_letters_to_read() >0) {
	  
		    emission_score   = this->model[next_state].get_emission_score(&subseq_x, last_rel_pos_x,
										  &subseq_y, last_rel_pos_y,
										  linear_index_emission);
		}	
		else {
		    emission_score   = 0.0;
		}
		
		if (indices) delete [] indices;
		indices = NULL;
		
		if ((emission_score > Logzero) && (transition_score > Logzero)) {
		    
		    steps++;  
		    sequence_of_states[steps] = next_state;
		    sequence_of_xsteps[steps] = rel_pos_x;
		    sequence_of_ysteps[steps] = rel_pos_y;
		    sequence_of_scores[steps] = emission_score + transition_score;
		    score += sequence_of_scores[steps];

	    // start new iteration
	      
		    last_state = next_state;

		}
		else {
		    cout << "ERROR : class Hmm::score_user_defined_state_path: emission_score ("	
			 << emission_score << ") or transition_score (" << transition_score 
			 << ") is Logzero (" << Logzero << ") for transition from state "
			 << last_state << " to " << next_state << " going from last_rel_pos_x ("
			 << last_rel_pos_x << ") and last_rel_pos_y (" << last_rel_pos_y 
			 << ") to rel_pos_x (" << rel_pos_x << ") and rel_pos_y (" << rel_pos_y << ").\n" << flush;
		    check++;
		    break;
		}
	  
		// **********************************************************************
		
	    } // loop over states in array states
	    
	    // go to End state
	    
	    if (check == 0) {
		
		next_state = number_of_states-1;
		
		transition_score = this->new_get_transition_score(NULL, 
								  last_state, next_state,
								  &subseq_x, last_rel_pos_x, 
								  &subseq_y, last_rel_pos_y);
		if (transition_score > Logzero) {
		    
		    steps++;  
		    sequence_of_states[steps] = next_state;
		    sequence_of_xsteps[steps] = rel_pos_x;
		    sequence_of_ysteps[steps] = rel_pos_y;
		    sequence_of_scores[steps] = transition_score;
		    score += sequence_of_scores[steps];

		}
		else {
		    cout << "ERROR : class Hmm::score_user_defined_state_path: transition_score ("
			 << transition_score << ") for transferring from last_state (" << last_state
			 << ") to End_state (" << next_state << ") <= Logzero. Cannot finish alignment.\n" << flush;
		    check++;
		}
	    } // if check == 0
	} // if check == 0
	
	// delete information about state path if error occurred
	
	if (check != 0) {
	    this->reset_variables_for_state_path();
	}
	else {
	    
	    (*sub_x) = subseq_x;
	    (*sub_y) = subseq_y;
	}
    }
    return(check);
}

void Hmm::reset_variables_for_state_path(void) 
{
    steps = 0;
    score = 0;
    
    if (sequence_of_scores) delete [] sequence_of_scores;
    if (sequence_of_states) delete [] sequence_of_states;
    if (sequence_of_xsteps) delete [] sequence_of_xsteps;
    if (sequence_of_ysteps) delete [] sequence_of_ysteps;

    sequence_of_scores = NULL;
    sequence_of_states = NULL;
    sequence_of_xsteps = NULL;
    sequence_of_ysteps = NULL;

    condensed_steps = 0;
    
    if (condensed_sequence_of_scores) delete [] condensed_sequence_of_scores;
    if (condensed_sequence_of_states) delete [] condensed_sequence_of_states;
    if (condensed_sequence_of_xsteps) delete [] condensed_sequence_of_xsteps;
    if (condensed_sequence_of_ysteps) delete [] condensed_sequence_of_ysteps;
    
    condensed_sequence_of_scores = NULL;
    condensed_sequence_of_states = NULL;
    condensed_sequence_of_xsteps = NULL;
    condensed_sequence_of_ysteps = NULL;
    
    return;
}

int Hmm::check_model_and_mirrored_model(const Hmm* mirror,
    //const int score_type,
    const Sequence *x, const int x_start, const int x_end,
    const Sequence *y, const int y_start, const int y_end)
{
    int check=0;
    
    // check Hmm
  
    if (mirrored != 0)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : Hmm : mirrored (" << mirrored 
	     << ") must be 0.\n" << flush;
	check+=1;
    }
    if (number_of_states<3)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : number of states = " 
	     << number_of_states << "<3.\n" << flush;
	check+=1;
    }
    if (model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : state 0 != Start \n" << flush;
	cout << "model[0].get_letters_to_read() : "<<model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : last state != End \n" << flush;
	cout<<"model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_mirrored()!=mirrored) 
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirrored << "\n" << flush;
		check+=1;
	    }
	    if (model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : number of states in state i = " << i 
		     << ", number of states = " 
		     << model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	   
	    if ((i!=0) && (i!=(number_of_states-1))&&(model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : state i = " << i << " should emit! \n" << flush;
		check+=1;
	    }
	    
	}
    }
    // check mirrored Hmm
    
    if (mirror->get_mirrored() != 1)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : mirrored (" 
	     << mirror->get_mirrored() << ") must be 1.\n" << flush;
	check+=1;
    }
    if (mirror->get_number_of_states()<3)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : number of states = " 
	     << mirror->get_number_of_states() << "<3.\n" << flush;
	check+=1;
    }
    if (mirror->model[0].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : state != Start \n" << flush;
	cout << "mirror->model[0].get_letters_to_read() : "
	     <<mirror->model[0].get_letters_to_read()<<endl;
	check+=1;
    }
    if (mirror->model[number_of_states-1].get_letters_to_read() != 0)
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : last state != End \n" << flush;
	cout<<"mirror->model["<<number_of_states-1<<"].get_letters_to_read() : "
	    <<mirror->model[number_of_states-1].get_letters_to_read()<<endl;
	check+=1;
    }  
    {
	for (int i=0; i<number_of_states; i++)
	{
	    if (mirror->model[i].get_alphabet()!=alphabet) 
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : alphabet of state i = " 
		     << i << ", alphabet = " 
		     << mirror->model[i].get_alphabet() << " != alphabet of Hmm alphabet = " 
		     << alphabet << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_mirrored()!=mirror->get_mirrored()) 
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : mirrored of state i = " 
		     << i << ", mirrored = " 
		     << mirror->model[i].get_mirrored() << " != mirrored of Hmm mirrored = " 
		     << mirror->get_mirrored() << "\n" << flush;
		check+=1;
	    }
	    if (mirror->model[i].get_number_of_states() != number_of_states)
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : number of states in state i = " 
		     << i << ", number of states = " 
		     << mirror->model[i].get_number_of_states() << " != number of states of Hmm = " 
		     << number_of_states << "\n" << flush;
		check+=1;
	    }
	    
	    if ((i!=0) && (i!=(number_of_states-1)) && (mirror->model[i].get_letters_to_read() == 0))
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : mirrored Hmm : state i = " 
		     << i << " should emit! \n" << flush;
		check+=1;
	    }
	   
	}
    }
    // check that this Hmm and the mirror Hmm belong together
    
    if (number_of_states != mirror->get_number_of_states())
    {
	cout << "ERROR: Hmm::check_model_and_mirrored_model : number_of_states (" << number_of_states 
	     << ") != number_of_states of mirrored Hmm (" << mirror->get_number_of_states() << ")\n" << flush;
	check+=1;
    }
    {
	for (int i=1; i<number_of_states-1; i++)
	{
	    if (model[i].get_letters_to_read() != mirror->model[number_of_states-1-i].get_letters_to_read())
	    {
		cout << "ERROR: Hmm::check_model_and_mirrored_model : letters_to_readof  model[" << i << "] (" 
		     << model[i].get_letters_to_read() << ") != letters_to_read of mirrored Hmm model["
		  << number_of_states-1-i << "] (" 
		     << mirror->model[number_of_states-1-i].get_letters_to_read() << ")\n" << flush;
		check+=1;
	    }
	}  
    }
    // check x_start, x_end and y_start, y_end values
 
    if ((x_end > x->length()) || (x_end<x_start))
    {
	cout << "ERROR Hmm::check_model_and_mirrored_model: x_end (" << x_end 
	     << ") out of range, may take values in [x_start (" << x_start
	     << "), x->length() (" << x->length() << ")].\n" << flush;
	check+=1;
    }
    if (x_start <0)
    {
	cout << "ERROR Hmm::check_model_and_mirrored_model: x_start (" << x_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    if ((y_end > y->length()) || (y_end<y_start))
    {
	cout << "ERROR Hmm::check_model_and_mirrored_model: y_end (" << y_end 
	     << ") out of range, may take values in [y_start (" << y_start
	     << "), y->length() (" << y->length() << ")].\n" << flush;
	check+=1;
    }
    if (y_start<0)
    {
	cout << "ERROR Hmm::check_model_and_mirrored_model: y_start (" << y_start 
	     << ") out of range (must be >0).\n" << flush;
	check+=1;
    }
    
    if (check == 0) {
	
	cout << "Hmm::check_model_and_mirrored_model\n" << flush;
	cout << "--------------------------------------------------\n" << flush;
	
    // parameters 
	
	const int direction   = 1;
	const int n_of_states = this->get_number_of_states();
	
	// variables which will only be needed to calculate alignment

	int shift_index = bitshift(alphabet);
	
	int i, j, p = 0;
	int xsteps, ysteps;
	int f_deltax, f_deltay = 0;
	int b_deltax, b_deltay = 0;
	
	Score f_emission       = Logzero;
	Score f_transition     = Logzero;
	Score f_transition_new = Logzero;
	
	Score b_emission       = Logzero;
	Score b_transition     = Logzero;
	Score b_transition_new = Logzero;
	
	int test_state, test_deltax, test_deltay;
	Score test_score;
	Score pre_test_score;
	
	Score transition_score_difference_1     = Logzero;
	Score max_transition_score_difference_1 = Logzero;
	int   max_new_state_f_1                 = 0;
	int   max_old_state_f_1                 = 0;
	int   max_x_position_f_1                = 0;
	int   max_y_position_f_1                = 0;
	
	Score transition_score_difference_2     = Logzero;
	Score max_transition_score_difference_2 = Logzero;
	int   max_new_state_f_2                 = 0;
	int   max_old_state_f_2                 = 0;
	int   max_x_position_f_2                = 0;
	int   max_y_position_f_2                = 0;
	int states=n_of_states;
	int start_x=0;
	int end_x=0;
	int start_y=0;
	int end_y=0;
	int offset=0;

	int f_new_state = 0;
	int b_new_state = 0;
	int f_old_state = 0;
	int b_old_state = 0;
	int f_number_of_old_states = 0;
	int f_old_state_number = 0;
	int f_x_position = 0;
	int f_y_position = 0;
	int b_x_position = 0;
	int b_y_position = 0;
	int f_linear_index_emission = 0;
	int b_linear_index_emission = 0;
	int* f_indices   = NULL;
	int* b_indices   = NULL;
	int* f_positions = NULL;
	int* b_positions = NULL;
	
	int y_length=y_end-y_start+1;
	
	if (direction==1) {
	    
	    offset=0;
	    start_x = x_start; 
	    end_x   = x_end+1; 
	    start_y = y_start;
	    end_y   = y_end+1; 
	}
	else if (direction==-1) { 
	    // permute start and end points if running backwards
	    
	    offset=1;	    
	    start_x = x_end-1;
	    end_x   = x_start-2;
	    start_y = y_end-1;
	    end_y   = y_start-2;
	}

	for ( xsteps=start_x; xsteps != end_x; xsteps+=direction) {

	    for ( ysteps=start_y; ysteps != end_y; ysteps+=direction) {

		for (f_new_state=1; f_new_state<states-1; f_new_state++) {
		    
		    f_deltax = this->model[f_new_state].get_letters_to_read_x();
		    f_deltay = this->model[f_new_state].get_letters_to_read_y();

		    if (((xsteps-start_x) >= f_deltax) && (ysteps-start_y   >= f_deltay)) {

			if (f_indices) delete [] f_indices;
			f_indices = NULL;
			f_indices = new int[this->model[f_new_state].get_number_of_dimensions_of_emission_scores()+1];
			
			if (b_indices) delete [] b_indices;
			b_indices = NULL;
			b_indices = new int[this->model[f_new_state].get_number_of_dimensions_of_emission_scores()+1];
			
			if (f_positions) delete [] f_positions;
			f_positions = NULL;
			f_positions = new int[this->model[f_new_state].get_number_of_dimensions_of_emission_scores()+1];
			
			if (b_positions) delete [] b_positions;
			b_positions = NULL;
			b_positions = new int[this->model[f_new_state].get_number_of_dimensions_of_emission_scores()+1];
			
			for (i=0; i<f_deltax; i++) {
			    
			    f_indices[i]   = x->letter(xsteps - direction*(f_deltax-i));
			    f_positions[i] = xsteps - direction*(f_deltax-i);

			}
			for (j=f_deltax; j<f_deltay+f_deltax; j++) {
			    
			    f_indices[j]   = y->letter(ysteps - direction*(f_deltay+f_deltax-j));
			    f_positions[j] = ysteps - direction*(f_deltay+f_deltax-j);

			}

			for (i=0; i<f_deltax+f_deltay; i++) {
			    
			    b_indices[i]   = f_indices[f_deltax+f_deltay-1-i];
			    b_positions[i] = f_positions[f_deltax+f_deltay-1-i];

			}

			f_number_of_old_states = this->model[f_new_state].get_number_of_previous_states();
			
			for (f_old_state_number=0; f_old_state_number<f_number_of_old_states; f_old_state_number++) {
			    
			    f_old_state = this->model[f_new_state].get_number_of_previous_state(f_old_state_number);

			    f_linear_index_emission = f_indices[0];
			    
			    for (p=1; p<this->model[f_new_state].get_number_of_dimensions_of_emission_scores(); p++) {
				f_linear_index_emission = (f_linear_index_emission << shift_index) | f_indices[p];
			    }

			    b_linear_index_emission = b_indices[0];
			    
			    for (p=1; p<this->model[f_new_state].get_number_of_dimensions_of_emission_scores(); p++) {
				b_linear_index_emission = (b_linear_index_emission << shift_index) | b_indices[p];
			    }
			    
			    f_x_position = xsteps - f_deltax;
			    f_y_position = ysteps - f_deltay;
			    
			    b_x_position = xsteps + 1;
			    b_y_position = ysteps + 1;
			    
			    b_new_state = this->get_number_of_states()-1 - f_old_state;
			    b_old_state = this->get_number_of_states()-1 - f_new_state;
			    
			    f_transition = this->get_transition_score(NULL, 
								      f_old_state, f_new_state, 
								      x, f_x_position, y, f_y_position);
			    f_emission   = this->model[f_new_state].get_emission_score(x, f_x_position, y, f_y_position,
										       f_linear_index_emission);

			    b_transition = mirror->get_transition_score(this, 
									b_old_state, b_new_state,
									x, b_x_position, y, b_y_position);
			    b_emission   = mirror->model[b_old_state].get_emission_score(x, b_x_position, y, b_y_position,
											 b_linear_index_emission);
			    
			    b_transition_new = mirror->new_get_transition_score(this, 
										b_old_state, b_new_state,
										x, b_x_position, y, b_y_position);
			    
			    f_transition_new = this->get_transition_score(NULL, 
									  f_old_state, f_new_state, 
									  x, f_x_position, y, f_y_position);
			    
			    if ((f_transition != b_transition)     ||
				(f_emission   != b_emission)       ||
				(f_transition != f_transition_new) ||
				(b_transition != b_transition_new) ||
				(f_transition != b_transition_new)) {
				
//		    if (f_transition != b_transition) {
//		      cout << "ERROR: f_transition != b_transition\n" << flush;
//		    }

//		      if (f_transition != b_transition) {
//			cout << "ERROR: f_transition != b_transition\n" << flush;
//		      }

//		      if (f_emission != b_emission) {
//			cout << "ERROR: f_emission != b_emission\n" << flush;
//		      }

//		    if (b_transition != b_transition_new) {
//		      cout << "ERROR: b_transition != b_transition_new\n" << flush;
//		    }

				if ((f_transition != b_transition_new) ||
				    (f_transition != f_transition_new)) {
				    
				    if (f_transition != b_transition_new) {
					cout << "ERROR: f_transition != b_transition_new\n" << flush;
					
					cout << "    f : [" << f_old_state << "] -> [" << f_new_state 
					     << "]\t (x,y) = (" << f_x_position << "," << f_y_position 
					     << ") : f_trans = " << f_transition << "\n" << flush;
					cout << "new b : [" << b_old_state << "] -> [" << b_new_state 
					     << "]\t (x,y) = (" << b_x_position << "," << b_y_position 
					     << ") : b_trans = " << b_transition_new << "\n" << flush;
				    }
				    if (f_transition != f_transition_new) {
					cout << "ERROR: f_transition != f_transition_new\n" << flush;
					
					cout << "    f : [" << f_old_state << "] -> [" << f_new_state 
			     << "]\t (x,y) = (" << f_x_position << "," << f_y_position 
					     << ") : f_trans = " << f_transition << "\n" << flush;
					cout << "new f : [" << f_old_state << "] -> [" << f_new_state 
					     << "]\t (x,y) = (" << f_x_position << "," << f_y_position 
					     << ") : f_trans = " << f_transition_new << "\n" << flush;
				    }
				}
				
				transition_score_difference_1 = abs(f_transition - f_transition_new);
				
				if ((transition_score_difference_1 > max_transition_score_difference_1) &&
				    (f_transition != Logzero) && (f_transition_new != Logzero)) {
				    
				    max_transition_score_difference_1 = transition_score_difference_1;
				    max_new_state_f_1                 = f_new_state;
				    max_old_state_f_1                 = f_old_state;
				    max_x_position_f_1                = f_x_position;
				    max_y_position_f_1                = f_y_position;
				    
				    cout << "    f : [" << f_old_state << "] -> [" << f_new_state 
					 << "]\t (x,y) = (" << f_x_position << "," << f_y_position 
					 << ") : f_trans = " << f_transition << "\n" << flush;
				    cout << "new f : [" << f_old_state << "] -> [" << f_new_state 
					 << "]\t (x,y) = (" << f_x_position << "," << f_y_position 
					 << ") : f_trans = " << f_transition_new << "\n" << flush;
				}
				
				transition_score_difference_2 = abs(f_transition - b_transition_new);
				
				if ((transition_score_difference_2 > max_transition_score_difference_2) &&
				    (f_transition != Logzero) && (b_transition_new != Logzero)) {
				    
				    max_transition_score_difference_2 = transition_score_difference_2;
				    max_new_state_f_2                 = f_new_state;
				    max_old_state_f_2                 = f_old_state;
				    max_x_position_f_2                = f_x_position;
				    max_y_position_f_2                = f_y_position;
				    
				    cout << "    f : [" << f_old_state << "] -> [" << f_new_state 
					 << "]\t (x,y) = (" << f_x_position << "," << f_y_position 
					 << ") : f_trans = " << f_transition << "\n" << flush;
				    cout << "new b : [" << b_old_state << "] -> [" << b_new_state 
					 << "]\t (x,y) = (" << b_x_position << "," << b_y_position 
					 << ") : b_trans = " << b_transition_new << "\n" << flush;
				}

			    }
			} // for loop f_old_state
		    } // for loop f_new_state
		} // if not ( (xsteps==0) && (ysteps==0) )
	    } // for loop ysteps
	} // for loop xsteps
	
	cout << "max_transition_score_difference_1 = " << max_transition_score_difference_1 << "\n" << flush;
	cout << "max_new_state_f_1                 = " << max_new_state_f_1                 << "\n" << flush;
	cout << "max_old_state_f_1                 = " << max_old_state_f_1                 << "\n" << flush;
	cout << "max_x_position_f_1                = " << max_x_position_f_1                << "\n" << flush;
	cout << "max_y_position_f_1                = " << max_y_position_f_1                << "\n" << flush;
	cout << "----------------------------------------------------------------------\n" << flush;
	cout << "max_transition_score_difference_2 = " << max_transition_score_difference_2 << "\n" << flush;
	cout << "max_new_state_f_2                 = " << max_new_state_f_2                 << "\n" << flush;
	cout << "max_old_state_f_2                 = " << max_old_state_f_2                 << "\n" << flush;
	cout << "max_x_position_f_2                = " << max_x_position_f_2                << "\n" << flush;
	cout << "max_y_position_f_2                = " << max_y_position_f_2                << "\n" << flush;
	
	if (f_indices) delete [] f_indices;
	f_indices = NULL;
	
	if (f_positions) delete [] f_positions;
	f_positions = NULL;
	
	if (b_indices) delete [] b_indices;
	b_indices = NULL;
	
	if (b_positions) delete [] b_positions;
	b_positions = NULL;
	
    } // if initial checks o.k.
    
    return(check);
}
int Hmm::plain_viterbi(const Hmm *s, 
			   const Sequence *x, 
			   const Sequence *y)
{
    int check = 0;

    if (check==0)
    {
	this->reset_variables_for_state_path();

	const int x_start = 0;
	const int x_end   = x->length();
	const int y_start = 0;
	const int y_end   = y->length();
	const int end_state = s->get_number_of_states()-1;

	int * start_states = new int[1];
	start_states[0] = 0;
	const int number_of_start_states = 1;
	
	Score*** viterbi_rectangle_array = NULL;
	
	check += this->allocate_memory_for_viterbi_rectangle(s,
							     x_start, x_end,
							     y_start, y_end,
							     1, // direction == 1
							     &viterbi_rectangle_array);
	if (check != 0) {
	    cout << "ERROR: Hmm::plain_viterbi: error occurred in function allocate_memory_for_viterbi_rectangle.\n";
	}
	if (check == 0) {
	    
	    check += this->calculate_viterbi_rectangle(viterbi_rectangle_array, 
						       s,
						       NULL,
						       x, x_start, x_end,
						       y, y_start, y_end, 
						       start_states,  
						       number_of_start_states,
						       0, 0, // x_margin, y_margin, 
						       1);   // direction

	    if (check != 0) {
		cout << "ERROR: Hmm::plain_viterbi: error occurred in function calculate_viterbi_rectangle.\n";
	    }
	}
	
	if (check == 0) {
	    
	    check += this->allocate_memory_for_state_path(x_start, x_end,
							  y_start, y_end,
							  &sequence_of_scores,
							  &sequence_of_states,
							  &sequence_of_xsteps,
							  &sequence_of_ysteps);
	    
	    if (check != 0) {
		cout << "ERROR: Hmm::plain_viterbi: error occurred in function allocate_memory_for_state_path.\n";
	    }
	}
	
	if (check == 0) {
	    
	    check += this->retrieve_state_path_from_viterbi_rectangle(viterbi_rectangle_array,
								      s, 
								      NULL,
								      x, x_start, x_end,
								      y, y_start, y_end, 
								      end_state, 
								      1, // direction
								      steps,
								      score,
								      sequence_of_scores,
								      sequence_of_states,
								      sequence_of_xsteps,
								      sequence_of_ysteps);
	    if (check != 0) {
		cout << "ERROR: Hmm::plain_viterbi: error occurred in function retrieve_state_path_from_viterbi_rectangle.\n";
	    }
	}
	
	if (check == 0) {
	    
	    cout << "score = " << score << "\n";

	    this->condense_solution(1); 
	}
	else { // if check != 0
	    
	    this->reset_variables_for_state_path();
	}
	
	this->delete_memory_for_viterbi_rectangle(s, 
						  x_start, x_end,
						  1, // direction
						  &viterbi_rectangle_array);
    }
    return(check);
}

#ifndef _INTEL

int Hmm::viterbi_tube(const Sequence &x, const Sequence &y, 
			  const vector<int> start_point, const vector<int> end_point, 
			  Tube<int> tube, StatePath &local_state_path) const
{
    // Note:
    //
    // - input tube can either be empty (than the full viterbi tube is allocated within this function) or
    //   not, in which case it has to comprise the entire length of sequence x
    // - start_point and end_point have to lie within the area defined by tube
    // - the state path starts at start_point (including emission at start and ends at end_point
    //   (including emission at the end state)
    // - at the start_point, start_state reads x letters [start_x - delta_x(start_state), start_x - 1] and
    //                                         y letters [start_y - delta_y(start_state), start_y - 1]
    // - at the end_point, end_state reads x letters [end_x - delta_x(end_state), end_x - 1] and
    //                                     y letters [end_y - delta_y(end_state), end_y - 1]
    
    
    int check = 0;
    
    if (this->get_mirrored() != 0) {
	
	cout << "ERROR: Hmm::viterbi_tube: this Hmm is mirrored. Cannot use it with this function.\n";
	check++;
    }

    // allocate default tube if input tube is empty
    // the default tube comprises the full Viterbi matrix
    
    if (tube.Empty()) {
	
	int max_y = y.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(x.length(), 2);
	
	for (int i = 0; i < x.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }
    if (tube.Lx() != x.length()) {
	
	cout << "ERROR: Hmm::viterbi_tube: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << x.length() << ").\n";
	check++;
    }
    // check that Start_x < End_x 
    
    if (start_point[0] > end_point[0]) {

	cout << "ERROR: Hmm::viterbi_tube: Start_x (" << start_point[0] << ") > End_x (" 
	     << end_point[0] << ").\n";
	check++;
    }

    if (check == 0) {
	
	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 
	
	// constants
	
	const int seq_start_x = 0;                 
	const int seq_end_x   = x.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = y.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = start_point[0];
	const int Start_y     = start_point[1];
	const int End_x       = end_point[0];
	const int End_y       = end_point[1];
	const int Start_state = start_point[2];
	const int End_state   = end_point[2];
	const int shift_index = bitshift(alphabet);
	
	cout << "Hmm::viterbi_tube: start_point[" << start_point[0] << "]";
	if(pair)
	{	  
	    cout<<"[" << start_point[1] << "]";
	}
	cout<<"["<< start_point[2]<<"] end_point["<< end_point[0]<<"]";
	if(pair)
	{
	    cout<<"[" << end_point[1] << "]";
	}
	cout<<"[" << end_point[2] << "]\n" << flush;
	
	// variables 

	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   new_state, old_state, deltax, deltay;
	int   max_state, index;
	int   end_y, old_start_y, old_end_y;
	Score test_score, new_score, max_score, old_score, emission_score, transition_score;
	
	long int l_index3, l_old_index3, l_this_index3, l_prev_index3; 
	
	vector<int > index3(3);       // x position, y position, state
	vector<int > old_index3(3);   // x position, y position, state
	vector<int > this_index3(3);  // x position, y position, state
	vector<int > prev_index3(3);  // x position, y position, state
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) {

	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate

	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }

	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));

	    }

	}

	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	int* offset_x = new int[End_x - Start_x + 1];
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	long int l = 0;

	for (i=(Start_x-1); i <= (End_x-1); i++) {

	    offset_x[i-(Start_x-1)] = l;

	    if (i < 0){
		if (tube.GetElement(0, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		}
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }

	    l += (int_tube.GetElement(i-(Start_x-1), 1) - int_tube.GetElement(i-(Start_x-1), 0) + 1) * states;

	}
	l+= (End_y-int_tube.GetElement(i-1-(Start_x-1),0))*states;
	
	const long int tube_area = l;

	// allocate memory for viterbi_tube

	ViterbiElement **viterbi_tube = NULL;
	viterbi_tube = new ViterbiElement*[tube_area];

	for (i = 0; i < tube_area; i++) 
	{
	    viterbi_tube[i] = NULL;
	}

	// calculate max_deltax value 

	int max_d = 0;

	for (i = 1; i < (states-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;

	// ----------------------------------------------------------------------
	// start : calculation of viterbi_tube
	// ----------------------------------------------------------------------
	
	last_state = states-2;

	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) {

	    index3[0] = xsteps;
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }

	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) {

		index3[1] = ysteps;

		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0)) 
		{
		    last_state = states-1;
		}

		if ((xsteps == Start_x) && (ysteps == Start_y)) { // initialise start of state path

		    index3[2] = Start_state;
		    l_index3  = offset_x[index3[0]-Start_x] + (index3[1] - start_y) * states + index3[2]; 
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();

		    // get emission score for Start_state
		    
		    if ((deltax + deltay) > 0) {
		  
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               {
			    linear_index_emission = (linear_index_emission << shift_index) | x.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) {
			    linear_index_emission = (linear_index_emission << shift_index) | y.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 
			
			emission_score = this->model[Start_state].get_emission_score(&x, xsteps - deltax, &y, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    viterbi_tube[l_index3]        = new ViterbiElement; 
		    viterbi_tube[l_index3]->score = emission_score;     

		}
		else {

		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			
			max_score      = Logzero;
			emission_score = Logzero; 
			
			index3[2]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();

			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) {                                                                      
			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_index3[0] = old_xsteps;
				old_index3[1] = old_ysteps;
				
				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | x.letter(old_xsteps + i);
				    }
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | y.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;} 
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;} 
				    
				    emission_score = this->model[new_state].get_emission_score(&x, old_xsteps, &y, old_ysteps, 
											       linear_index_emission);
				}
				else {emission_score = 0.0;}
				
				// loop over old states, if emission_score in new_state is larger than Logzero
			
				if (emission_score > Logzero) { 				 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) {
					
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index3[2] = old_state;
					l_old_index3  = offset_x[old_index3[0]-Start_x] + (old_index3[1] - old_start_y) * states + old_index3[2]; 
					
					// if viterbi element of old_state exists

					if (viterbi_tube[l_old_index3] != NULL) 
					{					 
					    					     
					    transition_score = this->new_get_transition_score(this, old_state, new_state, &x, old_xsteps, &y, old_ysteps);					   
					    new_score = transition_score + emission_score;
					    old_score = viterbi_tube[l_old_index3]->score;		

					    if ((transition_score > Logzero) && (old_score > Logzero)) 
					    {
			
						test_score = old_score + new_score;
						
						if (test_score > max_score) 
						{			  
						    max_score  = test_score;
						    max_state  = old_state;						   
						}
						
					    } 
					} // if viterbi element of old_state exists
					
				    } // loop over old states
				   
				    if (max_score > Logzero) {

					l_index3  = offset_x[index3[0]-Start_x] + (index3[1] - start_y) * states + index3[2]; 
					
					
					//cout<<l_index3<<endl;

					viterbi_tube[l_index3]            = new ViterbiElement; 
					viterbi_tube[l_index3]->score     = max_score;          
					viterbi_tube[l_index3]->prevState = max_state;          

					// update nParent number of old_state from which max was derived
		  
					old_index3[2] = max_state;
					l_old_index3  = offset_x[old_index3[0]-Start_x] + (old_index3[1] - old_start_y) * states + old_index3[2]; 
					viterbi_tube[l_old_index3]->nParent++; 
				
				    } // if max_score > Logzero					    
				} // if emission_score in new_state > Logzero				
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	} // loop over xsteps

	// ----------------------------------------------------------------------
	// end: calculation of viterbi_tube
	// ----------------------------------------------------------------------
    
	// ----------------------------------------------------------------------
	// start : traceback
	// ----------------------------------------------------------------------
    
	// start traceback at end_point

	index3[0] = End_x;
	index3[1] = End_y;
	index3[2] = End_state;    

	l_index3  = offset_x[index3[0]-Start_x] + (index3[1] - int_tube.GetElement(index3[0]-Start_x, 0)) * states + index3[2]; 

	if (viterbi_tube[l_index3] == NULL) {

	    cout << "ERROR: Hmm::viterbi_tube: end of state path [" << End_x << "][" << End_y 
		 << "][" << End_state << "] was not reached.\n";
	    check++;
	}
	else {

	    local_state_path.l_score = viterbi_tube[l_index3]->score; 

	    local_state_path.l_scores = new Score[x.length()+y.length()+2];
	    local_state_path.l_states = new   int[x.length()+y.length()+2];
	    local_state_path.l_xsteps = new   int[x.length()+y.length()+2];
	    local_state_path.l_ysteps = new   int[x.length()+y.length()+2];

	    this_index3   = index3;
	    l_this_index3 = l_index3; 
      
	    index = 0; 

	    while ((viterbi_tube[l_this_index3]->nParent   >= 0) &&  
		   (viterbi_tube[l_this_index3]->prevState >= 0)) { 

		prev_index3[0] = this_index3[0] - this->model[this_index3[2]].get_letters_to_read_x();
		prev_index3[1] = this_index3[1] - this->model[this_index3[2]].get_letters_to_read_y();
		prev_index3[2] = viterbi_tube[l_this_index3]->prevState; 
		l_prev_index3  = offset_x[prev_index3[0]-Start_x] + (prev_index3[1] - int_tube.GetElement(prev_index3[0]-Start_x, 0)) * states 
		    + prev_index3[2]; 
	
		local_state_path.l_xsteps[index] = this_index3[0];
		local_state_path.l_ysteps[index] = this_index3[1];
		local_state_path.l_states[index] = this_index3[2];
		local_state_path.l_scores[index] = viterbi_tube[l_this_index3]->score - viterbi_tube[l_prev_index3]->score;
		
		index++; 
		viterbi_tube[l_prev_index3]->nParent--; 
	
		this_index3   = prev_index3;
		l_this_index3 = l_prev_index3; 

	    } // while loop
    
	    // finish traceback at start_point

	    if ((viterbi_tube[l_this_index3]->nParent   >=  0) &&
		(viterbi_tube[l_this_index3]->prevState == -1) && 
		(this_index3[0] == Start_x) && 
		(this_index3[1] == Start_y) && 
		(this_index3[2] == Start_state)) {
      
		local_state_path.l_xsteps[index] = this_index3[0];
		local_state_path.l_ysteps[index] = this_index3[1];
		local_state_path.l_states[index] = this_index3[2];
		local_state_path.l_scores[index] = viterbi_tube[l_this_index3]->score;

		delete viterbi_tube[l_this_index3]; 
		viterbi_tube[l_this_index3] = NULL; 
	    }
	    local_state_path.l_steps = index; 

	    // reverse order of states etc. to get state path starting at start_point and finishing at end_point

	    for (i=0; i<((float)local_state_path.l_steps+1.0)/2.0; i++) {
		
		int   c_x = local_state_path.l_xsteps[i];
		int   c_y = local_state_path.l_ysteps[i];
		int   c_s = local_state_path.l_states[i];
		Score c_S = local_state_path.l_scores[i];
		
		local_state_path.l_xsteps[i] = local_state_path.l_xsteps[local_state_path.l_steps-i];
		local_state_path.l_ysteps[i] = local_state_path.l_ysteps[local_state_path.l_steps-i];
		local_state_path.l_states[i] = local_state_path.l_states[local_state_path.l_steps-i];
		local_state_path.l_scores[i] = local_state_path.l_scores[local_state_path.l_steps-i];
		
		local_state_path.l_xsteps[local_state_path.l_steps-i] = c_x;
		local_state_path.l_ysteps[local_state_path.l_steps-i] = c_y;
		local_state_path.l_states[local_state_path.l_steps-i] = c_s;
		local_state_path.l_scores[local_state_path.l_steps-i] = c_S;
	    }
	    
	    
	    Score sum_score = 0.0;
	    
	} // if end of state path could be reached

	// ----------------------------------------------------------------------
	// end: traceback
	// ----------------------------------------------------------------------
	
	// free memory
	
	if (offset_x) delete [] offset_x;   
	offset_x     = NULL;

	if (viterbi_tube) {
	    for (i = 0; i < tube_area; i++) {
		if (viterbi_tube[i]) delete viterbi_tube[i];
		viterbi_tube[i] = NULL;
	    }
	    if (viterbi_tube) delete [] viterbi_tube; 
	    viterbi_tube = NULL;   
	}
    }
    
    return(check);
}

#endif // #ifndef _INTEL

#ifndef _INTEL

int Hmm::viterbi_tube_strip(const Sequence &x, const Sequence &y, 
				const vector<int> start_point, const int End_x, 
				Tube<int> tube,
				ViterbiElement ****viterbi_strip) const {
				
    // Note:
    //
    // - this Hmm has to be unmirrored
    // - input tube can either be empty (than the full viterbi tube is allocated within this function) or
    //   not, in which case it has to comprise the entire length of sequence x
    // - start_point has to lie within the area defined by tube
    // - the state path starts at start_point (including emission at the start state)
    // - at the start_point, start_state reads x letters [start_x - delta_x(start_state), start_x - 1] and
    //                                         y letters [start_y - delta_y(start_state), start_y - 1]
    // - the last x letter read at End_x is    x letter  [End_x - 1]
    
    int check = 0;
    int direction = 0;
    
    if      (this->get_mirrored() == 0) {direction = 1;}
    else if (this->get_mirrored() == 1) {direction = -1;}

    if ((this->get_mirrored() != 0) and (this->get_mirrored() != 1)) {
	
	cout << "ERROR: Hmm::viterbi_tube_strip: this Hmm is neither mirrored nor unmirrored.\n";
	check++;
    }

    // allocate default tube if input tube is empty
    // the default tube comprises the full Viterbi matrix
    
    if (tube.Empty()) {
	
	const int max_y = y.length()-1;
	Tube<int> default_tube(x.length(), 2);
	
	for (int i = 0; i < x.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }
    if (tube.Lx() != x.length()) {
	
	cout << "ERROR: Hmm::viterbi_tube_strip: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << x.length() << ").\n";
	check++;
    }
  
    if (check == 0) {

	// calculate max_deltax value 
	
	int max_d = 0;
	for (int i = 1; i < (this->get_number_of_states()-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	
	int last_x_pos = 0; 
	int last_y_pos = 0; 
	
	// constants
	
	const int max_deltax  = max_d;
	const int seq_start_x = 0;                 
	const int seq_end_x   = x.length()-1;      
	const int seq_start_y = 0;                 
	const int seq_end_y   = y.length()-1;      
	const int states      = this->get_number_of_states(); 
	const int Start_x     = start_point[0];
	const int Start_y     = start_point[1];
	const int Start_state = start_point[2];
	const int shift_index = bitshift(alphabet);
	const int tube_width  = max_deltax+1;
	const int tube_height = y.length()+1;	

	cout << "Hmm::viterbi_tube_strip: Start_point[" << start_point[0] << "][" 
	     << start_point[1] << "][" << start_point[2] << "] End_x = " << End_x << "\n" << flush;
	
	// variables 
	
	int   j, i, k, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   new_state, old_state, deltax, deltay;
	int   max_state, index;
	int   end_y, old_start_y, old_end_y;
	Score test_score, new_score, max_score, old_score, emission_score, transition_score;
	
	vector<int > index3(3);       // x position, y position, state
	vector<int > old_index3(3);   // x position, y position, state
	vector<int > this_index3(3);  // x position, y position, state
	vector<int > prev_index3(3);  // x position, y position, state
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) {

	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0)));// modify tube according to start_y coordinate
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) {
		tube.SetElement(i, 0, tube.GetElement(i, 1));

	    }

	}
	
	// calculate int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	long int l = 0;
	
	for (i=(Start_x-1); i <= (End_x-1); i++) {
	    
	    if (i < 0){
		if (tube.GetElement(0, 0) == 0) {int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));}
		else                 {int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));}
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));}	
		else                 {int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }

	}
	
	// ----------------------------------------------------------------------
	// start : calculation of viterbi_strip
	// ----------------------------------------------------------------------
	
	last_state = states-2;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) {
	    
	    index3[0] = max_deltax;
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
      
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) {
		
		index3[1] = ysteps;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) {last_state = states-1;}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) { // initialise start of state path
		    
		    index3[2] = Start_state;
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state

		    if ((deltax + deltay) > 0) {
		  
			linear_index_emission = 0;

			for (i = 0; i < deltax; i++)  {
			    linear_index_emission = (linear_index_emission << shift_index) | x.letter(xsteps - deltax + i);

			}

			for (j = 0; j < deltay; j++) {
			    linear_index_emission = (linear_index_emission << shift_index) | y.letter(ysteps - deltay + j);

			}

			if ((xsteps-1) > last_x_pos) {last_x_pos = xsteps-1;} 
			if ((ysteps-1) > last_y_pos) {last_y_pos = ysteps-1;} 
			
			emission_score = this->model[Start_state].get_emission_score(&x, xsteps - deltax, &y, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    
		    viterbi_strip[index3[0]][index3[1]][index3[2]]        = new ViterbiElement; 
		    viterbi_strip[index3[0]][index3[1]][index3[2]]->score = emission_score;           

		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) {
			
			max_score      = Logzero;
			emission_score = Logzero; 
			
			index3[2]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) {                                                                      
			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_index3[0] = max_deltax - deltax;
				old_index3[1] = old_ysteps;
				
				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | x.letter(old_xsteps + i);
				    }
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | y.letter(old_ysteps + j);
				    }
				    if ((old_xsteps+deltax-1) > last_x_pos) {last_x_pos = old_xsteps+deltax-1;} 
				    if ((old_ysteps+deltay-1) > last_y_pos) {last_y_pos = old_ysteps+deltay-1;} 
				    
				    emission_score = this->model[new_state].get_emission_score(&x, old_xsteps, &y, old_ysteps, 
											       linear_index_emission);
				}
				else {emission_score = 0.0;}
				
				// loop over old states, if emission_score in new_state is larger than Logzero
				
				if (emission_score > Logzero) { 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) {
					
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index3[2] = old_state;
					
					// if viterbi element of old_state exists
					
					if (viterbi_strip[old_index3[0]][old_index3[1]][old_index3[2]] != NULL) {
					    
					    transition_score = this->new_get_transition_score(this, old_state, new_state, &x, old_xsteps, &y, old_ysteps);
					    new_score        = transition_score + emission_score;
					    old_score        = viterbi_strip[old_index3[0]][old_index3[1]][old_index3[2]]->score;
					    
					    if ((transition_score > Logzero) && (old_score > Logzero)) {
						
						test_score = old_score + new_score;
						
						if (test_score > max_score) {
						    
						    max_score  = test_score;
						    max_state  = old_state;
						}
					    } 
					} // if viterbi element of old_state exists
				    } // loop over old states
				    
				    if (max_score > Logzero) {
					
					viterbi_strip[index3[0]][index3[1]][index3[2]]            = new ViterbiElement; 
					viterbi_strip[index3[0]][index3[1]][index3[2]]->score     = max_score;          
					viterbi_strip[index3[0]][index3[1]][index3[2]]->prevState = max_state;          
					

					// update nParent number of old_state from which max was derived
					
					old_index3[2] = max_state;
					viterbi_strip[old_index3[0]][old_index3[1]][old_index3[2]]->nParent++; 
					
				    } // if max_score > Logzero		
				} // if emission_score in new_state > Logzero
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    
	    // copy columns in strip one to the left, initialise new right-most column
	    
	    if (xsteps != End_x) {
		
		for (i = 0; i < tube_width; i++) {
		    
		    if (i == 0) { // delete all elements in left-most column
			
			for (j = 0; j < tube_height; j++) {
			    for (k = 0; k < states; k++) {
				if (viterbi_strip[i][j][k]) {
				    
				    delete viterbi_strip[i][j][k];
				    viterbi_strip[i][j][k] = NULL;
				}
			    }
			}
		    }
		    else { // copy elements in column i to column i-1, then delete elements in column i
			
			for (j = 0; j < tube_height; j++) {
			    for (k = 0; k < states; k++) {
				
				if (viterbi_strip[i][j][k]) {
				    
				    if (viterbi_strip[i-1][j][k] == NULL) {
					
					viterbi_strip[i-1][j][k] = new ViterbiElement; 
		    
					viterbi_strip[i-1][j][k]->prevState = viterbi_strip[i][j][k]->prevState;
					viterbi_strip[i-1][j][k]->nParent   = viterbi_strip[i][j][k]->nParent;
					viterbi_strip[i-1][j][k]->score     = viterbi_strip[i][j][k]->score;
					
					delete viterbi_strip[i][j][k];
					viterbi_strip[i][j][k] = NULL;
				    }
				    else {
					cout << "ERROR: viterbi_strip[" << i-1 << "][" << j << "][" << k << "] should be NULL, but isn't.\n";
				    }
				}
			    }
			}
		    }
		}
	    }
	} // loop over xsteps
	
	// ----------------------------------------------------------------------
	// end: calculation of viterbi_strip
	// ----------------------------------------------------------------------
		
    }
    return(check);
}

#endif // #ifndef _INTEL

#ifndef _INTEL

int Hmm::viterbi_tube_strip_backwards(const Sequence &x, const Sequence &y, 
					  const int Start_x, const vector<int> end_point, 
					  Hmm *unmirrored_s, 					 
					  Tube<int> tube, 
					  ViterbiElement ****viterbi_strip) const {

    // Note:
    //
    // - function can only be called by a mirrored pairhmme
    // - input tube can either be empty (than the full viterbi tube is allocated within this function) or
    //   not, in which case it has to comprise the entire length of sequence x
    // - start_point has to lie within the area defined by tube
    // - the state path starts at end_point (including emission at the end state)
    // - at the end_point, end_state reads    x letters  [end_x + 1, end_x + delta_x(end_state)] and
    //                                        y letters  [end_y + 1, end_y + delta_y(end_state)]
    // - the 5-'most x letter read at Start_x is x letter [Start_x + 1]
    
    int check = 0;
    
    if (this->get_mirrored() == 0) {
	
	cout << "ERROR: Hmm::viterbi_tube_strip_backwards: this Hmm is unmirrored.\n";
	check++;
    }
    
    // allocate default tube if input tube is empty
    // the default tube comprises the full Viterbi matrix
    
    if (tube.Empty()) {
	
	const int max_y = y.length()-1;
	Tube<int> default_tube(x.length(), 2);
	
	for (int i = 0; i < x.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }
    if (tube.Lx() != x.length()) {
	
	cout << "ERROR: Hmm::viterbi_tube_strip_backwards: length of tube (" << tube.Lx() 
	     << ") != length of sequence X ("
	     << x.length() << ").\n";
	check++;
    }

    if (check == 0) {

	int max_d = 0;
	
	for (int i = 1; i < (this->get_number_of_states()-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	
	// constants

	const int max_deltax  = max_d;
	const int seq_start_x = 0;                 
	const int seq_end_x   = x.length()-1;      
	const int seq_start_y = 0;                 
	const int seq_end_y   = y.length()-1;      
	const int states      = this->get_number_of_states(); 
	const int End_x       = end_point[0]; 
	const int End_y       = end_point[1]; 
	const int End_state   = states - 1 - end_point[2];
	const int shift_index = bitshift(alphabet);
	const int tube_width  = max_deltax+1;
	const int tube_height = y.length()+1;
	
	cout << "Hmm::viterbi_tube_strip_backwards: Start_x = " << Start_x << " End_point[" 
	     << end_point[0] << "][" << end_point[1] << "][" << end_point[2] << "]\n" << flush;
	

	int last_x_pos = End_x; 
	int last_y_pos = End_y; 
	
	const long int tube_volume = tube_width * tube_height * states;

        // variables 

	int   j, i, k, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   new_state, old_state, deltax, deltay;
	int   max_state, index;
	int   start_x, start_y, old_start_y, old_end_y;
	Score test_score, new_score, max_score, old_score, emission_score, transition_score;

	vector<int > index3(3);       // x position, y position, state
	vector<int > old_index3(3);   // x position, y position, state
	vector<int > this_index3(3);  // x position, y position, state
	vector<int > prev_index3(3);  // x position, y position, state
	
	int end_y = End_y;                              // most 3' letter read from y
	int d_y   = this->model[End_state].get_letters_to_read_y();
	
	if (d_y > 0) {end_y = min(end_y + d_y, seq_end_y);}
	else         {end_y = min(end_y + 1, seq_end_y);}
	
	for (i=Start_x+1; i <= End_x; i++) {
	    
	    tube.SetElement(i, 1, min(end_y, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) {
		tube.SetElement(i, 1, tube.GetElement(i, 0));

	    }

	}

	// calculate int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]

	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	long int l = 0;
	
	for (i=(Start_x+1); i <= (End_x+1); i++) {
	    
	    if (i > seq_end_x){
		int_tube.SetElement(i-(Start_x+1), 0, max(tube.GetElement(seq_end_x, 0)-1, -1));
		if (tube.GetElement(seq_end_x, 1) == seq_end_y) {
		    int_tube.SetElement(i-(Start_x+1), 1, max(tube.GetElement(seq_end_x, 1), -1));
		}
		else {
		    int_tube.SetElement(i-(Start_x+1), 1, max(tube.GetElement(seq_end_x, 1)-1, -1));
		}
	    }
	    else {
		int_tube.SetElement(i-(Start_x+1), 0, max(tube.GetElement(i, 0)-1, -1));
		if (tube.GetElement(i, 1) == seq_end_y) {
		    int_tube.SetElement(i-(Start_x+1), 1, max(tube.GetElement(i, 1), -1));
		}
		else {
		    int_tube.SetElement(i-(Start_x+1), 1, max(tube.GetElement(i, 1)-1, -1));
		}
	    }

	}
	
	// ----------------------------------------------------------------------
	// start : calculation of viterbi_strip
	// ----------------------------------------------------------------------
    
	last_state  = states-2;
	
	for (xsteps = End_x; xsteps >= Start_x; xsteps -= 1) {
	    
	    index3[0] = 0;
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    	    
	    for (ysteps = end_y; ysteps >= start_y; ysteps -= 1) {

		index3[1] = ysteps+1;


		if (((xsteps+1) == seq_start_x) && ((ysteps+1) == seq_start_y)) { 

		    last_state = states - 1;

	}

		if ((xsteps == End_x) && (ysteps == End_y)) { // initialise start of state path

		    index3[2] = End_state;
		    deltax    = this->model[End_state].get_letters_to_read_x();
		    deltay    = this->model[End_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state

		    if ((deltax + deltay) > 0) {
		  
			linear_index_emission = 0;

			for (j = deltay; j > 0; j--) {
			    linear_index_emission = (linear_index_emission << shift_index) | y.letter(ysteps + j);

			}

			for (i = deltax; i > 0; i--)               {
			    linear_index_emission = (linear_index_emission << shift_index) | x.letter(xsteps + i);

			}

			if ((xsteps+1) < last_x_pos) {last_x_pos = xsteps+1;} 
			if ((ysteps+1) < last_y_pos) {last_y_pos = ysteps+1;} 

			emission_score = this->model[End_state].get_emission_score(&x, xsteps + 1, &y, ysteps + 1, 
										   linear_index_emission);
		    }
		    else {emission_score = 0.0;} 

		    viterbi_strip[index3[0]][index3[1]][index3[2]]        = new ViterbiElement; 
		    viterbi_strip[index3[0]][index3[1]][index3[2]]->score = emission_score;           

		}
		else {

		    for (new_state = 1; new_state <= last_state; new_state++) {

			max_score      = Logzero;
			emission_score = Logzero; 
			
			index3[2]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps + deltax;
			old_ysteps     = ysteps + deltay;

			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps + 1)) && ((xsteps + deltax) <= seq_end_x)) || 
			     (deltax == 0))) {                                   

			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);

			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps + 1)) && ((ysteps + deltay) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      

				old_index3[0] = deltax;
				old_index3[1] = old_ysteps+1;
				
				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;

				    for (j = deltay; j > 0; j--) {
					linear_index_emission = (linear_index_emission << shift_index) | y.letter(ysteps + j);

				    }

				    for (i = deltax; i > 0; i--)               {
					linear_index_emission = (linear_index_emission << shift_index) | x.letter(xsteps + i);

				    }

				    if ((xsteps+1) < last_x_pos) {last_x_pos = xsteps+1;} 
				    if ((ysteps+1) < last_y_pos) {last_y_pos = ysteps+1;} 

				    emission_score = this->model[new_state].get_emission_score(&x, xsteps + 1, &y, ysteps + 1, 
											       linear_index_emission);
				}
				else {emission_score = 0.0;}

				// loop over old states, if emission_score in new_state is larger than Logzero

				if (emission_score > Logzero) { 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) {
					
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index3[2] = old_state;

					// if viterbi element of old_state exists
					
					if (viterbi_strip[old_index3[0]][old_index3[1]][old_index3[2]] != NULL) {

					    transition_score = this->new_get_transition_score(unmirrored_s, old_state, new_state, 
											      &x, xsteps + 1, &y, ysteps + 1);
					    new_score        = transition_score + emission_score;
					    old_score        = viterbi_strip[old_index3[0]][old_index3[1]][old_index3[2]]->score;

					    if ((transition_score > Logzero) && (old_score > Logzero)) {
						
						test_score = old_score + new_score;
						
						if (test_score > max_score) {
						    
						    max_score  = test_score;
						    max_state  = old_state;

						}
					    } 
					} // if viterbi element of old_state exists
				    } // loop over old states
				    
				    if (max_score > Logzero) {
					
					viterbi_strip[index3[0]][index3[1]][index3[2]]            = new ViterbiElement; 
					viterbi_strip[index3[0]][index3[1]][index3[2]]->score     = max_score;          
					viterbi_strip[index3[0]][index3[1]][index3[2]]->prevState = max_state;          

					// update nParent number of old_state from which max was derived
					
					old_index3[2] = max_state;
					viterbi_strip[old_index3[0]][old_index3[1]][old_index3[2]]->nParent++; 
					
				    } // if max_score > Logzero		
				} // if emission_score in new_state > Logzero
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    

	    // copy columns in strip one to the right, initialise right-most column
	    
	    if (xsteps != Start_x) {
		
		for (i = tube_width-1; i >= 0; i--) {
		    
		    if (i == (tube_width-1)) { // delete all elements in right-most column
			
			for (j = 0; j < tube_height; j++) {
			    for (k = 0; k < states; k++) {
				if (viterbi_strip[i][j][k]) {
				    
				    delete viterbi_strip[i][j][k];
				    viterbi_strip[i][j][k] = NULL;
				}
			    }
			}
		    }
		    else { // copy elements in column i to column i+1, then delete elements in column i
			
			for (j = 0; j < tube_height; j++) {
			    for (k = 0; k < states; k++) {
				
				if (viterbi_strip[i][j][k]) {
				    
				    if (viterbi_strip[i+1][j][k] == NULL) {
					
					viterbi_strip[i+1][j][k] = new ViterbiElement; 
		    
					viterbi_strip[i+1][j][k]->prevState = viterbi_strip[i][j][k]->prevState;
					viterbi_strip[i+1][j][k]->nParent   = viterbi_strip[i][j][k]->nParent;
					viterbi_strip[i+1][j][k]->score     = viterbi_strip[i][j][k]->score;
					
					delete viterbi_strip[i][j][k];
					viterbi_strip[i][j][k] = NULL;
				    }
				    else {
					cout << "ERROR: viterbi_strip[" << i+1 << "][" << j << "][" << k << "] should be NULL, but isn't.\n";
				    }
				}
			    }
			}
		    }
		}
	    }
	} // loop over xsteps
	
	// ----------------------------------------------------------------------
	// end: calculation of viterbi_strip
	// ----------------------------------------------------------------------
	
	// rename the state numbers so they are the same as in the forward strip
	
	int middle_point = 0;
	
	if (states%2 == 0) {middle_point = (int)((float)(states)/2.)     - 1;}
	else               {middle_point = (int)((float)(states - 1)/2.) - 1;}
	
	ViterbiElement copy_element;
	
	for (j = 0; j < tube_height; j++) {
	    
	    for (i = 0; i < tube_width; i++) {

		for (k = 0; k < middle_point+1; k++) {
		    
		    copy_element.score = Logzero;
		    
		    if (viterbi_strip[i][j][states - 1 - k]) {
			
			copy_element.score     = viterbi_strip[i][j][states - 1 - k]->score;
			copy_element.nParent   = viterbi_strip[i][j][states - 1 - k]->nParent;
			copy_element.prevState = viterbi_strip[i][j][states - 1 - k]->prevState;
			
			delete viterbi_strip[i][j][states - 1 - k];
			viterbi_strip[i][j][states - 1 - k] = NULL;
		    }
		    if (viterbi_strip[i][j][k]) {

			viterbi_strip[i][j][states - 1 - k] = new ViterbiElement;
			viterbi_strip[i][j][states - 1 - k]->score     = viterbi_strip[i][j][k]->score;
			viterbi_strip[i][j][states - 1 - k]->nParent   = viterbi_strip[i][j][k]->nParent;
			viterbi_strip[i][j][states - 1 - k]->prevState = viterbi_strip[i][j][k]->prevState;
			
			delete viterbi_strip[i][j][k];
			viterbi_strip[i][j][k] = NULL;
		    }
		    if (copy_element.score >Logzero) {
			
			viterbi_strip[i][j][k] = new ViterbiElement;
			viterbi_strip[i][j][k]->score     = copy_element.score;
			viterbi_strip[i][j][k]->nParent   = copy_element.nParent;
			viterbi_strip[i][j][k]->prevState = copy_element.prevState;
		    }	    
		}

	    }

	}

    }
    
    return(check);
}

int Hmm::combine_strips(ViterbiElement ****f_strip, ViterbiElement ****b_strip,
			    const Sequence &x, const Sequence &y,
			    const int tube_width, const int tube_height,
			    const int x_middle,
			    vector<int> &left_point, vector<int> &right_point) const {
  
    int check = 0;

    if (this->get_mirrored() != 0) {
	cout << "ERROR: Hmm::combine_strips: this Hmm is mirrored.\n";
	check++;
    }
    if (f_strip == NULL) {
	cout << "ERROR: Hmm::combine_strips: f_strip is NULL.\n";
	check++;
    }
    if (b_strip == NULL) {
	cout << "ERROR: Hmm::combine_strips: b_strip is NULL.\n";
	check++;
    }
    if (check == 0) {

	cout << "Hmm::combine_strips\n" << flush;

	left_point[0] = -1;
	left_point[1] = -1;
	left_point[2] = -1;
    
	right_point[0] = -1;
	right_point[1] = -1;
	right_point[2] = -1;
	
	const int states = this->get_number_of_states();
    
	int i, j, left, right;
	Score max_score, max_l_score, max_r_score, max_t_score, score, t_score, l_score, r_score;
	
	max_score   = Logzero;
	max_l_score = Logzero;
	max_t_score = Logzero;
	max_r_score = Logzero;

	for (i=0; i<tube_width; i++) {

	    for (j=0; j<tube_height; j++) {
		for (left=1; left<states; left++) {
		    if (f_strip[i][j][left] != NULL) {
			for (right=1; right<states; right++) {
			    if (b_strip[i][j][right] != NULL) {

				// combine f_strip[i][j][left] witb b_strip[i][j][right]
		
				l_score = f_strip[i][j][left]->score;
				r_score = b_strip[i][j][right]->score;
				t_score = this->new_get_transition_score(this, left, right, 
									 &x, x_middle+1-(tube_width-1-i),
									 &y, j);		
				if (t_score > Logzero) {
		
				    score = l_score + t_score + r_score;

				    if (score > max_score) {

					max_score      = score;
					max_l_score    = l_score;
					max_r_score    = r_score;
					max_t_score    = t_score;
					
					left_point[0]  = x_middle+1-(tube_width-1-i);
					left_point[1]  = j;
					left_point[2]  = left;
					right_point[0] = left_point[0] + this->model[right].get_letters_to_read_x();
					right_point[1] = left_point[1] + this->model[right].get_letters_to_read_y();
					right_point[2] = right;

				    }
				}
			    }
			}
		    }
		}
	    }
	}
	if ((left_point[2] == -1) || (right_point[2] == -1)) {
	    
	    cout << "ERROR: Hmm::combine_strips: did not manage to merge strips: x_middle = " << x_middle 
		 << "\n" << flush;
	    check++;
	}
	else {
	    cout << "Hmm::combine_strips: left point [" << left_point[0] << "][" << left_point[1] << "][" << left_point[2] 
		 << "] right point [" << right_point[0] << "][" << right_point[1] << "][" << right_point[2] << "]\n" << flush;
	}
	
    }
    return(check);
}

int Hmm::GetTrainingParameters(const char* cIn_dir,
			       const char* XMLfile,
			       TransitionProb* TP,
			       EmissionProb* EP,
			       model_parameters* const MP,
			       char** seqfile,
			       char** seqannfile,
			       char** name_input_tube_file,
			       long long int& max_volume,
			       double& threshold,
			       int& Maxiter,
			       double& DirConstant
    )

{
    TiXmlDocument doc(XMLfile);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    int checktp = 0;
    int checkep = 0;
   
    if(!loadOkay){
	cout << "ERROR: GetTraninigParameters : cannot open file "<<XMLfile<< ".\n" <<flush;
	check++;
	return check;
    }
    TiXmlNode * SequenceNode = 0;
    TiXmlNode * TmpNode = 0; 
    TiXmlElement* Element = 0;

    SequenceNode = doc.FirstChild("HMMConverter");
    if(!SequenceNode)
    {
        check++;
        return check;
    }

    SequenceNode = SequenceNode->FirstChild("sequence_analysis");
    if(!SequenceNode)
    {
	cout<<"No sequence_analysis tag defined, training abort."<<endl;
	return check;
    }
    TmpNode = SequenceNode->FirstChild("Algorithm");    
    if(!TmpNode)
    {
	cout<<"No Algorithm tag defined, training abort."<<endl;
	return check;
    }
    
    Element = TmpNode->FirstChildElement("training_parameters");            
    if(!Element)  // no training is needed
    {
	cout<<"No training_parameters tag defined, training abort."<<endl;
	return check;
    }else{
	// read volume	
	if(!Element->Attribute("MaxVolume"))
	{
	    cout<<"Error : pairhmm class GetTrainingParameters : volume attribute "
		<<"at training_parameters tag is NULL!"<<endl;
	    check++;
	}else {
	    max_volume = atol(Element->Attribute("MaxVolume"));
	    if((max_volume==0)&&(strcmp(Element->Attribute("MaxVolume"),"0")))
	    {	       	    
		cout<<"ERROR : Hmm:: GetTrainingParameters : "
		    <<" MaxVolume attribute in training_parameters tag"
		    <<" is invalid("<<Element->Attribute("MaxVolume")<<endl;
		check++;
	    }
	}

	// read threshold	
	if(!Element->Attribute("threshold"))
	{
	    // set default
	    threshold = 0;
	}else {
	    threshold = atof(Element->Attribute("threshold"));
	    if((threshold==0)&&(strcmp(Element->Attribute("threshold"),"0")))
	    {	       	    
		cout<<"ERROR : Hmm:: GetTrainingParameters : "
		    <<" threshold attribute in training_parameters tag"
		    <<" is invalid("<<Element->Attribute("threshold")<<endl;
		check++;
	    }
	}

	// read Maxiter
	if(!Element->Attribute("Maxiter"))
	{
	    // set default
	    Maxiter = 10;
	}else {
	    Maxiter = atoi(Element->Attribute("Maxiter"));
	    if((Maxiter==0)&&(strcmp(Element->Attribute("Maxiter"),"0")))
	    {	       	    
		cout<<"ERROR : Hmm:: GetTrainingParameters : "
		    <<" Maxiter attribute in training_parameters tag"
		    <<" is invalid("<<Element->Attribute("Maxiter")<<endl;
		check++;
	    }
	}

	// read DirConstant
	if(!Element->Attribute("DirConstant"))
	{
	    // set default
	    DirConstant = 0;
	}else {
	    DirConstant = atof(Element->Attribute("DirConstant"));
	    if((DirConstant==0)&&(strcmp(Element->Attribute("DirConstant"),"0")))
	    {	       	    
		cout<<"ERROR : Hmm:: GetTrainingParameters : "
		    <<" DirConstant attribute in training_parameters tag"
		    <<" is invalid("<<Element->Attribute("DirConstant")<<endl;
		check++;
	    }
	}
	
	// read tube
	if(Element->Attribute("tube"))
	{
	    strcpy((*name_input_tube_file),Element->Attribute("tube"));
	    if(check)
	    {
		cout<<"ERROR: Hmm:: GetTrainingParameters: "
		    <<" error in Nullstrcpy for tube attribute "
		    <<" in training_parameters tag!"<<endl;
	    }else
	    {
		char tmptubefile[Max_word_length];
		strcpy(tmptubefile,cIn_dir);
		strcat(tmptubefile,(*name_input_tube_file));		    
		strcpy((*name_input_tube_file),tmptubefile);
	    }
	} 
    }
    
    checktp += TP->get_parameters_for_training(XMLfile,MP,this);
    if(checktp)
    {
	cout<<"ERROR: Hmm:: GetTrainingParameters: error in getting TP parameters";
    }
    checkep += EP->get_parameters_for_training(XMLfile,MP,this);
    if(checkep)
    {
	cout<<"ERROR: Hmm:: GetTrainingParameters: error in getting EP parameters";
    }

    check = check + checkep + checktp;   

    if(!check)
    {

	TmpNode = SequenceNode->FirstChild("Input_file");
	Element = TmpNode->FirstChildElement("training_file_parameters");
	if(!Element)
	{
	    cout<<"ERROR: Hmm:: GetTrainingParameters: "
		<<" training field tag is NULL under Input_file tag!"<<endl;	      
	    check++;
	}else{

	    // get training input files

	    if(!Element->Attribute("SeqFile"))
	    {
		cout<<"ERROR: Hmm:: GetTrainingParameters: "
		    <<" SeqFile attribute is NULL in training_file_parameters tag!"<<endl;      
	    }else
	    {
		strcpy((*seqfile),Element->Attribute("SeqFile"));
		if(check)
		{
		    cout<<"ERROR: Hmm:: GetTrainingParameters: "
			<<" error in Nullstrcpy for SeqFile attribute "
			<<" in training_file_parameters tag!"<<endl;
		}else
		{
		    char tmpseqfile[Max_word_length];
		    strcpy(tmpseqfile,cIn_dir);
		    strcat(tmpseqfile,(*seqfile));		    
		    strcpy((*seqfile),tmpseqfile);
		}
	    }

	    // get training annotation input file

	    if(Element->Attribute("AnnFile"))
	    {
		strcpy((*seqannfile),Element->Attribute("AnnFile"));
		if(check)
		{
		    cout<<"ERROR: Hmm:: GetTrainingParameters: "
			<<" error in Nullstrcpy for AnnFile attribute "
			<<" in training_file_parameters tag!"<<endl;
		}else
		{
		    char tmpannfile[Max_word_length];
		    strcpy(tmpannfile,cIn_dir);
		    strcat(tmpannfile,(*seqannfile));		   
		    strcpy((*seqannfile),tmpannfile);
		}  				
	    }       
	    
	}       
    }    
    return check;
}

int Hmm::get_sequence_and_tube_for_training(FILE* fIn_sequence,
					    const char* cIn_annotation,
					    const char* tube_file,
					    const int radius,
					    model_parameters* const MP,
					    Sequence * const sX, 
					    Sequence * const sY,
					    Tube<int> * const tTube,
					    vector<int>* Start_point,
					    vector<int>* End_point )
{
    int check = 0;
    
    // get sequence
 
    if(MP->is_PairHMM())
    {
	check+= get_sequence_pair(fIn_sequence, sX, sY, MP);
	if(check){
	    cout<<"Error: Parihmm class:: get_sequence_and_tube_for_training : "
		<<"get_sequence_pair failed. Abort.\n"<<flush;
	}
    }else{
	check+= get_single_sequence(fIn_sequence,sX,MP);
	if(check){
	    cout<<"Error: Hmm class:: get_sequence_and_tube_for_training : "
		<<"get_single_sequence failed. Abort.\n"<<flush;
	}
    }

    // if special emission is on and annotaiton file is provied,
    // set annotation to seq

    if((MP->is_SpecialEmit())&&(strcmp(cIn_annotation," ")))
    {	
	cout<<"set_annotation"<<endl;
	check += sX->set_annotation(cIn_annotation,MP);
	cout<<"out"<<endl;
	
	if (check != 0) 
	{
	    cout <<"ERROR: Hmm class:: get_sequence_and_tube_for_training "
		 <<"set_annotation for sequence "<<sX->get_ac() 
		 <<" and annotation file " 
		 << cIn_annotation << ".\n" << flush;
	}
	
	check += set_scores_for_annotated_sequence(sX,
						   0, // (0 = x, 1 = y)
						   //score_or_prob, // (0= score, 1 = prob)
						   this);						 	    	    
	if (check != 0) {
	    cout <<"ERROR: Hmm class:: get_sequence_and_tube_for_training "
		 <<"set_scores_for_annotated_sequence for "<<sX->get_ac() 
		 <<" and annotation file " 
		 << cIn_annotation << ".\n" << flush;
	}
	
	if(MP->is_PairHMM())
	{	    
	    check += sY->set_annotation(cIn_annotation, MP);
	    
	    if (check != 0) 
	    {
		cout <<"ERROR: Hmm class:: get_sequence_and_tube_for_training "
		     <<"set_annotation for sequence "<<sY->get_ac() 
		     <<" and annotation file " 
		     << cIn_annotation << ".\n" << flush;
	    }
	    
	    check += set_scores_for_annotated_sequence(sY,
						       1, // (0 = x, 1 = y)
						       //score_or_prob, // (0= score, 1 = prob)
						       this);						 
	    
	    if (check != 0) 
	    {
		cout <<"ERROR: Hmm class:: get_sequence_and_tube_for_training "
		     <<"set_scores_for_annotated_sequence for "<<sY->get_ac() 
		     <<" and annotation file " 
		     << cIn_annotation << ".\n" << flush;
	    }
	}	
    }
       
    // set Start point and End point
    (*Start_point)[0] = 0; 
    (*Start_point)[1] = 0;
    (*Start_point)[2] = 0;
    (*End_point)[0]   = sX->length();
    (*End_point)[1]   = sY->length();
    (*End_point)[2]   = this->get_number_of_states()-1;
    
    // set tube
    if((strcmp(tube_file," "))&&(!check))
    {
	// read tube from file 

        if(!cmp_nocase(tube_file,"Blastn")) // get tube from blastn
	{
	    Tube<int> tBlast_tube;
	    int iCount_matches, iCount_matches_new, iI;
	    unsigned long long int dTube_volume, dViterbi_volume;
	    Match* mMatches = NULL;
	    Match* mMatches_new = NULL;
	    check += get_blastn_results(sX,sY, &iCount_matches, &mMatches, MP);
	    if(check)
	    {
	        cout<<"ERROR: get_blastn_results for training: abort.\n"<<flush;     
	    }
	    check+= find_maximal_subset_of_matches(iCount_matches, mMatches, &iCount_matches_new, &mMatches_new);
	
	    if(check){
		cout<<"ERROR: find_maximal_subset_of_matches for training: abort.\n"<<flush;		
	    }	    

	    //mMatches_new->print(cout);
	      
	    tBlast_tube = convert_matches_to_tube(sX->length(),sY->length(),iCount_matches_new,mMatches_new, radius);

	    if(mMatches_new) delete [] mMatches_new;
	    mMatches_new=NULL;

	    if(mMatches) delete []mMatches;
	    mMatches = NULL;
	    
	    (*tTube) = tBlast_tube;

	}else if(!cmp_nocase(tube_file,"TBlastx")) // get tube from tblastx
	{
	    Tube<int> tBlast_tube;
	    int iCount_matches, iCount_matches_new, iI;
	    unsigned long long int dTube_volume, dViterbi_volume;
	    Match* mMatches = NULL;
	    Match* mMatches_new = NULL;
	    check += get_tblastx_results(sX,sY,&iCount_matches,&mMatches,MP);
	    if(check)
	    {
	      cout<<"ERROR: get_tblastx_results for training: abort.\n"<<flush;
	    }   
	    check+= find_maximal_subset_of_matches(iCount_matches, mMatches, &iCount_matches_new, &mMatches_new);

	    if(check){
		cout<<"ERROR: find_maximal_subset_of_matches for training: abort.\n"<<flush;
	    }	    

	    //mMatches_new->print(cout);
	      
	    tBlast_tube = convert_matches_to_tube(sX->length(),sY->length(),iCount_matches_new,mMatches_new, radius);

	    if(mMatches_new) delete [] mMatches_new;
	    mMatches_new=NULL;

	    if(mMatches) delete []mMatches;
	    mMatches = NULL;

	    (*tTube) = tBlast_tube;

	}else{

	    check += get_tube_from_file(tube_file,
				        sX->get_ac(),
				        sY->get_ac(),
				        sX->length(),
				        sY->length(),
				        tTube);
	    if(check)
	    {
	        cout <<"ERROR: Hmm class:: get_sequence_and_tube_for_training: "
		     <<"error in get_tube_from_file, file("<<tube_file<<").\n" << flush;	  
	    }		
	
        }
	
    }
       
    return check;
}

int Hmm::Viterbi_train_transition(const Sequence &sX, 
				  const Sequence &sY,
				  const int number_of_transitions_being_train,
				  int** num_of_from_list,
				  int*** from_list,
				  long** final_count,
				  Tube<int> tube,
				  const vector<int> Start_point,
				  const vector<int> End_point)								 
{
    int check = 0;
    
    if (tube.Empty()) 
    {	
	int max_y = sY.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(sX.length(), 2);
	
	for (int i = 0; i < sX.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }

    if (tube.Lx() != sX.length()) 
    {
	
	cout << "ERROR: Hmm::Viterbi_train_transition: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << sX.length() << ").\n";
	 check++;
    }
    // check that Start_x < End_x 
    
    if (Start_point[0] > End_point[0]) 
    {
	
	cout << "ERROR: Hmm::Viterbi_train_transition: Start_x (" << Start_point[0] << ") > End_x (" 
	     << End_point[0] << ").\n";
	check++;
    }
    
    if (check == 0) 
    {
	 
	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 

	// constants

	const int seq_start_x = 0;                 
	const int seq_end_x   = sX.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = sY.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = Start_point[0];
	const int Start_y     = Start_point[1];
	const int End_x       = End_point[0];
	const int End_y       = End_point[1];
	const int Start_state = Start_point[2];
	const int End_state   = End_point[2];
	const int shift_index = bitshift(alphabet);
	
	// variables 
	
	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   cur_offset, old_offset, temp_offset;
	int   new_state, old_state, deltax, deltay;
	int   max_state;
	int   end_y, old_start_y, old_end_y;
	Score emission_score, transition_score, max_score, test_score, new_score, old_score;	    
	
	long int cur_index, old_index, prev_index;
	
	vector<int> cur_index_vec(2);   // y position, state
	vector<int> old_index_vec(2);   // x position, y position, state
	
	//double final_score = Logzero;	
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) 
	{
	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    
	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }
	    
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) 
	    {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));

	    }
	}
	
	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	for (i=(Start_x-1); i <= (End_x-1); i++) 
	{
	    
	    if (i < 0)
	    {
		if (tube.GetElement(0, 0) == 0) 
		{
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		}
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }
	}

	// initialization
	for(i=0; i<number_of_transitions_being_train; i++)
	{
	    (*final_count)[i] = 0;
	}
		
	const long int tube_length = states*(seq_end_y+2);			
	
	// calculate max_deltax value 
	
	int max_d = 0;
	
	for (i = 1; i < (states-1); i++) 
	{
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;
	const int tube_width = max_deltax+1;

	// allocate memory for viterbi_tube
	
	long ***training_tube = NULL;	
	training_tube = new long**[number_of_transitions_being_train];	
	for(i=0; i<number_of_transitions_being_train; i++)
	{
	    training_tube[i] = new long*[tube_width];	    
	    for(j=0; j<tube_width; j++)
	    {
		training_tube[i][j] = new long[tube_length];
	    }
	}

	double **viterbi_tube = NULL;
	viterbi_tube = new double*[tube_width];	
	for(i=0; i<tube_width; i++)
	{	    
	    viterbi_tube[i]  = new double[tube_length];
	}
	
	int from = -1;
	int to = -1;

	// ----------------------------------------------------------------------
	// start : calculation of training_tube
	// ----------------------------------------------------------------------
	
	last_state = states-2;
	
	int count = 0;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) 
	{
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }
	    
	    if(xsteps==0)
	    {
		cur_offset = 0;
	    }else{
		cur_offset = (cur_offset+1)%tube_width;
	    }
	    
	    // initialize the current table

	    for(i=0; i<number_of_transitions_being_train;i++)
	    {
		for(j=0; j<tube_length; j++)
		{		  
		    training_tube[i][cur_offset][j] = 0;
		}
	    }
	    for(i=0; i<tube_length; i++)
	    {
		viterbi_tube[cur_offset][i] = Logzero;	       
	    }
	    
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) 
	    {
		
		cur_index_vec[0] = ysteps;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0)) 
		{
		    last_state = states-1;
		}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) 
		{ // initialise start of state path
		    
		    cur_index_vec[1] = Start_state;
		    cur_index  = (cur_index_vec[0] - Start_y) * states + cur_index_vec[1]; 
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state
		    
		    if ((deltax + deltay) > 0) 
		    {
			
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sX.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) 
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sY.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 
			
			emission_score = this->model[Start_state].get_emission_score(&sX, xsteps - deltax, &sY, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    viterbi_tube[cur_offset][cur_index] = emission_score; 		    
		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			max_score = Logzero;
			emission_score = Logzero; 
			
			cur_index_vec[1]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) 
			{                     			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_offset = cur_offset-(xsteps-old_xsteps);
				if(old_offset<0)
				{
				    old_offset+=tube_width;
				}
				
				old_index_vec[0] = old_ysteps;
				
				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | sX.letter(old_xsteps + i);
				    }
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | sY.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;}
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;}
				    
				    emission_score = this->model[new_state].get_emission_score(&sX, old_xsteps, &sY, old_ysteps, 
											       linear_index_emission);
				    
				}
				else // silent state
				{emission_score = 0.0;}
				
				// loop over old states, if emission_score in new_state is larger than Logzero
				cur_index  = (cur_index_vec[0] - Start_y) * states + cur_index_vec[1]; 							    				
				
				if (emission_score > Logzero) 
				{ 				 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();

				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) 
				    {
										
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index_vec[1]  = old_state;
					old_index = (old_index_vec[0] - Start_y) * states + old_index_vec[1];

					if(new_state==states-1)
					{
					    old_offset=cur_offset;
					}
					
					transition_score = this->new_get_transition_score(this, old_state, new_state, &sX, old_xsteps, &sY, old_ysteps);

					new_score = transition_score + emission_score;
					
					old_score = viterbi_tube[old_offset][old_index];
			
					if((transition_score > Logzero)&&(old_score > Logzero))
					{				
					    // find the maximum
					    test_score = old_score+new_score;
					    if(test_score > max_score)
					    {
						max_score = test_score;
						max_state = old_state;
					    }

					}				
				    } // loop over old states	    				   
				    
				    if(max_score > Logzero)
				    {
					viterbi_tube[cur_offset][cur_index] = max_score;
					old_index = (old_index_vec[0]-Start_y)*states + max_state;

					for(i=0; i<number_of_transitions_being_train; i++)
					{
					    if(num_of_from_list[i][new_state]>0)
					    {
						for(j=0; j<num_of_from_list[i][new_state];j++)
						{
						    from = from_list[i][new_state][j];
						    to = new_state;						  

						    if(from==max_state)
						    {							
							training_tube[i][cur_offset][cur_index] = training_tube[i][old_offset][old_index]+1;
						    }else{
							training_tube[i][cur_offset][cur_index] = training_tube[i][old_offset][old_index];
						    }
						    						  
						}
					    }else{
						training_tube[i][cur_offset][cur_index] = training_tube[i][old_offset][old_index];
					    }
					}
				    }// if max_score > Logzero
				    else{
					viterbi_tube[cur_offset][cur_index] = Logzero;
					for(i=0; i<number_of_transitions_being_train; i++)
					{
					    training_tube[i][cur_offset][cur_index] = 0;
					}
				    }			
				} // if emission_score in new_state > Logzero		      
				else{
				    viterbi_tube[cur_offset][cur_index] = Logzero;
				    for(i=0; i<number_of_transitions_being_train; i++)
				    {
					training_tube[i][cur_offset][cur_index] = 0;
				    }
				}
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    count++;
	} // loop over xsteps
	
	// ----------------------------------------------------------------------
	// end: calculation of training_tube
	// ----------------------------------------------------------------------
	
	// ----------------------------------------------------------------------
	// start : Obtained the result of the training
	//         set TTP[index].score = trainint_tube_score - forward_tube_score (log prob)
	// ----------------------------------------------------------------------
	
	cur_index_vec[0] = End_y;
	cur_index_vec[1] = End_state;    
       	cur_index = cur_index_vec[0]*states + cur_index_vec[1];
	
	if(viterbi_tube[cur_offset][cur_index]<Logzero)
	{
	    cout <<"ERROR: Hmm::Viteri_train_transition: the viterbi probability of sequence/sequence pair "
		 <<"is 0."<<endl;
	    check++;
	}
	
	if(!check) // set score
	{	    
	    for(int k =0; k<number_of_transitions_being_train; k++)
	    {
		if(training_tube[k][cur_offset][cur_index]<=0)
		{
		    cout <<"Warning: Hmm::Vitrerbi_train_transition: the count of transition ";
		    for(i=0; i<states; i++)
		    {		  
			for(j=0; j<num_of_from_list[k][i]; j++)
			{
			    cout<<from_list[k][i][j]<<"->"<<i<<" ";			 
			}
		    }
		    cout<<"is 0."<<endl;
		    (*final_count)[k] = 0;
		}else{
		    (*final_count)[k] = training_tube[k][cur_offset][cur_index];
		}	   
	    }
	}
	
	// free memory
	
	if (viterbi_tube) 
	{
	    for (i = 0; i < tube_width; i++) 
	    {
		if (viterbi_tube[i]) delete [] viterbi_tube[i];
		viterbi_tube[i] = NULL;
	    }
	    delete [] viterbi_tube; 
	}
	viterbi_tube = NULL;
	
	if (training_tube) 
	{
	    for(i=0; i< number_of_transitions_being_train; i++)
	    {
		if(training_tube[i])
		{
		    for (j = 0; j < tube_width; j++) 
		    {
			if (training_tube[i][j]) delete [] training_tube[i][j];
			training_tube[i][j] = NULL;
		    }
		    delete [] training_tube[i];
		}
		training_tube[i] = NULL;		
	    }
	    delete [] training_tube; 
	}
	training_tube = NULL;  

	
    }
    
    return check;
}

int Hmm::Posterior_train_transition(const Sequence &sX, 
				    const Sequence &sY,
				    const int number_of_transitions_being_train,
				    const int number_of_sample_paths,
				    int** num_of_from_list,
				    int*** from_list,
				    long** final_count,
				    Tube<int> tube,
				    const vector<int> Start_point,
				    const vector<int> End_point)					 
{
    int check = 0;
    
    if (tube.Empty()) 
    {	
	int max_y = sY.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(sX.length(), 2);
	
	for (int i = 0; i < sX.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }

    if (tube.Lx() != sX.length()) 
    {
	
	cout << "ERROR: Hmm::Posterior_train_transition: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << sX.length() << ").\n";
	 check++;
    }
    // check that Start_x < End_x 
    
    if (Start_point[0] > End_point[0]) 
    {
	
	cout << "ERROR: Hmm::Posterior_train_transition: Start_x (" << Start_point[0] << ") > End_x (" 
	     << End_point[0] << ").\n";
	check++;
    }
    
    if (check == 0) 
    {
	 
	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 

	// constants

	const int seq_start_x = 0;                 
	const int seq_end_x   = sX.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = sY.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = Start_point[0];
	const int Start_y     = Start_point[1];
	const int End_x       = End_point[0];
	const int End_y       = End_point[1];
	const int Start_state = Start_point[2];
	const int End_state   = End_point[2];
	const int shift_index = bitshift(alphabet);
	
	// variables 
	
	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   cur_offset, old_offset, temp_offset;
	int   deltax, deltay;
	int   sample_state, new_state,old_state;
	int   end_y, old_start_y, old_end_y;	   
	Score f_log_sum_score, f_sum_score, emission_score, transition_score;
	    
	Score first_forward_score, first_forward_transition_score;
	bool  f_get_first_score;

	long int cur_index, old_index, temp_index, prev_index;
	
	vector<int> cur_index_vec(2);   // y position, state
	vector<int> old_index_vec(2);   // x position, y position, state
	vector<int> temp_index_vec(2);  // for getting the sample pointer
	
	//double final_score = Logzero;	
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) 
	{		

	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    
	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }
	    
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) 
	    {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));

	    }
	}
	
	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	for (i=(Start_x-1); i <= (End_x-1); i++) 
	{
	    
	    if (i < 0)
	    {
		if (tube.GetElement(0, 0) == 0) 
		{
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		}
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }
	}

	// initialization
	for(i=0; i<number_of_transitions_being_train; i++)
	{
	    (*final_count)[i] = 0;
	}
		
	const long int tube_length = states*(seq_end_y+2);			
	
	// calculate max_deltax value 
	
	int max_d = 0;
	
	for (i = 1; i < (states-1); i++) 
	{
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;
	const int tube_width = max_deltax+1;

	// allocate memory for viterbi_tube
	
	long ****training_tube = NULL;	
	training_tube = new long***[number_of_sample_paths];
	for(i=0; i<number_of_sample_paths; i++)
	{
	    training_tube[i] = new long**[number_of_transitions_being_train];	
	    for(j=0; j<number_of_transitions_being_train; j++)
	    {
		training_tube[i][j] = new long*[tube_width];	    
		for(int k =0; k<tube_width; k++)
		{
		    training_tube[i][j][k] = new long[tube_length];		  
		}
	    }
	}

	double **forward_tube = NULL;
	forward_tube = new double*[tube_width];	
	for(i=0; i<tube_width; i++)
	{	    
	    forward_tube[i]  = new double[tube_length];
	}

	double* p = new double[tube_length];
	int* p_index = new int[tube_length];
	double* p_transition = new double[tube_length];
	int number_of_p_index = 0;
	
	int from = -1;
	int to = -1;

	// ----------------------------------------------------------------------
	// start : calculation of training_tube
	// ----------------------------------------------------------------------
	
	last_state = states-2;
	
	int count = 0;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) 
	{
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }
	    
	    if(xsteps==0)
	    {
		cur_offset = 0;
	    }else{
		// swap plate 0 and plate 1
		cur_offset = (cur_offset+1)%tube_width;
	    }
	    
	    // initialize the current table
	    
	    for(i=0; i<number_of_sample_paths; i++)
	    {
		for(j=0; j<number_of_transitions_being_train;j++)
		{
		    for(int k=0; k<tube_length; k++)
		    {		  
			training_tube[i][j][cur_offset][k] = 0;
		    }
		}
	    }

	    for(i=0; i<tube_length; i++)
	    {
		forward_tube[cur_offset][i] = Logzero;	       
	    }	    
	        
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) 
	    {
		cur_index_vec[0] = ysteps;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0))
		{
		    last_state = states-1;
		}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) 
		{ // initialise start of state path
		    
		    cur_index_vec[1] = Start_state;
		    cur_index  = (cur_index_vec[0] - Start_y) * states + cur_index_vec[1]; 
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state
		    
		    if ((deltax + deltay) > 0) 
		    {
			
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sX.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) 
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sY.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 
			
			emission_score = this->model[Start_state].get_emission_score(&sX, xsteps - deltax, &sY, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    forward_tube[cur_offset][cur_index] = emission_score; 		    
		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			number_of_p_index = 0;

			f_sum_score = Logzero;
			f_log_sum_score = Logzero;

			emission_score = Logzero;

			cur_index_vec[1]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) 
			{                     			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_offset = cur_offset-(xsteps-old_xsteps);
				if(old_offset<0)
				{
				    old_offset+=tube_width;
				}
				
				old_index_vec[0] = old_ysteps;
				
				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | sX.letter(old_xsteps + i);
				    }
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | sY.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;} 
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;} 
				    
				    emission_score = this->model[new_state].get_emission_score(&sX, old_xsteps, &sY, old_ysteps, 
											       linear_index_emission);
				    
				}
				else // silent state
				{emission_score = 0.0;}
				
				// loop over old states, if emission_score in new_state is larger than Logzero
				cur_index  = (cur_index_vec[0] - Start_y) * states + cur_index_vec[1]; 			

				if (emission_score > Logzero) 
				{ 				 				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();

				    f_get_first_score = false;

				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) 				
				    {
										
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);

					old_index_vec[1]  = old_state;
					old_index = (old_index_vec[0] - Start_y) * states + old_index_vec[1];
					if(new_state==states-1)
					{
					    old_offset=cur_offset;
					}
					
					transition_score = this->new_get_transition_score(this, old_state, new_state, &sX, old_xsteps, &sY, old_ysteps);
			
					if((transition_score > Logzero)&&(forward_tube[old_offset][old_index] > Logzero))
					{				
					    // calculate forward tube
					    p[old_index] = forward_tube[old_offset][old_index] + transition_score;
					    p_index[number_of_p_index] = old_index;
					    number_of_p_index++;       				

					    if(!f_get_first_score)
					    {					       
						f_sum_score = p[old_index];
						f_log_sum_score = 0;
						f_get_first_score = true;						
					    }else{
						f_log_sum_score = f_log_sum_score + exp(p[old_index] - f_sum_score);	
					    }
					}							
				    } // loop over old states			   		
				    
				    if (f_sum_score > Logzero)
				    {					
					f_log_sum_score = log(1+f_log_sum_score);
					
					f_sum_score = f_sum_score + f_log_sum_score+ emission_score;	   
					
					// assign sum to forward tube
					forward_tube[cur_offset][cur_index] = f_sum_score;	 				
				    } // if f_sum_score > Logzero					
				    // calculate the probability distribution vector p
				    				    
				    for(i=0; i<number_of_p_index; i++)
				    {					
					double arg = p[p_index[i]]+emission_score-forward_tube[cur_offset][cur_index];
					if(i==0)
					{
					    if(p[p_index[i]]>Logzero)
					    {
						p[p_index[i]] = exp(arg);
					    }else{
						p[p_index[i]] = 0;
					    }				
					}else{
					    if(p[p_index[i]]>Logzero)
					    {
						p[p_index[i]] = p[p_index[i-1]]+exp(arg);
					    }else{
						p[p_index[i]] = p[p_index[i-1]];
					    }
					    
					}			
				    }
				    
				    // choose number_of_sample_paths pointers and add to the
				    // corresponding entry of table t
				    for(i=0; i<number_of_sample_paths; i++)
				    {
					double rannum = static_cast<double>(rand())/
					    (static_cast<double>(RAND_MAX)+ static_cast<double>(1));

					for(j=0;j<number_of_p_index;j++)
					{
					    if(rannum<p[p_index[j]])
					    {
						sample_state = p_index[j]%states;
						break;
					    }
					}
				       			      
					old_index = (old_index_vec[0]-Start_y)*states + sample_state;
					for(j=0;j<number_of_transitions_being_train;j++)
					{
					    if(num_of_from_list[j][new_state]>0)
					    {
						for(int k=0; k<num_of_from_list[j][new_state];k++)
						{
						    from = from_list[j][new_state][k];
						    to = new_state;
						    if(from==sample_state)
						    {
							training_tube[i][j][cur_offset][cur_index] = training_tube[i][j][old_offset][old_index]+1;
						    }else{
							training_tube[i][j][cur_offset][cur_index] = training_tube[i][j][old_offset][old_index];
						    }
						}
					    }else{
						training_tube[i][j][cur_offset][cur_index] = training_tube[i][j][old_offset][old_index];
					    }
					}
				    } // finish sampling 			
				}
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    count++;
	} // loop over xsteps
	
	// ----------------------------------------------------------------------
	// end: calculation of training_tube
	// ----------------------------------------------------------------------
	
	// ----------------------------------------------------------------------
	// start : Obtained the result of the training
	//         set TTP[index].score = trainint_tube_score - forward_tube_score (log prob)
	// ----------------------------------------------------------------------
	
	cur_index_vec[0] = End_y;
	cur_index_vec[1] = End_state;    
	
	cur_index = cur_index_vec[0]*states + cur_index_vec[1];
	
	if(forward_tube[cur_offset][cur_index]<Logzero)
	{
	    cout <<"ERROR: Hmm::Posterior_train_transition: the viterbi probability of sequence/sequence pair "
		 <<"is 0."<<endl;
	    check++;
	}
	
	if(!check) // set score
	{	    
	    for(int k = 0; k<number_of_transitions_being_train; k++)
	    {
		(*final_count)[k] = 0;
		for(int l=0; l<number_of_sample_paths; l++)
		{
		    (*final_count)[k] += training_tube[l][k][cur_offset][cur_index];
		}
	   
		if((*final_count)[k]<=0)
		{
		    cout <<"Warning: Hmm:: Posterior_train_transition: the count of transition ";
		    for(i=0; i<states; i++)
		    {		  
			for(j=0; j<num_of_from_list[k][i]; j++)
			{
			    cout<<from_list[k][i][j]<<"->"<<i<<" ";			 
			}
		    }
		    cout<<"is 0."<<endl;
		}
	    }
	}
	
	// free memory
	
	if (forward_tube) 
	{
	    for (i = 0; i < tube_width; i++) 
	    {
		if (forward_tube[i]) delete [] forward_tube[i];
		forward_tube[i] = NULL;
	    }
	    delete [] forward_tube; 
	}
	forward_tube = NULL;
	
	if (training_tube) 
	{
	    for(i=0; i<number_of_sample_paths; i++)
	    {
		if(training_tube[i])
		{
		    for(j=0; j< number_of_transitions_being_train; j++)
		    {
			if(training_tube[i][j])
			{
			    for (int k  = 0; k < tube_width; k++) 
			    {
				if (training_tube[i][j][k]) delete [] training_tube[i][j][k];
				training_tube[i][j][k] = NULL;
			    }
			    delete [] training_tube[i][j];
			}
			training_tube[i][j] = NULL;		
		    }
		    delete [] training_tube[i];
		}
		training_tube[i] = NULL;
	    }
	    delete [] training_tube; 
	}
	training_tube = NULL;  

	if(p) delete [] p;
	p = NULL;

	if(p_index) delete [] p_index;
	p_index = NULL;
	
	if(p_transition) delete [] p_transition;
	p_transition = NULL;
	
    }
    
    return check;
}

int Hmm::BaumWelch_train_transition(const Sequence &sX, 
				    const Sequence &sY,
				    const int number_of_transitions_being_train,
				    int** num_of_from_list,
				    int*** from_list,
				    double** final_score,
				    double& forward_score,
				    Tube<int> tube,
				    const vector<int> Start_point,
				    const vector<int> End_point)								 
{
    int check = 0;
    
    if (tube.Empty()) 
    {	
	int max_y = sY.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(sX.length(), 2);
	
	for (int i = 0; i < sX.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }

    if (tube.Lx() != sX.length()) 
    {
	
	cout << "ERROR: Hmm::BaumWelch_train_transition: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << sX.length() << ").\n";
	 check++;
    }
    // check that Start_x < End_x 
    
    if (Start_point[0] > End_point[0]) 
    {
	
	cout << "ERROR: Hmm::BaumWelch_train_transition: Start_x (" << Start_point[0] << ") > End_x (" 
	     << End_point[0] << ").\n";
	check++;
    }
    
    

    if (check == 0) 
    {
	 

	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 

	// constants

	const int seq_start_x = 0;                 
	const int seq_end_x   = sX.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = sY.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = Start_point[0];
	const int Start_y     = Start_point[1];
	const int End_x       = End_point[0];
	const int End_y       = End_point[1];
	const int Start_state = Start_point[2];
	const int End_state   = End_point[2];
	const int shift_index = bitshift(alphabet);
	
	// variables 
	
	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   cur_offset, old_offset, temp_offset;
	int   new_state, old_state, deltax, deltay;
	int   end_y, old_start_y, old_end_y;
	Score f_log_sum_score, f_sum_score, emission_score, transition_score;
	Score* t_log_sum_score1 = new Score[number_of_transitions_being_train];
	Score* t_log_sum_score2 = new Score[number_of_transitions_being_train];
	Score* t_sum_score = new Score[number_of_transitions_being_train]; 
	    
	Score first_forward_score, first_forward_transition_score;
	bool  f_get_first_score; 
	bool* t_get_first_score = new bool[number_of_transitions_being_train];
	
	long int cur_index, old_index, prev_index;
	
	vector<int> cur_index_vec(2);   // y position, state
	vector<int> old_index_vec(2);   // x position, y position, state
	
	//double final_score = Logzero;	
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
		
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) 
	{
	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    
	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }
	    
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) 
	    {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));

	    }		
	}
	
	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	for (i=(Start_x-1); i <= (End_x-1); i++) 
	{
	    
	    if (i < 0)
	    {
		if (tube.GetElement(0, 0) == 0) 
		{
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		}
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }
	}

	// initialization
	
	for(i=0; i<number_of_transitions_being_train; i++)
	{
	    (*final_score)[i] = Logzero;
	}
		
	const long int tube_length = states*(seq_end_y+2);			
	
	// calculate max_deltax value 
	
	int max_d = 0;
	
	for (i = 1; i < (states-1); i++) 
	{
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;
	const int tube_width = max_deltax+1;

	// allocate memory for viterbi_tube
	
	double ***training_tube = NULL;	
	training_tube = new double**[number_of_transitions_being_train];	
	for(i=0; i<number_of_transitions_being_train; i++)
	{
	    training_tube[i] = new double*[tube_width];	    
	    for(j=0; j<tube_width; j++)
	    {
		training_tube[i][j] = new double[tube_length];
	    }
	}

	double **forward_tube = NULL;
	forward_tube = new double*[tube_width];	
	for(i=0; i<tube_width; i++)
	{	    
	    forward_tube[i]  = new double[tube_length];
	}
	
	int from = -1;
	int to = -1;

	// ----------------------------------------------------------------------
	// start : calculation of training_tube
	// ----------------------------------------------------------------------
	
	last_state = states-2;
	
	int count = 0;

	double min_exp_arg = 0;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) 
	{
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }
	    
	    if(xsteps==0)
	    {
		cur_offset = 0;
	    }else{
		// swap plate 0 and plate 1
		cur_offset = (cur_offset+1)%tube_width;		
	    }
	    
	    // initialize the current table

	    for(i=0; i<number_of_transitions_being_train;i++)
	    {
		for(j=0; j<tube_length; j++)
		{		  
		    training_tube[i][cur_offset][j] = Logzero;
		}
	    }
	    for(i=0; i<tube_length; i++)
	    {
		forward_tube[cur_offset][i] = Logzero;	       
	    }
	    
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) 
	    {
		
		cur_index_vec[0] = ysteps;
		long cur_temp_index = (cur_index_vec[0]-Start_y)*states;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0)) 
		{
		    last_state = states-1;
		}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) 
		{ // initialise start of state path
		    
		    cur_index_vec[1] = Start_state;
		    cur_index = cur_temp_index + cur_index_vec[1];
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state
		    
		    if ((deltax + deltay) > 0) 
		    {
			
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sX.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) 
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sY.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 
			
			emission_score = this->model[Start_state].get_emission_score(&sX, xsteps - deltax, &sY, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    forward_tube[cur_offset][cur_index] = emission_score; 		    
		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			f_sum_score      = Logzero;
			f_log_sum_score  = Logzero;
			
			for(i=0; i<number_of_transitions_being_train; i++)
			{
			    t_sum_score[i]      = Logzero;
			    t_log_sum_score1[i] = Logzero;
			    t_log_sum_score2[i] = Logzero;
			}
			
			emission_score = Logzero; 
			
			cur_index_vec[1]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) 
			{                     			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_offset = cur_offset-(xsteps-old_xsteps);
				if(old_offset<0)
				{
				    old_offset+=tube_width;
				}
				
				old_index_vec[0] = old_ysteps;
				
				long old_temp_index = (old_index_vec[0]-Start_y)*states;

				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | sX.letter(old_xsteps + i);
				    }
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | sY.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;} 
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;} 
				    
				    emission_score = this->model[new_state].get_emission_score(&sX, old_xsteps, &sY, old_ysteps, 
											       linear_index_emission);
				    
				}
				else // silent state
				{emission_score = 0.0;}
				
				// loop over old states, if emission_score in new_state is larger than Logzero
				
				cur_index = cur_temp_index + cur_index_vec[1];

				if (emission_score > Logzero) 
				{ 				 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    f_get_first_score = false;
				    for(i=0;i<number_of_transitions_being_train; i++)
				    {
					t_get_first_score[i] = false;
				    }
				    
				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) 
				    {
										
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);

					old_index_vec[1]  = old_state;
					old_index = old_temp_index + old_index_vec[1];

					if(new_state==states-1)
					{
					    old_offset=cur_offset;
					}
					
					transition_score = this->new_get_transition_score(this, old_state, new_state, &sX, old_xsteps, &sY, old_ysteps);
			
					if((transition_score > Logzero)&&(forward_tube[old_offset][old_index] > Logzero))
					{
					    // calculate forward tube
					    if(!f_get_first_score)
					    {					       
						first_forward_score = forward_tube[old_offset][old_index];	 				
						first_forward_transition_score = transition_score;				
						f_sum_score = first_forward_score + first_forward_transition_score;
						
						f_log_sum_score = 0;
						f_get_first_score = true;						
					    }else{
						f_log_sum_score = f_log_sum_score + exp(forward_tube[old_offset][old_index] + transition_score - f_sum_score);	
					    }
					    
					}
				
					// calculate the training tube

					if(transition_score > Logzero)
					{
					    for(i=0; i<number_of_transitions_being_train;i++)
					    {
						if(training_tube[i][old_offset][old_index] > Logzero)
						{
						    // calculate training tube
						    if(!t_get_first_score[i])
						    {
							//first_training_score[i] = training_tube[old_offset][old_index];
							//first_training_transition_score[i] = transition_score;
							t_sum_score[i] = training_tube[i][old_offset][old_index] + transition_score;
							t_log_sum_score1[i] = 0;
							t_get_first_score[i] = true;
						    }else{
							t_log_sum_score1[i] = t_log_sum_score1[i] + exp(training_tube[i][old_offset][old_index] + transition_score - t_sum_score[i]);
							
						    }					
						}				
					    }
					}						       	  
					
				    } // loop over old states
				        
				    if (f_sum_score > Logzero)
				    {
					
					f_log_sum_score = log(1+f_log_sum_score);
					
					f_sum_score = f_sum_score + f_log_sum_score+ emission_score;			
					
					// assign sum to forward tube
					forward_tube[cur_offset][cur_index] = f_sum_score;								
					
				    } // if f_sum_score > Logzero					    
				    
				    for(i=0; i<number_of_transitions_being_train; i++)
				    {
					if(num_of_from_list[i][new_state]>0) // need to add forward values
					{					
					    if(t_sum_score[i] <= Logzero) // i.e. the first column
					    {
						
						from = from_list[i][new_state][0];
						to   = new_state;					    
						transition_score =  this->new_get_transition_score(this, from, to, &sX, old_xsteps, &sY, old_ysteps);
						
						old_index = old_temp_index + from;
						if(forward_tube[old_offset][old_index] > Logzero)
						{
							t_log_sum_score1[i] = forward_tube[old_offset][old_index] + transition_score;
						}

						j = 1;					     
					    }else{
												
						t_log_sum_score1[i] = log(1+t_log_sum_score1[i]); 
						t_log_sum_score1[i] = t_sum_score[i] + t_log_sum_score1[i];
						// t_log_sum_score1[i] is the sum of the weighted probabilities for the i-th parameter
						// before adding the forward values of the current seq. position
						j = 0;

					    }
					    
					    t_log_sum_score2[i] = 0;
					    for(int k=j ; k<num_of_from_list[i][new_state]; k++)
					    {
						from = from_list[i][new_state][k];
						to   = new_state;
						transition_score =  this->new_get_transition_score(this, from, to, &sX, old_xsteps, &sY, old_ysteps);
						old_index = old_temp_index + from;
						if(forward_tube[old_offset][old_index]>Logzero)
						{
							if(t_log_sum_score1[i]>Logzero)
							{
								t_log_sum_score2[i] = t_log_sum_score2[i] + 
									exp(forward_tube[old_offset][old_index] + transition_score - t_log_sum_score1[i]);
							}else{		
								t_log_sum_score1[i] = forward_tube[old_offset][old_index] + transition_score;
							}
						}
					    }
					    
					    
					    t_log_sum_score2[i] = log(1+t_log_sum_score2[i]);
					    t_sum_score[i] = t_log_sum_score1[i] + t_log_sum_score2[i] + emission_score;      
					    training_tube[i][cur_offset][cur_index] = t_sum_score[i];	
					    
					}else{ // do not need to add the portion from the forward values
					    if(t_sum_score[i] > Logzero)
					    {
						t_log_sum_score1[i] = log(1+t_log_sum_score1[i]);
						t_sum_score[i] = t_sum_score[i] + t_log_sum_score1[i] + emission_score;
					    }
					    training_tube[i][cur_offset][cur_index] = t_sum_score[i];	    	
					}			        	    	
				    }
				} // if emission_score in new_state > Logzero
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states

		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    count++;
	} // loop over xsteps
	
	// ----------------------------------------------------------------------
	// end: calculation of training_tube
	// ----------------------------------------------------------------------
	
	// ----------------------------------------------------------------------
	// start : Obtained the result of the training
	//         set TTP[index].score = trainint_tube_score - forward_tube_score (log prob)
	// ----------------------------------------------------------------------
	
	cur_index_vec[0] = End_y;
	cur_index_vec[1] = End_state;    
	
	cur_index = cur_index_vec[0]*states + cur_index_vec[1];
	
	if(forward_tube[cur_offset][cur_index]<Logzero)
	{
	    cout <<"ERROR: Hmm::train_transition: the forward probability of sequence/sequence pair "
		 <<"is 0."<<endl;
	    check++;
	}
	
	if(!check) // set score
	{
	    forward_score = forward_tube[cur_offset][cur_index];
	    for(int k =0; k<number_of_transitions_being_train; k++)
	    {
		if(training_tube[k][cur_offset][cur_index]<=Logzero)
		{
		    cout <<"Warning: Hmm::BaumWelch_train_transition: the  probability of transition ";
		    for(i=0; i<states; i++)
		    {		  
			for(j=0; j<num_of_from_list[k][i]; j++)
			{
			    cout<<from_list[k][i][j]<<"->"<<i<<" ";			 
			}
		    }
		    cout<<"is 0."<<endl;
		    (*final_score)[k] = Logzero;
		}else{
		    (*final_score)[k] = training_tube[k][cur_offset][cur_index] - forward_score;
		}	   
	    }
	}
	
	// free memory
	
	if (forward_tube) 
	{
	    for (i = 0; i < tube_width; i++) 
	    {
		if (forward_tube[i]) delete [] forward_tube[i];
		forward_tube[i] = NULL;
	    }
	    delete [] forward_tube; 
	}
	forward_tube = NULL;
	
	if (training_tube) 
	{
	    for(i=0; i< number_of_transitions_being_train; i++)
	    {
		if(training_tube[i])
		{
		    for (j = 0; j < tube_width; j++) 
		    {
			if (training_tube[i][j]) delete [] training_tube[i][j];
			training_tube[i][j] = NULL;
		    }
		    delete [] training_tube[i];
		}
		training_tube[i] = NULL;		
	    }
	    delete [] training_tube; 
	}
	training_tube = NULL;   

	if(t_get_first_score) delete [] t_get_first_score;
	t_get_first_score = NULL;

	if(t_sum_score) delete [] t_sum_score;
	t_sum_score = NULL;

	if(t_log_sum_score1) delete [] t_log_sum_score1;
	t_log_sum_score1 = NULL;

	if(t_log_sum_score2) delete [] t_log_sum_score2;
	t_log_sum_score2 = NULL;
	
    }
    
    return check;
}

int Hmm::Viterbi_train_emission(const Sequence &sX, 
				const Sequence &sY,
				const long number_of_emissions_being_train,
				long* emission_index,
				int** from_list,		      
				long** final_count,
				Tube<int> tube,
				const vector<int> Start_point,
				const vector<int> End_point) 								 
{
    int check = 0;

    if (tube.Empty()) 
    {	
	int max_y = sY.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(sX.length(), 2);
	
	for (int i = 0; i < sX.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }

    if (tube.Lx() != sX.length()) 
    {
	
	cout << "ERROR: Hmm::Viterbi_train_emission: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << sX.length() << ").\n";
	 check++;
    }
    // check that Start_x < End_x 
    
    if (Start_point[0] > End_point[0]) 
    {
	
	cout << "ERROR: Hmm::Viterbi_train_emission: Start_x (" << Start_point[0] << ") > End_x (" 
	     << End_point[0] << ").\n";
	check++;
    }
    
    if (check == 0) 
    {
	 
	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 

	// constants

	const int seq_start_x = 0;                 
	const int seq_end_x   = sX.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = sY.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = Start_point[0];
	const int Start_y     = Start_point[1];
	const int End_x       = End_point[0];
	const int End_y       = End_point[1];
	const int Start_state = Start_point[2];
	const int End_state   = End_point[2];
	const int shift_index = bitshift(alphabet);
	
	// variables 
	
	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   cur_offset, old_offset, temp_offset;
	int   new_state, old_state, deltax, deltay;
	int   max_state, max_linear_index_emission;
	int   end_y, old_start_y, old_end_y;
	Score new_score, test_score, max_score, old_score, emission_score, transition_score;

	long int cur_index, old_index;
	
	vector<int> cur_index_vec(2);   // y position, state
	vector<int> old_index_vec(2);   // x position, y position, state       	
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) 
	{
	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    
	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }
	    
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) 
	    {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));
		
	    }
	}
	
	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	for (i=(Start_x-1); i <= (End_x-1); i++) {
	    
	    if (i < 0)
	    {
		if (tube.GetElement(0, 0) == 0) 
		{
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		 }
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }
	}
	
	// calculate max_deltax value 

	for(i=0; i<number_of_emissions_being_train; i++)
	{
	    (*final_count)[i] = 0;
	}

	const long int tube_length = states*(seq_end_y+2);
	
	// calculate max_deltax value

	int max_d = 0;
	
	for (i = 1; i < (states-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;
	const int tube_width = max_deltax+1;	
		 
	// allocate memory for viterbi_tube
	
	long ***training_tube = NULL;
	training_tube = new long**[number_of_emissions_being_train];
	for(i=0; i<number_of_emissions_being_train; i++)
	{
	    training_tube[i] = new long*[tube_width];
	    for(j=0; j<tube_width; j++)	 
	    {
		training_tube[i][j] = new long[tube_length];
	    }
	}
	

	double **viterbi_tube = NULL;
	viterbi_tube = new double*[tube_width];	
	for(i=0; i<tube_width; i++)
	{	   
	    viterbi_tube[i]  = new double[tube_length];
	}

	int** linear_emission_tube = NULL;
	linear_emission_tube = new int*[tube_width];
	for(i=0; i<tube_width;i++)
	{
	    linear_emission_tube[i] = new int[tube_length];
	}
	
	int from = -1;
	int to = -1;
	
	// ----------------------------------------------------------------------
	// start : calculation of training_tube
	// ----------------------------------------------------------------------
	
	int count = 0;
	
	last_state = states-2;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) 
	{
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }
	    
	    if(xsteps==0)
	    {
		cur_offset = 0;
	    }else{
		// swap plate 0 and plate 1
		cur_offset = (cur_offset+1)%tube_width;
	    }
	    
	    // initialize the current plate 
	    for(i=0; i<number_of_emissions_being_train;i++)
	    {
		for(j=0; j<tube_length; j++)
		{		  
		    training_tube[i][cur_offset][j] = 0;
		}
	    }
	    
	    for(i=0; i<tube_length; i++)
	    {
		viterbi_tube[cur_offset][i] = Logzero;
		linear_emission_tube[cur_offset][i] = -1;
	    }
	    
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) 
	    {
		
		cur_index_vec[0] = ysteps;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0)) 
		{
		    last_state = states-1;
		}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) 
		{ // initialise start of state path
		    
		    cur_index_vec[1] = Start_state;
		    cur_index  = (cur_index_vec[0] - start_y) * states + cur_index_vec[1]; 
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state
		    linear_index_emission = -1;
		    
		    if ((deltax + deltay) > 0) 
		    {
			
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sX.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) 
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sY.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 

			emission_score = this->model[Start_state].get_emission_score(&sX, xsteps - deltax, &sY, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    viterbi_tube[cur_offset][cur_index] = emission_score; 		    
		    linear_emission_tube[cur_offset][cur_index] = linear_index_emission;
		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			max_score = Logzero;
			emission_score = Logzero; 
			
			cur_index_vec[1]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) 
			{                     			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_offset = cur_offset-(xsteps-old_xsteps);
				if(old_offset<0)
				{
				    old_offset += tube_width;
				}
				old_index_vec[0] = old_ysteps;
				
				// get emission score for new state
				
				linear_index_emission = -1;
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | sX.letter(old_xsteps + i);
				    }
				    
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | sY.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;} 
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;} 

				    emission_score = this->model[new_state].get_emission_score(&sX, old_xsteps, &sY, old_ysteps, 
											       linear_index_emission);
				    
				}
				else // silent state
				{emission_score = 0.0;}
				
				cur_index  = (cur_index_vec[0] - start_y) * states + cur_index_vec[1]; 

				linear_emission_tube[cur_offset][cur_index] =linear_index_emission;

				 // loop over old states, if emission_score in new_state is larger than Logzero

				if (emission_score > Logzero) 
				{ 				 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) 
				    {
					
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index_vec[1]  = old_state;
					old_index = (old_index_vec[0] - old_start_y) * states + old_index_vec[1];
					
					if(new_state==states-1)
					{
					    old_offset=cur_offset;
					}

					transition_score = this->new_get_transition_score(this, old_state, new_state, &sX, old_xsteps, &sY, old_ysteps);
					
					new_score = transition_score + emission_score;
					
					old_score = viterbi_tube[old_offset][old_index];
					
					if((transition_score > Logzero)&&(old_score > Logzero))
					{
					     // find the maximum
					    test_score = old_score+new_score;
					    if(test_score > max_score)
					    {
						max_score = test_score;
						max_state = old_state;
						max_linear_index_emission = linear_emission_tube[old_offset][old_index];
					    }
					    
					}
					
				    } // loop over old states
				    
				    if(max_score>Logzero)
				    {
					viterbi_tube[cur_offset][cur_index] = max_score;
					
					old_index = (old_index_vec[0]-start_y)*states + max_state;
					
					for(i=0; i<number_of_emissions_being_train; i++)
					{
					    if(from_list[i][max_state]>0)
					    {						 
						if(emission_index[i]==max_linear_index_emission)
						{						     
						    training_tube[i][cur_offset][cur_index] = training_tube[i][old_offset][old_index]+1;	         
						}else{
						    training_tube[i][cur_offset][cur_index] = training_tube[i][old_offset][old_index];
						}
					    }else{
						training_tube[i][cur_offset][cur_index] = training_tube[i][old_offset][old_index];
					    }						
					}			
				    } // if max_score > Logzero				    
				    else
				    {
					viterbi_tube[cur_offset][cur_index] = Logzero;
					for(i=0; i<number_of_emissions_being_train; i++)
					{
					    training_tube[i][cur_offset][cur_index] = 0;
					}
				    }
				    
				} // if emission_score in new_state > Logzero		      				
				else
				{
				    viterbi_tube[cur_offset][cur_index] = Logzero;
				    for(i=0; i<number_of_emissions_being_train; i++)
				    {
					training_tube[i][cur_offset][cur_index] = 0;
				    }
				}				
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    count++;
	} // loop over xsteps

	// ----------------------------------------------------------------------
	// end: calculation of training_tube
	// ----------------------------------------------------------------------
	
	// ----------------------------------------------------------------------
	// start : Obtained the result of the training
	//       
	// ----------------------------------------------------------------------
	
	cur_index_vec[0] = End_y;
	cur_index_vec[1] = End_state;    
	
	cur_index = cur_index_vec[0]*states + cur_index_vec[1];
	
	if(viterbi_tube[cur_offset][cur_index]<Logzero)
	{
	    cout <<"ERROR: Hmm::Viterbi_train_emission: the forward probability of sequence/sequence pair "
		 <<"is 0."<<endl;
	    check++;
	}
	
	
	if(!check) // set score
	{
	    for(int k=0; k<number_of_emissions_being_train; k++)
	    {
		if(training_tube[k][cur_offset][cur_index]<=0)
		{
		    cout <<"Warning: Hmm::Viterbi_train_emission: the probability of emission ";
		    for(i=0; i<states; i++)
		    {
			if(from_list[k][i]!=-1)
			{			 
			    cout<<" state : "<<from_list[k][i]<<" emit : "<<emission_index[k];			 
			}
		    }
		    cout<<" is 0."<<endl;
		    (*final_count)[k] = 0;
		}else{
		    (*final_count)[k] = training_tube[k][cur_offset][cur_index];
		}	   
	    }
	}
	
	// free memory
	
	if (viterbi_tube) {
	    for (i = 0; i < tube_width; i++) {
		if (viterbi_tube[i]) delete viterbi_tube[i];
		viterbi_tube[i] = NULL;
	    }
	    delete [] viterbi_tube; 
	}
	viterbi_tube = NULL;
	
	if (linear_emission_tube) {
	    for (i = 0; i < tube_width; i++) {
		if (linear_emission_tube[i]) delete linear_emission_tube[i];
		linear_emission_tube[i] = NULL;
	    }
	    delete [] linear_emission_tube; 
	}
	linear_emission_tube = NULL;

	if (training_tube) 
	{
	    for(i=0; i< number_of_emissions_being_train; i++)
	    {
		if(training_tube[i])
		{
		    for (j = 0; j < tube_width; j++) 
		    {
			if (training_tube[i][j]) delete [] training_tube[i][j];
			training_tube[i][j] = NULL;
		    }
		    delete [] training_tube[i];
		 }
		training_tube[i] = NULL;		
	    }
	    delete [] training_tube; 
	}
	training_tube = NULL;   
		
    }
    
    return check;
}

int Hmm::Posterior_train_emission(const Sequence &sX, 
				  const Sequence &sY,
				  const long number_of_emissions_being_train,
				  const int number_of_sample_paths,
				  long* emission_index,				       
				  int** from_list,		      
				  long** final_count,
				  Tube<int> tube,
				  const vector<int> Start_point,
				  const vector<int> End_point) 								 
{
    int check = 0;

    if (tube.Empty()) 
    {	
	int max_y = sY.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(sX.length(), 2);
	
	for (int i = 0; i < sX.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }

    if (tube.Lx() != sX.length()) 
    {
	
	cout << "ERROR: Hmm::Posterior_train_emission: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << sX.length() << ").\n";
	 check++;
    }
    // check that Start_x < End_x 
    
    if (Start_point[0] > End_point[0]) 
    {
	
	cout << "ERROR: Hmm::Posterior_train_emission: Start_x (" << Start_point[0] << ") > End_x (" 
	     << End_point[0] << ").\n";
	check++;
    }
    
    if (check == 0) 
    {
	 
	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 

	// constants

	const int seq_start_x = 0;                 
	const int seq_end_x   = sX.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = sY.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = Start_point[0];
	const int Start_y     = Start_point[1];
	const int End_x       = End_point[0];
	const int End_y       = End_point[1];
	const int Start_state = Start_point[2];
	const int End_state   = End_point[2];
	const int shift_index = bitshift(alphabet);
	
	// variables 
	
	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   cur_offset, old_offset, temp_offset;
	int   sample_state;
	int   new_state, old_state, deltax, deltay;
	int   sample_linear_index_emission;
	int   end_y, old_start_y, old_end_y;
	Score f_log_sum_score, f_sum_score, emission_score, transition_score;
	    
	Score first_forward_score, first_forward_transition_score;
	bool  f_get_first_score;

	long int cur_index, old_index;
	
	vector<int> cur_index_vec(2);   // y position, state
	vector<int> old_index_vec(2);   // x position, y position, state       	
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) 
	{
	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    
	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }
	    
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) 
	    {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));
		
	    }
	}
	
	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	for (i=(Start_x-1); i <= (End_x-1); i++) {
	    
	    if (i < 0)
	    {
		if (tube.GetElement(0, 0) == 0) 
		{
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		 }
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }
	}
	
	// calculate max_deltax value 

	for(i=0; i<number_of_emissions_being_train; i++)
	{
	    (*final_count)[i] = 0;
	}

	const long int tube_length = states*(seq_end_y+2);
	
	// calculate max_deltax value

	int max_d = 0;
	
	for (i = 1; i < (states-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;
	const int tube_width = max_deltax+1;	
		 
	// allocate memory for viterbi_tube
	
	long ****training_tube = NULL;
	training_tube = new long***[number_of_sample_paths];
	for(i=0; i<number_of_sample_paths; i++)
	{
	    training_tube[i] = new long**[number_of_emissions_being_train];
	    for(j=0; j<number_of_emissions_being_train; j++)
	    {
		training_tube[i][j] = new long*[tube_width];
		for(int k =0; k<tube_width; k++)	 
		{
		    training_tube[i][j][k] = new long[tube_length];
		}
	    }
	}
	
	double **forward_tube = NULL;
	forward_tube = new double*[tube_width];	
	for(i=0; i<tube_width; i++)
	{	   
	    forward_tube[i]  = new double[tube_length];
	}

	int** linear_emission_tube = NULL;
	linear_emission_tube = new int*[tube_width];
	for(i=0; i<tube_width;i++)
	{
	    linear_emission_tube[i] = new int[tube_length];
	}

	double* p = new double[tube_length];
	int* p_index = new int[tube_length];
	double* p_transition = new double[tube_length];
	int number_of_p_index = 0;
	
	int from = -1;
	int to = -1;
	
	// ----------------------------------------------------------------------
	// start : calculation of training_tube
	// ----------------------------------------------------------------------
	
	int count = 0;
	
	last_state = states-2;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) 
	{
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }
	    
	    if(xsteps==0)
	    {
		cur_offset = 0;
	    }else{
		// swap plate 0 and plate 1
		cur_offset = (cur_offset+1)%tube_width;
	    }
	    
	    // initialize the current plate 
	    for(i=0; i<number_of_sample_paths; i++)
	    {
		for(j=0; j<number_of_emissions_being_train;j++)
		{
		    for(int k =0; k<tube_length; k++)
		    {		  
			training_tube[i][j][cur_offset][k] = 0;
		    }
		}
	    }
	    
	    for(i=0; i<tube_length; i++)
	    {
		forward_tube[cur_offset][i] = Logzero;
		linear_emission_tube[cur_offset][i] = -1;	       
	    }	    
	    
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) 
	    {
		
		cur_index_vec[0] = ysteps;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0)) 
		{
		    last_state = states-1;
		}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) 
		{ // initialise start of state path
		    
		    cur_index_vec[1] = Start_state;
		    cur_index  = (cur_index_vec[0] - start_y) * states + cur_index_vec[1]; 
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state
		    linear_index_emission = -1;
		    
		    if ((deltax + deltay) > 0) 
		    {
			
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sX.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) 
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sY.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 

			emission_score = this->model[Start_state].get_emission_score(&sX, xsteps - deltax, &sY, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    forward_tube[cur_offset][cur_index] = emission_score; 		    
		    linear_emission_tube[cur_offset][cur_index] = linear_index_emission;
		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			number_of_p_index = 0;

			f_sum_score = Logzero;
			f_log_sum_score = Logzero;
			
			emission_score = Logzero; 
			
			cur_index_vec[1]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) 
			{                     			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_offset = cur_offset-(xsteps-old_xsteps);
				if(old_offset<0)
				{
				    old_offset += tube_width;
				}
				old_index_vec[0] = old_ysteps;
				
				// get emission score for new state
				
				linear_index_emission = -1;
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | sX.letter(old_xsteps + i);
				    }
				    
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | sY.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;} 
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;} 

				    emission_score = this->model[new_state].get_emission_score(&sX, old_xsteps, &sY, old_ysteps, 
											       linear_index_emission);
				    
				}
				else // silent state
				{emission_score = 0.0;}
				
				cur_index  = (cur_index_vec[0] - start_y) * states + cur_index_vec[1]; 

				linear_emission_tube[cur_offset][cur_index] =linear_index_emission;

				 // loop over old states, if emission_score in new_state is larger than Logzero

				if (emission_score > Logzero) 
				{ 				 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    f_get_first_score = false; 

				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) 
				    {
					
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index_vec[1]  = old_state;
					old_index = (old_index_vec[0] - old_start_y) * states + old_index_vec[1];
					
					if(new_state==states-1)
					{
					    old_offset=cur_offset;
					}

					transition_score = this->new_get_transition_score(this, old_state, new_state, &sX, old_xsteps, &sY, old_ysteps);
				       
					
					if((transition_score > Logzero)&&(forward_tube[old_offset][old_index] > Logzero))
					{
					    // calculate forward tube
					    p[old_index] = forward_tube[old_offset][old_index];
					    p_index[number_of_p_index] = old_index;
					    p_transition[number_of_p_index] = transition_score;
					    number_of_p_index++;

					    if(!f_get_first_score)
					    {					       
						first_forward_score = forward_tube[old_offset][old_index];	 				
						first_forward_transition_score = transition_score;				
						f_sum_score = first_forward_score + first_forward_transition_score;
						
						f_log_sum_score = 0;
						f_get_first_score = true;						
					    }else{
						f_log_sum_score = f_log_sum_score + exp(forward_tube[old_offset][old_index] + transition_score - f_sum_score);	
					    }
					}
					
				    } // loop over old states

				    if (f_sum_score > Logzero)
				    {					
					f_log_sum_score = log(1+f_log_sum_score);
					
					f_sum_score = f_sum_score + f_log_sum_score+ emission_score;   			
					// assign sum to forward tube
					forward_tube[cur_offset][cur_index] = f_sum_score;	 
				    } // if f_sum_score > Logzero					
				    // calculate the probability distribution vector p
				    
				    for(i=0; i<number_of_p_index; i++)
				    {					
					if(i==0)
					{
					    if(p[p_index[i]]>Logzero)
					    {
						p[p_index[i]] = exp(p[p_index[i]]+p_transition[i]+emission_score-forward_tube[cur_offset][cur_index]);
					    }else{
						p[p_index[i]] = 0;
					    }				
					}else{
					    if(p[p_index[i]]>Logzero)
					    {
						p[p_index[i]]=p[p_index[i-1]]+exp(p[p_index[i]]+p_transition[i]+emission_score-forward_tube[cur_offset][cur_index]);
					    }else{
						p[p_index[i]] = p[p_index[i-1]];
					    }
					    
					}			
				    }
				    
				    // choose number_of_sample_paths pointers and add to the
				    // corresponding entry of table e
				    for(i=0; i<number_of_sample_paths; i++)
				    {
					double rannum = static_cast<double>(rand())/
					    (static_cast<double>(RAND_MAX)+ static_cast<double>(1));
					for(j=0;j<number_of_p_index;j++)
					{
					    if(rannum<p[p_index[j]])
					    {
						sample_state = p_index[j]%states;
						break;
					    }
					}
				
					old_index = (old_index_vec[0]-start_y)*states + sample_state;
					sample_linear_index_emission = linear_emission_tube[old_offset][old_index];

					for(j=0;j<number_of_emissions_being_train;j++)
					{
					    if(from_list[j][sample_state]>0)
					    {
						if(emission_index[j]==sample_linear_index_emission)
						{
						    training_tube[i][j][cur_offset][cur_index] = training_tube[i][j][old_offset][old_index]+1;
						}else
						    training_tube[i][j][cur_offset][cur_index] = training_tube[i][j][old_offset][old_index];
					    }else{
						training_tube[i][j][cur_offset][cur_index] = training_tube[i][j][old_offset][old_index];
					    }
					}
				    }				    	    		     
				} // if emission_score in new_state > Logzero		      
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    count++;
	} // loop over xsteps

	// ----------------------------------------------------------------------
	// end: calculation of training_tube
	// ----------------------------------------------------------------------
	
	// ----------------------------------------------------------------------
	// start : Obtained the result of the training
	//       
	// ----------------------------------------------------------------------
	
	cur_index_vec[0] = End_y;
	cur_index_vec[1] = End_state;    
	
	cur_index = cur_index_vec[0]*states + cur_index_vec[1];
	
	if(forward_tube[cur_offset][cur_index]<Logzero)
	{
	    cout <<"ERROR: Hmm::Posterior_train_emission: the forward probability of sequence/sequence pair "
		 <<"is 0."<<endl;
	    check++;
	}
	
	
	if(!check) // set score
	{
	    for(int k=0; k<number_of_emissions_being_train; k++)
	    {
		(*final_count)[k] = 0;
		for(int l =0; l<number_of_sample_paths; l++)
		{
		    (*final_count)[k] += training_tube[l][k][cur_offset][cur_index];		    
		}

		if((*final_count)[k]<=0)
		{
		    cout <<"Warning: Hmm::Posterior_train_emission: the probability of emission ";
		    for(i=0; i<states; i++)
		    {
			if(from_list[k][i]!=-1)
			{			 
			    cout<<" state : "<<from_list[k][i]<<" emit : "<<emission_index[k];		 
			}
		    }
		    cout<<" is 0."<<endl;
		}   
	    }
	}
	
	// free memory
	
	if (forward_tube) {
	    for (i = 0; i < tube_width; i++) {
		if (forward_tube[i]) delete forward_tube[i];
		forward_tube[i] = NULL;
	    }
	    delete [] forward_tube; 
	}
	forward_tube = NULL;
	
	if (linear_emission_tube) {
	    for (i = 0; i < tube_width; i++) {
		if (linear_emission_tube[i]) delete linear_emission_tube[i];
		linear_emission_tube[i] = NULL;
	    }
	    delete [] linear_emission_tube; 
	}
	linear_emission_tube = NULL;

	if (training_tube) 
	{
	    for(i=0; i<number_of_sample_paths; i++)
	    {
		if(training_tube[i])
		{
		    for(j=0; j< number_of_emissions_being_train; j++)
		    {
			if(training_tube[i][j])
			{
			    for (int k = 0; k < tube_width; k++) 
			    {
				if (training_tube[i][j][k]) delete [] training_tube[i][j][k];
				training_tube[i][j][k] = NULL;
			    }
			    delete [] training_tube[i][j];
			}
			training_tube[i][j] = NULL;		
		    }
		    delete [] training_tube[i];
		}
		training_tube[i] = NULL;
	    }
	    delete [] training_tube;
	}
	training_tube = NULL;

	if(p) delete[] p;
	p = NULL;

	if(p_index) delete [] p_index;
	p_index = NULL;

	if(p_transition) delete [] p_transition;
	p_transition = NULL;
		
    }
    
    return check;
}

int Hmm::BaumWelch_train_emission(const Sequence &sX, 
				  const Sequence &sY,
				  const int number_of_emissions_being_train,
				  int* emission_index,
				  int** from_list,		      
				  double** final_score,
				  double& forward_score,
				  Tube<int> tube,
				  const vector<int> Start_point,
				  const vector<int> End_point) 								 
{
    int check = 0;
    
    if (tube.Empty()) 
    {	
	int max_y = sY.length()-1;
	if(max_y < 0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(sX.length(), 2);
	
	for (int i = 0; i < sX.length(); i++) {
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);
	}
	tube = default_tube;
    }

    if (tube.Lx() != sX.length()) 
    {
	
	cout << "ERROR: Hmm::BaumWelch_train_emission: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << sX.length() << ").\n";
	 check++;
    }
    // check that Start_x < End_x 
    
    if (Start_point[0] > End_point[0]) 
    {
	
	cout << "ERROR: Hmm::BaumWelch_train_emission: Start_x (" << Start_point[0] << ") > End_x (" 
	     << End_point[0] << ").\n";
	check++;
    }
    
    if (check == 0) 
    {
	 
	int first_x_pos = -1; 
	int first_y_pos = -1; 
	int last_x_pos  = -1; 
	int last_y_pos  = -1; 

	// constants

	const int seq_start_x = 0;                 
	const int seq_end_x   = sX.length()-1;
	const int seq_start_y = 0;                 
	int seq_end_y   = sY.length()-1;
	if(seq_end_y < 0)
	{
	    seq_end_y = 0;
	}
	const int states      = this->get_number_of_states(); 
	const int Start_x     = Start_point[0];
	const int Start_y     = Start_point[1];
	const int End_x       = End_point[0];
	const int End_y       = End_point[1];
	const int Start_state = Start_point[2];
	const int End_state   = End_point[2];
	const int shift_index = bitshift(alphabet);
	
	// variables 
	
	int   j, i, linear_index_emission, number_of_old_states, old_state_number;
	int   xsteps, ysteps, old_xsteps, old_ysteps, last_state;
	int   cur_offset, old_offset, temp_offset;
	int   new_state, old_state, deltax, deltay;
	int   end_y, old_start_y, old_end_y;
	Score f_log_sum_score, f_sum_score, emission_score, transition_score;
	Score* e_log_sum_score1 = new Score[number_of_emissions_being_train];
	Score* e_log_sum_score2 = new Score[number_of_emissions_being_train];
	Score* e_sum_score = new Score[number_of_emissions_being_train];
	
	Score first_forward_score, first_forward_transition_score;
	//first_forward_transition_score, first_training_transition_score;
	bool  f_get_first_score;
	bool* e_get_first_score = new bool[number_of_emissions_being_train];
	
	long int cur_index, old_index;
	
	vector<int> cur_index_vec(2);   // y position, state
	vector<int> old_index_vec(2);   // x position, y position, state       	
	
	int start_x = Start_x; // most 5' letter read from x
	int start_y = Start_y; // most 5' letter read from y
	int d_x     = this->model[Start_state].get_letters_to_read_x();
	int d_y     = this->model[Start_state].get_letters_to_read_y();
	
	if (d_x > 0) {start_x = max(0, Start_x - d_x);}
	else         {start_x = max(0, Start_x-1);}
	if (d_y > 0) {start_y = max(0, Start_y - d_y);}
	else         {start_y = max(0, Start_y-1);}
	
	for (i=start_x; i < End_x; i++) 
	{
	    tube.SetElement(i, 0, max(start_y, tube.GetElement(i, 0))); // modify tube according to start_y coordinate
	    
	    tube.SetElement(i, 1, min(End_y-1, tube.GetElement(i, 1))); // modify tube according to end_y coordinate
	    
	    if(End_y==0)
	    {
		tube.SetElement(i,1,min(0,tube.GetElement(i,1)));
	    }
	    
	    if (tube.GetElement(i, 0) > tube.GetElement(i, 1)) 
	    {
		
		tube.SetElement(i, 0, tube.GetElement(i, 1));
		
	    }
	}
	
	// calculate offsets for Viterbi and int_tube
	// int_tube[x][0] = lower_y(x) and int_tube[x][1] = upper_y(x)
	// where x is the xsteps value and the corresponding ysteps values are contained in [lower_y(x), upper_y(x)]
	
	Tube<int> int_tube(End_x - Start_x + 1, 2);
	
	for (i=(Start_x-1); i <= (End_x-1); i++) {
	    
	    if (i < 0)
	    {
		if (tube.GetElement(0, 0) == 0) 
		{
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(0, 0)+1, seq_end_y+1));
		 }
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(0, 1)+1, seq_end_y+1));
	    }
	    else {
		if (tube.GetElement(i, 0) == 0) {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0), seq_end_y+1));
		}
		else {
		    int_tube.SetElement(i-(Start_x-1), 0, min(tube.GetElement(i, 0)+1, seq_end_y+1));
		}	
		int_tube.SetElement(i-(Start_x-1), 1, min(tube.GetElement(i, 1)+1, seq_end_y+1));
	    }
	}
	
	// calculate max_deltax value 

	for(i=0; i<number_of_emissions_being_train; i++)
	{
	    (*final_score)[i] = Logzero;
	}

	const long int tube_length = states*(seq_end_y+2);
	
	// calculate max_deltax value

	int max_d = 0;
	
	for (i = 1; i < (states-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
	}
	const int max_deltax = max_d;
	const int tube_width = max_deltax+1;	
		 
	// allocate memory for viterbi_tube
	
	double ***training_tube = NULL;
	training_tube = new double**[number_of_emissions_being_train];
	for(i=0; i<number_of_emissions_being_train; i++)
	{
	    training_tube[i] = new double*[tube_width];
	    for(j=0; j<tube_width; j++)	 
	    {
		training_tube[i][j] = new double[tube_length];
	    }
	}
	

	double **forward_tube = NULL;
	forward_tube = new double*[tube_width];	
	for(i=0; i<tube_width; i++)
	{	   
	    forward_tube[i]  = new double[tube_length];
	}
	
	int from = -1;
	int to = -1;
	
	// ----------------------------------------------------------------------
	// start : calculation of training_tube
	// ----------------------------------------------------------------------
	
	int count = 0;
	
	last_state = states-2;
	
	for (xsteps = Start_x; xsteps <= End_x; xsteps += 1) 
	{
	    start_y   = int_tube.GetElement(xsteps-Start_x, 0);
	    end_y     = int_tube.GetElement(xsteps-Start_x, 1);
	    
	    if (xsteps == End_x) {
		end_y = End_y;
	    }
	    
	    if(xsteps==0)
	    {
		cur_offset = 0;
	    }else{
		// swap plate 0 and plate 1
		cur_offset = (cur_offset+1)%tube_width;
	    }
	    
	    // initialize the current plate 
	    for(i=0; i<number_of_emissions_being_train;i++)
	    {
		for(j=0; j<tube_length; j++)
		{		  
		    training_tube[i][cur_offset][j] = Logzero;
		}
	    }
	    
	    for(i=0; i<tube_length; i++)
	    {
		forward_tube[cur_offset][i] = Logzero;
	    }
	    
	    for (ysteps = start_y; ysteps <= end_y; ysteps += 1) 
	    {
		
		cur_index_vec[0] = ysteps;
		
		if ((xsteps == (seq_end_x+1)) && (ysteps == (seq_end_y+1))) 
		{
		    last_state = states-1;
		}else if((xsteps==(seq_end_x+1))&&(seq_end_y==0)) 
		{
		    last_state = states-1;
		}
		
		if ((xsteps == Start_x) && (ysteps == Start_y)) 
		{ // initialise start of state path
		    
		    cur_index_vec[1] = Start_state;
		    cur_index  = (cur_index_vec[0] - start_y) * states + cur_index_vec[1]; 
		    deltax    = this->model[Start_state].get_letters_to_read_x();
		    deltay    = this->model[Start_state].get_letters_to_read_y();
		    
		    // get emission score for Start_state
		    
		    if ((deltax + deltay) > 0) 
		    {
			
			linear_index_emission = 0;
			for (i = 0; i < deltax; i++)               
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sX.letter(xsteps - deltax + i);
			}
			for (j = 0; j < deltay; j++) 
			{
			    linear_index_emission = (linear_index_emission << shift_index) | sY.letter(ysteps - deltay + j);
			}
			if (deltax > 0) {first_x_pos = xsteps-deltax;} 
			if (deltay > 0) {first_y_pos = ysteps-deltay;} 

			emission_score = this->model[Start_state].get_emission_score(&sX, xsteps - deltax, &sY, ysteps - deltay, 
										     linear_index_emission);
		    }
		    else {emission_score = 0.0;} 
		    forward_tube[cur_offset][cur_index] = emission_score; 		    
		}
		else {
		    
		    for (new_state = 1; new_state <= last_state; new_state++) 
		    {
			
			f_sum_score      = Logzero;
			f_log_sum_score  = Logzero;
			for(i=0; i<number_of_emissions_being_train; i++)
			{
			    e_sum_score[i]      = Logzero;
			    e_log_sum_score1[i] = Logzero;
			    e_log_sum_score2[i] = Logzero;
			}
			emission_score = Logzero; 
			
			cur_index_vec[1]      = new_state;
			deltax         = this->model[new_state].get_letters_to_read_x();
			deltay         = this->model[new_state].get_letters_to_read_y();
			
			old_xsteps     = xsteps - deltax;
			old_ysteps     = ysteps - deltay;
			
			if ((Start_x <= old_xsteps) && (old_xsteps <= End_x) &&                                     
			    (((deltax > 0) && (seq_start_x <= (xsteps - deltax)) && ((xsteps - 1) <= seq_end_x)) || 
			     (deltax == 0))) 
			{                     			    
			    old_start_y = int_tube.GetElement(old_xsteps-Start_x, 0);
			    old_end_y   = int_tube.GetElement(old_xsteps-Start_x, 1);
			    
			    if (old_xsteps == End_x)   {old_end_y   = End_y;}
			    
			    if ((old_start_y <= old_ysteps) && (old_ysteps <= old_end_y) &&                             
				(((deltay > 0) && (seq_start_y <= (ysteps - deltay)) && ((ysteps - 1) <= seq_end_y)) || 
				 (deltay == 0))) {                                                                      
				
				old_offset = cur_offset-(xsteps-old_xsteps);
				if(old_offset<0)
				{
				    old_offset += tube_width;
				}
				old_index_vec[0] = old_ysteps;
				
				// get emission score for new state
				
				if ((deltax + deltay) > 0) {
				    
				    linear_index_emission = 0;
				    for (i = 0; i < deltax; i++)               {
					linear_index_emission = (linear_index_emission << shift_index) | sX.letter(old_xsteps + i);
				    }
				    
				    for (j = 0; j < deltay; j++) {
					linear_index_emission = (linear_index_emission << shift_index) | sY.letter(old_ysteps + j);
				    }
				    if (deltax > 0) {last_x_pos = old_xsteps+deltax-1;} 
				    if (deltay > 0) {last_y_pos = old_ysteps+deltay-1;} 

				    emission_score = this->model[new_state].get_emission_score(&sX, old_xsteps, &sY, old_ysteps, 
											       linear_index_emission);
				    
				}
				else // silent state
				{emission_score = 0.0;}

				
				cur_index  = (cur_index_vec[0] - start_y) * states + cur_index_vec[1]; 

				 // loop over old states, if emission_score in new_state is larger than Logzero

				if (emission_score > Logzero) 
				{ 				 
				    
				    number_of_old_states = this->model[new_state].get_number_of_previous_states();
				    
				    f_get_first_score = false;
				    for(i=0; i<number_of_emissions_being_train; i++)
				    {
					e_get_first_score[i] = false;
				    }
				    
				    for (old_state_number = 0; old_state_number < number_of_old_states; old_state_number++) 
				    {
					
					old_state     = this->model[new_state].get_number_of_previous_state(old_state_number);
					old_index_vec[1]  = old_state;
					old_index = (old_index_vec[0] - old_start_y) * states + old_index_vec[1];
					
					if(new_state==states-1)
					{
					    old_offset=cur_offset;
					}
					
					transition_score = this->new_get_transition_score(this, old_state, new_state, &sX, old_xsteps, &sY, old_ysteps);
					
					if((transition_score > Logzero)&&(forward_tube[old_offset][old_index] > Logzero))
					{
					    // calculate forward tube
					    if(!f_get_first_score)
					    {
						first_forward_score = forward_tube[old_offset][old_index];
						first_forward_transition_score = transition_score;
						f_sum_score = first_forward_score + first_forward_transition_score;
						f_log_sum_score = 0;
						f_get_first_score = true;
					    }else{
						f_log_sum_score = f_log_sum_score + exp(forward_tube[old_offset][old_index] + transition_score - first_forward_score-first_forward_transition_score);
					    }
					    
					}
					if(transition_score > Logzero)
					{
					    for(i=0; i<number_of_emissions_being_train; i++)
					    {
						if(training_tube[i][old_offset][old_index] > Logzero)
						{
						    // calculate training tube
						    if(!e_get_first_score[i])
						    {
							e_sum_score[i] = training_tube[i][old_offset][old_index] + transition_score;
							e_log_sum_score1[i] = 0;
							e_get_first_score[i] = true;
						    }else{
							e_log_sum_score1[i] = e_log_sum_score1[i] + exp(training_tube[i][old_offset][old_index] + transition_score - e_sum_score[i]);
						    }	
						}
					    }				
					}														  
				    } // loop over old states
				    
				    if (f_sum_score > Logzero)
				    {
					f_log_sum_score = log(1+f_log_sum_score);
					
					f_sum_score = f_sum_score + f_log_sum_score+ emission_score;			
					
					// assign sum to forward tube
					forward_tube[cur_offset][cur_index] = f_sum_score;
				    } // if f_sum_score > Logzero			
				     
				    for(i=0; i<number_of_emissions_being_train; i++)
				    {
					if((from_list[i][new_state]>0)&&
					   (emission_index[i]==linear_index_emission))
					{				
					    
					    if(e_sum_score[i] <= Logzero)
					    {
						e_sum_score[i] = forward_tube[cur_offset][cur_index];
						
					    }else{
						e_log_sum_score1[i] = log(1+e_log_sum_score1[i]);
						e_log_sum_score1[i] = e_sum_score[i] + e_log_sum_score1[i] + emission_score; // the log sum of the training table part
						e_log_sum_score1[i] = exp(e_log_sum_score1[i]-forward_tube[cur_offset][cur_index]);
						e_log_sum_score1[i] = log(1+e_log_sum_score1[i]);
						e_sum_score[i] = forward_tube[cur_offset][cur_index] + e_log_sum_score1[i];					  
					    }
					    
					    training_tube[i][cur_offset][cur_index] = e_sum_score[i];	
					    
					}else{
					    if(e_sum_score[i] > Logzero)
					    {
						e_log_sum_score1[i] = log(1+e_log_sum_score1[i]);
						e_sum_score[i] = e_sum_score[i] + e_log_sum_score1[i] + emission_score;
					    }
					    training_tube[i][cur_offset][cur_index] = e_sum_score[i];	
					}			
				    }					    
				}else // if emission_score in new_state > Logzero	       
				{
				    forward_tube[cur_offset][cur_index] = Logzero;
				    for(i=0; i<number_of_emissions_being_train; i++)
				    {
					training_tube[i][cur_offset][cur_index] = Logzero;
				    }
				}
			    } // if old_ysteps within tube
			} // if old_xsteps within tube
		    } // loop over new states
		} // if not ((xsteps == Start_x) && (ysteps == Start_y))
	    } // loop over ysteps
	    count++;
	} // loop over xsteps
	
	// ----------------------------------------------------------------------
	// end: calculation of training_tube
	// ----------------------------------------------------------------------
	
	// ----------------------------------------------------------------------
	// start : Obtained the result of the training
	//       
	// ----------------------------------------------------------------------
	
	cur_index_vec[0] = End_y;
	cur_index_vec[1] = End_state;    
	
	cur_index = cur_index_vec[0]*states + cur_index_vec[1];
	
	if(forward_tube[cur_offset][cur_index]<Logzero)
	{
	    cout <<"ERROR: Hmm::BaumWelch_train_emission: the forward probability of sequence/sequence pair "
		 <<"is 0."<<endl;
	    check++;
	}
	
	if(!check) // set score
	{
	    forward_score = forward_tube[cur_offset][cur_index];

	    for(int k=0; k<number_of_emissions_being_train; k++)
	    {
		if(training_tube[k][cur_offset][cur_index]<=Logzero)
		{
		    cout <<"Warning: Hmm::BaumWelch_train_emission: the probability of emission ";
		    for(i=0; i<states; i++)
		    {
			if(from_list[k][i]!=-1)
			{			 
			    cout<<" state : "<<from_list[k][i]<<" emit : "<<emission_index[k];			 
			}
		    }
		    cout<<" is 0."<<endl;
		    (*final_score)[k] = Logzero;
		}else{
		    (*final_score)[k] = training_tube[k][cur_offset][cur_index] - forward_score;
		}	   
	    }
	}
	
	// free memory
	
	if (forward_tube) {
	    for (i = 0; i < tube_width; i++) {
		if (forward_tube[i]) delete forward_tube[i];
		forward_tube[i] = NULL;
	    }
	    delete [] forward_tube; 
	}
	forward_tube = NULL;
	
	if (training_tube) 
	{
	    for(i=0; i< number_of_emissions_being_train; i++)
	    {
		if(training_tube[i])
		{
		    for (j = 0; j < tube_width; j++) 
		    {
			if (training_tube[i][j]) delete [] training_tube[i][j];
			training_tube[i][j] = NULL;
		    }
		    delete [] training_tube[i];
		 }
		training_tube[i] = NULL;		
	    }
	     delete [] training_tube; 
	}
	training_tube = NULL;   
	
	if(e_get_first_score) delete [] e_get_first_score;
	e_get_first_score = NULL;
	
	if(e_sum_score) delete [] e_sum_score;
	e_sum_score = NULL;
	 
	if(e_log_sum_score1) delete [] e_log_sum_score1;
	e_log_sum_score1 = NULL;
	
	if(e_log_sum_score2) delete [] e_log_sum_score2;
	e_log_sum_score2 = NULL;
	
    }    
    return check;
}

int Hmm::Viterbi_train_FTP(const char* seqfile,
			   const char* seqannfile,
			   const char* name_input_tube_file,
			   const long int max_volume, 
			   const int radius,
			   int& SeqCount,
			   long*** prev_GTP_count,
			   long*** cur_GTP_count, 
			   model_parameters* const MP,
			   TransitionProb* TP)
{


    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: Viterbi_train_TTP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
       
    int i = 0;
    int j = 0;

    int GTPsize = 2*TP->get_GTPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;   
     
    long* GTPcount = new long[GTPsize];
    long* total_GTP_count = new long[GTPsize];
    
    int** num_of_from_list = new int*[GTPsize];
    int*** from_list = new int**[GTPsize];
    for(i=0; i<GTPsize; i++)
    {
	GTPcount[i] = 0;
	total_GTP_count[i] = 0;
	num_of_from_list[i] = new int[number_of_states];      
	from_list[i] = new int*[number_of_states];
	for(j=0; j<number_of_states; j++)
	{
	    num_of_from_list[i][j] = 0;
	    from_list[i][j] = new int[number_of_states];
	    for(int k=0; k<number_of_states; k++)
	    {
		from_list[i][j][k] = -1;
	    }
	}		
    }
    int index = 0;
    int number_of_GTP_being_train = 0;
    
    // ---------------------------------------------------------------------
    // Start: training of FTP, train all the GTP
    // ---------------------------------------------------------------------
    
    // for each pair of sequences		    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: Viterbi_train_FTP: "
	    <<"training sequence file : "<<seqfile
	    <<" can not be read, training terminated."<<endl;
	check++;
    }
	
    // train transition_scores(to,from)		    
    SeqCount = 0;

    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);

	for(i=0; i<GTPsize; i++)
	{
	    GTPcount[i] = 0;
	}

	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point);	
	
	int length_y = sY.length();
	
	long long int tube_length = number_of_states*(length_y+2);
	    
	long long int tube_area = tube_width*tube_length;
	long long int total_volume = 2*tube_area*8;
	
	
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: Viterbi_train_FTP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" is "<<total_volume<<endl;	       
	    break;
	}

	// calculate the max number of GTP can be trained in one train_transition function
	number_of_GTP_being_train = 1;
	
	while(number_of_GTP_being_train<GTPsize)
	{
	    total_volume += tube_area;
	    if(total_volume <= max_volume)
	    {
		number_of_GTP_being_train++;
	    }else{
		break;
	    }
	}

	cout<<"Number of group transition parameters being trained for the "<<SeqCount
	    <<"-th input sequence (pairs) in parallel : "<<number_of_GTP_being_train<<endl;
	
	// for each GTP
	i = 0;
	while(i<GTPsize)  
	{		    
	    // calculate number of GTP can be trained in for this pair of sequences
	    
	    if(i+number_of_GTP_being_train>GTPsize)
	    {
		number_of_GTP_being_train = GTPsize-i;
	    }

	    for(j=0;j<number_of_GTP_being_train;j++)
	    {
		GTPcount[j] = 0;
	    }
    	    
	    // initialize
	    for(j=0; j<number_of_GTP_being_train; j++)
	    {	
		for(int k=0; k<number_of_states; k++)
		{
		    num_of_from_list[j][k] = 0;	   		  
		    for(int l=0; l<number_of_states; l++)
		    {
			if(from_list[j][k][l] == -1)
			{
			    break;
			}
			from_list[j][k][l] = -1;
		    }
		}
	    }
	    
	    for(j=0; j<number_of_GTP_being_train; j++)
	    {
		if(j%2==0)
		{
		    for(int k=0; k<TP->get_GTP_NumOffromto((j+i)/2); k++)
		    {
			index = TP->get_GTP_to((j+i)/2,k);
			from_list[j][index][num_of_from_list[j][index]] = TP->get_GTP_from((j+i)/2,k);
			num_of_from_list[j][index]++;
		    }
		}else{
		    for(int k=0; k<TP->get_GTP_NumOfOverfromto((j+i)/2); k++)
		    {
			index = TP->get_GTP_Overto((j+i)/2,k);
			from_list[j][index][num_of_from_list[j][index]] = TP->get_GTP_Overfrom((j+i)/2,k);
			num_of_from_list[j][index]++;
		    }
		}
	    }
	   	 	   
	    check+= Viterbi_train_transition(sX,sY,
					     number_of_GTP_being_train,
					     num_of_from_list,
					     from_list,
					     &GTPcount,
					     tTube,
					     Start_point,
					     End_point);
	    
	    if(check)
	    {
		cout<<"Warning: Hmm class:: error in Viterbi_train_transition of the "
		    <<SeqCount<<"-th sequence (pair): Viterbi_train_FTP";		
		for(j=i; j<i+number_of_GTP_being_train; j+=2)
		{
		    cout<<"FTP: "<<(j/2)<<endl;
		}
	    }
	    
	    // set the counts of all sequences to TTP[i].score ( a_i,j)
	    	    
	    for(j=0; j<number_of_GTP_being_train; j++)
	    {
		total_GTP_count[j+i]+=GTPcount[j];
	    }	
	    
	    i+= number_of_GTP_being_train;	    	    	 	    

	}// loop over the TTP
	for(i=0; i<GTPsize; i++)
	{
	    (*prev_GTP_count)[SeqCount][i] = (*cur_GTP_count)[SeqCount][i];
	    (*cur_GTP_count)[SeqCount][i] = total_GTP_count[i];
	}
	SeqCount++;	
    } // after loop of sequence		          
    fclose(fIn_sequence);
	    

    // normalization
    // the i-th GTP score is calculated by total_GTP_count[2*i] / total_GTP_count[2*i+1]

    for(i=0; i<TP->get_GTPsize(); i++)
    {
	TP->set_GTP_score(i,static_cast<double>(static_cast<double>(total_GTP_count[2*i])
						/static_cast<double>(total_GTP_count[2*i+1])));
    }

    //derive the FTPs from the GTPs
    check += TP->derive_FTPs_from_GTPs();
    if(check)
    {
	cout<<"ERROR: Hmm class:: Viterbi_train_FTP : "
	    <<"error in derive_FTPs_from_GTPs"<<endl;
    }
    
    // add pseudocount for FTPs
    for(i=0; i<TP->get_GTPsize(); i++)
    {
	TP->set_FTP_prob(i,TP->get_FTP_prob(i)+TP->get_FTP_pseudocount(i));
    }
    //-----------------------------------------------------------------------------
    // End : Training of TTP in this loop
    //-----------------------------------------------------------------------------

    // free memory

    if(total_GTP_count) delete [] total_GTP_count;
    total_GTP_count = NULL;

    if(GTPcount) delete [] GTPcount;
    GTPcount = NULL;

    if(num_of_from_list) 
    {
	for(i=0; i<GTPsize; i++)
	{
	    if(num_of_from_list[i]) delete[] num_of_from_list[i];
	    num_of_from_list[i] = NULL;
	}
	delete [] num_of_from_list;
    }
    num_of_from_list = NULL;

    if(from_list)
    {
	for(i=0; i<GTPsize; i++)
	{
	    if(from_list[i])
	    {
		for(j=0; j<number_of_states; j++)
		{
		    if(from_list[i][j]) delete [] from_list[i][j];
		    from_list[i][j] = NULL;
		}
		delete [] from_list[i];
	    }
	    from_list[i] = NULL;	   
	}
	delete [] from_list;
    }
    from_list = NULL;   

    return check;
}

int Hmm::Posterior_train_FTP(const char* seqfile,
			     const char* seqannfile,
			     const char* name_input_tube_file,
			     const long int max_volume, 
			     const int radius,
			     int& SeqCount,
			     model_parameters* const MP,
			     TransitionProb* TP,
			     const int SamplePaths)
{


    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: Posterior_train_FTP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
       
    int i = 0;
    int j = 0;

    int GTPsize = 2*TP->get_GTPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;   
     
    long* GTPcount = new long[GTPsize];
    long* total_GTP_count = new long[GTPsize];
    
    int** num_of_from_list = new int*[GTPsize];
    int*** from_list = new int**[GTPsize];
    for(i=0; i<GTPsize; i++)
    {
	GTPcount[i] = 0;
	total_GTP_count[i] = 0;
	num_of_from_list[i] = new int[number_of_states];      
	from_list[i] = new int*[number_of_states];
	for(j=0; j<number_of_states; j++)
	{
	    num_of_from_list[i][j] = 0;
	    from_list[i][j] = new int[number_of_states];
	    for(int k=0; k<number_of_states; k++)
	    {
		from_list[i][j][k] = -1;
	    }
	}		
    }

    int index = 0;
    int number_of_GTP_being_train = 0;
    int number_of_sample_paths_being_train = 0;
    
    // ---------------------------------------------------------------------
    // Start: training of FTP, train all the GTP
    // ---------------------------------------------------------------------
    
    // for each pair of sequences		    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: Posterior_train_FTP: "
	    <<"training sequence file : "<<seqfile
	    <<" can not be read, training terminated."<<endl;
	check++;
    }
	
    // train transition_scores(to,from)		    
    SeqCount = 0;

    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);

	for(i=0; i<GTPsize; i++)
	{
	    GTPcount[i] = 0;
	}

	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point);	
	
	int length_y = sY.length();
	
	long long int tube_length = number_of_states*(length_y+2);
	    
	long long int tube_area = tube_width*tube_length;
	long long int total_volume = 2*tube_area*8;
	
	
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: Posterior_train_FTP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" is "<<total_volume<<endl;	       
	    break;
	}


	number_of_sample_paths_being_train = 0;
	// check whether all the GTP can be train with one sample path
	if(total_volume+(GTPsize-1)*tube_area<=max_volume) 
	{
	    // train all TTP with a sample path
	    number_of_GTP_being_train = GTPsize;
	    total_volume += (GTPsize-1)*tube_area;
	    tube_area = total_volume - tube_area;
	    // calculate number of sample paths can be sampled in one iteration
	    while(number_of_sample_paths_being_train<SamplePaths)
	    {
		total_volume += tube_area;
		if(total_volume <= max_volume)
		{
		    number_of_sample_paths_being_train++;
		}else{
		    break;
		}
	    }
	}else{
	    // can not train all GTP with a sample path
	    number_of_sample_paths_being_train =1;
	    // calculate number of GTP can be train with one sample path
	    number_of_GTP_being_train = 1;
	    while(number_of_GTP_being_train<GTPsize)
	    {
		total_volume += tube_area;
		if(total_volume <= max_volume)
		{
		    number_of_GTP_being_train++;
		}else{
		    break;
		}
	    }
	}

				    
	cout<<"Number of group transition parameters being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "
	    <<number_of_GTP_being_train<<endl;
	cout<<"Number of Sample paths being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "
	    <<number_of_sample_paths_being_train<<endl;

	int s = 0;
	while(s<SamplePaths) // for each Sample path
	{
	    
	    // for each GTP
	    i = 0;

	    while(i<GTPsize)  
	    {		    
		// calculate number of GTP can be trained in for this pair of sequences
		
		if(i+number_of_GTP_being_train>GTPsize)
		{
		    number_of_GTP_being_train = GTPsize-i;
		}

		for(j=0;j<number_of_GTP_being_train;j++)
		{
		    GTPcount[j] = 0;
		}
		
		// initialize
		for(j=0; j<number_of_GTP_being_train; j++)
		{	
		    for(int k=0; k<number_of_states; k++)
		    {
			num_of_from_list[j][k] = 0;	   		  
			for(int l=0; l<number_of_states; l++)
			{
			    if(from_list[j][k][l] == -1)
			    {
				break;
			    }
			    from_list[j][k][l] = -1;
			}
		    }
		}
		
		for(j=0; j<number_of_GTP_being_train; j++)
		{
		    if(j%2==0)
		    {
			for(int k=0; k<TP->get_GTP_NumOffromto((j+i)/2); k++)
			{
			    index = TP->get_GTP_to((j+i)/2,k);
			    from_list[j][index][num_of_from_list[j][index]] = TP->get_GTP_from((j+i)/2,k);
			    num_of_from_list[j][index]++;
			}
		    }else{
			for(int k=0; k<TP->get_GTP_NumOfOverfromto((j+i)/2); k++)
			{
			    index = TP->get_GTP_Overto((j+i)/2,k);
			    from_list[j][index][num_of_from_list[j][index]] = TP->get_GTP_Overfrom((j+i)/2,k);
			    num_of_from_list[j][index]++;
			}
		    }
		    GTPcount[j] = 0;
		}

		check+= Posterior_train_transition(sX,sY,
						   number_of_GTP_being_train,
						   number_of_sample_paths_being_train,
						   num_of_from_list,
						   from_list,
						   &GTPcount,
						   tTube,
						   Start_point,
						   End_point);

		cout<<"after train"<<endl;
	    
		if(check)
		{
		    cout<<"Warning: Hmm class:: error in train_transition of the "
			<<SeqCount<<"-th sequence : Posterior_train_FTP";		
		    for(j=i; j<i+number_of_GTP_being_train; j+=2)
		    {
			cout<<"FTP: "<<(j/2)<<endl;
		    }
		}
		
		// set the counts of all sequences to TTP[i].score ( a_i,j)
		
		for(j=0; j<number_of_GTP_being_train; j++)
		{
		    total_GTP_count[j+i]+=GTPcount[j];
		}	
		
		i+= number_of_GTP_being_train;	    	    	 	    
		
	    }// loop over the FTP
	    s+= number_of_sample_paths_being_train;
	}
	SeqCount++;	
    } // after loop of sequence		          
    fclose(fIn_sequence);
	    
    // normalization
    // the i-th GTP score is calculated by total_GTP_count[2*i] / total_GTP_count[2*i+1]

    for(i=0; i<TP->get_GTPsize(); i++)
    {
	TP->set_GTP_score(i,static_cast<double>(static_cast<double>(total_GTP_count[2*i])
						/static_cast<double>(total_GTP_count[2*i+1])));
    }

    //derive the FTPs from the GTPs
    check += TP->derive_FTPs_from_GTPs();
    if(check)
    {
	cout<<"ERROR: Hmm class:: Posterior_train_FTP : "
	    <<"error in derive_FTPs_from_GTPs"<<endl;
    }
    
    // add pseudocount for FTPs
    for(i=0; i<TP->get_GTPsize(); i++)
    {
	TP->set_FTP_prob(i,TP->get_FTP_prob(i)+TP->get_FTP_pseudocount(i));
    }
    //-----------------------------------------------------------------------------
    // End : Training of TTP in this loop
    //-----------------------------------------------------------------------------

    // free memory

    if(total_GTP_count) delete [] total_GTP_count;
    total_GTP_count = NULL;

    if(GTPcount) delete [] GTPcount;
    GTPcount = NULL;

    if(num_of_from_list) 
    {
	for(i=0; i<GTPsize; i++)
	{
	    if(num_of_from_list[i]) delete[] num_of_from_list[i];
	    num_of_from_list[i] = NULL;
	}
	delete [] num_of_from_list;
    }
    num_of_from_list = NULL;

    if(from_list)
    {
	for(i=0; i<GTPsize; i++)
	{
	    if(from_list[i])
	    {
		for(j=0; j<number_of_states; j++)
		{
		    if(from_list[i][j]) delete [] from_list[i][j];
		    from_list[i][j] = NULL;
		}
		delete [] from_list[i];
	    }
	    from_list[i] = NULL;	   
	}
	delete [] from_list;
    }
    from_list = NULL;   

    return check;
}

int Hmm::Viterbi_train_TTP(const char* seqfile,
			   const char* seqannfile,
			   const char* name_input_tube_file,
			   const long int max_volume, 
			   const int radius,
			   int& SeqCount,
			   long*** prev_TTP_count,
			   long*** cur_TTP_count, 
			   model_parameters* const MP,
			   TransitionProb* TP)
{


    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: Viterbi_train_TTP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
    
    
    int i = 0;
    int j = 0;

    int TTPsize = TP->get_TTPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;   
     
    long* trancount = new long[TTPsize];
    double* tmp_TTP_score = new double[TTPsize];
    double bij = Logzero;
    
    int** num_of_from_list = new int*[TTPsize];
    int*** from_list = new int**[TTPsize];
    for(i=0; i<TTPsize; i++)
    {
	num_of_from_list[i] = new int[number_of_states];      
	from_list[i] = new int*[number_of_states];
	for(j=0; j<number_of_states; j++)
	{
	    num_of_from_list[i][j] = 0;
	    from_list[i][j] = new int[number_of_states];
	    for(int k=0; k<number_of_states; k++)
	    {
		from_list[i][j][k] = -1;
	    }
	}		
    }

    int index = 0;
    int number_of_TTP_being_train = 0;
    long* total_TTP_count = new long[TTPsize];
    
    // ---------------------------------------------------------------------
    // Start: training of TTP
    // ---------------------------------------------------------------------
    
    // for each pair of sequences		    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: Viterbi_train_TTP: "
	    <<"training sequence file : "<<seqfile
	    <<" can not be read, training terminated."<<endl;
	check++;
    }
	
    // train transition_scores(to,from)		    
    SeqCount = 0;
    
    for(i=0;i<TTPsize;i++)
    {
	total_TTP_count[i] = 0;
    }	
    
    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);

	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point);	
	
	int length_y = sY.length();
	
	long long int tube_length = number_of_states*(length_y+2);
	    
	long long int tube_area = tube_width*tube_length;
	long long int total_volume = 2*tube_area*8;
	
	
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: Viterbi_train_TTP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" is "<<total_volume<<endl;	       
	    break;
	}

	// calculate the max number of TTP can be trained in one train_transition function
	number_of_TTP_being_train = 1;
	
	while(number_of_TTP_being_train<TTPsize)
	{
	    total_volume += tube_area;
	    if(total_volume <= max_volume)
	    {
		number_of_TTP_being_train++;
	    }else{
		break;
	    }
	}

	cout<<"Number of transition probabilities being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "<<number_of_TTP_being_train<<endl;
	

	i = 0;
	while(i<TTPsize)  // for each TTP
	{		    
	    // calculate number of TTP can be trained in for this pair of sequences
	    
	    if(i+number_of_TTP_being_train>TTPsize)
	    {
		number_of_TTP_being_train = TTPsize-i;
	    }

	    for(j=0;j<number_of_TTP_being_train;j++)
	    {
		trancount[j] = 0;
	    }
    	    
	    // initialize
	    for(j=0; j<number_of_TTP_being_train; j++)
	    {	
		for(int k=0; k<number_of_states; k++)
		{
		    num_of_from_list[j][k] = 0;	   		  
		    for(int l=0; l<number_of_states; l++)
		    {
			if(from_list[j][k][l] == -1)
			{
			    break;
			}
			from_list[j][k][l] = -1;
		    }
		}
	    }
	    
	    for(j=0; j<number_of_TTP_being_train; j++)
	    {
		index = TP->get_TTP_to(j+i);
		from_list[j][index][num_of_from_list[j][index]] = TP->get_TTP_from(j+i);
		num_of_from_list[j][index]++;
	    }
	   	 	   
	    check+= Viterbi_train_transition(sX,sY,
					     number_of_TTP_being_train,
					     num_of_from_list,
					     from_list,
					     &trancount,
					     tTube,
					     Start_point,
					     End_point);
	    
	    if(check)
	    {
		cout<<"Warning: Hmm class:: error in Viterbi_train_transition of the "
		    <<SeqCount<<"-th sequence : Viterbi_train_TTP";		
		for(j=i; j<i+number_of_TTP_being_train; j++)
		{
		    cout<<TP->get_TTP_from(j)<<"->"<<TP->get_TTP_to(j)<<endl;
		}
	    }
	    
	    // set the counts of all sequences to TTP[i].score ( a_i,j)
	    	    
	    for(j=0; j<number_of_TTP_being_train; j++)
	    {
		total_TTP_count[j+i]+=trancount[j];
	    }	
	    
	    i+= number_of_TTP_being_train;	    	    	 	    

	}// loop over the TTP
	for(i=0; i<TTPsize; i++)
	{
	    (*prev_TTP_count)[SeqCount][i] = (*cur_TTP_count)[SeqCount][i];
	    (*cur_TTP_count)[SeqCount][i] = total_TTP_count[i];
	}
	SeqCount++;	
    } // after loop of sequence		          
    fclose(fIn_sequence);
	    
    // calculate the final score for each TTP set to be trained
    // i.e. a_i,j / b_i,j
    

    int fromindex = -1;
    // calculate the ratio from the counts
    long sumcount = 0;
    for(i=0; i<TTPsize ; i++)
    {		  
	sumcount = 0;
	fromindex = TP->get_TTP_from(i);		      		    
	j = 0;	
	while(j<TTPsize)
	{
	    if(TP->get_TTP_from(j)==fromindex)
	    {	
		sumcount+= total_TTP_count[j];

		if(j+1<TTPsize)
		{
		    if(TP->get_TTP_from(j+1)!=fromindex)
		    {
			break;
		    }
		}

	    }
	    j++;
	}
	
	tmp_TTP_score[i] = static_cast<double>(static_cast<double>(total_TTP_count[i])/static_cast<double>(sumcount));
    }
    
    // Add pseudoprob, now tmp_TTP_score store the probability
    for(i=0; i<TTPsize; i++)
    {		    
	tmp_TTP_score[i]+= TP->get_TTP_pseudoprob(i);
    }			
    // normalizing the prob again, fix transition to end

    for(i=0; i<TTPsize ; i++)
    {		  
	bij = 0;
	fromindex = TP->get_TTP_from(i);		      		    
	j = 0;	
	while(j<TTPsize)
	{
	    if(TP->get_TTP_from(j)==fromindex)
	    {	
		bij+=tmp_TTP_score[j];	
		if(j+1<TTPsize)
		{
		    if(TP->get_TTP_from(j+1)!=fromindex)
		    {
			break;
		    }
		}

	    }
	    j++;
	}	
	double to_end_prob = model[number_of_states-1].get_transition_prob(fromindex);
	check+=TP->set_TTP_score(i,(tmp_TTP_score[i]/bij)*(1-to_end_prob));	
    }
 
    //-----------------------------------------------------------------------------
    // End : Training of TTP in this loop
    //-----------------------------------------------------------------------------

    // free memory

    if(total_TTP_count) delete [] total_TTP_count;
    total_TTP_count = NULL;

    if(trancount) delete [] trancount;
    trancount = NULL;

    if(num_of_from_list) 
    {
	for(i=0; i<TTPsize; i++)
	{
	    if(num_of_from_list[i]) delete[] num_of_from_list[i];
	    num_of_from_list[i] = NULL;
	}
	delete [] num_of_from_list;
    }
    num_of_from_list = NULL;

    if(from_list)
    {
	for(i=0; i<TTPsize; i++)
	{
	    if(from_list[i])
	    {
		for(j=0; j<number_of_states; j++)
		{
		    if(from_list[i][j]) delete [] from_list[i][j];
		    from_list[i][j] = NULL;
		}
		delete [] from_list[i];
	    }
	    from_list[i] = NULL;	   
	}
	delete [] from_list;
    }
    from_list = NULL;   

    if(tmp_TTP_score) delete[] tmp_TTP_score;
    tmp_TTP_score = NULL;
  
    return check;
}

int Hmm::Posterior_train_TTP(const char* seqfile,
			     const char* seqannfile,
			     const char* name_input_tube_file,
			     const long int max_volume, 
			     const int radius,
			     int& SeqCount,
			     model_parameters* const MP,
			     TransitionProb* TP,
			     const int SamplePaths)
{


    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: Posterior_train_TTP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
    
    
    int i = 0;
    int j = 0;

    int TTPsize = TP->get_TTPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;   
     
    long* trancount = new long[TTPsize];
    double* tmp_TTP_score = new double[TTPsize];
    double bij = Logzero;
    
    int** num_of_from_list = new int*[TTPsize];
    int*** from_list = new int**[TTPsize];
    for(i=0; i<TTPsize; i++)
    {
	num_of_from_list[i] = new int[number_of_states];      
	from_list[i] = new int*[number_of_states];
	for(j=0; j<number_of_states; j++)
	{
	    num_of_from_list[i][j] = 0;
	    from_list[i][j] = new int[number_of_states];
	    for(int k=0; k<number_of_states; k++)
	    {
		from_list[i][j][k] = -1;
	    }
	}		
    }
    long* total_TTP_count = new long[TTPsize];

    int index = 0;
    int number_of_TTP_being_train = 0;
    int number_of_sample_paths_being_train = 0;
    
    // ---------------------------------------------------------------------
    // Start: training of TTP
    // ---------------------------------------------------------------------

    // for each pair of sequences		    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: Posterior_train_TTP: "
	    <<"training sequence file : "<<seqfile
	    <<" can not be read, training terminated."<<endl;
	check++;
    }
       	    
    // initialization
    SeqCount = 0;    
    for(j=0;j<TTPsize;j++)
    {
	total_TTP_count[j] = 0;
    }	
    
    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);

	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point);	
	
	int length_y = sY.length();
	
	int tube_length = number_of_states*(length_y+2);
	    
	int tube_area = tube_width*tube_length*8;
	int total_volume = 2*tube_area;
		
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: train_transitions : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" is "<<total_volume<<endl;	       
	    break;
	}

	number_of_sample_paths_being_train = 0;
	// check whether all the TTP can be train with one sample path
	if(total_volume+(TTPsize-1)*tube_area<=max_volume) 
	{
	    // train all TTP with a sample path
	    number_of_TTP_being_train = TTPsize;
	    total_volume += (TTPsize-1)*tube_area;
	    tube_area = total_volume - tube_area;
	    // calculate number of sample paths can be sampled in one iteration
	    while(number_of_sample_paths_being_train<SamplePaths)
	    {
		total_volume += tube_area;
		if(total_volume <= max_volume)
		{
		    number_of_sample_paths_being_train++;
		}else{
		    break;
		}
	    }
	}else{
	    // can not train all TTP with a sample path
	    number_of_sample_paths_being_train =1;
	    // calculate number of TTP can be train with one sample path
	    number_of_TTP_being_train = 1;
	    while(number_of_TTP_being_train<TTPsize)
	    {
		total_volume += tube_area;
		if(total_volume <= max_volume)
		{
		    number_of_TTP_being_train++;
		}else{
		    break;
		}
	    }
	}

	cout<<"Number of transition parameters being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "
	    <<number_of_TTP_being_train<<endl;
	cout<<"Number of Sample paths being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "
	    <<number_of_sample_paths_being_train<<endl;

	int s = 0;
	while(s<SamplePaths) // for each Sample path
	{	    
	    i = 0;
	    while(i<TTPsize)  // for each TTP
	    {		    
		// not enough memory to train all TTP with a sample path
		// calculate number of TTP can be trained for this pair of sequences	    

		if(i+number_of_TTP_being_train>TTPsize)
		{
		    number_of_TTP_being_train = TTPsize-i;
		}

		// initialize
		for(j=0; j<number_of_TTP_being_train; j++)
		{	
		    for(int k=0; k<number_of_states; k++)
		    {
			num_of_from_list[j][k] = 0;	   		  
			for(int l=0; l<number_of_states; l++)
			{
			    if(from_list[j][k][l] == -1)
			    {
				break;
			    }
			    from_list[j][k][l] = -1;
			}
		    }
		}
		
		for(j=0; j<number_of_TTP_being_train; j++)
		{
		    index = TP->get_TTP_to(j+i);
		    from_list[j][index][num_of_from_list[j][index]] = TP->get_TTP_from(j+i);
		    num_of_from_list[j][index]++;
		    trancount[j] = 0;
		}
		
		check+= Posterior_train_transition(sX,sY,
						   number_of_TTP_being_train,
						   number_of_sample_paths_being_train,
						   num_of_from_list,
						   from_list,
						   &trancount,
						   tTube,
						   Start_point,
						   End_point);	
		    
		if(check)
		{
		    cout<<"Warning: Hmm class:: error in Posterior_train_transition of the "
			<<SeqCount<<"-th sequence : Posterior_train_TTP";		
		    for(j=i; j<i+number_of_TTP_being_train; j++)
		    {
			cout<<TP->get_TTP_from(j)<<"->"<<TP->get_TTP_to(j)<<endl;
		    }
		}
		
		// set the counts of all sequences to TTP[i].score ( a_i,j)
	    	
		for(j=0; j<number_of_TTP_being_train; j++)
		{
		    total_TTP_count[j+i]+=trancount[j];
		}	
				
		i+= number_of_TTP_being_train;
		
	    }// loop over the TTP
	    s += number_of_sample_paths_being_train;
	}//loop over number of sample paths
	SeqCount++;	
    } // after loop of sequence		          
    fclose(fIn_sequence);
      
    // calculate the final score for each TTP set to be trained
    // i.e. a_i,j / b_i,j   

    int fromindex = -1;
    // calculate the ratio from the counts
    long sumcount = 0;
    for(i=0; i<TTPsize ; i++)
    {		  
	sumcount = 0;
	fromindex = TP->get_TTP_from(i);		      		    
	j = 0;	
	while(j<TTPsize)
	{
	    if(TP->get_TTP_from(j)==fromindex)
	    {	
		sumcount+= total_TTP_count[j];

		if(j+1<TTPsize)
		{
		    if(TP->get_TTP_from(j+1)!=fromindex)
		    {
			break;
		    }
		}

	    }
	    j++;
	}
	
	tmp_TTP_score[i] = static_cast<double>(static_cast<double>(total_TTP_count[i])/static_cast<double>(sumcount));
    }
    
    // Add pseudoprob, now tmp_TTP_score store the probability
    for(i=0; i<TTPsize; i++)
    {		    
	tmp_TTP_score[i]+= TP->get_TTP_pseudoprob(i);
    }			
    // normalizing the prob again

    for(i=0; i<TTPsize ; i++)
    {		  
	bij = 0;
	fromindex = TP->get_TTP_from(i);		      		    
	j = 0;	
	while(j<TTPsize)
	{
	    if(TP->get_TTP_from(j)==fromindex)
	    {	
		bij+=tmp_TTP_score[j];	
		if(j+1<TTPsize)
		{
		    if(TP->get_TTP_from(j+1)!=fromindex)
		    {
			break;
		    }
		}

	    }
	    j++;
	}	
	double to_end_prob = model[number_of_states-1].get_transition_prob(fromindex);
	check+=TP->set_TTP_score(i,(tmp_TTP_score[i]/bij)*(1-to_end_prob));	
    }
    
    // put the TTP back to Original TP for this iteration		

    //-----------------------------------------------------------------------------
    // End : Training of TTP in this loop
    //-----------------------------------------------------------------------------

    // free memory

    if(trancount) delete [] trancount;
    trancount = NULL;

    if(total_TTP_count) delete [] total_TTP_count;
    total_TTP_count = NULL;

    if(num_of_from_list) 
    {
	for(i=0; i<TTPsize; i++)
	{
	    if(num_of_from_list[i]) delete[] num_of_from_list[i];
	    num_of_from_list[i] = NULL;
	}
	delete [] num_of_from_list;
    }
    num_of_from_list = NULL;

    if(from_list)
    {
	for(i=0; i<TTPsize; i++)
	{
	    if(from_list[i])
	    {
		for(j=0; j<number_of_states; j++)
		{
		    if(from_list[i][j]) delete [] from_list[i][j];
		    from_list[i][j] = NULL;
		}
		delete [] from_list[i];
	    }
	    from_list[i] = NULL;	   
	}
	delete [] from_list;
    }
    from_list = NULL;   

    if(tmp_TTP_score) delete[] tmp_TTP_score;
    tmp_TTP_score = NULL;
  
    return check;
}

int Hmm::BaumWelch_train_TTP(const char* seqfile,
			     const char* seqannfile,
			     const char* name_input_tube_file,
			     const long int max_volume,
			     const int radius,
			     bool& get_forward_score,
			     double& cur_forward_score,
			     int& SeqCount,
			     model_parameters* const MP,
			     TransitionProb* TP)
{
    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: BaumWelch_train_TTP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
    
    
    int i = 0;
    int j = 0;

    int TTPsize = TP->get_TTPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;   
	    
    double* aij = new double[TTPsize];
    double bij=Logzero;
    double* sum_aij = new double[TTPsize];
    double sum_bij = Logzero;
    double* tmp_aij = new double[TTPsize];
    double tmp_bij = Logzero;
    double* tmp_TTP_score = new double[TTPsize];
    double forward_score = Logzero;
    int** num_of_from_list = new int*[TTPsize];
    int*** from_list = new int**[TTPsize];
    for(i=0; i<TTPsize; i++)
    {
	aij[i] = Logzero;
	sum_aij[i] = Logzero;
	tmp_aij[i] = Logzero;
	num_of_from_list[i] = new int[number_of_states];      
	from_list[i] = new int*[number_of_states];
	for(j=0; j<number_of_states; j++)
	{
	    num_of_from_list[i][j] = 0;
	    from_list[i][j] = new int[number_of_states];
	    for(int k=0; k<number_of_states; k++)
	    {
		from_list[i][j][k] = -1;
	    }
	}		
    }

    int index = 0;
    int number_of_TTP_being_train = 0;
    
    
    // ---------------------------------------------------------------------
    // Start: training of TTP
    // ---------------------------------------------------------------------
    
    // for each pair of sequences		    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: BaumWelch_train_TTP: "
	    <<"training sequence file : "<<seqfile
	    <<" can not be read, training terminated."<<endl;
	check++;
    }
    
    // train transition_scores(to,from)		    
    SeqCount = 0;
    for(i=0;i<TTPsize;i++)
    {
	aij[i] = Logzero;
    }

    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);
	for(i=0; i<TTPsize; i++)
	{
	    tmp_aij[i] = Logzero;
	}
	
	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point);	
	
	int length_y = sY.length();
	
	int tube_length = number_of_states*(length_y+2);
	
	int tube_area = tube_width*tube_length;
	int total_volume = 2*tube_area*8;
	
	
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: BaumWelch_train_TTP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" is "<<total_volume<<endl;	       
	    break;
	}

	// calculate max number of TTP can be trained in one train_transition function

	number_of_TTP_being_train = 1;

	while(number_of_TTP_being_train<TTPsize)
	{		
	    total_volume += tube_area;
	    if(total_volume <= max_volume)
	    {
		number_of_TTP_being_train++;
	    }else{
		break;
	    }
	}
     
	cout<<"Number of transition probabilities being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "<<number_of_TTP_being_train<<endl;

	// loop over TTP
	i =0 ;
	while(i<TTPsize)
	{	
	    // calculate number of TTP can be trained in for this pair of sequences
	    
	    if(i+number_of_TTP_being_train>=TTPsize)
	    {
		number_of_TTP_being_train = TTPsize-i;
	    }	    	    	    	 	    
	    // initialize
	    for(j=0; j<number_of_TTP_being_train; j++)
	    {	
		for(int k=0; k<number_of_states; k++)
		{
		    num_of_from_list[j][k] = 0;	   		  
		    for(int l=0; l<number_of_states; l++)
		    {
			if(from_list[j][k][l] == -1)
			{
			    break;
			}
			from_list[j][k][l] = -1;
		    }
		}
	    }
	    
	    for(j=0; j<number_of_TTP_being_train; j++)
	    {
		index = TP->get_TTP_to(j+i);
		from_list[j][index][num_of_from_list[j][index]] = TP->get_TTP_from(j+i);
		num_of_from_list[j][index]++;
	    }
	   	 	   
	    check+= BaumWelch_train_transition(sX,sY,
					       number_of_TTP_being_train,
					       num_of_from_list,
					       from_list,
					       &tmp_aij,
					       forward_score,
					       tTube,
					       Start_point,
					       End_point);
	    
	    if(check)
	    {
		cout<<"Warning: Hmm class:: error in BaumWelch_train_transition of the "
		    <<SeqCount<<"-th sequence : BaumWelch_train_TTP";		
		for(j=i; j<i+number_of_TTP_being_train; j++)
		{
		    cout<<TP->get_TTP_from(j)<<"->"<<TP->get_TTP_to(j)<<endl;
		}
	    }
	    
	    for(j=i; j<i+number_of_TTP_being_train; j++)
	    {
		if(aij[j]<=Logzero)
		{
		    aij[j] = tmp_aij[j];
		    sum_aij[j] = 0;
		}else{
		    if(tmp_aij[j] > Logzero)
		    {
			sum_aij[j] += exp(tmp_aij[j] - aij[j]);
		    }
		}  		
	    }      
	    
	    // get forward score
	    if(!get_forward_score)
	    {
		if(cur_forward_score <= Logzero)
		{
		    if(forward_score > Logzero)
		    {
			cur_forward_score = forward_score;
			get_forward_score = true;
		    }
		}else{
		    if(forward_score > Logzero)
		    {
			cur_forward_score += forward_score;
		    }				      	     
		}		
	    }
	    
	    i+= number_of_TTP_being_train;
	    
	}	
	SeqCount++;
	
    } // after loop of sequence		   
	    	
    fclose(fIn_sequence);
   
    // set the counts of all sequences to TTP[i].score ( a_i,j)
    for(i=0; i<TTPsize; i++)
    {
	sum_aij[i] = log(1+sum_aij[i]);
	aij[i] += sum_aij[i];					 
	
	TP->set_TTP_score(i,aij[i]);
	
	//the forward score of current sequence is got
	
    }// loop over the TTP
   
    // calculate the final score for each TTP set to be trained
    // i.e. a_i,j / b_i,j
    

    int fromindex = -1;
    
    for(i=0; i<TTPsize ; i++)
    {		  
	bij = Logzero;
	sum_bij = Logzero;
	fromindex = TP->get_TTP_from(i);		      		    
	j = 0;	
	while(j<TTPsize)
	{
	    if((TP->get_TTP_from(j)==fromindex)&&(TP->get_TTP_to(j)!=number_of_states-1))
	    {			
		if(bij<=Logzero)
		{
		    bij = TP->get_TTP_score(j);
		    sum_bij = 0;		       
		}else
		{
		    sum_bij = sum_bij + exp(TP->get_TTP_score(j)-bij);
		}
		
		if(j+1<TTPsize)
		{
		    if(TP->get_TTP_from(j+1)!=fromindex)
		    {
			break;
		    }
		}

	    }
	    j++;
	}
	
	// get a_i,j / b_i,j and put in TTP[i].score
	
	sum_bij = log(1+sum_bij);
	bij += sum_bij;
	
	tmp_TTP_score[i] = TP->get_TTP_score(i) - bij; // aij := TPP[i].score
	
    }
    
    // Add pseudoprob, now tmp_TTP_score store the probability
    for(i=0; i<TTPsize; i++)
    {		    
	tmp_TTP_score[i] = exp(tmp_TTP_score[i])+TP->get_TTP_pseudoprob(i);
    }			
    // normalizing the prob again, fix transition to end  

    for(i=0; i<TTPsize ; i++)
    {		  
	bij = 0;
	fromindex = TP->get_TTP_from(i);		      		    
	j = 0;	
	while(j<TTPsize)
	{
	    if(TP->get_TTP_from(j)==fromindex)
	    {	
		bij+=tmp_TTP_score[j];	
		if(j+1<TTPsize)
		{
		    if(TP->get_TTP_from(j+1)!=fromindex)
		    {
			break;
		    }
		}

	    }
	    j++;
	}	
	double to_end_prob = model[number_of_states-1].get_transition_prob(fromindex);
	check+=TP->set_TTP_score(i,(tmp_TTP_score[i]/bij)*(1-to_end_prob));	
    }

    //-----------------------------------------------------------------------------
    // End : Training of TTP in this loop
    //-----------------------------------------------------------------------------

    // free memory
    
    if(aij) delete [] aij;
    aij = NULL;

    if(sum_aij) delete [] sum_aij;
    sum_aij = NULL;
    
    if(tmp_aij) delete [] tmp_aij;
    tmp_aij = NULL;

    if(num_of_from_list) 
    {
	for(i=0; i<TTPsize; i++)
	{
	    if(num_of_from_list[i]) delete[] num_of_from_list[i];
	    num_of_from_list[i] = NULL;
	}
	delete [] num_of_from_list;
    }
    num_of_from_list = NULL;

    if(from_list)
    {
	for(i=0; i<TTPsize; i++)
	{
	    if(from_list[i])
	    {
		for(j=0; j<number_of_states; j++)
		{
		    if(from_list[i][j]) delete [] from_list[i][j];
		    from_list[i][j] = NULL;
		}
		delete [] from_list[i];
	    }
	    from_list[i] = NULL;	   
	}
	delete [] from_list;
    }
    from_list = NULL;   

    if(tmp_TTP_score) delete[] tmp_TTP_score;
    tmp_TTP_score = NULL;
  
    return check;
}

int Hmm::BaumWelch_train_FTP(const char* seqfile,
			     const char* seqannfile,
			     const char* name_input_tube_file,
			     const long int max_volume,
			     const int radius,
			     bool& get_forward_score,
			     double& cur_forward_score,
			     int& SeqCount,
			     model_parameters* const MP,
			     TransitionProb* TP)
{
    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: BaumWelch_train_FTP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
    
    
    int i = 0;
    int j = 0;

     // store both the counts of from and over-from probabilities of a GTP
    int GTPsize = 2*TP->get_GTPsize(); 
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;   
	    
    double* aij = new double[GTPsize];
    double bij=Logzero;
    double* sum_aij = new double[GTPsize];
    double sum_bij = Logzero;
    double* tmp_aij = new double[GTPsize];
    double tmp_bij = Logzero;
    double forward_score = Logzero;
    int** num_of_from_list = new int*[GTPsize];
    int*** from_list = new int**[GTPsize];
    for(i=0; i<GTPsize; i++)
    {
	aij[i] = Logzero;
	sum_aij[i] = Logzero;
	tmp_aij[i] = Logzero;
	num_of_from_list[i] = new int[number_of_states];      
	from_list[i] = new int*[number_of_states];
	for(j=0; j<number_of_states; j++)
	{
	    num_of_from_list[i][j] = 0;
	    from_list[i][j] = new int[number_of_states];
	    for(int k=0; k<number_of_states; k++)
	    {
		from_list[i][j][k] = -1;
	    }
	}		
    }

    int index = 0;
    int number_of_GTP_being_train = 0;
        
    // ---------------------------------------------------------------------
    // Start: training of FTP, train all the GTP
    // ---------------------------------------------------------------------
    	
    // for each pair of sequences		    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: BaumWelch_train_FTP: "
	    <<"training sequence file : "<<seqfile
	    <<" can not be read, training terminated."<<endl;
	check++;
    }
    
    // train transition_scores(to,from)		    
    SeqCount = 0;

    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);
	for(i=0; i<GTPsize; i++)
	{
	    tmp_aij[i] = Logzero;
	}
	
	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point);	
	
	int length_y = sY.length();
	
	int tube_length = number_of_states*(length_y+2);
	
	int tube_area = tube_width*tube_length;
	int total_volume = 2*tube_area*8;
	
	
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: BaumWelch_train_FTP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" is "<<total_volume<<endl;	       
	    break;
	}

	// calculate max number of GTP can be trained in one train_transition function

	number_of_GTP_being_train = 1;

	while(number_of_GTP_being_train<GTPsize)
	{		
	    total_volume += tube_area;
	    if(total_volume <= max_volume)
	    {
		number_of_GTP_being_train++;
	    }else{
		break;
	    }
	}
	
	cout<<"Number of group transition parameters being trained for the "<<SeqCount
	    <<"-th input sequnce in parallel : "<<number_of_GTP_being_train<<endl;
	
	// loop over GTP
	i =0 ;
	while(i<GTPsize)
	{		
	    // calculate number of GTP can be trained in for this pair of sequences
	    
	    if(i+number_of_GTP_being_train>=GTPsize)
	    {
		number_of_GTP_being_train = GTPsize-i;
	    }	    	    	    	 	    
	    // initialize
	    for(j=0; j<number_of_GTP_being_train; j++)
	    {	
		for(int k=0; k<number_of_states; k++)
		{
		    num_of_from_list[j][k] = 0;	   		  
		    for(int l=0; l<number_of_states; l++)
		    {
			if(from_list[j][k][l] == -1)
			{
			    break;
			}
			from_list[j][k][l] = -1;
		    }
		}
	    }
	    
	    for(j=0; j<number_of_GTP_being_train; j++)
	    {
		if(j%2==0)
		{
		    for(int k=0; k<TP->get_GTP_NumOffromto((j+i)/2); k++)
		    {
			index = TP->get_GTP_to((j+i)/2,k);
			from_list[j][index][num_of_from_list[j][index]] = TP->get_GTP_from((j+i)/2,k);
			num_of_from_list[j][index]++;
		    }
		}else{
		    for(int k=0; k<TP->get_GTP_NumOfOverfromto((j+i)/2); k++)
		    {
			index = TP->get_GTP_Overto((j+i)/2,k);
			from_list[j][index][num_of_from_list[j][index]] = TP->get_GTP_Overfrom((j+i)/2,k);
			num_of_from_list[j][index]++;
		    }
		}
	    }
	    	   	 	   
	    check+= BaumWelch_train_transition(sX,sY,
					       number_of_GTP_being_train,
					       num_of_from_list,
					       from_list,
					       &tmp_aij,
					       forward_score,
					       tTube,
					       Start_point,
					       End_point);	    
	    if(check)
	    {
		cout<<"Warning: Hmm class:: error in BaumWelch_train_transition of the "
		    <<SeqCount<<"-th sequence : train_FTP";		
		for(j=i; j<i+number_of_GTP_being_train; j+=2)
		{
		    cout<<"FTP: "<<(j/2)<<endl;
		}
	    }	
	    
	    // tmp_aij[j] store the counts of the GTP j for this sequence
	    
	    for(j=i; j<i+number_of_GTP_being_train; j++)
	    {
		if(aij[j]<=Logzero)
		{
		    aij[j] = tmp_aij[j];
		    sum_aij[j] = 0;
		}else{
		    if(tmp_aij[j] > Logzero)
		    {
			sum_aij[j] += exp(tmp_aij[j] - aij[j]); // sum the aij for each sequence
		    }
		}  		
	    }      
	    
	    // get forward score
	    if(!get_forward_score)
	    {
		if(cur_forward_score <= Logzero)
		{
		    if(forward_score > Logzero)
		    {
			cur_forward_score = forward_score;
			get_forward_score = true;
		    }
		}else{
		    if(forward_score > Logzero)
		    {
			cur_forward_score += forward_score;
		    }				      	     
		}		
	    }
	    
	    i+= number_of_GTP_being_train;
	    
	}	
	SeqCount++;
	
    } // after loop of sequence		   
	    	
    fclose(fIn_sequence);
   
    // calculate the final counts for each GTP
    for(i=0; i<GTPsize; i++)
    {
	sum_aij[i] = log(1+sum_aij[i]);
	aij[i] += sum_aij[i];					 
    }

    // normalization 
    // the i-th GTP score is calculated by exp(aij[2*i] - aij[2*i+1]); 
    
    for(i=0; i<TP->get_GTPsize() ; i++)
    {		
	TP->set_GTP_score(i,exp(aij[2*i]-aij[2*i+1]));
    }

    // derive the FTPs from the GTPs
    check += TP->derive_FTPs_from_GTPs();
    if(check)
    {
	cout<<"ERROR: Hmm class:: BaumWelch_train_FTP : "
	    <<"error in derive_FTPs_from_GTPs"<<endl;      
    }

    // add psuedocount for FTPs
    for(i=0; i<TP->get_GTPsize(); i++)
    {
     	TP->set_FTP_prob(i,TP->get_FTP_prob(i)+ TP->get_FTP_pseudocount(i));
    }

    //-----------------------------------------------------------------------------
    // End : Training of GTP in this loop
    //-----------------------------------------------------------------------------

    // free memory
    
    if(aij) delete [] aij;
    aij = NULL;

    if(sum_aij) delete [] sum_aij;
    sum_aij = NULL;
    
    if(tmp_aij) delete [] tmp_aij;
    tmp_aij = NULL;

    if(num_of_from_list) 
    {
	for(i=0; i<GTPsize; i++)
	{
	    if(num_of_from_list[i]) delete[] num_of_from_list[i];
	    num_of_from_list[i] = NULL;
	}
	delete [] num_of_from_list;
    }
    num_of_from_list = NULL;

    if(from_list)
    {
	for(i=0; i<GTPsize; i++)
	{
	    if(from_list[i])
	    {
		for(j=0; j<number_of_states; j++)
		{
		    if(from_list[i][j]) delete [] from_list[i][j];
		    from_list[i][j] = NULL;
		}
		delete [] from_list[i];
	    }
	    from_list[i] = NULL;	   
	}
	delete [] from_list;
    }
    from_list = NULL;   
  
    return check;
}

int Hmm::Viterbi_train_EP(const char* seqfile,
			  const char* seqannfile,
			  const char* name_input_tube_file,
			  const long int max_volume,
			  const int radius,
			  int& SeqCount,			  
			  long**** prev_EP_count,
			  long**** cur_EP_count,
			  model_parameters* const MP,
			  EmissionProb* EP)
{
    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: Viterbi_train_EP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }

    int i = 0;
    int j = 0;

    int FEPsize = EP->get_FEPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }

    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;
    
    double bij=Logzero; 
    int index = -1;
    long number_of_emission = 0;
    long max_number_of_parameters = 0;
    bool needtrain;

    // calculate max_number_of_parameters
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index>=0)
	    {
		number_of_emission = static_cast<long>(pow(alphabet,model[EP->get_GEP_from(index,0)].get_emission_prob().GetNumberofDimensions()));
		for(j=0; j<number_of_emission; j++)
		{
		    needtrain = false;
		    for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		    {		  		
			if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
			   ||(EP->get_FEP_pseudoprob(i,j)!=0))
			{
			    needtrain = true;
			    break;
			}						
		    }		    
		    if(needtrain)
		    {			   
			max_number_of_parameters++;
		    }
		}		
	    }
	}
    }
       
    int** emit_from_list = new int*[max_number_of_parameters];
    for(i=0; i<max_number_of_parameters; i++)
    {
	emit_from_list[i] = new int[number_of_states];
    }
    array<Prob>* tmp_FEP_prob = new array<Prob>[FEPsize];

    int* FEP_index = new int[max_number_of_parameters];
    int* GEP_index = new int[max_number_of_parameters];
    long* emit_index = new long[max_number_of_parameters];
    long* emission_index = new long[max_number_of_parameters];
    long* emit_count = new long[max_number_of_parameters];
    long* EP_count = new long[max_number_of_parameters];
    		
    // initialize variables
    
    long number_of_emissions_being_train = 0;   
    long parameters_count = 0;

    for(i=0; i<FEPsize; i++)
    {
	tmp_FEP_prob[i] = EP->get_FEP_prob(i);
    }

    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index>=0)
	    {
		number_of_emission = static_cast<long>(pow(alphabet,model[EP->get_GEP_from(index,0)].get_emission_prob().GetNumberofDimensions()));
		for(j=0; j<number_of_emission; j++)
		{
		    needtrain = false;
		    for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		    {		  		
			if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
			   ||(EP->get_FEP_pseudoprob(i,j)!=0))
			{
			    needtrain = true;
			    break;
			}							
		    }		    
		    if(needtrain)
		    {			
			FEP_index[parameters_count] = i;
			emit_index[parameters_count] = j;		
			GEP_index[parameters_count] = index;
			EP_count[parameters_count] = 0;
			parameters_count++;
		    }		
		}		
	    }
	}
    }

    // ---------------------------------------------------------------------
    // Start: training of EP
    // ---------------------------------------------------------------------    
   
    SeqCount = 0;

    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: Viterbi_train_EP: "
	    <<"training sequence file : "<<seqfile
	    <<" does not exist"<<endl;
	check++;
    }		
    
    // for each sequences
    

    while((!feof(fIn_sequence))&&(!check))
    {			
	//number_of_emissions_being_train = 0;					    
	
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);

	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point); 
	
	int length_y = sY.length();		       
	long long tube_length = number_of_states*(length_y+1);
	long long tube_area = tube_width*tube_length;
	long long total_volume = 2*tube_area*8;
	long long volume_for_one_emission = total_volume;
	
	if(volume_for_one_emission>max_volume)
	{
	    cout<<"ERROR:: Viterbi_train_EP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" : "<<volume_for_one_emission<<endl;	       
	    break;
	}

	// calculate the maximum number of parameters can be trained in 1 train_emission function
	number_of_emissions_being_train = 1;
	while(number_of_emissions_being_train<max_number_of_parameters)
	{
	    total_volume += tube_area;
	    if(total_volume <= max_volume)
	    {
		number_of_emissions_being_train++;
	    }else{
		break;
	    }
	}

	cout<<"Number of emission parameters being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel: "<<number_of_emissions_being_train<<endl;

	i =0 ;
	while(i<max_number_of_parameters)
	{
	    if(i+number_of_emissions_being_train > max_number_of_parameters)
	    {
		number_of_emissions_being_train = max_number_of_parameters-i;
	    }

	    for(j=0; j<number_of_emissions_being_train; j++)
	    {		    	
		for(int k =0; k<number_of_states; k++)
		{
		    emit_from_list[j][k] = -1;
		}
		for(int k =0; k<EP->get_GEP_NumOffrom(GEP_index[i+j]); k++)
		{		  
		    emit_from_list[j][EP->get_GEP_from(GEP_index[i+j],k)] = EP->get_GEP_from(GEP_index[i+j],k);				
		}
		emit_count[j] = 0;
		emission_index[j] = emit_index[i+j];
	    }
	    	  	
	    check+= Viterbi_train_emission(sX,sY,
					   number_of_emissions_being_train,
					   emission_index,
					   emit_from_list,
					   &emit_count,
					   tTube,
					   Start_point,
					   End_point);
	    
	    if(check)
	    {
		cout<<"Warning: Hmm class:: error in Viterbi_train_emission of the "
		    <<SeqCount<<"-th sequence : Viterbi_train_EP: ";
		for(j = i; j<i+number_of_emissions_being_train; j++)
		{
		    cout<<emission_index[j]<<" ";
		}
		cout<<"of FEP("<<FEP_index[i]<<")"<<endl;    
	    }

	    // set the counts of all sequences to 
	    
	    for(j = 0; j<number_of_emissions_being_train; j++)
	    {
		EP_count[i+j] += emit_count[j];	
	    }			
	    
	    i+= number_of_emissions_being_train;		 

	}// loop over i-th parameter
	SeqCount++;
	
    } // after loop of sequence		   
  
    fclose(fIn_sequence);
    
    int paramcount = 0;
    int prevparamcount = 0;


    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index < 0)
	    {
		cout<<"ERROR: class Hmm:: Viterbi_train_EP : "
		    <<"GEP id of the "<<i<<"-th FEP is not valid : "
		    <<EP->get_FEP_exp(i)<<endl;
		check++;
		break;
	    }	

	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));
	    long sumcount = 0;
	    prevparamcount = paramcount;
	    for(j=0; j<number_of_emission; j++)
	    {
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;
			break;
		    }						
		}		    
		if(needtrain)
		{			   
		    sumcount+= EP_count[paramcount];		
		    paramcount++;
		}
	    }	
	    
	    paramcount = prevparamcount;
	    for(j=0; j<number_of_emission; j++)
	    {
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;
			break;
		    }						
		}		    
		if(needtrain)
		{
		    tmp_FEP_prob[i].SetElement(j,
					       static_cast<double>(EP_count[paramcount]/
								   static_cast<double>(sumcount)) 
					       + EP->get_FEP_pseudoprob(i,j));	
		    paramcount++;
		}else{
		    tmp_FEP_prob[i].SetElement(j,0.0);
		}
	    }	   
	}
    }


    // put the tmp_GEP back to FEP and normalization
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));
	    bij = 0;
	    for(j=0; j<number_of_emission; j++)		
	    {
		bij += tmp_FEP_prob[i].GetElement(j);	
	    }
	    for(j=0; j<number_of_emission; j++)
	    {	    
		EP->set_FEP_prob(i,j,tmp_FEP_prob[i].GetElement(j)/bij);
	    }
	}
    }

    // free memory

    if(FEP_index) delete [] FEP_index;
    FEP_index = NULL;
    if(GEP_index) delete [] GEP_index;
    GEP_index = NULL;
    if(emit_index) delete [] emit_index;
    emit_index = NULL;
    if(emit_count) delete [] emit_count;
    emit_count = NULL;
    if(EP_count) delete [] EP_count;
    EP_count = NULL;
     
    if(tmp_FEP_prob) delete [] tmp_FEP_prob;
    tmp_FEP_prob = NULL;
    
    if(emit_from_list)
    {
	for(i=0; i<max_number_of_parameters; i++)
	{
	    if(emit_from_list[i]) delete [] emit_from_list[i];
	    emit_from_list[i] = NULL;
	}
	delete[] emit_from_list;
    }
    emit_from_list = NULL;

    return check;
}

int Hmm::Posterior_train_EP(const char* seqfile,
			    const char* seqannfile,
			    const char* name_input_tube_file,
			    const long int max_volume,
			    const int radius,
			    int& SeqCount,				 
			    model_parameters* const MP,
			    EmissionProb* EP,
			    const int SamplePaths)
{
    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: Posterior_train_EP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }

    int i = 0;
    int j = 0;

    int FEPsize = EP->get_FEPsize();
    int max_d = 0;
	
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }

    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;
    
    double bij=Logzero; 
    int index = -1;
    long max_number_of_parameters = 0;
    long number_of_emission = 0;
    bool needtrain;

    // calculate max_number_of_parameters 
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index>=0)
	    {
		number_of_emission = static_cast<long>(pow(alphabet,model[EP->get_GEP_from(index,0)].get_emission_prob().GetNumberofDimensions()));
		for(j=0; j<number_of_emission; j++)
		{
		    needtrain = false;
		    for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		    {		  		
			if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
			   ||(EP->get_FEP_pseudoprob(i,j)!=0))
			{
			    needtrain = true;
			    break;
			}						
		    }		    
		    if(needtrain)
		    {			   
			max_number_of_parameters++;
		    }
		}		
	    }
	}
    }

    int** emit_from_list = new int*[max_number_of_parameters];
    for(i=0; i<max_number_of_parameters; i++)
    {
	emit_from_list[i] = new int[number_of_states];
    }
    array<Prob>* tmp_FEP_prob = new array<Prob>[FEPsize];
    for(i=0; i<FEPsize; i++)
    {
	tmp_FEP_prob[i] = EP->get_FEP_prob(i);
    }

    int* FEP_index = new int[max_number_of_parameters];
    int* GEP_index = new int[max_number_of_parameters];
    long* emit_index = new long[max_number_of_parameters];
    long* emission_index = new long[max_number_of_parameters];
    long* emit_count = new long[max_number_of_parameters];
    long* EP_count = new long[max_number_of_parameters];

    // initialize variables
    
    int number_of_emissions_being_train = 0;
    int number_of_sample_paths_being_train = 0;
    int parameters_count = 0;
    
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index>=0)
	    {
		number_of_emission = static_cast<long>(pow(alphabet,model[EP->get_GEP_from(index,0)].get_emission_prob().GetNumberofDimensions()));
		for(j=0; j<number_of_emission; j++)
		{
		    needtrain = false;
		    for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		    {		  		
			if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
			   ||(EP->get_FEP_pseudoprob(i,j)!=0))
			{
			    needtrain = true;
			    break;
			}							
		    }		    
		    if(needtrain)
		    {			
			FEP_index[parameters_count] = i;
			emit_index[parameters_count] = j;		
			GEP_index[parameters_count] = index;
			EP_count[parameters_count] = 0;
			parameters_count++;
		    }		
		}		
	    }
	}
    }

    // ---------------------------------------------------------------------
    // Start: training of EP
    // ---------------------------------------------------------------------    
   
    SeqCount = 0;

    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: Posterior_train_EP: "
	    <<"training sequence file : "<<seqfile
	    <<" does not exist"<<endl;
	check++;
    }		
    
    // for each sequences
    
    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);

	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point); 
	
	long long int length_y = sY.length();		       
	long long int tube_length = number_of_states*(length_y+1);
	long long int tube_area = tube_width*tube_length;
	long long int total_volume = 2*tube_area*8;
	
	if(total_volume>max_volume)
	{
	    cout<<"ERROR:: Posterior_train_EP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" : "<<total_volume<<endl;	       
	    break;
	}
	
	number_of_sample_paths_being_train = 0;
	// calculate the maximum nunmber of parameters can be trained in 1 train_emission function
	if(total_volume+(max_number_of_parameters-1)*tube_area<=max_volume)
	{
	    // train all emissions in one iteration
	    number_of_emissions_being_train = max_number_of_parameters;
	    total_volume += (max_number_of_parameters-1)*tube_area;
	    tube_area = total_volume - tube_area;
	    while(number_of_sample_paths_being_train<SamplePaths)
	    {
		total_volume += tube_area;
		if(total_volume <= max_volume)
		{
		    number_of_sample_paths_being_train++;
		}else{
		    break;
		}
	    }
	}else{
	    // can not train all emissions with a sample path
	    number_of_sample_paths_being_train = 1;
	    number_of_emissions_being_train = 1;
	    while(number_of_emissions_being_train<max_number_of_parameters)
	    {
		total_volume+= tube_area;
		if(total_volume <= max_volume)
		{
		    number_of_emissions_being_train++;
		}else{
		    break;
		}
	    }
	}
	
	cout<<"Number of emission parameters being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "
	    <<number_of_emissions_being_train<<endl;
	cout<<"Number of Sample paths being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "
	    <<number_of_sample_paths_being_train<<endl;

	int s = 0;
	while(s<SamplePaths)
	{
	    i=0;

	    while(i<max_number_of_parameters)
	    {
		// not enough memory to train all emissinos with a sample path
		// calculate number of emissions can be trained for this pair of sequences	    
		if(i+number_of_emissions_being_train > max_number_of_parameters)
		{
		    number_of_emissions_being_train = max_number_of_parameters - i;
		}
			
		for(j=0; j<number_of_emissions_being_train; j++)
		{	
		    for(int k = 0; k<number_of_states; k++)
		    {
			emit_from_list[j][k] = -1;
		    }
		    for(int k =0; k<EP->get_GEP_NumOffrom(GEP_index[i+j]); k++)
		    {		  
			emit_from_list[j][EP->get_GEP_from(GEP_index[i+j],k)] = EP->get_GEP_from(GEP_index[i+j],k);		
		    }
		    emit_count[j] = 0;
		    emission_index[j] = emit_index[i+j];
		}

		check+= Posterior_train_emission(sX,sY,
						 number_of_emissions_being_train,
						 number_of_sample_paths_being_train,
						 emission_index,
						 emit_from_list,
						 &emit_count,				
						 tTube,
						 Start_point,
						 End_point);
		
		if(check)
		{
		    cout<<"Warning: Hmm class:: error in Posterior_train_emission of the "
			<<SeqCount<<"-th sequence : Posterior_train_EP: ";
		    for(j = i; j<i+number_of_emissions_being_train; j++)
		    {
			cout<<emission_index[j]<<" ";
		    }
		    cout<<"of FEP("<<FEP_index[j]<<")"<<endl;    
		}

		// set the counts of all sequences to 
		
		for(j = 0; j<number_of_emissions_being_train; j++)
		{
		    EP_count[i+j] += emit_count[j];	
		}			
		
		i+= number_of_emissions_being_train;

	    }
	    s += number_of_sample_paths_being_train;
	}
	
	SeqCount++;
	
    } // after loop of sequence		   

    fclose(fIn_sequence);
    // convert counts to prob
    
    int paramcount = 0;
    int prevparamcount = 0;
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index < 0)
	    {
		cout<<"ERROR: class Hmm:: Posterior_train_EP : "
		    <<"GEP id of the "<<i<<"-th FEP is not valid : "
		    <<EP->get_FEP_exp(i)<<endl;
		check++;
		break;
	    }	

	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));

	    long sumcount = 0;
	    prevparamcount = paramcount;
	    for(j=0; j<number_of_emission; j++)		
	    {	
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;
			break;
		    }						
		}
		if(needtrain)
		{
		    sumcount+= EP_count[paramcount];
		    paramcount++;
		}
	    }

	    paramcount = prevparamcount;
	    for(j=0; j<number_of_emission; j++)
	    {		
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;
			break;
		    }						
		}	
		if(needtrain)
		{
		    tmp_FEP_prob[i].SetElement(j,
					       static_cast<double>(EP_count[paramcount]/
								   static_cast<double>(sumcount)) 
					       + EP->get_FEP_pseudoprob(i,j));
		    paramcount++;
		}else{
		    tmp_FEP_prob[i].SetElement(j,0.0);
		}
	    }	   
	}
    }
    
    // put the tmp_GEP back to FEP and normalization
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));
	    bij = 0;
	    for(j=0; j<number_of_emission; j++)		
	    {
		 bij += tmp_FEP_prob[i].GetElement(j);	
	    }
	    for(j=0; j<number_of_emission; j++)
	    {	    
		EP->set_FEP_prob(i,j,tmp_FEP_prob[i].GetElement(j)/bij);
	    }
	}
    }
    // free memory

    if(FEP_index) delete [] FEP_index;
    FEP_index = NULL;
    if(GEP_index) delete [] GEP_index;
    GEP_index = NULL;
    if(emit_index) delete [] emit_index;
    emit_index = NULL;
    if(emit_count) delete [] emit_count;
    emit_count = NULL;
    if(EP_count) delete [] EP_count;
    EP_count = NULL;
     
    if(tmp_FEP_prob) delete [] tmp_FEP_prob;
    tmp_FEP_prob = NULL;
    
    if(emit_from_list)
    {
	for(i=0; i<max_number_of_parameters; i++)
	{
	    if(emit_from_list[i]) delete [] emit_from_list[i];
	    emit_from_list[i] = NULL;
	}
	delete[] emit_from_list;
    }
    emit_from_list = NULL;
    
    return check;
}

int Hmm::BaumWelch_train_EP(const char* seqfile,
			    const char* seqannfile,
			    const char* name_input_tube_file,
			    const long int max_volume,
			    const int radius,
			    bool& get_forward_score,
			    double& cur_forward_score,
			    int& SeqCount,
			    model_parameters* const MP,
			    EmissionProb* EP)
{
    int check = 0;
    if(!seqfile)
    {
	cout<<"ERROR:: class hmm: BaumWelch_train_EP: seqfile is NULL!"<<endl;
	check++;
	return check;
    }
    
    int i = 0;
    int j = 0;
    
    int FEPsize = EP->get_FEPsize();
    int GEPsize = EP->get_GEPsize();
    int max_d = 0;
    
    for (i = 1; i < (number_of_states-1); i++) 
    {
	if (this->model[i].get_letters_to_read_x() > max_d) {max_d = this->model[i].get_letters_to_read_x();}
    }
    
    const int max_deltax = max_d;
    const int tube_width = max_deltax+1;
    
    double bij=Logzero; 
    double sum_bij = Logzero;
    double tmp_bij = Logzero;
    double forward_score = Logzero;
    int index = -1;
    int number_of_emission = 0;
    int max_number_of_parameters = 0;
    bool needtrain;

    // calculate max_number_of_parameters
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index>=0)
	    {
		number_of_emission = static_cast<long>(pow(alphabet,model[EP->get_GEP_from(index,0)].get_emission_prob().GetNumberofDimensions()));
		for(j=0; j<number_of_emission; j++)
		{
		    needtrain = false;
		    for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		    {		  		
			if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
			   ||(EP->get_FEP_pseudoprob(i,j)!=0))
			{
			    needtrain = true;
			    break;
			}						
		    }		    
		    if(needtrain)
		    {			   
			max_number_of_parameters++;
		    }
		}		
	    }
	}
    }

    
    int** emit_from_list = new int*[max_number_of_parameters];
    for(i=0; i<max_number_of_parameters; i++)
    {
	emit_from_list[i] = new int[number_of_states];
    }
    array<Prob>* tmp_FEP_prob = new array<Prob>[FEPsize];
    
    int* FEP_index = new int[max_number_of_parameters];
    int* GEP_index = new int[max_number_of_parameters];
    int* emit_index = new int[max_number_of_parameters];
    int* emission_index = new int[max_number_of_parameters];
    double* aij = new double[max_number_of_parameters];
    double* sum_aij = new double[max_number_of_parameters];
    double* tmp_aij = new double[max_number_of_parameters];

    // initialize variables
    
    array<Prob>* tmp_GEP = new array<Prob>[GEPsize];

    long number_of_emissions_being_train = 0;
    long parameters_count = 0;

    for(i=0; i<FEPsize; i++)
    {
	tmp_FEP_prob[i] = EP->get_FEP_prob(i);
    }

    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index>=0)
	    {
		number_of_emission = static_cast<long>(pow(alphabet,model[EP->get_GEP_from(index,0)].get_emission_prob().GetNumberofDimensions()));
		for(j=0; j<number_of_emission; j++)
		{
		    needtrain = false;
		    for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		    {		  		
			if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
			   ||(EP->get_FEP_pseudoprob(i,j)!=0))
			{
			    needtrain = true;
			    break;
			}							
		    }		    
		    if(needtrain)
		    {			
			FEP_index[parameters_count] = i;
			emit_index[parameters_count] = j;		
			GEP_index[parameters_count] = index;
			aij[parameters_count] = Logzero;
			sum_aij[parameters_count] = Logzero;
			parameters_count++;
		    }		
		}		
	    }
	}
    }
    
    for(i=0; i<FEPsize; i++)
    {	
	if(EP->is_FEP_train(i))
	{
	    // tmp_GEP for temporary storage of the score
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index < 0)
	    {
		cout<<"ERROR: class Hmm:: BaumWelch_train_EP : "
		    <<"GEP id of the "<<i<<"-th FEP is not valid : "
		    <<EP->get_FEP_exp(i)<<endl;
		check++;
		break;
	    }
	    tmp_GEP[index] = model[EP->get_GEP_from(index,0)].get_emission_prob();
	  	    
	}
    }
    
    SeqCount = 0;
    
    FILE* fIn_sequence = fopen(seqfile,"rt");
    
    if(!fIn_sequence)
    {
	cout<<"Error:: in hmm class: BaumWelch_train_EP: "
	    <<"training sequence file : "<<seqfile
	    <<" does not exist"<<endl;
	check++;
    }		
    
    // for each sequences
    
    
    while((!feof(fIn_sequence))&&(!check))
    {			
	Sequence sX, sY;					
	Tube<int> tTube;
	vector<int> Start_point(3), End_point(3);
	
	check+=get_sequence_and_tube_for_training(fIn_sequence,
						  seqannfile,
						  name_input_tube_file,
						  radius,
						  MP,
						  &sX, &sY,
						  &tTube,
						  &Start_point,
						  &End_point); 

	int length_y = sY.length();		       
	long long int tube_length = number_of_states*(length_y+1);
	long long int tube_area = tube_width*tube_length;
	long long int total_volume = 2*tube_area*8;
	long long int volume_for_one_emission = total_volume;
	
	if(volume_for_one_emission>max_volume)
	{
	    cout<<"ERROR:: BaumWelch_train_EP : "
		<<"Maximum volume allowed : "<<max_volume
		<<" not enough "<<endl;
	    cout<<"total memory needed for training "
		<<"Seq : "<<sX.get_ac();
	    if(length_y>0)
	    {
		cout<<" and Seq : "<<sY.get_ac();
	    }
	    cout<<" : "<<volume_for_one_emission<<endl;	       
	    break;
	}
	
	// calculate the maximum number of parameters can be trained in 1 train_emission function
	number_of_emissions_being_train = 1;
	while(number_of_emissions_being_train<max_number_of_parameters)
	{
	    total_volume += tube_area;
	    if(total_volume <= max_volume)
	    {
		number_of_emissions_being_train++;
	    }else{
		break;
	    }
	}

	cout<<"Number of emission parameters being trained for the "<<SeqCount
	    <<"-th input sequence (pair) in parallel : "<<number_of_emissions_being_train<<endl;

	i =0 ;

	while(i<max_number_of_parameters)
	{
	    if(i+number_of_emissions_being_train > max_number_of_parameters)
	    {
		number_of_emissions_being_train = max_number_of_parameters-i;
	    }

	    for(j=0; j<number_of_emissions_being_train; j++)
	    {		    	
		for(int k =0; k<number_of_states; k++)
		{
		    emit_from_list[j][k] = -1;
		}
		for(int k =0; k<EP->get_GEP_NumOffrom(GEP_index[i+j]); k++)
		{		  
		    emit_from_list[j][EP->get_GEP_from(GEP_index[i+j],k)] = EP->get_GEP_from(GEP_index[i+j],k);				
		}
		aij[j+i] = Logzero;
		sum_aij[j+i] = Logzero;
		tmp_aij[j] = Logzero;
		emission_index[j] = emit_index[i+j];
	    }

	    // calculate 
	
	    check+= BaumWelch_train_emission(sX,sY,
					     number_of_emissions_being_train,
					     emission_index,
					     emit_from_list,
					     &tmp_aij,
					     forward_score,
					     tTube,
					     Start_point,
					     End_point);

		
	    if(check)
	    {
		cout<<"Warning: Hmm class:: error in BaumWelch_train_emission of the "
		    <<SeqCount<<"-th sequence : BaumWelch_train_EP: ";
		for(j = i; j<i+number_of_emissions_being_train; j++)
		{
		    cout<<emission_index[j]<<" ";
		}
		cout<<"of FEP("<<FEP_index[i]<<")"<<endl;    
	    }
			
	    for(j= 0; j<number_of_emissions_being_train; j++)
	    {
		if(aij[j+i]<=Logzero)
		{
		    aij[i+j] = tmp_aij[j];
		    sum_aij[i+j] = 0;
		}else{
		    if(tmp_aij[j] > Logzero)
		    {
			sum_aij[i+j] += exp(tmp_aij[j] - aij[i+j]);
		    }
		}  		      
	    }
	    // get forward score
	    if(!get_forward_score)
	    {
		if(cur_forward_score <= Logzero)
		{
		    if(forward_score > Logzero)
		    {
			cur_forward_score = forward_score;
			get_forward_score = true;
		    }
		}else{
		    if(forward_score > Logzero)
		    {
			cur_forward_score += forward_score;
		    }				      	     
		}				    
	    }			  		
	    i+=number_of_emissions_being_train;
	} // loop over number of parameters		
	SeqCount++;	
    } // after loop of sequence		   
    
    fclose(fIn_sequence);
    
    // sum up the counts for all FEPs and all sequences 

    long paramcount = 0;
    long prevparamcount = 0;
    
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{

	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));	   
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index < 0)
	    {
		cout<<"ERROR: class Hmm:: train_emission : "
		    <<"GEP id of the "<<i<<"-th FEP is not valid : "
		    <<EP->get_FEP_exp(i)<<endl;
		check++;
		break;
	    }	    
	    
	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));	   
	    for(j = 0; j<number_of_emission; j++)
	    {
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;			
		    }			
		}
		if(needtrain)
		{
		    if(aij[paramcount]>Logzero)
		    {
			sum_aij[paramcount] = log(1+sum_aij[paramcount]);
			aij[paramcount] += sum_aij[paramcount];			  
		    }
		    tmp_GEP[index].SetElement(j,aij[paramcount]); // some may set to be Logzero		
		    paramcount++;
		}
	    }	
	}
    }

    
    // calculate the final score for each GEP set to be trained
    // i.e. a_i,j / b_i,j
    
    paramcount = 0;
    for(i=0; i<FEPsize; i++)
    {	
	if(EP->is_FEP_train(i))
	{
	    index = convert_GEP_id_to_int(EP,EP->get_FEP_exp(i));
	    if(index < 0)
	    {
		cout<<"ERROR: class Hmm:: BaumWelch_train_EP : "
		    <<"GEP id of the "<<i<<"-th FEP is not valid : "
		    <<EP->get_FEP_exp(i)<<endl;
		check++;
		break;
	    }
	    
	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));
	    
	    bij = Logzero;
	    sum_bij = Logzero;
	    
	    for(j=0; j<number_of_emission ; j++)
	    {			       
		
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;			
			break;
		    }			
		}

	    }
	    if(!needtrain)
	    {
		continue;
	    }
	    
	    for(j=0; j<number_of_emission; j++)
	    {
		if(bij<=Logzero)
		{		       
		    bij = tmp_GEP[index].GetElement(j);
		    sum_bij = 0;		       
		}else
		{
		    if(tmp_GEP[index].GetElement(j)>Logzero)
		    {
			sum_bij = sum_bij + exp(tmp_GEP[index].GetElement(j)-bij);
		    }
		}		
	    }		
	    // get a_i,j / b_i,j and put in EP[i].score
	    if(sum_bij>Logzero)
	    {
		sum_bij = log(1+sum_bij);
		bij += sum_bij; // the sum of the counts of all emission probs in j-th GEP
	    }
	    
	    for(j=0; j<number_of_emission; j++)
	    {		
		needtrain = false;
		for(int k=0; k<EP->get_GEP_NumOffrom(index); k++)
		{		  		
		    if((model[EP->get_GEP_from(index,k)].get_emission_prob(j)!=0)
		       ||(EP->get_FEP_pseudoprob(i,j)!=0))
		    {
			needtrain = true;			
			break;
		    }			
		}
		if(needtrain)
		{
		    aij[paramcount] = tmp_GEP[index].GetElement(j);		   
		    if(aij[paramcount]>Logzero)
		    {
			tmp_GEP[index].SetElement(j,exp(aij[paramcount]-bij)+EP->get_FEP_pseudoprob(i,j));  
		    }else{
			tmp_GEP[index].SetElement(j,EP->get_FEP_pseudoprob(i,j));
		    }
		}
	    }
	    tmp_FEP_prob[i] = tmp_GEP[index] ;
	}
    }
    
    // put the tmp_GEP back to FEP and normalization
    for(i=0; i<FEPsize; i++)
    {
	if(EP->is_FEP_train(i))
	{
	    number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(i).GetNumberofDimensions()));
	    bij = 0;
	    for(j=0; j<number_of_emission; j++)		
	    {
		bij += tmp_FEP_prob[i].GetElement(j);	
	    }
	    for(j=0; j<number_of_emission; j++)
	    {	    
		EP->set_FEP_prob(i,j,tmp_FEP_prob[i].GetElement(j)/bij);
	    }
	}
    }
    
    
    // free memory
    
    if(aij) delete [] aij;
    aij = NULL;
    
    if(tmp_aij) delete [] tmp_aij;
    tmp_aij = NULL;
   
    if(sum_aij) delete [] sum_aij;
    sum_aij = NULL; 

    if(FEP_index) delete [] FEP_index;
    FEP_index = NULL;

    if(GEP_index) delete [] GEP_index;
    GEP_index = NULL;
       
    if(emit_index) delete [] emit_index;
    emit_index = NULL;

    if(emission_index) delete [] emission_index;
    emission_index = NULL;
        
    if(tmp_GEP) delete[] tmp_GEP;
    tmp_GEP = NULL;    

    if(tmp_FEP_prob) delete [] tmp_FEP_prob;
    tmp_FEP_prob = NULL;
    
    if(emit_from_list) 
    {
	for(i=0; i<max_number_of_parameters; i++)
	{
	    if(emit_from_list[i]) delete [] emit_from_list[i];
	    emit_from_list[i] = NULL;
	}
	delete[] emit_from_list;
    }
    emit_from_list = NULL;
       
    return check;
}

int Hmm::ViterbiTraining(const char* XMLfile,
			 TransitionProb* TP,
			 EmissionProb* EP,
			 model_parameters* const MP,
			 const char* seqfile,
			 const char* seqannfile,
			 const char* name_input_tube_file,
			 const long long int max_volume,
			 const int radius,
			 const int Maxiter)
{
    int check = 0;
    
    if(!XMLfile)
    {
	cout<<"ERROR:: class Hmm: ViterbiTraining: input XMLfile is NULL!"<<endl;
	check++;
    }

    if(!check)
    {	
	
	int checktp = 0;
	int checkep = 0;

	int i=0;
	int j=0;
	
	checktp += TP->get_parameters_for_training(XMLfile,MP,this);
	if(checktp)
	{
	    cout<<"ERROR: Hmm:: ViterbiTraining: error in getting TP parameters for training";
	}
	
	checkep += EP->get_parameters_for_training(XMLfile,MP,this);
	if(checkep)
	{
	    cout<<"ERROR: Hmm:: ViterbiTraining: error in getting EP parameters for training";
	}
	
	check += checktp + checkep;
	
	if(!check)
	{
	    // start the training here
	    
	    // allocate memory for temporary storage of set of transition and emission probs
	    // emission: always declared in the EmissionProb Class
	    // transition: two cases:
	    //             - stored in TransitionProb Class
	    //             - no TransitionProb Class created

	    int index = -1;
	    int iCount = 0;	    	    
	    const int max_num_of_seq = 5000;
	    int SeqCount = 0;

	    int GTPsize = TP->get_GTPsize();
	    int FTPsize = GTPsize*2;
	    
	    long** prev_FTP_count = NULL;
	    long** cur_FTP_count = NULL;
	    bool trainFTP = false;
	    
	    if(GTPsize>0)
	    {
		for(i=0; i<FTPsize; i++)
		{
		    if(TP->is_FTP_train(i))
		    {
			trainFTP = true;
			break;
		    }
		}
	    }
	    if(trainFTP)
	    {
		prev_FTP_count = new long*[max_num_of_seq];
		cur_FTP_count = new long*[max_num_of_seq];
		
		for(i=0; i<max_num_of_seq; i++)
		{
		    prev_FTP_count[i] = new long[FTPsize];
		    cur_FTP_count[i] = new long[FTPsize];
		    for(j=0; j<FTPsize; j++)		    
		    {
			prev_FTP_count[i][j] = 0;
			cur_FTP_count[i][j] = 0;
		    }
		}
	    }
	    int FEPsize = EP->get_FEPsize();	    
	    long*** prev_EP_count = NULL;
	    long*** cur_EP_count = NULL;
	    bool trainEP = false;
	    for(i=0; i<FEPsize; i++)
	    {
		if(EP->is_FEP_train(i))
		{
		    trainEP = true;
		    break;
		}
	    }
	    if(trainEP)
	    {
		prev_EP_count = new long**[max_num_of_seq];
		cur_EP_count = new long**[max_num_of_seq];
	    
		for(i=0; i<max_num_of_seq; i++)
		{
		    prev_EP_count[i] = new long*[FEPsize];
		    cur_EP_count[i] = new long*[FEPsize];
		    
		    for(j=0; j<FEPsize; j++)		    
		    {
			int number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(j).GetNumberofDimensions()));
			prev_EP_count[i][j] = new long[number_of_emission];
			cur_EP_count[i][j] = new long[number_of_emission];
			for(int k=0; k<number_of_emission; k++)
			{
			    prev_EP_count[i][j][k] = 0;
			    cur_EP_count[i][j][k] = 0;
			}
		    }
		}
	    }

	    int TTPsize = TP->get_TTPsize();
	    long** prev_TTP_count = NULL;
	    long** cur_TTP_count = NULL;
	    bool trainTTP = false;
	    if(TTPsize>0)
	    {
		trainTTP = true;
	    }
	    if(trainTTP)
	    {
		prev_TTP_count = new long*[max_num_of_seq];
		cur_TTP_count = new long*[max_num_of_seq];
		
		for(i=0; i<max_num_of_seq; i++)
		{
		    prev_TTP_count[i] = new long[TTPsize];
		    cur_TTP_count[i] = new long[TTPsize];
		    for(j=0; j<TTPsize; j++)		    
		    {
			prev_TTP_count[i][j] = 0;
			cur_TTP_count[i][j] = 0;
		    }
		}
	    }

	    bool change = trainFTP|trainEP|trainTTP;
	    bool FTPchange = false;
	    bool tranchange = false;
	    bool emitchange = false;
	    
	    srand(time(NULL));

	    // read the superfastexptable for sample paths
	    //readsuperfastexptable();

	    cout<<"Start Viterbi training"<<endl;
	    time_t t = 0;
	    
	    while((change)&&(iCount<Maxiter))
	    {
		change = false;
		FTPchange = false;
		tranchange = false;
		emitchange = false;
		
		// train FTP
		cout<<"iteration : "<<iCount<<endl;							    

		if(trainFTP)
		{
		    time_t t1 = time(NULL);
			
		    check+= Viterbi_train_FTP(seqfile,
					      seqannfile,
					      name_input_tube_file,
					      max_volume,
					      radius,
					      SeqCount,
					      &prev_FTP_count,
					      &cur_FTP_count,
					      MP,
					      TP);
		    
		    time_t t2 = time(NULL);
		    
		    if(check)
		    {
		        cout<<"ERROR:: hmm class: ViterbiTraining : "
			    <<"error in Viterbi_train_FTP "<<endl;
			break;
		    }else{			
			t += t2-t1;
		        cout<<endl<<"time for training free transition parameters: "<<t2-t1<<"s"<<endl<<endl; 				
		    }
			
		    for(i=0; i<SeqCount;i++)
		    {
		        for(j=0; j<FTPsize; j++)		    
			{
			    if(TP->is_FTP_train(j))
			    {
			        if(prev_FTP_count[i][j]!=cur_FTP_count[i][j])
				{
				    FTPchange = true;
				    break;
				}
			    }
			}
		    }

		}
		    
		if(trainTTP)
		{
		    time_t t1=time(NULL);				    		    
		    
		    check+= Viterbi_train_TTP(seqfile,
					      seqannfile,
					      name_input_tube_file,
					      max_volume,
					      radius,
					      SeqCount,
					      &prev_TTP_count,
					      &cur_TTP_count,
					      MP,
					      TP);
			
		    time_t t2 = time(NULL);
			
		    if(check)
		    {
		        cout<<"ERROR:: hmm class: ViterbiTraining : "
			    <<"error in Viterbi_train_TTP "<<endl;
			break;
		    }else{			
				t+= t2-t1;
		        cout<<endl<<"time for training transition probabilities: "<<t2-t1<<"s"<<endl<<endl; 
		    }

		    for(i=0; i<SeqCount;i++)
		    {
		        for(j=0; j<TTPsize; j++)
			{
			    if(prev_TTP_count[i][j]!=cur_TTP_count[i][j])
			    {
				tranchange = true;
				break;
			    }
			}
		    }	 
		    
		}

		if(trainEP)
		{
		    time_t t1 = time(NULL);
		    check+= Viterbi_train_EP(seqfile,
					     seqannfile,
					     name_input_tube_file,
					     max_volume,
					     radius,
					     SeqCount,
					     &prev_EP_count,
					     &cur_EP_count,
					     MP,
					     EP);
		    time_t t2 = time(NULL);

		    if(check)
		    {
		        cout<<"ERROR:: hmm class: ViterbiTraining : "
			    <<"error in Viterbi_train_EP "<<endl;
			break;
		    }else{			    		    				    
				t+=t2-t1;
		        cout<<endl<<"time for training emission parameters: "<<t2-t1<<"s"<<endl<<endl;
		    }

		    for(i=0; i<SeqCount;i++)
		    {
		        for(j=0; j<FEPsize; j++)
			{
			    if(EP->is_FEP_train(j))
			    {
				int number_of_emission = static_cast<long>(pow(alphabet,EP->get_FEP_prob(j).GetNumberofDimensions()));
				for(int k=0;k<number_of_emission;k++)
				{
				    if(prev_EP_count[i][j][k]!=cur_EP_count[i][j][k])
				    {
				        emitchange = true;
					break;
				    }
				}
			    }
			    if(emitchange)
			    {
			        break;
			    }
			}
		    }			
		}
		change = tranchange|FTPchange|emitchange;

		if(trainFTP)
		{
			check+= derive_transition_probs_from_transition_parameters(this,MP,TP);
				for(i=0; i<TP->get_GTPsize(); i++)
    		{
      		cout<<"FTP_prob("<<i<<") : "<<TP->get_FTP_prob(i)<<endl;
    		}
		    if(check)
		    {
			cout<<"Error: hmm class:: ViterbiTraining : "
			    <<"error in derive_transition_probs_from_transition_parameters()"<<endl;
		    }	
		}
       
		if(trainTTP)
		{
			for(i=0; i<TTPsize; i++)
			{
				this->model[TP->get_TTP_from(i)].set_transition_prob(TP->get_TTP_to(i),TP->get_TTP_score(i));
			}
		}

		if(trainEP)
		{
			check+= derive_emission_probs(this,MP,EP);
			if(check)
			{
				cout<<"Error: hmm class:: ViterbiTraining : " 
					<<"error in training emission prob: derive_emission_probs()"<<endl;
			}
		}

		// print transition probs		
		if((trainTTP)||(trainFTP))
		{
			for(i=0; i<number_of_states-1; i++)
			{
				cout<<"From ["<<i<<"] to "<<endl;
				model[i].print_transition_probs(cout);
			}
		}

		//print emission
		if(trainEP)
		{
			for(i=1; i<number_of_states-1; i++)
			{
				cout<<"State : "<<i<<endl;
				model[i].print_emission_probs(cout);
			}
		}
		//modify the score array after an interation
		check += this->calculate_scores_from_probs();
		if(check)
		{
		    cout<<"Error: hmm class:: ViterbiTraining : "
			<<"error in  calculate_scores_from_probs() "
			<<"in loop("<<iCount<<endl;
		}	  

		check += this->check_consistency_of_probs();

		if(check)
		{
		    cout<<"Error: hmm class:: ViterbiTraining : "
			<<"error in check_consistency_of_probs() "
			<<"in loop("<<iCount<<endl;
		}	
		
		iCount++;
	    } // finished training

		cout<<"Total time used for the Viterbi training : "<<t<<" s."<<endl;
	    
	    // free memeory

	    if(prev_TTP_count)
	    {
		for(i=0; i<max_num_of_seq; i++)
		{
		    if(prev_TTP_count[i]) delete [] prev_TTP_count[i];
		    prev_TTP_count[i] = NULL;
		}
		delete [] prev_TTP_count;
	    }
	    prev_TTP_count = NULL;
	    
	    if(cur_TTP_count) 
	    {
		for(i=0; i<max_num_of_seq; i++)
		{
		    if(cur_TTP_count[i]) delete [] cur_TTP_count[i];
		    cur_TTP_count[i] = NULL;
		}
		delete [] cur_TTP_count;
	    }
	    cur_TTP_count = NULL;


	    if(prev_EP_count)
	    {
		for(i=0; i<max_num_of_seq; i++)
		{
		    if(prev_EP_count[i])
		    { 
			for(j=0; j<FEPsize; j++)
			{
			    if(prev_EP_count[i][j]) delete [] prev_EP_count[i][j];
			    prev_EP_count[i][j] = NULL;
			}
			delete [] prev_EP_count[i];
		    }
		    prev_EP_count[i] = NULL;
		}
		delete [] prev_EP_count;
	    }
	    prev_EP_count = NULL;
	    
	    if(cur_EP_count) 
	    {
		for(i=0; i<max_num_of_seq; i++)
		{
		    if(cur_EP_count[i]) 
		    {
			for(j=0; j<FEPsize; j++)
			{
			    if(cur_EP_count[i][j]) delete [] cur_EP_count[i][j];
			    cur_EP_count[i][j] = NULL;
			}
			delete [] cur_EP_count[i];
		    }
		    cur_EP_count[i] = NULL;
		}
		delete [] cur_EP_count;
	    }
	    cur_EP_count = NULL;	    
	}	
    }
    return check;
}

int Hmm::PosteriorTraining(const char* XMLfile,
			   TransitionProb* TP,
			   EmissionProb* EP,
			   model_parameters* const MP,
			   const char* seqfile,
			   const char* seqannfile,
			   const char* name_input_tube_file,
			   const long long int max_volume,
			   const int radius,
			   const int Maxiter,
			   const int SamplePaths) 
{
    int check = 0;
    
    if(!XMLfile)
    {
	cout<<"ERROR:: class Hmm: PosteriorTraining: input XMLfile is NULL!"<<endl;
	check++;
    }

    if(!check)
    {	
	
	int checktp = 0;
	int checkep = 0;

	int i=0;
	int j=0;
	
	checktp += TP->get_parameters_for_training(XMLfile,MP,this);
	if(checktp)
	{
	    cout<<"ERROR: Hmm:: PosteriorTraining: error in getting TP parameters for training";
	}
	
	checkep += EP->get_parameters_for_training(XMLfile,MP,this);
	if(checkep)
	{
	    cout<<"ERROR: Hmm:: PosteriorTraining: error in getting EP parameters for training";
	}
	
	check += checktp + checkep;
	
	if(!check)
	{
	    // start the training here
	    
	    // allocate memory for temporary storage of set of transition and emission probs
	    // emission: always declared in the EmissionProb Class
	    // transition: two cases:
	    //             - stored in TransitionProb Class
	    //             - no TransitionProb Class created

	    int index = -1;
	    int iCount = 0;	    	    
	    const int max_num_of_seq = 5000;
	    int SeqCount = 0;

	    int GTPsize = TP->get_GTPsize();
	    int FTPsize = TP->get_FTPsize();
	    bool trainFTP = false;
	    
	    if(GTPsize>0)
	    {
		for(i=0; i<FTPsize; i++)
		{
		    if(TP->is_FTP_train(i))
		    {
			trainFTP = true;
			break;
		    }
		}
	    }

	    int FEPsize = EP->get_FEPsize();	    
	    bool trainEP = false;
	    for(i=0; i<FEPsize; i++)
	    {
		if(EP->is_FEP_train(i))
		{
		    trainEP = true;
		    break;
		}
	    }

	    int TTPsize = TP->get_TTPsize();
	    bool trainTTP = false;
	    if(TTPsize>0)
	    {
		trainTTP = true;
	    }
	    
	    srand(time(NULL));

		cout<<"Start Posterior training "<<endl;
		time_t t = 0;

	    // read the superfastexptable for sample paths
	    //readsuperfastexptable();
	    
	    while(iCount<Maxiter)
	    {
		
		// train FTP
		cout<<"iteration : "<<iCount<<endl;				

		if(trainFTP)
		{
		    time_t t1 = time(NULL);

		    check+= Posterior_train_FTP(seqfile,
						seqannfile,
						name_input_tube_file,
						max_volume,
						radius,
						SeqCount,
						MP,
						TP,
						SamplePaths);
		    time_t t2 = time(NULL);

		    if(check)
		    {
		        cout<<"ERROR:: hmm class: PosteriorTraining : "
			    <<"error in Posterior_train_FTP "<<endl;
			break;
		    }else{
				t+= t2-t1;
		        cout<<endl<<"time for training free transition parameters: "<<t2-t1<<"s"<<endl<<endl;
		    }

		}
		    
		if(trainTTP)
		{ 
		    time_t t1=time(NULL);				    		    
		    
		    check+= Posterior_train_TTP(seqfile,
						     seqannfile,
						     name_input_tube_file,
						     max_volume,
						     radius,
						     SeqCount,
						     MP,
						     TP,
						     SamplePaths);
			
		    time_t t2 = time(NULL);
		    
		    if(check)
		    {
		        cout<<"ERROR:: hmm class: PosteriorTraining : "
			    <<"error in Posterior_train_TTP "<<endl;
			    break;
		    }else{
				t+= t2-t1;
		        cout<<endl<<"time for training transition probabilities: "<<t2-t1<<"s"<<endl<<endl;
		    }
		}
		    
		if(trainEP)
		{
		    time_t t1=time(NULL);
		    check+= Posterior_train_EP(seqfile,
						    seqannfile,
						    name_input_tube_file,
						    max_volume,
						    radius,
						    SeqCount,
						    MP,
						    EP,
						    SamplePaths);

		    time_t t2 = time(NULL);
			
		    if(check)
		    {
		        cout<<"ERROR:: hmm class: PosteriorTraining : "
			    <<"error in Posterior_train_EP "<<endl;
			break;
		    }else{
			t+= t2-t1;
		        cout<<endl<<"time for training emission parameters: "<<t2-t1<<"s"<<endl<<endl;
		    }
		}		

		if(trainFTP)
		{
		    check+= derive_transition_probs_from_transition_parameters(this,MP,TP);
		    for(i=0; i<TP->get_GTPsize(); i++)
		    {
			cout<<"FTP_prob("<<i<<") : "<<TP->get_FTP_prob(i)<<endl;
		    }
		    if(check)
		    {
			cout<<"Error: hmm class:: PosteriorTraining : "
			    <<"error in derive_transition_probs_from_transition_parameters()"<<endl;
		    }	
		}
		
		if(trainTTP)
		{
		    for(i=0; i<TTPsize; i++)
		    {
			this->model[TP->get_TTP_from(i)].set_transition_prob(TP->get_TTP_to(i),TP->get_TTP_score(i));
		    }
		}
	
		if(trainEP)
		{
		    check+= derive_emission_probs(this,MP,EP);
		    if(check)
		    {
			cout<<"Error: hmm class:: PosteriorTraining : " 
			    <<"error in training emission prob: derive_emission_probs()"<<endl;
		    }
		}

		// print transition probs
		if((trainTTP)||(trainFTP))
		{
		    for(i=0; i<number_of_states-1; i++)
		    {
			cout<<"From ["<<i<<"] to "<<endl;
			model[i].print_transition_probs(cout);
		    }
		}

		
		//print emission
		if(trainEP)
		{
		    for(i=1; i<number_of_states-1; i++)
		    {
			cout<<"State : "<<i<<endl;
			model[i].print_emission_probs(cout);
		    }
		}

		//modify the score array after an interation
		check += this->calculate_scores_from_probs();
		if(check)
		{
		    cout<<"Error: hmm class:: PosteriorTraining : "
			<<"error in  calculate_scores_from_probs() "
			<<"in loop("<<iCount<<endl;
		}	  

		check += this->check_consistency_of_probs();

		if(check)
		{
		    cout<<"Error: hmm class:: PosteriorTraining : "
			<<"error in check_consistency_of_probs() "
			<<"in loop("<<iCount<<endl;
		}	
		
		iCount++;
	    } // finished training
	    cout<<"Total time used for the Posterior training : "<<t<<" s."<<endl;
	}	    
    }
    return check;
}

int Hmm::BaumWelchTraining(const char* XMLfile,
			   TransitionProb* TP,
			   EmissionProb* EP,
			   model_parameters* const MP,
			   const char* seqfile,
			   const char* seqannfile,
			   const char* name_input_tube_file,
			   const long long int max_volume,
			   const int radius,
			   const int Maxiter,
			   const double threshold) 
{
    int check = 0;
    
    if(!XMLfile)
    {
	cout<<"ERROR:: class Hmm: BaumWelchTraining: input XMLfile is NULL!"<<endl;
	check++;
    }
    
    if(!check)
    {	

	bool   get_forward_score = false;
	double prev_forward_score = Logzero;
	double cur_forward_score = Logzero;
	double forward_score = Logzero;
	
	double log_odd_ratio = 1000;	

	int checktp = 0;
	int checkep = 0;
       
	int i=0;
	int j=0;

	checktp += TP->get_parameters_for_training(XMLfile,MP,this);
	if(checktp)
	{
	    cout<<"ERROR: Hmm:: BaumWelchTraining: error in getting TP parameters for training";
	}
		
	checkep += EP->get_parameters_for_training(XMLfile,MP,this);
	if(checkep)
	{
	    cout<<"ERROR: Hmm:: BaumWelchTraining: error in getting EP parameters for trianing";
	}

	check += checktp + checkep;
	
	if(!check)
	{
	    // start the training here
	    
	    // allocate memory for temporary storage of set of transition and emission probs
	    // emission: always declared in the EmissionProb Class
	    // transition: two cases:
	    //             - stored in TransitionProb Class
	    //             - no TransitionProb Class created
	    
	    int index = -1;
	    int iCount = 0;
	    
	    int GTPsize = TP->get_GTPsize();
	    int FTPsize = TP->get_FTPsize();
	    bool trainFTP = false;
	    
	    if(GTPsize>0)
	    {
		for(i=0; i<FTPsize; i++)
		{
		    if(TP->is_FTP_train(i))
		    {
			trainFTP = true;
			break;
		    }
		}
	    }
	    
	    int FEPsize = EP->get_FEPsize();	    
	    bool trainEP = false;
	    for(i=0; i<FEPsize; i++)
	    {
		if(EP->is_FEP_train(i))
		{
		    trainEP = true;
		    break;
		}
	    }

	    int TTPsize = TP->get_TTPsize();
	    bool trainTTP = false;
	    if(TTPsize>0)
	    {
		trainTTP = true;
	    }

	    bool train = trainEP|trainFTP|trainTTP;
	    int SeqCount = 0;

	    //check if there is parameters set to be trained
	    cout<<"Start Baum-Welch training"<<endl;
	    
	    time_t t = 0;
	    
	    while((log_odd_ratio>threshold)&&(iCount<Maxiter)&&(train))
	    {
		get_forward_score = false;
		prev_forward_score = cur_forward_score;		
		cur_forward_score = Logzero;
				
		// train FTP
		cout<<"iteration : "<<iCount<<endl;		
		// train FTP		
		if(trainFTP)
		{
		    time_t t1=time(NULL);
		    
		    check+=BaumWelch_train_FTP(seqfile,
					       seqannfile,
					       name_input_tube_file,
					       max_volume,
					       radius,
					       get_forward_score,
					       cur_forward_score,
					       SeqCount,
					       MP,
					       TP);
		    time_t t2=time(NULL);

		    if(check)
		    {
			cout<<"ERROR:: hmm class: BaumWelchTraining : "
			    <<"error in BaumWelch_train_FTP "<<endl;
			break;
		    }else{
			cout<<endl<<"time for training free transition paramters: "<<t2-t1<<"s"<<endl<<endl;
			t += t2 - t1;
		    }		
		}
		
		// train TTP
		if(trainTTP)
		{
		    time_t t1 = time(NULL);

		    check+= BaumWelch_train_TTP(seqfile,
						seqannfile,
						name_input_tube_file,
						max_volume,
						radius,
						get_forward_score,
						cur_forward_score,
						SeqCount,
						MP,
						TP);
		
		    time_t t2 = time(NULL);
		
		    if(check)
		    {
			cout<<"ERROR:: hmm class: BaumWelchTraining : "
			    <<"error in BaumWelch_train_TTP "<<endl;
			break;
		    }else{
			 cout<<endl<<"time for training transition probabilities: "<<t2-t1<<"s"<<endl<<endl;	
			 t += t2 - t1;
		    }
		}
				
		// train EP
		
		if(trainEP)
		{
		    time_t t1 = time(NULL);
		    
		    check+= BaumWelch_train_EP(seqfile,
					       seqannfile,
					       name_input_tube_file,
					       max_volume,
					       radius,
					       get_forward_score,
					       cur_forward_score,
					       SeqCount,
					       MP,
					       EP);
		    time_t t2 = time(NULL);

		    if(check)
		    {
			cout<<"ERROR:: hmm class: BaumWelchTraining : "
			    <<"error in BaumWelch_train_EP "<<endl;
			break;
		    }else{		
			cout<<endl<<"time for training emission parameters: "<<t2-t1<<"s"<<endl<<endl;
			t += t2 - t1;
		    } 
		}

		// Put trained parameters back to original place
		
		if(trainFTP)
		{
		    check+= derive_transition_probs_from_transition_parameters(this,MP,TP);
		    for(i=0; i<TP->get_GTPsize(); i++)
    		{
      		cout<<"FTP_prob("<<i<<") : "<<TP->get_FTP_prob(i)<<endl;
    		}
		    if(check)
		    {
			cout<<"Error: hmm class:: BaumWelchTraining : "
			    <<"error in derive_transition_probs_from_transition_parameters()"<<endl;
		    }
		}

		if(trainTTP)
		{
		    // put the TTP back to Original TP for this iteratio
		    for(i=0; i<TTPsize; i++)
		    {
			this->model[TP->get_TTP_from(i)].set_transition_prob(TP->get_TTP_to(i),TP->get_TTP_score(i));
		    }
		}
	 
		
		if(trainEP)
		{
		    // put the trained emission prob back to Original EP for this iteration	     
		    check += derive_emission_probs(this,MP,EP);
		    if(check)
		    {
			cout<<"Error: hmm class:: BaumWelchTraining : "
			    <<"error in derive_emission_probs()"<<endl;
		    }
		}
		// print transition probs
		
		if((trainTTP)||(trainFTP))
		{
		    for(i=0; i<number_of_states-1; i++)
		    {
			cout<<"From ["<<i<<"] to "<<endl;
			model[i].print_transition_probs(cout);
		    }		
		}

		//print emission
		if(trainEP)
		{
		    for(i=1; i<number_of_states-1; i++)
		    {
			cout<<"State : "<<i<<endl;
			model[i].print_emission_probs(cout);
		    }
		}

		//modify the score array after an interation
		check += this->calculate_scores_from_probs();
		if(check)
		{
		    cout<<"Error: hmm class:: BaumWelchTraining : "
			<<"error in training FEP: calculate_scores_from_probs() "
			<<"in loop("<<iCount<<endl;
		}

		check += this->check_consistency_of_probs();

		if(check)
		{
		    cout<<"Error: hmm class:: BaumWelchTraining : "
			<<"error in training TTP: check_consistency_of_probs() "
			<<"in loop("<<iCount<<endl;
		}

		
		cout<<"prev_forward_score : "<<prev_forward_score<<endl;
		cout<<"cur_forward_score : "<<cur_forward_score<<endl;
		cout<<"Number of sequences : "<<SeqCount<<endl;
		log_odd_ratio = abs((cur_forward_score - prev_forward_score)/SeqCount);

		cout<<"log_odd_ratio : "<<log_odd_ratio<<endl<<endl;
		
		iCount++;
	    } // finished training	    

		cout<<"Total time used for the Baum-Welch training : "<<t<<" s."<<endl;
	}	
    }
    return check;
}

int Hmm::Hirschberg_tube(const Hmm &mirror,
			 const Sequence &x, const Sequence &y,
			 const unsigned long long int max_volume, 
			 Tube<int> &tube) {

    int check = 0;

    if (this->get_mirrored() != 0) {
	
	cout << "ERROR: Hmm::Hirschberg_tube: this Hmm is not unmirrored.\n";
	check++;
    }
    if (mirror.get_mirrored() != 1) {
	
	cout << "ERROR: Hmm::Hirschberg_tube: mirror Hmm is unmirrored.\n";
	check++;
    }
    if (tube.Empty()) {
	
	int max_y = y.length()-1;
	if(max_y<0)
	{
	    max_y = 0;
	}
	Tube<int> default_tube(x.length(), 2);
	
	cout << "Hirschberg_tube: setting default tube.\n" << flush;
	
	for (int i = 0; i < x.length(); i++) 
	{
	    
	    default_tube.SetElement(i, 0, 0);
	    default_tube.SetElement(i, 1, max_y);

	
	    tube = default_tube;
	}
    }
    if (tube.Lx() != x.length()) {
	
	cout << "ERROR: Hmm::Hirschberg_tube: length of tube (" << tube.Lx() << ") != length of sequence X ("
	     << x.length() << ").\n";
	check++;
    }
    if (check == 0) {
	
	int max_deltax, i;
	vector<int> Start_point(3), End_point(3);
      
	Start_point[0] = 0; 
	Start_point[1] = 0;
	Start_point[2] = 0;
	End_point[0]   = x.length();
	End_point[1]   = y.length();
	End_point[2]   = this->get_number_of_states()-1;
	
	max_deltax = 0;
	for (i = 1; i < (this->get_number_of_states()-1); i++) {
	    if (this->model[i].get_letters_to_read_x() > max_deltax) {max_deltax = this->model[i].get_letters_to_read_x();}
	}
      
	// allocate memory for global state path
      
	this->allocate_memory_for_state_path(x.length()+y.length()+2, -1, 0.0, -1, Logzero);

	// run Hirschberg
	
	check += this->Hirschberg_tube_internal(mirror, x, y, Start_point, End_point, max_deltax, max_volume, tube);
	
	// if Hirschberg ok, condense solution discarding long one, otherwise delete memory for state path
	
	if (check == 0) {
	    
	    
	    check += check_consistency_of_solution(&x, &y);
	    
	    cout << "Hmm::Hirschberg_tube: after function check_consistency_of_solution check = " << check << "\n" << flush;
	    
	    if (check == 0) {
		condense_solution(0);
	    }
	    else {
		cout << "ERROR: Hmm::Hirschberg_tube: check_consistency_of_solution failed.\n";      
	    }
	}
	else {
	    cout << "ERROR: Hmm::Hirschberg_tube: Hirschberg_tube_internal did not terminate successfully.\n";      
	    this->reset_variables_for_state_path();
	}
    }
    return(check);
}

int Hmm::Hirschberg_tube_internal(const Hmm &mirror,
				      const Sequence &x, const Sequence &y,
				      const vector<int> Start_point, const vector<int> End_point,
				      const int max_deltax,
				      const unsigned long long int max_volume, 
				      Tube<int> &tube) {

    cout << "Hmm::Hirschberg_tube_internal: Start_point[" << Start_point[0] << "]";
    if(pair)
    {
	cout<<"[" << Start_point[1] << "]";
    } 
    cout <<"["<< Start_point[2] << "] End_point[" << End_point[0] << "]";
    if(pair)
    {
	cout<<"[" << End_point[1] << "]";
    }
    cout<<"[" << End_point[2] << "] max_deltax = " << max_deltax 
	<< "\t max_volume = " << max_volume << "\n" << flush;

    int i,j,k, check;
    unsigned long long int tube_volume;
  
    check = 0;

    tube_volume = 0;
    for (i=Start_point[0]-this->model[Start_point[2]].get_letters_to_read_x(); i<End_point[0]; i++) {
	tube_volume += (tube.GetElement(i, 1) - tube.GetElement(i, 0) + 1) * this->get_number_of_states();
    }
    
    cout<<"tube_volume : "<<tube_volume<<endl;
    
    if (tube_volume <= max_volume) {

	StatePath state_path;

	// calculate viterbi_tube with fixed start and end points

	check += this->viterbi_tube(x, y, Start_point, End_point, tube, state_path);

	// add local state path to global solution
	
	if (check == 0) {

	    this->add_local_state_path(&x, &y, &state_path);
      
	    if (check != 0) {
		cout << "ERROR: Hirschberg_tube_internal: add_local_solution did not complete successfully.\n" << flush;
	    }
	}
	else {
	    cout << "ERROR: Hirschberg_tube_internal: viterbi_tube did not complete successfully.\n" << flush;
	}
    }
    else {
	
	const int tube_width        = max_deltax+1;
	const int tube_height       = y.length()+1;
	const int states            = this->get_number_of_states();
	const unsigned long long int strip_volume = tube_width * tube_height * states;
	
	if ((2 * strip_volume) > max_volume) {
	    
	    cout << "ERROR: Hirschberg_tube_internal: 2 * strip_volume = " << 2 * strip_volume << " > max_volume = "
		 << max_volume << "\n" << flush;
	    check++;
	}
	else {

	    int x_middle_volume;
	    unsigned long long int test_volume;
	    const int plus = static_cast<int>(ceil(static_cast<float>(max_deltax)/2.));    
	    
	    // get x point of middle volume = x_middle_volume
	    
	    x_middle_volume = 0;
	    test_volume     = 0;
	    for (i=Start_point[0]-this->model[Start_point[2]].get_letters_to_read_x(); i<End_point[0]; i++) {
		test_volume += (tube.GetElement(i, 1) - tube.GetElement(i, 0) + 1) * this->get_number_of_states();
		if (test_volume*2 < tube_volume) {x_middle_volume = i;}
		else {break;}
	    }
	    x_middle_volume += plus;
	    
	    // check that x_middle_volume falls within lower and upper bound
	    
	    if (x_middle_volume > (End_point[0]-1+this->model[End_point[2]].get_letters_to_read_x())) {
		x_middle_volume = End_point[0]-1+this->model[End_point[2]].get_letters_to_read_x();

	    }
	    if ((x_middle_volume-max_deltax) < (Start_point[0]-this->model[Start_point[2]].get_letters_to_read_x())) {
		x_middle_volume = Start_point[0]-this->model[Start_point[2]].get_letters_to_read_x()+max_deltax;

	    }

	    // get forward and backward strip
      
	    if ((x_middle_volume <= (End_point[0]-1+this->model[End_point[2]].get_letters_to_read_x())) &&
		((x_middle_volume-max_deltax) >= (Start_point[0]-this->model[Start_point[2]].get_letters_to_read_x()))) {
	
		vector<int> end_point_b(3), left_point(3), right_point(3);
	
		// allocate memory for strips
	
		ViterbiElement ****f_strip = new ViterbiElement***[tube_width];
		ViterbiElement ****b_strip = new ViterbiElement***[tube_width];
		
		for (i = 0; i < tube_width; i++) {    
		    f_strip[i] = new ViterbiElement**[tube_height];
		    b_strip[i] = new ViterbiElement**[tube_height];
	  
		    for (j = 0; j < tube_height; j++) {
			f_strip[i][j] = new ViterbiElement*[states];
			b_strip[i][j] = new ViterbiElement*[states];
	    
			for (k = 0; k < states; k++) {
			    f_strip[i][j][k] = NULL;
			    b_strip[i][j][k] = NULL;
			}
		    }
		}
	
		// calculate forward and backward strip
	
		end_point_b[0] = End_point[0]-1-this->model[End_point[2]].get_letters_to_read_x();
		end_point_b[1] = End_point[1]-1-this->model[End_point[2]].get_letters_to_read_y();
		end_point_b[2] = End_point[2];
	
		if (check == 0) {
		    check += this->viterbi_tube_strip(x, y, Start_point, x_middle_volume+1, tube, f_strip);
		    if (check != 0) {
			cout << "ERROR: Hmm::Hirschberg_tube_internal:viterbi_tube_strip Start_point[" << Start_point[0] << "]["
		 << Start_point[1] << "][" << Start_point[2] << "] end_x = " << x_middle_volume+1 << "\n" << flush;
		    }
		}
		if (check == 0) {
		    check += mirror.viterbi_tube_strip_backwards(x, y, x_middle_volume-max_deltax, end_point_b, this, 
						       tube, b_strip);
		    if (check != 0) {
			cout << "ERROR: Hmm::Hirschberg_tube_internal:viterbi_tube_strip_backwards Start_point[" 
			     << Start_point[0] << "][" << Start_point[1] << "][" << Start_point[2] << "] end_x = " 
			     << x_middle_volume+1 << "\n" << flush;
		    }
		}

		// merge strips
	
		if (check == 0) {
		    check += this->combine_strips(f_strip, b_strip, x, y, tube_width, tube_height, x_middle_volume, 
						  left_point, right_point);
		    if (check != 0) {
			cout << "ERROR: Hmm::Hirschberg_tube_internal:combine_strips Start_point[" 
			     << Start_point[0] << "][" << Start_point[1] << "][" << Start_point[2] 
			     << "] End_point[" << End_point[0] << "][" << End_point[1] << "][" << End_point[2]  
			     << "]\n" << flush;
		    }
		}

		// free memory
	
		for (i = 0; i < tube_width; i++) {
		    for (j = 0; j < tube_height; j++) {
			for (k = 0; k < states; k++) {
			    
			    if (f_strip[i][j][k]) delete f_strip[i][j][k];
			    if (b_strip[i][j][k]) delete b_strip[i][j][k];
			    f_strip[i][j][k] = NULL;
			    b_strip[i][j][k] = NULL;
			}
			if (f_strip[i][j]) delete [] f_strip[i][j];
			if (b_strip[i][j]) delete [] b_strip[i][j];
			f_strip[i][j] = NULL;
			b_strip[i][j] = NULL;
		    }
		    if (f_strip[i]) delete [] f_strip[i];
		    if (b_strip[i]) delete [] b_strip[i];
		    f_strip[i] = NULL;
		    b_strip[i] = NULL;
		}
		if (f_strip) delete [] f_strip;
		if (b_strip) delete [] b_strip;
		f_strip = NULL;
		b_strip = NULL;
		
		// start new iteration

		// left rectangle
		
		if (check == 0) {
		    
		    cout << "Hmm::Hirschberg_tube_internal: calculating left rectangle\n" << flush;
		    check += this->Hirschberg_tube_internal(mirror, x, y, Start_point, left_point, max_deltax, 
							    max_volume, tube);
		    
		    if (check != 0) {
			cout << "ERROR: Hmm::Hirschberg_tube_internal:Hirschberg_tube_internal of left rectangle Start_point["
			     << Start_point[0] << "][" << Start_point[1] << "][" << Start_point[2] << "] End_point["
			     << left_point[0] << "][" << left_point[1] << "][" << left_point[2] << "]\n" << flush;
		    }
		}

		// right rectangle

		if (check == 0) {

		    cout << "Hmm::Hirschberg_tube_internal: calculating right rectangle\n" << flush;
		    check += this->Hirschberg_tube_internal(mirror, x, y, right_point, End_point, max_deltax, 
							    max_volume, tube);
		    
		    if (check != 0) {
			cout << "ERROR: Hmm::Hirschberg_tube_internal:Hirschberg_tube_internal of right rectangle Start_point["
			     << right_point[0] << "][" << right_point[1] << "][" << right_point[2] << "] End_point["
			     << End_point[0] << "][" << End_point[1] << "][" << End_point[2] << "]\n" << flush;
		    }
		}
		
	    }
	    else {
		cout << "ERROR: Hmm::Hirschberg_tube_internal: x_middle_volume (" << x_middle_volume << ") not in ["
		     << Start_point[0]-this->model[Start_point[2]].get_letters_to_read_x() << ", " 
		     << End_point[0]-1+this->model[End_point[2]].get_letters_to_read_x() << "]\n" << flush;
		check++;
	    }
	} // else: 2 * strip_volume < max_volume
    } // else: tube_volume > max_volume
    
    return(check);
}

#endif // #ifndef _INTEL

