 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare scoring functions with which sequences can be assinged scores

   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/scoring_functions.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#include <fstream.h>
#include "scoring_functions.h"


int set_scores_for_annotated_sequence(// input and output
				      Sequence* const sequence,
				      // input
				      const int            x_or_y, // 0 = x, 1 = y
				      Hmm* const pairhmm)
{
    int check = 0;
    
    if (sequence == NULL) 
    {
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : sequence is NULL.\n";
	check++;
    }
    else {
	if(sequence->get_constraint()==1)
	{
	    if (sequence->annotation_labels_exist() != 1) 
	    {
		
		cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : Sequence sequence has no labels "
		     << "associated with it.\n";
		check++;
	    }	    
	}else{
	    return check;
	}
    }
    if (pairhmm == NULL) {
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : Hmm pairhmm is NULL.\n";
	check++;
    }
    if ((x_or_y != 0) && (x_or_y != 1)) {
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : x_or_y (" << x_or_y 
	     << ") has to be 0 or 1.\n";
	check++;
    }

    if (check == 0) 
    {
	
	int i, j = 0;

	const int n_of_states = pairhmm->get_number_of_states();
	
	int* numbers_of_child_states_of_emitxy_states     = new int[n_of_states];
	int* numbers_of_states_with_virtual_child_state   = new int[n_of_states];
	int* numbers_of_states_with_own_special_emissions = new int[n_of_states];
	int* numbers_of_child_states_x_or_y               = new int[n_of_states];
	
	int sequence_delta_x = 0;
	int sequence_delta_y = 0;
	
	if(x_or_y==0){
	    sequence_delta_x++;
	}else if(x_or_y==1){
	    sequence_delta_y++;
	}


	int count_chi = 0; // counts entries in array numbers_of_child_states_of_emitxy_states
	int count_vir = 0; // counts entries in array numbers_of_states_with_virtual_child_state
	int count_own = 0; // counts entries in array numbers_of_states_with_own_special_emissions
	int count_all = 0; // counts entries in array numbers_of_child_states_x_or_y
	
	int   child        = 0;  // number of child_x or child_y state
	int   child_x      = 0;
	int   child_y      = 0;
	int   delta_x      = 0;
	int   delta_y      = 0;
	int   state_number = 0;
	
	for (i=0; i<n_of_states; i++) {
	    
	    state_number = (*pairhmm)[i]->get_number_of_state();
	    delta_x = (*pairhmm)[i]->get_letters_to_read_x();
	    delta_y = (*pairhmm)[i]->get_letters_to_read_y();
	    	    
	    if ((*pairhmm)[i]->get_special_emission() == 1) { // if state has special emissions	

		child_x = (*pairhmm)[i]->get_number_of_child_state_x();
		child_y = (*pairhmm)[i]->get_number_of_child_state_y();
		
		if (compare_state_type(sequence_delta_x,sequence_delta_y,1,0)==1) // sequence_type = EmitX
		{
		    child = child_x;
		}
		else if (compare_state_type(sequence_delta_x,sequence_delta_y,0,1)==1) // sequence_type = EmitY
		{
		    child = child_y;
		}
		
		// if this is a EmitXY and the child state is virtual
		
		if ((delta_x>0)&&(delta_y>0)) { 
		    
		    numbers_of_child_states_of_emitxy_states[count_chi] = child;
		    count_chi++;

		    if (child > n_of_states) {
			
			numbers_of_states_with_virtual_child_state[count_vir] = state_number;
			count_vir++;
			
		    }	
		}
		// if state type is EmitX or EmitY (depending on the sequence)
		else if (compare_state_type(sequence_delta_x,sequence_delta_y,delta_x,delta_y)==1) 
		{  
		    
		    numbers_of_child_states_x_or_y[count_all] = child;
		    count_all++;

		    if ((state_number == child_x) &&  // if state has got its own special emissions
			(state_number == child_y)) {
			
			numbers_of_states_with_own_special_emissions[count_own] = state_number;
			count_own++;
			
		    } // if state has got its own special emissions
		} // else if state_type == sequence_type
	    } // if state has special emissions
	} // loop over all states

	// make checks:
	//
	// - 1.) every entry in array numbers_of_child_states_x_or_y has to be found in array 
	//       numbers_of_states_with_own_special_emissions
	// - 2.) every entry in numbers_of_states_with_own_special_emissions must be unique
	// - 3.) every entry in array numbers_of_child_states_of_emitxy_states has to be found
	//       in array numbers_of_states_with_own_special_emissions
	
	int entry       = 0;
	int found_entry = 0;

	// 1.)

	if (check == 0) {
          
	    for (i=0; i<count_all; i++) {
		if (check == 0) {
		    
		    entry       = numbers_of_child_states_x_or_y[i];
		    found_entry = 0;

		    for (j=0; j<count_own; j++) {
			if (numbers_of_states_with_own_special_emissions[j] == entry) {

			    found_entry++;
			    break;
			}
		    }
	  
		    if (found_entry == 0) {
			cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : (check 1) "
			     << "child state with number " << entry 
			     << " has no associated state with special emission probs.\n";
			check++;
			break;
		    }
		} // if check == 0
	    }
	} // if check == 0
	
	// 2.)
	
	if (check == 0) {
	   
	    int new_count_own = count_own;

	    check += make_elements_of_array_unique(// input and output
		&new_count_own,
		numbers_of_states_with_own_special_emissions);
	    
	    if (new_count_own != count_own) {
		cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : "
		     << "the numbers of states with their own special emissions is not unique.\n";
		check++;	
	    }
	    if (check != 0) {
		cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : error occurred in function "
		     << "make_elements_of_array_unique.\n";
	    }
	} // if check == 0
	
	// 3.)
	
	if (check == 0) {

	    entry       = 0;
	    found_entry = 0;
	    
	    for (i=0; i<count_chi; i++) {
		if (check == 0) {
		    
		    entry       = numbers_of_child_states_of_emitxy_states[i];
		    found_entry = 0;

		    if (entry > n_of_states) { // if this is a virtual child state
			found_entry++;
		    }
		    else { // if this is not a virtual child state
			
			for (j=0; j<count_own; j++) 
			{
			    if (numbers_of_states_with_own_special_emissions[j] == entry) 
			    {				
				found_entry++;
				break;
			    }
			}
		    }
		    
		    if (found_entry == 0) {
			cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : (check 3) "
			     << "child state with number " << entry 
			     << " has no associated state with special emission probs.\n";
			check++;
			break;
		    }
		} // if check == 0
	    }
	} // if check == 0
	
	// loop over states which have their own special emissions and score the sequence
	// ------------------------------------------------------------------------------
	
	if (check == 0) {
	    
	    for (i=0; i<count_own; i++) {
		if (check == 0) {	
		    
		    state_number = numbers_of_states_with_own_special_emissions[i];

		    check += set_scores_for_annotated_sequence_for_given_state(// input
			sequence,
			x_or_y,		
			(*pairhmm)[state_number]);

		    if (check != 0) {
			cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : error occurred in function "
			     << "set_scores_for_annotated_sequence for state " << state_number << ".\n";
		    }
		} // if check == 0
	    } // loop over states which have their own special emissions
	    
	    // loop over states with virtual child states
	    // ------------------------------------------------------------------------------
	    
	    int  already_implemented = 0;
	    int  count_implemented   = 0;
	    int* numbers_of_already_implemented_virtual_states = new int[n_of_states];      
	    
	    for (i=0; i<count_vir; i++) { // loop over entries of array numbers_of_states_with_virtual_child_state
		if (check == 0) {	
		    
		    state_number = numbers_of_states_with_virtual_child_state[i];
		    
		    child_x = (*pairhmm)[state_number]->get_number_of_child_state_x();
		    child_y = (*pairhmm)[state_number]->get_number_of_child_state_y();
		    if      (compare_state_type(sequence_delta_x,sequence_delta_y,1,0) == 1)  // sequence_type == EmitX
		    {
			child = child_x;
		    }
		    else if (compare_state_type(sequence_delta_x,sequence_delta_y,0,1) == 1)  // sequence_type == EmitY
		    {
			child = child_y;
		    }
		    
		    // check if virtual state was already implemented
		    
		    for (j=0; j<count_implemented; j++) {
			if (numbers_of_already_implemented_virtual_states[j] == child) {
			    already_implemented++;
			    break;
			}
		    }
		    
		    if (already_implemented == 0) { // if virtual state was not already implemented
			
			numbers_of_already_implemented_virtual_states[count_implemented] = child;
			count_implemented++;
			
			check += set_scores_for_annotated_sequence_for_given_state(// input
			    sequence,
			    x_or_y,
			    (*pairhmm)[state_number]);
		
			if (check != 0) {
			    cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence : "
				 << "error occurred in function set_scores_for_annotated_sequence for state " 
				 << state_number << ".\n";
			}
		    } // if virtual state was not already implemented
		} // if check == 0
	    } // for loop
	    
	    if (numbers_of_already_implemented_virtual_states) delete [] numbers_of_already_implemented_virtual_states;
	    numbers_of_already_implemented_virtual_states = NULL;
	    
	} // if check == 0
	
	if (numbers_of_child_states_of_emitxy_states) delete [] numbers_of_child_states_of_emitxy_states;
	numbers_of_child_states_of_emitxy_states = NULL;
	
	if (numbers_of_states_with_own_special_emissions) delete [] numbers_of_states_with_own_special_emissions;
	numbers_of_states_with_own_special_emissions = NULL;

	if (numbers_of_child_states_x_or_y) delete [] numbers_of_child_states_x_or_y;
	numbers_of_child_states_x_or_y = NULL;
	
	if (numbers_of_states_with_virtual_child_state) delete [] numbers_of_states_with_virtual_child_state;
	numbers_of_states_with_virtual_child_state = NULL;
	
    } // if check == 0
    
    return(check);
}

int set_scores_for_annotated_sequence_for_given_state(// input
						      Sequence*      const sequence,
						      const int            x_or_y, // 0 = x, 1 = y
						      const Hmm_State* const state)
{
    int check = 0;

    bool***     labels = NULL;
    double***    label_scores = NULL;
   
    if (sequence == NULL) {
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : "
	     << "Sequence sequence is NULL.\n";
	check++;
    }
    else {
	if (sequence->annotation_labels_exist() != 1) {
	    cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : "
		 << "Sequence sequence has no labels associated with it.\n";
	    check++;
	}
	if (sequence->length() < 1) {
	    cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : "
		 << "Sequence sequence length (" << sequence->length() << ") too short (has to be > 0).\n";
	    check++;
	}    
	if (check == 0) 
	{
	    check += sequence->get_combine_annotation_labels(&labels,
							     &label_scores);

	    if (check != 0) 
	    {
		cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : error "
		     << "occurred in function Sequence::get_annotation_labels.\n";
	    }
	}	  
    }
    if ((x_or_y != 0) && (x_or_y != 1)) {
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : x_or_y (" << x_or_y 
	     << ") has to be 0 or 1.\n";
	check++;
    }

    if (state == NULL) {
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : "
	     << "Hmm_State state is NULL.\n";
	check++;
    }
    if ((compare_state_type(state->get_letters_to_read_x(),state->get_letters_to_read_y(),1,0)!=1)&&
	(compare_state_type(state->get_letters_to_read_x(),state->get_letters_to_read_y(),0,1)!=1)&&
	(compare_state_type(state->get_letters_to_read_x(),state->get_letters_to_read_y(),1,1)!=1))
    {	
	cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : Hmm_State state "
	     << "is not of type EmitX or EmitY or EmitXY.\n"
	     << "number_of_letters_to_read_x("<<state->get_letters_to_read_x()
	     <<") number_of_letters_to_read_y("<<state->get_letters_to_read_y()
	     <<").\n"<<flush;
	check++;
    }
    int i = 0; 
    int j = 0;
    int k = 0;

    if (check == 0) {

	const int delta_x = state->get_letters_to_read_x();
	const int delta_y = state->get_letters_to_read_y();
	int delta  = 0;
	int offset = 0;
	int child  = 0;
	
	const int length = sequence->length();

	int number_of_sources = sequence->get_number_of_sources();
	int number_of_annotation_labels = sequence->get_number_of_annotation_labels();
	int* number_of_each_annotation_label = sequence->get_number_of_each_annotation_label();

		
	int sequence_delta_x = 0;
	int sequence_delta_y = 0;
	
	if(x_or_y==0){
	    sequence_delta_x++;
	}else if(x_or_y==1){
	    sequence_delta_y++;
	}
	
	if (compare_state_type(delta_x,delta_y,1,1)==1){ // if state is EmitXY

	    if (compare_state_type(sequence_delta_x,sequence_delta_y,1,0)==1)
	    {		
		delta  = state->get_letters_to_read_x();
		offset = 0; 
		child  = state->get_number_of_child_state_x();
	    }
	    else if ((compare_state_type(sequence_delta_x,sequence_delta_y,0,1)==1)) 
	    { 	
		delta  = state->get_letters_to_read_y();
		offset = state->get_letters_to_read_x();
		child  = state->get_number_of_child_state_y();
	    }
	}
	else if (compare_state_type(delta_x,delta_y,1,0)==1)  
	{	    
	    delta  = state->get_letters_to_read_x();
	    offset = 0;
	    child  = state->get_number_of_child_state_x();
	}
	else if (compare_state_type(delta_x,delta_y,0,1)==1)  
	{	    
	    delta  = state->get_letters_to_read_y();
	    offset = 0;
	    child  = state->get_number_of_child_state_y();
	}

	Score* scores_for_sequence = new Score[length];
	
	int j_in_state = 0;
	
	int position_not_ok = 0;

	bool* score_labels = new bool[sequence->get_number_of_annotation_labels()];
	for(i=0; i<sequence->get_number_of_annotation_labels(); i++)
	{	    
	    score_labels[i] = sequence->get_score_labels(i);
	}
			
	// loop over every position in sequence

	for (i=0; i<length; i++) 
	{	    
	    if(check)
	    {
		break;
	    }
	    position_not_ok = 0;
	    
	    scores_for_sequence[i] = 0;
	    
	    // loop over every position in state

	    for (j=0; j<delta; j++) 
	    {
		
		if(check)
		{
		    break;
		}
		
		j_in_state = offset+j;

		if(position_not_ok)
		{
		    break;
		}
		
		if ((i+j) < length) 
		{

		    for(k=0;k<number_of_annotation_labels;k++)
		    {

			if(check)
			{
			    break;
			}
					
			if(!score_labels[k])
			{
			    if((!check_undefine(labels[i+j][k],
					      number_of_each_annotation_label[k]))&&
			       
			       (!labels[i+j][k][state->get_state_labels(k,j_in_state)]))
			    {
				position_not_ok++;
				break;
			    }
			}
		    }

		    if(!position_not_ok)
		    {
			if(check)
			{
			    break;
			}

			for(k=0; k<number_of_annotation_labels; k++)
			{
			    if(score_labels[k])
			    {
				if( check_undefine(labels[i+j][k],
						   number_of_each_annotation_label[k]))
				{
				    scores_for_sequence[i] += 0;
				}
				else if(!labels[i+j][k][state->get_state_labels(k,j_in_state)])
				{
				    position_not_ok++;
				    break;
				}else{		
				 
				    if((label_scores[i+j][k][state->get_state_labels(k,j_in_state)]<0)||
				       (label_scores[i+j][k][state->get_state_labels(k,j_in_state)]>1))
				    {
					cout << "ERROR : scoring_functions : label_scores["<<i+j<<"]["
					     <<k<<"]["<<state->get_state_labels(k,j_in_state)<<"]("
					     <<label_scores[i+j][k][state->get_state_labels(k,j_in_state)]<<") out of range [0,1]"<<endl;
					check++;
					break;
				    }
				    
				    scores_for_sequence[i]+= static_cast<Score>(log(label_scores[i+j][k][state->get_state_labels(k,j_in_state)]));

				}
			    }
			}
		    }			  	      		    
		    
		} // if (i+j) < sequence length
	    }
	    
	    
	    if(position_not_ok)
	    {
		scores_for_sequence[i] = Logzero;
	    }
	} // loop over positions in sequence
	
	// set scores for state

	if (compare_state_type(delta_x,delta_y,1,1)!=1)
	{	   
	    check += sequence->set_scores(state->get_number_of_state(),
					  scores_for_sequence);
	}else{
	    check += sequence->set_scores(child,
					  scores_for_sequence); 
	}

	if (check != 0) {
	    cout << "ERROR : scoring_functions : set_scores_for_annotated_sequence_for_given_state : "
		 << "error occurred in function Sequence::set_scores for state number " 
		 << state->get_number_of_state() << ".\n";
	}      

	if (scores_for_sequence) delete [] scores_for_sequence;
	scores_for_sequence = NULL;
	
    } // if check == 0

    const int length = sequence->length();
    int number_of_annotation_labels = sequence->get_number_of_annotation_labels();
    
    if(labels){
	for(i=0; i<length; i++){
	    if(labels[i]){
		for(j=0; j<number_of_annotation_labels; j++){
		    if(labels[i][j]) delete [] labels[i][j];
		    labels[i][j] = NULL;
		}
		delete [] labels[i];
	    }
	    labels[i] = NULL;
	}
	delete [] labels;
    }
    labels = NULL;
    
    if(label_scores){
	for(i=0; i<length; i++){
	    if(label_scores[i]){
		for(j=0; j<number_of_annotation_labels; j++){
		    if(label_scores[i][j]) delete [] label_scores[i][j];
		    label_scores[i][j] = NULL;
		}
		delete [] label_scores[i];
	    }
	    label_scores[i] = NULL;
	}
	delete [] label_scores;
    }
    label_scores = NULL;

    return(check);
}
