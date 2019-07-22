/* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare sequence-class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/sequence.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
*/

#include <malloc.h>
#include <fstream>
#include <unistd.h>
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
#include <iostream>
#include "sequence.h"
#include "evaluation.h"

/* define constructors */

Sequence::Sequence()
{
    sequence_type=0;
    length_of_sequence=0;
    start_position = 0;
    end_position = 0;
    
    sequence=NULL;
    ac=NULL;

    // annotation
    number_of_sources = 0;
    number_of_annotation_labels = 0;
    number_of_each_annotation_label = NULL;
    annotation_labels = NULL;
    annotation_label_scores = NULL;

    score_labels = NULL;

    constraint = 0;

    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;
    
    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    implementation_status_of_special_transitions.SetNumberofDimensions(0);
    indices_of_special_transitions.SetNumberofDimensions(0);
    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
    priors_of_special_transitions.SetNumberofDimensions(0);
    scores_for_special_transitions.SetNumberofDimensions(0);
    
    implementation_status_of_special_emissions.SetNumberofDimensions(0);
    indices_of_special_emissions.SetNumberofDimensions(0);
    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
    scores_for_special_emissions.SetNumberofDimensions(0);
    
    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 
    
    sp_steps = 0;
    
}

Sequence::Sequence(const char* const char_seq,
		   model_parameters* const MP)
{
    int check = 0;

    sequence_type=0;
    length_of_sequence=0;
    start_position = 0;
    end_position = 0;
    
    sequence=NULL;
    ac=NULL;

    // annotation
    number_of_sources  = 0;
    number_of_annotation_labels = 0;
    number_of_each_annotation_label = NULL;
    annotation_labels = NULL;
    annotation_label_scores = NULL;

    score_labels = NULL;
    constraint = 0;

    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;

    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    implementation_status_of_special_transitions.SetNumberofDimensions(0);
    indices_of_special_transitions.SetNumberofDimensions(0);
    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
    priors_of_special_transitions.SetNumberofDimensions(0);
    scores_for_special_transitions.SetNumberofDimensions(0);
    
    implementation_status_of_special_emissions.SetNumberofDimensions(0);
    indices_of_special_emissions.SetNumberofDimensions(0);
    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
    scores_for_special_emissions.SetNumberofDimensions(0);
    
    if (!char_seq) {
	cout << "ERROR class Sequence::constructor : char_seq NULL\n" << flush;
	check++;
    }
    if (strlen(char_seq)<1) {
	cout << "ERROR class Sequence::constructor : length of char_seq (" << strlen(char_seq) << ") < 1\n" << flush;
	check++;
    }
    if (check == 0) {
	
	int i = 0;
	length_of_sequence=strlen(char_seq);
	sequence=new int[length_of_sequence];
	for (i=0; i<length_of_sequence; i++) {
	    sequence[i]=convert_alphabet_to_int(MP,char_seq[i]);
	}      
	start_position = 1;
	end_position   = length_of_sequence;
    }
    
    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 
    
    sp_steps = 0;
    			    
}

Sequence::Sequence(const int seq_type, 
		   const char* const char_seq,
		   model_parameters* const MP)
{
    // note: - probs_or_scores == 1 : use probs
    //         probs_or_scores == 0 : use scores
    
    int check = 0;

    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 
    
    sp_steps = 0;
    
    sequence_type=0;
    length_of_sequence=0;
    start_position = 0;
    end_position = 0;
    
    sequence=NULL;
    ac=NULL;

    number_of_sources = 0;
    number_of_annotation_labels = 0;
    number_of_each_annotation_label = NULL;
    annotation_labels = NULL;
    annotation_label_scores = NULL;

    score_labels = NULL;
    
    constraint = 0;

    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;
    
    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    implementation_status_of_special_transitions.SetNumberofDimensions(0);
    indices_of_special_transitions.SetNumberofDimensions(0);
    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
    priors_of_special_transitions.SetNumberofDimensions(0);
    scores_for_special_transitions.SetNumberofDimensions(0);
    
    implementation_status_of_special_emissions.SetNumberofDimensions(0);
    indices_of_special_emissions.SetNumberofDimensions(0);
    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
    scores_for_special_emissions.SetNumberofDimensions(0);
    
    if (!char_seq) {
	cout << "ERROR class Sequence::constructor : char_seq NULL\n" << flush;
	check++;
    }
    if (strlen(char_seq)<1) {
	cout << "ERROR class Sequence::constructor : length of char_seq (" << strlen(char_seq) << ") < 1\n" << flush;
	check++;
    }
   
    if (check == 0) {
	
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	
	sequence_type=seq_type;
	length_of_sequence=strlen(char_seq);

	sequence=new int[length_of_sequence];
	for (int i=0; i<length_of_sequence; i++) {
	    sequence[i]=convert_alphabet_to_int(MP,char_seq[i]);
	}   
	score_labels = new bool[MP->get_Total_Number_of_Annotation_Labels()];
	for(i=0; i<MP->get_Total_Number_of_Annotation_Labels(); i++)
	{
	    score_labels[i] = MP->get_score_of_Annotation_Label(i);	    	    
	}
	int n_of_special_transitions = 0;
	int n_of_special_emissions = MP->get_Number_of_Special_Emissions();
	int max_n_of_child_state = MP->get_Max_n_of_Child_State() + 1;
	
	constraint = 0;

	start_position = 1;
	end_position   = length_of_sequence;
	
	if ((n_of_special_transitions > 0) && (max_n_of_child_state > 0)) {
	    
	    number_of_special_transitions = n_of_special_transitions;
	    
	    implementation_status_of_special_transitions.SetNumberofDimensions(2);
	    implementation_status_of_special_transitions.SetDimension(0, max_n_of_child_state, 0);
	    implementation_status_of_special_transitions.SetDimension(1, max_n_of_child_state, 0);
	    
	    indices_of_special_transitions.SetNumberofDimensions(2);
	    indices_of_special_transitions.SetDimension(0, max_n_of_child_state, 0);
	    indices_of_special_transitions.SetDimension(1, max_n_of_child_state, 0);
	
	    priors_of_special_transitions.SetNumberofDimensions(1);
	    priors_of_special_transitions.SetDimension(0, number_of_special_transitions, 0.0);
	    
	    scores_for_special_transitions.SetNumberofDimensions(2);
	    scores_for_special_transitions.SetDimension(0, number_of_special_transitions, Logzero);
	    scores_for_special_transitions.SetDimension(1, length_of_sequence, Logzero);
	}
	
	if ((n_of_special_emissions > 0) && (max_n_of_child_state > 0)) {
	    
	    number_of_special_emissions = n_of_special_emissions;
	    
	    implementation_status_of_special_emissions.SetNumberofDimensions(1);
	    implementation_status_of_special_emissions.SetDimension(0, max_n_of_child_state, 0);
	    
	    indices_of_special_emissions.SetNumberofDimensions(1);
	    indices_of_special_emissions.SetDimension(0, max_n_of_child_state, 0);
	    
	    scores_for_special_emissions.SetNumberofDimensions(2);
	    scores_for_special_emissions.SetDimension(0, number_of_special_emissions, Logzero);
	    scores_for_special_emissions.SetDimension(1, length_of_sequence, Logzero);
	}
    }
} 

Sequence::Sequence(const int* const seq, const int length)
{
    int check = 0;			    

    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 
    
    sp_steps = 0;
    
    sequence_type=0;
    length_of_sequence=0;
    start_position = 0;
    end_position = 0;
    
    sequence=NULL;
    ac=NULL;

    // annotation
    number_of_sources = 0;
    number_of_annotation_labels = 0;
    number_of_each_annotation_label = NULL;
    annotation_labels = NULL;
    annotation_label_scores = NULL;

    score_labels = NULL;

    constraint = 0;

    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;
    
    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    implementation_status_of_special_transitions.SetNumberofDimensions(0);
    indices_of_special_transitions.SetNumberofDimensions(0);
    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
    priors_of_special_transitions.SetNumberofDimensions(0);
    scores_for_special_transitions.SetNumberofDimensions(0);
    
    implementation_status_of_special_emissions.SetNumberofDimensions(0);
    indices_of_special_emissions.SetNumberofDimensions(0);
    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
    scores_for_special_emissions.SetNumberofDimensions(0);

    if (!seq) {
	cout << "ERROR class Sequence::constructor : sequence NULL\n" << flush;
	check++;
    }
    if (length<1) {
	cout << "ERROR class Sequence::constructor : length of sequence (" << length << ") < 1\n" << flush;
	check++;
    }
    if (!check) {
	
	int i = 0;
	length_of_sequence=length;
	sequence=new int[length];
	for (i=0; i<length; i++) {
	    sequence[i]=seq[i];
	}
	start_position = 1;
	end_position = length_of_sequence;
    }
}

Sequence::Sequence(const int seq_type, 
		   const int* const seq, 
		   const int length, 
		   const bool probs_or_scores,
		   model_parameters* const MP)
{
    // note: - probs_or_scores == 1 : use probs
    //         probs_or_scores == 0 : use scores

    int check = 0;
    
    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 
    
    sp_steps = 0;

    sequence_type = 0;
    length_of_sequence=0;
    start_position = 0;
    end_position = 0;
    
    sequence=NULL;
    ac=NULL;
    
    // annotation
    number_of_sources = 0;
    number_of_annotation_labels = 0;
    number_of_each_annotation_label = NULL;
    annotation_labels = NULL;
    annotation_label_scores = NULL;

    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;
  
    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    implementation_status_of_special_transitions.SetNumberofDimensions(0);
    indices_of_special_transitions.SetNumberofDimensions(0);
    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
    priors_of_special_transitions.SetNumberofDimensions(0);
    scores_for_special_transitions.SetNumberofDimensions(0);

    implementation_status_of_special_emissions.SetNumberofDimensions(0);
    indices_of_special_emissions.SetNumberofDimensions(0);
    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
    scores_for_special_emissions.SetNumberofDimensions(0);
    
    if (!seq) {
	cout << "ERROR class Sequence::constructor : sequence NULL\n" << flush;
	check++;
    }
    if (length<1) {
	cout << "ERROR class Sequence::constructor : length of sequence (" << length << ") < 1\n" << flush;
	check++;
    }
  
    if (check == 0)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	
	sequence_type=seq_type;
	length_of_sequence=length;
	sequence=new int[length];
	for (i=0; i<length; i++) {
	    sequence[i]=seq[i];
	}

	number_of_sources = 0;
	
	number_of_annotation_labels = 0;
	number_of_each_annotation_label = NULL;
	annotation_labels = NULL;
	annotation_label_scores = NULL;

	score_labels = new bool[MP->get_Total_Number_of_Annotation_Labels()];
	for(i=0; i<MP->get_Total_Number_of_Annotation_Labels(); i++)
	{
	    score_labels[i] = MP->get_score_of_Annotation_Label(i);
	}

	constraint = 0;
	
	start_position = 1;
	end_position = length_of_sequence;

	int n_of_special_transitions = 0;
	int n_of_special_emissions = MP->get_Number_of_Special_Emissions();
	int max_n_of_child_state = MP->get_Max_n_of_Child_State() + 1;
	
	if ((n_of_special_transitions > 0) && (max_n_of_child_state > 0)) {

	    number_of_special_transitions = n_of_special_transitions;
	    
	    implementation_status_of_special_transitions.SetNumberofDimensions(2);
	    implementation_status_of_special_transitions.SetDimension(0, max_n_of_child_state, 0);
	    implementation_status_of_special_transitions.SetDimension(1, max_n_of_child_state, 0);
	    
	    indices_of_special_transitions.SetNumberofDimensions(2);
	    indices_of_special_transitions.SetDimension(0, max_n_of_child_state, 0);
	    indices_of_special_transitions.SetDimension(1, max_n_of_child_state, 0);
	    
	    if (probs_or_scores) {
		
		posterior_probs_for_special_transitions.SetNumberofDimensions(2);
		posterior_probs_for_special_transitions.SetDimension(0, number_of_special_transitions, 
								     static_cast<Prob>(0.));
		posterior_probs_for_special_transitions.SetDimension(1, length_of_sequence, static_cast<Prob>(0.));
	    }
	    else {
		priors_of_special_transitions.SetNumberofDimensions(1);
		priors_of_special_transitions.SetDimension(0, number_of_special_transitions, 0.0);
		scores_for_special_transitions.SetNumberofDimensions(2);
		scores_for_special_transitions.SetDimension(0, number_of_special_transitions, Logzero);
		scores_for_special_transitions.SetDimension(1, length_of_sequence, Logzero);
	    }
	}
	if ((n_of_special_emissions > 0) && (max_n_of_child_state > 0)) {
	    
	    number_of_special_emissions = n_of_special_emissions;
	    
	    implementation_status_of_special_emissions.SetNumberofDimensions(1);
	    implementation_status_of_special_emissions.SetDimension(0, max_n_of_child_state, 0);
	    
	    indices_of_special_emissions.SetNumberofDimensions(1);
	    indices_of_special_emissions.SetDimension(0, max_n_of_child_state, 0);
	    
	    if (probs_or_scores) {
		
		posterior_probs_for_special_emissions.SetNumberofDimensions(2);
		posterior_probs_for_special_emissions.SetDimension(0, number_of_special_emissions, 
							     static_cast<Prob>(0.));
		posterior_probs_for_special_emissions.SetDimension(1, length_of_sequence, static_cast<Prob>(0.));
	    }
	    else {
		scores_for_special_emissions.SetNumberofDimensions(2);
		scores_for_special_emissions.SetDimension(0, number_of_special_emissions, Logzero);
		scores_for_special_emissions.SetDimension(1, length_of_sequence, Logzero);
	    }
	}
    }
}

Sequence::Sequence(const Sequence &s)
{
    int check = 0;

    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 

    sp_steps = 0;

    sequence_type=0; 
    length_of_sequence=0;
    start_position = 0;
    end_position = 0;

    sequence=NULL;     

    ac=NULL;

    // annotation
    number_of_sources = 0;
    number_of_annotation_labels = 0;
    number_of_each_annotation_label = NULL;
    annotation_labels = NULL;
    annotation_label_scores = NULL;
    
    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;
    
    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    implementation_status_of_special_transitions.SetNumberofDimensions(0);
    indices_of_special_transitions.SetNumberofDimensions(0);
    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
    priors_of_special_transitions.SetNumberofDimensions(0);
    scores_for_special_transitions.SetNumberofDimensions(0);

    implementation_status_of_special_emissions.SetNumberofDimensions(0);
    indices_of_special_emissions.SetNumberofDimensions(0);
    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
    scores_for_special_emissions.SetNumberofDimensions(0);
    
    if (!s.sequence) {
	cout << "ERROR class Sequence::copy constructor : sequence NULL\n" << flush;
	check++;
    }
    if (s.length_of_sequence<1) {
	cout << "ERROR class Sequence::copy constructor : length_of_sequence (" 
	     << s.length_of_sequence << ") < 1\n" << flush;
	check++;
    }
    
    if (check == 0)
    {
	*this=s;
    }
}

Sequence::~Sequence()
{  
    int i = 0;
    int j = 0;
    int k = 0;

    sequence_type=0;
    start_position = 0;
    end_position = 0;
    
    if (sequence) delete [] sequence;
    sequence=NULL;

    if (ac) delete [] ac;
    ac=NULL;

    length_of_sequence=0;

    // remove annotation_labels
    if(annotation_labels){
	for(i=0; i<length_of_sequence; i++){
	    if(annotation_labels[i]){
		for(j=0; j<number_of_sources; j++){	
		    if(annotation_labels[i][j]){
			for(k=0; k<number_of_annotation_labels; k++){
			    if(annotation_labels[i][j][k]) delete[] annotation_labels[i][j][k];
			    annotation_labels[i][j][k] = NULL;
			}
			delete[] annotation_labels[i][j];
		    }
		    annotation_labels[i][j] = NULL;
		}
		delete [] annotation_labels[i];
	    }
	    annotation_labels[i] = NULL;
	}
	delete[] annotation_labels;
    }
    annotation_labels= NULL;    

    // remove annotation_label_scores
    if(annotation_label_scores){
	for(i=0; i<length_of_sequence; i++){
	    if(annotation_label_scores[i]){
		for(j=0; j<number_of_sources; j++){	
		    if(annotation_label_scores[i][j]){
			for(k=0; k<number_of_annotation_labels; k++){
			    if(annotation_label_scores[i][j][k]) delete[] annotation_label_scores[i][j][k];
			    annotation_label_scores[i][j][k] = NULL;
			}
			delete[] annotation_label_scores[i][j];
		    }
		    annotation_label_scores[i][j] = NULL;
		}
		delete [] annotation_label_scores[i];
	    }
	    annotation_label_scores[i] = NULL;
	}
	delete[] annotation_label_scores;
    }
    annotation_label_scores = NULL;   

    number_of_sources = 0;
    number_of_annotation_labels = 0;
    if(number_of_each_annotation_label) delete[] number_of_each_annotation_label;
    number_of_each_annotation_label = NULL;
    if(score_labels) delete[] score_labels;
    score_labels = NULL;
    
    constraint = 0;

    number_of_implemented_special_transitions=0;
    number_of_special_transitions=0;
    
    number_of_implemented_special_emissions=0;
    number_of_special_emissions=0;
    
    
    if (sp_xsteps) delete [] sp_xsteps;
    if (sp_ysteps) delete [] sp_ysteps;
    if (sp_states) delete [] sp_states;
    if (sp_scores) delete [] sp_scores;
    
    sp_xsteps = NULL; 
    sp_ysteps = NULL; 
    sp_states = NULL; 
    sp_scores = NULL; 
    sp_steps = 0;
    
}

// private functions

int Sequence::set_implementation_status_of_special_transition(const int number_of_from_state,
							      const int number_of_to_state,
							      const int status)
{
    int check = 0;
    
   if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::set_implementation_status_of_special_transition : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	if ((status != 0) && (status != 1))
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_transition : status (" << status
		 << ") must be 0 or 1.\n" << flush;
	    check++;
	}
	if (number_of_from_state < 0)
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") < 0.\n" << flush;
	    check++;
	}
	if (number_of_to_state < 0)
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") < 0.\n" << flush;
	    check++;
	}
	
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }
    
    
    if (check == 0)
    {
	if (((this->get_implementation_status_of_special_transition(number_of_from_state, 
								    number_of_to_state) == 0) &&
	     (status == 0)) ||
	    ((this->get_implementation_status_of_special_transition(number_of_from_state, 
								    number_of_to_state) == 1) &&
	     (status == 1)))
	{ 
	    // there is nothing to do as old_status = new_status
	}
	else 
	{
	    int max_n_of_child_state  = implementation_status_of_special_transitions.GetDimension(0);
	    int linear_index = max_n_of_child_state*number_of_from_state + number_of_to_state;
	    
	    implementation_status_of_special_transitions.SetElement(linear_index, status);
	}
    }
    return(check);
}

int Sequence::set_index_of_special_transition(const int number_of_from_state,
					      const int number_of_to_state,
					      const int index)
{
    int check = 0;
    
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::set_index_of_special_transition : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	if (index != (number_of_implemented_special_transitions-1))
	{
	    cout << "ERROR class Sequence::set_index_of_special_transition : index (" << index 
		 << ") must be = number_of_implemented_special_transitions-1 ("
		 << number_of_implemented_special_transitions-1 << ").\n" << flush;
	    check++;
	}
	if (number_of_from_state < 0)
	{
	    cout << "ERROR class Sequence::set_index_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") < 0.\n" << flush;
	    check++;
	}
	if (number_of_to_state < 0)
	{
	    cout << "ERROR class Sequence::set_index_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") < 0.\n" << flush;
	    check++;
	}
	
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::set_index_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::set_index_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) != 1)
	{
	    cout << "ERROR class Sequence::set_index_of_special_transition : implementation \n"
		 << "status of special transition number_of_from_state (" << number_of_from_state
		 << ") -> number_of_to_state (" << number_of_to_state 
		 << ") must be 1 in order to set the index for it.\n" << flush;
	    check++;
	}
    }
    
    if (check == 0)
    {
	int max_n_of_child_state  = implementation_status_of_special_transitions.GetDimension(0);
	int linear_index = max_n_of_child_state*number_of_from_state + number_of_to_state;
	
	indices_of_special_transitions.SetElement(linear_index, index);
    }
    return(check);
}

int Sequence::set_implementation_status_of_special_emission(const int number_of_state,
							    const int status)
{
    int check = 0;
    
    if (number_of_special_emissions == 0)
    {
	cout << "ERROR class Sequence::set_implementation_status_of_special_emission : this sequence has not \n"
	     << "been set up to implement special emissions.\n" << flush;
	check++;
    }
    else
    {
	if ((status != 0) && (status != 1))
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_emission : status (" << status
		 << ") must be 0 or 1.\n" << flush;
	    check++;
	}
	if (number_of_state < 0)
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_emission : number_of_state ("
		 << number_of_state << ") < 0.\n" << flush;
	    check++;
	}
	
	int max_n_of_child_state = implementation_status_of_special_emissions.GetDimension(0);
	
	if (number_of_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::set_implementation_status_of_special_emission : number_of_state ("
		 << number_of_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }
    if (check == 0)
    {
	if (((this->get_implementation_status_of_special_emission(number_of_state) == 0) && (status == 0)) ||
	    ((this->get_implementation_status_of_special_emission(number_of_state) == 1) && (status == 1)))
	{ 
	    // there is nothing to do as old_status = new_status
	}
	else 
	{
	    implementation_status_of_special_emissions.SetElement(number_of_state, status);
	}
    }
    return(check);
}

int Sequence::set_index_of_special_emission(const int number_of_state,
					    const int index)
{
    int check = 0;
    
    if (number_of_special_emissions == 0)
    {
	cout << "ERROR class Sequence::set_index_of_special_emission : this sequence has not \n"
	     << "been set up to implement special emissions.\n" << flush;
	check++;
    }
    else
    {
	if (index != (number_of_implemented_special_emissions-1))
	{
	    cout << "ERROR class Sequence::set_index_of_special_emission : index (" << index 
		 << ") must be = number_of_implemented_special_emissions-1 ("
		 << number_of_implemented_special_emissions-1 << ").\n" << flush;
	    check++;
	}
	if (number_of_state < 0)
	{
	    cout << "ERROR class Sequence::set_index_of_special_emission : number_of_state ("
		 << number_of_state << ") < 0.\n" << flush;
	    check++;
	}
	
	int max_n_of_child_state = implementation_status_of_special_emissions.GetDimension(0);
	
	if (number_of_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::set_index_of_special_emission : number_of_state ("
		 << number_of_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (this->get_implementation_status_of_special_emission(number_of_state) != 1)
	{
	    cout << "ERROR class Sequence::set_index_of_special_emission : implementation \n"
		 << "status of special emission number_of_state (" << number_of_state
		 << ") must be 1 in order to set the index for it.\n" << flush;
	    check++;
	}
    }
    
    if (check == 0)
    {
	indices_of_special_emissions.SetElement(number_of_state, index);
    }
    return(check);
}

// access functions

int Sequence::get_subsequence(const int start_pos, const int end_pos,
			      Sequence* seq) const
{
    int check = 0;
  
    // check start_pos and end_pos
    
    if (start_pos < start_position) {
	cout << "ERROR : class Sequence::get_subsequence : start_pos (" << start_pos
	     << ") < start_position (" << start_position << ").\n" << flush;
	check++;
    }
    if (end_pos > end_position) {
	cout << "ERROR : class Sequence::get_subsequence : end_pos (" << end_pos
	     << ") > end_position (" << end_position << ").\n" << flush;
	check++;
    }
    if (end_pos < start_pos) {
	cout << "ERROR : class Sequence::get_subsequence : end_pos (" << end_pos
	     << ") < start_pos (" << start_pos << ").\n" << flush;
	check++;
    }
    
    // check this sequence

    if ((sequence == NULL) || (length_of_sequence == 0)){
	cout << "ERROR: class Sequence::get_subsequence: sequence is NULL and length_of_sequence = 0.\n" << flush;
	check++;
    }
    if (this->get_ac() == NULL){
	cout << "ERROR: class Sequence::get_subsequence: ac is NULL.\n" << flush;
	check++;
    }
    if (this->get_start_position() > this->get_end_position()){
	cout << "ERROR: class Sequence::get_subsequence: start_position (" << this->get_start_position()
	     << ") > end_position (" << this->get_end_position() << ").\n" << flush;
	check++;
    }

    if (check == 0){
	// copy this sequence to subsequence

	Sequence sub((*this));
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
      
	// unchanged w.r.t. this sequence:
	//
	// id
	// ac
	// number_of_implemented_special_transitions
	// number_of_special_transitions
	// implementation_status_of_special_transitions
	// indices_of_special_transitions
	// 

	// modify variables of sub

	if(sub.annotation_labels){
	    for(i=0; i<sub.length_of_sequence; i++){
		if(sub.annotation_labels[i]){
		    for(j=0; j<sub.number_of_sources;j++){
			if(sub.annotation_labels[i][j]){
			    for(k=0; k<sub.number_of_annotation_labels; k++){
				if(sub.annotation_labels[i][j][k]) delete[] sub.annotation_labels[i][j][k];
				sub.annotation_labels[i][j][k] = NULL;
			    }
			    delete[] sub.annotation_labels[i][j];
			}
			sub.annotation_labels[i][j] = NULL;
		    }
		    delete[] sub.annotation_labels[i];
		}
		sub.annotation_labels = NULL;
	    }
	    delete [] sub.annotation_labels;
	}
	sub.annotation_labels = NULL;
    
	if(sub.annotation_label_scores){
	    for(i=0; i<sub.length_of_sequence; i++){
		if(sub.annotation_label_scores[i]){
		    for(j=0; j<sub.number_of_sources;j++){
			if(sub.annotation_label_scores[i][j]){
			    for(k=0; k<sub.number_of_annotation_labels; k++){
				if(sub.annotation_label_scores[i][j][k]) delete[] sub.annotation_label_scores[i][j][k];
				sub.annotation_label_scores[i][j][k] = NULL;
			    }
			    delete[] sub.annotation_label_scores[i][j];
			}
			sub.annotation_label_scores[i][j] = NULL;
		    }
		    delete[] sub.annotation_label_scores[i];
		}
		sub.annotation_label_scores = NULL;
	    }
	    delete [] sub.annotation_label_scores;
	}
	sub.annotation_label_scores = NULL;
      
  
	// length_of_sequence
	sub.length_of_sequence = abs(start_pos-end_pos)+1; 
	// start_position
	sub.start_position     = start_pos;                 
	// end_position
	sub.end_position       = end_pos;

	// truncate arrays
	
	int rel_start = 0;  // first index of subsequence in sequence array
	int rel_end   = 0;  // last index of subsequence in sequence array
	
	rel_start = start_pos - start_position;
	rel_end   = end_pos   - start_position;
	
	// sequence
	if (sub.sequence) delete [] sub.sequence;
	sub.sequence = NULL;
	sub.sequence = new int[sub.length_of_sequence];
	for (i=rel_start; i<rel_end+1; i++) {
	    sub.sequence[i-rel_start] = sequence[i];
	}

	// annotation labels
	if(annotation_labels)
	{	    
	    sub.annotation_labels = new bool***[sub.length_of_sequence];
	    for(i=rel_start; i<rel_end+1; i++){
		sub.annotation_labels[i] = new bool**[sub.number_of_sources];
		for(j=0; j<sub.number_of_sources; j++){
		    sub.annotation_labels[i-rel_start][j] = new bool*[sub.number_of_annotation_labels];
		    for(k=0; k<sub.number_of_annotation_labels; k++){
			sub.annotation_labels[i-rel_start][j][k] = new bool[sub.number_of_each_annotation_label[k]];
			for(l=0; l<sub.number_of_each_annotation_label[k]; l++){
			    sub.annotation_labels[i-rel_start][j][k][l] = annotation_labels[i][j][k][l];
			}
		    }
		}
	    }
	    
	    sub.annotation_label_scores = new double***[sub.length_of_sequence];
	    for(i=rel_start; i<rel_end+1; i++){
		sub.annotation_label_scores[i] = new double**[sub.number_of_sources];
		for(j=0; j<sub.number_of_sources; j++){
		    sub.annotation_label_scores[i][j-rel_start] = new double*[sub.number_of_annotation_labels];
		    for(k=0; k<sub.number_of_annotation_labels; k++){
			sub.annotation_label_scores[i][j-rel_start][k] = new double[sub.number_of_each_annotation_label[k]];
			for(l=0; l<sub.number_of_each_annotation_label[k]; l++){
			    sub.annotation_label_scores[i][j-rel_start][k][l] = annotation_label_scores[i][j][k][l];
			}
		    }
		}
	    }
	}
    
	array<int> index_1(1);
	index_1.SetDimension(0, 2);
      
	array<int> index_2(1);
	index_2.SetDimension(0, 2);
	
	// posterior_probs_for_special_transitions
	
	if (posterior_probs_for_special_transitions.GetNumberofDimensions() == 2) {
	    
	    sub.posterior_probs_for_special_transitions.Reset();
	    sub.posterior_probs_for_special_transitions.SetNumberofDimensions(2);
	    sub.posterior_probs_for_special_transitions.SetDimension(0, sub.number_of_special_transitions, 
								     static_cast<Prob>(0.));
	    sub.posterior_probs_for_special_transitions.SetDimension(1, sub.length_of_sequence, 
								     static_cast<Prob>(0.));
	    
	    for (j=0; j<sub.number_of_special_transitions; j++) {
		for (i=rel_start; i<rel_end+1; i++) {
		    
		    index_1.SetElement(0, j); 
		    index_1.SetElement(1, i-rel_start);
		    
		    index_2.SetElement(0, j); 
		    index_2.SetElement(1, i);
		    
		    sub.posterior_probs_for_special_transitions.
			SetElement(index_1,posterior_probs_for_special_transitions.GetElement(index_2));
		}
	    }
	}
	
	// scores_for_special_transitions
	
	if (scores_for_special_transitions.GetNumberofDimensions() == 2) {
	    
	    sub.scores_for_special_transitions.Reset();
	    sub.scores_for_special_transitions.SetNumberofDimensions(2);
	    sub.scores_for_special_transitions.SetDimension(0, sub.number_of_special_transitions, 
							    Logzero);
	    sub.scores_for_special_transitions.SetDimension(1, sub.length_of_sequence, 
							    Logzero);
	    
	    for (j=0; j<sub.number_of_special_transitions; j++) {
		for (i=rel_start; i<rel_end+1; i++) {
		    
		    index_1.SetElement(0, j); 
		    index_1.SetElement(1, i-rel_start);
		    
		    index_2.SetElement(0, j); 
		    index_2.SetElement(1, i);
		    
		    sub.scores_for_special_transitions.
			SetElement(index_1,scores_for_special_transitions.GetElement(index_2));
		}
	    }
	}
	
	// posterior_probs_for_special_emissions
	
	if (posterior_probs_for_special_emissions.GetNumberofDimensions() == 2) {
	    
	    sub.posterior_probs_for_special_emissions.Reset();
	    sub.posterior_probs_for_special_emissions.SetNumberofDimensions(2);
	    sub.posterior_probs_for_special_emissions.SetDimension(0, sub.number_of_special_emissions, 
								   static_cast<Prob>(0.));
	    sub.posterior_probs_for_special_emissions.SetDimension(1, sub.length_of_sequence, 
								   static_cast<Prob>(0.));
	    
	    for (j=0; j<sub.number_of_special_emissions; j++) {
		for (i=rel_start; i<rel_end+1; i++) {
		    
		    index_1.SetElement(0, j); 
		    index_1.SetElement(1, i-rel_start);
		    
		    index_2.SetElement(0, j); 
		    index_2.SetElement(1, i);
		    
		    sub.posterior_probs_for_special_emissions.
			SetElement(index_1,posterior_probs_for_special_emissions.GetElement(index_2));
		}
	    }
	}
	
	// scores_for_special_emissions
	
	if (scores_for_special_emissions.GetNumberofDimensions() == 2) {
	    
	    sub.scores_for_special_emissions.Reset();
	    sub.scores_for_special_emissions.SetNumberofDimensions(2);
	    sub.scores_for_special_emissions.SetDimension(0, sub.number_of_special_emissions, 
							  Logzero);
	    sub.scores_for_special_emissions.SetDimension(1, sub.length_of_sequence, 
							  Logzero);
	    
	    for (j=0; j<sub.number_of_special_emissions; j++) {
		for (i=rel_start; i<rel_end+1; i++) {
		    
		    index_1.SetElement(0, j); 
		    index_1.SetElement(1, i-rel_start);
		    
		    index_2.SetElement(0, j); 
		    index_2.SetElement(1, i);
		    
		    sub.scores_for_special_emissions.
			SetElement(index_1,scores_for_special_emissions.GetElement(index_2));
		}
	    }
	}
	
	// copy sub to seq output
	
	(*seq) = sub;
	
    }
    return(check);
}

int Sequence::get_subsequence(const int start_pos, 
			      const int end_pos, 
			      int** subsequence)
{
    int check = 0;

    if ((*subsequence) != NULL) {
	cout << "ERROR : class Sequence::get_subsequence : subsequence is not NULL.\n" << flush;
	check++;
    }
    
    if (length_of_sequence == 0) {
	cout << "ERROR : class Sequence::get_subsequence : length_of_sequence is 0.\n" << flush;
	check++;
    }

    if (start_position > end_position) {
	cout << "ERROR : class Sequence::get_subsequence : start_position (" << start_position
	     << ") > end_position (" << end_position << ").\n" << flush;
	check++;
    }

    // check input positions

    if (start_pos < start_position) {
	cout << "ERROR : class Sequence::get_subsequence : start_pos (" << start_pos
	 << ") < start_position (" << start_position << ").\n" << flush;
	check++;
    }
    if (end_pos > end_position) {
	cout << "ERROR : class Sequence::get_subsequence : end_pos (" << end_pos
	     << ") > end_position (" << end_position << ").\n" << flush;
	check++;
    }
    if (end_pos < start_pos) {
	cout << "ERROR : class Sequence::get_subsequence : end_pos (" << end_pos
	     << ") < start_pos (" << start_pos << ").\n" << flush;
	check++;
    }
    
    if (check == 0) {
	
	(*subsequence) = new int[end_pos-start_pos+1];

	int i = 0;
	int rel_start = 0;  // first index of subsequence in sequence array
	int rel_end   = 0;  // last index of subsequence in sequence array
	
	rel_start = start_pos - start_position;
	rel_end   = end_pos   - start_position;

	for (i=rel_start; i<rel_end+1; i++) {
	    (*subsequence)[i-rel_start] = sequence[i];
	}
    }
    return(check);
}

int Sequence::get_subsequence_rel(const int rel_start_pos, const int rel_end_pos,
				  Sequence* seq) const
{
    int check = 0;
  
    // check start_pos and end_pos
    
    if (rel_start_pos < 0) {
	cout << "ERROR : class Sequence::get_subsequence_rel : rel_start_pos (" << rel_start_pos
	     << ") < 0.\n" << flush;
	check++;
    }
    if (rel_end_pos > (length_of_sequence-1)) {
	cout << "ERROR : class Sequence::get_subsequence_rel : rel_end_pos (" << rel_end_pos
	     << ") > length_of_sequence-1 (" << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    if (rel_start_pos > rel_end_pos) {
	cout << "ERROR : class Sequence::get_subsequence_rel : rel_start_pos (" << rel_start_pos
	     << ") > rel_end_pos (" << rel_end_pos << ").\n" << flush;
	check++;
    }
    
    // check this sequence

    if ((sequence == NULL) || (length_of_sequence == 0)){
	cout << "ERROR: class Sequence::get_subsequence_rel: sequence is NULL and length_of_sequence = 0.\n" << flush;
	check++;
    }
    if (this->get_ac() == NULL){
	cout << "ERROR: class Sequence::get_subsequence_rel: ac is NULL.\n" << flush;
	check++;
    }
    if (this->get_start_position() > this->get_end_position()){
	cout << "ERROR: class Sequence::get_subsequence_rel: start_position (" << this->get_start_position()
	     << ") > end_position (" << this->get_end_position() << ").\n" << flush;
	check++;
    }
    
    if (check == 0) {
	
	Sequence sub;
	
	int start_pos = 0;
	int end_pos   = 0;

	start_pos = this->start_position + rel_start_pos;
	end_pos   = this->start_position + rel_end_pos;
	
	check += this->get_subsequence(start_pos, end_pos, &sub);
	
	if (check == 0) {
	    (*seq) = sub;
	}
    }
    return(check);
}

int Sequence::get_subsequence_rel(const int rel_start_pos, 
				  const int rel_end_pos, 
				  int** subsequence)
{
    int check = 0;
  
    if ((*subsequence) != NULL) {
	cout << "ERROR : class Sequence::get_subsequence_rel : subsequence is not NULL.\n" << flush;
	check++;
    }
    
    if (length_of_sequence == 0) {
	cout << "ERROR : class Sequence::get_subsequence_rel : length_of_sequence is 0.\n" << flush;
	check++;
    }

    if (start_position > end_position) {
	cout << "ERROR : class Sequence::get_subsequence_rel : start_position (" << start_position
	     << ") > end_position (" << end_position << ").\n" << flush;
	check++;
    }

    // check input positions

    if (rel_start_pos < 0) {
	cout << "ERROR : class Sequence::get_subsequence_rel : rel_start_pos (" << rel_start_pos
	     << ") < 0.\n" << flush;
	check++;
    }
    
    if (rel_end_pos > (length_of_sequence-1)) {
	cout << "ERROR : class Sequence::get_subsequence_rel : rel_end_pos (" << rel_end_pos
	     << ") > length_of_sequence-1 (" << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    
    if (rel_end_pos < rel_start_pos) {
	cout << "ERROR : class Sequence::get_subsequence_rel : rel_end_pos (" << rel_end_pos
	     << ") < rel_start_pos (" << rel_start_pos << ").\n" << flush;
	check++;
    }
  
    if (check == 0) {

	(*subsequence) = new int[rel_end_pos-rel_start_pos+1];
    
	int i = 0;
	
	for (i=rel_start_pos; i<rel_end_pos+1; i++) {
	    (*subsequence)[i-rel_start_pos] = sequence[i];
	}
    }
    return(check);
}

char* Sequence::get_ac(void) const 
{
    return(ac);
}

int Sequence::get_start_position(void) const
{
    return(start_position);
}

int Sequence::get_end_position(void) const
{
    return(end_position);
}

int Sequence::get_number_of_sources(void) const
{
    return(number_of_sources);
}

int Sequence::get_number_of_annotation_labels(void) const
{
    return(number_of_annotation_labels);
}

int* Sequence::get_number_of_each_annotation_label(void) const
{
    return(number_of_each_annotation_label);
}

int* Sequence::get_sequence(void) const
{
    return(sequence);
}

char* Sequence::get_char_sequence(model_parameters* const MP) const
{
  // note : calling program is responsible for deleting memory 

    char* return_sequence = NULL;

    if (length_of_sequence > 0) {
     
	return_sequence = new char[length_of_sequence+1];
    
	int i = 0;
	for (i=0; i<length_of_sequence; i++) {
	    return_sequence[i] = MP->get_Alphabet_name(sequence[i]);
	}
	return_sequence[length_of_sequence] = '\0';
    }
    return(return_sequence);
}

int Sequence::length(void) const
{
    return(length_of_sequence);
}

int Sequence::letter(const int i) const {
    int check = 0;
    if ((i>length_of_sequence-1) || (i<0)) {
	cout << "ERROR class Sequence::letter : index i (" << i << ") out of range (0.." 
	     << length_of_sequence-1 << ")\n" << flush;
	check++;
    }
    return(sequence[i]);
}

const char* Sequence::get_sequence_type(void) const
{
    return(convert_int_to_sequence_type[sequence_type]);
}

int Sequence::get_sequence_type_int(void) const  
{
    return(sequence_type);
}

int Sequence::get_number_of_special_transitions(void) const
{
    return(number_of_special_transitions);
}

int Sequence::get_number_of_implemented_special_transitions(void) const
{
    return(number_of_implemented_special_transitions);
}

int Sequence::get_index_of_special_transition(const int number_of_from_state,
					      const int number_of_to_state) const
{
    int check = 0;

    if (number_of_from_state < 0)
    {
	cout << "ERROR class Sequence::get_index_of_special_transition : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0)
    {
	cout << "ERROR class Sequence::get_index_of_special_transition : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::get_index_of_special_transition : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::get_index_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::get_index_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }


    if (check == 0)
    {
	if (number_of_implemented_special_transitions == 0)
	{
	    return(0);
	}
	else
	{
	    int max_n_of_child_state  = implementation_status_of_special_transitions.GetDimension(0);
	    int linear_index = max_n_of_child_state*number_of_from_state + number_of_to_state;
	    
	    return(indices_of_special_transitions.GetElement(linear_index));
	}
    }
    else
    {
	return(-1);
    }
}

int Sequence::get_implementation_status_of_special_transition(const int number_of_from_state,
							      const int number_of_to_state) const
{

    int check = 0;
    if (number_of_from_state < 0) {
	cout << "ERROR class Sequence::get_implementation_status_of_special_transition : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0) {
	cout << "ERROR class Sequence::get_implementation_status_of_special_transition : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0) {
	cout << "ERROR class Sequence::get_implementation_status_of_special_transition : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_implementation_status_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (number_of_to_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_implementation_status_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }
    if (check == 0) {
	if (number_of_implemented_special_transitions == 0)
	{
	    return(0);
	}
	else {
	    int max_n_of_child_state  = implementation_status_of_special_transitions.GetDimension(0);
	    int linear_index          = max_n_of_child_state*number_of_from_state + number_of_to_state;
	    
	    return(implementation_status_of_special_transitions.GetElement(linear_index));
	}
    }
    else {
	return(-1);
    }
}

Prob Sequence::get_prior_of_special_transition(const int number_of_from_state,
					       const int number_of_to_state) const
{

    int check = 0;
    
    if ((posterior_probs_for_special_transitions.GetNumberofDimensions() == 2) &&
	(scores_for_special_transitions.GetNumberofDimensions() == 0)) {
	
	cout << "ERROR: Sequence::get_prior_of_special_transition: sequence object is ment to not be used with scores, \n"
	     << "but with posterior probs.\n" << flush;
	check++;
    }
    if (number_of_from_state < 0)
    {
	cout << "ERROR class Sequence::get_prior_of_special_transition : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0)
    {
	cout << "ERROR class Sequence::get_prior_of_special_transition : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::get_prior_of_special_transition : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::get_prior_of_special_transition : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::get_prior_of_special_transition : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }
    
    
    if (check == 0)
    {
	if (number_of_implemented_special_transitions == 0)
	{
	    return(0);
	}
	else
	{
	    
	    const int index = this->get_index_of_special_transition(number_of_from_state,
								    number_of_to_state);
	    return(priors_of_special_transitions.GetElement(index));
	    
	}
    }
    else
    {
	return(-1.0);
    }
}

Prob Sequence::get_posterior_prob(const int number_of_from_state,
				  const int number_of_to_state,
				  const int position_in_sequence) const
{

    int check = 0;
    
    if (number_of_from_state < 0) {
	cout << "ERROR class Sequence::get_posterior_prob : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0) {
	cout << "ERROR class Sequence::get_posterior_prob : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0) {
	cout << "ERROR class Sequence::get_posterior_prob : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_posterior_prob : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (number_of_to_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_posterior_prob : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) == 0) {
	    cout << "ERROR class Sequence::get_posterior_prob : special transition number_of_from_state ("
		 << number_of_from_state << ") -> number_of_to_state (" << number_of_to_state
		 << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }
    if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1))) {
	cout << "ERROR class Sequence::get_posterior_prob : position_in_sequence ("
	     << position_in_sequence << ") out of range (0.." 
	     << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    if (check == 0) {
	int index_of_special_transition = this->get_index_of_special_transition(number_of_from_state,
										number_of_to_state);
	int linear_index = index_of_special_transition * length_of_sequence + position_in_sequence;
	Prob posterior_prob = posterior_probs_for_special_transitions.GetElement(linear_index);
	
	return(posterior_prob);
    }
    else {
	return(static_cast<Prob>(0.0));
    }
}

Score Sequence::get_score(const int number_of_from_state,
			 const int number_of_to_state,
			 const int position_in_sequence) const
{

    int check = 0;

    if ((posterior_probs_for_special_transitions.GetNumberofDimensions() == 2) &&
	(scores_for_special_transitions.GetNumberofDimensions() == 0))
    {
	cout << "ERROR: Sequence::get_score: sequence object is not ment to be used with scores, but \n"
	     << "with posterior probs.\n" << flush;
	check++;
    }
    if (number_of_from_state < 0) 
    {
	cout << "ERROR class Sequence::get_score : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0) 
    {
	cout << "ERROR class Sequence::get_score : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0) 
    {
	cout << "ERROR class Sequence::get_score : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1)) 
	{
	    cout << "ERROR class Sequence::get_score : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (number_of_to_state > (max_n_of_child_state-1)) 
	{
	    cout << "ERROR class Sequence::get_score : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) == 0) 
	{
	    cout << "ERROR class Sequence::get_score : special transition number_of_from_state ("
		 << number_of_from_state << ") -> number_of_to_state (" << number_of_to_state
		 << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }
    if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1))) 
    {
	cout << "ERROR class Sequence::get_score : position_in_sequence ("
	     << position_in_sequence << ") out of range (0.." 
	     << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    if (check == 0) 
    {
	int index_of_special_transition = this->get_index_of_special_transition(number_of_from_state,
										number_of_to_state);
	int linear_index   = index_of_special_transition * length_of_sequence + position_in_sequence;
	Score score        = scores_for_special_transitions.GetElement(linear_index);
	
	return(score);
    }
    else {
	return(0);
    }
}

int Sequence::get_number_of_special_emissions(void) const
{
    return(number_of_special_emissions);
}

int Sequence::get_number_of_implemented_special_emissions(void) const
{
    return(number_of_implemented_special_emissions);
}

int Sequence::get_index_of_special_emission(const int number_of_state) const
{
    int check = 0;
    
    if (number_of_state < 0) {
	cout << "ERROR class Sequence::get_index_of_special_emission : number_of_state ("
	     << number_of_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_emissions == 0) {
	cout << "ERROR class Sequence::get_index_of_special_emission : this sequence has not \n"
	     << "been set up to implement special emissions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_emissions.GetDimension(0);
	
	if (number_of_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_index_of_special_emission : number_of_state ("
		 << number_of_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }
    
    if (check == 0) {
	
	if (number_of_implemented_special_emissions == 0)
	{
	    return(0);
	}
	else
	{
	    return(indices_of_special_emissions.GetElement(number_of_state));
	}
    }
    else {
	return(-1);
    }
}

int Sequence::get_implementation_status_of_special_emission(const int number_of_state) const
{
    int check = 0;
    if (number_of_state < 0) {
	cout << "ERROR class Sequence::get_implementation_status_of_special_emission : number_of_state ("
	     << number_of_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_emissions == 0) {
	cout << "ERROR class Sequence::get_implementation_status_of_special_emission : this sequence has not \n"
	     << "been set up to implement special emissions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_emissions.GetDimension(0);
	
	if (number_of_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_implementation_status_of_special_emission : number_of_state ("
		 << number_of_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
    }
    if (check == 0) {
	if (number_of_implemented_special_emissions == 0)
	{
	    return(0);
	}
	else {
	    return(implementation_status_of_special_emissions.GetElement(number_of_state));
	}
    }
    else {
	return(-1);
    }
}


Prob Sequence::get_posterior_prob(const int number_of_state,
				  const int position_in_sequence) const
{

    int check = 0;
    
    if (number_of_state < 0) {
	cout << "ERROR class Sequence::get_posterior_prob : number_of_state ("
	     << number_of_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_emissions == 0) {
	cout << "ERROR class Sequence::get_posterior_prob : this sequence has not \n"
	     << "been set up to implement special emissions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_emissions.GetDimension(0);
	
	if (number_of_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_posterior_prob : number_of_state ("
		 << number_of_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (this->get_implementation_status_of_special_emission(number_of_state) == 0) {
	    cout << "ERROR class Sequence::get_posterior_prob : special emission number_of_from_state ("
		 << number_of_state << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }
    if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1))) {
	cout << "ERROR class Sequence::get_posterior_prob : position_in_sequence ("
	     << position_in_sequence << ") out of range (0.." 
	     << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    if (check == 0) {

	int index_of_special_emission = this->get_index_of_special_emission(number_of_state);
	int linear_index = index_of_special_emission * length_of_sequence + position_in_sequence;
	Prob posterior_prob = posterior_probs_for_special_emissions.GetElement(linear_index);
	
	return(posterior_prob);

    }
    else {
	return(static_cast<Prob>(0.0));
    }
}

Prob Sequence::get_score(const int number_of_state,
			 const int position_in_sequence) const
{

    int check = 0;
    
    if (number_of_state < 0) {
	cout << "ERROR class Sequence::get_score : number_of_state ("
	     << number_of_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_emissions == 0) {
	cout << "ERROR class Sequence::get_score : this sequence has not \n"
	     << "been set up to implement special emissions.\n" << flush;
	check++;
    }
    else {
	int max_n_of_child_state = implementation_status_of_special_emissions.GetDimension(0);
	
	if (number_of_state > (max_n_of_child_state-1)) {
	    cout << "ERROR class Sequence::get_score : number_of_state ("
		 << number_of_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	if (this->get_implementation_status_of_special_emission(number_of_state) == 0) {
	    cout << "ERROR class Sequence::get_score : special emission number_of_state ("
		 << number_of_state<< ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }
    if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1))) {
	cout << "ERROR class Sequence::get_score : position_in_sequence ("
	 << position_in_sequence << ") out of range (0.." 
	     << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    if (check == 0) {

	Score score = Logzero;
	
	if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1))) {
	    
	    cout << "NOTE class Sequence::get_score : position_in_sequence ("
		 << position_in_sequence << ") out of range (0.." 
		 << length_of_sequence-1 << ") => set score to Logzero.\n" << flush;
	    
	    score = Logzero;
	}
	else {
	    
	    int index_of_special_emission = this->get_index_of_special_emission(number_of_state);
	    int linear_index   = index_of_special_emission * length_of_sequence + position_in_sequence;
	    
	    score = scores_for_special_emissions.GetElement(linear_index);
	   
	}

	return(score);
	
    }
    else {
	return(0); 
    }
}

int Sequence::get_next_annotation_labels_rel(// input
    const int        rel_position, 
    bool***    const label,
    double***  const label_score,
    // output
    int*             next_rel_position,
    bool****   const next_label,
    double**** const next_label_score) const
{
    // note: this function will return the input values if no next_label is different
    //       from the input label
    
    int check = 0;
  
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_next_annotation_labels_rel: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_next_annotation_labels_rel: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }

    if ((rel_position < 0) || (rel_position > (length_of_sequence-1))) 
    {
	cout << "ERROR: class Sequence::get_next_annotation_labels_rel: rel_position ("
	     << rel_position << ") out of range [0," << length_of_sequence-1 << "].\n" << flush;
	check++;
    }
 
    if (check == 0) 
    {
	
	// initialise output variables and set them to input variables		
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	
	for(i=0; i<number_of_sources; i++){
	    for(j=0; j<number_of_annotation_labels; j++){
		for(k=0; k<number_of_each_annotation_label[j]; k++){
		    (*next_label)[i][j][k]        = label[i][j][k];
		    (*next_label_score)[i][j][k]  = label_score[i][j][k];
		}
	    }
	}

	(*next_rel_position) = rel_position;
	bool found = false;
	
	for(i=rel_position; ((i<length_of_sequence)&&(!found)); i++){
	    for (j=0; ((j<number_of_sources) && (!found)); j++) {
		for(k=0; ((k<number_of_annotation_labels)&&(!found)); k++){
		    for(l=0; ((l<number_of_each_annotation_label[k])&&(!found));l++){
			if ((label[j][k][l]!=annotation_labels[i][j][k][l])||
			    (label_score[j][k][l] != annotation_label_scores[i][j][k][l]))
			{
			    found = true;
			    (*next_rel_position) = i;
			    break;
			}
		    }
		}
	    }
	}
	// assign value for next_label
	for(i=0; i<number_of_sources; i++){
	    for(j=0; j<number_of_annotation_labels; j++){
		for(k=0; k<number_of_each_annotation_label[j];k++){
		    (*next_label)[i][j][k] = annotation_labels[(*next_rel_position)][i][j][k];
		    (*next_label_score)[i][j][k] = annotation_label_scores[(*next_rel_position)][i][j][k];
		}
	    }
	}
    }
    return(check);
}

int Sequence::get_annotation_labels(//output
    bool*****   const label,
    double***** const label_score) const
    
{
    int check = 0;
    
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }

    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
       
    if (check == 0) {
	
	if(!(*label))
	{
	    (*label)  = new bool***[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label)[i] = new bool**[number_of_sources];
		for(j=0; j<number_of_sources; j++)
		{
		    (*label)[i][j] = new bool*[number_of_annotation_labels];
		    for(k=0; k<number_of_annotation_labels; k++)
		    {
			(*label)[i][j][k] = new bool[number_of_each_annotation_label[k]];	
			for(l=0; l<number_of_each_annotation_label[k]; l++)
			{
			    (*label)[i][j][k][l] = false;		
			}
		    }
		}
	    }
	}
	
	if(!(*label_score))
	{
	    (*label_score) = new double***[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label_score)[i] = new double**[number_of_sources];
		for(j=0; j<number_of_sources; j++)
		{
		    (*label_score)[i][j] = new double*[number_of_annotation_labels];
		    for(k=0; k<number_of_annotation_labels; k++)
		    {
			(*label_score)[i][j][k] = new double[number_of_each_annotation_label[k]];
			for(l=0; l<number_of_each_annotation_label[k]; l++)
			{
			    (*label_score)[i][j][k][l] = 0;
			}
		    }
		}
	    }
	}

	for(i=0; i<length_of_sequence; i++)
	{
	    for(j=0; j<number_of_sources; j++)
	    {
		for(k=0; k<number_of_annotation_labels; k++)
		{
		    for(l=0; l<number_of_each_annotation_label[j]; l++)
		    {
			(*label)[i][j][k][l] = annotation_labels[i][j][k][l];
			(*label_score)[i][j][k][l] = annotation_label_scores[i][j][k][l];
		    }
		}
	    }
	}
    }

    // release memory for label and label scores if there is error
    if(check)
    {
	if(*label)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label)[i])
		{
		    for(j=0; j<number_of_sources; j++)
		    {
			if((*label)[i][j]){
			    for(k=0; k<number_of_annotation_labels; k++)
			    {
				if((*label)[i][j][k]) delete[] (*label)[i][j][k];
				(*label)[i][j][k] = NULL;
			    }
			    delete [] (*label)[i][j];
			}			
			(*label)[i][j] = NULL;
		    }
		    delete [] (*label)[i];
		}
		(*label)[i] = NULL;
	    }
	    delete [] (*label);
	}
	(*label) = NULL;

	if(*label_score)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label_score)[i])
		{
		    for(j=0; j<number_of_sources; j++)
		    {
			if((*label_score)[i][j]){ 
			    for(k=0; k<number_of_annotation_labels; k++)
			    {
				if((*label_score)[i][j][k]) delete[] (*label_score)[i][j][k];
				(*label_score)[i][j][k] = NULL;
			    }			 
			    delete [] (*label_score)[i][j];
			}
			(*label_score)[i][j] = NULL;
		    }
		    delete [] (*label_score)[i];
		}
		(*label_score)[i] = NULL;
	    }
	    delete [] (*label_score);
	}
	(*label_score) = NULL;

    }

    return(check);
}

int Sequence::get_combine_annotation_labels(bool****   const label,
					    double**** const label_score)
{
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_combine_annotation_labels: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_combine_annotation_labels: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }
    
    if(!check)
    {	
	if(!(*label))
	{
	    (*label) = new bool**[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label)[i] = new bool*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label)[i][j] = new bool[number_of_each_annotation_label[j]];	
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label)[i][j][k] = false;		
		    }
		}
	    }
	}

	if(!(*label_score))
	{
	    (*label_score) = new double**[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label_score)[i] = new double*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label_score)[i][j] = new double[number_of_each_annotation_label[j]];
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label_score)[i][j][k] = 0;
		    }
		}
	    }
	}

	bool* label_set = new bool[number_of_annotation_labels];

	for(i=0; i<length_of_sequence; i++)
	{
	    if(check)
	    {
		break;
	    }
	    for(j=0; j<number_of_annotation_labels; j++)
	    {
		if(check)
		{
		    break;
		}
		label_set[j] = false;

		if(!score_labels[j])
		{
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			if(check)
			{
			    break;
			}
			for(l=0; l<number_of_sources; l++)
			{
			    if(check)
			    {
				break;
			    }
			    if(annotation_labels[i][l][j][k])
			    {			
				(*label)[i][j][k] = annotation_labels[i][l][j][k];
				(*label_score)[i][j][k] = annotation_label_scores[i][l][j][k];
			    }
			}
		    }							  
		}
		else{
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			if(check)
			{
			    break;
			}
			for(l=0; l<number_of_sources; l++)
			{
			    if(check)
			    {
				break;
			    }
			    if(annotation_labels[i][l][j][k])
			    {			
				if(!(*label)[i][j][k])
				{
				    (*label)[i][j][k] = annotation_labels[i][l][j][k];
				    (*label_score)[i][j][k] = annotation_label_scores[i][l][j][k];
				}else{
				    if((*label_score)[i][j][k]!=annotation_label_scores[i][l][j][k])
				    {
					cout << "ERROR: class Sequence::get_combine_annotation_labels: "
					     <<"contradicting label("<<k<<") of label set("<<j<<") "
					     <<"from source : "<<l<<" at position "<<i<<endl;
					
					cout <<"original score : "<<(*label_score)[i][j][k]<<", new score : "<<annotation_label_scores[i][l][j][k]<<endl;
					check++;
					break;
				    }
				}
			    }
			}
		    }		    
		}
		
	    }
	}

	if(label_set) delete [] label_set;
	label_set = NULL;

	
    }

    // release memory for label and label scores if there is error
    if(check)
    {
	if(*label)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label)[i][j]) delete [] (*label)[i][j];
			(*label)[i][j] = NULL;
		    }
		    delete [] (*label)[i];
		}
		(*label)[i] = NULL;
	    }
	    delete [] (*label);
	}
	(*label) = NULL;

	if(*label_score)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label_score)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label_score)[i][j]) delete [] (*label_score)[i][j];
			(*label_score)[i][j] = NULL;
		    }
		    delete [] (*label_score)[i];
		}
		(*label_score)[i] = NULL;
	    }
	    delete [] (*label_score);
	}
	(*label_score) = NULL;

    }
	    
    return check;
}


int Sequence::get_annotation_labels_for_position(const int abs_position, 
						 //output
						 bool****   const label,
						 double**** const label_score) const
    
{
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }

    if ((abs_position < start_position) || (abs_position > end_position)) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position: abs_position ("
	     << abs_position << ") out of range [" << start_position << ","
	     << end_position << "].\n" << flush;
	check++;
    }
   
    if (check == 0) {
	
	int rel_position = 0;
	rel_position = abs_position - start_position;

	if(!(*label))
	{
	    (*label) = new bool**[number_of_sources];
	    for(i=0; i<number_of_sources; i++)
	    {
		(*label)[i] = new bool*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label)[i][j] = new bool[number_of_each_annotation_label[j]];	
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label)[i][j][k] = false;		
		    }
		}
	    }
	}

	if(!(*label_score))
	{
	    (*label_score) = new double**[number_of_sources];
	    for(i=0; i<number_of_sources; i++)
	    {
		(*label_score)[i] = new double*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label_score)[i][j] = new double[number_of_each_annotation_label[j]];
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label_score)[i][j][k] = 0;
		    }
		}
	    }
	}

	for(i=0; i<number_of_sources; i++)
	{
	    for(j=0; j<number_of_annotation_labels; j++)
	    {
		for(k=0; k<number_of_each_annotation_label[j]; k++)
		{
		    (*label)[i][j][k] = annotation_labels[rel_position][i][j][k];
		    (*label_score)[i][j][k] = annotation_label_scores[rel_position][i][j][k];
		}
	    }
	}

    }

    // release memory for label and label scores if there is error
    if(check)
    {
	if(*label)
	{
	    for(i=0; i<number_of_sources; i++)
	    {
		if((*label)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label)[i][j]) delete [] (*label)[i][j];
			(*label)[i][j] = NULL;
		    }
		    delete [] (*label)[i];
		}
		(*label)[i] = NULL;
	    }
	    delete [] (*label);
	}
	(*label) = NULL;

	if(*label_score)
	{
	    for(i=0; i<number_of_sources; i++)
	    {
		if((*label_score)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label_score)[i][j]) delete [] (*label_score)[i][j];
			(*label_score)[i][j] = NULL;
		    }
		    delete [] (*label_score)[i];
		}
		(*label_score)[i] = NULL;
	    }
	    delete [] (*label_score);
	}
	(*label_score) = NULL;

    }

    return(check);
}

int Sequence::get_annotation_labels_for_position_rel(const int rel_position, 
						     // output
						     bool****   const label,
						     double**** const label_score) const
{
    int check = 0;

    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position_rel: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position_rel: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }

    if ((rel_position < 0) || (rel_position > length_of_sequence-1)) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_position_rel: rel_position ("
	     << rel_position << ") out of range [0," << length_of_sequence-1 << "].\n" << flush;
	check++;
    }

    int i = 0;
    int j = 0;
    int k = 0;


    if (check == 0) 
    {
      
	if(!(*label))
	{
	    (*label) = new bool**[number_of_sources];
	    for(i=0; i<number_of_sources; i++)
	    {
		(*label)[i] = new bool*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label)[i][j] = new bool[number_of_each_annotation_label[j]];	
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label)[i][j][k] = false;		
		    }
		}
	    }
	}

	if(!(*label_score))
	{
	    (*label_score) = new double**[number_of_sources];
	    for(i=0; i<number_of_sources; i++)
	    {
		(*label_score)[i] = new double*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label_score)[i][j] = new double[number_of_each_annotation_label[j]];
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label_score)[i][j][k] = 0;
		    }
		}
	    }
	}
       
	for(i=0; i<number_of_sources; i++)
	{
	    for(j=0; j<number_of_annotation_labels; j++)
	    {
		for(k=0; k<number_of_each_annotation_label[j]; k++)
		{
		    (*label)[i][j][k] = annotation_labels[rel_position][i][j][k];
		    (*label_score)[i][j][k] = annotation_label_scores[rel_position][i][j][k];
		}
	    }
	}
    }

    // release memory for label and label scores if there is error
    if(check)
    {
	if(*label)
	{
	    for(i=0; i<number_of_sources; i++)
	    {
		if((*label)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label)[i][j]) delete [] (*label)[i][j];
			(*label)[i][j] = NULL;
		    }
		    delete [] (*label)[i];
		}
		(*label)[i] = NULL;
	    }
	    delete [] (*label);
	}
	(*label) = NULL;

	if(*label_score)
	{
	    for(i=0; i<number_of_sources; i++)
	    {
		if((*label_score)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label_score)[i][j]) delete [] (*label_score)[i][j];
			(*label_score)[i][j] = NULL;
		    }
		    delete [] (*label_score)[i];
		}
		(*label_score)[i] = NULL;
	    }
	    delete [] (*label_score);
	}
	(*label_score) = NULL;
    }
    
    return(check);
}

int Sequence::get_annotation_labels_for_set(const int        label_set, 
					    bool****   const label,
					    double**** const label_score) const
    		
{
    int check = 0;

    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_set: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_set: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }
    
    if((label_set<0)||(label_set>=number_of_annotation_labels))
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_set: label_set("
	     << label_set << ") out of range [" << number_of_annotation_labels-1<<".\n";
	check++;
    }
    	
    int i = 0;
    int j = 0;
    int k = 0;
	
    if (check == 0) {

	if(!(*label))
	{
	    (*label) = new bool**[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label)[i] = new bool*[number_of_sources];
		for(j=0; j<number_of_sources; j++)
		{
		    (*label)[i][j] = new bool[number_of_each_annotation_label[label_set]];	
		    for(k=0; k<number_of_each_annotation_label[label_set]; k++)
		    {
			(*label)[i][j][k] = false;		
		    }
		}
	    }
	}

	if(!(*label_score))
	{
	    (*label_score) = new double**[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label_score)[i] = new double*[number_of_sources];
		for(j=0; j<number_of_sources; j++)
		{
		    (*label_score)[i][j] = new double[number_of_each_annotation_label[label_set]];
		    for(k=0; k<number_of_each_annotation_label[label_set]; k++)
		    {
			(*label_score)[i][j][k] = 0;
		    }
		}
	    }
	}
	for(i=0; i<length_of_sequence; i++){
	    for(j=0; j<number_of_sources; j++){
		for(k=0; k<number_of_each_annotation_label[label_set]; k++){
		    (*label)[i][j][k] = annotation_labels[i][j][label_set][k];
		    (*label_score)[i][j][k] = annotation_label_scores[i][j][label_set][k];
		}
	    }
	}       
    }

    // release memory for label and label scores if there is error
    if(check)
    {
	if(*label)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label)[i])
		{
		    for(j=0; j<number_of_sources; j++)
		    {
			if((*label)[i][j]) delete [] (*label)[i][j];
			(*label)[i][j] = NULL;
		    }
		    delete [] (*label)[i];
		}
		(*label)[i] = NULL;
	    }
	    delete [] (*label);
	}
	(*label) = NULL;

	if(*label_score)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label_score)[i])
		{
		    for(j=0; j<number_of_sources; j++)
		    {
			if((*label_score)[i][j]) delete [] (*label_score)[i][j];
			(*label_score)[i][j] = NULL;
		    }
		    delete [] (*label_score)[i];
		}
		(*label_score)[i] = NULL;
	    }
	    delete [] (*label_score);
	}
	(*label_score) = NULL;
    }

    return(check);
}

int Sequence::get_annotation_labels_for_source(const int        source,
		                               bool****   const label,
					       double**** const label_score) const
{
    int check = 0;

    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_source: annotation_labels array is NULL.\n" << flush;
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_source: annotation_label_scores array is NULL.\n" << flush;
	check++;
    }
    
    if((source<0)||(source>=number_of_sources))
    {
	cout << "ERROR: class Sequence::get_annotation_labels_for_source: source("
	     << source << ") out of range [" << number_of_sources-1<<".\n";
	check++;
    }

    int i = 0;
    int j = 0;
    int k = 0;

    if (check == 0) {
	
	if(!(*label))
	{
	    (*label) = new bool**[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label)[i] = new bool*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label)[i][j] = new bool[number_of_each_annotation_label[j]];	
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label)[i][j][k] = false;		
		    }
		}
	    }
	}

	if(!(*label_score))
	{
	    (*label_score) = new double**[length_of_sequence];
	    for(i=0; i<length_of_sequence; i++)
	    {
		(*label_score)[i] = new double*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++)
		{
		    (*label_score)[i][j] = new double[number_of_each_annotation_label[j]];
		    for(k=0; k<number_of_each_annotation_label[j]; k++)
		    {
			(*label_score)[i][j][k] = 0;
		    }
		}
	    }
	}
	
	for(i=0; i<length_of_sequence; i++){
	    for(j=0; j<number_of_annotation_labels; j++){
		for(k=0; k<number_of_each_annotation_label[j]; k++){
		    (*label)[i][j][k] = annotation_labels[i][source][j][k];
		    (*label_score)[i][j][k] = annotation_label_scores[i][source][j][k];
		}
	    }
	}     
    }

    // release memory for label and label scores if there is error
    if(check)
    {
	if(*label)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label)[i][j]) delete [] (*label)[i][j];
			(*label)[i][j] = NULL;
		    }
		    delete [] (*label)[i];
		}
		(*label)[i] = NULL;
	    }
	    delete [] (*label);
	}
	(*label) = NULL;

	if(*label_score)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if((*label_score)[i])
		{
		    for(j=0; j<number_of_annotation_labels; j++)
		    {
			if((*label_score)[i][j]) delete [] (*label_score)[i][j];
			(*label_score)[i][j] = NULL;
		    }
		    delete [] (*label_score)[i];
		}
		(*label_score)[i] = NULL;
	    }
	    delete [] (*label_score);
	}
	(*label_score) = NULL;
    }

    return(check);    
    
}

bool Sequence:: get_score_labels(const int i) const
{   
    if((i<0)||(i>=number_of_annotation_labels))
    {
	cout<<"Error: Sequence class:: get_score_labels "
	    <<"index("<<i<<") out of range[0.."
	    <<number_of_annotation_labels<<"].\n";
	return false;
    }
    if(!score_labels)
    {
	cout<<"Error: Sequence class:: get_score_labels "
	    <<"score_labels array is NULL.\n";
	return false;
    }
    return score_labels[i];
}

int Sequence:: get_constraint() const
{
    return constraint;
}

int Sequence::print_non_alphabet_contents_of_sequence(model_parameters* const MP,
						      std::ostream &o) const
{
    int check = 0;

    if (this->sequence == NULL) {
	cout << "ERROR: class Sequence::print_non_alphabet_contents_of_sequence: sequence array is NULL.\n" << flush;
	check++;
    }
    if (this->annotation_labels == NULL) {
	cout << "ERROR: class Sequence::print_non_alphabet_contents_of_sequence: labels array is NULL.\n" << flush;
	check++;
    }
    if (check == 0) {
	
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;
	int size = MP->get_Alphabet_size();

	for (i=0; i<length_of_sequence; i++) {
	    if ((sequence[i]<0)||(sequence[i]>=size)){
		int range_j = min(length_of_sequence,i+3);
		for (j= max(0,i-2); j<range_j; j++) {
		    o << "sequence[" << j << "] (" << MP->get_Alphabet_name(sequence[j]) << ")"; 
		    for(k=0; k<number_of_sources; k++){
			for(l=0;l<number_of_annotation_labels;l++){
			    o<<"( "<<MP->get_Annotation_Label_setname(l)<<" : "<<endl;
			    for(m=0;m<number_of_each_annotation_label[l]; m++){
				if((annotation_labels[i][k][l][m])&&(annotation_label_scores[i][k][l][m]!=0)){
				    o<<MP->get_Annotation_Label_name(l,m)<<" : "<<annotation_label_scores[i][k][l][m]<<".\n";
				}
			    }
			}
		    }
		    if (j == i) {
			o << " <= not in alphabet.";
		    }
		    o << "\n" << flush;
		}
	    }
	}
    }
    return(check);
}


int Sequence::replace_non_alphabet_contents_of_sequence(model_parameters* const MP)
{
    // note : replaces letters in sequence array which are not in alphabet by randomly generated letter in alphabet
  
    int check = 0;

    if (this->sequence == NULL) {
	cout << "ERROR: class Sequence::replace_non_alphabet_contents_of_sequence: sequence array is NULL.\n" << flush;
	check++;
    }

    if (check == 0) {
	
	int i = 0;    
	int size = MP->get_Alphabet_size();
	for (i=0; i<length_of_sequence; i++) {
	    if ((sequence[i]<0)||(sequence[i]>=size)) {
		sequence[i] = generate_random_letter(size);
	    }
	}
	
    }
    return(check);
}			

int Sequence::annotation_labels_exist() const
{
    if((!annotation_labels)||(!annotation_label_scores)){
	return (0);
    }
    return (1);
}

int Sequence::set_sp_xsteps(const int steps, const int* const xsteps_array) {

    int check = 0;

    if ((steps != sp_steps) || (steps <= 0)) {
	cout << "ERROR: Sequence::set_sp_xsteps : steps (" << steps << ") != sp_steps ("
	     << sp_steps << ") or steps (" << steps << ") <= 0.\n" << flush;
	check++;
    }
    if (xsteps_array == NULL) {
	cout << "ERROR: Sequence::set_sp_xsteps : array xsteps_array is NULL.\n" << flush;
	check++;
    }
    if (check == 0) {
	
	if (sp_xsteps) delete [] sp_xsteps;
	sp_xsteps = NULL;
	
	sp_xsteps = new int[sp_steps];
    
	int i = 0;

	for (i = 0; i < sp_steps; i++) {
	    sp_xsteps[i] = xsteps_array[i];

	}
    }
    return(check);
}

int Sequence::set_sp_ysteps(const int steps, const int* const ysteps_array) {

    int check = 0;

    if ((steps != sp_steps) || (steps <= 0)) {
	cout << "ERROR: Sequence::set_sp_ysteps : steps (" << steps << ") != sp_steps ("
	     << sp_steps << ") or steps (" << steps << ") <= 0.\n" << flush;
	check++;
    }
    if (ysteps_array == NULL) {
	cout << "ERROR: Sequence::set_sp_ysteps : array ysteps_array is NULL.\n" << flush;
	check++;
    }
    if (check == 0) {
	
	if (sp_ysteps) delete [] sp_ysteps;
	sp_ysteps = NULL;
	
	sp_ysteps = new int[sp_steps];
	
	int i = 0;
	
	for (i = 0; i < sp_steps; i++) {
	    sp_ysteps[i] = ysteps_array[i];
	    
	}
    }
    return(check);
}

int Sequence::set_sp_scores(const int steps, const Score* const scores_array) {

    int check = 0;

    if ((steps != sp_steps) || (steps <= 0)) {
	cout << "ERROR: Sequence::set_sp_scores : steps (" << steps << ") != sp_steps ("
	     << sp_steps << ") or steps (" << steps << ") <= 0.\n" << flush;
	check++;
    }
    if (scores_array == NULL) {
	cout << "ERROR: Sequence::set_sp_scores : array scores_array is NULL.\n" << flush;
	check++;
    }
    if (check == 0) {
	
	if (sp_scores) delete [] sp_scores;
	sp_scores = NULL;
	
	sp_scores = new Score[sp_steps];
	
	int i = 0;
	
	for (i = 0; i < sp_steps; i++) {
	    sp_scores[i] = scores_array[i];
	    
	}
    }
    return(check);
}

int Sequence::set_sp_states(const int steps, const int* const states_array) {

    int check = 0;

    if ((steps != sp_steps) || (steps <= 0)) {
	cout << "ERROR: Sequence::set_sp_states : steps (" << steps << ") != sp_steps ("
	     << sp_steps << ") or steps (" << steps << ") <= 0.\n" << flush;
	check++;
    }
    if (states_array == NULL) {
	cout << "ERROR: Sequence::set_sp_states : array states_array is NULL.\n" << flush;
	check++;
    }
    if (check == 0) {

	if (sp_states) delete [] sp_states;
	sp_states = NULL;
	
	sp_states = new int[sp_steps];
    
	int i = 0;
	
	for (i = 0; i < sp_steps; i++) {
	    sp_states[i] = states_array[i];
	}
    }
    return(check);
}

int Sequence::set_sp_steps(const int steps) {

    int check = 0;

    if (steps <= 0) {
	cout << "ERROR: Sequence::set_sp_steps : steps (" << steps << ") <= 0.\n" << flush;
	check++;
    }
    if ((sp_xsteps != NULL) ||
	(sp_ysteps != NULL) ||
	(sp_states != NULL) ||
	(sp_scores != NULL)) {
	
	cout << "ERROR: Sequence::set_sp_steps : arrays sp_xsteps, sp_ysteps, sp_states and/or sp_scores have already been filled."
	     << " Cannot set sp_steps.\n" << flush;
	check++;
    }
    if (check == 0) {
	
	sp_steps = steps;
	
    }
    return(check);
}

Score Sequence::get_sp_score(const int xsteps, const int ysteps, const int state) const {
  
    Score return_score = Logzero;
  
    int check = 0;

    if ((sp_xsteps == NULL) ||
	(sp_ysteps == NULL) ||
	(sp_states == NULL) ||
	(sp_scores == NULL)) {
	
	cout << "ERROR: Sequence::get_sp_score : arrays sp_xsteps, sp_ysteps, sp_states and/or sp_scores have not yet been filled.";
	check++;
    }
    if (check == 0) {
	
	int i = 0;
	
	for (i=0; i<sp_steps; i++) {
	    
	    if ((sp_xsteps[i] == xsteps) &&
		(sp_ysteps[i] == ysteps) &&
		(sp_states[i] == state)) {

		return_score = sp_scores[i];
		break;
	    }
	}
    }
    
    return(return_score);
}

int Sequence::get_sp_steps(void) {
    return(sp_steps);
}

int Sequence::set_annotation(const char* file_name,
			     model_parameters* const MP)

    // note: use this function on special-file (gff-file option has not yet
    //       been implemented)
    //       old annotation is discarded if input values o.k.
    
{
    int check=0;
    
    // check input

    if (file_name == NULL)
    {
	cout << "ERROR: class Sequence::set_annotation: file_name is NULL.\n" << flush;
	check++;
    }
    if ((sequence == NULL) || (length_of_sequence == 0))
    {
	cout << "ERROR: class Sequence::set_annotation: sequence is NULL and length_of_sequence = 0.\n" << flush;
	check++;
    }
    if (this->get_ac() == NULL)
    {
	cout << "ERROR: class Sequence::set_annotation: ac is NULL.\n" << flush;
	check++;
    }
    if (this->get_start_position() > this->get_end_position())
    {
	cout << "ERROR: class Sequence::set_annotation: start_position (" << this->get_start_position()
	     << ") > end_position (" << this->get_end_position() << ").\n" << flush;
	check++;
    }
    
    if (check == 0)
    {
	// delete any existing information
	
	this->delete_annotation();
      
	// read annotation from special_file
	
	Info annotation;

	check += this->read_special_file(file_name,
					 MP,
					 &annotation);

	if (check != 0) {
	    cout << "ERROR: class Sequence::set_annotation: error occurred in function Sequence::read_special_file.\n" << flush;
	}
	
	// get arrays to set annotation
	
	bool****         new_annotation_labels        = NULL;
	double****       new_annotation_label_scores  = NULL;    

	int              new_number_of_sources = 0;
      
	if (check == 0) 
	{
	    new_number_of_sources = annotation.get_info_number_of_sources();

	    if(new_number_of_sources == 0)
	    {
		constraint = 0;
		return check;
	    }else{
		constraint = 1;
	    }
	  
	    check+=annotation.get_annotation_of_sequence(// input
		this,
		MP,
		// output
		&new_annotation_labels,
		&new_annotation_label_scores
		);
	  
	    if (check != 0) {
		cout << "ERROR: class Sequence::set_annotation: error occurred in function "
		     << "Info::get_annotation_of_sequence.\n" << flush;
	    }
	}
	// set annotation
	
	if (check == 0) {

	    int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();

	    check += this->set_number_of_sources(new_number_of_sources);
	    if (check != 0) {
		cout << "ERROR: class Sequence::set_annotation: error occurred in function "
		     << "Sequence::set_number_of_sources.\n" << flush;
	    }

	    check += this->set_number_of_annotation_labels(number_of_annotation_labels);
	    if (check != 0) {
		cout << "ERROR: class Sequence::set_annotation: error occurred in function "
		     << "Sequence::set_number_of_annotation_labels.\n" << flush;
	    }

	    int* tmp_number_of_each_annotation_label = new int[number_of_annotation_labels];
	    for(int i=0; i<number_of_annotation_labels; i++)
	    {
		tmp_number_of_each_annotation_label[i] = MP->get_Annotation_Label_size(i);
	    }

	    check += this->set_number_of_each_annotation_label(number_of_annotation_labels,
							       tmp_number_of_each_annotation_label);
	    if (check != 0) {
		cout << "ERROR: class Sequence::set_annotation: error occurred in function "
		     << "Sequence::set_number_of_each_annotation_label.\n" << flush;
	    }

	    if(tmp_number_of_each_annotation_label) delete[] tmp_number_of_each_annotation_label;
	    tmp_number_of_each_annotation_label = NULL;

	    check += this->set_annotation_labels(new_annotation_labels,
						 new_annotation_label_scores);
  	    
	    if (check != 0) {
		cout << "ERROR: class Sequence::set_annotation: error occurred in function "
		     << "Sequence::set_annotation_labels.\n" << flush;
	    }

	    // set score_labels
	    if(score_labels) delete [] score_labels;
	    score_labels = NULL;
	    
	    score_labels = new bool[number_of_annotation_labels];
	    for(int i =0; i<number_of_annotation_labels; i++)
	    {
		score_labels[i] = MP->get_score_of_Annotation_Label(i);
	    }	   
	}
	    
	// delete memory
	int i = 0;
	int j = 0;
	int k = 0;
	
	
	// int number_of_other_labels = MP->get_Total_Number_of_Other_Labels();

	if(new_annotation_labels)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if(new_annotation_labels[i]){
		    for(j=0; j<number_of_sources; j++)
		    {
			if(new_annotation_labels[i][j]){ 
			    for(k=0; k<number_of_annotation_labels; k++){
				if(new_annotation_labels[i][j][k]) delete[] new_annotation_labels[i][j][k];
				new_annotation_labels[i][j][k] = NULL;
			    }
			    delete [] new_annotation_labels[i][j];			   
			}
			new_annotation_labels[i][j] = NULL;
		    }
		    delete [] new_annotation_labels[i];
		}
		new_annotation_labels[i] = NULL;
	    }
	    delete[] new_annotation_labels;
	}
	new_annotation_labels = NULL;

	if(new_annotation_label_scores)
	{
	    for(i=0; i<length_of_sequence; i++)
	    {
		if(new_annotation_label_scores[i]){
		    for(j=0; j<number_of_sources; j++)
		    {
			if(new_annotation_label_scores[i][j]){ 
			    for(k=0; k<number_of_annotation_labels; k++){
				if(new_annotation_label_scores[i][j][k]) delete[] new_annotation_label_scores[i][j][k];
				new_annotation_label_scores[i][j][k] = NULL;
			    }
			    delete [] new_annotation_label_scores[i][j];			   
			}
			new_annotation_label_scores[i][j] = NULL;
		    }
		    delete [] new_annotation_label_scores[i];
		}
		new_annotation_label_scores[i] = NULL;
	    }
	    delete[] new_annotation_label_scores;
	}
	new_annotation_label_scores = NULL;
	 
	if (check != 0) {
	    this->delete_annotation();
	}
    } // if check == 0
    return(check);
}

int Sequence::change_annotation_labels(const int source,
				       const int label_set,
				       const int old_label,
				       const int new_label)
{

    int check = 0;
    
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::change_annotation_label: annotation_labels array is NULL.\n" << flush; 
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::change_annotation_label: annotation_label_scores array is NULL.\n" << flush; 
	check++;
    }
    
    if((source<0)||(source>=number_of_sources))
    {
	cout << "ERROR: class Sequence::change_annotation_label: source("
	     <<source<<") out of range [0,"<<number_of_sources-1
	     <<"].\n" << flush; 
	check++;
    }

    if ((label_set<0)||(label_set>=number_of_annotation_labels))
    {
	cout << "ERROR: class Sequence::change_annotation_label: label_set("
	     <<label_set<<") out of range [0,"<<number_of_annotation_labels-1
	     <<"].\n" << flush; 
	check++;
    }
    
    if(!number_of_each_annotation_label){
	cout << "ERROR: class Sequence::change_annotation_label: numebr_of_each_annotation_label array is NULL.\n" << flush; 
	check++;
    }else{
	if((old_label<0)||(old_label>=number_of_each_annotation_label[label_set]))
	{
	    cout << "ERROR: class Sequence::change_annotation_label: old_label("
		 <<old_label<<") out of range [0,"<<number_of_each_annotation_label[label_set]-1
		 <<"].\n" << flush; 
	    check++;
	}
	if((new_label<0)||(new_label>=number_of_each_annotation_label[label_set]))
	{
	    cout << "ERROR: class Sequence::change_annotation_label: old_label("
		 <<new_label<<") out of range [0,"<<number_of_each_annotation_label[label_set]-1
		 <<"].\n" << flush; 
	    check++;
	}
    }
    

    if (check == 0) 
    {	
	int i = 0;
	
	for (i=0; i<length_of_sequence; i++) 
	{
	    if (annotation_labels[i][source][label_set][old_label])		    
	    {
		if(!score_labels[label_set])
		{
		    annotation_labels[i][source][label_set][new_label] = true;
		}else if(annotation_label_scores[i][source][label_set][old_label]!=0){
		    annotation_labels[i][source][label_set][new_label] = true;
		    annotation_label_scores[i][source][label_set][new_label] = annotation_label_scores[i][source][label_set][old_label];
		    annotation_label_scores[i][source][label_set][old_label] = 0;
		}
	    }
	}
    }
    return(check);
}

int Sequence::change_annotation_labels(const int label_set,
				       const int old_label,
				       const int new_label)
{
    int check = 0;
    
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::change_annotation_label: annotation_labels array is NULL.\n" << flush; 
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::change_annotation_label: annotation_label_scores array is NULL.\n" << flush; 
	check++;
    }

    if ((label_set<0)||(label_set>=number_of_annotation_labels))
    {
	cout << "ERROR: class Sequence::change_annotation_label: label_set("
	     <<label_set<<") out of range [0,"<<number_of_annotation_labels-1
	     <<"].\n" << flush; 
	check++;
    }
    
    if(!number_of_each_annotation_label){
	cout << "ERROR: class Sequence::change_annotation_label: numebr_of_each_annotation_label array is NULL.\n" << flush; 
	check++;
    }else{
	if((old_label<0)||(old_label>=number_of_each_annotation_label[label_set]))
	{
	    cout << "ERROR: class Sequence::change_annotation_label: old_label("
		 <<old_label<<") out of range [0,"<<number_of_each_annotation_label[label_set]-1
		 <<"].\n" << flush; 
	    check++;
	}
	if((new_label<0)||(new_label>=number_of_each_annotation_label[label_set]))
	{
	    cout << "ERROR: class Sequence::change_annotation_label: old_label("
		 <<new_label<<") out of range [0,"<<number_of_each_annotation_label[label_set]-1
		 <<"].\n" << flush; 
	    check++;
	}
    }
    

    if (check == 0) 
    {	
	int i = 0;
	int j = 0;
	for(i=0; i<length_of_sequence; i++)
	{
	    for (j=0; j<number_of_sources; j++) 
	    {
		if (annotation_labels[i][j][label_set][old_label])		    
		{
		    if(!score_labels[label_set])
		    {
			annotation_labels[i][j][label_set][new_label] = true;
		    }else if(annotation_label_scores[i][j][label_set][old_label]!=0){
			annotation_labels[i][j][label_set][new_label] = true;
			annotation_label_scores[i][j][label_set][new_label] = annotation_label_scores[i][j][label_set][old_label];
			annotation_label_scores[i][j][label_set][old_label] = 0;
		    }
		}
	    }
	}
    }
    return(check);
}

int Sequence::change_all_annotation_labels(const int source,
					   const int label_set,
					   const int new_label)
{
    int check = 0;
    
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::change_all_annotation_labels: annotation_labels array is NULL.\n" << flush; 
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::change_all_annotation_labels: annotation_label_scores array is NULL.\n" << flush; 
	check++;
    }

    if((source<0)||(source>=number_of_sources))
    {
	cout << "ERROR: class Sequence::change_annotation_label: source("
	     <<source<<") out of range [0,"<<number_of_sources-1
	     <<"].\n" << flush; 
	check++;
    }
    
    if ((label_set<0)||(label_set>=number_of_annotation_labels))
    {
	cout << "ERROR: class Sequence::change_annotation_label: label_set("
	     <<label_set<<") out of range [0,"<<number_of_annotation_labels-1
	     <<"].\n" << flush; 
	check++;
    }
    if (check == 0) {
	
	int i = 0;
	int j = 0;
	for (i=0; i<length_of_sequence; i++) 
	{
	    if(score_labels[label_set])
	    {
		for(j=0; j<number_of_each_annotation_label[label_set]; j++){
		    if(j!=new_label){
			annotation_label_scores[i][source][label_set][j]  = 0;
		    }
		}
	    }else{
		for(j=0; j<number_of_each_annotation_label[label_set]; j++){
		    annotation_label_scores[i][source][label_set][j]  = 0;	
		    annotation_labels[i][source][label_set][j] = false;
		}
		annotation_labels[i][source][label_set][new_label] = true;		    
	    }	    
	}
    }
    return(check);
}

int Sequence::change_all_annotation_labels(const int label_set,
					   const int new_label)
{
    int check = 0;
    
    if (annotation_labels == NULL) 
    {
	cout << "ERROR: class Sequence::change_all_annotation_labels: annotation_labels array is NULL.\n" << flush; 
	check++;
    }
    if (annotation_label_scores == NULL) 
    {
	cout << "ERROR: class Sequence::change_all_annotation_labels: annotation_label_scores array is NULL.\n" << flush; 
	check++;
    }

    if ((label_set<0)||(label_set>=number_of_annotation_labels))
    {
	cout << "ERROR: class Sequence::change_annotation_label: label_set("
	     <<label_set<<") out of range [0,"<<number_of_annotation_labels-1
	     <<"].\n" << flush; 
	check++;
    }
    if (check == 0) {
	
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0; i<length_of_sequence; i++)
	{
	    for (j=0; j<number_of_sources; j++) 
	    {
		if(score_labels[label_set])
		{
		    for(k=0; k<number_of_each_annotation_label[label_set]; k++){
			if(k!=new_label){
			    annotation_label_scores[i][j][label_set][k]  = 0;
			}
		    }
		}else{
		    for(k=0; k<number_of_each_annotation_label[label_set]; k++){
			annotation_label_scores[i][j][label_set][k]  = 0;	
			annotation_labels[i][j][label_set][k] = false;
		    }
		    annotation_labels[i][j][label_set][new_label] = true;		    
		}
	    }
	}
    }
    return(check);
}

void Sequence::delete_annotation(void)
{
    int i = 0;
    int j = 0;
    int k = 0;

     if(annotation_labels){
	for(i=0; i<length_of_sequence; i++){
	    if(annotation_labels[i]){
		for(j=0; j<number_of_sources; j++){	
		    if(annotation_labels[i][j]){
			for(k=0; k<number_of_annotation_labels; k++){
			    if(annotation_labels[i][j][k]) delete[] annotation_labels[i][j][k];
			    annotation_labels[i][j][k] = NULL;
			}
			delete[] annotation_labels[i][j];
		    }
		    annotation_labels[i][j] = NULL;
		}
		delete [] annotation_labels[i];
	    }
	    annotation_labels[i] = NULL;
	}
	delete[] annotation_labels;
    }
    annotation_labels= NULL;    

    // remove annotation_label_scores
    if(annotation_label_scores){
	for(i=0; i<length_of_sequence; i++){
	    if(annotation_label_scores[i]){
		for(j=0; j<number_of_sources; j++){	
		    if(annotation_label_scores[i][j]){
			for(k=0; k<number_of_annotation_labels; k++){
			    if(annotation_label_scores[i][j][k]) delete[] annotation_label_scores[i][j][k];
			    annotation_label_scores[i][j][k] = NULL;
			}
			delete[] annotation_label_scores[i][j];
		    }
		    annotation_label_scores[i][j] = NULL;
		}
		delete [] annotation_label_scores[i];
	    }
	    annotation_label_scores[i] = NULL;
	}
	delete[] annotation_label_scores;
    }
    annotation_label_scores = NULL;
    
    return;
}

int Sequence::set_annotation_labels(const int       source,
				    bool***   const ann_labels,
				    double*** const ann_label_scores)
{
    // note: this function replaces any existing labels array by the ann_labels array

    int check = 0;

    if (ann_labels == NULL) 
    {
	cout << "ERROR : class Sequence::set_annotation_labels: ann_labels are NULL.\n" << flush;
	check++;
    }
    if (ann_label_scores == NULL) 
    {
	cout << "ERROR : class Sequence::set_annotation_labels: ann_label_scores are NULL.\n" << flush;
	check++;
    }
    if ((source<0)||(source<=number_of_sources))
    {
	cout<<"ERROR : class Sequence::set_annotation_labels: source("<<source
	    <<") out of range[0"<<number_of_sources-1<<"].\n"<<flush;
	check++;
    }

    if ((length_of_sequence == 0) || (sequence == NULL)) 
    {
	cout << "ERROR : class Sequence::set_annotation_labels: this sequence has not yet been initialise"
	     << "(length_of_sequence = 0 and sequence = NULL).\n" << flush;
	check++;
    }
  
    if (check == 0) 
    {
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0; i<length_of_sequence; i++)
	{
	    for (j=0; j<number_of_annotation_labels; j++) 
	    {
		for(k=0; k<number_of_each_annotation_label[j]; k++)
		{			
		    annotation_labels[source][i][j][k] = ann_labels[i][j][k];
		    annotation_label_scores[source][i][j][k] = ann_label_scores[i][j][k];
		}
	    }
	}
    }
    return(check);
}

int Sequence::set_annotation_labels(bool****   const ann_labels,
				    double**** const ann_label_scores)
{
    // note: this function replaces any existing labels array by the ann_labels array

    int check = 0;

    if (ann_labels == NULL) 
    {
	cout << "ERROR : class Sequence::set_annotation_labels: ann_labels are NULL.\n" << flush;
	check++;
    }
    if (ann_label_scores == NULL) 
    {
	cout << "ERROR : class Sequence::set_annotation_labels: ann_label_scores are NULL.\n" << flush;
	check++;
    }
    if ((length_of_sequence == 0) || (sequence == NULL)) 
    {
	cout << "ERROR : class Sequence::set_annotation_labels: this sequence has not yet been initialise"
	     << "(length_of_sequence = 0 and sequence = NULL).\n" << flush;
	check++;
    }
  
    if (check == 0) 
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	this->delete_annotation();

	// allocate memory for annotation_labels
	annotation_labels = new bool ***[length_of_sequence];
	for(i=0; i<length_of_sequence; i++){
	    annotation_labels[i] = new bool**[number_of_sources];
	    for(j=0; j<number_of_sources; j++){
		annotation_labels[i][j] = new bool*[number_of_annotation_labels];
		for(k=0; k<number_of_annotation_labels; k++){
		    annotation_labels[i][j][k] = new bool[number_of_each_annotation_label[k]];
		    for(l=0; l<number_of_each_annotation_label[k]; l++){
			annotation_labels[i][j][k][l] = false;
		    }
		}
	    }
	}

	// allocate memory for annotation_label_scores
	annotation_label_scores = new double***[length_of_sequence];
	for(i=0; i<length_of_sequence; i++){
	    annotation_label_scores[i] = new double**[number_of_sources];
	    for(j=0; j<number_of_sources; j++){
		annotation_label_scores[i][j] = new double*[number_of_annotation_labels];
		for(k=0; k<number_of_annotation_labels; k++){
		    annotation_label_scores[i][j][k] = new double[number_of_each_annotation_label[k]];
		    for(l=0; l<number_of_each_annotation_label[k]; l++){
			annotation_label_scores[i][j][k][l] = 0;
		    }
		}
	    }
	}	
	
	for(i=0; i<length_of_sequence; i++)
	{
	    for(j=0; j<number_of_sources; j++)
	    {
		for (k=0; k<number_of_annotation_labels; k++) 
		{
		    for(l=0; l<number_of_each_annotation_label[k]; l++)
		    {			
			annotation_labels[i][j][k][l] = ann_labels[i][j][k][l];
			annotation_label_scores[i][j][k][l] = ann_label_scores[i][j][k][l];
		    }
		}
	    }
	}
    }
    return(check);
}

int Sequence:: set_score_labels(const int i, const bool s)
{
    int check = 0;
    if(!score_labels)
    {
	cout<<"Error: Sequence class:: set_score_labels "
	    <<"score_labels array is NULL.\n";
	check++;
    }
    if((i<0)||(i>=number_of_annotation_labels))
    {
	cout<<"Error: Sequence class:: set_score_labels "
	    <<"index("<<i<<") out of range[0.."
	    <<number_of_annotation_labels<<"].\n";
	check++;
    }
    if(!check)
    {
	score_labels[i] = s;
    }
    return check;
}

int Sequence:: set_constraint(const int c)
{
    int check = 0;
    if((c!=0)&&(c!=1))
    {
	check++;
    }
    constraint = c;
    return check;
}

int Sequence::set_ac(const char* const n)
{
    if (ac) delete [] ac;
    ac=NULL;
    int check = 0;
  
    if(!n){
	check++;
    }else{
	int length=strlen(n);
	ac=new char[length+1]; 
	for (int i=0; i<length+1; i++)
	{
	    ac[i]=n[i];
	}
    }
    return check;
}
 
int Sequence::set_sequence_start_and_end(const int sequence_start,
					 const int sequence_end)
{
    int check = 0;

    if (sequence_end < sequence_start)
    {
	cout << "ERROR class Sequence::set_sequence_start_and_end : sequence_start (" << sequence_start
	     << ") must be <= sequence_end (" << sequence_end << ").\n" << flush;
	check++;
    }
    if ((sequence_end - sequence_start + 1) != length_of_sequence)
    {
	cout << "ERROR class Sequence::set_sequence_start_and_end : sequence_end (" << sequence_end
	   << ") - sequence_start (" << sequence_start << ") + 1 != length_of_sequence ("
	     << length_of_sequence << ").\n" << flush;
	check++;
    }
    if (check == 0)
    {
	start_position = sequence_start;
	end_position   = sequence_end;
    }
    return(check);
}

int Sequence:: set_number_of_special_emissions(const int n_of_special_emissions)
{
    int check=0;
    if(n_of_special_emissions<0)
    {
	check++;
    }
    number_of_special_emissions = n_of_special_emissions;
    return check;
}

int Sequence:: set_number_of_sources(const int n_of_source)
{
    int check = 0;
    if(n_of_source<0){
	cout<<"ERROR: class Sequence: set_number_of_sources, n_of_source("
	    <<n_of_source<<") out of range.\n"<<flush;
	check++;
    }
    if(!check){
	number_of_sources = n_of_source;
    }
    return check;
}

int Sequence:: set_number_of_annotation_labels(const int n_of_annotation_labels)
{
    int check = 0;
    if(n_of_annotation_labels<0){
	cout<<"ERROR: class Sequence: set_number_of_annotation_labels, n_of_annotation_labels("
	    <<n_of_annotation_labels<<") out of range.\n"<<flush;
	check++;
    }
    if(!check){
	number_of_annotation_labels = n_of_annotation_labels;
    }
    return check;
}

int Sequence::set_number_of_each_annotation_label(const int n_of_annotation_labels,
						   int* const n_of_each_annotation_label)
{
    int check = 0;
    int i = 0;
    if(n_of_annotation_labels!=number_of_annotation_labels){
	cout<<"ERROR: class Sequence: set_number_of_annotation_labels," 
	    <<"n_of_annotation_labels("<<n_of_annotation_labels<<") != "
	    <<"number_of_annotation_labels("<<number_of_annotation_labels<<").\n"
	    <<flush;
	check++;
    }
    
    if(!check){
	if(!number_of_each_annotation_label)
	{
	    number_of_each_annotation_label = new int[number_of_annotation_labels];
	}
	for(i=0; i<number_of_annotation_labels; i++){
	    number_of_each_annotation_label[i] = n_of_each_annotation_label[i];
	}
    }
    return check;    

}

int Sequence::set_posterior_probs(const int   number_of_state,
				  const Prob* probs)
{
    int check = 0;

    if (probs == NULL) {
	cout << "ERROR: Sequence::set_posterior_probs: array probs is NULL.\n" << flush;
	check++;
    }
    else {
	int i = 0;
	for (i=0; i<length_of_sequence; i++) {
	    if ((probs[i] < 0.0) || (probs[i] > 1.0)) {
		cout << "ERROR: Sequence::set_posterior_probs: probs[" << i << "] ("
		     << probs[i] << ") out of range [0,1].\n" << flush;
		check++;
		break;
	    }
	}
    }
    
    // first check if sequence object was declared to have special transitions
    
    if (number_of_special_emissions == 0) {
	cout << "ERROR: Sequence::set_posterior_probs: sequence object was declared to have zero special \n"
	     << "emissions. Cannot set any posterior probs.\n" << flush;
	check++;
    }
    else { // if sequence object was declared to have > 0 special emissions
	if ((posterior_probs_for_special_emissions.GetNumberofDimensions() == 0) &&
	    (scores_for_special_emissions.GetNumberofDimensions() == 2)) {
	    cout << "ERROR: Sequence::set_posterior_probs: sequence object is ment to be used with scores \n"
		 << "not posterior probs. Cannot set any posterior probs.\n" << flush;
	    check++;
	}
	
	// make sure that number_of_state is in the correct range
	
	if ((number_of_state < 0) || 
	    (number_of_state > (indices_of_special_emissions.GetDimension(0)-1))) {
	    cout << "ERROR: Sequence::set_posterior_probs: number_of_state ("
		 << number_of_state << ") out of range (0.."
		 << indices_of_special_emissions.GetDimension(0) << ").\n" << flush;
	    check++;
	}
	// make sure that posterior probs of this special emission were not yet implemented
	
	if (check == 0) {
	    if (this->get_implementation_status_of_special_emission(number_of_state) == 1) {
		cout << "ERROR: Sequence::set_posterior_probs: posterior probs for special emission \n"
		     << "number_of_state (" << number_of_state << ") were already implemented.\n" << flush;
		check++;
	    }
	}    
	// make sure that there is space left to implement the new information for this special state
	
	if (number_of_implemented_special_emissions >= number_of_special_emissions) {
	    cout << "ERROR: Sequence::set_posterior_probs: cannot implement information for special emission \n"
		 << "number_of_state (" << number_of_state 
		 << ") as the number_of_implemented_special_emissions ("
		 << number_of_implemented_special_emissions << ") > number_of_special_emissions ("
		 << number_of_special_emissions << ").\n" << flush;
	    check++;
	}
	// check that Sequence object is in the right state
	
	if (start_position >= end_position) {
	    cout << "ERROR: Sequence::set_posterior_probs: start_position (" << start_position
		 << ") >= end_position (" << end_position 
		 << ") for this Sequence instance. \n"
		 << "Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
	
	if ((end_position - start_position +1) != length_of_sequence) {
	    cout << "ERROR: Sequence::set_posterior_probs: end_position (" << end_position
		 << ") - start_position (" << start_position << ") + 1 != length_of_sequence ("
		 << length_of_sequence << "). Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
    } // if sequence object was declared to have special emissions
    
    if (check == 0) {
	
	int values_set = 0;
	
	int i = 0;
	array<int> Index(1);
	const int index_of_special_emission = number_of_implemented_special_emissions;    
	
	// store info about the implementation in sequence
	
	number_of_implemented_special_emissions++;
	this->set_implementation_status_of_special_emission(number_of_state, 1);    
	this->set_index_of_special_emission(number_of_state,
					    index_of_special_emission);
	
#ifdef _DEBUG
	cout << "implementation_status_of_special_emission(" << number_of_state << ") = "
	     << implementation_status_of_special_emissions.GetElement(number_of_state) << "\n" << flush;
#endif

	Index.SetDimension(0, 2);
	Index.SetElement(0, index_of_special_emission); 
	
	for (i=0; i<length_of_sequence; i++) {
	    Index.SetElement(1, i);
	    posterior_probs_for_special_emissions.SetElement(Index, probs[i]);
	    values_set++;
	    
#ifdef _DEBUG
	    cout << "probs[" << Index.GetElement(0) << "][" << Index.GetElement(1) << "] = "
	   << posterior_probs_for_special_emissions.GetElement(Index) << "\n" << flush;
#endif
	}
	if (values_set <= 0) {
	    cout << "NOTE : Sequence::set_posterior_probs: state (" << number_of_state << ") number of values set = "
		 << values_set << " is zero.\n" << flush;
	}
    }
    return(check);
}

int Sequence::set_posterior_probs(const int number_of_new_state, 
				  const int number_of_already_implemented_state)
{
    int check = 0;

    if (number_of_special_emissions == 0) {
	cout << "ERROR: Sequence::set_posterior_probs: sequence object was declared to have zero special \n"
	     << "emissions. Cannot set any posterior probs.\n" << flush;
	check++;
    }
    else { // if sequence object was declared to have > 0 special emissions
	if ((posterior_probs_for_special_emissions.GetNumberofDimensions() == 0) &&
	    (scores_for_special_emissions.GetNumberofDimensions() == 2)) {
	    cout << "ERROR: Sequence::set_posterior_probs: sequence object is ment to be used with scores \n"
		 << "not posterior probs. Cannot set any posterior probs.\n" << flush;
	    check++;
	}
	
	// make sure that number_of_from_state and number_of_to_state are in the correct range
	
	if ((number_of_new_state < 0) || 
	    (number_of_new_state > (indices_of_special_emissions.GetDimension(0)-1))) {
	    cout << "ERROR: Sequence::set_posterior_probs: number_of_new_state ("
		 << number_of_new_state << ") out of range (0.."
		 << indices_of_special_emissions.GetDimension(0) << ").\n" << flush;
	    check++;
	}
	if ((number_of_already_implemented_state < 0) || 
	    (number_of_already_implemented_state > (indices_of_special_emissions.GetDimension(0)-1))) {
	    cout << "ERROR: Sequence::set_posterior_probs: number_of_already_implemented_state ("
		 << number_of_already_implemented_state << ") out of range (0.."
		 << indices_of_special_emissions.GetDimension(0) << ").\n" << flush;
	    check++;
	}
	// make sure that posterior probs of new special emission were not yet implemented
	
	if (check == 0) {
	    if (this->get_implementation_status_of_special_emission(number_of_new_state) == 1) {
		cout << "ERROR: Sequence::set_posterior_probs: posterior probs for new special emission \n"
		     << "number_of_new_state (" << number_of_new_state 
		     << ") were already implemented.\n" << flush;
		check++;
	    }
	}
	// make sure that posterior probs of already implemented special emission were really already implemented
	
	if (check == 0) {
	    if (this->get_implementation_status_of_special_emission(number_of_already_implemented_state) == 0) {
		cout << "ERROR: Sequence::set_posterior_probs: posterior probs for supposedly \n"
		     << "already implemented special emission number_of_already_implemented_state (" 
		     << number_of_already_implemented_state << ") were not yet implemented.\n" << flush;
		check++;
	    }
	}
	// make sure that there is space left to implement the new information for this special state
	
	if (number_of_implemented_special_emissions >= number_of_special_emissions) {
	    cout << "ERROR: Sequence::set_posterior_probs: cannot implement information for special emission \n"
		 << "number_of_new_state (" << number_of_new_state 
		 << ") as the number_of_implemented_special_emissions ("
		 << number_of_implemented_special_emissions << ") > number_of_special_emissions ("
		 << number_of_special_emissions << ").\n" << flush;
	    check++;
	}
	// check that Sequence object is in the right state
	
	if (start_position >= end_position) {
	    cout << "ERROR: Sequence::set_posterior_probs: start_position (" << start_position
		 << ") >= end_position (" << end_position 
		 << ") for this Sequence instance. \n"
		 << "Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
	if ((end_position - start_position +1) != length_of_sequence) {
	    cout << "ERROR: Sequence::set_posterior_probs: end_position (" << end_position
		 << ") - start_position (" << start_position << ") + 1 != length_of_sequence ("
		 << length_of_sequence << "). Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
    } // if sequence object was declared to have special emissions
    
    if (check == 0) {

	int values_set = 0;
	
	int index_of_already_implemented_special_emission = 
	    this->get_index_of_special_emission(number_of_already_implemented_state);
	
	// store info about the implementation in sequence
	
	int index_of_new_special_emission = number_of_implemented_special_emissions;
	number_of_implemented_special_emissions++;
	this->set_implementation_status_of_special_emission(number_of_new_state, 1);
	
	this->set_index_of_special_emission(number_of_new_state,
					    index_of_new_special_emission);
	
	// copy entries for special emission number_of_already_implemented_state
	// into number_of_new_state
	
	array<int> index(1);
	index.SetDimension(0, 2);
	index.SetElement(0, index_of_already_implemented_special_emission);
	
	array<int> new_index(1);
	new_index.SetDimension(0, 2);
	new_index.SetElement(0, index_of_new_special_emission); 
	
	Prob prob = static_cast<Prob>(0.0);
	
	int i = 0;
	for (i=0; i<length_of_sequence; i++) {
	    index.SetElement(1, i);
	    new_index.SetElement(1, i);
	    
	    prob = static_cast<Prob>(0.0);
	    prob = posterior_probs_for_special_emissions.GetElement(index);
	    posterior_probs_for_special_emissions.SetElement(new_index, prob);
	    values_set++;

	}      

	if (values_set <= 0) {
	    cout << "NOTE : Sequence::set_posterior_probs: state (" << number_of_new_state 
		 << ") number of values set = "
		 << values_set << " is zero.\n" << flush;
	}
    }
    return(check);
}

int Sequence::set_scores(const int    number_of_state,
			 const Score* scores)
{
    int check = 0;

    if (scores == NULL) 
    {
	cout << "ERROR: Sequence::set_scores: array scores is NULL.\n" << flush;
	check++;
    }

    // first check if sequence object was declared to have special transitions
    
    if (number_of_special_emissions == 0) 
    {
	cout << "ERROR: Sequence::set_scores: sequence object was declared to have zero special \n"
	     << "emissions. Cannot set any scores.\n" << flush;
	check++;
    }
    else 
    {
	
        // if sequence object was declared to have > 0 special emissions
	
	if ((posterior_probs_for_special_emissions.GetNumberofDimensions() == 2) &&
	    (scores_for_special_emissions.GetNumberofDimensions() == 0)) {
	    cout << "ERROR: Sequence::set_scores: sequence object is ment to be used with posterior probs \n"
		 << "not scores. Cannot set any scores.\n" << flush;
	    check++;
	}

	// make sure that number_of_state is in the correct range
	
	if ((number_of_state < 0) || 
	    (number_of_state > (indices_of_special_emissions.GetDimension(0)-1))) {
	    cout << "ERROR: Sequence::set_scores: number_of_state ("
		 << number_of_state << ") out of range (0.."
		 << indices_of_special_emissions.GetDimension(0) << ").\n" << flush;
	    check++;
	}

	// make sure that posterior probs of this special emission were not yet implemented	
	
	if (check == 0) 
	{
	    if (this->get_implementation_status_of_special_emission(number_of_state) == 1) 
	    {
		cout << "ERROR: Sequence::set_scores: scores for special emission \n"
		     << "number_of_state (" << number_of_state << ") were already implemented.\n" << flush;
		check++;
	    }
	}    
	
	// make sure that there is space left to implement the new information for this special state
	
	if (number_of_implemented_special_emissions >= number_of_special_emissions) 
	{
	    cout << "ERROR: Sequence::set_scores: cannot implement information for special emission \n"
		 << "number_of_state (" << number_of_state 
		 << ") as the number_of_implemented_special_emissions ("
		 << number_of_implemented_special_emissions << ") > number_of_special_emissions ("
		 << number_of_special_emissions << ").\n" << flush;
	    check++;
	}

	// check that Sequence object is in the right state
	
	if (start_position >= end_position) 
	{
	    cout << "ERROR: Sequence::set_scores: start_position (" << start_position
		 << ") >= end_position (" << end_position 
		 << ") for this Sequence instance. \n"
		 << "Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}

	if ((end_position - start_position +1) != length_of_sequence) 
	{
	    cout << "ERROR: Sequence::set_scores: end_position (" << end_position
		 << ") - start_position (" << start_position << ") + 1 != length_of_sequence ("
		 << length_of_sequence << "). Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
	
    } // if sequence object was declared to have special emission
    
    if (check == 0) {
	int values_set = 0;
	
	int i = 0;
	array<int> Index(1);
	const int index_of_special_emission = number_of_implemented_special_emissions;    
	number_of_implemented_special_emissions++;
	this->set_implementation_status_of_special_emission(number_of_state, 1);    
	this->set_index_of_special_emission(number_of_state,
					    index_of_special_emission);
	
#ifdef _DEBUG
	cout << "implementation_status_of_special_emission(" << number_of_state << ") = "
	     << implementation_status_of_special_emissions.GetElement(number_of_state) << "\n" << flush;
#endif

	Index.SetDimension(0, 2);
	Index.SetElement(0, index_of_special_emission); 
	
	for (i=0; i<length_of_sequence; i++) {
	    
	    Index.SetElement(1, i);
	    scores_for_special_emissions.SetElement(Index, scores[i]);
	    values_set++;
	    
#ifdef _DEBUG
	    cout << "scores[" << Index.GetElement(0) << "][" << Index.GetElement(1) << "] = "
		 << scores_for_special_emissions.GetElement(Index) << "\n" << flush;
#endif
	}
	
	if (values_set <= 0) {
	    cout << "NOTE : Sequence::set_scores: state (" << number_of_state << ") number of values set = "
		 << values_set << " is zero.\n" << flush;
	}
    }
    return(check);
}
 

int Sequence::set_scores(const int number_of_new_state,
			 const int number_of_already_implemented_state)
{
    int check = 0;
    
    if (number_of_special_emissions == 0) 
    {
	cout << "ERROR: Sequence::set_scores: sequence object was declared to have zero special \n"
	     << "emissions. Cannot set any scores.\n" << flush;
	check++;
    }
    else 
    { // if sequence object was declared to have > 0 special emissions
	if ((posterior_probs_for_special_emissions.GetNumberofDimensions() == 2) &&
	    (scores_for_special_emissions.GetNumberofDimensions() == 0)) {
	    cout << "ERROR: Sequence::set_scores: sequence object is ment to be used with posterior probs \n"
		 << "not scores. Cannot set any scores.\n" << flush;
	    check++;
	}
	
	// make sure that number_of_from_state and number_of_to_state are in the correct range
	
	if ((number_of_new_state < 0) || 
	    (number_of_new_state > (indices_of_special_emissions.GetDimension(0)-1))) {
	    cout << "ERROR: Sequence::set_scores: number_of_new_state ("
		 << number_of_new_state << ") out of range (0.."
		 << indices_of_special_emissions.GetDimension(0) << ").\n" << flush;
	    check++;
	}
	if ((number_of_already_implemented_state < 0) || 
	    (number_of_already_implemented_state > (indices_of_special_emissions.GetDimension(0)-1))) {
	    cout << "ERROR: Sequence::set_scores: number_of_already_implemented_state ("
		 << number_of_already_implemented_state << ") out of range (0.."
		 << indices_of_special_emissions.GetDimension(0) << ").\n" << flush;
	    check++;
	}
	// make sure that scores of new special emission were not yet implemented
	
	if (check == 0) {
	    if (this->get_implementation_status_of_special_emission(number_of_new_state) == 1) {
		cout << "ERROR: Sequence::set_scores: scores for new special emission \n"
		     << "number_of_new_state (" << number_of_new_state 
		     << ") were already implemented.\n" << flush;
		check++;
	    }
	}
	// make sure that scores of already implemented special emission were really already implemented
	
	if (check == 0) {
	    if (this->get_implementation_status_of_special_emission(number_of_already_implemented_state) == 0) {
		cout << "ERROR: Sequence::set_scores: scores for supposedly \n"
		     << "already implemented special emission number_of_already_implemented_state (" 
		     << number_of_already_implemented_state << ") were not yet implemented.\n" << flush;
		check++;
	    }
	}
	// make sure that there is space left to implement the new information for this special state
	
	if (number_of_implemented_special_emissions >= number_of_special_emissions) {
	    cout << "ERROR: Sequence::set_scores: cannot implement information for special emission \n"
		 << "number_of_new_state (" << number_of_new_state 
		 << ") as the number_of_implemented_special_emissions ("
		 << number_of_implemented_special_emissions << ") > number_of_special_emissions ("
		 << number_of_special_emissions << ").\n" << flush;
	    check++;
	}
	// check that Sequence object is in the right state

	if (start_position >= end_position) {
	    cout << "ERROR: Sequence::set_scores: start_position (" << start_position
		 << ") >= end_position (" << end_position 
		 << ") for this Sequence instance. \n"
		 << "Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
	if ((end_position - start_position +1) != length_of_sequence) {
	    cout << "ERROR: Sequence::set_scores: end_position (" << end_position
		 << ") - start_position (" << start_position << ") + 1 != length_of_sequence ("
		 << length_of_sequence << "). Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
    } // if sequence object was declared to have special emissions
    
    if (check == 0) {

	int values_set = 0;
	
	int index_of_already_implemented_special_emission = 
	    this->get_index_of_special_emission(number_of_already_implemented_state);
	
	// store info about the implementation in sequence
	
	int index_of_new_special_emission = number_of_implemented_special_emissions;
	number_of_implemented_special_emissions++;
	this->set_implementation_status_of_special_emission(number_of_new_state, 1);
	
#ifdef _DEBUG
	cout << "implementation_status_of_special_emission(" << number_of_new_state << ") = "
	     << implementation_status_of_special_emissions.GetElement(number_of_new_state) 
	     << " = " << this->get_implementation_status_of_special_emission(number_of_new_state) 
	 << "\n" << flush;
#endif    
	
	this->set_index_of_special_emission(number_of_new_state,
					    index_of_new_special_emission);
	
	// copy entries for special emission number_of_already_implemented_state
	// into number_of_new_state
	
	array<int> index(1);
	index.SetDimension(0, 2);
	index.SetElement(0, index_of_already_implemented_special_emission);
	
	array<int> new_index(1);
	new_index.SetDimension(0, 2);
	new_index.SetElement(0, index_of_new_special_emission); 
    
	Score score = static_cast<Score>(0.0);
	
	int i = 0;
	for (i=0; i<length_of_sequence; i++) {
	    index.SetElement(1, i);
	    new_index.SetElement(1, i);
	    
	    score = static_cast<Score>(0.0);
	    score = scores_for_special_emissions.GetElement(index);
	    scores_for_special_emissions.SetElement(new_index, score);
	    
	    values_set++;

	}      

	if (values_set <= 0) {
	    cout << "NOTE : Sequence::set_scores: state (" << number_of_new_state << ") number of values set = "
		 << values_set << " is zero.\n" << flush;
	}
    }
    return(check);
}

int Sequence::change_scores(const Score new_default_score)
{
    int check = 0;

    if (number_of_special_transitions == 0)
    {
	cout << "ERROR: Sequence::change_scores: sequence object was declared to have zero special \n"
	     << "transitions. Cannot set any scores.\n" << flush;
	check++;
    }
    else // if sequence object was declared to have > 0 special transitions
    {
	if ((posterior_probs_for_special_transitions.GetNumberofDimensions() == 2) &&
	    (scores_for_special_transitions.GetNumberofDimensions() == 0))
	{
	    cout << "ERROR: Sequence::change_scores: sequence object is ment to be used with posterior probs \n"
		 << "not scores. Cannot set any scores.\n" << flush;
	    check++;
	}
	
	// check that Sequence object is in the right state
	
	if (start_position >= end_position)
	{
	    cout << "ERROR: Sequence::change_scores: start_position (" << start_position
		 << ") >= end_position (" << end_position 
		 << ") for this Sequence instance. \n"
		 << "Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
	
	if ((end_position - start_position +1) != length_of_sequence)
	{
	    cout << "ERROR: Sequence::change_scores: end_position (" << end_position
		 << ") - start_position (" << start_position << ") + 1 != length_of_sequence ("
		 << length_of_sequence << "). Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
	if (ac == NULL)
	{
	    cout << "ERROR: Sequence::change_scores: ac of sequence is NULL. "
		 << "Sequence instance lacks information to implement new information.\n" << flush;
	    check++;
	}
    } // if sequence object was declared to have special transitions
    
    if (check == 0)
    {
	int i, j = 0;
	array<int> index(1);
	index.SetDimension(0, 2);
	Score score = 0;
	
	const int number_of_special_transitions = scores_for_special_transitions.GetDimension(0);
	const int length                        = scores_for_special_transitions.GetDimension(1);
	
	cout << "length                        = " << length << "\n" << flush;
	cout << "number_of_special_transitions = " << number_of_special_transitions << "\n" << flush;
	
	for (i=0; i<number_of_special_transitions; i++) {
	    
	    index.SetElement(0, i);
	    
	    for (j=0; j<length; j++) {
		
		index.SetElement(1, j);
		
		score = scores_for_special_transitions.GetElement(index);
		
		if (score != Logzero) {

		    cout << "scores_for_special_transitions[" << i << "][" << j << "] = "
			 << scores_for_special_transitions.GetElement(index);
		    
		    scores_for_special_transitions.SetElement(index, new_default_score);
		    
		    cout << " => " << scores_for_special_transitions.GetElement(index) << "\n" << flush;
		}
	    }
	}

    }
    return(check);
}

int Sequence::print_sequence_to_fasta_file(model_parameters* const MP, 
					   std::ostream &o) const
{
    // note : keep_orientation can be either 0 or 1, if 1 the sequence will be printed as indicated by
    //        its orientation (the plus strand, if orientation = 1 and the minus strand, if orientation = -1),
    //        if 0 the plus strand of the sequence will be printed whatever the orientation of the sequence
    
    int check = 0;

    if (sequence == NULL) {
	cout << "ERROR: class Sequence::print_sequence_to_fasta_file: sequence array is NULL.\n" << flush;
	check++;
    }
    if (check == 0) {

	// header 

	o << ">" << ac << "\t" << start_position << "-" << end_position;
	o << "\t" << "forward\n" << flush;

	// sequence
	
	const int max_fasta_line_length = 60;
	int i = 0;
	int j = 0;
	int alphabet = MP->get_Alphabet_size();
	
	int* new_sequence = NULL;
	new_sequence = new int[length_of_sequence];
    
	for (i=0; i<length_of_sequence; i++) {
	    new_sequence[i] = sequence[i];
	}

	if (check == 0) {
	    i = 0;
	    while (i<length_of_sequence) {
		for (j=i; j<min(max_fasta_line_length+i, length_of_sequence); j++) {
		    o << MP->get_Alphabet_name(new_sequence[j]);
		}
		o << "\n" << flush;
		i += min(max_fasta_line_length+i, length_of_sequence) - i;
	    }
	}
	o << flush;
	
	if (new_sequence) delete [] new_sequence;
	new_sequence = NULL;
    }
    return(check);
}

void Sequence::print_emission_information_linewise(std::ostream &o,
						   model_parameters* const MP) const
{
    if ((length_of_sequence > 0) && (sequence != NULL)) 
    {
	const int print_width = 3;
	int i, j, k, l, count = 0;
	int state_number = 0;
	array<int> index(1);
	index.SetDimension(0, 2);
	
	int abs_position = 0;
	
	abs_position = start_position;
	
	for (i=0; i<length_of_sequence; i++) {
	    
	    o << abs_position << "\t" << i << "\t";
	    if (sequence)     
	    {
		o << MP->get_Alphabet_name(sequence[i]) << "\t";
	    }

	    if(annotation_labels)
	    {
		for(j=0; j<number_of_sources; j++)
		{
		    o <<"source : "<<j<<"\t";
		    for(k=0; k<number_of_annotation_labels; k++)
		    {				
			for(l=0; l<number_of_each_annotation_label[k]; l++){
			    if(annotation_labels[i][j][k][l]){
				o<<MP->get_Annotation_Label_setname(k)<<":";
				o<<MP->get_Annotation_Label_name(k,l)<<"(";
				o<<annotation_label_scores[i][j][k][l]<<")\t";
			    }
			}
		    }
		}
	    }
	 
	    if (this->get_number_of_special_emissions() > 0) {
		for (l=0; l<this->get_number_of_implemented_special_emissions(); l++) {
		    
		    index.SetElement(0, l);
		    index.SetElement(1, i);
		    
		    for (j=0; j<implementation_status_of_special_emissions.GetDimension(0); j++) {
			if ((implementation_status_of_special_emissions.GetElement(j) == 1) &&
			    (indices_of_special_emissions.GetElement(j) == l)) {
			    
			    state_number = j;
			    break;
			}
		    }
		    if (scores_for_special_emissions.GetElement(index) != Logzero) {
			o << state_number << "\t";
		    }
		}
	    }
	    
	    o << "\n" << flush;
	    
	    abs_position += 1;
	}
    }
    return;
}

void Sequence::print_emission_and_transition_information_linewise(std::ostream &o,
								  model_parameters* const MP) const
{
    if ((length_of_sequence > 0) && (sequence != NULL)) {
	
	const int print_width = 3;
	int i, j, k, l, count = 0;
	int state_number = 0;
	array<int> index(1);
	index.SetDimension(0, 2);
	
	int abs_position = 0;
	
	abs_position = start_position;
	
	
	for (i=0; i<length_of_sequence; i++) {
	    
	    o << abs_position << "\t" << i << "\t";
	    if (sequence)     
	    {
		o << MP->get_Alphabet_name(sequence[i]) << "\t";
	    }
	    
	    if(annotation_labels)
	    {
		for(j=0; j<number_of_sources; j++)
		{
		    o <<"source : "<<j<<"\t";
		    for(k=0; k<number_of_annotation_labels; k++)
		    {				
			for(l=0; l<number_of_each_annotation_label[k]; l++){
			    if(annotation_labels[i][j][k][l]){
				o<<MP->get_Annotation_Label_setname(k)<<":";
				o<<MP->get_Annotation_Label_name(k,l)<<"(";
				o<<annotation_label_scores[i][j][k][l]<<")\t";
			    }
			}
		    }
		}
	    }
	    
	    if (this->get_number_of_special_emissions() > 0) {
		for (l=0; l<this->get_number_of_implemented_special_emissions(); l++) {
		    
		    index.SetElement(0, l);
		    index.SetElement(1, i);
		    
		    for (j=0; j<implementation_status_of_special_emissions.GetDimension(0); j++) {
			if ((implementation_status_of_special_emissions.GetElement(j) == 1) &&
			    (indices_of_special_emissions.GetElement(j) == l)) {
			    
			    state_number = j;
			    break;
			}
		    }
		    if (scores_for_special_emissions.GetElement(index) != Logzero) {
			o << state_number << "\t";
		    }
		}
	    }
	    
	    if ((this->get_number_of_special_transitions() > 0) &&
		(scores_for_special_transitions.GetNumberofDimensions() > 0)) {
		
		int index_i, j_1, j_2, k, count, linear_index = 0;
		int state_number_1 = 0;
		int state_number_2 = 0;
		Score score = Logzero;
		
		array<int> index(1);
		index.SetDimension(0, 2);
		
		for (index_i=0; index_i<this->get_number_of_implemented_special_transitions(); index_i++) {
		    
		    score = Logzero;
		    
		    for (j_1=0; j_1<implementation_status_of_special_transitions.GetDimension(0); j_1++) {
			
			index.SetElement(0, j_1);
			
			for (j_2=0; j_2<implementation_status_of_special_transitions.GetDimension(0); j_2++) {
			    
			    index.SetElement(1, j_2);
			    
			    if ((implementation_status_of_special_transitions.GetElement(index) == 1) &&
				(indices_of_special_transitions.GetElement(index) == index_i)) {
	    
				state_number_1 = j_1;
				state_number_2 = j_2;
				
				break;
			    }
			}
		    }
		    
		    linear_index = index_i * length_of_sequence + i;
		    score        = scores_for_special_transitions.GetElement(linear_index);
	  
		    if (score != Logzero) {
			o << state_number_1 << " -> " << state_number_2 << "\t";
		    }
		}
	    }
	    
	    o << "\n" << flush;
	    
	    abs_position += 1;
	}
    }
    return;
}

void Sequence::print_annotation_linewise(std::ostream &o,
					 model_parameters* const MP) const 
{
    if ((length_of_sequence > 0) && (sequence != NULL)) 
    {
	
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int abs_position = 0;

	abs_position = start_position;
      
	for(i=0; i<length_of_sequence; i++)
	{	    	    
	    if(annotation_labels)
	    {
		for(j=0; j<number_of_sources; j++)
		{	
		    for(k=0; k<number_of_annotation_labels; k++)
		    {				
			for(l=0; l<number_of_each_annotation_label[k]; l++){
			    if(annotation_labels[i][j][k][l]){
				o<<i<<"["
				 <<j<<"]"<<":"
				 <<"("<<MP->get_Annotation_Label_setname(k)<<","
				 <<MP->get_Annotation_Label_name(k,l)<<")"<<endl;
			    }
			}
		    }
		}
	    }
		
	    abs_position += 1;
	}
    }
    return;
}

void Sequence::print_annotation_linewise(std::ostream &o, 
					 const int start_rel_pos, 
					 const int end_rel_pos,
					 model_parameters* const MP) const
{
    int check = 0;

    if (end_rel_pos < start_rel_pos) 
    {
	cout << "ERROR: Sequence::print_annotation_linewise: start_rel_pos ("
	     << start_rel_pos << ") has to be <= end_rel_pos ("
	     << end_rel_pos << ").\n" << flush;
	check++;
    }
    if (start_rel_pos < 0) 
    {
	cout << "ERROR: Sequence::print_annotation_linewise: start_rel_pos ("
	     << start_rel_pos << ") has to be >= 0.\n" << flush;
	check++;
    }
    if (end_rel_pos > length_of_sequence-1) 
    {
	cout << "ERROR: Sequence::print_annotation_linewise: end_rel_pos ("
	     << end_rel_pos << ") has to be < length_of_sequence (" << length_of_sequence << ").\n" << flush;
	check++;
    }
    
    if (check == 0) 
    {
	if ((length_of_sequence > 0) && (sequence != NULL)) {
	    
	    int i = 0;
	    int j = 0;
	    int k = 0;
	    int l = 0;
	    int abs_position = 0;
	    
	    abs_position = start_position;
	    
	    for (i=start_rel_pos; i<end_rel_pos+1; i++) 
	    {

		o << abs_position << "\t" << i << "\t";
		if (sequence)     
		{
		    o << MP->get_Alphabet_name(sequence[i]) << "\t";
		}

		if(annotation_labels)
		{
		    for(j=0; j<number_of_sources; j++)
		    {
			o <<"source : "<<j<<"\t";
			for(k=0; k<number_of_annotation_labels; k++)
			{				
			    for(l=0; l<number_of_each_annotation_label[k]; l++){
				if(annotation_labels[i][j][k][l]){
				    o<<MP->get_Annotation_Label_setname(k)<<":";
				    o<<MP->get_Annotation_Label_name(k,l)<<"(";
				    o<<annotation_label_scores[i][j][k][l]<<")\t";
				}
			    }
			}
		    }
		}
		
		o << "\n" << flush;
		
		abs_position += 1;
	    }
	}
    }
    return;
}

void Sequence::print_char_sequence(model_parameters* const MP, std::ostream &o) const
{
    for (int i=0; i<length_of_sequence; i++)
    {
	o << MP->get_Alphabet_name(sequence[i]);
    }
    o << '\n';
    return;
}

void Sequence::print_sequence(std::ostream &o) const
{
    for (int i=0; i<length_of_sequence; i++) 
    {
	o << sequence[i];
    }
    o << '\n';
    return;
}

void Sequence::print_indices_of_implemented_special_transitions(std::ostream &o) const
{
    // check that sequence object has special transitions

    if (number_of_special_transitions > 0)
    {
	int i;
	int j;
	
	int max_n_of_child_state  = indices_of_special_transitions.GetDimension(0);
	int linear_index = 0;
	
	for (i=0; i<max_n_of_child_state; i++) // first index (from-state number)
	{
	    for (j=0; j<max_n_of_child_state; j++) // second index (to-state number)
	    {
		linear_index = max_n_of_child_state * i + j;
		
		if (implementation_status_of_special_transitions.GetElement(linear_index) == 1)
		{
		    o << "[" << i << ", " << j << "] " 
		      << indices_of_special_transitions.GetElement(linear_index) << "\n" << flush;
		}
	    }
	}
    }
    else
    {
	cout << "WARNING class Sequence::print_indices_of_implemented_special_transitions : \n"
	     << "sequence object has not been set up for any special transitions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_indices_of_implemented_special_emissions(std::ostream &o) const
{
    // check that sequence object has special emissions
    
    if (number_of_special_emissions > 0)
    {
	int i;
	int max_n_of_child_state  = indices_of_special_emissions.GetDimension(0);
	
	for (i=0; i<max_n_of_child_state; i++) {
	    if (implementation_status_of_special_emissions.GetElement(i) == 1) {
		o << "[" << i << "] " << indices_of_special_emissions.GetElement(i) << "\n" << flush;
	    }
	}
    }
    else
    {
	cout << "WARNING class Sequence::print_indices_of_implemented_special_emissions : \n"
	     << "sequence object has not been set up for any special emissions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_non_zero_posterior_probs(std::ostream &o) const
{
    if (number_of_special_transitions > 0)
    {
	posterior_probs_for_special_transitions.PrintonlyNonZerowithIndices(o);
    }
    else
    {
	cout << "WARNING class Sequence::print_non_zero_posterior_probs : \n"
	     << "sequence object has not been set up for any special transitions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_non_Logzero_scores(std::ostream &o) const
{
    if (number_of_special_transitions > 0)
    {
	scores_for_special_transitions.PrintonlyNonLogzerowithIndices(o);
    }
    else
    {
	cout << "WARNING class Sequence::print_non_Logzero_scores : \n"
	     << "sequence object has not been set up for any special transitions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_non_Logzero_emission_scores(std::ostream &o) const
{
    if (number_of_special_emissions > 0)
    {
	scores_for_special_emissions.PrintonlyNonLogzerowithIndices(o);
    }
    else
    {
	cout << "WARNING class Sequence::print_non_Logzero_emission_scores : \n"
	     << "sequence object has not been set up for any special emissions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_emission_scores(std::ostream &o) const
{
    if (number_of_special_emissions > 0)
    {
	const int length = scores_for_special_emissions.GetDimension(0) * 
	    scores_for_special_emissions.GetDimension(1);
	int i = 0;
	
	for (i=0; i<length; i++) {
	    cout << "[" << i << "]\t" << scores_for_special_emissions.GetElement(i) << "\n" << flush;
	}
    }
    else
    {
	cout << "WARNING class Sequence::print_emission_scores : \n"
	     << "sequence object has not been set up for any special emissions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_emission_scores_with_sequence_info(std::ostream &o,
							model_parameters* const MP) const
{

    if (this->get_number_of_special_emissions() > 0) 
    {

	const int print_width = 3;
	
	int i, j, k, count = 0;
	int state_number = 0;
	array<int> index(1);
	index.SetDimension(0, 2);

	for (i=0; i<this->get_number_of_implemented_special_emissions(); i++) {
	    
	    index.SetElement(0, i);
	    for (j=0; j<implementation_status_of_special_emissions.GetDimension(0); j++) {
		
		if ((implementation_status_of_special_emissions.GetElement(j) == 1) &&
		    (indices_of_special_emissions.GetElement(j) == i)) {
		    
		    state_number = j;

		    break;
		}
	    }
	    
	    count = 0;
	    
	    for (j=0; j<this->length(); j++) {
		
		index.SetElement(1, j);
		
		if (scores_for_special_emissions.GetElement(index) != Logzero) {
		    
		    count++;
		    
		    o << "[" << i << " (state " << state_number << "), " << j << "] = ";
		    
		    for (k = -print_width; k<print_width+1; k++) {
			
			if (((k+j) >= 0) && ((k+j) < this->length())) {
			    if (k == 0) {o << "(";}
			    o << MP->get_Alphabet_name(this->letter(k+j));
			    if (k == 0) {o << ")";}
			}
		    }
		    
		    o << "\t" << scores_for_special_emissions.GetElement(index) << "\n" << flush;
		}
	    }
	    
	    if (count == 0) {
		cout << "NOTE: special emission for state " << state_number 
		     << " has been implemented, but this state may never be used.\n" << flush;
	    }
	}
    }
    else {
	cout << "WARNING class Sequence::print_emission_scores_with_sequence_info : \n"
	     << "sequence object has not been set up for any special emissions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_implemented_special_transitions(std::ostream &o) const
{
    if (number_of_special_transitions > 0)
    {
	implementation_status_of_special_transitions.PrintonlyNonZerowithIndices(o);
    }
    else
    {
	cout << "WARNING class Sequence::print_implemented_special_transitions : \n"
	     << "sequence object has not been set up for any special transitions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_implemented_special_emissions(std::ostream &o) const
{
    if (number_of_special_emissions > 0)
    {
	implementation_status_of_special_emissions.PrintonlyNonZerowithIndices(o);
    }
    else
    {
	cout << "WARNING class Sequence::print_implemented_special_emissions : \n"
	     << "sequence object has not been set up for any special emissions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_non_zero_posterior_probs_with_sequence_info(std::ostream &o,
								 model_parameters* const MP) const
{
    if ((this->get_number_of_special_transitions() > 0) &&
	(posterior_probs_for_special_transitions.GetNumberofDimensions() > 0)) {
	
	int i, j, j_1, j_2, k, count = 0;
	int state_number_1 = 0;
	int state_number_2 = 0;
	
	array<int> index(1);
	index.SetDimension(0, 2);
	
	cout << "number_of_implemented_special_transitions = "
	     << this->get_number_of_implemented_special_transitions() << "\n" << flush;
      
	for (i=0; i<this->get_number_of_implemented_special_transitions(); i++) 
	{

	    for (j_1=0; j_1<implementation_status_of_special_transitions.GetDimension(0); j_1++) 
	    {

		index.SetElement(0, j_1);
		
		for (j_2=0; j_2<implementation_status_of_special_transitions.GetDimension(0); j_2++) 
		{
		    
		    index.SetElement(1, j_2);

		    if ((implementation_status_of_special_transitions.GetElement(index) == 1) &&
			(indices_of_special_transitions.GetElement(index) == i)) {
			
			state_number_1 = j_1;
			state_number_2 = j_2;

			break;
		    }
		}
	    }
	    
	    this->print_non_zero_posterior_probs_with_sequence_info(state_number_1, state_number_2, o,MP);
	}
    }
    else {
	o << "Sequence::print_non_zero_posterior_probs_with_sequence_info: this sequence has not been set up to have "
	  << "posterior probs and/or special transitions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_non_Logzero_scores_with_sequence_info(model_parameters* const MP,
							   std::ostream &o) const
{
    if ((this->get_number_of_special_transitions() > 0) &&
	(scores_for_special_transitions.GetNumberofDimensions() > 0)) {
	
	int i, j, j_1, j_2, k, count = 0;
	int state_number_1 = 0;
	int state_number_2 = 0;
	
	array<int> index(1);
	index.SetDimension(0, 2);

	for (i=0; i<this->get_number_of_implemented_special_transitions(); i++) {

	    for (j_1=0; j_1<implementation_status_of_special_transitions.GetDimension(0); j_1++) {
		
		index.SetElement(0, j_1);
		
		for (j_2=0; j_2<implementation_status_of_special_transitions.GetDimension(0); j_2++) {
		    
		    index.SetElement(1, j_2);

		    if ((implementation_status_of_special_transitions.GetElement(index) == 1) &&
			(indices_of_special_transitions.GetElement(index) == i)) {
			
			state_number_1 = j_1;
			state_number_2 = j_2;
			break;
		    }
		}
	    }
	    
	    this->print_non_Logzero_scores_with_sequence_info(state_number_1, state_number_2, o, MP);
	}
    }
    else {
	o << "Sequence::print_non_Logzero_scores_with_sequence_info: this sequence has not been set up to have "
	  << "scores and/or special transitions. There is nothing to print.\n" << flush;
    }
    return;
}

void Sequence::print_non_zero_posterior_probs_with_sequence_info(const int number_of_from_state,
								 const int number_of_to_state,
								 std::ostream &o,
								 model_parameters* const MP) const
{
    int check = 0;

    if (number_of_from_state < 0)
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0)
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    if (check == 0) 
    {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) == 0)
	{
	    cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : \n"
		 << "special transition number_of_from_state ("
		 << number_of_from_state << ") -> number_of_to_state (" << number_of_to_state
		 << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }

    if ((check == 0) && (posterior_probs_for_special_transitions.GetNumberofDimensions() > 0)) 
    {
	
	int  count = 0;
	int  i;
	int  position_in_sequence = 0;
	Prob posterior_prob       = static_cast<Prob>(0.0);
	
	int linear_index = 0;
	int index_of_special_transition = this->get_index_of_special_transition(number_of_from_state,
										number_of_to_state);
	
	for (position_in_sequence=0; position_in_sequence<length_of_sequence; position_in_sequence++)
	{
	    linear_index = index_of_special_transition * length_of_sequence + position_in_sequence;
	    posterior_prob = posterior_probs_for_special_transitions.GetElement(linear_index);
	    
	    if (posterior_prob != static_cast<Prob>(0.0))
	    {
		count++;
		
		o << number_of_from_state << " -> " << number_of_to_state << " : ";
		o << " +";

		o << " [" << position_in_sequence << "] ";
		
		for (i=-2; i<3; i++)
		{
		    if (((position_in_sequence+i) >= 0) &&
			((position_in_sequence+i) < length_of_sequence))
		    {
			if (i != 0)
			{o << "[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]";}
			else
			{o << " ->[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]<- ";}
		    }
		}
		o << " \t" << posterior_prob << "\n" << flush;
	    }
	}
	if (count == 0) {
	    o << "NOTE: special transition from state " << number_of_from_state << " -> "
	      << number_of_to_state << " has been implemented, but may never be used.\n" << flush;
	}
    }
    return;
}

int Sequence::get_number_of_features_in_sequence(const int source,
						 const char* const file_name,
						 const int feature_label_index,
						 const int feature,
						 model_parameters* const MP)
{
    int check = 0;
    
    if((source<0)||(source>=number_of_sources))
    {
	cout<<" ERROR class Sequence::get_number_of_features_in_sequence : number_of_sources("
	    <<number_of_sources<<") out of range[0,"<<number_of_sources-1<<"].\n"<<flush;
	check++;
    }
    if (feature < 0) 
    {
	cout << "ERROR class Sequence::get_number_of_features_in_sequence : feature is Undefined.\n" << flush;
	check++;
    }
    if((feature_label_index<0)||(feature_label_index >= MP->get_Total_Number_of_Annotation_Labels()))
    {
	cout << "ERROR class Sequence::get_number_of_features_in_sequence : feature_label_index("
	     <<feature_label_index<<") out of range[0,"<<MP->get_Total_Number_of_Annotation_Labels()-1
	     <<"].\n"<<flush;
	check++;
    }    
    if (file_name == NULL) 
    {
	cout << "ERROR class Sequence::get_number_of_features_in_sequence : gtf_file_name is NULL.\n" << flush;
	check++;
    }
    
    if (check == 0)
    {
	int i = 0;
	int    number_of_features             = 0; // incl. partial features
	int    tmp_number_of_sources          = 0;
	int**  start_positions                = NULL;
	int**  end_positions                  = NULL;
	int*   number_of_start_positions      = NULL;
	int*   number_of_end_positions        = NULL;
      
	// get feature start and end positions
	// ===================================
	
	check += get_positions_of_feature_start_and_end(// input
	    this,
	    file_name,
	    feature_label_index,
	    feature,
	    MP,
	    // output
	    &tmp_number_of_sources,
	    &number_of_start_positions,
	    &number_of_end_positions,
	    &start_positions,
	    &end_positions);
	
	if (check != 0) {
	    cout << "ERROR: Sequence::get_number_of_features_in_sequence: error occurred in function "
		 << "get_positions_of_feature_start_and_end.\n" << flush;}
	else {
	    
	    number_of_features = max(number_of_start_positions[source], number_of_end_positions[source]);

	}
	
	// release memory
	if(number_of_start_positions) delete[] number_of_start_positions;
	number_of_start_positions = NULL;
	if(number_of_end_positions) delete[] number_of_end_positions;
	number_of_end_positions = NULL;
	
	if(start_positions){
	    for(i =0; i<tmp_number_of_sources; i++){
		if(start_positions[i]) delete[] start_positions[i];
		start_positions[i] = NULL;
	    }
	    delete [] start_positions;
	}
	start_positions = NULL;

	if(end_positions){
	    for(i =0; i<tmp_number_of_sources; i++){
		if(end_positions[i]) delete[] end_positions[i];
		end_positions[i] = NULL;
	    }
	    delete[] end_positions;
	}
	end_positions = NULL;

	return(number_of_features);
    }
    else {// if check != 0
	return(0);
    }
}

void Sequence::print_non_Logzero_scores_with_sequence_info(const int number_of_from_state,
							   const int number_of_to_state,
							   std::ostream &o,
							   model_parameters* const MP) const
{
    int check = 0;

    if (number_of_from_state < 0)
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0)
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) == 0)
	{
	    cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : \n"
		 << "special transition number_of_from_state ("
		 << number_of_from_state << ") -> number_of_to_state (" << number_of_to_state
		 << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }

    if ((check == 0) && (scores_for_special_transitions.GetNumberofDimensions() > 0)) 
    {
	
	int  count = 0;
	int  i;
	int  position_in_sequence = 0;
	Prob score                = static_cast<Score>(0.0);
	
	int linear_index = 0;
	int index_of_special_transition = this->get_index_of_special_transition(number_of_from_state,
										number_of_to_state);
	
	for (position_in_sequence=0; position_in_sequence<length_of_sequence; position_in_sequence++)
	{
	    linear_index = index_of_special_transition * length_of_sequence + position_in_sequence;
	    score        = scores_for_special_transitions.GetElement(linear_index);
	    
	    if (score != Logzero)
	    {
		count++;
		
		o << number_of_from_state << " -> " << number_of_to_state << " : ";
		o << " +";		// forward
		
		o << " [" << position_in_sequence << "] ";
		
		for (i=-2; i<3; i++)
		{
		    if (((position_in_sequence+i) >= 0) &&
			((position_in_sequence+i) < length_of_sequence))
		    {
			if (i != 0)
			{
			    o << "[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]";
			}
			else
			{
			    o << " ->[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]<- ";
			}
		    }
		}
		o << " \t" << score << "\n" << flush;
	    }
	}
	if (count == 0) 
	{
	    o << "NOTE: special transition from state " << number_of_from_state << " -> "
	      << number_of_to_state << " has been implemented, but may never be used.\n" << flush;
	}
    }
    return;
}

void Sequence::print_non_zero_posterior_probs_with_sequence_info(const int number_of_from_state,
								 const int number_of_to_state,
								 const int position_in_sequence,
								 std::ostream &o,
								 model_parameters* const MP) const
{
    int check = 0;
    
    if (number_of_from_state < 0)
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0)
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_from_state ("
	      << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) == 0)
	{
	    cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info : \n"
		 << "special transition number_of_from_state ("
		 << number_of_from_state << ") -> number_of_to_state (" << number_of_to_state
		 << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }
    
    if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1)))
    {
	cout << "ERROR class Sequence::print_non_zero_posterior_probs_with_sequence_info: \n"
	     << "position_in_sequence ("
	     << position_in_sequence << ") out of range (0.." 
	     << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    
    if ((check == 0) && (posterior_probs_for_special_transitions.GetNumberofDimensions() > 0)) {
	
	int  i;
	
	int index_of_special_transition = this->get_index_of_special_transition(number_of_from_state,
										number_of_to_state);
	int linear_index = index_of_special_transition * length_of_sequence + position_in_sequence;
	
	Prob posterior_prob = posterior_probs_for_special_transitions.GetElement(linear_index); 
	
	o << number_of_from_state << " -> " << number_of_to_state << " : ";
	o << " +"; //forward
	
	o << " [" << position_in_sequence << "] ";
	
	for (i=-2; i<3; i++)
	{
	    if (((position_in_sequence+i) >= 0) &&
		((position_in_sequence+i) < length_of_sequence))
	    {
		if (i != 0)
		{o << "[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]";}
		else
		{o << " ->[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]<- ";}
	    }
	}
	o << " \t" << posterior_prob << "\n" << flush;
      
    }
    return;
}

void Sequence::print_non_Logzero_scores_with_sequence_info(const int number_of_from_state,
							   const int number_of_to_state,
							   const int position_in_sequence,
							   std::ostream &o,
							   model_parameters* const MP) const
{
    int check = 0;

    if (number_of_from_state < 0)
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_from_state ("
	     << number_of_from_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_to_state < 0)
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_to_state ("
	     << number_of_to_state << ") < 0.\n" << flush;
	check++;
    }
    if (number_of_special_transitions == 0)
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : this sequence has not \n"
	     << "been set up to implement special transitions.\n" << flush;
	check++;
    }
    else
    {
	int max_n_of_child_state = implementation_status_of_special_transitions.GetDimension(0);
	
	if (number_of_from_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_from_state ("
		 << number_of_from_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (number_of_to_state > (max_n_of_child_state-1))
	{
	    cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : number_of_to_state ("
		 << number_of_to_state << ") > max_n_of_child_state-1 (" << max_n_of_child_state-1 << ").\n" << flush;
	    check++;
	}
	
	if (this->get_implementation_status_of_special_transition(number_of_from_state,
								  number_of_to_state) == 0)
	{
	    cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info : \n"
		 << "special transition number_of_from_state ("
		 << number_of_from_state << ") -> number_of_to_state (" << number_of_to_state
		 << ") has not yet been implemented.\n" << flush;
	    check++;
	}
    }
    
    if ((position_in_sequence < 0) || (position_in_sequence > (length_of_sequence-1)))
    {
	cout << "ERROR class Sequence::print_non_Logzero_scores_with_sequence_info: \n"
	     << "position_in_sequence ("
	     << position_in_sequence << ") out of range (0.." 
	     << length_of_sequence-1 << ").\n" << flush;
	check++;
    }
    
    if ((check == 0) && (scores_for_special_transitions.GetNumberofDimensions() > 0)) {
	
	int  i;
	
	int index_of_special_transition = this->get_index_of_special_transition(number_of_from_state,
										number_of_to_state);
	int linear_index = index_of_special_transition * length_of_sequence + position_in_sequence;
	Prob score       = scores_for_special_transitions.GetElement(linear_index); 
	
	o << number_of_from_state << " -> " << number_of_to_state << " : ";
	o << " +";
	
	o << " [" << position_in_sequence << "] ";
	
	for (i=-2; i<3; i++)
	{
	    if (((position_in_sequence+i) >= 0) &&
		((position_in_sequence+i) < length_of_sequence))
	    {
		if (i != 0)
		{o << "[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]";}
		else
		{o << " ->[" << MP->get_Alphabet_name(this->letter(position_in_sequence+i)) << "]<- ";}
	    }
	}
	o << " \t" << score << "\n" << flush;
	
    }
    return;
}

int Sequence::calculate_abs_pos(const int             rel_start,
				const int             rel_end,
				int*      const abs_start,
				int*      const abs_end) const
{
    int check = 0;

    if (this == NULL){
	cout << "ERROR : Sequence::calculate_abs_pos : Sequence is NULL.\n" << flush;
	check++;
    }

    if (this->get_start_position() > this->get_end_position()){
	cout << "ERROR : Sequence::calculate_abs_pos : start_position of Sequence ("
	     << this->get_start_position() << ") > end_position ("
	     << this->get_end_position() << ").\n" << flush;
	check++;
    }
    if (rel_start > rel_end){
	cout << "ERROR : Sequence::calculate_abs_pos : rel_start (" << rel_start
	     << ") > rel_end (" << rel_end << ").\n" << flush;
	check++;
    }
    if ((rel_start < 0) || (rel_start > this->length()-1)){
	cout << "ERROR : Sequence::calculate_abs_pos : rel_start (" << rel_start
	     << ") is out of range [0," << this->length()-1 << "].\n" << flush;
	check++;
    }
    if ((rel_end < 0) || (rel_end > this->length()-1)){
	cout << "ERROR : Sequence::calculate_abs_pos : rel_end (" << rel_end
	     << ") is out of range [0," << this->length()-1 << "].\n" << flush;
	check++;
    }

    if (check == 0){	
	(*abs_start) = this->get_start_position() + rel_start;
	(*abs_end)   = this->get_start_position() + rel_end;
    }

    return(check);
}

int Sequence::read_special_file(// input
			    const char*     const file_name,
			    model_parameters*  const MP,
			    // output
			    Info* const info)
{
    int check=0;

    // check input 
    
    if (file_name == NULL)
    {
	cout << "ERROR: class Sequence::read_special_file: file_name is NULL.\n" << flush;
	check++;
    }
    if (this->get_ac() == NULL)
    {
	cout << "ERROR: class Sequence::read_special_file: ac is NULL.\n" << flush;
	check++;
    }
    if (this->get_start_position() > this->get_end_position())
    {
	cout << "ERROR: class Sequence::read_special_file: start_position (" << this->get_start_position()
	     << ") > end_position (" << this->get_end_position() << ").\n" << flush;
	check++;
    }
  
    if (info == NULL)
    {
	cout << "ERROR: class Sequence::read_special_file: info is NULL.\n" << flush;
	check++;
    }
    
    // check that variables for output are NULL
    
    if (check == 0)
    {

	int        file_number_of_allocated_lines  = 0;
	int        file_number_of_sources          = 0;
	int*       file_number_of_lines            = NULL;
	char***    file_seq_names                  = NULL;  
	char**     file_source_names               = NULL;  
        bool****   file_annotation_labels          = NULL;  
	double**** file_scores                     = NULL;
	int**      file_start_positions            = NULL;
	int**      file_end_positions              = NULL; 
	int position_offset                = 0;
	
	check += read_special_file_for_sequence(// input
	    file_name,
	    ac,
	    start_position,
	    end_position,
	    position_offset,
	    MP,
	    // output
	    &file_number_of_allocated_lines,
	    &file_number_of_sources,
	    &file_number_of_lines,
	    &file_seq_names,
	    &file_source_names,
	    &file_annotation_labels,
	    &file_scores,
	    &file_start_positions,
	    &file_end_positions);
	
	if (check != 0)
	{
	    cout << "ERROR: class Sequence::read_special_file: occurred in function read_special_file_for_sequence.\n" << flush;
	}
	
	if (check == 0)
	{
	    if (file_number_of_sources > 0) 
	    {

		Info new_info(file_number_of_sources,
			      file_number_of_lines,
			      MP,
			      check);
	
		if(check)
		{
		    cout << "ERROR: class Sequence::read_special_file: occurred in Info constructor.\n" << flush;
		}
	
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;

		int tmp_number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
		
		for(i=0; i<file_number_of_sources; i++)
		{				   
		    check+=new_info.set_info_source_names(i, file_source_names[i]); 
		    for (j=0; j<file_number_of_lines[i]; j++)
		    {
			check+=new_info.set_info_seq_names(i, j, file_seq_names[i][j]);    	      
			
			for(k=0; k<tmp_number_of_annotation_labels; k++){
			    for(l=0; l<MP->get_Annotation_Label_size(k); l++){
				check+=new_info.set_info_annotation_labels(i,j,k,l,file_annotation_labels[i][j][k][l]);
				check+=new_info.set_info_scores(i,j,k,l,file_scores[i][j][k][l]);       							
			    }
			}		    
			check+=new_info.set_start_positions    (i,j,file_start_positions[i][j]);
			check+=new_info.set_end_positions      (i,j,file_end_positions[i][j]);  
		    }
		
		}
		(*info) = new_info;

	    }
	    else // if number_of_sources = 0
	    {	
		Info new_info;
	      
		(*info) = new_info;
		
		cout << "WARNING: class Sequence::read_special_file: did not find any lines for sequence " 
		     << this->get_ac() << ".\n" << flush;
	    }

	}
			
	// delete memory
	
	check+= release_memory_for_labels(&file_number_of_sources,
					 &file_number_of_lines,
					 MP,
					 &file_seq_names,
					 &file_source_names,
					 &file_annotation_labels,
					 &file_scores,
					 &file_start_positions,
					 &file_end_positions);
	
	if(check)
	{
	    cout << "ERROR: class Sequence::read_special_file: occurred in release_memory_for_labels.\n" << flush;
	}       
    }
    return(check);
}

void Sequence::print(model_parameters* const MP,
		     std::ostream &o) const 
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    
    o << "sequence_type      = " << convert_int_to_sequence_type[sequence_type]
      << "\n" << flush;
    o << "ac                 = " << ac << "\n" << flush;
    o << "start_position     = " << start_position << "\n" << flush;
    o << "end_position       = " << end_position << "\n" << flush;
    o << "length_of_sequence = " << length_of_sequence << "\n" << flush;
    if (sequence) {
	o << "sequence = \n" << flush;
	for (i=0; i < length_of_sequence; i++) {
	    o << sequence[i];
	}
	o << "\n" << flush;
    }
    o << "constraint = "<<constraint<<"\n"<<flush;

    if(annotation_labels)
    {
	for(i=0;i<length_of_sequence;i++)
	{
	    for(j=0; j<number_of_sources; j++)
	    {
		o <<"source : "<<j<<"\t";
		for(k=0; k<number_of_annotation_labels; k++)
		{				
		    for(l=0; l<number_of_each_annotation_label[k]; l++){
			if(annotation_labels[i][j][k][l]){
			    o<<MP->get_Annotation_Label_setname(k)<<":";
			    o<<MP->get_Annotation_Label_name(k,l)<<"(";
			    o<<annotation_label_scores[i][j][k][l]<<")\t";
			}
		    }
		}
	    }
	}
    }
   
    o << "number_of_implemented_special_emissions = " 
	 << number_of_implemented_special_emissions << "\n" << flush;
    o << "number_of_special_emissions             = " 
	 << number_of_special_emissions << "\n" << flush;
    o << "implementation_status_of_special_emissions = \n" << flush;
    implementation_status_of_special_emissions.PrintwithIndicesRestrict(cout, 0);
    o << "indices_of_special_emissions = \n" << flush;
    o << "the special emission with index 0 is missing from this printout!\n" << flush;
    indices_of_special_emissions.PrintwithIndicesRestrict(cout, 0);
    o << "posterior_probs_for_special_emissions = \n" << flush;
    posterior_probs_for_special_emissions.PrintwithIndicesRestrict(cout, 0.0);
    o << "scores_for_special_emissions = \n" << flush;
    this->print_emission_scores_with_sequence_info(cout,
						   MP);
    
    o << "number_of_implemented_special_transitions = " 
	 << number_of_implemented_special_transitions << "\n" << flush;
    o << "number_of_special_transitions = " 
	 << number_of_special_transitions << "\n" << flush;
    o << "implementation_status_of_special_transitions = \n" << flush;
    implementation_status_of_special_transitions.PrintwithIndicesRestrict(cout, 0);
    o << "indices_of_special_transitions = \n" << flush;
    o << "the special transition with index 0 is missing from this printout!\n" << flush;
    indices_of_special_transitions.PrintwithIndicesRestrict(cout, 0);
    
    o << "posterior_probs_for_special_transitions = \n" << flush;
    this->print_non_zero_posterior_probs_with_sequence_info(cout,
							    MP);
    o << "scores_for_special_transitions = \n" << flush;
    this->print_non_Logzero_scores_with_sequence_info(MP,
						      cout);
    o << "priors_of_special_transitions                = \n" << flush;
    priors_of_special_transitions.PrintwithIndicesRestrict(cout, 0.0);
    return;
}

void Sequence::print_short(void) const 
{
    cout << "type " << convert_int_to_sequence_type[sequence_type]
	 << " | ac " << ac
	 << " | start " << start_position
	 << " | end " << end_position
	 << " | length " << length_of_sequence << "\n" << flush;
    return;
}

Sequence & Sequence::operator = (const Sequence &s)
{
    int check = 0;

    if (!s.sequence)
    {
	cout << "ERROR class Sequence::operator = : sequence NULL\n" << flush;
	check++;
    }
    if (s.length_of_sequence<1)
    {
	cout << "ERROR class Sequence::operator = : length_of_sequence (" << s.length_of_sequence << ") < 1\n" << flush;
	check++;
    }
    
    if (check == 0)
    {
	if ( this != &s)
	{
	    int i = 0;
	    int j = 0;
	    int k = 0;
	    int l = 0;
	    
	    sequence_type=0;   // Nosequence
	    start_position = 0;
	    end_position = 0;
	    
	    if (sequence) delete [] sequence;
	    sequence=NULL;
      
	    if (ac) delete [] ac;
	    ac=NULL;
	    if(number_of_each_annotation_label) delete[] number_of_each_annotation_label;
	    number_of_each_annotation_label = NULL;

	    if(annotation_labels){
		for(i=0; i<length_of_sequence; i++){
		    if(annotation_labels[i]){
			for(j=0;j<number_of_sources;j++){
			    if(annotation_labels[i][j]){
				for(k=0; k<number_of_annotation_labels; k++){
				    if(annotation_labels[i][j][k]) delete [] annotation_labels[i][j][k];
				    annotation_labels[i][j][k] = NULL;
				}
				delete [] annotation_labels[i][j];
			    }
			    annotation_labels[i][j] = NULL;
			}
			delete[] annotation_labels[i];
		    }
		    annotation_labels[i] = NULL;
		}
		delete [] annotation_labels;
	    }
	    annotation_labels = NULL;
	    
	    if(annotation_label_scores){
		for(i=0; i<length_of_sequence; i++){
		    if(annotation_label_scores[i]){
			for(j=0;j<number_of_sources;j++){
			    if(annotation_label_scores[i][j]){
				for(k=0; k<number_of_annotation_labels; k++){
				    if(annotation_label_scores[i][j][k]) delete [] annotation_label_scores[i][j][k];
				    annotation_label_scores[i][j][k] = NULL;
				}
				delete [] annotation_label_scores[i][j];
			    }
			    annotation_label_scores[i][j] = NULL;
			}
			delete[] annotation_label_scores[i];
		    }
		    annotation_label_scores[i] = NULL;
		}
		delete [] annotation_label_scores;
	    }
	    annotation_label_scores = NULL;

	    number_of_sources = 0;
	    length_of_sequence=0;	  
	    number_of_annotation_labels = 0;
	    
	    if(score_labels) delete [] score_labels;
	    score_labels = NULL;
	    
	    constraint = 0;
	    
	    number_of_implemented_special_transitions=0;
	    number_of_special_transitions=0;
	  
	    number_of_implemented_special_emissions=0;
	    number_of_special_emissions=0;
	    
	    implementation_status_of_special_transitions.SetNumberofDimensions(0);
	    indices_of_special_transitions.SetNumberofDimensions(0);
	    posterior_probs_for_special_transitions.SetNumberofDimensions(0);
	    priors_of_special_transitions.SetNumberofDimensions(0);
	    scores_for_special_transitions.SetNumberofDimensions(0);
	    
	    implementation_status_of_special_emissions.SetNumberofDimensions(0);
	    indices_of_special_emissions.SetNumberofDimensions(0);
	    posterior_probs_for_special_emissions.SetNumberofDimensions(0);
	    scores_for_special_emissions.SetNumberofDimensions(0);
	    
	    if (s.sequence      &&              
		s.length_of_sequence>0 )
	    {
		// get sequence, its length, start and end positions and orientation

		sequence_type = s.sequence_type;
		start_position = s.start_position;
		end_position = s.end_position;
	      
		number_of_sources = s.number_of_sources;
		length_of_sequence=s.length_of_sequence;      
		sequence=new int[length_of_sequence];
		for (i=0; i<length_of_sequence; i++)
		{
		    sequence[i]=s.sequence[i];
		}
	      
		i=0;

		number_of_implemented_special_transitions=s.number_of_implemented_special_transitions;
		number_of_special_transitions=s.number_of_special_transitions;
		number_of_implemented_special_emissions=s.number_of_implemented_special_emissions;
		number_of_special_emissions=s.number_of_special_emissions;
		
		implementation_status_of_special_transitions=
		    s.implementation_status_of_special_transitions;
		indices_of_special_transitions=s.indices_of_special_transitions;
		priors_of_special_transitions=s.priors_of_special_transitions;
		posterior_probs_for_special_transitions=s.posterior_probs_for_special_transitions;
		scores_for_special_transitions=s.scores_for_special_transitions;
		
		implementation_status_of_special_emissions=
		    s.implementation_status_of_special_emissions;
		indices_of_special_emissions=s.indices_of_special_emissions;
		posterior_probs_for_special_emissions=s.posterior_probs_for_special_emissions;
		scores_for_special_emissions=s.scores_for_special_emissions;
		
		if (s.ac && (strlen(s.ac)>0)) 
		{
		    int length=0;
		    length=strlen(s.ac);
		    ac=new char[length+1];              
		    for (i=0; i<length+1; i++) {
			ac[i]=s.ac[i];
		    }
		}
		number_of_annotation_labels = s.number_of_annotation_labels;
		
		if(s.number_of_each_annotation_label){
		    number_of_each_annotation_label = new int[number_of_annotation_labels];
		    for(i=-0; i<number_of_annotation_labels; i++){
			number_of_each_annotation_label[i] = s.number_of_each_annotation_label[i];
		    }
		}
		
		if(s.annotation_labels)
		{
		    annotation_labels       = new bool   ***[length_of_sequence];
		    annotation_label_scores = new double ***[length_of_sequence]; 
		    for(i=0; i<length_of_sequence; i++)
		    {
			annotation_labels[i]       = new bool   **[number_of_sources];
			annotation_label_scores[i] = new double **[number_of_sources];
			for(j=0;j<number_of_sources; j++)
			{
			    annotation_labels[i][j]       = new bool   *[number_of_annotation_labels];
			    annotation_label_scores[i][j] = new double *[number_of_annotation_labels];
			    for(k=0; k<number_of_annotation_labels; k++){
				annotation_labels[i][j][k]       = new bool   [number_of_each_annotation_label[k]];
				annotation_label_scores[i][j][k] = new double [number_of_each_annotation_label[k]];
				for(l=0; l<number_of_each_annotation_label[k]; l++){
				    annotation_labels[i][j][k][l] = s.annotation_labels[i][j][k][l];
				    annotation_label_scores[i][j][k][l] = s.annotation_label_scores[i][j][k][l];
				}
			    }
			}
		    }
		}
		
		if(s.score_labels!=NULL)
		{
		    score_labels = new bool[number_of_annotation_labels];
		    
		    for(i=0; i<number_of_annotation_labels;i++)
		    {
			score_labels[i] = s.score_labels[i];
		    }
		}
		constraint = s.constraint;
				
	    }
	}
    }
    return(*this);
}
