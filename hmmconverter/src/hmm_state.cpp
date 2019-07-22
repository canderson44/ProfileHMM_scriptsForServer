/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: header-file for state-class
*/

#include "hmm_state.h"

// constructors and destructor
// ----------------------------------------------------------------------

Hmm_State::Hmm_State()
{
  // this constructor shall not be used

    special = 0;
    special_emission = 0;
    mirrored=0;
    alphabet=0;
    number_of_state=0;
    number_of_states=0;
    number_of_state_labels = 0;
    state_labels = NULL;
    state_labels_type = NULL;
    
    number_of_letters_to_read_x=0;
    number_of_letters_to_read_y=0;
    number_of_letters_to_read = 0;
    number_of_previous_states=0;
    number_of_special_transitions_to_previous_states=0;
    number_of_next_states=0;
    number_of_special_transitions_to_next_states=0;
  
    emission_probs_expression = NULL;
    emission_get_from = -1;
    emission_sum_over = false;
    emission_product = false;
    num_sum_over=0;
    sum_over_FromState = -1;
    sum_over_FromParam = -1;
    number_of_products = 0;
    product_FromState = NULL;
    product_FromParam = NULL;
    product_ThisPos = NULL;
    product_FromPos = NULL;    
}

Hmm_State::Hmm_State(const int n_of_letter_to_read_x, 
		     const int n_of_letter_to_read_y,
		     const int alph, 
		     const int n_of_state,
		     const int n_of_states,
		     const int n_of_previous_states,
		     const array<int> ns_of_previous_states,
		     const int n_of_state_labels,
		     array<int>* const s_labels,
		     char** const s_labels_type,
		     model_parameters* const MP)
{
    int check = 0; 
    int size = -1;
    int i = 0;
    // check alphabet
   
    if(alph!=MP->get_Alphabet_size()) {
	cout << "ERROR class Hmm_State::constructor : alphabet was not recognized.\n";
	check++;
    }
    // check number of states
    if (n_of_states<2)
    {cout << "ERROR class Hmm_State::constructor : number of states n = " << n_of_states
	  << "<2 (need at least 2 for Start and End state)\n";
    check++;}
    // check number of state
    if ((n_of_state <0) || (n_of_state > (n_of_states-1)))
    {cout << "ERROR class Hmm_State::constructor : number of state (" << n_of_state 
	  << ") out of range (0.." << n_of_states-1 << ").\n";
    check++;}
    // check number of previous states to this state
    if ((n_of_previous_states>(n_of_states-1)) || (n_of_previous_states<0))
    {cout << "ERROR class Hmm_State::constructor : number of previous states (" << n_of_previous_states
	  << ") has to < number_of_states (" << n_of_states << ") and >= 0.\n";
    check++;}
    // check the vector of number of previous states to this state
    if (ns_of_previous_states.GetNumberofDimensions() != 1)
    {cout << "ERROR class Hmm_State::constructor : number of dimensions of numbers_of_previous_states (" 
	  << ns_of_previous_states.GetNumberofDimensions() << ") != 1.\n";
    check++;}
    else
    {
	if ((n_of_previous_states>0) && (ns_of_previous_states.GetDimension(0) != n_of_previous_states))
	{cout << "ERROR class Hmm_State::constructor : length of array numbers_of_previous_states (" 
	      << ns_of_previous_states.GetDimension(0) << ") != number_of_previous_states (" 
	      << n_of_previous_states << ")\n";
	check++;}
	else
	{
	    // check that numbers of previous states are in the allowed range (0...number_of_states-1)
	    
	    for (int j=0; j<n_of_previous_states; j++)
	    {
		if ((ns_of_previous_states.GetElement(j) < 0) || 
		    (ns_of_previous_states.GetElement(j) > (n_of_states-1)))
		{cout << "ERROR class Hmm_State::constructor : number of previous state " << j << " ("
		      << ns_of_previous_states.GetElement(j) << ") out of range (0..."
		      << n_of_states-1 << ").\n";
		check++;}
	    }
	    
	    if (check == 0)
	    {
		// check that no previous state appears more than once
		
		for (int j=0; j<n_of_previous_states; j++)
		{
		    int state = ns_of_previous_states.GetElement(j);
		    int number_of_implementations = 0;
		    
		    for (int k=0; k<n_of_previous_states; k++)
		    {if (ns_of_previous_states.GetElement(k) == state) {number_of_implementations++;}}
		    
		    if (number_of_implementations != 1)
		    {cout << "ERROR class Hmm_State::constructor : number of previous state " << j << " ("
			  << ns_of_previous_states.GetElement(j) << ") is implemented more than once ("
			  << number_of_implementations << ").\n";
		    check++;}
		}
	    }
	}
    }
    //check state labels
    //number_of_dimensions is to calculate the total number of dimensions of the labels, xdim + ydim
    //add xread, yread to the state may be
    size=MP->get_Total_Number_of_Annotation_Labels();
    
    for(i=0;i<size;i++){
	if(s_labels[i].GetNumberofDimensions()!=1)
	{cout << "ERROR class Hmm_State::constructor : number of dimensions in "<<state_labels_type[i]<<"(" 
	      << s_labels[i].GetNumberofDimensions() << ") != 1.\n";
	check++;}
	
    }

    if (check == 0)
    {
   
	mirrored=0;
	
	alphabet=alph;	
	number_of_state=n_of_state;
	number_of_states=n_of_states;
	
	number_of_letters_to_read_x=n_of_letter_to_read_x;
	number_of_letters_to_read_y=n_of_letter_to_read_y;
	number_of_letters_to_read = number_of_letters_to_read_x + number_of_letters_to_read_y;
      
	number_of_previous_states=n_of_previous_states;
	numbers_of_previous_states=ns_of_previous_states;
	
	number_of_next_states = 0;
	number_of_state_labels = n_of_state_labels;

	//create state_labels
	state_labels = new array<int>[number_of_state_labels];
 
	//create state_labels_type      
	state_labels_type = new char*[number_of_state_labels];

	// set state_labels, state_labels_type
	for(i=0;i<size;i++){

	    state_labels[i] = s_labels[i];
	    state_labels_type[i]=Nullstrcpy(s_labels_type[i],check);
	    if(check!=0){
		cout<<"ERROR: Hmm_State constructor : "<<endl;
		cout<<"state_labels_type[i] : "<<s_labels_type[i]<<endl;
		break;
	    }
	}
	// set number of dimensions of probs (initialised with 0) and scores (initialised with Logzero)

	// For special model

	special = 0;
	special_emission = 0;
	number_of_child_state_x=0;
	number_of_child_state_y=0;
      
	number_of_previous_states=n_of_previous_states;
	number_of_special_transitions_to_previous_states=0;
	numbers_of_previous_states=ns_of_previous_states;
	special_flags_of_transitions_to_previous_states.SetNumberofDimensions(1);
	special_flags_of_transitions_to_previous_states.SetDimension(0, n_of_previous_states, 0);
	offset_for_previous_states_x.SetNumberofDimensions(0);
	offset_for_previous_states_y.SetNumberofDimensions(0);
	
	numbers_of_from_child_states_x_previous.SetNumberofDimensions(0);
	numbers_of_from_child_states_y_previous.SetNumberofDimensions(0);
	numbers_of_to_child_states_x_previous.SetNumberofDimensions(0);
	numbers_of_to_child_states_y_previous.SetNumberofDimensions(0); 
	
	if (number_of_letters_to_read_x!=0)
	{
	    offset_for_previous_states_x.SetNumberofDimensions(1);
	    offset_for_previous_states_x.SetDimension(0, n_of_previous_states, 0);      
	    
	    numbers_of_from_child_states_x_previous.SetNumberofDimensions(1);
	    numbers_of_to_child_states_x_previous.SetNumberofDimensions(1);
	    numbers_of_from_child_states_x_previous.SetDimension(0, n_of_previous_states, 0);      
	    numbers_of_to_child_states_x_previous.SetDimension(0, n_of_previous_states, 0);      
	}
	if (number_of_letters_to_read_y!=0)
	{
	    offset_for_previous_states_y.SetNumberofDimensions(1);
	    offset_for_previous_states_y.SetDimension(0, n_of_previous_states, 0);      
	    
	    numbers_of_from_child_states_y_previous.SetNumberofDimensions(1);
	    numbers_of_to_child_states_y_previous.SetNumberofDimensions(1);
	    numbers_of_from_child_states_y_previous.SetDimension(0, n_of_previous_states, 0);      
	    numbers_of_to_child_states_y_previous.SetDimension(0, n_of_previous_states, 0);      
	}
	number_of_special_transitions_to_next_states=0;
	
	transition_probs.SetNumberofDimensions(1);                             
	transition_probs.SetDimension(0,number_of_states, static_cast<Prob>(0.0));

	transition_probs_expressions.SetNumberofDimensions(1);                             
	transition_probs_expressions.SetDimension(0,number_of_states, static_cast<char*>(NULL));
	
	emission_probs.SetNumberofDimensions(number_of_letters_to_read); 
	{
	    for (int j=0; j<number_of_letters_to_read; j++)
	    {
		emission_probs.SetDimension(j,alphabet, static_cast<Prob>(0.0));
	    }
	}
       
	emission_probs_expression = NULL;
	emission_get_from = -1;
	emission_sum_over = false;
	emission_product = false;
	num_sum_over = 0 ;
	sum_over_FromState = -1;
	sum_over_FromParam = -1;
	number_of_products = 0;
	product_FromState = NULL;
	product_FromParam = NULL;
	product_ThisPos = NULL;
	product_FromPos = NULL;  
      
	transition_scores.SetNumberofDimensions(1);                            
	transition_scores.SetDimension(0, number_of_states, Logzero);  
      
	emission_scores.SetNumberofDimensions(number_of_letters_to_read);      
	{
	    for (int j=0; j<number_of_letters_to_read; j++){
		emission_scores.SetDimension(j, alphabet, Logzero);
	    }
	}
	
	viterbi_strip.SetNumberofDimensions(2);    
	viterbi_rectangle.SetNumberofDimensions(2);    
	viterbi_scores.SetNumberofDimensions(2);    
	
	this->order_next_and_previous_states();

    }
    if(check){
	this->~Hmm_State();
    }
    
}

Hmm_State::Hmm_State(const int n_of_letter_to_read_x, 
			     const int n_of_letter_to_read_y,
			     const int   alph,  
			     const int   n_of_state,
			     const int   n_of_states,
			     const int   n_of_previous_states,
			     const array<int> ns_of_previous_states,
			     const array<int> special_flags_to_previous_states,
			     const array<int> offset_previous_states_x,
			     const array<int> offset_previous_states_y,
			     const int n_of_from_child_state_x, 
			     const int n_of_from_child_state_y, 
			     const int n_of_to_child_state_x,   
			     const int n_of_to_child_state_y,   
			     const int   n_of_state_labels,
			     array<int>* const s_labels,
			     char** const s_labels_type,
			     model_parameters* const MP)
{
    int check = 0; 
    int size = -1;
    int i = 0;
    // check alphabet
   
    if(alph!=MP->get_Alphabet_size()) {
	cout << "ERROR class Hmm_State::constructor : alphabet was not recognized.\n";
	check++;
    }
    // check number of states
    if (n_of_states<2)
    {cout << "ERROR class Hmm_State::constructor : number of states n = " << n_of_states
	  << "<2 (need at least 2 for Start and End state)\n";
    check++;}
    // check number of state
    if ((n_of_state <0) || (n_of_state > (n_of_states-1)))
    {cout << "ERROR class Hmm_State::constructor : number of state (" << n_of_state 
	  << ") out of range (0.." << n_of_states-1 << ").\n";
    check++;}
    // check number of previous states to this state
    if ((n_of_previous_states>(n_of_states-1)) || (n_of_previous_states<0))
    {cout << "ERROR class Hmm_State::constructor : number of previous states (" << n_of_previous_states
	  << ") has to < number_of_states (" << n_of_states << ") and >= 0.\n";
    check++;}
    // check the vector of number of previous states to this state
    if (ns_of_previous_states.GetNumberofDimensions() != 1)
    {cout << "ERROR class Hmm_State::constructor : number of dimensions of numbers_of_previous_states (" 
	  << ns_of_previous_states.GetNumberofDimensions() << ") != 1.\n";
    check++;}
    else
    {
	if ((n_of_previous_states>0) && (ns_of_previous_states.GetDimension(0) != n_of_previous_states))
	{cout << "ERROR class Hmm_State::constructor : length of array numbers_of_previous_states (" 
	      << ns_of_previous_states.GetDimension(0) << ") != number_of_previous_states (" 
	      << n_of_previous_states << ")\n";
	check++;}
	else
	{
	    // check that numbers of previous states are in the allowed range (0...number_of_states-1)
	    
	    for (int j=0; j<n_of_previous_states; j++)
	    {
		if ((ns_of_previous_states.GetElement(j) < 0) || 
		    (ns_of_previous_states.GetElement(j) > (n_of_states-1)))
		{cout << "ERROR class Hmm_State::constructor : number of previous state " << j << " ("
		      << ns_of_previous_states.GetElement(j) << ") out of range (0..."
		      << n_of_states-1 << ").\n";
		check++;}
	    }
	    
	    if (check == 0)
	    {
		// check that no previous state appears more than once
		
		for (int j=0; j<n_of_previous_states; j++)
		{
		    int state = ns_of_previous_states.GetElement(j);
		    int number_of_implementations = 0;
		    
		    for (int k=0; k<n_of_previous_states; k++)
		    {if (ns_of_previous_states.GetElement(k) == state) {number_of_implementations++;}}
		    
		    if (number_of_implementations != 1)
		    {cout << "ERROR class Hmm_State::constructor : number of previous state " << j << " ("
			  << ns_of_previous_states.GetElement(j) << ") is implemented more than once ("
			  << number_of_implementations << ").\n";
		    check++;}
		}
	    }
	}
    }
    //check state labels
    //number_of_dimensions is to calculate the total number of dimensions of the labels, xdim + ydim
    //add xread, yread to the state may be
    size=MP->get_Total_Number_of_Annotation_Labels();
    
    for(i=0;i<size;i++){
	if(s_labels[i].GetNumberofDimensions()!=1)
	{cout << "ERROR class Hmm_State::constructor : number of dimensions in "<<state_labels_type[i]<<"(" 
	      << s_labels[i].GetNumberofDimensions() << ") != 1.\n";
	check++;}
	
    }

    // check input for special model

    if (special_flags_to_previous_states.GetNumberofDimensions() != 1)
    {
	cout << "ERROR class Hmm_State::constructor : number of dimensions in special_flags_to_previous_states (" 
	     << special_flags_to_previous_states.GetNumberofDimensions() << ") != 1.\n";
	check++;
    }
    if (special_flags_to_previous_states.GetDimension(0) != ns_of_previous_states.GetDimension(0))
    {
	cout << "ERROR class Hmm_State::constructor : length of array special_flags_to_previous_states (" 
	     << special_flags_to_previous_states.GetDimension(0) << ") != length of array ns_of_previous_states (" 
	     << ns_of_previous_states.GetDimension(0) << ")\n";
	check++;
    }
    if (n_of_letter_to_read_x!=0)
    {
	if (offset_previous_states_x.GetNumberofDimensions() != 1)
	{
	    cout << "ERROR class Hmm_State::constructor : number of dimensions in offset_previous_states_x (" 
		 << offset_previous_states_x.GetNumberofDimensions() << ") != 0.\n";
	    check++;
	}
	if (offset_previous_states_x.GetDimension(0) != ns_of_previous_states.GetDimension(0))
	{
	    cout << "ERROR class Hmm_State::constructor : length of array offset_previous_states_x (" 
		 << offset_previous_states_x.GetDimension(0) << ") != length of array ns_of_previous_states (" 
		 << ns_of_previous_states.GetDimension(0) << ")\n";
	    check++;
	}
    }
    else
    {
	if (offset_previous_states_x.GetNumberofDimensions() != 0)
	{
	    cout << "ERROR class Hmm_State::constructor : number of dimensions in offset_previous_states_x (" 
		 << offset_previous_states_x.GetNumberofDimensions() 
		 << ") should be 0 as this state does not read any letters from sequence x.\n";
	    check++;
	}
    }
    if (n_of_letter_to_read_y!=0)
    {
	if (offset_previous_states_y.GetNumberofDimensions() != 1)
	{
	    cout << "ERROR class Hmm_State::constructor : number of dimensions in offset_previous_states_y (" 
		 << offset_previous_states_y.GetNumberofDimensions() << ") != 1.\n";
	    check++;
	}
	if (offset_previous_states_y.GetDimension(0) != ns_of_previous_states.GetDimension(0))
	{
	    cout << "ERROR class Hmm_State::constructor : length of array offset_previous_states_y (" 
		 << offset_previous_states_y.GetDimension(0) << ") != length of array ns_of_previous_states (" 
		 << ns_of_previous_states.GetDimension(0) << ")\n";
	    check++;
	}
    }
    else
    {
	if (offset_previous_states_y.GetNumberofDimensions() != 0)
	{
	    cout << "ERROR class Hmm_State::constructor : number of dimensions in offset_previous_states_y (" 
		 << offset_previous_states_y.GetNumberofDimensions() 
		 << ") should be 0 as this state does not read any letters from sequence y.\n";
	    check++;
	}
    }

  // start new_trans

    if ((0 > n_of_from_child_state_x) || (n_of_from_child_state_x > (n_of_states-1))) {
	cout << "ERROR class Hmm_State::constructor : n_of_from_child_state_x (" << n_of_from_child_state_x
	     << ") out of range [0, " << n_of_states-1 << "].\n";
	check++;
    }
    if ((0 > n_of_from_child_state_y) || (n_of_from_child_state_y > (n_of_states-1))) {
	cout << "ERROR class Hmm_State::constructor : n_of_from_child_state_y (" << n_of_from_child_state_y
	     << ") out of range [0, " << n_of_states-1 << "].\n";
	check++;
    }
    if ((0 > n_of_to_child_state_x) || (n_of_to_child_state_x > (n_of_states-1))) {
	cout << "ERROR class Hmm_State::constructor : n_of_to_child_state_x (" << n_of_to_child_state_x
	     << ") out of range [0, " << n_of_states-1 << "].\n";
	check++;
    }
    if ((0 > n_of_to_child_state_y) || (n_of_to_child_state_y > (n_of_states-1))) {
	cout << "ERROR class Hmm_State::constructor : n_of_to_child_state_y (" << n_of_to_child_state_y
	     << ") out of range [0, " << n_of_states-1 << "].\n";
	check++;
    }
    if (((n_of_letter_to_read_x!=0)&&(n_of_letter_to_read_y==0)) && (n_of_from_child_state_y != 0)) {
	cout << "ERROR class Hmm_State::constructor : state is of type EmitX, n_of_from_child_state_y (" 
	     << n_of_from_child_state_y << ") has to be 0.\n";
	check++;
    }
    if (((n_of_letter_to_read_x!=0)&&(n_of_letter_to_read_y==0)) && (n_of_to_child_state_y != 0)) {
	cout << "ERROR class Hmm_State::constructor : state is of type EmitX, n_of_to_child_state_y (" 
	     << n_of_to_child_state_y << ") has to be 0.\n";
	check++;
    }
    if (((n_of_letter_to_read_y!=0)&&(n_of_letter_to_read_x==0)) && (n_of_from_child_state_x != 0)) {
	cout << "ERROR class Hmm_State::constructor : state is of type EmitY, n_of_from_child_state_x (" 
	     << n_of_from_child_state_x << ") has to be 0.\n";
	check++;
    }
    if (((n_of_letter_to_read_y!=0)&&(n_of_letter_to_read_x==0)) && (n_of_to_child_state_x != 0)) {
	cout << "ERROR class Hmm_State::constructor : state is of type EmitY, n_of_to_child_state_x (" 
	     << n_of_to_child_state_x << ") has to be 0.\n";
	check++;
    }

    // check vectors for previous special states
    if (check == 0)
    {
	int g;
	
	for (g=0; g<special_flags_to_previous_states.GetDimension(0); g++)
	{
	    // check that the entries in array special_flags_to_previous_states are 0 or 1
	    
	    if ((special_flags_to_previous_states.GetElement(g) != 0) &&
		(special_flags_to_previous_states.GetElement(g) != 1))
	    {
		cout << "ERROR class Hmm_State::constructor : element g (" << g 
		     << ") of array special_flags_to_previous_states ("
		     << special_flags_to_previous_states.GetElement(g) 
		     << ") is neither 0 nor 1.\n";
		check++;
	    }
	    else
	    {
		// check that if entry in array special_flags_to_previous_states is 0
		// that the corresponding entry or entries in offset_previous_states_x/y
		// array(s) are 0 
		
		if (special_flags_to_previous_states.GetElement(g) == 0)
		{
		    if  ((n_of_letter_to_read_x!=0) &&
			 (offset_previous_states_x.GetElement(g) != 0))
		    {
			cout << "ERROR class Hmm_State::constructor : element g (" << g 
			     << ") of array offset_previous_states_x ("
			     << offset_previous_states_x.GetElement(g) 
			     << ") must be 0 as the corresponding entry of array special_flags_to_previous_states ("
			     << special_flags_to_previous_states.GetElement(g) 
			     << ") is 0.\n";
			check++;
		    }
		    if ((n_of_letter_to_read_y!=0) &&
			(offset_previous_states_y.GetElement(g) != 0))
		    {
			cout << "ERROR class Hmm_State::constructor : element g (" << g 
			     << ") of array offset_previous_states_y ("
			     << offset_previous_states_y.GetElement(g) 
			     << ") must be 0 as the corresponding entry of array special_flags_to_previous_states ("
			     << special_flags_to_previous_states.GetElement(g) 
			     << ") is 0.\n";
			check++;
		    }
		}

		// check that if entry in array special_flags_to_previous_states is 1
		// that the corresponding entry or entries in offset_previous_states_x/y
		// array(s) is within the allowed range (0-delta_x-1) and (0-delta_y-1)
		
		else if (special_flags_to_previous_states.GetElement(g) == 1)
		{
		    if ((n_of_letter_to_read_x!=0) &&
			((offset_previous_states_x.GetElement(g) < 0) ||
			 (offset_previous_states_x.GetElement(g) > (n_of_letter_to_read_x-1))))
		    {
			cout << "ERROR class Hmm_State::constructor : element g (" << g 
			     << ") of array offset_previous_states_x ("
			     << offset_previous_states_x.GetElement(g) 
			     << ") out of range (0..." << (n_of_letter_to_read_x-1) << ").\n";
			check++;
		    }
		    if ((n_of_letter_to_read_y!=0) &&
			((offset_previous_states_y.GetElement(g) < 0) ||
			 (offset_previous_states_y.GetElement(g) > (n_of_letter_to_read_y-1))))
		    {
			cout << "ERROR class Hmm_State::constructor : element g (" << g 
			     << ") of array offset_previous_states_y ("
			     << offset_previous_states_y.GetElement(g) 
			     << ") out of range (0..." << (n_of_letter_to_read_y-1) << ").\n";
			check++;
		    }
		}
	    }
	}
    }

    if (check == 0)
    {
	special = 1 ;
	special_emission= 0;
	mirrored=0;
	
	alphabet=alph;
	number_of_state=n_of_state;
	number_of_child_state_x=0;
	number_of_child_state_y=0;
	number_of_states=n_of_states;

	number_of_letters_to_read_x=n_of_letter_to_read_x;
	number_of_letters_to_read_y=n_of_letter_to_read_y;
	number_of_letters_to_read = number_of_letters_to_read_x + number_of_letters_to_read_y;

	// start new_trans

	numbers_of_from_child_states_x_previous.SetNumberofDimensions(0);
	numbers_of_from_child_states_y_previous.SetNumberofDimensions(0);
	numbers_of_to_child_states_x_previous.SetNumberofDimensions(0);
	numbers_of_to_child_states_y_previous.SetNumberofDimensions(0); 
	
	if (number_of_letters_to_read_x!=0) {
	    
	    numbers_of_from_child_states_x_previous.SetNumberofDimensions(1); 
	    numbers_of_to_child_states_x_previous.SetNumberofDimensions(1);  
	    numbers_of_from_child_states_x_previous.SetDimension(0, n_of_previous_states, 0); 
	    numbers_of_to_child_states_x_previous.SetDimension(0, n_of_previous_states, 0);  
	}

	if (number_of_letters_to_read_y!=0) {
	    
	    numbers_of_from_child_states_y_previous.SetNumberofDimensions(1);  
	    numbers_of_to_child_states_y_previous.SetNumberofDimensions(1);   
	    numbers_of_from_child_states_y_previous.SetDimension(0, n_of_previous_states, 0);  
	    numbers_of_to_child_states_y_previous.SetDimension(0, n_of_previous_states, 0);   
	}
	
	// end new_trans
	
	number_of_previous_states=n_of_previous_states;
	numbers_of_previous_states=ns_of_previous_states;

	int g=0;
	for (g=0; g<special_flags_to_previous_states.GetDimension(0); g++) {
	    if (special_flags_to_previous_states.GetElement(g) == 1) {
		
		if (number_of_letters_to_read_x!=0) {
		    numbers_of_from_child_states_x_previous.SetElement(g, n_of_from_child_state_x);
		    numbers_of_to_child_states_x_previous.SetElement(g, n_of_to_child_state_x);
		}
		if (number_of_letters_to_read_y!=0) {
		    numbers_of_from_child_states_y_previous.SetElement(g, n_of_from_child_state_y);
		    numbers_of_to_child_states_y_previous.SetElement(g, n_of_to_child_state_y);
		}
		number_of_special_transitions_to_previous_states++;
	    }
	}
	
	numbers_of_previous_states=ns_of_previous_states;
	special_flags_of_transitions_to_previous_states=special_flags_to_previous_states;
	offset_for_previous_states_x = offset_previous_states_x;
	offset_for_previous_states_y = offset_previous_states_y;

	number_of_special_transitions_to_next_states=0;
	number_of_next_states = 0;
	number_of_state_labels = n_of_state_labels;

	//create state_labels
	state_labels = new array<int>[number_of_state_labels];
 
	//create state_labels_type      
	state_labels_type = new char*[number_of_state_labels];

	// set state_labels, state_labels_type
	for(i=0;i<size;i++){

	    state_labels[i] = s_labels[i];
	    state_labels_type[i]=Nullstrcpy(s_labels_type[i],check);
	    if(check!=0){
		cout<<"ERROR: Hmm_State constructor : "<<endl;
		cout<<"state_labels_type[i] : "<<s_labels_type[i]<<endl;
		break;
	    }
	}
	// set number of dimensions of probs (initialised with 0) and scores (initialised with Logzero)
	
	transition_probs.SetNumberofDimensions(1);                             
	transition_probs.SetDimension(0,number_of_states, static_cast<Prob>(0.0));

	transition_probs_expressions.SetNumberofDimensions(1);                             
	transition_probs_expressions.SetDimension(0,number_of_states, static_cast<char*>(NULL));
	
	emission_probs.SetNumberofDimensions(number_of_letters_to_read); 
	{
	    for (int j=0; j<number_of_letters_to_read; j++)
	    {
		emission_probs.SetDimension(j,alphabet, static_cast<Prob>(0.0));
	    }
	}
       
	emission_probs_expression = NULL;
	emission_get_from = -1;
	emission_sum_over = false;
	emission_product = false;
	num_sum_over = 0 ;
	sum_over_FromState = -1;
	sum_over_FromParam = -1;
	number_of_products = 0;
	product_FromState = NULL;
	product_FromParam = NULL;
	product_ThisPos = NULL;
	product_FromPos = NULL;  
      
	transition_scores.SetNumberofDimensions(1);                            
	transition_scores.SetDimension(0, number_of_states, Logzero);  
      
	emission_scores.SetNumberofDimensions(number_of_letters_to_read);      
	{
	    for (int j=0; j<number_of_letters_to_read; j++){
		emission_scores.SetDimension(j, alphabet, Logzero);
	    }
	}
	
	viterbi_strip.SetNumberofDimensions(2);    
	viterbi_rectangle.SetNumberofDimensions(2);    
	viterbi_scores.SetNumberofDimensions(2);    
	
	this->order_next_and_previous_states();
    }
    if(check){
	this->~Hmm_State();
    }
}
    
Hmm_State::Hmm_State(const Hmm_State &p)
{
    *this=p;
    
    // the elements of the following arrays which hold result of a particular run 
    // of an algorithm will not be copied

    viterbi_strip.SetNumberofDimensions(2);    
    viterbi_rectangle.SetNumberofDimensions(2);   
    viterbi_scores.SetNumberofDimensions(2);    
}

Hmm_State::~Hmm_State(void)
{
    int i = 0;
    int j = 0;
    special = 0;
    special_emission = 0;
    mirrored=0;
    
    alphabet=0;
    number_of_state=0;
    number_of_child_state_x=0;
    number_of_child_state_y=0;
    number_of_states=0;
    
    number_of_letters_to_read=0;
    number_of_letters_to_read_x = 0;
    number_of_letters_to_read_y = 0;
    number_of_previous_states=0;
    number_of_special_transitions_to_previous_states=0;
    number_of_next_states=0;
    number_of_special_transitions_to_next_states=0;

    emission_probs_expression = NULL;
    emission_get_from = -1;
    emission_sum_over = false;
    emission_product = false;
    num_sum_over=0;
    sum_over_FromState = -1;
    sum_over_FromParam = -1;
    number_of_products = 0;
    if(product_FromState) delete [] product_FromState;
    product_FromState = NULL;
    if(product_FromParam) delete [] product_FromParam;
    product_FromParam = NULL;
    if(product_ThisPos) delete [] product_ThisPos;
    product_ThisPos = NULL;
    if(product_FromPos) delete [] product_FromPos;
    product_FromPos = NULL;
 
}

// private member functions
// ----------------------------------------------------------------------

int Hmm_State::add_connection_to_new_state(const int number_of_new_state,
					       const int new_next_or_new_previous, // next = 1, previous = -1
					       const int transition_to_new_state_special, 
					       const int number_of_from_child_state_x_new_state,
					       const int number_of_from_child_state_y_new_state,
					       const int number_of_to_child_state_x_new_state,
					       const int number_of_to_child_state_y_new_state,
					       const int offset_x_new_state,
					       const int offset_y_new_state)
{
    // - this function can be used to add a new previous or new next state to this state
    // - if this state is mirrored = 0:
    //      - can add special or non-special transition to new previous state
    //      - can add non-special transition to new next state
    // - if this state is mirrored = 1:
    //      - can add special or non-special transition to new next state
    //      - can add non-special transition to new previous state
    //
    // - this function only updates the following scalar or array privat variables of this state
    //
    //         int        special;  
    //        
    //         and (
    // 
    //         int number_of_next_states;                            
    //         int number_of_special_transitions_to_next_states;     
    //         array<int> numbers_of_next_states;                    
    //         array<int> special_flags_of_transitions_to_next_states
    //         array<int> offset_for_next_states_x;                  
    //         array<int> offset_for_next_states_y;                  
    //         array<int> numbers_of_from_child_states_x_next;       
    //         array<int> numbers_of_from_child_states_y_next;       
    //         array<int> numbers_of_to_child_states_x_next;         
    //         array<int> numbers_of_to_child_states_y_next;         
    //
    //         or
    //
    //         int        number_of_previous_states;                       
    //         int        number_of_special_transitions_to_previous_states;
    //         array<int> numbers_of_previous_states;
    //         array<int> special_flags_of_transitions_to_previous_states; 
    //         array<int> offset_for_previous_states_x;              
    //         array<int> offset_for_previous_states_y;              
    //         array<int> numbers_of_from_child_states_x_previous;   
    //         array<int> numbers_of_from_child_states_y_previous;   
    //         array<int> numbers_of_to_child_states_x_previous;     
    //         array<int> numbers_of_to_child_states_y_previous;     
    //
    //         )
    //        
    // - NOTE: this function does not update the array of transition_probs or transition_scores
    //         and emission_probs or emission_scores
    
    int check = 0;
    
    // check if number_of_new_state within range 
    
    if ((number_of_new_state < 0) || (number_of_new_state > (number_of_states-1))) {
	cout << "ERROR: Hmm_State::add_connection_to_new_state : number_of_new_state("
	     << number_of_new_state << ") out of range [0, " << number_of_state-1 << "].\n";
	check++;
    }
    
    // check if value of new_next_or_new_previous is valid
    
    if ((new_next_or_new_previous != 1) && (new_next_or_new_previous != -1)) {
	cout << "ERROR: Hmm_State::add_connection_to_new_state : new_next_or_new_previous ("
	     << new_next_or_new_previous << ") has to be either 1 (next) or -1 (previous).\n";
	check++;
    }
    
    // check if new_state is already implement as next or previous state
    
    if ((new_next_or_new_previous ==  1) && (this->is_state_next_state(number_of_new_state) == 1)) {
	cout << "ERROR: Hmm_State::add_connection_to_new_state : number_of_new_state ("
	     << number_of_new_state << ") is already implemented as next state.\n";
	check++;
    }
    if ((new_next_or_new_previous == -1) && (this->is_state_previous_state(number_of_new_state) == 1)) {
	cout << "ERROR: Hmm_State::add_connection_to_new_state : number_of_new_state ("
	     << number_of_new_state << ") is already implemented as previous state.\n";
	check++;
    }
    
    // check if request for implement new special next or previous state is compatible with
    // mirrored-ness of this state
    
    if ((this->mirrored == 1) && (new_next_or_new_previous == -1) && (transition_to_new_state_special == 1)) {
	
	cout << "ERROR: Hmm_State::add_connection_to_new_state : new special transition to new previous "
	     << "state cannot be implemented as this state is mirrored.\n";
	check++;
    }
    if ((this->mirrored == 0) && (new_next_or_new_previous == 1) && (transition_to_new_state_special == 1)) {
	
	cout << "ERROR: Hmm_State::add_connection_to_new_state : new special transition to new next "
	     << "state cannot be implemented as this state is not mirrored.\n";
	check++;
    }

    if (check == 0) {
	
	int i      = 0;
	
	// make copy of this state
	// ------------------------------------------------------------
	
	Hmm_State copy = *this;

	// update scalars
	// ------------------------------------------------------------
	
	// * udpate : special 
	
	if ((transition_to_new_state_special == 1) && (special == 0)) {
	    
	    this->special = 1;

	}
	
	// update arrays
	// ------------------------------------------------------------
	
	if (new_next_or_new_previous == 1) { // if new next state is to be implemented
	    
	    // * update : number_of_next_states
	    // * update : array numbers_of_next_states
	    
	    this->number_of_next_states++;
	    
	    this->numbers_of_next_states.SetNumberofDimensions(1); 
	    this->special_flags_of_transitions_to_next_states.SetNumberofDimensions(1);
	    
	    this->numbers_of_next_states.SetDimension(0, this->number_of_next_states, 0); 
	    this->special_flags_of_transitions_to_next_states.SetDimension(0, this->number_of_next_states, 0);

	    if (this->mirrored == 1) {
		if (this->get_letters_to_read_x() !=0){
		    
		    this->offset_for_next_states_x.SetNumberofDimensions(1);                  
		    this->numbers_of_from_child_states_x_next.SetNumberofDimensions(1);       
		    this->numbers_of_to_child_states_x_next.SetNumberofDimensions(1);         
		    
		    this->offset_for_next_states_x.SetDimension(0, this->number_of_next_states, 0);                  
		    this->numbers_of_from_child_states_x_next.SetDimension(0, this->number_of_next_states, 0);       
		    this->numbers_of_to_child_states_x_next.SetDimension(0, this->number_of_next_states, 0);         
		}
		if (this->get_letters_to_read_y() !=0){
		    
		    this->offset_for_next_states_y.SetNumberofDimensions(1);                  
		    this->numbers_of_from_child_states_y_next.SetNumberofDimensions(1);       
		    this->numbers_of_to_child_states_y_next.SetNumberofDimensions(1);         
		    
		    this->offset_for_next_states_y.SetDimension(0, this->number_of_next_states, 0);                  
		    this->numbers_of_from_child_states_y_next.SetDimension(0, this->number_of_next_states, 0);       
		    this->numbers_of_to_child_states_y_next.SetDimension(0, this->number_of_next_states, 0);         
		}
	    }
	    
	    for (i=0; i<this->number_of_next_states; i++) {
		
		if (i<this->number_of_next_states-1) {
		    this->numbers_of_next_states.SetElement(i, copy.numbers_of_next_states.GetElement(i));
		}
		else {
		    this->numbers_of_next_states.SetElement(i, number_of_new_state);
		}
		
	    }
	    
	    // * update : number_of_special_transitions_to_next_states
	    // * update : array special_flags_of_transitions_to_next_states
	    // * update : array offset_for_next_states_x
	    // * update : array offset_for_next_states_y
	    // * update : array numbers_of_from_child_states_x_next
	    // * update : array numbers_of_from_child_states_y_next
	    // * update : array numbers_of_to_child_states_x_next
	    // * update : array numbers_of_to_child_states_y_next             
	    
	    this->number_of_special_transitions_to_next_states += transition_to_new_state_special; 
	    
	    for (i=0; i<this->number_of_next_states; i++) {
		
		if (i<this->number_of_next_states-1) {
		    
		    this->special_flags_of_transitions_to_next_states.
			SetElement(i, copy.special_flags_of_transitions_to_next_states.GetElement(i));
		    
		    if (this->mirrored == 1) {
			if (this->get_letters_to_read_x()!=0) {
			    
			    this->offset_for_next_states_x.
				SetElement(i, copy.offset_for_next_states_x.GetElement(i));	      
			    this->numbers_of_from_child_states_x_next.
				SetElement(i, copy.numbers_of_from_child_states_x_next.GetElement(i));			 
			    this->numbers_of_to_child_states_x_next.
				SetElement(i, copy.numbers_of_to_child_states_x_next.GetElement(i));
			}
			if (this->get_letters_to_read_y()!=0) {
			    
			    this->offset_for_next_states_y.
				SetElement(i, copy.offset_for_next_states_y.GetElement(i));
			    this->numbers_of_from_child_states_y_next.
				SetElement(i, copy.numbers_of_from_child_states_y_next.GetElement(i));			 
			    this->numbers_of_to_child_states_y_next.
				SetElement(i, copy.numbers_of_to_child_states_y_next.GetElement(i));
			}
		    }
		}
		else {
		    
		    this->special_flags_of_transitions_to_next_states.
			SetElement(i, transition_to_new_state_special);
		    
		    if (this->mirrored == 1) {
			if (this->get_letters_to_read_x()!=0) {
			    
			    this->offset_for_next_states_x.
				SetElement(i, offset_x_new_state);
			    this->numbers_of_from_child_states_x_next.
				SetElement(i, number_of_from_child_state_x_new_state);
			    this->numbers_of_to_child_states_x_next.
				SetElement(i, number_of_to_child_state_x_new_state);
			}
			if (this->get_letters_to_read_y()!=0) {
	      
			    this->offset_for_next_states_y.
				SetElement(i, offset_y_new_state);
			    this->numbers_of_from_child_states_y_next.
				SetElement(i, number_of_from_child_state_y_new_state);
			    this->numbers_of_to_child_states_y_next.
				SetElement(i, number_of_to_child_state_y_new_state);
			}
		    }
		}
	    }
	} // if new_next_or_new_previous == 1
	
	else if (new_next_or_new_previous == -1) {
	    
	    // * update : number_of_previous_states
	    // * update : array numbers_of_previous_states
	    
	    this->number_of_previous_states++;
	    
	    this->numbers_of_previous_states.SetNumberofDimensions(1);                   
	    this->special_flags_of_transitions_to_previous_states.SetNumberofDimensions(1);
	    
	    this->numbers_of_previous_states.SetDimension(0, this->number_of_previous_states, 0);  
	    this->special_flags_of_transitions_to_previous_states.SetDimension(0, this->number_of_previous_states, 0);

	    if (this->mirrored == 0) {
		if (this->get_letters_to_read_x()!=0){
		    
		    this->offset_for_previous_states_x.SetNumberofDimensions(1);                  
		    this->numbers_of_from_child_states_x_previous.SetNumberofDimensions(1);       
		    this->numbers_of_to_child_states_x_previous.SetNumberofDimensions(1);         
		    
		    this->offset_for_previous_states_x.SetDimension(0, this->number_of_previous_states, 0); 
		    this->numbers_of_from_child_states_x_previous.SetDimension(0, this->number_of_previous_states, 0); 
		    this->numbers_of_to_child_states_x_previous.SetDimension(0, this->number_of_previous_states, 0);
		}
		if (this->get_letters_to_read_y()!=0){
		    
		    this->offset_for_previous_states_y.SetNumberofDimensions(1);                  
		    this->numbers_of_from_child_states_y_previous.SetNumberofDimensions(1);       
		    this->numbers_of_to_child_states_y_previous.SetNumberofDimensions(1);         
		    
		    this->offset_for_previous_states_y.SetDimension(0, this->number_of_previous_states, 0);
		    this->numbers_of_from_child_states_y_previous.SetDimension(0, this->number_of_previous_states, 0); 
		    this->numbers_of_to_child_states_y_previous.SetDimension(0, this->number_of_previous_states, 0);
		}
	    }
	    
	    for (i=0; i<this->number_of_previous_states; i++) {
		
		if (i<this->number_of_previous_states-1) {
		    this->numbers_of_previous_states.SetElement(i, copy.numbers_of_previous_states.GetElement(i));
		}
		else {
		    this->numbers_of_previous_states.SetElement(i, number_of_new_state);
		}
	    }

	    // * update : number_of_special_transitions_to_previous_states
	    // * update : array special_flags_of_transitions_to_previous_states
	    // * update : array offset_for_previous_states_x
	    // * update : array offset_for_previous_states_y
	    // * update : array numbers_of_from_child_states_x_previous
	    // * update : array numbers_of_from_child_states_y_previous
	    // * update : array numbers_of_to_child_states_x_previous
	    // * update : array numbers_of_to_child_states_y_previous             
	    
	    this->number_of_special_transitions_to_previous_states += transition_to_new_state_special; 
      
	    for (i=0; i<this->number_of_previous_states; i++) {
		
		if (i<this->number_of_previous_states-1) {
		    
		    this->special_flags_of_transitions_to_previous_states.
			SetElement(i, copy.special_flags_of_transitions_to_previous_states.GetElement(i));
		    
		    if (this->mirrored == 0) {
			if (this->get_letters_to_read_x()!=0) {
			    
			    this->offset_for_previous_states_x.
				SetElement(i, copy.offset_for_previous_states_x.GetElement(i));	      
			    this->numbers_of_from_child_states_x_previous.
				SetElement(i, copy.numbers_of_from_child_states_x_previous.GetElement(i));			 
			    this->numbers_of_to_child_states_x_previous.
				SetElement(i, copy.numbers_of_to_child_states_x_previous.GetElement(i));
			}
			if (this->get_letters_to_read_y()!=0){
			    
			    this->offset_for_previous_states_y.
				SetElement(i, copy.offset_for_previous_states_y.GetElement(i));
			    this->numbers_of_from_child_states_y_previous.
				SetElement(i, copy.numbers_of_from_child_states_y_previous.GetElement(i));			 
			    this->numbers_of_to_child_states_y_previous.
				SetElement(i, copy.numbers_of_to_child_states_y_previous.GetElement(i));
			}
		    }
		}
		else {
		    
		    this->special_flags_of_transitions_to_previous_states.
			SetElement(i, transition_to_new_state_special);
		    
		    if (this->mirrored == 0) {
			if (this->get_letters_to_read_x()!=0) {
			    
			    this->offset_for_previous_states_x.
				SetElement(i, offset_x_new_state);
			    this->numbers_of_from_child_states_x_previous.
				SetElement(i, number_of_from_child_state_x_new_state);
			    this->numbers_of_to_child_states_x_previous.
				SetElement(i, number_of_to_child_state_x_new_state);
			}
			if (this->get_letters_to_read_y()!=0) {
	      
			    this->offset_for_previous_states_y.
				SetElement(i, offset_y_new_state);
			    this->numbers_of_from_child_states_y_previous.
				SetElement(i, number_of_from_child_state_y_new_state);
			    this->numbers_of_to_child_states_y_previous.
				SetElement(i, number_of_to_child_state_y_new_state);
			}
		    }
		}
	    }	       
	    
	} // else if new_next_or_new_previous == -1
	
	// finalise state
	// ------------------------------------------------------------
	
	this->order_next_and_previous_states();
	
    }
    
    return(check);
}

void Hmm_State::order_next_and_previous_states()
{

    int i = 0;
    int k = 0;
   
    if (number_of_previous_states > 1){
	
	// copy numbers of states into array states
       
	int* states        = new int[number_of_previous_states];
	int* special_flags = new int[number_of_previous_states];
	int* offset_x      = new int[number_of_previous_states];
	int* offset_y      = new int[number_of_previous_states];
	
	int* from_child_x  = new int[number_of_previous_states];
	int* from_child_y  = new int[number_of_previous_states];
	int* to_child_x    = new int[number_of_previous_states];
	int* to_child_y    = new int[number_of_previous_states];
	
	for (i=0; i<number_of_previous_states; i++)
	{
	    states[i]        = numbers_of_previous_states.GetElement(i);
	    special_flags[i] = special_flags_of_transitions_to_previous_states.GetElement(i);

	    if ((mirrored == 0) && (number_of_letters_to_read_x!=0)) {
		offset_x[i]     = offset_for_previous_states_x.GetElement(i);
		from_child_x[i] = numbers_of_from_child_states_x_previous.GetElement(i);
		to_child_x[i]   = numbers_of_to_child_states_x_previous.GetElement(i);
	    }
	    if ((mirrored == 0) && (number_of_letters_to_read_y!=0)) {
		offset_y[i] = offset_for_previous_states_y.GetElement(i);
		from_child_y[i] = numbers_of_from_child_states_y_previous.GetElement(i);
		to_child_y[i]   = numbers_of_to_child_states_y_previous.GetElement(i);
	    }
	}
	// order numbers of states in array states by increasing number

	int copy_state_number = 0;
	int copy_special_flag = 0;
	int copy_offset_x     = 0;
	int copy_offset_y     = 0;

	int copy_from_child_x = 0;
	int copy_from_child_y = 0;
	int copy_to_child_x   = 0;
	int copy_to_child_y   = 0;
	
	int min_state_number  = 0;
	
	int min_state_number_index = 0;
	
	for (i=0; i<number_of_previous_states-1; i++)
	{
	    // modified to selection sort
	    // look in all states with index [i+1] -> [final_index]
	    // for the smallest state number which is smaller than 
	    // states[i] 
	  
	    min_state_number = states[i];
	    min_state_number_index = i;
	    
	    for (k=i+1; k<number_of_previous_states; k++) {
		if (states[k] < min_state_number) {		  
		    min_state_number = states[k];
		    min_state_number_index = k;
		}
	    }
	    copy_state_number = states[i];
	    copy_special_flag = special_flags[i];

	    if ((mirrored == 0) && (number_of_letters_to_read_x!=0))
	    {
		copy_offset_x  = offset_x[i];
		copy_from_child_x = from_child_x[i];
		copy_to_child_x   = to_child_x[i];
	    }
	    if ((mirrored == 0) && (number_of_letters_to_read_y!=0))
	    {
		copy_offset_y  = offset_y[i];
		copy_from_child_y = from_child_y[i];
		copy_to_child_y   = to_child_y[i];
	    }
	    
	    states[i] =  states[min_state_number_index];
	    special_flags[i]  = special_flags[min_state_number_index];

	    if ((mirrored == 0) && (number_of_letters_to_read_x!=0)) {
		offset_x[i]    = offset_x[min_state_number_index];
		from_child_x[i] = from_child_x[min_state_number_index];
		to_child_x[i]   = to_child_x[min_state_number_index];
	    }
	    if ((mirrored == 0) && (number_of_letters_to_read_y!=0)) {
		offset_y[i]    = offset_y[min_state_number_index];
		from_child_y[i] = from_child_y[min_state_number_index];
		to_child_y[i]   = to_child_y[min_state_number_index];
	    }
	  
	    states[min_state_number_index] = copy_state_number;

	    special_flags[min_state_number_index]  = copy_special_flag;
	    if ((mirrored == 0) && (number_of_letters_to_read_x!=0)) {
		offset_x[min_state_number_index]    = copy_offset_x;
		from_child_x[min_state_number_index] = copy_from_child_x;
		to_child_x[min_state_number_index]   = copy_to_child_x;
	    }
	    if ((mirrored == 0) && (number_of_letters_to_read_y!=0)) {
		offset_y[min_state_number_index]    = copy_offset_y;
		from_child_y[min_state_number_index] = copy_from_child_y;
		to_child_y[min_state_number_index]   = copy_to_child_y;
	    }
	    
	}  
		
	// copy ordered numbers of states from array states back into array numbers_of_previous_states
	
	for (i=0; i<number_of_previous_states; i++)
	{
	    numbers_of_previous_states.SetElement(i, states[i]);
	    
	    special_flags_of_transitions_to_previous_states.SetElement(i, special_flags[i]);

	    if ((mirrored == 0) && (number_of_letters_to_read_x!=0)) {
		offset_for_previous_states_x.SetElement(i, offset_x[i]);
		numbers_of_from_child_states_x_previous.SetElement(i, from_child_x[i]);
		numbers_of_to_child_states_x_previous.SetElement(i, to_child_x[i]);
	    }
	    if ((mirrored == 0) && (number_of_letters_to_read_y!=0)) {
		offset_for_previous_states_y.SetElement(i, offset_y[i]);
		numbers_of_from_child_states_y_previous.SetElement(i, from_child_y[i]);
		numbers_of_to_child_states_y_previous.SetElement(i, to_child_y[i]);
	    }
	    
	}
	
	if (states) delete [] states;
	states = NULL;

	if (special_flags) delete [] special_flags;
	special_flags = NULL;
	if (offset_x) delete [] offset_x;
	offset_x = NULL;
	if (offset_y) delete [] offset_y;
	offset_y = NULL;
	
	if (from_child_x) delete [] from_child_x;
	from_child_x = NULL;
	if (from_child_y) delete [] from_child_y;
	from_child_y= NULL;
	if (to_child_x) delete [] to_child_x;
	to_child_x = NULL;
	if (to_child_y) delete [] to_child_y;
	to_child_y = NULL;
    }

    if (number_of_next_states > 1)
    {
	
	// copy numbers of states into array states
	
	int* states        = new int[number_of_next_states];
	int* special_flags = new int[number_of_next_states];
	int* offset_x      = new int[number_of_next_states];
	int* offset_y      = new int[number_of_next_states];
	
	int* from_child_x  = new int[number_of_next_states];
	int* from_child_y  = new int[number_of_next_states];
	int* to_child_x    = new int[number_of_next_states];
	int* to_child_y    = new int[number_of_next_states];
	
	for (i=0; i<number_of_next_states; i++)
	{
	    states[i]        = numbers_of_next_states.GetElement(i);
	    if(special_flags_of_transitions_to_next_states.GetNumberofDimensions()!=0){
		special_flags[i] = special_flags_of_transitions_to_next_states.GetElement(i);
		if ((mirrored == 1) && (number_of_letters_to_read_x!=0)) {
		    offset_x[i] = offset_for_next_states_x.GetElement(i);
		    from_child_x[i] = numbers_of_from_child_states_x_next.GetElement(i);
		    to_child_x[i]   = numbers_of_to_child_states_x_next.GetElement(i);
		}
		if ((mirrored == 1) && (number_of_letters_to_read_y!=0)) {
		    offset_y[i] = offset_for_next_states_y.GetElement(i);
		    
		    from_child_y[i] = numbers_of_from_child_states_y_next.GetElement(i);
		    to_child_y[i]   = numbers_of_to_child_states_y_next.GetElement(i);
		}
	    }
	}
	
	// order numbers of states in array states by increasing number
	
	int copy_state_number = 0;
	int copy_special_flag = 0;
	int copy_offset_x     = 0;
	int copy_offset_y     = 0;
	
	int copy_from_child_x = 0;
	int copy_from_child_y = 0;
	int copy_to_child_x   = 0;
	int copy_to_child_y   = 0;
	
	int min_state_number  = 0;

	int min_state_number_index = 0;
	
	for (i=0; i<number_of_next_states-1; i++)
	{
	    // look in all states with index [i+1] -> [final_index]
	    // for the smallest state number which is smaller than 
	    // states[i]
	    // modified to selection sort

	    min_state_number = states[i];
	    min_state_number_index = i;
	    
	    for (k=i+1; k<number_of_next_states; k++)
	    {
		if (states[k] < min_state_number)
		{
		    min_state_number = states[k];
		    min_state_number_index = k;
		}
		
	    }	  
	    copy_state_number = states[i];
	    if(special_flags_of_transitions_to_next_states.GetNumberofDimensions()!=0){
		copy_special_flag = special_flags[i];
		if ((mirrored == 1) && (number_of_letters_to_read_x!=0)) {
		    copy_offset_x  = offset_x[i];
		    
		    copy_from_child_x = from_child_x[i];
		    copy_to_child_x   = to_child_x[i];
		}
		if ((mirrored == 1) && (number_of_letters_to_read_y!=0)) {
		    copy_offset_y  = offset_y[i];
		    
		    copy_from_child_y = from_child_y[i];
		    copy_to_child_y   = to_child_y[i];
		}
	    }

	    states[i]         = states[min_state_number_index];
	    if(special_flags_of_transitions_to_next_states.GetNumberofDimensions()!=0){
		special_flags[i]  = special_flags[min_state_number_index];
		if ((mirrored == 1) && (number_of_letters_to_read_x!=0 )) {
		    offset_x[i]    = offset_x[min_state_number_index];
		    from_child_x[i] = from_child_x[min_state_number_index];
		    to_child_x[i]   = to_child_x[min_state_number_index];
		}
		if ((mirrored == 1) && (number_of_letters_to_read_y!=0)) {
		    offset_y[i]    = offset_y[min_state_number_index];
		    from_child_y[i] = from_child_y[min_state_number_index];
		    to_child_y[i]   = to_child_y[min_state_number_index];
		}
	    }
	    states[min_state_number_index]         = copy_state_number;
	    if(special_flags_of_transitions_to_next_states.GetNumberofDimensions()!=0){
		special_flags[min_state_number_index]  = copy_special_flag;
		if ((mirrored == 1) && (number_of_letters_to_read_x!=0)){
		    offset_x[min_state_number_index]    = copy_offset_x;
		    from_child_x[min_state_number_index] = copy_from_child_x;
		    to_child_x[min_state_number_index]   = copy_to_child_x;
		}
		if ((mirrored == 1) && (number_of_letters_to_read_y!=0)) {
		    offset_y[min_state_number_index]    = copy_offset_y;
		    from_child_y[min_state_number_index] = copy_from_child_y;
		    to_child_y[min_state_number_index]   = copy_to_child_y;
		}
	    }
	}
	
	// copy ordered numbers of states from array states back into array numbers_of_next_states
	
	for (i=0; i<number_of_next_states; i++) {
	    numbers_of_next_states.SetElement(i, states[i]);
	    if(special_flags_of_transitions_to_next_states.GetNumberofDimensions()!=0){
		special_flags_of_transitions_to_next_states.SetElement(i, special_flags[i]);
		
		if ((mirrored == 1) && (number_of_letters_to_read_x!=0)){
		    offset_for_next_states_x.SetElement(i, offset_x[i]);
		    numbers_of_from_child_states_x_next.SetElement(i, from_child_x[i]);
		    numbers_of_to_child_states_x_next.SetElement(i, to_child_x[i]);
		}
		if ((mirrored == 1) && (number_of_letters_to_read_y!=0)) {
		    offset_for_next_states_y.SetElement(i, offset_y[i]);
		    numbers_of_from_child_states_y_next.SetElement(i, from_child_y[i]);
		    numbers_of_to_child_states_y_next.SetElement(i, to_child_y[i]);
		}
	    }
	}
	if (states) delete [] states;
	states = NULL;

	if (special_flags) delete [] special_flags;
	special_flags = NULL;
	if (offset_x) delete [] offset_x;
	offset_x = NULL;
	if (offset_y) delete [] offset_y;
	offset_y = NULL;
	
	if (from_child_x) delete [] from_child_x;
	from_child_x = NULL;
	if (from_child_y) delete [] from_child_y;
	from_child_y= NULL;
	if (to_child_x) delete [] to_child_x;
	to_child_x = NULL;
	if (to_child_y) delete [] to_child_y;
	to_child_y = NULL;
    }
    return;
}

int Hmm_State::permute_x_and_y_indices_of_emission_probs(const int x_first_or_y_first) 
{

    int check = 0;
    
    if ((x_first_or_y_first != 1) && (x_first_or_y_first != -1)) {
	cout << "Hmm_State::permute_x_and_y_indices_of_emission_probs: x_first_or_y_first (" << x_first_or_y_first 
	     << ") has to be 1 (if the x indices come first in the old indices) or -1 (if the y indices come first "
	     << "in the old indices).\n";
	check++;
    }
    
    if (check == 0) {
      
	// permute x-indices <-> y-indices 
	// emission indices: (x_1,x_2,x_3,y_1,y_2,y_3) => (y_1,y_2,y_3,x_1,x_2,x_3)

	int i   = 0;
	int j   = 0;
    
	int d   = this->emission_probs.GetNumberofDimensions();
	int max = 0;
    
	if (d > 0) {
      
	    max= static_cast<int>(pow( static_cast<float>(alphabet), static_cast<float>(d)));
      
	    array<int> old_index(1);
	    old_index.SetDimension(0, d);
	    array<int> new_index(1);
	    new_index.SetDimension(0, d);
	    
	    Prob prob = 0;
	    
	    // copy array of emission probs before modifying it
      
	    array<Prob> copy_emission_probs = emission_probs;
	    
	    for (j=0; j<max; j++) { // loop over all possible indices
		
		prob=0;
		old_index.ResetData();
		convert_to_base(j, alphabet, &old_index);
		prob= copy_emission_probs.GetElement(old_index);
		
		new_index.ResetData();
		new_index.SetDimension(0, d);
		
		// if the x indices come first in the old indices
		
		if (x_first_or_y_first == 1) { 
		    
		    for (i=0; i<this->get_letters_to_read_x(); i++) {
			new_index.SetElement(this->get_letters_to_read_y()+i, old_index.GetElement(i));
		    }
		    for (i=0; i<this->get_letters_to_read_y(); i++) {
			new_index.SetElement(i, old_index.GetElement(this->get_letters_to_read_x()+i));

		    }
		}
		// if the y indices come first in the old indices
		
		else if (x_first_or_y_first == -1) { 
		    
		    for (i=0; i<this->get_letters_to_read_y(); i++) {
			new_index.SetElement(this->get_letters_to_read_x()+i, old_index.GetElement(i));
		    }	
		    for (i=0; i<this->get_letters_to_read_x(); i++) {
			new_index.SetElement(i, old_index.GetElement(this->get_letters_to_read_y()+i));
			
		    }	
		}
		
		this->set_emission_prob(new_index, prob);
		
	    }
	}
    }  // if check == 0
    
    return(check);
}

int Hmm_State::set_info_on_transitions_to_next_states(const array<int> ns_of_next_states,
							  const array<int> special_flags_to_next_states)
{
    
  int check = 0;

  if (mirrored == 1)
  {
      cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
	   << "this is a mirrored state. The function can be used only on un-mirrored states.\n";
      check++;
  }

  // check the number of dimensions of array ns_of_next_states 

  if (ns_of_next_states.GetNumberofDimensions() != 1)
  {
      cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
	   << "number of dimensions in ns_of_next_states (" 
	   << ns_of_next_states.GetNumberofDimensions() << ") != 1.\n";
      check++;
  }
  else
  {
      if (ns_of_next_states.GetDimension(0) < 0)
      {
	  cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
	       << "array ns_of_next_states must have at least >= 0 next states, but it has (" 
	       << ns_of_next_states.GetDimension(0) << ") next states. \n";
	  check++;
      }
  }
  
  if (check == 0)
  {
      // check that numbers of next states are in the allowed range (0...number_of_states-1)

      int j;

      for (j=0; j<ns_of_next_states.GetDimension(0); j++)
      {
	  if ((ns_of_next_states.GetElement(j) < 0) || 
	      (ns_of_next_states.GetElement(j) > (number_of_states-1)))
	  {
	      cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
		   << "number of next state " << j << " (" << ns_of_next_states.GetElement(j) 
		   << ") out of range (0..." << number_of_states-1 << ").\n";
	      check++;
	  }
      }

      if (check == 0)
      {
	  // check that no next state appears more than once

	  for (j=0; j<ns_of_next_states.GetDimension(0); j++)
	  {
	      int state = ns_of_next_states.GetElement(j);
	      int number_of_implementations = 0;
	      
	      for (int k=0; k<ns_of_next_states.GetDimension(0); k++)
	      {
		  if (ns_of_next_states.GetElement(k) == state)
		  {
		      number_of_implementations++;
		  }
	      }
	      
	      if (number_of_implementations != 1)
	      {
		  cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
		       << "number of next state " << j << " (" << ns_of_next_states.GetElement(j) 
		       << ") is implemented more than once (" << number_of_implementations << ").\n";
		  check++;
	      }
	  }
      }
  }

  // check the number of dimensions of array special_flags_to_next_states

  if (special_flags_to_next_states.GetNumberofDimensions() != 1)
    {cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
	  << "number of dimensions in special_flags_to_next_states (" 
	  << special_flags_to_next_states.GetNumberofDimensions() << ") != 1.\n";
    check++;}
  if (special_flags_to_next_states.GetDimension(0) != ns_of_next_states.GetDimension(0))
    {cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
	  << "length of array special_flags_to_next_states (" 
	  << special_flags_to_next_states.GetDimension(0) << ") != length of array ns_of_next_states (" 
	  << ns_of_next_states.GetDimension(0) << ")\n";
    check++;}

  // check that each entry in array special_flags_to_next_states is either 0 or 1

  if (check == 0)
  {
      int j;
      
      for (j=0; j<special_flags_to_next_states.GetDimension(0); j++)
      {
	  if ((special_flags_to_next_states.GetElement(j) != 0) &&
	      (special_flags_to_next_states.GetElement(j) != 1))
	  {cout << "ERROR class Hmm_State::set_info_on_transitions_to_next_states : \n"
		<< "element j (" << j << ") of array special_flags_to_next_states (" 
		<< special_flags_to_next_states.GetElement(j) << ") has to be either 0 or 1.\n";
	  check++;}
	}
  }
  
  if (check == 0)
  {
      number_of_next_states = ns_of_next_states.GetDimension(0);
      numbers_of_next_states = ns_of_next_states;
      
      number_of_special_transitions_to_next_states = 0;
      int i;
      for (i=0; i<special_flags_to_next_states.GetDimension(0); i++)
      {
	  if (special_flags_to_next_states.GetElement(i) == 1)
	  {number_of_special_transitions_to_next_states++;}
      }
      special_flags_of_transitions_to_next_states = special_flags_to_next_states;
      
      this->order_next_and_previous_states();
  }
  
  return(check);
  
}

// member functions
// ----------------------------------------------------------------------

int Hmm_State::modify_state(const int new_n_of_states,
			    const int shift_of_state_numbers,
			    const int n_of_fixed_state,
			    const int new_n_of_fixed_state)
{
    // NOTE : - if n_of_fixed_state == new_n_of_fixed_state, jump is set to 0
    //          and to 1 otherwise
    //        - for jump = 1 the state of number n_of_fixed_state is removed
    //          from the list of state and the references to this state are
    //          replaced by references to state of number new_n_of_fixed state
    //        - for jump = 0 the state n_of_fixed_state is treated the same as
    //          any other (non start, non end) state
    //        - for n_of_fixed_state != new_n_of_fixed_state neither of the
    //          numbers should refer to start or end states (this cannot be 
    //          checked within this function as has to be done by the user)
    
    int check = 0;

    if (new_n_of_states <= number_of_states) {
	cout << "ERROR: Hmm_State::modify_state: new_n_of_states (" << new_n_of_states
	     << ") <= number_of_states (" << number_of_states << ").\n";
	check++;
    }
    if (shift_of_state_numbers < 0) {
	cout << "ERROR: Hmm_State::modify_state: shift_of_state_number (" << shift_of_state_numbers
	     << ") < 0.\n";
	check++;
    }
    if ((n_of_fixed_state < 0) || (n_of_fixed_state > (number_of_states-1))) {
	cout << "ERROR: Hmm_State::modify_state: n_of_fixed_state (" << n_of_fixed_state
	     << ") out of range [0, " << number_of_states-1 << "].\n";
	check++;
    }
    if ((new_n_of_fixed_state < 0) || (new_n_of_fixed_state > (number_of_states + shift_of_state_numbers - 1))) {
	cout << "ERROR: Hmm_State::modify_state: new_n_of_fixed_state (" << new_n_of_fixed_state
	     << ") out of range [0, " << number_of_states + shift_of_state_numbers - 1 << "].\n";
	check++;
    }
    if ((shift_of_state_numbers > 0) &&
	(n_of_fixed_state != new_n_of_fixed_state) &&  // jump = 1
	((number_of_states + shift_of_state_numbers - 1) != new_n_of_states)) {
	cout << "ERROR: Hmm_State::modify_state: new_n_of_states (" << new_n_of_states
	     << ") != number_of_states (" << number_of_states << ") + shift_of_state_numbers ("
	     << shift_of_state_numbers << ") - 1 = " << number_of_states + shift_of_state_numbers - 1 << ".\n";
	check++;
    }
    if ((shift_of_state_numbers > 0) &&
	(n_of_fixed_state == new_n_of_fixed_state) &&  // jump = 0
	((number_of_states + shift_of_state_numbers) != new_n_of_states)) {
	cout << "ERROR: Hmm_State::modify_state: new_n_of_states (" << new_n_of_states
	     << ") != number_of_states (" << number_of_states << ") + shift_of_state_numbers ("
	 << shift_of_state_numbers << ") = " << number_of_states + shift_of_state_numbers - 1 << ".\n";
	check++;
    }
    
    if (check == 0) {
	
	int i      = 0;
	int j      = 0;
	int n      = 0;
	int k      = 0;
	int length = 0;
	int jump   = 0;  // decides whether to jump (1) the fixed state or not (0)

	// if jump == 1 : fixed state is removed from list of states
	//                and all references to fixed state are turned into
	//                references to state new_n_of_fixed_state
	
	
	if (n_of_fixed_state != new_n_of_fixed_state) {
	    jump = 1;
	}
	
	// ------------------------------------------------------------
	
	Hmm_State copy = *this;
	this->~Hmm_State();
    
	// new = old
	// ------------------------------------------------------------
	// scalars

	special                      = copy.special;
	special                      = copy.special_emission;
	mirrored                     = copy.mirrored;
	alphabet                     = copy.alphabet;
	
	number_of_state_labels = copy.number_of_state_labels;
	
	// copy the state_labels and state_labels_type
	
	if(state_labels) delete[] state_labels;
	state_labels = NULL;
	state_labels = new array<int>[number_of_state_labels];

	for(i=0;i<number_of_state_labels;i++){
	    state_labels[i] = copy.state_labels[i];
	}
	
	if(state_labels_type) delete[] state_labels_type;
	state_labels_type = NULL;
	state_labels_type = new char*[number_of_state_labels];
	
	for(i=0;i<number_of_state_labels;i++){
	    state_labels_type[i]=Nullstrcpy(copy.state_labels_type[i],check);
	}
	    
	number_of_letters_to_read    = copy.number_of_letters_to_read;
	number_of_letters_to_read_x = copy.number_of_letters_to_read_x;
	number_of_letters_to_read_y = copy.number_of_letters_to_read_y;
	
	number_of_previous_states    = copy.number_of_previous_states;    
	number_of_next_states        = copy.number_of_next_states;

	// arrays

	offset_for_previous_states_x = copy.offset_for_previous_states_x;
	offset_for_previous_states_y = copy.offset_for_previous_states_y;
	offset_for_next_states_x     = copy.offset_for_next_states_x;
	offset_for_next_states_y     = copy.offset_for_next_states_y;
	emission_probs               = copy.emission_probs;

	special_flags_of_transitions_to_previous_states  = 
	    copy.special_flags_of_transitions_to_previous_states;
	special_flags_of_transitions_to_next_states      =
	    copy.special_flags_of_transitions_to_next_states;
	number_of_special_transitions_to_previous_states = 
	    copy.number_of_special_transitions_to_previous_states;
	number_of_special_transitions_to_next_states     =
	    copy.number_of_special_transitions_to_next_states;
	
	// new
	// ------------------------------------------------------------
	// scalars
	
	number_of_states = new_n_of_states;
	
	if (copy.number_of_state == 0) {                                     // this == start state
	    
	    number_of_state = copy.number_of_state;
	    
	}
	else if (copy.number_of_state == (copy.number_of_states-1)) {        // this == end state
	    
	    number_of_state = new_n_of_states-1;
	}
	else if ((copy.number_of_state == n_of_fixed_state) && (jump == 1)) { // this == fixed state 
	    
	    number_of_state = new_n_of_fixed_state;
	}
	else {
	    
	    number_of_state = copy.number_of_state + shift_of_state_numbers;
	    
	    if ((jump == 1) && (copy.number_of_state > n_of_fixed_state)) {
		
		number_of_state -= 1;

	    }
	}
	
	if ((copy.number_of_child_state_x == n_of_fixed_state) &&
	    (jump == 1))                                      { // child_state == fixed state
	    
	    number_of_child_state_x = new_n_of_fixed_state;
	}
	else if (copy.number_of_child_state_x == 0) {           // child_state == default
	    
	    number_of_child_state_x = 0;
	}
	else if (number_of_child_state_x == (copy.number_of_states-1)) {   // this == end state
	    
	    number_of_child_state_x = new_n_of_states-1;
	}
	else {
	    
	    number_of_child_state_x = copy.number_of_child_state_x + shift_of_state_numbers;
	    
	    if ((jump == 1) && (copy.number_of_child_state_x > n_of_fixed_state)) {
		number_of_child_state_x -= 1;
	    }
	}
	
	if ((copy.number_of_child_state_y == n_of_fixed_state) &&
	    (jump == 1))                                      { // child_state == fixed state
      
	    number_of_child_state_y = new_n_of_fixed_state;
	}
	else if (copy.number_of_child_state_y == 0) {           // child_state == default
	    
	    number_of_child_state_y = 0;
	}
	else if (number_of_child_state_y == (copy.number_of_states-1)) {   // this == end state
	    
	    number_of_child_state_y = new_n_of_states-1;
	}
	else {
	    
	    number_of_child_state_y = copy.number_of_child_state_y + shift_of_state_numbers;
	    
	    if ((jump == 1) && (copy.number_of_child_state_y > n_of_fixed_state)) {
		number_of_child_state_y -= 1;
	    }
	}
	
	// 1-d arrays
	
	numbers_of_previous_states = copy.numbers_of_previous_states;
	
	if (copy.numbers_of_previous_states.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_previous_states.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_previous_states.GetElement(i);	
		k = n;
		
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
			
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		    
		}
		else if (n == (copy.number_of_states-1)) {

		    n = new_n_of_states-1;

		}
		numbers_of_previous_states.SetElement(i, n);
		
	    }
	}
	
	numbers_of_next_states = copy.numbers_of_next_states;
	
	if (copy.numbers_of_next_states.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_next_states.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_next_states.GetElement(i);	
		k = n;
		
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_next_states.SetElement(i, n);
		
	    }
	}

	numbers_of_from_child_states_x_previous = copy.numbers_of_from_child_states_x_previous;

	if (copy.numbers_of_from_child_states_x_previous.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_from_child_states_x_previous.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_from_child_states_x_previous.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_from_child_states_x_previous.SetElement(i, n);
	    }
	}
	
	numbers_of_from_child_states_y_previous = copy.numbers_of_from_child_states_y_previous;
	
	if (copy.numbers_of_from_child_states_y_previous.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_from_child_states_y_previous.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_from_child_states_y_previous.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_from_child_states_y_previous.SetElement(i, n);
	    }
	}
	
	numbers_of_to_child_states_x_previous = copy.numbers_of_to_child_states_x_previous;
	
	if (copy.numbers_of_to_child_states_x_previous.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_to_child_states_x_previous.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_to_child_states_x_previous.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_to_child_states_x_previous.SetElement(i, n);
	    }
	}
	
	numbers_of_to_child_states_y_previous = copy.numbers_of_to_child_states_y_previous;
	
	if (copy.numbers_of_to_child_states_y_previous.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_to_child_states_y_previous.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_to_child_states_y_previous.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_to_child_states_y_previous.SetElement(i, n);
	    }
	}
	
	numbers_of_from_child_states_x_next = copy.numbers_of_from_child_states_x_next;
	
	if (copy.numbers_of_from_child_states_x_next.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_from_child_states_x_next.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_from_child_states_x_next.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_from_child_states_x_next.SetElement(i, n);
	    }
	}
	
	numbers_of_from_child_states_y_next = copy.numbers_of_from_child_states_y_next;
	
	if (copy.numbers_of_from_child_states_y_next.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_from_child_states_y_next.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_from_child_states_y_next.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_from_child_states_y_next.SetElement(i, n);
	    }
	}
	
	numbers_of_to_child_states_x_next = copy.numbers_of_to_child_states_x_next;
	
	if (copy.numbers_of_to_child_states_x_next.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_to_child_states_x_next.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_to_child_states_x_next.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		
		    n = new_n_of_states-1;
		}
		numbers_of_to_child_states_x_next.SetElement(i, n);
	    }
	}
	
	numbers_of_to_child_states_y_next = copy.numbers_of_to_child_states_y_next;
	
	if (copy.numbers_of_to_child_states_y_next.GetNumberofDimensions() > 0) {
	    length = copy.numbers_of_to_child_states_y_next.GetDimension(0);
	    for (i=0; i<length; i++) {
		n = copy.numbers_of_to_child_states_y_next.GetElement(i);	
		k = n;
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (k > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {

		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		numbers_of_to_child_states_y_next.SetElement(i, n);
	    }
	}
	
	transition_probs = copy.transition_probs;
	if (copy.transition_probs.GetNumberofDimensions() > 0) {
	    
	    const int d   = copy.transition_probs.GetNumberofDimensions();
	    const int max = copy.transition_probs.GetDimension(0);
	    
	    transition_probs.SetDimension(0, new_n_of_states, 0);
	    
	    Prob prob = 0.0;
	    
	    for (j=0; j<max; j++) {
	    
		prob = copy.transition_probs.GetElement(j);
		
		n = j;
		
		if ((n != 0)                                                && // if state != start and
		    (n != (copy.number_of_states-1))                        && //    state != end   and
		    (! ((n == n_of_fixed_state) && (jump == 1)))) {            //    state != fixed state
		    
		    n += shift_of_state_numbers;
		    
		    if ((jump == 1) && (j > n_of_fixed_state)) {
			n -= 1;
		    }
		}
		else if ((n == n_of_fixed_state) && (jump == 1)) {
		    
		    n = new_n_of_fixed_state;
		}
		else if (n == (copy.number_of_states-1)) {
		    
		    n = new_n_of_states-1;
		}
		transition_probs.SetElement(n, prob);
				
	    }
	}
	
	// not copied but initialised
	// ------------------------------------------------------------
	//
	viterbi_strip.SetNumberofDimensions(2);    
	viterbi_rectangle.SetNumberofDimensions(2);    
	viterbi_scores.SetNumberofDimensions(2);    
	
	transition_scores.SetNumberofDimensions(1);                            
	transition_scores.SetDimension(0, number_of_states, Logzero);  
	
	const int n_of_dimensions = number_of_letters_to_read;
	emission_scores.SetNumberofDimensions(n_of_dimensions);      
	for (j=0; j<n_of_dimensions; j++) {
	    emission_scores.SetDimension(j, alphabet, Logzero);
	}
	
	// finalise state
	// ------------------------------------------------------------
	
	this->order_next_and_previous_states();
	
    }
    return(check);
}

// return 0 if no error
int Hmm_State:: set_alphabet(int const alph)
{
    alphabet = alph;
    return 0;
}

int Hmm_State:: set_mirrored(int const mir)
{
    mirrored = mir;
    return 0;
}

int Hmm_State:: set_number_of_state_labels(int const n_s_labels)
{
    number_of_state_labels = n_s_labels;
    return 0;
}

int Hmm_State:: set_state_labels_type(char** const s_label_type)
{
    int check = 0;
    if(number_of_state_labels<=0){
	cout<<"Error : Hmm_State:: set_state_labels_type, number_of_state_labels : "<<number_of_state_labels<<endl;
	check++;
	return check;
    }
    state_labels_type = new char*[number_of_state_labels];
    for(int i = 0; i<number_of_state_labels; i++){
	if(!s_label_type[i]){
	    cout<<"Error : Hmm_State:: set_state_labels_type, s_label_type["<<i<<"] is NULL "<<endl;
	    check++;
	    break;
	}else{
	    state_labels_type[i]=Nullstrcpy(s_label_type[i],check);
	    if(check){
		cout<<"Error : Hmm_State:: set_state_labels_type, input string invalid : "<<s_label_type<<endl;
		break;
	    }
	}
    }
    if(check){
	if(state_labels_type) delete[] state_labels_type;
	state_labels_type = NULL;
    }
    return check;
}
	
int Hmm_State:: set_state_labels_type(char* const s_label_type, const int i)
{
    int check = 0;
    if(!state_labels_type){
	state_labels_type = new char*[number_of_state_labels];
    }
    state_labels_type[i]=Nullstrcpy(s_label_type,check);
    if(check){
	cout<<"Error : Hmm_State:: set_state_labels_type, input string invalid : "<<s_label_type<<endl;
	if(state_labels_type) delete[] state_labels_type;
	state_labels_type = NULL;
    }
    return check;
}      

int Hmm_State::set_number_of_state(int ns) 
{
    int check = 0;
    if((ns<0)||(ns>number_of_states))
    {
	check++;
    }else{
	number_of_state = ns;
    }
    return check;
}

int Hmm_State::set_number_of_states(int n) 
{
    int check = 0;
    if(n<0)
    {
	check++;
    }else{
	number_of_states = n;
    }
    return check;
}


int Hmm_State::set_transition_probs_expression(const int j, char* const exp)
{
    int check = 0;
    
  if ((j<0) || (j>number_of_states))
    {cout << "ERROR class Hmm_State::set_transition_prob_expression index : j = " << j << " out of range. \n";
    check++;}

  if (transition_probs_expressions.GetNumberofDimensions()!=1 ||
      transition_probs_expressions.GetDimension(0)!=number_of_states)
    {cout << "ERROR class Hmm_State::set_transition_probs_expression : array transition_probs_expression does not exist. \n";
    check++;}

  if (check == 0)
  {
      int state_is_next_state = 0;

      for (int i=0; i<number_of_next_states; i++)
      {
	  if (numbers_of_next_states.GetElement(i) == j)
	  {
	      state_is_next_state++;
	      break;
	  }
      }

      if (state_is_next_state != 1) 
      {cout << "ERROR class Hmm_State::set_transition_probs_expression : cannot set transition prob expression to state (" << j 
	    << ") to expression = " << exp <<" as this state is no next state.\n";
      cout << "ERROR class Hmm_State::set_transition_probs_expression : number of next states = "
	   << number_of_next_states << "\n";
      cout << "ERROR class Hmm_State::set_transition_probs_expression : next states are : ";
      numbers_of_next_states.Print(cout);
      cout << "\n";
      check++;}
    }
  
  if (check == 0) {transition_probs_expressions.SetElement(j,exp);}
  return(check);
}

int Hmm_State::set_transition_prob(const int j, const Prob x)
{
  int check = 0;

  if ((j<0) || (j>number_of_states))
    {cout << "ERROR class Hmm_State::set_transition_prob index : j = " << j << " out of range. \n";
    check++;}
#ifndef _SPECIAL_TRANS
  if ((x<0) || (x>1))
    {cout << "ERROR class Hmm_State::set_transition_prob Prob : x = " << x << " out of range. \n";
    check++;}
#endif
  if (transition_probs.GetNumberofDimensions()!=1 ||
      transition_probs.GetDimension(0)!=number_of_states)
    {cout << "ERROR class Hmm_State::set_transition_prob : array transition_probs does not exist. \n";
    check++;}
  if (check == 0)
  {
      int state_is_next_state = 0;

      for (int i=0; i<number_of_next_states; i++)
      {
	  if (numbers_of_next_states.GetElement(i) == j)
	  {
	      state_is_next_state++;
	      break;
	  }
      }
      if (state_is_next_state != 1) 
      {
	  cout << "ERROR class Hmm_State::set_transition_prob : cannot set transition prob to state (" << j 
	       << ") to prob = " << x <<" as this state is no next state.\n";
	  cout << "ERROR class Hmm_State::set_transition_prob : number of next states = "
	       << number_of_next_states << "\n";
	  cout << "ERROR class Hmm_State::set_transition_prob : next states are : ";
	  numbers_of_next_states.Print(cout);
	  cout << "\n";
	  check++;
      }
  }
  
  if (check == 0) {transition_probs.SetElement(j,x);}
  return(check);
}

int Hmm_State::set_special_emissions(void) 
{
    int check = 0;
    if ((number_of_letters_to_read_x!=0) || (number_of_letters_to_read_y != 0)) {
	
	special_emission=1;
	number_of_child_state_x = number_of_state;
	number_of_child_state_y = number_of_state;
    }
    else {
	cout << "ERROR: class Hmm_State::set_special_emissions: this state is of wrong type, "
	     <<"it has to be EmitX or EmitY. Cannot enable special emissions.\n";
	check++;
    }
    return(check);
}

int Hmm_State::set_special_emissions(const Hmm_State* const child_state_x,
				     const Hmm_State* const child_state_y) 
{
    int check = 0;
    
    if (((this->get_letters_to_read_x() != 0)&&(this->get_letters_to_read_y() != 0))&&
	(child_state_x->get_letters_to_read_x() != 0) &&
	(child_state_y->get_letters_to_read_y() != 0)) 
    {	
	special_emission=1;
	number_of_child_state_x = child_state_x->get_number_of_state();
	number_of_child_state_y = child_state_y->get_number_of_state();
    }
    else {
	if ((this->get_letters_to_read_x() == 0)||(this->get_letters_to_read_y() == 0)) {
	    cout << "ERROR: class Hmm_State::set_special_emissions: this state is of wrong type "
		 <<"(it has to be EmitXY). Cannot enable special emissions.\n";
	    check++;
	}
	if (child_state_x->get_letters_to_read_x() == 0) {
	    cout << "ERROR: class Hmm_State::set_special_emissions: child_state_x is of wrong type "
		 << "(it has to be EmitX). Cannot enable special emissions for this state.\n";
	    check++;
	}
	if (child_state_y->get_letters_to_read_y() == 0) {
	    cout << "ERROR: class Hmm_State::set_special_emissions: child_state_y is of wrong type "
		 << "(it has to be EmitY). Cannot enable special emissions for this state.\n";
	    check++;
	}
    }
    return(check);
}

int Hmm_State::set_special_emissions(const Hmm_State* const child_state) 
{
    int check = 0;
    
    if ((((this->get_letters_to_read_x() != 0) && (this->get_letters_to_read_y() == 0)) // EmitX
	 || ((this->get_letters_to_read_y() != 0) && (this->get_letters_to_read_x() == 0)))// EmitY
	&&(((child_state->get_letters_to_read_x() != 0) && (child_state->get_letters_to_read_y() == 0)) //EmitX
	   ||((child_state->get_letters_to_read_y() != 0) && (child_state->get_letters_to_read_x() == 0)))) //EmitY
    {
	
	special_emission=1;
	number_of_child_state_x = child_state->get_number_of_state();
	number_of_child_state_y = child_state->get_number_of_state();
    }
    else {
	if ( ! (((this->get_letters_to_read_x() != 0)&&(this->get_letters_to_read_y()==0))
		|| ((this->get_letters_to_read_y() !=0)&&(this->get_letters_to_read_x()==0)))) {
	    cout << "ERROR: class Hmm_State::set_special_emissions: this state is of wrong type "
		 << "(it has to be EmitX or EmitY). Cannot enable special emissions.\n";
	    check++;
	}
	if ( ! (((child_state->get_letters_to_read_x() != 0)&&(child_state->get_letters_to_read_y()==0))
		|| ((child_state->get_letters_to_read_y() !=0)&&(child_state->get_letters_to_read_x()==0)))) {
	    cout << "ERROR: class Hmm_State::set_special_emissions: child_state is of wrong type "
		 << "(it has to be EmitX or EmitY). Cannot enable special emissions for this state.\n";
	    check++;
	}
    }
    return(check);
}

int  Hmm_State::set_special_emissions(const int n_of_child_state_x,
					  const int n_of_child_state_y) 
{
    int check = 0;
    
    if ((this->get_letters_to_read_x() != 0)&& (this->get_letters_to_read_y() != 0)){
	
	special_emission=1;
	number_of_child_state_x = n_of_child_state_x;
	number_of_child_state_y = n_of_child_state_y;
    }
    else {
	if (!((this->get_letters_to_read_x() != 0)&& (this->get_letters_to_read_y() != 0))) {
	    cout << "ERROR: class Hmm_State::set_special_emissions: this state is of wrong type "
		 << "(it has to be EmitXY). Cannot enable special emissions.\n";
	    check++;
	}
    }
    return(check);
}

int  Hmm_State::set_special_emissions(const int n_of_child_state) 
{
    int check = 0;
    
    if (((this->get_letters_to_read_x() != 0) && (this->get_letters_to_read_y() == 0)) // EmitX
	|| ((this->get_letters_to_read_y() != 0) && (this->get_letters_to_read_x() == 0))) //EmitY
    {
	special_emission=1;
	number_of_child_state_x = n_of_child_state;
	number_of_child_state_y = n_of_child_state;
    }
    else {
	if ( ! (((this->get_letters_to_read_x() != 0) && (this->get_letters_to_read_y() == 0)) // EmitX
		|| ((this->get_letters_to_read_y() != 0) && (this->get_letters_to_read_x() != 0)))) //EmitY
	{
	    cout << "ERROR: class Hmm_State::set_special_emissions: this state is of wrong type "
		 << "(it has to be EmitXY). Cannot enable special emissions.\n";
	    check++;
	}
    }
    return(check);
}

int Hmm_State::set_emission_probs_expression(const char* exp)
{
  int check = 0;

  if (check == 0) 
  {
      if(emission_probs_expression) delete[] emission_probs_expression;
      emission_probs_expression = NULL;
      emission_probs_expression=Nullstrcpy(exp,check);
      if(check){
	  cout<<"Error class Hmm_State:: set_emission_probs_expression, input expression : "<<exp<< " invalid "<<endl;
      }
  }
  return(check);
}

int Hmm_State::set_emission_prob(const array<int> &indices, const Prob x)
{
  int check = 0;

  if ( (x<0) || (x>1) )
  {
      cout << "ERROR class Hmm_State::set_emission_prob Prob : x = " << x << " out of range. \n";
      check++;
  }
  if (indices.GetNumberofDimensions()!=1)
  {
      cout << "ERROR class Hmm_State::set_emission_prob : indices is no vector\n";
      check++;
  }
  if (indices.GetDimension(0)!=emission_probs.GetNumberofDimensions())
  {
      cout << "ERROR class Hmm_State::set_emission_prob : length of vector indices, " 
	   << indices.GetDimension(0)
	   << ", does not match number of dimensions of array emission_probs, " 
	   << emission_probs.GetNumberofDimensions() << "\n";
      check++;
  }
  if (emission_probs.GetNumberofDimensions()!= number_of_letters_to_read )
  {
      cout << "ERROR class Hmm_State::set_emission_prob : array emission_probs does not exist. \n";
      check++;
  }
  if (check == 0)
  {
      // check that entries of indices are in correct range
      for (int i=0; i<emission_probs.GetNumberofDimensions(); i++)
      {
	  if ((indices.GetElement(i) < 0) || (indices.GetElement(i) > (emission_probs.GetDimension(i)-1)))
	  {
	      cout << "ERROR class Hmm_State::set_emission_prob : element " << i << " of indices (" 
		   << indices.GetElement(i) << ") out of range (0..."
		   << emission_probs.GetDimension(i)-1 << ").\n";
	      check++;
	  }
      }
  }
  
  if (check == 0) 
  {
      emission_probs.SetElement(indices,x);
  }
  return(check);
}

int Hmm_State::set_emission_prob(const int linear_index, const Prob x)
{
  int check = 0;
  
  if ( (x<0) || (x>1) )
    {cout << "ERROR class Hmm_State::set_emission_prob Prob : x = " << x << " out of range. \n";
    check++;}
  if (emission_probs.GetNumberofDimensions()!= number_of_letters_to_read )
    {cout << "ERROR class Hmm_State::set_emission_prob : array emission_probs does not exist. \n";
    check++;}
  int length=emission_probs.GetDimension(0);
  for (int i=1; i<emission_probs.GetNumberofDimensions(); i++)
    {length*=emission_probs.GetDimension(i);}
  if ((linear_index<0) || (linear_index>(length-1)))
    {cout << "ERROR class Hmm_State::set_emission_prob : linear_index " << linear_index
          << " out of range (0.." << length-1 << ").\n";
    check++;}

  if (check == 0) {emission_probs.SetElement(linear_index, x);}
  return(check);
}

int Hmm_State::set_emission_prob(const array<Prob> prob)
{
  int check = 0;
  int NumDim = prob.GetNumberofDimensions();
  int* Dim = NULL;
  long i =0;
  
  if (NumDim!= number_of_letters_to_read )
  {
      cout << "ERROR class Hmm_State::set_emission_prob : array emission_probs does not exist. \n";
      check++;
  }
  
  Dim = new int[NumDim];
  double Tprob=0;
  
  for(i=0; i<NumDim; i++){
      Dim[i] = prob.GetDimension(i);
  }

  long TNum = 1;
  for(i=0;i<NumDim;i++){
      TNum*=Dim[i];
  }
  for(i=0;i<TNum;i++){
      if((prob.GetElement(i)<0)||(prob.GetElement(i)>1)){
	  cout << "ERROR class Hmm_State::set_emission_prob : "
	       <<" prob.GetElement("<<i<<") : "<<prob.GetElement(i)<<endl;
	  check++;
	  break;
      }
      Tprob+=prob.GetElement(i);
  }
  
  if(abs(Tprob-1.0)>Max_deviation){
      cout << "ERROR class Hmm_State::set_emission_prob : "
	   <<" abs(Tprob-1.0) > Max_deviation, Tprob : "<<Tprob<<endl;
      check++;
  }

  if(Dim) delete[] Dim;
  Dim = NULL;

  if (check == 0){
      emission_probs = prob;
  }
  return(check);
}


void Hmm_State::set_emission_sum_over(bool input_sum_over){
  emission_sum_over=input_sum_over;
}

void Hmm_State::set_emission_product(bool input_product){
  emission_product = input_product;
}

void Hmm_State::set_num_sum_over(int input_num_sum){
  num_sum_over=input_num_sum;
}

int Hmm_State::set_emission_get_from(int from){
    int check = 0;
    if((from<0)||(from>=number_of_states)){
	cout<<"Error : Hmm_State:: set_emission_get_from, input not in range 0.."<<number_of_states<<endl;
	check++;
    }else{
	emission_get_from = from;
    }
    return check;
}

int Hmm_State::set_sum_over_ThisPos(const char* cur_pos)
{
    int i=0;
    int check =0;

    if(!check){
	
	int NumOfItems = 0;
	
	char** Items = new char*[Max_number_of_items];	
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	check+= splitstring(cur_pos,NumOfItems,&Items,' ');
	sum_over_ThisPos.SetNumberofDimensions(1);
	sum_over_ThisPos.SetDimension(0,NumOfItems);
	for(i=0; i<NumOfItems; i++)
	{
	    sum_over_ThisPos.SetElement(i,atoi(Items[i]));
	}
	
	for(i=0; i< Max_number_of_items; i++){
	    if(Items[i]) delete[] Items[i];
	    Items[i] = NULL;
	}
	if(Items) delete[] Items;
	Items = NULL;
    }else{
	cout<<"Error : Hmm_State:: set_sum_over_ThisPos, input cur_pos : "<<cur_pos<< " invalid "<<endl;
    }
    
    return(check);
}
  
int Hmm_State::set_sum_over_FromPos(const char* corr_pos){
    int i=0;
    int check =0;
    
    if(check==0){
	
	int NumOfItems = 0;
	char** Items = new char*[Max_number_of_items];	
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	check+= splitstring(corr_pos,NumOfItems,&Items,' ');
	sum_over_FromPos.SetNumberofDimensions(1);
	sum_over_FromPos.SetDimension(0,NumOfItems);
	for(i=0; i<NumOfItems; i++)
	{
	    sum_over_FromPos.SetElement(i,atoi(Items[i]));
	}
      
	for(i=0; i< Max_number_of_items; i++){
	    if(Items[i]) delete[] Items[i];
	    Items[i] = NULL;
	}
	if(Items) delete[] Items;
	Items = NULL;
	
    }else{
	cout<<"Error : Hmm_State:: set_sum_over_FromPos, input corr_pos : "<<corr_pos<< " invalid "<<endl;
    }
    
    return(check);
}

int Hmm_State::set_sum_over_pos(const char* sum_pos)
{
    int i=0;
    int check =0;
    
    if(check==0){
	
	int NumOfItems = 0;
	char** Items = new char*[Max_number_of_items];	
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	check+= splitstring(sum_pos,NumOfItems,&Items,' ');
	sum_over_pos.SetNumberofDimensions(1);
	sum_over_pos.SetDimension(0,NumOfItems);
	for(i=0; i<NumOfItems; i++)
	{
	    sum_over_pos.SetElement(i,atoi(Items[i]));
	}
	
	for(i=0; i< Max_number_of_items; i++){
	    if(Items[i]) delete[] Items[i];
	    Items[i] = NULL;
	}
	if(Items) delete[] Items;
	Items = NULL;
      
    }else{
	cout<<"Error : Hmm_State:: set_sum_over_pos, input sum_pos : "<<sum_pos<< " invalid "<<endl;
    } 
    
    return(check);
}

int Hmm_State::set_sum_over_FromState(int fromstate)
{
    int check = 0;
    if((fromstate<0)||(fromstate>=number_of_states)){
	cout<<"Error : Hmm_State:: set_sum_over_FromState, input not in range 0.."
	    <<number_of_states<<endl;
	check++;
    }else{
	sum_over_FromState= fromstate;
    }
    return check;
}

int Hmm_State::set_sum_over_FromParam(int fromparam)
{
    int check = 0;
    if(fromparam<0){
	cout<<"Error : Hmm_State:: set_sum_over_FromState, input out of range."<<endl;
	check++;
    }else{
	sum_over_FromParam= fromparam;
    }
    return check;
}

int Hmm_State::set_number_of_products(int n_of_products)
{
    int check = 0;
    if(n_of_products <= 0)
    {
	cout<<"Error : Hmm_State : set_number_of_products: n_of_products("<<n_of_products<<") out of range."<<endl;
	check++;
    }
    number_of_products = n_of_products;
    return check;
}

int Hmm_State::set_product_ThisPos(int index, const char* cur_pos)
{
    int i=0;
    int check =0;

    if(number_of_products<0)
    {
	cout<<"Hmm_State :: set_product_ThisPos, number_of_products("
	    <<number_of_products<<") <0."<<endl;
	check++;
    }  

    if(check==0)
    {

	if(!product_ThisPos)
	{
	    product_ThisPos = new array<int>[number_of_products];
	}
	
	int NumOfItems = 0;
	char** Items = new char*[Max_number_of_items];	
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	check+= splitstring(cur_pos,NumOfItems,&Items,' ');
	product_ThisPos[index].SetNumberofDimensions(1);
	product_ThisPos[index].SetDimension(0,NumOfItems);
	for(i=0; i<NumOfItems; i++)
	{
	    product_ThisPos[index].SetElement(i,atoi(Items[i]));
	}
	
	for(i=0; i< Max_number_of_items; i++){
	    if(Items[i]) delete[] Items[i];
	    Items[i] = NULL;
	}
	if(Items) delete[] Items;
	Items = NULL;
	
    }else{
	cout<<"Error : Hmm_State:: set_product_ThisPos, input cur_pos : "<<cur_pos<< " invalid "<<endl;
    }
  
    return(check);
}

int Hmm_State::set_product_FromPos(int index,const char* corr_pos)
{
    int i=0;
    int check =0;

    if(number_of_products<0)
    {
	cout<<"Hmm_State :: set_product_FromPos, number_of_products("
	    <<number_of_products<<") <0."<<endl;
	check++;
    }    

    if(check==0)
    {

	if(!product_FromPos)
	{	
	    product_FromPos = new array<int>[number_of_products];
	}
	
      	int NumOfItems = 0;
	
	char** Items = new char*[Max_number_of_items];	
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	check+= splitstring(corr_pos,NumOfItems,&Items,' ');
	product_FromPos[index].SetNumberofDimensions(1);
	product_FromPos[index].SetDimension(0,NumOfItems);
	for(i=0; i<NumOfItems; i++)
	{
	    product_FromPos[index].SetElement(i,atoi(Items[i]));
	}
	
	for(i=0; i< Max_number_of_items; i++){
	    if(Items[i]) delete[] Items[i];
	    Items[i] = NULL;
	}
	if(Items) delete[] Items;
	Items = NULL;
	
    }else{
	cout<<"Error : Hmm_State:: set_product_FromPos, input corr_pos : "<<corr_pos<< " invalid "<<endl;
    }
    
    return(check);	
}

int Hmm_State::set_product_FromState(int index, int fromstate){
    int check = 0;

    if(number_of_products<0)
    {
	cout<<"Hmm_State :: set_product_From, number_of_products("
	    <<number_of_products<<") <0."<<endl;
	check++;
    } 

    if(!check)
    {
	if(!product_FromState)
	{
	    product_FromState = new int[number_of_products];
	    for(int i=0; i<number_of_products; i++)
	    {
		product_FromState[i] = -1;
	    }
	}
	if((fromstate<0)||(fromstate>=number_of_states)){
	    cout<<"Error : Hmm_State:: set_product_FromState, input not in range 0.."<<number_of_states<<endl;
	    check++;
	}else{
	    product_FromState[index] = fromstate;
	}
    }
    return check;
}

int Hmm_State::set_product_FromParam(int index, int fromparam){
    int check = 0;

    if(number_of_products<0)
    {
	cout<<"Hmm_State :: set_product_From, number_of_products("
	    <<number_of_products<<") <0."<<endl;
	check++;
    } 

    if(!check)
    {
	if(!product_FromParam)
	{
	    product_FromParam = new int[number_of_products];
	    for(int i=0; i<number_of_products; i++)
	    {
		product_FromParam[i] = -1;
	    }
	}
	if(fromparam<0){
	    cout<<"Error : Hmm_State:: set_product_FromParam, input out of range."<<endl;
	    check++;
	}else{
	    product_FromParam[index] = fromparam;
	}
    }
    return check;
}

int Hmm_State::set_transition_score(const int j, const Score x)
{
    int check = 0;

    if ( (j<0) || (j>number_of_states) )
    {
	cout << "ERROR class Hmm_State::set_transition_score index : j = " << j << " out of range. \n";
	check++;
    }
    if ((transition_scores.GetNumberofDimensions()!=1)  ||
	(transition_scores.GetDimension(0)!=number_of_states))
    {
	cout << "ERROR class Hmm_State::set_transition_score : array transition_scores does not exist. \n";
	check++;
    }
    if (!check)
    {
	int state_is_previous_state = 0;

	for (int i=0; i<number_of_previous_states; i++)
	{
	    if (numbers_of_previous_states.GetElement(i) == j)
	    {
		state_is_previous_state++;
		break;
	    }
	}
	if (state_is_previous_state != 1) 
	{
	    cout << "ERROR class Hmm_State::set_transition_score: state (" << j << ") is no previous state.\n";
	    cout << "ERROR class Hmm_State::set_transition_score : number of previous states = "
		 << number_of_previous_states << "\n";
	    cout << "ERROR class Hmm_State::set_transition_score : previous states are : ";
	    numbers_of_previous_states.Print(cout);
	    cout << "\n";
	    check++;
	}
    }
    
    if (!check) 
    {
	transition_scores.SetElement(j,x);
    }
    return(check);
}

int Hmm_State::set_emission_score(const array<int> &indices, const Score x)
{
    int check = 0;

    if (emission_scores.GetNumberofDimensions()!= number_of_letters_to_read)
    {
	cout << "ERROR class Hmm_State::set_emission_score : array emission_scores does not exist. \n";
	check++;
    }
    if (indices.GetNumberofDimensions()!=1)
    {
	cout << "ERROR class Hmm_State::set_emission_score : indices is no vector\n";
	check++;
    }
    if (indices.GetDimension(0)!=emission_scores.GetNumberofDimensions())
    {
	cout << "ERROR class Hmm_State::set_emission_score : length of vector indices, " << indices.GetDimension(0)
	     << ", does not match number of dimensions of array emission_scores, " 
	     << emission_scores.GetNumberofDimensions() << "\n";
	check++;
    }
    if (!check)
    {
	// check that entries of indices are in correct range
	for (int i=0; i<emission_scores.GetNumberofDimensions(); i++)
	{
	    if ((indices.GetElement(i) < 0) || (indices.GetElement(i) > (emission_scores.GetDimension(i)-1)))
	    {
		cout << "ERROR class Hmm_State::set_emission_score : element " << i << " of indices (" 
		     << indices.GetElement(i) << ") out of range (0..."
		     << emission_scores.GetDimension(i)-1 << ").\n";
		check++;
	    }
	}
    }
    
    if (!check) 
    {
	emission_scores.SetElement (indices,x);
    }
    return(check);
}

int Hmm_State::set_viterbi_strip_score(const int j, const int k, const Score x)
{
  int check = 0;

  if (viterbi_strip.GetNumberofDimensions()!=2)
    {cout << "ERROR class Hmm_State::set_viterbi_strip_score : array viterbi_strip has not 2 "
	  << "dimensions or does not exist\n";
    check++;}
  if ((j<0) || (j>(viterbi_strip.GetDimension(0)-1)))
    {cout << "ERROR class Hmm_State::set_viterbi_strip_score : first index j = " << j 
	  << " out of range (0.." << viterbi_strip.GetDimension(0)-1 << ")\n";
    check++;}
  if ((k<0) || (k>(viterbi_strip.GetDimension(1)-1)))
    {cout << "ERROR class Hmm_State::set_viterbi_strip_score : second index k = " << k 
	  << " out of range (0.." << viterbi_strip.GetDimension(1)-1 << ")\n";
    check++;}

  if (check == 0)
    {
      array<int> indices(1);
      indices.SetDimension(0,2);
      indices.SetElement(0,j);
      indices.SetElement(1,k);
      viterbi_strip.SetElement(indices, x);
    }
  return(check);
}

int Hmm_State::set_viterbi_rectangle_score(const int j, const int k, const Score x)
{
  int check = 0 ;

  if (viterbi_rectangle.GetNumberofDimensions()!=2)
    {cout << "ERROR class Hmm_State::set_viterbi_rectangle_score : array viterbi_rectangle has "
	  << "not 2 dimensions or does not exist\n";
    check++;}
  if ((j<0) || (j>(viterbi_rectangle.GetDimension(0)-1)))
    {cout << "ERROR class Hmm_State::set_viterbi_rectangle_score : first index j = " << j 
	  << " out of range (0.." << viterbi_rectangle.GetDimension(0)-1 << ")\n";
    check++;}
  if ((k<0) || (k>(viterbi_rectangle.GetDimension(1)-1)))
    {cout << "ERROR class Hmm_State::set_viterbi_rectangle_score : second index k = " << k 
	  << " out of range (0.." << viterbi_rectangle.GetDimension(1)-1 << ")\n";
    check++;}

  if (check == 0)
    {
      array<int> indices(1);
      indices.SetDimension(0,2);
      indices.SetElement(0,j);
      indices.SetElement(1,k);
      viterbi_rectangle.SetElement(indices, x);
    }
  return(check);
}

int Hmm_State::set_viterbi_score(const int j, const int k, const Score x)
{
  int check = 0;

  if (viterbi_scores.GetNumberofDimensions()!=2)
    {cout << "ERROR class Hmm_State::set_viterbi_score : array viterbi_scores has not "
	  << "2 dimensions or does not exist\n";
    check++;}
  if ((j<0) || (j>(viterbi_scores.GetDimension(0)-1)))
    {cout << "ERROR class Hmm_State::set_viterbi_score : first index j = " << j 
	  << " out of range (0.." << viterbi_scores.GetDimension(0)-1 << ")\n";
    check++;}
  if ((k<0) || (k>(viterbi_scores.GetDimension(1)-1)))
    {cout << "ERROR class Hmm_State::set_viterbi_score : second index k = " << k 
	  << " out of range (0.." << viterbi_scores.GetDimension(1)-1 << ")\n";
    check++;}

  if (check == 0)
    {
      array<int> indices(1);
      indices.SetDimension(0,2);
      indices.SetElement(0,j);
      indices.SetElement(1,k);
      viterbi_scores.SetElement(indices, x);
    }
  return(check);
}

int Hmm_State::set_dimensions_viterbi_strip(const int j, const int k)
{
  int check = 0;

  if (j<1) 
    {cout << "ERROR class Hmm_State::set_dimensions_viterbi_strip : length of first dimension ("
	  << j << ") invalid.\n";
    check++;}
  if (k<1) 
    {cout << "ERROR class Hmm_State::set_dimensions_viterbi_strip : length of second dimension ("
	  << k << ") invalid.\n";
    check++;}

  if (check == 0)
    {
      viterbi_strip.SetDimension(0,j,Logzero);
      viterbi_strip.SetDimension(1,k,Logzero);
    }
  return(check);
}

int Hmm_State::set_dimensions_viterbi_rectangle(const int j, const int k)
{
  int check = 0;

  if (j<1) 
    {cout << "ERROR class Hmm_State::set_dimensions_viterbi_rectangle: length of first dimension ("
	  << j << ") invalid.\n";
    check++;}
  if (k<1) 
    {cout << "ERROR class Hmm_State::set_dimensions_viterbi_rectangle : length of second dimension ("
	  << k << ") invalid.\n";
    check++;}

  if (check == 0)
    {
      viterbi_rectangle.SetDimension(0,j,Logzero);
      viterbi_rectangle.SetDimension(1,k,Logzero);
    }
  return(check);
}

int Hmm_State::set_dimensions_viterbi_scores(const int j, const int k)
{
  int check = 0;

  if (j<1) 
    {cout << "ERROR class Hmm_State::set_dimensions_viterbi_scores: length of first dimension ("
	  << j << ") invalid.\n";
    check++;}
  if (k<1) 
    {cout << "ERROR class Hmm_State::set_dimensions_viterbi_scores : length of second dimension ("
	  << k << ") invalid.\n";
    check++;}

  if (check == 0)
    {
      viterbi_scores.SetDimension(0,j,Logzero);
      viterbi_scores.SetDimension(1,k,Logzero);
    }
  return(check);
}

int Hmm_State::set_emission_probs_below_threshold_to_threshold(const Prob threshold) 
{
  int check = 0;  
  int i = 0;
  int length = emission_probs.GetDimension(0);
  for (i=1; i<emission_probs.GetNumberofDimensions(); i++)
    {length*=emission_probs.GetDimension(i);}

  if (length == 0) {
    cout << "ERROR: class Hmm_State::set_emission_probs_below_threshold_to_threshold: "
	 << "array of emission_probs does not exist.\n";
    check++;
  }
  if ((threshold <= static_cast<Prob>(0.0)) || (threshold >= static_cast<Prob>(1.0))) {
    cout << "ERROR: class Hmm_State::set_emission_probs_below_threshold_to_threshold: "
	 << "threshold (" << threshold << ") out of range ]0,1[.\n";
    check++;
  }

  if (check == 0) {

#ifdef _PRINT
    array<int> indices;
    const int d   = emission_probs.GetNumberofDimensions(); 
    int*      dim = new int[d];
    for (i=0; i<d; i++) {
      dim[i] = emission_probs.GetDimension(i);
    }
#endif
    
    for (i=0; i<length; i++) {
      if ((emission_probs.GetElement(i) != static_cast<Prob>(0.0)) &&
	  (emission_probs.GetElement(i) <  threshold)) {

#ifdef _PRINT
	indices = get_indices(i, dim, d);
	indices.PrintwithoutNewLine(cout);
	cout << " = " << emission_probs.GetElement(i);
#endif	
	emission_probs.SetElement(i, threshold);
#ifdef _PRINT
	cout << " => " << emission_probs.GetElement(i) << "\n";
#endif
      }
    }
  }
  return(check);
}


int Hmm_State::get_mirrored_copy_of(const Hmm_State &p)
{
    // NOTE : new state has empty arrays for emission and transition scores
    //        and transition probs, only array of emission probs is filled
    //        use Hmm::calculate_scores_from_probs on complete pairhmm to 
    //        derive them contents of empty arrays
    
    int check = 0;
    int i = 0;
    
    if (this == &p) {
	cout << "ERROR class Hmm_State::get_mirrored_copy_of : cannot be applied to state itself.\n"
	     << "Calling deconstructor for this state.\n";
	
	this->~Hmm_State();
    }
    else {
	
	special=p.special;
	special_emission = p.special_emission;
	mirrored=(p.mirrored+1)%2; // make sure that mirrored is 0 or 1
	
	alphabet=p.alphabet;
	number_of_state=p.number_of_states-1-p.number_of_state;
	number_of_child_state_x=p.number_of_child_state_x;
	number_of_child_state_y=p.number_of_child_state_y;
	
	number_of_states=p.number_of_states;
	//set state labels
	number_of_state_labels = p.number_of_state_labels;
	if(state_labels) delete[] state_labels;
	state_labels = NULL;
	state_labels = new array<int>[number_of_state_labels];
	for(i =0; i<number_of_state_labels; i++){
	    state_labels[i] = p.state_labels[i];
	}
	if(state_labels_type) delete[] state_labels_type;
	state_labels_type = NULL;
	state_labels_type = new char*[number_of_state_labels];
	for(i=0;i<number_of_state_labels;i++){
	    state_labels_type[i] = Nullstrcpy(p.state_labels_type[i],check);
	    if(check){
		cout<<"Error: get_mirrored copy of, error in copy state_labels_type["<<i<<"] : "<<state_labels_type[i]<<endl;
		break;
	    }
	}
	
	// Do not have to care about the start and end state, as they both do not emit
	number_of_letters_to_read=p.number_of_letters_to_read;
	number_of_letters_to_read_x = p.number_of_letters_to_read_x;
	number_of_letters_to_read_y = p.number_of_letters_to_read_y;

	// set info on previous and next states 
	// ----------------------------------------------------------------------
	
	number_of_previous_states=p.number_of_next_states;
	number_of_special_transitions_to_previous_states=p.number_of_special_transitions_to_next_states;
      
	numbers_of_previous_states.SetNumberofDimensions(1);
	numbers_of_previous_states.SetDimension(0, number_of_previous_states);
	special_flags_of_transitions_to_previous_states.SetNumberofDimensions(1);
	special_flags_of_transitions_to_previous_states.SetDimension(0, number_of_previous_states);

	if (mirrored == 1) {

	    // offset_for_previous_states_x/y empty
	    
	    offset_for_previous_states_x=p.offset_for_next_states_x;
	    offset_for_previous_states_y=p.offset_for_next_states_y;
	    
	    numbers_of_from_child_states_x_previous = p.numbers_of_from_child_states_x_next;
	    numbers_of_from_child_states_y_previous = p.numbers_of_from_child_states_y_next;
	    numbers_of_to_child_states_x_previous   = p.numbers_of_to_child_states_x_next;  
	    numbers_of_to_child_states_y_previous   = p.numbers_of_to_child_states_y_next;  
	}
	else // if mirrored == 0
	{
	    if ((number_of_letters_to_read_x!=0)&&(number_of_letters_to_read_y!=0))
	    {
		offset_for_previous_states_x.SetNumberofDimensions(1);
		offset_for_previous_states_x.SetDimension(0, number_of_previous_states);
		offset_for_previous_states_y.SetNumberofDimensions(1);
		offset_for_previous_states_y.SetDimension(0, number_of_previous_states);
		
		numbers_of_from_child_states_x_previous.SetNumberofDimensions(1);
		numbers_of_from_child_states_x_previous.SetDimension(0, number_of_previous_states);
		numbers_of_from_child_states_y_previous.SetNumberofDimensions(1);
		numbers_of_from_child_states_y_previous.SetDimension(0, number_of_previous_states);
		numbers_of_to_child_states_x_previous.SetNumberofDimensions(1);
		numbers_of_to_child_states_x_previous.SetDimension(0, number_of_previous_states);
		numbers_of_to_child_states_y_previous.SetNumberofDimensions(1);
		numbers_of_to_child_states_y_previous.SetDimension(0, number_of_previous_states);
	    }
	    else if (number_of_letters_to_read_x!=0)
	    {
		offset_for_previous_states_x.SetNumberofDimensions(1);
		offset_for_previous_states_x.SetDimension(0, number_of_previous_states);
		offset_for_previous_states_y.SetNumberofDimensions(0);
		
		numbers_of_from_child_states_x_previous.SetNumberofDimensions(1);
		numbers_of_from_child_states_x_previous.SetDimension(0, number_of_previous_states);
		numbers_of_from_child_states_y_previous.SetNumberofDimensions(0);
		numbers_of_to_child_states_x_previous.SetNumberofDimensions(1);
		numbers_of_to_child_states_x_previous.SetDimension(0, number_of_previous_states);
		numbers_of_to_child_states_y_previous.SetNumberofDimensions(0);
	    }
	    else if (number_of_letters_to_read_y!=0)
	    {
		offset_for_previous_states_x.SetNumberofDimensions(0);
		offset_for_previous_states_y.SetNumberofDimensions(1);
		offset_for_previous_states_y.SetDimension(0, number_of_previous_states);
		
		numbers_of_from_child_states_x_previous.SetNumberofDimensions(0);
		numbers_of_from_child_states_y_previous.SetNumberofDimensions(1);
		numbers_of_from_child_states_y_previous.SetDimension(0, number_of_previous_states);
		numbers_of_to_child_states_x_previous.SetNumberofDimensions(0);
		numbers_of_to_child_states_y_previous.SetNumberofDimensions(1);
		numbers_of_to_child_states_y_previous.SetDimension(0, number_of_previous_states);
	    }
	}
	
	{
	    for (int i=0; i<number_of_previous_states; i++)
	    {
		numbers_of_previous_states.SetElement(number_of_previous_states-1-i, number_of_states - 1
						      - p.numbers_of_next_states.GetElement(i));
		special_flags_of_transitions_to_previous_states.SetElement(number_of_previous_states-1-i,
									   p.special_flags_of_transitions_to_next_states.GetElement(i));
		if (mirrored == 0)
		{
		    if (number_of_letters_to_read_x!=0)
		    {
			offset_for_previous_states_x.SetElement(number_of_previous_states-1-i,
								p.offset_for_next_states_x.GetElement(i));
			numbers_of_from_child_states_x_previous.
			    SetElement(number_of_previous_states-1-i, 
				       p.numbers_of_from_child_states_x_next.GetElement(i));     
			numbers_of_to_child_states_x_previous.
			    SetElement(number_of_previous_states-1-i, 
				       p.numbers_of_to_child_states_x_next.GetElement(i));     
		    }
		    if (number_of_letters_to_read_y!=0)
		    {
			offset_for_previous_states_y.SetElement(number_of_previous_states-1-i,
								p.offset_for_next_states_y.GetElement(i));
			numbers_of_from_child_states_y_previous.
			    SetElement(number_of_previous_states-1-i, 
				       p.numbers_of_from_child_states_y_next.GetElement(i));     
			numbers_of_to_child_states_y_previous.
			    SetElement(number_of_previous_states-1-i, 
				       p.numbers_of_to_child_states_y_next.GetElement(i));     
		    }
		}
	    }
	}
	number_of_next_states=p.number_of_previous_states;
	number_of_special_transitions_to_next_states=p.number_of_special_transitions_to_previous_states;
	
	numbers_of_next_states.SetNumberofDimensions(1);
	numbers_of_next_states.SetDimension(0, number_of_next_states);
	special_flags_of_transitions_to_next_states.SetNumberofDimensions(1);
	special_flags_of_transitions_to_next_states.SetDimension(0, number_of_next_states);
	
	if (mirrored == 0)
	{
	    // offset_for_next_states_x/y empty
	    
	    offset_for_next_states_x=p.offset_for_previous_states_x;
	    offset_for_next_states_y=p.offset_for_previous_states_y;
	    
	    numbers_of_from_child_states_x_next = p.numbers_of_from_child_states_x_previous;
	    numbers_of_from_child_states_y_next = p.numbers_of_from_child_states_y_previous;
	    numbers_of_to_child_states_x_next   = p.numbers_of_to_child_states_x_previous;  
	    numbers_of_to_child_states_y_next   = p.numbers_of_to_child_states_y_previous;  
	}
	else // if mirrored == 1
	{
	    if ((number_of_letters_to_read_x!=0)&&(number_of_letters_to_read_y!=0))
	    {
		offset_for_next_states_x.SetNumberofDimensions(1);
		offset_for_next_states_x.SetDimension(0, number_of_next_states);
		offset_for_next_states_y.SetNumberofDimensions(1);
		offset_for_next_states_y.SetDimension(0, number_of_next_states);
		
		numbers_of_from_child_states_x_next.SetNumberofDimensions(1);
		numbers_of_from_child_states_x_next.SetDimension(0, number_of_next_states);
		numbers_of_from_child_states_y_next.SetNumberofDimensions(1);
		numbers_of_from_child_states_y_next.SetDimension(0, number_of_next_states);
		numbers_of_to_child_states_x_next.SetNumberofDimensions(1);
		numbers_of_to_child_states_x_next.SetDimension(0, number_of_next_states);
		numbers_of_to_child_states_y_next.SetNumberofDimensions(1);
		numbers_of_to_child_states_y_next.SetDimension(0, number_of_next_states);
	    }
	    else if (number_of_letters_to_read_x!=0)
	    {
		offset_for_next_states_x.SetNumberofDimensions(1);
		offset_for_next_states_x.SetDimension(0, number_of_next_states);
		offset_for_next_states_y.SetNumberofDimensions(0);
		
		numbers_of_from_child_states_x_next.SetNumberofDimensions(1);
		numbers_of_from_child_states_x_next.SetDimension(0, number_of_next_states);
		numbers_of_from_child_states_y_next.SetNumberofDimensions(0);
		numbers_of_to_child_states_x_next.SetNumberofDimensions(1);
		numbers_of_to_child_states_x_next.SetDimension(0, number_of_next_states);
		numbers_of_to_child_states_y_next.SetNumberofDimensions(0);
	    }
	    else if (number_of_letters_to_read_y!=0)
	    {
		offset_for_next_states_x.SetNumberofDimensions(0);
		offset_for_next_states_y.SetNumberofDimensions(1);
		offset_for_next_states_y.SetDimension(0, number_of_next_states);
		
		numbers_of_from_child_states_x_next.SetNumberofDimensions(0);
		numbers_of_from_child_states_y_next.SetNumberofDimensions(1);
		numbers_of_from_child_states_y_next.SetDimension(0, number_of_next_states);
		numbers_of_to_child_states_x_next.SetNumberofDimensions(0);
		numbers_of_to_child_states_y_next.SetNumberofDimensions(1);
		numbers_of_to_child_states_y_next.SetDimension(0, number_of_next_states);
	    }
	}
	
	{
	    for (int i=0; i<number_of_next_states; i++)
	    {
		numbers_of_next_states.SetElement(number_of_next_states-1-i, number_of_states - 1
						  - p.numbers_of_previous_states.GetElement(i));
		special_flags_of_transitions_to_next_states.SetElement(number_of_next_states-1-i,
								       p.special_flags_of_transitions_to_previous_states.GetElement(i));
		
		if (mirrored == 1)
		{
		    if (number_of_letters_to_read_x!=0)
		    {
			offset_for_next_states_x.SetElement(number_of_next_states-1-i,
							    p.offset_for_previous_states_x.GetElement(i));
			
			numbers_of_from_child_states_x_next.
			    SetElement(number_of_next_states-1-i, 
				       p.numbers_of_from_child_states_x_previous.GetElement(i));     
			numbers_of_to_child_states_x_next.
			    SetElement(number_of_next_states-1-i, 
				       p.numbers_of_to_child_states_x_previous.GetElement(i));     
		    }
		    if (number_of_letters_to_read_y!=0)
		    {
			offset_for_next_states_y.SetElement(number_of_next_states-1-i,
							    p.offset_for_previous_states_y.GetElement(i));
			
			numbers_of_from_child_states_y_next.
			    SetElement(number_of_next_states-1-i, 
				       p.numbers_of_from_child_states_y_previous.GetElement(i));     
			numbers_of_to_child_states_y_next.
			    SetElement(number_of_next_states-1-i, 
				       p.numbers_of_to_child_states_y_previous.GetElement(i));     
		    }
		}
	    }
	}
	
	// set number of dimensions of probs (initialised with 0) and scores (initialised with Logzero)
	
	transition_probs.SetNumberofDimensions(1);                             
	transition_probs.SetDimension(0,number_of_states, static_cast<Prob>(0.0));
	
	transition_probs_expressions.SetNumberofDimensions(1);                             
	transition_probs_expressions.SetDimension(0,number_of_states, static_cast<char*>(NULL));
	
	int n_of_dimensions = p.number_of_letters_to_read; 
	
	if(emission_probs_expression) delete[] emission_probs_expression;
	emission_probs_expression = NULL;
	emission_probs_expression = Nullstrcpy(p.emission_probs_expression,check);
	if(check){
	    cout<<"Error: get_mirrored_copy_of, error in emission_probs_expression : "<<p.emission_probs_expression<<endl;
	}
	emission_get_from = p.emission_get_from;
	emission_sum_over = p.emission_sum_over;
	emission_product = p.emission_product;
	num_sum_over = p.num_sum_over;
	sum_over_FromState = p.sum_over_FromState;
	sum_over_FromParam = p.sum_over_FromParam;
	product_FromState = p.product_FromState;
	product_FromParam = p.product_FromParam;
      
	emission_probs.SetNumberofDimensions(n_of_dimensions); 
	{
	    for (int j=0; j<n_of_dimensions; j++)
	    {emission_probs.SetDimension(j,alphabet, static_cast<Prob>(0.0));}
	}
	
	if(strcmp(emission_probs_expression,"none"))
	{
	    int max=0;
	    int d=0;
	  
	    d= this->get_letters_to_read();
	    if (d > 0) {
		max= static_cast<int>(pow( static_cast<float>(alphabet), static_cast<float>(d)));
	    }
	    else {
		max=0;
	    }

	    if (d>0) 
	    { // if this is not a silent state
		
		array<int> index(1);
		index.SetDimension(0, d);
		
		Prob emission_prob=0;
	      
		for (int number=0; number<max; number++) { // loop over all possible indices
		  
		    // transform linear index into index 
		  
		    emission_prob=0;
		    
		    convert_to_base(number, alphabet, &index);		  
		    emission_prob=p.get_emission_prob(index);
			
		    index.ResetData();
		    convert_to_base_and_invert(number, alphabet, &index);
		    this->set_emission_prob(index, emission_prob);
		
		    index.ResetData();
		}	     
	    }
	}
	
	transition_scores.SetNumberofDimensions(1);                            
	transition_scores.SetDimension(0, number_of_states, Logzero);  
	
	emission_scores.SetNumberofDimensions(n_of_dimensions);      
	{
	    for (int j=0; j<n_of_dimensions; j++)
	    {emission_scores.SetDimension(j, alphabet, Logzero);}
	}
      
	// the elements of the following arrays which hold result of a particular run 
	// of an algorithm will not be copied
	
	viterbi_strip.SetNumberofDimensions(2);    
	viterbi_rectangle.SetNumberofDimensions(2);    
	viterbi_scores.SetNumberofDimensions(2);    
      
	this->order_next_and_previous_states();
    }
    return(check);
}


int Hmm_State::get_complemented_copy_of(const Hmm_State &p)

    // NOTE : this function does not mirror the input state
    //        new state has empty arrays for emission and transition scores
    //        and transition probs, only array of emission probs is filled
    //        use Hmm::calculate_scores_from_probs on complete pairhmm to 
    //        derive them contents of empty arrays
{
    int check = 0;
    
    if (this == &p) {
	cout << "ERROR class Hmm_State::get_complemented_copy_of : cannot be applied to state itself.\n"
	     << "Calling deconstructor for this state.\n";
	
	this->~Hmm_State();
    }
    if (check == 0) {
	
	// emission indices: (x_3,x_2,x_1,y_3,y_2,y_1) => (x_3',x_2',x_1',y_3',y_2',y_1')
	//
	// where x' = complement_alphabet(x) i.e. for DNA : A <-> T etc.
	
	// first get mirrored copy of input state p
	// NOTE: this state then has empty array transition_scores and emission_scores

	*this = p;
	
	int i   = 0;
	int j   = 0;
	
	int d   = this->emission_probs.GetNumberofDimensions();
	int max = 0;
	
	if (d > 0) {
	    
	    max= static_cast<int>(pow( static_cast<float>(alphabet), static_cast<float>(d)));
	    
	    
	    array<int> old_index(1);
	    old_index.SetDimension(0, d);
	    array<int> new_index(1);
	    new_index.SetDimension(0, d);
	    
	    Prob prob = 0;
	    
	    for (j=0; j<max; j++) { // loop over all possible indices
		
		// order of x-indices reverses, order of y-indices reverses
		
		prob=0;
		old_index.ResetData();
		convert_to_base(j, alphabet, &old_index);	  
		prob= p.get_emission_prob(old_index);
       		
		new_index.ResetData();		
		new_index.SetDimension(0, d);
	
	    }	
	    
	    this->set_emission_prob(new_index, prob);
		
	}
    }
    return(check);
}

int Hmm_State::get_alphabet(void) const
{
    return(alphabet);
}

int Hmm_State::get_special(void) const
{
    return(special);
}

int Hmm_State::get_special_emission(void) const
{
    return(special_emission);
}

int Hmm_State::get_mirrored(void) const
{
    return(mirrored);
}

int Hmm_State::get_number_of_state_labels(void) const
{
    return (number_of_state_labels);
}

char* Hmm_State::get_state_labels_type(const int i) const
{
    if ((i<0) || (i > number_of_state_labels))
    {
      cout << "ERROR class Hmm_State::get_state_labels_type: i (" << i 
	   << ") out of range (0.." << number_of_state_labels << ").\n";
      return NULL;
    }
    return state_labels_type[i];
}

int Hmm_State:: get_state_labels_dim(const int i) const
{
    if ((i<0) || (i > number_of_state_labels))
    {
	cout << "ERROR class Hmm_State::get_state_labels_dim: i (" << i 
	     << ") out of range (0.." << number_of_state_labels << ").\n";
	return -1 ;
    }
    return state_labels[i].GetDimension(0);
}

int Hmm_State:: get_state_labels(const int i, const int j) const
{
    if ((i<0) || (i > number_of_state_labels))
    {
	cout << "ERROR class Hmm_State::get_state_labels: i (" << i 
	     << ") out of range (0.." << number_of_state_labels -1<< ").\n";
	return -1 ;
    }
  
    if((j<0) || (j > state_labels[i].GetDimension(0)))
    {
	cout<<"ERROR class Hmm_State::get_state_labels: j(" << j
	  << ") out of range(0.." <<state_labels[i].GetDimension(0) <<").\n";
	return -1;
    }
    return state_labels[i].GetElement(j);
    
}

int Hmm_State::get_number_of_state(void) const
{
    return(number_of_state);
}



int Hmm_State::get_number_of_child_state_x(void) const
{
    return(number_of_child_state_x);
}

int Hmm_State::get_number_of_child_state_y(void) const
{
    return(number_of_child_state_y);
}

int Hmm_State::get_number_of_states(void) const
{
    return(number_of_states);
}

int Hmm_State::get_number_of_previous_states(void) const
{
    return(number_of_previous_states);
}

int Hmm_State::get_number_of_special_transitions_to_previous_states(void) const
{
    return(number_of_special_transitions_to_previous_states);
}

int Hmm_State::get_number_of_previous_state(const int i) const
{

  if ((i<0) || (i> (number_of_previous_states-1)))
  {
      cout << "ERROR class Hmm_State::get_number_of_previous_state : index i ("
	   << i << ") out of range (0.." << number_of_previous_states-1 << ").\n";

      return(0);
  }
  else
  {
      return(numbers_of_previous_states.GetElement(i));
  }
  
}

int Hmm_State::get_index_of_previous_state(const int number_of_previous_state) const
{

    if ((number_of_previous_state < 0) || (number_of_previous_state > (number_of_states-1)))
    {
	cout << "ERROR class Hmm_State::get_index_of_previous_state : number_of_previous_state ("
	     << number_of_previous_state << ") out of range (0.." << number_of_states-1 << ").\n";
	
	return(0);
    }
    else
    {

	int i=0;
	int index=-1;
	
	for (i=0; i<number_of_previous_states; i++)
	{
	    if (numbers_of_previous_states.GetElement(i) == number_of_previous_state)
	    {
		index = i;
		break;
	    }
	}
	return(index);
	
    }

}

int Hmm_State::get_special_flag_of_transition_to_previous_state(const int i) const
{
    if ((i<0) || (i> (number_of_previous_states-1)))
    {
	cout << "ERROR class Hmm_State::get_special_flag_of_transition_to_previous_state : index i ("
	     << i << ") out of range (0.." << number_of_previous_states-1 << ").\n";
	return(0);
    }
    else
    {
	return(special_flags_of_transitions_to_previous_states.GetElement(i));
    }
}

int Hmm_State::get_offset_for_previous_state_x(const int i) const
{
    int check = 0;
    
    if ((i<0) || (i> (number_of_previous_states-1)))
    {
	cout << "ERROR class Hmm_State::get_offset_for_previous_state_x : index i ("
	     << i << ") out of range (0.." << number_of_previous_states-1 << ").\n";
	check++;
    }
    if (mirrored == 1)
    {
	cout << "ERROR class Hmm_State::get_offset_for_previous_state_x : array offset_for_previous_states_x \n"
	     << "is empty as this state is mirrored.\n";
	check++;
    }
    
    if (check == 0)
    {
	if (number_of_letters_to_read_x!=0)
	{
	    return(offset_for_previous_states_x.GetElement(i));
	}
	else
	{
	    return(0);
	}
    }
    else
    {
	return(0);
    }
}

int Hmm_State::get_offset_for_previous_state_y(const int i) const
{
    int check = 0;
    
    if ((i<0) || (i> (number_of_previous_states-1)))
    {
	cout << "ERROR class Hmm_State::get_offset_for_previous_state_y : index i ("
	     << i << ") out of range (0.." << number_of_previous_states-1 << ").\n";
	check++;
    }
    if (mirrored == 1)
    {
	cout << "ERROR class Hmm_State::get_offset_for_previous_state_y : array offset_for_previous_states_y \n"
	     << "is empty as this state is mirrored.\n";
	check++;
    }
    
    if (check == 0)
    {
	if (number_of_letters_to_read_y!=0)
	{
	    return(offset_for_previous_states_y.GetElement(i));
	}
	else
	{
	    return(0);
	}
    }
    else
    {
	return(0);
    }
}

int Hmm_State::get_number_of_next_states(void) const
{
    return(number_of_next_states);
}

int Hmm_State::get_number_of_special_transitions_to_next_states(void) const
{
    return(number_of_special_transitions_to_next_states);
}

int Hmm_State::get_number_of_next_state(const int i) const
{
    if ((i<0) || (i> (number_of_next_states-1)))
    {
      cout << "ERROR class Hmm_State::get_number_of_next_state : index i ("
	   << i << ") out of range (0.." << number_of_next_states-1 << ").\n";
      
      return(0);
    }
    else
    {
	return(numbers_of_next_states.GetElement(i));

    }
    
}

int Hmm_State::get_index_of_next_state(const int number_of_next_state) const
{
    
    if ((number_of_next_state < 0) || (number_of_next_state > (number_of_states-1)))
    {
	cout << "ERROR class Hmm_State::get_index_of_next_state : number_of_next_state ("
	     << number_of_next_state << ") out of range (0.." << number_of_states-1 << ").\n";

	return(-1);
    }
    else
    {

	int i=0;
	int index=-1;

	for (i=0; i<number_of_next_states; i++)
	{
	    if (numbers_of_next_states.GetElement(i) == number_of_next_state)
	    {
		index = i;
		break;
	    }
	}
	return(index);
    }
    
}

int Hmm_State::get_special_flag_of_transition_to_next_state(const int i) const
{
    if ((i<0) || (i> (number_of_next_states-1)))
    {
	cout << "ERROR class Hmm_State::get_special_flag_of_transition_to_next_state : index i ("
	     << i << ") out of range (0.." << number_of_next_states-1 << ").\n";
	return(0);
    }
    else
    {
	return(special_flags_of_transitions_to_next_states.GetElement(i));
    }
}

int Hmm_State::get_offset_for_next_state_x(const int i) const
{
    int check = 0;
  
    if ((i<0) || (i> (number_of_next_states-1)))
    {
	cout << "ERROR class Hmm_State::get_offset_for_next_state_x : index i ("
	     << i << ") out of range (0.." << number_of_next_states-1 << ").\n";
	check++;
    }
    if (mirrored == 0)
    {
	cout << "ERROR class Hmm_State::get_offset_for_next_state_x : array offset_for_next_states_x \n"
	     << "is empty as this state is un-mirrored.\n";
	check++;
    }
    
    if (check == 0)
    {
	if (number_of_letters_to_read_x!=0)
	{
	    return(offset_for_next_states_x.GetElement(i));
	}
	else
	{
	    return(0);
	}
    }
    else
    {
	return(0);
    }
}

int Hmm_State::get_offset_for_next_state_y(const int i) const
{
    int check = 0;
    
    if ((i<0) || (i> (number_of_next_states-1)))
    {
	cout << "ERROR class Hmm_State::get_offset_for_next_state_y : index i ("
	     << i << ") out of range (0.." << number_of_next_states-1 << ").\n";
	check++;
    }
    if (mirrored == 0)
    {
	cout << "ERROR class Hmm_State::get_offset_for_next_state_y : array offset_for_next_states_y \n"
	     << "is empty as this state is un-mirrored.\n";
	check++;
    }
    
    if (check == 0)
    {
	if (number_of_letters_to_read_y!=0)
	{
	    return(offset_for_next_states_y.GetElement(i));
	}
	else
	{
	    return(0);
	}
    }
    else
    {
	return(0);
    }
}

int Hmm_State::get_letters_to_read(void) const
{
  return(number_of_letters_to_read);
}

int Hmm_State::get_letters_to_read_x(void) const
{
    return (number_of_letters_to_read_x);
}

int Hmm_State::get_letters_to_read_y(void) const
{
    return (number_of_letters_to_read_y);
}

int Hmm_State::is_state_previous_state(const int i) const
{
    int found_previous_state = 0;
    
    if ((i<0) || (i>(number_of_states-1))) 
    {
	cout << "ERROR class Hmm_State::is_state_previous_state: \n"
	     << "state i (" << i << ") out of range (0..." << number_of_states-1 << ").\n";
    }
    else
    {
	
	for (int j=0; j<number_of_previous_states; j++)
	{
	    if (numbers_of_previous_states.GetElement(j) == i) 
	    {
	      found_previous_state++;
	      break;
	    }
	}
	
    }
    
    return(found_previous_state);
}

int Hmm_State::is_transition_to_previous_state_special(const int i) const
{
    int transition_to_previous_state_i_is_special = 0;
    
    if ((i<0) || (i>(number_of_states-1))) 
    {
	cout << "ERROR class Hmm_State::is_transition_to_previous_state_special: \n"
	     << "state i (" << i << ") out of range (0..." << number_of_states-1 << ").\n";
    }
    else
    {
	if (number_of_special_transitions_to_previous_states > 0)
	{
	    int j;
	    for (j=0; j<number_of_previous_states; j++)
	    {
		if ((special_flags_of_transitions_to_previous_states.GetElement(j) == 1) &&
		    (numbers_of_previous_states.GetElement(j) == i))
		{
		    transition_to_previous_state_i_is_special++;
		    break;
		}
	    }
	}
    }
    return(transition_to_previous_state_i_is_special);
}

int Hmm_State::get_number_of_from_child_state_x_previous(const int number_of_previous_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_previous_state<0) || (number_of_previous_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_previous: \n"
	     << "state i (" << number_of_previous_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 1) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_previous: \n"
	     << "this state is mirrored and has no previous child states.\n";
	check++;
    }
    if (number_of_letters_to_read_x == 0){
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_previous: \n"
	     << "this state is not of type EmitX or EmitXY, which has no previous child states x.\n";
	check++;
    }
    
    if (check == 0) {
	const int index = get_index_of_previous_state(number_of_previous_state);

	if (index > -1) {
	    return_number = numbers_of_from_child_states_x_previous.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_previous: \n"
		 << "input state with number (" << number_of_previous_state 
		 << ") is no previous state. Returning value -1.\n";
	}
    }
    return(return_number);
}

int Hmm_State::get_number_of_from_child_state_y_previous(const int number_of_previous_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_previous_state<0) || (number_of_previous_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_previous: \n"
	     << "state i (" << number_of_previous_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 1) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_previous: \n"
	     << "this state is mirrored and has no previous child states.\n";
	check++;
    }
    if (number_of_letters_to_read_y==0){
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_previous: \n"
	     << "this state is not of type EmitY or EmitXY, which has no previous child states y.\n";
	check++;
    }
    
    if (check == 0) {

	const int index = get_index_of_previous_state(number_of_previous_state);
	
	if (index > -1) {
	    return_number = numbers_of_from_child_states_y_previous.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_previous: \n"
		 << "input state with number (" << number_of_previous_state 
		 << ") is no previous state. Returning value -1.\n";
	}
  
    }
 
    return(return_number);
 }

int Hmm_State::get_number_of_from_child_state_x_next(const int number_of_next_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_next_state<0) || (number_of_next_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_next: \n"
	     << "state i (" << number_of_next_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 0) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_next: \n"
	     << "this state is unmirrored and has no next child states.\n";
	check++;
    }
    if (number_of_letters_to_read_x==0) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_next: \n"
	     << "this state is not of type EmitX or EmitXY, which has no next child states x.\n";
	check++;
    }
    
    if (check == 0) {
	const int index = get_index_of_next_state(number_of_next_state);
	
	if (index > -1) {
	    return_number = numbers_of_from_child_states_x_next.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_from_child_state_x_next: \n"
		 << "input state with number (" << number_of_next_state 
		 << ") is no next state. Returning value -1.\n";
	}
    }
    return(return_number);
}

int Hmm_State::get_number_of_from_child_state_y_next(const int number_of_next_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_next_state<0) || (number_of_next_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_next: \n"
	 << "state i (" << number_of_next_state << ") out of range (0..." 
	 << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 0) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_next: \n"
	     << "this state is unmirrored and has no next child states.\n";
	check++;
    }
    if (number_of_letters_to_read_y==0) {
	cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_next: \n"
	     << "this state is not of type EmitY or EmitXY, but of type which has no next child states y.\n";
	check++;
    }
    
    if (check == 0) {

	const int index = get_index_of_next_state(number_of_next_state);
	
	if (index > -1) {
	    return_number = numbers_of_from_child_states_y_next.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_from_child_state_y_next: \n"
		 << "input state with number (" << number_of_next_state 
		 << ") is no next state. Returning value -1.\n";
	}
     
    }
    
    return(return_number);
}

int Hmm_State::get_number_of_to_child_state_x_previous(const int number_of_previous_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_previous_state<0) || (number_of_previous_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_previous: \n"
	     << "state i (" << number_of_previous_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 1) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_previous: \n"
	     << "this state is mirrored and has no previous child states.\n";
	check++;
    }
    if (number_of_letters_to_read_x==0) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_previous: \n"
	     << "this state is not of type EmitX or EmitXY, but of type which has no previous child states x.\n";
	check++;
    }
    
    if (check == 0) {

	const int index = get_index_of_previous_state(number_of_previous_state);
	
	if (index > -1) {
	    return_number = numbers_of_to_child_states_x_previous.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_previous: \n"
		 << "input state with number (" << number_of_previous_state 
		 << ") is no previous state. Returning value -1.\n";
	}
    }
    return(return_number);
 }

int Hmm_State::get_number_of_to_child_state_y_previous(const int number_of_previous_state) const
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_previous_state<0) || (number_of_previous_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_previous: \n"
	     << "state i (" << number_of_previous_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 1) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_previous: \n"
	     << "this state is mirrored and has no previous child states.\n";
	check++;
    }
    if (number_of_letters_to_read_y==0) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_previous: \n"
	     << "this state is not of type EmitY or EmitXY, but of type which has no previous child states y.\n";
	check++;
    }
    
    if (check == 0) {
	
	const int index = get_index_of_previous_state(number_of_previous_state);
	
	if (index > -1) {
	    return_number = numbers_of_to_child_states_y_previous.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_previous: \n"
		 << "input state with number (" << number_of_previous_state 
		 << ") is no previous state. Returning value -1.\n";
	}
       
    }
    
    return(return_number);
}

int Hmm_State::get_number_of_to_child_state_x_next(const int number_of_next_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_next_state<0) || (number_of_next_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_next: \n"
	     << "state i (" << number_of_next_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 0) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_next: \n"
	     << "this state is unmirrored and has no next child states.\n";
	check++;
    }
    if (number_of_letters_to_read_x==0) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_next: \n"
	     << "this state is not of type EmitX or EmitXY, but of type which has no next child states.\n";
	check++;
    }
    
    if (check == 0) {
	
	const int index = get_index_of_next_state(number_of_next_state);
	
	if (index > -1) {
	    return_number = numbers_of_to_child_states_x_next.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_to_child_state_x_next: \n"
		 << "input state with number (" << number_of_next_state 
		 << ") is no next state. Returning value -1.\n";
	}
    }

    return(return_number);
}

int Hmm_State::get_number_of_to_child_state_y_next(const int number_of_next_state) const 
{
    int check  = 0;
    
    int return_number = -1;
    
    if ((number_of_next_state<0) || (number_of_next_state>(number_of_states-1))) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_next: \n"
	     << "state i (" << number_of_next_state << ") out of range (0..." 
	     << number_of_states-1 << ").\n";
	check++;
    }
    if (mirrored == 0) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_next: \n"
	 << "this state is unmirrored and has no next child states.\n";
	check++;
    }
    if (number_of_letters_to_read_y==0) {
	cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_next: \n"
	     << "this state is not of type EmitX or EmitXY, but of type which has no next child states.\n";
	check++;
    }
    
    if (check == 0) {
	const int index = get_index_of_next_state(number_of_next_state);
	
	if (index > -1) {
	    return_number = numbers_of_to_child_states_y_next.GetElement(index);
	}
	else {
	    cout << "ERROR class Hmm_State::get_number_of_to_child_state_y_next: \n"
		 << "input state with number (" << number_of_next_state 
		 << ") is no next state. Returning value -1.\n";
	}
    }
    
    return(return_number);
}

int Hmm_State::is_state_next_state(const int i) const
{
    int found_next_state = 0;
    

    if ((i<0) || (i>(number_of_states-1))) 
    {
	cout << "ERROR class Hmm_State::is_state_next_state: \n"
	     << "state i (" << i << ") out of range (0..." << number_of_states-1 << ").\n";
    }
    else
    {
      
	for (int j=0; j<number_of_next_states; j++)
	{
	    if (numbers_of_next_states.GetElement(j) == i) 
	    {
		found_next_state++;
		break;
	    }
	}
    }
    
    return(found_next_state);
}

int Hmm_State::is_transition_to_next_state_special(const int i) const
{
    int transition_to_next_state_i_is_special = 0;
    
    if ((i<0) || (i>(number_of_states-1))) 
    {
	cout << "ERROR class Hmm_State::is_transition_to_next_state_special: \n"
	   << "state i (" << i << ") out of range (0..." << number_of_states-1 << ").\n";
    }
    else
    {
	if (number_of_special_transitions_to_next_states > 0)
	{
	    int j;
	    for (j=0; j<number_of_next_states; j++)
	    {
		if ((special_flags_of_transitions_to_next_states.GetElement(j) == 1) &&
		    (numbers_of_next_states.GetElement(j) == i))
		{
		    transition_to_next_state_i_is_special++;
		    break;
		}
	    }
	}
    }
    return(transition_to_next_state_i_is_special);
}

int Hmm_State::is_state_special(void) const
{
    return(special);
}

int Hmm_State::has_state_special_emissions(void) const
{
    return(special_emission);
}

int Hmm_State::is_state_mirrored(void) const
{
    return(mirrored);
}

int Hmm_State::exponentiate_matrix_of_emission_probs(const int exponent)
{
  // note : - integer part of exponent is used to exponentiate matrix (i.e. 2.5 => 2)
  //        - exponent must be >= 1

  int check=0;

  if (exponent < 1)
    {
      cout << "ERROR class Hmm_State:: exponentiate_matrix_of_emission_probs : exponent (" << exponent
	   << ") < 1.\n";
      check++;
    }
  if (emission_probs.GetNumberofDimensions()%2 != 0)
    {
      cout << "ERROR class Hmm_State:: exponentiate_matrix_of_emission_probs : matrix of emission_probs "
	   << "has not even number of dimensions, but " << emission_probs.GetNumberofDimensions() 
	   << " dimensions.\n";
      check++;
    }
  
  if (check == 0)
    {
      if (exponent > 1) // for exponent == 1 there is nothing to do
	{
	  array<Prob> old_matrix = emission_probs;
	  array<Prob> new_matrix = emission_probs;
	  
	  for (int i=1; i<exponent; i++)
	    {new_matrix = new_matrix.SpecialMultiply(old_matrix);}
	  
	  emission_probs = new_matrix;
	}
#ifndef _EXPERIMENT 
      check += this->check_matrix_of_emission_probs();
#endif 
    }
  return(check);
}

int Hmm_State::check_matrix_of_emission_probs() const
{
    int check=0;

    if (number_of_letters_to_read>0)
    {
	const int number_of_dimensions = emission_probs.GetNumberofDimensions();
	const int l   = alphabet;
	const int max = static_cast<int>(pow( static_cast<float>(l), static_cast<float>(number_of_dimensions)));
	
	// test if sum of all elements in matrix is 1 and that all elements are >= 0
	
	Prob element             = 0.0;
	Prob sum_of_all_elements = 0.0;
	
	for (int l_index=0; l_index<max; l_index++) // loop over all possible indices
	{
	    element = 0.0;
	    element = emission_probs.GetElement(l_index);
	    
	    if (element < 0) 
	    {
		cout << "ERROR: class Hmm_State:: check_matrix_of_emission_probs : found element < 0\n";
		check++;
		break;
	    }
	    else
	    {
		sum_of_all_elements += element;
	    }
	}
	
	if (abs(sum_of_all_elements - 1.) > Max_deviation)
	{
	    cout << "ERROR: class Hmm_State:: check_matrix_of_emission_probs : matrix is not normalised "
		 << "i.e. the sum of all its elements is not 1.\n";
	    cout << "ERROR: class Hmm_State:: check_matrix_of_emission_probs : |sum_of_all_elements - 1| = "
		 << abs(sum_of_all_elements - 1.) << " > Max_deviation = " << Max_deviation << "\n";
	    check++;
	}
    }
    return(check);
}

char* Hmm_State::get_transition_probs_expressions(const int j) const
{

    int check = 0;
    if ((j<0) || (j>number_of_states-1))
    {cout << "ERROR class Hmm_State::get_transition_probs_expressions : index j (" << j << ") out of range (0.." 
	  << number_of_states-1 << ")\n";
    check++;}
    if (transition_probs.GetNumberofDimensions()!=1 ||
	transition_probs.GetDimension(0)!=number_of_states)
    {cout << "ERROR class Hmm_State::get_transition_probs_expressions : array transition_probs does not exist. \n";
    check++;}

    return (transition_probs_expressions.GetElement(j));
}

Prob Hmm_State::get_transition_prob(const int j) const
{

    int check = 0;
    if ((j<0) || (j>number_of_states-1))
    {cout << "ERROR class Hmm_State::get_transition_prob : index j (" << j << ") out of range (0.." 
	  << number_of_states-1 << ")\n";
    check++;}
    if (transition_probs.GetNumberofDimensions()!=1 ||
	transition_probs.GetDimension(0)!=number_of_states)
    {cout << "ERROR class Hmm_State::get_transition_prob : array transition_probs does not exist. \n";
    check++;}

    return(transition_probs.GetElement(j));
}

#ifdef _PROB
Prob Hmm_State::get_special_transition_prob(const Hmm_State* const to,
						const Sequence*      const x,
						const int                  x_position,
						const Sequence*      const y,
						const int                  y_position) const
{

    Prob return_prob = static_cast<Prob>(0.0);

    int  check = 0;
    
    if (to == NULL)
    {
	cout << "ERROR class Hmm_State::get_special_transition_prob : Hmm_State to is NULL.\n";
	check++;
    }
    else
    {
	if ((to->get_number_of_state() < 0) || (to->get_number_of_state() > (number_of_states-1)))
	{
	    cout << "ERROR class Hmm_State::get_special_transition_prob : number of to state ("
		 << to->get_number_of_state() << ") out of range (0.." << number_of_states-1 << ").\n";
	    check++;
	}
	if (to->get_mirrored() != this->get_mirrored())
	{
	    cout << "ERROR class Hmm_State::get_special_transition_prob : Hmm_State to has mirrored flag ("
		 << to->get_mirrored() << ") != mirrored flag of this state (" << this->get_mirrored() << ").\n";
	    check++;
	}
    }

    if (check == 0)
    {
	int   number_of_to_state     = to->get_number_of_state();
	int   number_of_from_state   = number_of_state;
	int   x_offset_mirrored      = 0;
	int   y_offset_mirrored      = 0;
	int   x_offset_un_mirrored   = 0;
	int   y_offset_un_mirrored   = 0;
	
	if (mirrored == 1)
	{
	    int index = this->get_index_of_next_state(number_of_to_state);
	    x_offset_mirrored      = - this->get_offset_for_next_state_x(index) -1;
	    y_offset_mirrored      = - this->get_offset_for_next_state_y(index) -1;
	}
	if (to->get_mirrored() == 0)
	{
	    int index = to->get_index_of_previous_state(number_of_from_state);
	    x_offset_un_mirrored   = to->get_offset_for_previous_state_x(index);
	    y_offset_un_mirrored   = to->get_offset_for_previous_state_y(index);
	}
      
	int   x_position_mirrored    = x_position + x_offset_mirrored;
	int   y_position_mirrored    = y_position + y_offset_mirrored;
	int   x_position_un_mirrored = x_position + x_offset_un_mirrored;
	int   y_position_un_mirrored = y_position + y_offset_un_mirrored;
	
	if (x == NULL)
	{
	    cout << "ERROR class Hmm_State::get_special_transition_prob : Sequence x is NULL.\n";
	    check++;
	}
	if (y == NULL)
	{
	    cout << "ERROR class Hmm_State::get_special_transition_prob : Sequence y is NULL.\n";
	    check++;
	}
	
	if (check == 0)
	{
	    int n_of_letters_to_read_x_to = to->get_letters_to_read_x();
	    int n_of_letters_to_read_y_to = to->get_letters_to_read_y();

	    int n_of_letters_to_read_x_from = number_of_letters_to_read_x;
	    int n_of_letters_to_read_y_from = number_of_letters_to_read_y;
	    
	    // if transition from -> to is special
	    
	    if (this->is_transition_to_next_state_special(number_of_to_state) == 1)
	    {
		return_prob = this->get_transition_prob(number_of_to_state);
		
		// if there is no sequence to read the special transition prob is set to 0
		
		if (to->get_mirrored() == 0)
		{
		    if ((n_of_letters_to_read_x_to!=0) &&
			((x_position_un_mirrored < 0) || (x_position_un_mirrored > (x->length()-1))))
		    {
			cout << "ERROR class Hmm_State::get_special_transition_prob : x_position ("
			     << x_position << ") + x_offset_un_mirrored ("
			     << x_offset_un_mirrored << ") out of range (0.." << x->length()-1 << ").\n";
			check++;
		    }
		    if ((n_of_letters_to_read_y_to!=0) &&
			((y_position_un_mirrored < 0) || (y_position_un_mirrored > (y->length()-1))))
		    {
			cout << "ERROR class Hmm_State::get_special_transition_prob : y_position ("
			     << y_position << ") + y_offset_un_mirrored ("
			     << y_offset_un_mirrored << ") out of range (0.." << y->length()-1 << ").\n";
			check++;
		    }
		}
		if (to->get_mirrored() == 1) 
		{
		    if ((n_of_letters_to_read_x_from!=0) &&
			((x_position_mirrored < 0) || (x_position_mirrored > (x->length()-1))))
		    {
			cout << "ERROR class Hmm_State::get_special_transition_prob : x_position ("
			     << x_position << ") + x_offset_mirrored ("
			     << x_offset_mirrored << ") out of range (0.." << x->length()-1 << ").\n";
			check++;
		    }
		    if ((n_of_letters_to_read_y_from!=0) &&
			((y_position_mirrored < 0) || (y_position_mirrored > (y->length()-1))))
		    {
			cout << "ERROR class Hmm_State::get_special_transition_prob : y_position ("
			     << y_position << ") + y_offset_mirrored ("
			     << y_offset_mirrored << ") out of range (0.." << y->length()-1 << ").\n";
			check++;
		    }
		}
		
		if (check == 0) // if x_position and y_position are in the correct range
		{
		    if (mirrored == 0)
		    {
			if (n_of_letters_to_read_x_to!=0)
			{
			    return_prob *= x->get_posterior_prob(number_of_from_state,
								 number_of_to_state,
								 x_position + x_offset_un_mirrored);
			}
			if (n_of_letters_to_read_y_to!=0)
			{
			    return_prob *= y->get_posterior_prob(number_of_from_state,
								 number_of_to_state,
								 y_position + y_offset_un_mirrored);
			}
		    }
		    else // if from and to state are mirrored
		    {
			if (n_of_letters_to_read_x_from!=0)
			{
			    return_prob *= x->get_posterior_prob(number_of_states-1-number_of_to_state,
								 number_of_states-1-number_of_from_state,
								 x_position + x_offset_mirrored);
			}
			if (n_of_letters_to_read_y_from!=0)
			{
			    return_prob *= y->get_posterior_prob(number_of_states-1-number_of_to_state,
								 number_of_states-1-number_of_from_state,
								 y_position + y_offset_mirrored);
			    
			}
		    }
		} // if check == 0 i.e. if x_position and y_position are in the correct range
		else
		{
		    return_prob = static_cast<Prob>(0.0);
		}
	    } // if transition is special
	} // if check == 0
    } // if check == 0
    return(return_prob);
}

#endif

Prob Hmm_State::get_emission_prob(const array<int> &indices) const
{

  int check = 0;
  if (emission_probs.GetNumberofDimensions()!=number_of_letters_to_read)
  {cout << "ERROR class Hmm_State::get_emission_prob : array emission_probs does not exist. \n";
  check++;}
  if (indices.GetNumberofDimensions()!=1)
  {cout << "ERROR class Hmm_State::get_emission_prob : indices is no vector\n";
  check++;}
  if (indices.GetDimension(0)!=emission_probs.GetNumberofDimensions())
  {cout << "ERROR class Hmm_State::get_emission_prob : length of vector indices, " << indices.GetDimension(0)
	<< ", does not match number of dimensions of array emission_probs, " 
	<< emission_probs.GetNumberofDimensions() << "\n";
  check++;}
  if (check == 0)
  {
      // check that entries of indices are in correct range
      for (int i=0; i<emission_probs.GetNumberofDimensions(); i++)
      {
	  if ((indices.GetElement(i) < 0) || (indices.GetElement(i) > (emission_probs.GetDimension(i)-1)))
	  {
	      cout << "ERROR class Hmm_State::get_emission_prob : element " << i << " of indices (" 
		   << indices.GetElement(i) << ") out of range (0..."
		   << emission_probs.GetDimension(i)-1 << ").\n";
	      check++;
	  }
      }
  }
  return(emission_probs.GetElement(indices));
}

Prob Hmm_State::get_emission_prob(const int linear_index) const
{

    if (emission_probs.GetNumberofDimensions()!= number_of_letters_to_read)
    {cout << "ERROR class Hmm_State::get_emission_prob : array emission_probs does not exist. \n";}
    int length=emission_probs.GetDimension(0);
    for (int i=1; i<emission_probs.GetNumberofDimensions(); i++)
    {
	length*=emission_probs.GetDimension(i);
    }
    if ((linear_index<0) || (linear_index>(length-1)))
    {cout << "ERROR class Hmm_State::get_emission_prob : linear_index " << linear_index
          << " out of range (0.." << length-1 << ").\n";}
    
    return(emission_probs.GetElement(linear_index));
}

array<Prob> Hmm_State:: get_emission_prob() const
{
	return emission_probs;
}

Prob* Hmm_State::get_emission_prob_array() const
{
    Prob* return_array = NULL;
  
    if (emission_probs.GetNumberofDimensions()== number_of_letters_to_read) {
	
	int i = 0;
	int length=emission_probs.GetDimension(0);
	for (i=1; i<emission_probs.GetNumberofDimensions(); i++) {
	    length*=emission_probs.GetDimension(i);
	}
	return_array = new Prob[length];
	
	for (i=0; i<length; i++) {
	    return_array[i] = emission_probs.GetElement(i);
	}
    }
    return(return_array);
}

char* Hmm_State::get_emission_probs_expression() const{
    return emission_probs_expression;
}

int Hmm_State::get_emission_get_from() const{
    return emission_get_from;
}

bool Hmm_State::get_emission_sum_over() const{
    return emission_sum_over;
}
 
bool Hmm_State::get_emission_product() const{
    return emission_product;
}

int Hmm_State::get_num_sum_over() const{
    return num_sum_over;
}

array<int> Hmm_State::get_sum_over_ThisPos() const{
    return sum_over_ThisPos;
}

int Hmm_State::get_sum_over_ThisPos(int i) const{
    return sum_over_ThisPos.GetElement(i);
}

int Hmm_State::get_sum_over_ThisPos_dim(int i) const{
    return sum_over_ThisPos.GetDimension(i);
}

array<int> Hmm_State::get_sum_over_FromPos() const{
    return sum_over_FromPos;
}

int Hmm_State::get_sum_over_FromPos(int i) const{
    return sum_over_FromPos.GetElement(i);
}
 
int Hmm_State::get_sum_over_FromPos_dim(int i) const{
    return sum_over_FromPos.GetDimension(i);
}

array<int> Hmm_State::get_sum_over_pos() const{
    return sum_over_pos;
}

int Hmm_State::get_sum_over_pos(int i) const{
    return sum_over_pos.GetElement(i);
}

int Hmm_State::get_sum_over_pos_dim(int i)const{
    return sum_over_pos.GetDimension(i);
}

int Hmm_State::get_sum_over_FromState() const{
    return sum_over_FromState;
}

int Hmm_State::get_sum_over_FromParam() const{
    return sum_over_FromParam;
}

int Hmm_State::get_number_of_products() const{
    return number_of_products;
}

array<int> Hmm_State::get_product_ThisPos(int index) const{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_ThisPos: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	array<int> temp(0);
	return temp;
    }
    return product_ThisPos[index];
}

int Hmm_State::get_product_ThisPos(int index, int i) const{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_ThisPos: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	return -1;
    }else if((i<0)||(i>=product_ThisPos[index].GetDimension(0)))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_ThisPos: "
	    <<"i("<<i<<") out of range[0.."
	    <<product_ThisPos[index].GetDimension(0)-1<<"]."<<endl;
	return -1;
    }
    return product_ThisPos[index].GetElement(i);
}

int Hmm_State::get_product_ThisPos_dim(int index, int i) const
{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_ThisPos_dim: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	return -1;
    }else if((i<0)||(i>=product_ThisPos[index].GetNumberofDimensions()))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_ThisPos_dim: "
	    <<"i("<<i<<") out of range[0.."
	    <<product_ThisPos[index].GetNumberofDimensions()-1<<"]."<<endl;
	return -1;
    }
    return product_ThisPos[index].GetDimension(i);
}

array<int> Hmm_State::get_product_FromPos(int index) const
{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromPos: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	array<int> temp(0);
	return temp;
    }
    return product_FromPos[index];
}

int Hmm_State::get_product_FromPos(int index, int i) const
{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromPos: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	return -1;
    }else if((i<0)||(i>=product_FromPos[index].GetDimension(0)))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromPos: "
	    <<"i("<<i<<") out of range[0.."
	    <<product_FromPos[index].GetDimension(0)-1<<"]."<<endl;
	return -1;
    }
    return product_FromPos[index].GetElement(i);
}

int Hmm_State::get_product_FromPos_dim(int index, int i)const
{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromPos_dim: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	return -1;
    }else if((i<0)||(i>=product_FromPos[index].GetNumberofDimensions()))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromPos_dim: "
	    <<"i("<<i<<") out of range[0.."
	    <<product_FromPos[index].GetNumberofDimensions()-1<<"]."<<endl;
	return -1;
    }
    return product_FromPos[index].GetDimension(i);
}

int Hmm_State::get_product_FromState(int index) const
{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromState: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	return -1;
    }
    if(!product_FromState)
    {
	return -1;
    }
    return product_FromState[index];
}

int Hmm_State::get_product_FromParam(int index) const
{
    if((index<0)||(index>=number_of_products))
    {
	cout<<"Error: Hmm_State class:: get_prodcut_FromParam: "
	    <<"index("<<index<<") out of range[0.."
	    <<number_of_products-1<<"]."<<endl;
	return -1;
    }
    if(!product_FromParam)
    {
	return -1;
    }
    return product_FromParam[index];
}
                       
Score Hmm_State::get_transition_score(const int j) const
{

  int check = 0;
  if ((j<0) || (j>number_of_states-1))
    {cout << "ERROR class Hmm_State::get_transition_score : index j (" << j << ") out of range (0.." 
	  << number_of_states-1 << ")\n";
    check++;}
  if ((transition_scores.GetNumberofDimensions()!=1)  && 
      (transition_scores.GetDimension(0)!=number_of_states))
    {cout << "ERROR class Hmm_State::get_transition_score : array transition_scores does not exist. \n";
    check++;}

  return(transition_scores.GetElement(j));
}

Score Hmm_State::get_emission_score(const array<int> &indices) const
{

    int check = 0;
    if (emission_scores.GetNumberofDimensions()!= number_of_letters_to_read)
    {
	cout << "ERRO`R class Hmm_State::get_emission_score : array emission_scores does not exist. \n";
	check++;
    }
    if (indices.GetNumberofDimensions()!=1)
    {
	cout << "ERROR class Hmm_State::get_emission_score : indices is no vector\n";
	check++;
    }
    if (indices.GetDimension(0)!=emission_scores.GetNumberofDimensions())
    {
	cout << "ERROR class Hmm_State::get_emission_score : length of vector indices, " << indices.GetDimension(0)
	     << ", does not match number of dimensions of array emission_scores, " 
	     << emission_scores.GetNumberofDimensions() << "\n";
	check++;
    }
    if (check == 0)
    {
	// check that entries of indices are in correct range
	for (int i=0; i<emission_scores.GetNumberofDimensions(); i++)
	{
	    if ((indices.GetElement(i) < 0) || (indices.GetElement(i) > (emission_scores.GetDimension(i)-1)))
	    {
		cout << "ERROR class Hmm_State::get_emission_score : element " << i << " of indices (" 
		     << indices.GetElement(i) << ") out of range (0..."
		     << emission_scores.GetDimension(i)-1 << ").\n";
		check++;
	    }
	}
    }
    return(emission_scores.GetElement(indices));
}

Score Hmm_State::get_emission_score(const int linear_index) const
{
    if (emission_scores.GetNumberofDimensions()!= number_of_letters_to_read )
    {
	cout << "ERROR class Hmm_State::get_emission_score : array emission_scores does not exist. \n";
    }
    int length=emission_scores.GetDimension(0);
    for (int i=1; i<emission_scores.GetNumberofDimensions(); i++)
    {
	length*=emission_scores.GetDimension(i);
    }
    if ((linear_index<0) || (linear_index>(length-1)))
    {
	cout << "ERROR class Hmm_State::get_emission_score : linear_index " << linear_index
	     << " out of range (0.." << length-1 << ").\n";
    }
    return(emission_scores.GetElement(linear_index));
}

Score Hmm_State::get_emission_score(const Sequence* const x,        
					const int             x_position,
					const Sequence* const y,        
					const int             y_position,
					const int linear_index) const 
{
    Score return_score = Logzero;
    
    int check = 0;
    
    if (x == NULL) {
	cout << "ERROR class Hmm_State::get_emission_score : Sequence x is NULL.\n";
	check++;
    }
    if (y == NULL) {
	cout << "ERROR class Hmm_State::get_emission_score : Sequence y is NULL.\n";
	check++;
    }
    
    if (emission_scores.GetNumberofDimensions()!= number_of_letters_to_read) {
	cout << "ERROR class Hmm_State::get_emission_score : array emission_scores does not exist. \n";
    }
    
    int length=emission_scores.GetDimension(0);
    for (int i=1; i<emission_scores.GetNumberofDimensions(); i++) {
	length*=emission_scores.GetDimension(i);
    }
    if ((linear_index<0) || (linear_index>(length-1))) {
	cout << "ERROR class Hmm_State::get_emission_score : linear_index " << linear_index
	     << " out of range (0.." << length-1 << ").\n";
    }
    if (check == 0) {
	return_score = emission_scores.GetElement(linear_index);

	if (special_emission == 1) {
	    
	    Score special_score     = 0; // NB: not Logzero
	    int internal_x_position = x_position;
	    int internal_y_position = y_position;
	    int constraint_x = x->get_constraint();
	    int constraint_y = y->get_constraint();
	    
	    if ((number_of_letters_to_read_x!=0)&&(number_of_letters_to_read_y!=0)) {
		
		if(constraint_x){
		    special_score  = x->get_score(this->get_number_of_child_state_x(),
						  internal_x_position);
		}
		if(constraint_y){
		    special_score += y->get_score(this->get_number_of_child_state_y(),
						  internal_y_position);
		}
	    }else if (number_of_letters_to_read_x != 0) {
		if(constraint_x){
		    special_score = x->get_score(this->get_number_of_child_state_x(),
						 internal_x_position);
		}
	    }else if (number_of_letters_to_read_y !=0) {
		if(constraint_y){
		    special_score = y->get_score(this->get_number_of_child_state_y(),
						 internal_y_position);
		}
	    }

	    //cout<<special_score<<endl;
	    
	    if (special_score <= Logzero) {
		return_score = Logzero;
	    }
	    else {
		if (return_score <= Logzero) {
		    return_score = Logzero;
		}
		else {
		    return_score += special_score;
		}
	    }
	} // if state has special emission
    } // if check == 0
    return(return_score);
}

void Hmm_State:: reset_transition_probs_expressions(int number_of_states)
{
    transition_probs_expressions.SetDimension(0,number_of_states, static_cast<char*>(NULL));
    return;
}

void Hmm_State::reset_transition_probs(void)
{
    transition_probs.ResetData();
    return;
}

void Hmm_State::reset_emission_probs_expression(void)
{
    if(emission_probs_expression) delete[] emission_probs_expression;
    emission_probs_expression=NULL;
    return;
}

void Hmm_State::reset_emission_probs(void)
{
    emission_probs.ResetData();
    return;
}

void Hmm_State::completely_reset_transition_probs_expressions(void)
{
    transition_probs_expressions.ResetDataandDimensions();
    return;
}

void Hmm_State::completely_reset_transition_probs(void)
{
    transition_probs.ResetDataandDimensions();
    return;
}

void Hmm_State::completely_reset_emission_probs_expression(void)
{
    if(emission_probs_expression) delete[] emission_probs_expression;
    emission_probs_expression=NULL;
    return;
}

void Hmm_State::completely_reset_emission_probs(void)
{
    emission_probs.ResetDataandDimensions();
    return;
}

void Hmm_State::reset_transition_scores(void)
{
    transition_scores.ResetData(Logzero);
    return;
}

void Hmm_State::reset_emission_scores(void)
{
    emission_scores.ResetData(Logzero);
    return;
}

void Hmm_State::completely_reset_transition_scores(void)
{
    transition_scores.ResetDataandDimensions();
    return;
}

void Hmm_State::completely_reset_emission_scores(void)
{
    emission_scores.ResetDataandDimensions();
    return;
}

Score Hmm_State::get_viterbi_strip_score(const int j, const int k)
{

    array<int> indices(1);
    indices.SetDimension(0,2);
    indices.SetElement(0,j);
    indices.SetElement(1,k);
    
    return(viterbi_strip.GetElement(indices));
}

Score Hmm_State::get_viterbi_rectangle_score(const int j, const int k)
{

    array<int> indices(1);
    indices.SetDimension(0,2);
    indices.SetElement(0,j);
    indices.SetElement(1,k);
    
    return(viterbi_rectangle.GetElement(indices));
}

Score Hmm_State::get_viterbi_score(const int j, const int k)
{

    array<int> indices(1);
    indices.SetDimension(0,2);
    indices.SetElement(0,j);
    indices.SetElement(1,k);
  
    return(viterbi_scores.GetElement(indices));
}

int Hmm_State::get_number_of_dimensions_of_emission_probs(void)
{
    return(emission_probs.GetNumberofDimensions());          
}

int Hmm_State::get_number_of_dimensions_of_transition_probs(void)
{
    return(transition_probs.GetNumberofDimensions());          
}

int Hmm_State::get_number_of_dimensions_of_transition_scores(void)
{
    return(transition_scores.GetNumberofDimensions());          
}

int Hmm_State::get_number_of_dimensions_of_emission_scores(void)
{
    return(emission_scores.GetNumberofDimensions());          
}

int  Hmm_State::get_dimension_viterbi_strip(const int j)
{
    return(viterbi_strip.GetDimension(j));
}

int  Hmm_State::get_dimension_viterbi_rectangle(const int j)
{
    return(viterbi_rectangle.GetDimension(j));
}

int  Hmm_State::get_dimension_viterbi_scores(const int j)
{
    return(viterbi_scores.GetDimension(j));

}

void Hmm_State::reset_viterbi_strip(void)
{
    viterbi_strip.ResetDataandDimensions();
    return;
}

void Hmm_State::reset_viterbi_rectangle(void)
{
    viterbi_rectangle.ResetDataandDimensions();
    return;
}

void Hmm_State::reset_viterbi_scores(void)
{
    viterbi_scores.ResetDataandDimensions();
    return;
}

void Hmm_State::print( std::ostream &o ) const
{
    o << "special                               : " << special << "\n" << flush;
    o << "special_emission                      : " << special_emission << "\n" << flush;
    o << "mirrored                              : " << mirrored << "\n" << flush;   
    o << "letters to read                       : " << number_of_letters_to_read << "\n" << flush;
    o << "letters to read x                     : " << number_of_letters_to_read_x <<"\n" <<flush;
    o << "letters to read y                     : " << number_of_letters_to_read_y <<"\n"<<flush;
    o << "number of state labels                : " << number_of_state_labels<< "\n"<<flush;

    for(int i=0; i<number_of_state_labels; i++){
	o<<state_labels_type[i]<<"                            : ";
	state_labels[i].Print(o);
    }
    o << "alphabet                              : " << alphabet << "\n" << flush;
    o << "number of state                       : " << number_of_state << "\n" << flush;
    o << "number_of_child_state_x               : " << number_of_child_state_x << "\n" << flush;
    o << "number_of_child_state_y               : " << number_of_child_state_y << "\n" << flush;
    o << "number of states                      : " << number_of_states << "\n" << flush;
    o << "number_of_previous_states             : " << number_of_previous_states << "\n" << flush;
    o << "numbers_of_previous_states            : ";
    numbers_of_previous_states.Print(o);
    o << "number_special_transitions_previous   : " << number_of_special_transitions_to_previous_states << "\n" << flush;
    o << "special_flags_previous_states         : ";
    special_flags_of_transitions_to_previous_states.Print(o);
//  // start new_trans
    o << "numbers_of_from_child_states_x_previous : ";
    numbers_of_from_child_states_x_previous.Print(o);
    o << "numbers_of_from_child_states_y_previous : ";
    numbers_of_from_child_states_y_previous.Print(o);
    o << "numbers_of_to_child_states_x_previous : ";
    numbers_of_to_child_states_x_previous.Print(o);
    o << "numbers_of_to_child_states_y_previous : ";
    numbers_of_to_child_states_y_previous.Print(o);
//  // end new_trans
    o << "offset_for_previous_states_x          : ";
    offset_for_previous_states_x.Print(o);
    o << "offset_for_previous_states_y          : ";
    offset_for_previous_states_y.Print(o);   
    o << "number_of_next_states                 : " << number_of_next_states << "\n" << flush;
    o << "numbers_of_next_states                : ";
    numbers_of_next_states.Print(o);
    o << "\n" << flush;
    o << "number_special_transitions_next       : " << number_of_special_transitions_to_next_states << "\n" << flush;
    o << "special_flags_next_states             : ";
    special_flags_of_transitions_to_next_states.Print(o);
//  // start new_trans
    o << "numbers_of_from_child_states_x_next : ";
    numbers_of_from_child_states_x_next.Print(o);
    o << "numbers_of_from_child_states_y_next : ";
    numbers_of_from_child_states_y_next.Print(o);
    o << "numbers_of_to_child_states_x_next : ";
    numbers_of_to_child_states_x_next.Print(o);
    o << "numbers_of_to_child_states_y_next : ";
    numbers_of_to_child_states_y_next.Print(o);
//  // end new_trans
    o << "offset_for_next_states_x              : ";
    offset_for_next_states_x.Print(o);
    o << "offset_for_next_states_y              : ";
    offset_for_next_states_y.Print(o);
    
    o << "dimensions(transition_probs)          : ";
    transition_probs.PrintDimensions(o);
    o << "transition_probs                      : \n" << flush;
    transition_probs.PrintonlyNonZerowithIndices(o);
    o << "\n" << flush;
    o << "dimensions(emission_probs)            : ";
    emission_probs.PrintDimensions(o);
    o << "emission_probs                        : \n" << flush;
    emission_probs.PrintonlyNonZerowithIndices(o);
    o << "\n" << flush;
    
    o << "dimensions(transition_scores)         : ";
    transition_scores.PrintDimensions(o);
    
    o << "transition_scores                     : \n" << flush;
    transition_scores.PrintonlyNonLogzerowithIndices(cout);
    o << "\n" << flush;
    
    o << "dimensions(emission_scores)           : ";
    emission_scores.PrintDimensions(o);
    
    o << "emission_scores                       : \n" << flush;
    emission_scores.PrintonlyNonLogzerowithIndices(cout);
    o << "\n" << flush;
    return;
}

void Hmm_State::print_transition_probs( std::ostream &o ) const
{
  transition_probs.PrintonlyNonZerowithIndices(o);
  return;
}

void Hmm_State::print_emission_probs( std::ostream &o ) const
{
  emission_probs.PrintonlyNonZerowithIndices(o);
  return;
}

void Hmm_State::print_transition_scores( std::ostream &o ) const
{
  transition_scores.PrintwithIndices(o);
  return;
}

void Hmm_State::print_emission_scores( std::ostream &o ) const
{
  emission_scores.PrintwithIndices(o);
  return;
}

void Hmm_State::print_previous_states( std::ostream &o ) const
{
  numbers_of_previous_states.PrintwithIndices(o);
  return;
}

void Hmm_State::print_special_flags_of_transitions_to_previous_states( std::ostream &o ) const
{
    special_flags_of_transitions_to_previous_states.PrintwithIndices(o);
    return;
}

void Hmm_State::print_offset_for_previous_states_x( std::ostream &o ) const
{
    offset_for_previous_states_x.PrintwithIndices(o);
    return;
}

void Hmm_State::print_offset_for_previous_states_y( std::ostream &o ) const
{
    offset_for_previous_states_y.PrintwithIndices(o);
    return;
}

void Hmm_State::print_next_states( std::ostream &o ) const
{
    numbers_of_next_states.PrintwithIndices(o);
    return;
}

void Hmm_State::print_special_flags_of_transitions_to_next_states( std::ostream &o ) const
{
    special_flags_of_transitions_to_next_states.PrintwithIndices(o);
    return;
}

void Hmm_State::print_offset_for_next_states_x( std::ostream &o ) const
{
    offset_for_next_states_x.PrintwithIndices(o);
    return;
}

void Hmm_State::print_offset_for_next_states_y( std::ostream &o ) const
{
    offset_for_next_states_y.PrintwithIndices(o);
    return;
}

// operators

// ======================================================================

Hmm_State & Hmm_State::operator = (const Hmm_State &p)
{
    int check = 0;
    if ( this != &p)
    {
	special=p.special;
	special_emission=p.special_emission;
	mirrored=p.mirrored;	

	alphabet=p.alphabet;
	number_of_state=p.number_of_state;	
	number_of_states=p.number_of_states;
	state_labels = p.state_labels;
       
	number_of_letters_to_read=p.number_of_letters_to_read;
	number_of_letters_to_read_x=p.number_of_letters_to_read_x;
	number_of_letters_to_read_y=p.number_of_letters_to_read_y; 

	number_of_child_state_x=p.number_of_child_state_x;
	number_of_child_state_y=p.number_of_child_state_y;

	// start new_trans
	numbers_of_from_child_states_x_previous=p.numbers_of_from_child_states_x_previous;
	numbers_of_from_child_states_y_previous=p.numbers_of_from_child_states_y_previous;
	numbers_of_to_child_states_x_previous  =p.numbers_of_to_child_states_x_previous;  
	numbers_of_to_child_states_y_previous  =p.numbers_of_to_child_states_y_previous;  
	numbers_of_from_child_states_x_next    =p.numbers_of_from_child_states_x_next;
	numbers_of_from_child_states_y_next    =p.numbers_of_from_child_states_y_next;
	numbers_of_to_child_states_x_next      =p.numbers_of_to_child_states_x_next;  
	numbers_of_to_child_states_y_next      =p.numbers_of_to_child_states_y_next;  
	// end new_trans

	number_of_previous_states=p.number_of_previous_states;
	special_flags_of_transitions_to_previous_states=p.special_flags_of_transitions_to_previous_states;
	numbers_of_previous_states=p.numbers_of_previous_states;
	number_of_special_transitions_to_previous_states=p.number_of_special_transitions_to_previous_states;
	offset_for_previous_states_x=p.offset_for_previous_states_x;
	offset_for_previous_states_y=p.offset_for_previous_states_y;
	
	number_of_next_states=p.number_of_next_states;
	number_of_special_transitions_to_next_states=p.number_of_special_transitions_to_next_states;
	numbers_of_next_states=p.numbers_of_next_states;
	special_flags_of_transitions_to_next_states=p.special_flags_of_transitions_to_next_states;
	offset_for_next_states_x=p.offset_for_next_states_x;
	offset_for_next_states_y=p.offset_for_next_states_y;
	
	number_of_state_labels = p.number_of_state_labels;
	for(int i=0; i<number_of_state_labels; i++){
	    state_labels[i]=p.state_labels[i];
	}

	state_labels_type = p.state_labels_type;
      
	transition_probs_expressions=p.transition_probs_expressions;
	transition_probs=p.transition_probs;

	emission_probs_expression=p.emission_probs_expression;
	emission_get_from=p.emission_get_from;
	emission_sum_over=p.emission_sum_over;
	emission_product=p.emission_product;
	num_sum_over=p.num_sum_over;
	sum_over_FromState=p.sum_over_FromState;
	number_of_products = p.number_of_products;
	if(product_FromState) delete [] product_FromState;
	product_FromState = NULL;
	if(product_FromParam) delete [] product_FromParam;
	product_FromParam = NULL;
	if(product_ThisPos) delete [] product_ThisPos;
	product_ThisPos = NULL;
	if(product_FromPos) delete [] product_FromPos;
	product_FromPos = NULL;  
	if(p.product_FromState){
	    product_FromState = new int[number_of_products];
	}
	if(p.product_FromParam){
	    product_FromParam = new int[number_of_products];
	}

	if(p.product_ThisPos){
	    product_ThisPos = new array<int>[number_of_products];
	}
	if(p.product_FromPos){
	    product_FromPos = new array<int>[number_of_products];
	}
	for(int i=0; i<number_of_products; i++)
	{
	    product_FromState[i]=p.product_FromState[i];
	    product_FromParam[i]=p.product_FromParam[i];
	    product_ThisPos[i] = p.product_ThisPos[i];
	    product_FromPos[i] = p.product_FromPos[i];
	}
	
	emission_probs=p.emission_probs;
	transition_scores=p.transition_scores;
	emission_scores=p.emission_scores;
	
	// the elements of the following arrays which hold result of a particular run 
	// of an algorithm will not be copied

	viterbi_strip.SetNumberofDimensions(2);    
	viterbi_rectangle.SetNumberofDimensions(2);    
	viterbi_scores.SetNumberofDimensions(2);    
    }
    return(*this);
}
