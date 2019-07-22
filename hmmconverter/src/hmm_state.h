/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: header-file for state-class
*/

#ifndef _hmm_state_h
#define _hmm_state_h

#include "multi_array_2.h"
#include "define.h"
#include "model_parameters.h"
#include "sequence.h"

class Hmm;

class Hmm_State
{

 private: friend class Hmm;
 
 int special;                   // = 1 (state special), = 0 (else) 
 int special_emission;          // = 1 (state has special emissions), = 0 (else) is this neccessary?
 
 int mirrored;                  // = 1 (state mirrored), = 0 (else)
 
 int alphabet;         		// # of letters in alphabet
 int number_of_state;           // number/ id of this state in pairhmm
 int number_of_states; 		// # of states in pairhmm

 int number_of_child_state_x;                          // number of child_state_x in pairhmm (emission scores)
  int number_of_child_state_y;                          // number of child_state_y in pairhmm (emission scores)
 
 int number_of_state_labels;
 array<int>* state_labels;
 char** state_labels_type;

 int number_of_letters_to_read;    // # of letters to read of the state
 int number_of_letters_to_read_x;  // # of letters to read of the state from sequence x
 int number_of_letters_to_read_y;  // # of letters to read of the state from sequence y
  
 // information on transitions to previous states

 int number_of_previous_states;   // number of states which can connect to this state

 int number_of_special_transitions_to_previous_states; // number of previous states for which the transition 
 
 array<int> numbers_of_previous_states;  // numbers of states which can connect to this state

 array<int> special_flags_of_transitions_to_previous_states; 
                                                        // special flags of the transitions from each previous 
                                                        // state to this state
  array<int> offset_for_previous_states_x;              // (for transition scores) empty if this state is mirrored
  array<int> offset_for_previous_states_y;              // (for transition scores) empty if this state is mirrored
  array<int> numbers_of_from_child_states_x_previous;   // (for transition scores) empty if state is mirrored 
  array<int> numbers_of_from_child_states_y_previous;   // (for transition scores) empty if state is mirrored 
  array<int> numbers_of_to_child_states_x_previous;     // (for transition scores) empty if state is mirrored 
  array<int> numbers_of_to_child_states_y_previous;     // (for transition scores) empty if state is mirrored 
  
 // information on transitions to next states
 
 int number_of_next_states;       // number of states to which this state can connect

 int number_of_special_transitions_to_next_states;     // number of next states for which the transition


 array<int> numbers_of_next_states;  // numbers of states to which this state can connect

 array<int> special_flags_of_transitions_to_next_states; 
                                                        // special flags of the transitions from this state 

 array<int> offset_for_next_states_x;                  // (for transition scores) empty if this state is not mirrored
 array<int> offset_for_next_states_y;                  // (for transition scores) empty if this state is not mirrored
 array<int> numbers_of_from_child_states_x_next;       // (for transition scores) empty if this state is not mirrored 
 array<int> numbers_of_from_child_states_y_next;       // (for transition scores) empty if this state is not mirrored 
 array<int> numbers_of_to_child_states_x_next;         // (for transition scores) empty if this state is not mirrored 
 array<int> numbers_of_to_child_states_y_next;         // (for transition scores) empty if this state is not mirrored 

  
 array<char*> transition_probs_expressions;	// array with transition probs expressions 
 array<Prob> transition_probs;       // array with transition probs
 array<Prob> emission_probs;         // array with emission probs
  
 // variables specially defined for this kind of emission prob
 char* emission_probs_expression;         
 int emission_get_from;
 bool emission_sum_over;
 bool emission_product;
 int num_sum_over;
 array<int> sum_over_ThisPos;
 array<int> sum_over_FromPos;
 array<int> sum_over_pos;
 int sum_over_FromState;
 int sum_over_FromParam;

 int number_of_products;
 array<int>* product_ThisPos;
 array<int>* product_FromPos;
 int* product_FromState;
 int* product_FromParam;

 array<Score> transition_scores;          // array with transition scores
 array<Score> emission_scores;            // array with emission scores
 
 // for memory_efficient_viterbi algorithm
 
 array<Score> viterbi_strip;
 array<Score> viterbi_rectangle;
 
 // for viterbi algorithm
 
 array<Score> viterbi_scores;            

 int add_connection_to_new_state(const int number_of_new_state,
				 const int new_next_or_new_previous, // next = 1, previous = -1
				 const int transition_to_new_state_special, 
				 const int number_of_from_child_state_x_new_state,
				 const int number_of_from_child_state_y_new_state,
				 const int number_of_to_child_state_x_new_state,
				 const int number_of_to_child_state_y_new_state,
				 const int offset_x_new_state,
				 const int offset_y_new_state);
 
 void order_next_and_previous_states(); 
 int  permute_x_and_y_indices_of_emission_probs(const int x_first_or_y_first);
 int  set_info_on_transitions_to_next_states(const array<int> ns_of_next_states,
					     const array<int> special_flags_to_next_states);
  
 public:
 
 // constructors and destructor
  
 Hmm_State();
 
 Hmm_State(const int n_of_letter_to_read_x, 
	       const int n_of_letter_to_read_y,
	       const int   alph,  
	       const int   n_of_state,
	       const int   n_of_states,
	       const int   n_of_previous_states,
	       const array<int> ns_of_previous_states,
	       const int   n_of_state_labels,
	       array<int>* const s_labels,   
	       char** const s_labels_type,
	       model_parameters* const MP);  

 Hmm_State(const int n_of_letter_to_read_x, 
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
	       model_parameters* const MP);  
 
 Hmm_State(const Hmm_State &p);
 ~Hmm_State(void);
  
  // member_functions
  
  int modify_state(const int new_n_of_states,
		   const int shift_state_numbers,
		   const int n_of_fixed_state,
		   const int new_n_of_fixed_state);

  int set_alphabet(int const alph);
  int set_mirrored(int const mir);
  int set_number_of_state_labels(int const n_s_labels);
  int set_state_labels_type(char** const s_label_type);
  int set_state_labels_type(char* const s_label_type, const int i);

  int set_number_of_state(int ns);
  int set_number_of_states(int n);

  int  set_transition_probs_expression(const int j, char* const exp);
  int  set_transition_prob(const int j, const Prob x);

  int  set_special_emissions(void);
  int  set_special_emissions(const Hmm_State* const child_state_x,
			    const Hmm_State* const child_state_y);
  int  set_special_emissions(const Hmm_State* const child_state);
  int  set_special_emissions(const int n_of_child_state_x,
			    const int n_of_child_state_y);
  int  set_special_emissions(const int n_of_child_state);

  int  set_emission_probs_expression(const char* exp);
  int  set_emission_prob(const array<int> &indices, const Prob x);
  int  set_emission_prob(const int linear_index, const Prob x);
  int  set_emission_prob(const array<Prob> prob);
 
  void set_emission_sum_over(bool input_sum_over);
  void set_emission_product(bool input_product);
  void set_num_sum_over(int input_num_sum);
  int set_emission_get_from(int from);
  int set_sum_over_ThisPos(const char* cur_pos); 
  int set_sum_over_FromPos(const char* corr_pos);
  int set_sum_over_pos(const char* sum_pos);
  int set_sum_over_FromState(int fromstate);
  int set_sum_over_FromParam(int fromparam);

  int set_number_of_products(int n_of_products);
  int set_product_ThisPos(int index, const char* cur_pos);
  int set_product_FromPos(int index, const char* corr_pos);
  int set_product_FromState(int index, int fromstate);
  int set_product_FromParam(int index, int fromparam);

  int  set_transition_score(const int j, const Score x);
  int  set_emission_score(const array<int> &indices, const Score x);

  int  set_viterbi_strip_score(const int j, const int k, const Score x);
  int  set_viterbi_rectangle_score(const int j, const int k, const Score x);
  int  set_viterbi_score(const int j, const int k, const Score x);

  int  set_dimensions_viterbi_strip(const int j, const int k);
  int  set_dimensions_viterbi_rectangle(const int j, const int k);
  int  set_dimensions_viterbi_scores(const int j, const int k); 

  int  set_emission_probs_below_threshold_to_threshold(const Prob threshold);
  
  int  get_mirrored_copy_of(const Hmm_State &p);
  int  get_complemented_copy_of(const Hmm_State &p); 

  int get_alphabet(void) const;
  int get_special(void) const;
  int get_special_emission(void) const;
  int get_mirrored(void) const;
  int get_number_of_state_labels(void) const;
  char* get_state_labels_type(const int i) const; // returns -1 when input out of range
  int get_state_labels_dim(const int i)const;
  int get_state_labels(const int i, const int j) const; // returns -1 when input out of range
  
  int get_number_of_state(void) const;
  int get_number_of_child_state_x(void) const; 
  int get_number_of_child_state_y(void) const; 
  int get_number_of_states(void) const;

  int get_number_of_previous_states(void) const;
  int get_number_of_special_transitions_to_previous_states(void) const;
  int get_number_of_previous_state(const int i) const; 
  int get_index_of_previous_state(const int number_of_previous_state) const; // returns -1 if number_of_previous_state is no previous state
  int get_special_flag_of_transition_to_previous_state(const int i) const;
  int get_offset_for_previous_state_x(const int i) const;
  int get_offset_for_previous_state_y(const int i) const;

  int get_number_of_next_states(void) const;
  int get_number_of_special_transitions_to_next_states(void) const; 
  int get_number_of_next_state(const int i) const; 
  int get_index_of_next_state(const int number_of_next_state) const;  // returns -1 if number_of_next_state is no next state
  int get_special_flag_of_transition_to_next_state(const int i) const;
  int get_offset_for_next_state_x(const int i) const; 
  int get_offset_for_next_state_y(const int i) const; 

  int get_letters_to_read(void) const;
  int get_letters_to_read_x(void) const;
  int get_letters_to_read_y(void) const;

  int is_state_previous_state(const int i) const;  // 1 (yes), 0 (no)
  int is_transition_to_previous_state_special(const int i) const; // 1 (yes), 0 (no)
  int get_number_of_from_child_state_x_previous(const int number_of_previous_state) const; // returns -1 if number_of_previous_state is no previous state
  int   get_number_of_from_child_state_y_previous(const int number_of_previous_state) const; // returns -1 if number_of_previous_state is no previous state
  int   get_number_of_from_child_state_x_next(const int number_of_next_state) const; // returns -1 if number_of_next_state is no next state
  int   get_number_of_from_child_state_y_next(const int number_of_next_state) const; // returns -1 if number_of_next_state is no next state
  int   get_number_of_to_child_state_x_previous(const int number_of_previous_state) const; // returns -1 if number_of_previous_state is no previous state
  int   get_number_of_to_child_state_y_previous(const int number_of_previous_state) const; // returns -1 if number_of_previous_state is no previous state
  int   get_number_of_to_child_state_x_next(const int number_of_next_state) const; // returns -1 if number_of_next_state is no next state
  int   get_number_of_to_child_state_y_next(const int number_of_next_state) const; // returns -1 if number_of_next_state is no next state

  int   is_state_next_state(const int i) const;  // 1 (yes), 0 (no)
  int   is_transition_to_next_state_special(const int i) const;  // 1 (yes), 0 (no)
  int   is_state_special(void) const; // 1 (yes), 0 (no) 
  int   has_state_special_emissions(void) const; // 1 (yes), 0 (no)  // checked
  int   is_state_mirrored(void) const; // 1 (yes), 0 (no) 

  int   exponentiate_matrix_of_emission_probs(const int exponent);
  int   check_matrix_of_emission_probs() const; 

  char* get_transition_probs_expressions(const int j) const;
  Prob  get_transition_prob(const int j) const;
  Prob  get_special_transition_prob(const Hmm_State* const to,
				    const Sequence*      const x,
				    const int                  x_position,
				    const Sequence*      const y,
				    const int                  y_position) const;
  int   get_scores_for_special_transition(// input
					  const Hmm_State* const to,
					  const Sequence*      const x,
					  const int                  x_position,
					  const Sequence*      const y,
					  const int                  y_position,
					  // output 
					  Score* x_sequence_score,
					  Score* y_sequence_score) const;
  Prob  get_special_transition_prob_with_print(const Hmm_State* const to,
					       const Sequence*      const x,
					       const int                  x_position,
					       const Sequence*      const y,
					       const int                  y_position) const;

  Prob  get_emission_prob(const array<int> &indices) const;  
  Prob  get_emission_prob(const int linear_index) const;   
  array<Prob> get_emission_prob() const;
	
  Prob* get_emission_prob_array() const;

  char* get_emission_probs_expression() const;
  int get_emission_get_from() const;
  bool get_emission_sum_over() const;
  bool get_emission_product() const;
  int get_num_sum_over() const;
  array<int> get_sum_over_ThisPos() const; 
  int get_sum_over_ThisPos(int i) const;
  int get_sum_over_ThisPos_dim(int i) const;
  array<int> get_sum_over_FromPos() const;
  int get_sum_over_FromPos(int i) const;
  int get_sum_over_FromPos_dim(int i) const;
  array<int> get_sum_over_pos() const;
  int get_sum_over_pos(int i) const;
  int get_sum_over_pos_dim(int i) const;
  int get_sum_over_FromState() const;
  int get_sum_over_FromParam() const;
  int get_number_of_products() const;
  array<int> get_product_ThisPos(int index) const;
  int get_product_ThisPos(int index, int i) const;
  int get_product_ThisPos_dim(int index, int i) const;
  array<int> get_product_FromPos(int index) const;
  int get_product_FromPos(int index, int i) const;
  int get_product_FromPos_dim(int index, int i) const;
  int get_product_FromState(int index) const;
  int get_product_FromParam(int index) const;

  Score get_transition_score(const int j) const;
  Score get_emission_score(const array<int> &indices) const; 
  Score get_emission_score(const int linear_index) const;

  Score get_emission_score(const Sequence* const x,        
			   const int             x_position,
			   const Sequence* const y,        
			   const int             y_position,
			   const int linear_index) const; 

  void  reset_transition_probs_expressions(int number_of_states);
  void  reset_transition_probs(void);
  void  reset_emission_probs_expression(void);
  void  reset_emission_probs(void);
  void  completely_reset_transition_probs_expressions(void);
  void  completely_reset_transition_probs(void);
  void  completely_reset_emission_probs_expression(void);
  void  completely_reset_emission_probs(void);

  void  reset_transition_scores(void);
  void  reset_emission_scores(void);
  void  completely_reset_transition_scores(void);
  void  completely_reset_emission_scores(void);

  Score get_viterbi_strip_score(const int j, const int k);
  Score get_viterbi_rectangle_score(const int j, const int k);
  Score get_viterbi_score(const int j, const int k);

  int   get_number_of_dimensions_of_transition_probs(void);
  int   get_number_of_dimensions_of_emission_probs(void);
  int   get_number_of_dimensions_of_transition_scores(void);
  int   get_number_of_dimensions_of_emission_scores(void);

  int   get_dimension_viterbi_strip(const int j);
  int   get_dimension_viterbi_rectangle(const int j);
  int   get_dimension_viterbi_scores(const int j);

  void  reset_viterbi_strip(void);
  void  reset_viterbi_rectangle(void);
  void  reset_viterbi_scores(void);

  void print( std::ostream &o ) const;
  void print_transition_probs( std::ostream &o ) const;
  void print_emission_probs( std::ostream &o ) const;
  void print_emission_probs_with_letter_indices( std::ostream &o ) const; 
  void print_transition_scores( std::ostream &o ) const;
  void print_emission_scores( std::ostream &o ) const;
  void print_previous_states( std::ostream &o ) const;
  void print_special_flags_of_transitions_to_previous_states( std::ostream &o ) const;
  void print_offset_for_previous_states_x( std::ostream &o ) const;
  void print_offset_for_previous_states_y( std::ostream &o ) const;

  void print_next_states( std::ostream &o ) const;
  void print_special_flags_of_transitions_to_next_states( std::ostream &o ) const;
  void print_offset_for_next_states_x( std::ostream &o ) const;
  void print_offset_for_next_states_y( std::ostream &o ) const;
  
  // operators

  Hmm_State & operator = (const Hmm_State &p);

};

#endif 

