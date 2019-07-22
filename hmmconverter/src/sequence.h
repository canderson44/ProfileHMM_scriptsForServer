/* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare sequence-class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/sequence.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _sequence_h
#define _sequence_h

#include <stdio.h>  
#include "multi_array_2.h"
#include "define.h"
#include "model_parameters.h"
#include "match.h"

class Info;

class Sequence
{
 private:
    
    int sequence_type; 
    int *sequence;                 
 
    int start_position;               /* start position of sequence
					 (always <= end_position) */
    int end_position;                 /* end position of sequence */
    int length_of_sequence;           /* length of sequence */
    char *ac;                         /* AC of the sequence */

    int    constraint;  // whether the sequence is the one with constraint states

    int*   sp_xsteps;
    int*   sp_ysteps;
    int*   sp_states;
    Score* sp_scores;
    
    int    sp_steps;

    //annotation

    int number_of_sources;
    int number_of_annotation_labels;
    int* number_of_each_annotation_label;
    bool**** annotation_labels;   // (source, pos, label_set, label)
    double**** annotation_label_scores;
    bool* score_labels;
        
    int  number_of_implemented_special_emissions;
    int  number_of_special_emissions;
    
    array<int> implementation_status_of_special_emissions; 
    // 1d-array: index 1: state number (un-mirrored hmm)
    //           entry  : 0 (info on special emission not yet implemented), 1 (implemented)
    
    array<int> indices_of_special_emissions;
    // 1d-array: index 1: state number (un-mirrored hmm)
    //           entry  : index of this special emission to be used as index 1
    //                    for array posterior_probs_for_special_emission
    
    array<Prob> posterior_probs_for_special_emissions;
    // 2d-array: index 1: index of special emission
    //           index 2: position within this sequence
    //           entry  : posterior-prob for emission index 1 at position
    //                    index 2 in the sequence

    array<Score> scores_for_special_emissions;
    // 2d-array: index 1: index of special emission
    //           index 2: position within this sequence
    //           entry  : score for emission index 1 at position
    //                    index 2 in the sequence
    // note:     score = log_2(prob(model)/prob(random_model))
    
    int         number_of_implemented_special_transitions;
    int         number_of_special_transitions; 

    array<int> implementation_status_of_special_transitions; 
    // 2d-array: index 1: from-state number (un-mirrored hmm)
    //           index 2: to-state number   (un-mirrored hmm)
    //           entry  : 0 (info on special transition from->to not yet implemented), 1 (implemented)
    
    array<int> indices_of_special_transitions;
    // 2d-array: index 1: from-state number (un-mirrored hmm)
    //           index 2: to-state number   (un-mirrored hmm)
    //           entry  : index of this special transition to be used as index 1
    //                    for array posterior_probs_for_special_states

    
    array<Prob> posterior_probs_for_special_transitions;
    // 2d-array: index 1: index of special transition
    //           index 2: position within this sequence
    //           entry  : posterior-prob for transition index 1 at position
    //                    index 2 in the sequence
    
    array<Prob> priors_of_special_transitions; // buh new_trans
    // 1d-array: index 1: index of special transition
    //           entry  : prior for transition index 1

    array<Score> scores_for_special_transitions;
    // 2d-array: index 1: index of special transition
    //           index 2: position within this sequence
    //           entry  : score for transition index 1 at position
    //                    index 2 in the sequence
    // note:     score = log_2(prob(model)/prob(random_model))

    int set_implementation_status_of_special_transition(const int number_of_from_state,
							const int number_of_to_state,
							const int status); 
    int set_index_of_special_transition(const int number_of_from_state,
					const int number_of_to_state,
					const int index); 
    int set_implementation_status_of_special_emission(const int number_of_state,
						      const int status);
    int set_index_of_special_emission(const int number_of_state,
				      const int index); 
    
 public:
    
    // constructors and destructor

    Sequence();            
    Sequence(const char* const char_seq,
	     model_parameters* const MP);   
    Sequence(const int seq_type, 
	     const char* const char_seq,
	     model_parameters* const MP); 
    Sequence(const int* const seq, 
	     const int length);   
    Sequence(const int seq_type,
	     const int* const seq, 
	     const int length, 
	     const bool probs_or_scores,
	     model_parameters* const MP);      
    Sequence(const Sequence &s); 
  

    ~Sequence();                              
  
    // All accessor functions
    
    int     get_subsequence(const int start_pos, const int end_pos,
			    Sequence* seq) const; 
    int     get_subsequence(const int start_pos, const int end_pos,
			    int** subsequence); 
    int     get_subsequence_rel(const int rel_start_pos, const int rel_end_pos,
				Sequence* seq) const; 
    int     get_subsequence_rel(const int rel_start_pos, const int rel_end_pos,
				int** subsequence); 
    char*   get_ac(void) const;                              
    int     get_start_position(void) const;                 
    int     get_end_position(void) const;    
    int     get_number_of_sources(void) const;
    int     get_number_of_annotation_labels(void) const;
    int*    get_number_of_each_annotation_label(void) const;
    int*    get_sequence(void) const;                
    char*   get_char_sequence(model_parameters* const MP) const;
    int     length(void) const;  
    int     letter(const int i) const;
    const char* get_sequence_type(void) const;
    int     get_sequence_type_int(void) const; 
    int     get_number_of_special_transitions(void) const;
    int     get_number_of_implemented_special_transitions(void) const;
    int     get_index_of_special_transition(const int number_of_from_state,
					    const int number_of_to_state) const;
    // note: the above function returns value >= 0 upon successful completion
    //       and value = -1 if error occurred
    int     get_implementation_status_of_special_transition(const int number_of_from_state,
							    const int number_of_to_state) const;
    // note: the above function returns value >= 0 upon successful completion
    //       and value = -1 if error occurred
    Prob     get_prior_of_special_transition(const int number_of_from_state,
					     const int number_of_to_state) const;
    // note: the above function returns value >= 0 upon successful completion
    //       and value = -1 if error occurred
    Prob    get_posterior_prob(const int number_of_from_state,
			       const int number_of_to_state,
			       const int position_in_sequence) const;
    Score   get_score(const int number_of_from_state,
		      const int number_of_to_state,
		      const int position_in_sequence) const;
  
    int     get_number_of_special_emissions(void) const;
    int     get_number_of_implemented_special_emissions(void) const;
    int     get_index_of_special_emission(const int number_of_state) const;
    // note: the above function returns value >= 0 upon successful completion
    //       and value = -1 if error occurred
    int     get_implementation_status_of_special_emission(const int number_of_state) const; // checked
    // note: the above function returns value >= 0 upon successful completion
    //       and value = -1 if error occurred
    Prob    get_posterior_prob(const int number_of_state,
			       const int position_in_sequence) const;
    Score   get_score(const int number_of_state,
		      const int position_in_sequence) const;
    int     get_next_annotation_labels_rel(//input
	const int        rel_position, 
	bool***    const label,
	double***  const label_score,
	//output
	int*             next_rel_position,
	bool****   const next_label,
	double**** const next_label_score) const;
    int     get_annotation_labels(bool*****   const label,
				  double***** const label_score) const;
    int     get_combine_annotation_labels(bool****   const label,
					  double**** const score);
    int     get_annotation_labels_for_position(const int abs_position, 
				 // output
				 bool****   const label,
				 double**** const label_score) const;
    int     get_annotation_labels_for_position_rel(const int rel_position, 
				     // output 
				     bool****   const label,
				     double**** const label_score) const;
    int     get_annotation_labels_for_set(const int        label_set, 
					  bool****   const label,
					  double**** const label_score) const;
    int     get_annotation_labels_for_source(const int        source,
					     bool****   const label,
					     double**** const label_score) const;
    bool    get_score_labels(const int i) const;
    int     get_constraint() const;
    
    int     print_non_alphabet_contents_of_sequence(model_parameters* const MP,
						    std::ostream &o) const; 
    int     replace_non_alphabet_contents_of_sequence(model_parameters* MP); 
    int     annotation_labels_exist() const; 

    // All mutators

    int     set_sp_xsteps(const int steps, const int* const xsteps_array);
    int     set_sp_ysteps(const int steps, const int* const ysteps_array);
    int     set_sp_scores(const int steps, const Score* const scores_array);
    int     set_sp_states(const int steps, const int* const states_array);  
    int     set_sp_steps(const int steps); 
    Score   get_sp_score(const int xsteps, const int ysteps, const int state) const;
    int     get_sp_steps(void);
    int     set_annotation(const char* file_name,
			   model_parameters* const MP);           
    int     change_annotation_labels(const int source,
				     const int label_set,
				     const int old_label,
				     const int new_label);
    int     change_annotation_labels(const int label_set,
				     const int old_label,
				     const int new_label);
    int     change_all_annotation_labels(const int source,
					 const int label_set,
					 const int new_label);
    int     change_all_annotation_labels(const int label_set,
					 const int new_label);
    void    delete_annotation(void);                             
    int     set_annotation_labels(const int   source,
				  bool***   const ann_labels,
				  double*** const ann_label_scores); 
    int     set_annotation_labels(bool****   const ann_labels,
				  double**** const ann_label_scores);
    int     set_score_labels(const int i, const bool s); 
    int     set_constraint(const int c);
    int     set_ac(const char* const n);
    int     set_sequence_start_and_end(const int sequence_start,
				       const int sequence_end);       
    int     set_number_of_special_emissions(const int n_of_special_emissions);
    int     set_number_of_sources(const int n_of_source) ;
    int     set_number_of_annotation_labels(const int n_of_annotation_labels) ;
    int     set_number_of_each_annotation_label(const int n_of_annotation_labels,
						int* const n_of_each_annotation_label);

    int     set_posterior_probs(const int   number_of_state,
				const Prob* probs);     
    int     set_posterior_probs(const int number_of_new_state,
				const int number_of_already_implemented_state);
    int     set_scores(const int    number_of_state, 
		       const Score* scores); 
    int     set_scores(const int number_of_new_state, 
		       const int number_of_already_implemented_state);
    int  change_scores(const Score new_default_score);
    int print_sequence_to_fasta_file(model_parameters* const MP,
				     std::ostream &o) const; 

    void print_emission_information_linewise(std::ostream &o,
					     model_parameters* const MP) const;
    void print_emission_and_transition_information_linewise(std::ostream &o,
							    model_parameters* const MP) const;
    void print_annotation_linewise(std::ostream &o,
				   model_parameters* const MP) const;
    void print_annotation_linewise(std::ostream &o, 
				   const int start_rel_pos, 
				   const int end_rel_pos,
				   model_parameters* const MP) const;
    void print_char_sequence(model_parameters* const MP, std::ostream &o) const;
    void print_sequence(std::ostream &o) const;
    void print_indices_of_implemented_special_transitions(std::ostream &o) const;
    void print_indices_of_implemented_special_emissions(std::ostream &o) const;
    void print_non_zero_posterior_probs(std::ostream &o) const;
    void print_non_Logzero_scores(std::ostream &o) const;
    void print_non_Logzero_emission_scores(std::ostream &o) const;
    void print_emission_scores(std::ostream &o) const;
    void print_emission_scores_with_sequence_info(std::ostream &o,
						  model_parameters* const MP) const;
    void print_implemented_special_transitions(std::ostream &o) const;
    void print_implemented_special_emissions(std::ostream &o) const;
    void print_non_zero_posterior_probs_with_sequence_info(std::ostream &o,
							   model_parameters* const MP) const;
    void print_non_Logzero_scores_with_sequence_info(model_parameters* const MP,
						     std::ostream &o ) const;
    void print_non_zero_posterior_probs_with_sequence_info(const int number_of_from_state,
							   const int number_of_to_state,
							   std::ostream &o,
							   model_parameters* const MP) const;
    int  get_number_of_features_in_sequence(const int source,
					    const char* const file_name,
					    const int feature_label_index,
					    const int feature,
					    model_parameters* const MP);
    void print_non_Logzero_scores_with_sequence_info(const int number_of_from_state,
						     const int number_of_to_state,
						     std::ostream &o,
						     model_parameters* const MP) const;
    void print_non_zero_posterior_probs_with_sequence_info(const int number_of_from_state,
							   const int number_of_to_state,
							   const int position_in_sequence,
							   std::ostream &o,
							   model_parameters* const MP) const;
    void print_non_Logzero_scores_with_sequence_info(const int number_of_from_state,
						     const int number_of_to_state,
						     const int position_in_sequence,
						     std::ostream &o,
						     model_parameters* const MP) const;
    int calculate_abs_pos(const int       rel_start,
			  const int       rel_end,
			  int*      const abs_start,
			  int*      const abs_end) const;
    int read_special_file(// input
	const char*     const file_name,
	model_parameters*  const MP,
	// output
	Info* const info);

    void print(model_parameters* const MP,
	       std::ostream &o) const;
    void print_short(void) const;
  
    // operators

    Sequence & operator = (const Sequence &p);
};

#include "info.h"

#endif
