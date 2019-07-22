 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: header-file for pairhmm class

   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/hmm.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _hmm_h
#define _hmm_h

#include <stdio.h>
#ifndef _INTEL
#include <vector>
#include <vector.h>
#include <map.h>
#endif // #ifndef _INTEL

#include "tube.h"
#include "multi_array_2.h"
#include "define.h"
#include "model_parameters.h"
#include "sequence.h"
#include "hmm_state.h"
#include "input_from_files.h"
#include "scoring_functions.h"

using namespace std;

class TransitionProb;

class EmissionProb;

class ViterbiElement {
  
public:

  int     prevState;
  int     nParent;
  Score   score;

  // constructor and destructor

  ViterbiElement() : prevState(-1), nParent(0), score(Logzero) {;}
  ~ViterbiElement() {;}

};

class StatePath {

 public:

    int    l_steps;
    Score  l_score;
    int*   l_states; 
    int*   l_xsteps;
    int*   l_ysteps;
    Score* l_scores;
    
    // constructor and destructor
    
    StatePath() : l_steps(0), l_score(0.0), l_states(NULL), l_xsteps(NULL), l_ysteps(NULL), l_scores(NULL) {;}
    
    ~StatePath() {

	if (l_states) delete [] l_states;
	if (l_xsteps) delete [] l_xsteps;
	if (l_ysteps) delete [] l_ysteps;      
	if (l_scores) delete [] l_scores;
	
	l_states = NULL;
	l_xsteps = NULL;
	l_ysteps = NULL;
	l_scores = NULL;
	l_steps  = 0;
	l_score  = 0;
    }
    
};

class ForBackElement {
  
public:

    Score   score;
  
    // constructor and destructor
    
    ForBackElement() : score(Logzero) {;}
    ~ForBackElement() {;}
    
};

class Hmm
{
 private:
    
    int mirrored;                                         // = 1 (pairhmm mirrored), = 0 (else)
    int alphabet;
    int number_of_states;
    bool pair;
    Hmm_State* model; 
    
    // variables for results of alignment algorithms
    
    int steps;                           /* counts the # of steps in
					    real pairhmm for the
					    alignment */
    Score score;                         // score of the alignment 
    Score* sequence_of_scores;           /* score for every step in
					    alignment (non cumulative) */
    int* sequence_of_states;             // state for every step in alignment
    int* sequence_of_xsteps;
    int* sequence_of_ysteps;

    int condensed_steps;

    Score* condensed_sequence_of_scores;
    int*   condensed_sequence_of_states;    
    int*   condensed_sequence_of_xsteps;
    int*   condensed_sequence_of_ysteps;

    // variables needed for the standard random pairhmm

    Prob* standard_random_probs;         /* [letter_x] is the prob for
					    randomly emitting letter x */
    Prob eta;                            /* standard prob in standard
					  random pairhmm (R. Durbin,
					  p. 83) */

    double* superfastexptable;


    // private functions

    int copy_rectangle_from_strip_to_next_strip(Score**** const strip,
						const int x_start, const int x_end,
						const int y_start, const int y_end,
						const int x_start_copy, const int x_end_copy,
						const int y_start_copy, const int y_end_copy,
						Score**** const next_strip,
						const int x_start_next, const int x_end_next,
						const int y_start_next, const int y_end_next,
						const int x_start_copy_next, const int x_end_copy_next,
						const int y_start_copy_next, const int y_end_copy_next,
						const int direction);

    int allocate_memory_for_strip(const Hmm *s, 
				  const int x_margin,
				  const int y_start, const int y_end,
				  Score**** const strip);
    
    int check_that_there_are_valid_values_in_strip(Score**** const strip,
						   const Hmm *s, 
						   const int x_margin,
						   const int y_start, const int y_end);	
    int delete_memory_for_strip(const Hmm *s, 
				const int x_margin,
				Score**** const strip);
    int allocate_memory_for_viterbi_rectangle(// input values       
	const Hmm *s, 
	const int x_start, const int x_end,
	const int y_start, const int y_end,
	const int direction,
	// output values;
	Score**** const viterbi_rectangle);

    int allocate_memory_for_state_path(const int x_start, const int x_end,
				       const int y_start, const int y_end,
				       Score** const local_sequence_of_scores,
				       int** const local_sequence_of_states,
				       int** const local_sequence_of_xsteps,
				       int** const local_sequence_of_ysteps);

    int allocate_memory_for_state_path(const int max_length_of_state_path,
				       const int default_steps,
				       const Score default_overall_score,
				       const int default_int,
				       const Score default_score); // buh

    int delete_memory_for_state_path(Score** const local_sequence_of_scores,
				     int** const local_sequence_of_states,
				     int** const local_sequence_of_xsteps,
				     int** const local_sequence_of_ysteps);

    int delete_memory_for_viterbi_rectangle(const Hmm *s, 
					    const int x_start, const int x_end,
					    const int direction,
					    Score**** const viterbi_rectangle);
    

    int calculate_viterbi_rectangle(Score*** const viterbi_rectangle,
				    const Hmm *s, 
				    const Hmm *s_unmirrored,
				    const Sequence *x, const int x_start, const int x_end,
				    const Sequence *y, const int y_start, const int y_end,
				    const int* start_states, 
				    const int number_of_start_states, 
				    int x_margin, 
				    int y_margin, 
				    const int direction);
    
    int retrieve_state_path_from_viterbi_rectangle(Score*** const viterbi_rectangle,
						   const Hmm *s, 
						   const Hmm *s_unmirrored,
						   const Sequence *x,
						   const int x_start,
						   const int x_end,
						   const Sequence *y,
						   const int y_start, 
						   const int y_end,
						   const int end_state, 
						   const int direction,
						   int& local_steps,
						   Score& local_score,
						   Score* const local_sequence_of_scores,
						   int* const local_sequence_of_states,
						   int* const local_sequence_of_xsteps,
						   int* const local_sequence_of_ysteps);
    
    int get_strips_and_merge(// input values
	const Hmm *mirror, 
	const Sequence *x, const int x_start, const int x_end,
	const Sequence *y, const int y_start, const int y_end,
	const int start_state, const int end_state,
	const int x_midpoint, const int min_strip_width,
	// output values			    
	int* max_i, int* max_j, int* max_state_forward, int* max_state_backward);
    
    int new_get_strips_and_merge(// input values
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
	int* max_i, int* max_j, int* max_state_forward, int* max_state_backward);

    
    int get_strip(Score**** const viterbi_strip, 
		  const Hmm *s, 
		  const Hmm *unmirrored_s,
		  const Sequence *x, const int x_start, const int x_end,
		  const Sequence *y, const int y_start, const int y_end, 
		  const int strip_width,
		  const int x_margin,
		  const int y_margin,
		  const int direction);
    
    int memory_viterbi(const Hmm* mirror,		
		       const Sequence *x, const int x_start, const int x_end,
		       const Sequence *y, const int y_start, const int y_end,
		       const int start_state, const int end_state,
		       const int min_strip_width, const int max_area);
    
    int new_memory_viterbi(// input values
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
	int* start_state_traceback);
    
    int add_local_solution(
	const Sequence *x, 
	const Sequence *y, 
	const int& local_steps,
	const Score& local_score,
	const Score* const local_sequence_of_scores,
	const int* const local_sequence_of_states,			 
	const int* const local_sequence_of_xsteps,
	const int* const local_sequence_of_ysteps,
	const int direction);

    int add_local_state_path(const Sequence *x, 
			     const Sequence *y, 
			     const StatePath *local_state_path);

    int check_consistency_of_solution(const Sequence *x, const Sequence *y);
    
    int viterbi_rectangle(Score*** viterbi_rectangle_array,
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
			  int* start_state_traceback);

    int get_sequence_labels_for_state_path(
	int*** const x_labels, int*** const y_labels,
	const Sequence &x, const Sequence &y,
	model_parameters* const MP); 

    void readsuperfastexptable();

    double superfastinterpexp(const double x);
    
    // requires non-condensed state path solution to exist
    // memory for x_labels and y_labels is allocated within the function and has to be released outside the function
    
 public:
    
    // constructors and destructors
  
    Hmm();
    Hmm(model_parameters* const MP);
    Hmm(const Hmm &p);
    Hmm(const Hmm &p, int mirror, int change_strand);
    Hmm(const Hmm &p_1, const int number_of_merge_state_1,
	    const Hmm &p_2, const int number_of_merge_state_2,
	    const Prob merge_prob);
    ~Hmm();
    
    // member functions
    
    int   get_alphabet(void) const;
    int   get_mirrored(void) const;
    int   get_number_of_states(void) const;
    
    bool  get_pair(void) const;

    // accessors
    int   set_alphabet(const int a);
    int   set_mirrored(const int m);
    int   set_number_of_states(const int n_s);
    int   set_pair(const bool p);

    
    Score get_transition_score(
	const Hmm *unmirrored_s,
	const int             n_of_from_state,
	const int             n_of_to_state,
	const Sequence* const x,
	const int             x_position,
	const Sequence* const y,
	const int             y_position) const;

    Score new_get_transition_score(
	const Hmm *unmirrored_s,
	const int             n_of_from_state,
	const int             n_of_to_state,
	const Sequence* const x,
	const int             x_position,
	const Sequence* const y,
	const int             y_position) const; 
    
    int   implement_state(int i, Hmm_State* state);
    int   build_connected_pairhmm() const;
    
    int   check_consistency() const;
    int   check_consistency_of_probs() const;
    
    void  print(std::ostream &o ) const;
    
    int   get_steps(void) const;
    Score get_score(void) const;
    int   get_state_in_alignment(int i) const;
    Score get_local_score_in_alignment(int i) const;
    
    void  print_results(std::ostream &o) const;
    
    void  output1(std::ostream &o, 
		  model_parameters * const MP,
		  const char* seq_name_x,
		  const char* seq_name_y) const;

    void  output2(std::ostream &o, 
		  model_parameters * const MP,
		  const char* seq_name_x,
		  const int seq_length_x,
		  const char* seq_name_y,
		  const int seq_length_y) const;

    void transition_prob_output(std::ostream &o,
				TransitionProb* const TP) const;

    void XML_output(const char* xmlfile,
		    std::ostream &o,
		    model_parameters* const MP,
		    TransitionProb* const TP) const;
    
    void emission_prob_output(std::ostream &o,
			      model_parameters* const MP,
			      EmissionProb* const EP) const;

    int   condense_solution(int discard_long_solution);
    
    int viterbi(
	const Sequence *x, const Sequence *y); 
    
    int new_viterbi(
	const Sequence *x, const Sequence *y);
    
    int Hirschberg_viterbi(
	const Hmm* mirror,
	const Sequence *x, const Sequence *y,
	const int max_area);
    
    int Hirschberg_viterbi_rectangle(Score*** viterbi_rectangle_array,
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
				     int* start_state_traceback);

    int heuristic_viterbi(const Hmm*   mirror,
			  const Sequence *x, const int* x_coordinates,
			  const Sequence *y, const int* y_coordinates,
			  const int number_of_xy_pairs,
			  const int x_margin,
			  const int y_margin,
			  const int max_area);

    int new_heuristic_viterbi(const Hmm*   mirror,
			      const Sequence *x, const int* x_coordinates,
			      const Sequence *y, const int* y_coordinates,
			      const int number_of_xy_pairs,
			      const int x_margin,
			      const int y_margin,
			      const int max_area);
    
    int calculate_scores_from_probs();

    // operators
    
    Hmm & operator = (const Hmm &p);
    Hmm_State* operator [] (const int i); 
    
    // **********************************************************************
    
    void set_information_of_fake_alignment(const Score  c_score,
					   const int    c_steps,
					   const int*   c_sequence_of_states,
					   const int*   c_sequence_of_xsteps,
					   const int*   c_sequence_of_ysteps,
					   const Score* c_sequence_of_scores);
    
    void delete_information_of_fake_alignment(void);
    
    // **********************************************************************
    
    Score get_log_odds_emission_score(const int reference_state,
				      const int state,
				      const int linear_index) const;

    int   score_user_defined_state_path(const Sequence* x,
					const Sequence* y,
					const int  abs_start_x,
					const int  abs_start_y,
					const int  number_of_states,
					const int* states,
					// output
					Sequence* sub_x,
					Sequence* sub_y);
   
    void   reset_variables_for_state_path(void);

    int    check_model_and_mirrored_model(const Hmm* mirror,
					  const Sequence *x, const int x_start, const int x_end,
					  const Sequence *y, const int y_start, const int y_end);

    int    plain_viterbi(const Hmm *s, // unmirrored: delete and use this
			 const Sequence *x, 
			 const Sequence *y);


#ifndef _INTEL

  int viterbi_tube(const Sequence &x, const Sequence &y, 
		   const vector<int> start_point, const vector<int> end_point, 
		   Tube<int> tube,
		   StatePath &local_state_path) const;
		   	  
  int viterbi_tube_strip(const Sequence &x, const Sequence &y, 
			 const vector<int> start_point, const int End_x, Tube<int> tube,
			 ViterbiElement ****viterbi_strip) const;

  int viterbi_tube_strip_backwards(const Sequence &x, const Sequence &y, 
				   const int Input_Start_x, const vector<int> end_point,
				   Hmm *unmirrored_s, 
				   Tube<int> tube, 
				   ViterbiElement ****viterbi_strip) const;

  int combine_strips(ViterbiElement ****f_strip, ViterbiElement ****b_strip, 
		     const Sequence &x, const Sequence &y,
		     const int tube_width, const int tube_height,
		     const int x_middle,
		     vector<int> &left_point, vector<int> &right_point) const;

  int GetTrainingParameters(const char* cIn_dir,
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
      );

  int get_sequence_and_tube_for_training(FILE* fIn_sequence,
					 const char* cIn_annotation,
					 const char* tube_file,
					 const int radius,
					 model_parameters* const MP,
					 Sequence * const sX, 
					 Sequence * const sY,
					 Tube<int> * const tTube,
					 vector<int>* Start_point,
					 vector<int>* End_point );
  
  int Viterbi_train_emission(const Sequence &sX, 
			     const Sequence &sY,
			     const long number_of_emissions_being_train,
			     long* emission_index,
			     int** from_list,		      
			     long** final_count,
			     Tube<int> tube,
			     const vector<int> Start_point,
			     const vector<int> End_point);

  int Posterior_train_emission(const Sequence &sX, 
			       const Sequence &sY,
			       const long number_of_emissions_being_train,
			       const int number_of_sample_paths,
			       long* emission_index,				       
			       int** from_list,		      
			       long** final_count,
			       Tube<int> tube,
			       const vector<int> Start_point,
			       const vector<int> End_point);
  
  int Posterior_train_transition(const Sequence &sX, 
				 const Sequence &sY,
				 const int number_of_transitions_being_train,
				 const int number_of_sample_paths,
				 int** num_of_from_list,
				 int*** from_list,
				 long** final_count,
				 Tube<int> tube,
				 const vector<int> Start_point,
				 const vector<int> End_point);

      
  int Viterbi_train_transition(const Sequence &sX, 
			       const Sequence &sY,
			       const int number_of_transitions_being_train,
			       int** num_of_from_list,
			       int*** from_list,
			       long** final_count,
			       Tube<int> tube,
			       const vector<int> Start_point,
			       const vector<int> End_point);

  int BaumWelch_train_transition(const Sequence &sX,
				 const Sequence &sY,
				 const int number_of_transitions_being_train,
				 int** num_of_from_list,
				 int*** from_list,
				 double** final_score,
				 double& forward_score,
				 Tube<int> tube,
				 const vector<int> Start_point, 
				 const vector<int> End_point);

  int BaumWelch_train_emission(const Sequence &sX, 
			       const Sequence &sY,
			       const int number_of_emissions_being_train,
			       int* emission_index,
			       int** from_list,
			       double** final_score,
			       double& forward_score,
			       Tube<int> tube,
			       const vector<int> Start_point,
			       const vector<int> End_point);

  int Viterbi_train_TTP(const char* seqfile,
			const char* seqannfile,
			const char* name_input_tube_file,
			const long int max_volume,  
			const int radius,
			int& SeqCount,
			long*** prev_TTP_count,
			long*** cur_TTP_count,		
			model_parameters* const MP,
			TransitionProb* TP);

  int Posterior_train_TTP(const char* seqfile,
			  const char* seqannfile,
			  const char* name_input_tube_file,
			  const long int max_volume,  
			  const int radius,
			  int& SeqCount,
			  model_parameters* const MP,
			  TransitionProb* TP,
			  const int SamplePaths);

  int BaumWelch_train_TTP(const char* seqfile,
		const char* seqannfile,
		const char* name_input_tube_file,
		const long int max_volume,
		const int radius,
		bool& get_forward_score,
		double& cur_forward_score,
		int& SeqCount,
		model_parameters* const MP,
		TransitionProb* TP);

  int Viterbi_train_FTP(const char* seqfile,
			const char* seqannfile,
			const char* name_input_tube_file,
			const long int max_volume, 
			const int radius,
			int& SeqCount,
			long*** prev_GTP_count,
			long*** cur_GTP_count, 
			model_parameters* const MP,
			TransitionProb* TP);

  int Posterior_train_FTP(const char* seqfile,
			  const char* seqannfile,
			  const char* name_input_tube_file,
			  const long int max_volume, 
			  const int radius,
			  int& SeqCount,
			  model_parameters* const MP,
			  TransitionProb* TP,
			  const int SamplePaths);
  
  int BaumWelch_train_FTP(const char* seqfile,
		const char* seqannfile,
		const char* name_input_tube_file,
		const long int max_volume,
		const int radius,
		bool& get_forward_score,
		double& cur_forward_score,
		int& SeqCount,
		model_parameters* const MP,
		TransitionProb* TP);
  
  int Viterbi_train_EP(const char* seqfile,
		       const char* seqannfile,
		       const char* name_input_tube_file,
		       const long int max_volume,
		       const int radius,
		       int& SeqCount,
		       long**** prev_EP_count,
		       long**** cur_EP_count,
		       model_parameters* const MP,
		       EmissionProb* EP);

  int Posterior_train_EP(const char* seqfile,
			 const char* seqannfile,
			 const char* name_input_tube_file,
			 const long int max_volume,
			 const int radius,
			 int& SeqCount,
			 model_parameters* const MP,
			 EmissionProb* EP,
			 const int SamplePaths);

  int BaumWelch_train_EP(const char* seqfile,
	       const char* seqannfile,
	       const char* name_input_tube_file,
	       const long int max_volume,
	       const int radius,
	       bool& get_forward_score,
	       double& cur_forward_score,
	       int& SeqCount,
	       model_parameters* const MP,
	       EmissionProb* EP);

  int PosteriorTraining(const char* XMLfile,
		        TransitionProb* TP,
		        EmissionProb* EP,
		        model_parameters* const MP,
		        const char* seqfile,
		        const char* seqannfile,
		        const char* name_input_tube_file,
		        const long long int max_volume,
		        const int radius,
		        const int Maxiter,
		        const int SamplePaths);

  int ViterbiTraining(const char* XMLfile,
		      TransitionProb* TP,
		      EmissionProb* EP,
		      model_parameters* const MP,
		      const char* seqfile,
		      const char* seqannfile,
		      const char* name_input_tube_file,
		      const long long int max_volume,
		      const int radius,
		      const int Maxiter);
  
  int BaumWelchTraining(const char* XMLfile,
			TransitionProb* TP,
			EmissionProb* EP,
			model_parameters* const MP,
			const char* seqfile,
			const char* seqannfile,
			const char* name_input_tube_file,
			const long long int max_volume,
			const int radius,
			const int Maxiter,
			const double threshold);

  int Hirschberg_tube(const Hmm &mirror,
		      const Sequence &x, const Sequence &y,
		      const unsigned long long int max_volume, 
		      Tube<int> &tube);

  int Hirschberg_tube_internal(const Hmm &mirror,
			       const Sequence &x, const Sequence &y,
			       const vector<int> Start_point, const vector<int> End_point,	
			       const int max_deltax,
			       const unsigned long long int max_volume, 
			       Tube<int> &tube);

#endif // #ifndef _INTEL

};

#include "transitionprob.h"
#include "emissionprob.h"

#endif




