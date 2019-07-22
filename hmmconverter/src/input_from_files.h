/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define functions which read in information on genes or emission probs
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/input_from_files.h,v 1.3 2008/08/17 07:22:37 natural Exp $
*/

#ifndef _input_from_files_h
#define _input_from_files_h

#include <stdio.h>  
#include "multi_array_2.h" 
#include "define.h"
#include "model_parameters.h"
#include "sequence.h"
#include "blastn.h"
#include "tube.h"

int get_next_sequence(FILE* sequence_file_with_pairs_of_sequences,
		      char** name_of_sequence_x,
		      char** sequence_x,
		      int* length_of_sequence_x,
		      int* start_of_sequence_x,
		      int* end_of_sequence_x,
		      model_parameters* const MP);

// note: - release memory for ac_of_gene and sequence in program which calls this function
//       - o.k. if return-value 0, else not o.k.
//       - the input sequence-file is assumed to come in the following format:
//
//                         - there are no blank lines (if there are any, they will be ignored)
//                         - the header line looks like '>seqname start-end (length:'forward' or 'reverse':number)'
//                         - each sequence directly follows its header line
//                         - the sequences come in pairs

int get_single_sequence(FILE* sequence_files, 
			Sequence* const x, 
			model_parameters* const MP);

int get_next_pair_of_sequences(FILE* sequence_file_with_pairs_of_sequences,
			       char** name_of_sequence_x,
			       char** sequence_x,
			       int* length_of_sequence_x,
			       int* start_of_sequence_x,
			       int* end_of_sequence_x,
			       char** name_of_sequence_y,
			       char** sequence_y,
			       int* length_of_sequence_y,
			       int* start_of_sequence_y,
			       int* end_of_sequence_y);

// note: - release memory for ac_of_gene and sequence in program which calls this function
//       - o.k. if return-value 0, else not o.k.
//       - the input sequence-file is assumed to come in the following format:
//
//                         - there are no blank lines (if there are any, they will be ignored)
//                         - the header line looks like '>seqname start-end (length:'forward' or 'reverse':number)'
//                         - each sequence directly follows its header line
//                         - the sequences come in pairs


int get_sequence_pair(FILE* sequence_file_with_pairs_of_sequences,
		      Sequence* const x,
		      Sequence* const y,
		      model_parameters* MP);

int get_sequence_decoding_parameters(char* const input_filename,
				     char** name_input_annotation_file,	  
				     char** name_input_sequence_file,
				     char** name_input_tube_file,				 
				     int& algorithm,
				     unsigned long long int& volume,
				     int& radius,
				     char** name_output_file);	

int get_parameter_training_parameters(char* const input_filename,
				      char** name_input_annotation_file,	  
				      char** name_input_sequence_file,
				      char** name_input_tube_file,				 
				      int& algorithm,
				      unsigned long long int& volume,
				      int& radius,
				      double& threshold,
				      int& Maxiter,
				      int& Samples,
				      char** name_transition_prob_file,
                                      char** name_emission_prob_file,
                                      char** name_XMLfile);

int get_tube_from_file(const char* filename,
		       const char* seq_x_name,
		       const char* seq_y_name,
		       const int seq_x_length,
		       const int seq_y_length,
		       Tube<int> * const tTube);

#endif
