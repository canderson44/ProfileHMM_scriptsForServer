 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#ifndef _evaluation_h
#define _evaluation_h

/*
  RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/evaluation.h,v 1.3 2008/12/14 10:39:23 natural Exp $
*/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include "multi_array_2.h"
#include "define.h"
#include "sequence.h"
#include "hmm_state.h"
#include "info.h"

// note: - none of the functions below assumes that their input (gtf or
//         fasta header lines) be ordered
//       - of those functions whose name does not explicitly state that
//         it orders their output only function
//         apply_restriction_set_to_gtf_lines returns ordered ouput

int get_details_for_position(// input
    const char*   const seq_ac,
    const int           seq_position,
    model_parameters* const MP,
    //
    const int           number_of_sources,         
    int*          const number_of_lines,
    char***       const seq_names,
    char**        const source_names,
    bool****      const annotation_labels,
    double****    const scores,
    int**         const start_positions,
    int**         const end_positions,
    // output
    bool****      const annotation_label,
    double****    const score,
    char****      const other_label
    );

int preselect_lines(// input
    const int             number_of_labels,
    const int             label_type,
    int* const            preselect_labels,
    model_parameters* const  MP,
    // input and output
    int            number_of_sources,
    int*     const number_of_lines,
    char***  const seq_names,
    char**   const source_names,
    bool**** const annotation_labels,
    int**    const start_positions,
    int**    const end_positions
    ); // optional

int preselect_lines(// input
    const int             number_of_labels,
    char** const          preselect_labels,
    model_parameters* const  MP,
    // input and output
    int            number_of_sources,
    int*     const number_of_lines,
    char***  const seq_names,
    char**   const source_names,
    bool**** const annotation_labels,
    int**    const start_positions,
    int**    const end_positions  
    ); // optional

int read_special_file_for_sequence(// input         // checked
    const char*     const file_name,
    const char*     const seq_ac,
    const int             seq_start_position,
    const int             seq_end_position,
    const int             position_offset,
    model_parameters*  const MP,
    // output
    int*         const number_of_allocated_lines,
    int*         const number_of_sources,
    int**        const number_of_lines,
    char****     const seq_names,
    char***      const source_names,
    bool*****    const annotation_labels,
    double*****  const scores,
    int***       const start_positions,
    int***       const end_positions
    );

int initialize_labels(
    // input
    const int            number_of_sources,
    int*           const number_of_lines,
    model_parameters* const MP,
    // output
    char****     const seq_names,
    char***      const source_names,
    bool*****    const annotation_labels,
    double*****  const scores,
    int***       const start_positions,
    int***       const end_positions
    );

int initialize_labels(
    // input
    const int            number_of_lines,
    model_parameters* const MP,
    // output
    char***      const seq_names,
    bool****     const annotation_labels,
    double****   const scores,
    int**        const start_positions,
    int**        const end_positions
    );

int release_memory_for_labels(
    // input
    int*           const number_of_sources,
    int**          const number_of_lines,
    model_parameters* const MP,
    // output
    char****     const seq_names,
    char***      const source_names,
    bool*****    const annotation_labels,
    double*****  const scores,
    int***       const start_positions,
    int***       const end_positions
    );

int release_memory_for_labels(
    // input
    int*           const number_of_lines,
    model_parameters* const MP,
    //output
    char***      const seq_names,
    bool****     const annotation_labels,
    double****   const scores,
    int**        const start_positions,
    int**        const end_positions
    );

int sort_lines_by_sources(//input
    const int            number_of_line,
    model_parameters* const MP,
    char***        const line_seq_names,
    char***        const line_source_names,
    bool****       const line_annotation_labels,
    double****     const line_scores,
    int**          const line_start_positions,
    int**          const line_end_positions,
    // output
    int*                 number_of_sources,
    int**                number_of_lines);

int read_special_file_for_double_sequences(// input        // checked
    const char*     const gtf_file_name,
    const char*     const seq_ac_1,
    const int             seq_start_position_1,
    const int             seq_end_position_1,
    const int             seq_orientation_1,
    const char*     const seq_ac_2,
    const int             seq_start_position_2,
    const int             seq_end_position_2,
    const int             seq_orientation_2,
    const int             feature_label_index,
    model_parameters*  const MP,
    // output
    int*       const number_of_allocated_lines,
    int*       const number_of_lines,
    // sequence 1
    char***    const seq_names_1,
    char***    const source_names_1,
    char****   const annotation_labels_1,
    char****   const other_labels_1,
    int**      const start_positions_1,
    int**      const end_positions_1,
    int**      const strands_1,
    // sequence 2
    char***     const seq_names_2,
    char***     const source_names_2,
    char****    const annotation_labels_2,
    char****    const other_labels_2,
    int**       const start_positions_2,
    int**       const end_positions_2,
    int**       const strands_2
    );

int get_positions_of_feature_start_and_end(// input // checked
    const Sequence* x,
    const char* const file_name,
    const int feature_label_index,
    const int feature_of_interest,
    model_parameters* const MP,
    // output
    int*   output_number_of_sources,
    int**  number_of_start_flags,
    int**  number_of_end_flags,
    int*** flags_start, // absolute start_positions of feature
    int*** flags_end);   // absolute end_positions of feature

int get_feature_label(// input
    const int              position,
    model_parameters* const   MP,
    // interval
    char**     const seq_name,
    const int*       const start_position,
    const int*       const end_position,
    // 
    const int**      const start_line,
    const int**      const end_line,
    // lines
    const int*  const number_of_sources,
    int**       const number_of_lines,
    char****    const seq_names,        
    char***     const source_names,     
    bool*****   const annotation_labels,
    double***** const scores,
    int***      const start_positions,  
    int***      const end_positions,    
    // output
    bool****     const feature,
    double****   const score);
   
int read_special_file(// input                     // checked
    FILE* file_ptr,
    const char* const file_name,
    int               max_number_of_lines,
    model_parameters* const MP,
    // output
    int*        number_of_sources,
    int**       number_of_lines,
    char****    seq_names,
    char***     source_names,
    bool*****   annotation_labels,
    double***** scores,
    int***      start_positions,
    int***      end_positions
    );

int apply_restriction_set_to_lines(// input             // purify checked 
    const int max_number_of_sources,
    const int max_number_of_lines,
    model_parameters* const MP,
    // restriction set
    int*      const res_number_of_sources,
    int**     const res_number_of_lines,
    char****  const res_seq_names,        
    int***    const res_start_positions,  
    int***    const res_end_positions,    
    // gtf-set
    // variables are used as input, are modified and then used as output
    int*        const number_of_sources,
    int**       const number_of_lines,
    char****    const seq_names,        
    char***     const source_names,     
    bool*****   const annotation_labels,
    double***** const scores,
    int***      const start_positions,  
    int***      const end_positions
    );

#endif				    
				  
