/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare info-class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/info.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _info_h
#define _info_h

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
#include <stdio.h>  
#include "multi_array_2.h"
#include "define.h"
#include "sequence.h"

class Info
{
 private:

    int  info_number_of_sources;
    int* info_number_of_lines;
  
    int  info_number_of_annotation_labels;
    int* info_number_of_each_annotation_label;

    char**        info_source_names;	    
    char***       info_seq_names;                  

    bool****      info_annotation_labels;  
    double****    info_scores;

    int**         start_positions;        
    int**         end_positions;        
    
    char*     output;
    
 public:
    
    // constructors and destructor
    
    Info();
    Info(const int n_of_sources,
	 int* const  n_of_lines,
	 model_parameters* const MP,
	 int& check
	);
    Info(const Info &s);

    ~Info();
    
    // access functions

    int set_info_number_of_sources(const int n_of_sources);
    int set_info_number_of_lines(const int source, const int n_of_lines);
    int set_info_seq_names       (const int source, const int line, const char* seq_name);
    int set_info_source_names    (const int index, const char* source_name);
    int set_info_number_of_annotation_labels(const int n_of_annotation_labels);
    int set_info_number_of_each_annotation_label(const int type, const int n_of_label);
    int set_info_annotation_labels (const int source, const int line, const int type, const int label, const bool on);
    int set_info_scores          (const int source, const int line, const int type, const int label,  const double s);

    int set_start_positions (const int source, const int line, const int start_pos);
    int set_end_positions   (const int source, const int line, const int end_pos); 
  
    int set_output(const char* o);

    int      get_info_number_of_sources(void) const;
    int      get_info_number_of_lines(const int source) const; 

    char*    get_info_seq_names       (const int source, const int line) const;
    char*    get_info_source_names    (const int source) const;
    
    int      get_info_number_of_annotation_labels() const;
    int      get_info_number_of_each_annotation_label(const int type) const;
    bool     get_info_annotation_labels(const int source, const int line, const int type, const int label) const;
    double   get_info_scores           (const int source, const int line, const int type, const int label) const; 

    int      get_start_positions (const int source, const int line) const;      // returns -1 if errors on input
    int      get_end_positions   (const int source, const int line) const;      // returns -1 if errors on input 
 
    char*    get_output (void) const;


    // member functions
    
    int  remove_source(const int source);
    int  remove_line (const int source, const int line);
    int  preselect_lines(const int type, const int number_of_labels, int* const preselect_labels);
    int  preselect_lines(const int number_of_labels, char** const preselect_labels, model_parameters* const MP);    
    
    int  get_boundary_positions(// output
	int**  const number_of_boundary_positions,
	int*** const boundary_positions) const;

    int  get_details_for_position(// input
	const Sequence* const seq,
	const int             seq_position, //abs coordinates 
	model_parameters* const MP,
	// output
	bool****              annotation_label,
	double****            score
	) const;
   
    int get_details_for_next_label(// input
	const Sequence* const seq,				     
	const int             direction,
	const int             seq_position,
	bool***   const       seq_annotation_label,
	double*** const       seq_score,
	model_parameters* const  MP,
	// output
	int*                  position,
	bool****              annotation_label,
	double****            score
	) const;

    int get_annotation_of_sequence(// input 
	const Sequence* const seq,
	model_parameters* const MP,
	// output
	bool*****              seq_annotation_labels,
	double*****            seq_scores
	) const;

    int  get_length_of_output(void) const;
    int  get_line_of_output(// input
	const int start,
	const int length_of_line,
	// output
	char* line) const;

    int print_source (const int source,
		      model_parameters* const MP,
		      std::ostream &o) const;

    int print_line  (const int source, 
		     const int line,
		     model_parameters* const MP,
		     std::ostream &o) const;
    
    void print (model_parameters* const MP, std::ostream &o) const;
    
    // operators

    Info & operator = (const Info &p);
};

#endif
