/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define constants, typedefs and enums used in the pairhmm, sequence classes 
           
  RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/define.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _define_h
#define _define_h

#include <iostream>
#include <math.h>

#include <vector.h>
#include "parameters.h"
#include "tube.h"
#include "multi_array_2.h"
#include "tinyxml.h"


#define abs(x)    (((x) >= 0) ? (x) : -(x))
#define nint(x)   (((x) > 0) ? ((int)(x + .5)) : -((int)(.5 - x)) )
#define min(a,b)  (((a) <= (b)) ? (a) : (b))
#define max(a,b)  (((a) >= (b)) ? (a) : (b))
#define dabs(x)   (double)abs(x)
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)



//tested, return the number of charactors of the input char*
int NumOfChar(const char* const inchar); 

char* Nullstrcpy(const char* str2, int& check);

// For converting tname to int for Annotation_Label, Alphabets, TransitionProb, EmissionProb
template <typename label>
int convert_typename_to_int(label* label_type, char* type_name, int size)
{
    for(int i=0; i<size; i++){
	if(!cmp_nocase(label_type[i].setname,type_name))
	    return i;
    }    return -1;
}

// For converting labelname to int for Annoation_Label, Alphabets, TransitionProb, EmissionProb
template <typename label>
int convert_labelname_to_int(label* label_type, int type, char* name)
{
    int size = label_type[type].size;
    for(int i = 0; i<size ; i++){
	if(!cmp_nocase(label_type[type].name[i],name)){
	    return i;
	}
    }
    return -1;
}

// For converting labelname to int for Annotation_Label, Alphabets, TransitionProb, EmissionProb
template <typename label>
int convert_labelname_to_int(label* label_type, char* type, char* name, int size)
{
    int stype = -1;
    int i = 0;
    for(i=0;i<size;i++){
	if(!cmp_nocase(label_type[i].tname,type)){
	    stype = i;
	    break;
	}
    }
    if(stype == -1){
	return -1;
    }
    int ssize = label_type[stype].size;
    for(i=0;i<ssize;i++){
	if(!cmp_nocase(label_type[stype].name[i],name)){
	    return i;
	}
    }
    return -1;
}

// To break the id in the form of : ID.x.y into tokens of {ID,x,y}
int break_id_into_tokens(const char* const id, char*** tokens, int& num_of_token);

void convert_to_base_and_invert(int number, const int base, array<int>* index);

array<int> get_indices(int linear_index, int* dim, int number_of_dimensions, int& check);

int generate_random_letter(const int alphabet);

int get_substr(const int start, const int length, 
	       const char* const string, char* const substring);

Tube<int> convert_boundaries_to_tube(const int length_x,
				     const int length_y,
				     vector<vector<int> > & upper_bounds, 
				     vector<vector<int> > & lower_bounds);

void convert_to_variable_base(int number, const int number_of_bases,
			      const int* const bases, int* const index);

double log_Stirling(const int n); // Stirling approximation of log(n!), see Bronstein page 103

int bitshift(const int value);
// note: function returns positive value (>0) if successful, else returns negative value (<0)
//       the return value is the minimum integer i for which alphabet =< 2^i

double calculate_number_from_score(const Score      score_value);

Score calculate_score(const double number);

int set_executable_path(const char* const new_executable_path); 

int set_temporary_path(const char* const new_temporary_path); 

void convert_to_base(int number, const int base, array<int>* index);

int compare_state_type(int delta_x, int delta_y, int delta_x2, int delta_y2);

int get_n_chars(const int n, const char letter, char** word);

int splitstring(const char* str1, int& NoOfItems, char*** Items, char splitchar);

int cmp_nocase(const char* s1, const char* s2);

int sort_elements_of_array_by_increasing_order(const int length_array,
					       int* const array); 

int sort_elements_of_two_arrays_and_create_new_array(// input
						     const int        length_array_1,
						     const int* const array_1,
						     const int        length_array_2,
						     const int* const array_2,
						     // output
						     int*             length_array,
						     int**      const array);

int make_elements_of_array_unique(// input and output
				  int* length_array,
				  int* array);

template <typename label>
void swap(label* a, label* b)
{
    label temp;
    temp = (*a);
    (*a) = (*b);
    (*b) = temp;	
}

int swap_string(char** a, char** b);

int check_undefine(bool* const labels,
		   int   length);

bool get_tag(const char* line,
	     const char* tag);
	    

#endif
