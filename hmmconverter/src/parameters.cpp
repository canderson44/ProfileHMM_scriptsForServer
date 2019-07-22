/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
*/

#include <math.h>
#include <stdlib.h>
#include "parameters.h"

const int MaxTokens = 10;
const int Max_number_of_items = 100;

const int Max_number_of_gff_lines               = 400000;
const int Max_number_of_gtf_lines               =  10000;
const int Max_number_of_sequences_in_fasta_file =  20000;

const int Max_word_length = 1000;
const int Max_line_length = 10000;

const double Base=2.718281828459;   // natural log
const Score Logzero= - pow(10,20);  

const double Max_deviation = powf(10,-4);
const double dPi          = 3.141592654;

const char* convert_int_to_sequence_type[3]={"Nosequence","SequenceX","SequenceY"};

const int Fasta_line_length = 60;

char*  executable_path = NULL;
char*  temporary_path  = NULL;

