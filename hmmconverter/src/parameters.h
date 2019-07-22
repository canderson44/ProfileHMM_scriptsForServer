/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
*/

#ifndef _parameters_h
#define _parameters_h

typedef double Prob;
typedef double Score;

extern const int MaxTokens;
extern const int Max_number_of_items;

extern const int Max_word_length;
extern const int Max_line_length;
extern const int Max_number_of_gff_lines;
extern const int Max_number_of_gtf_lines;
extern const int Max_number_of_sequences_in_fasta_file;

extern const double Base;
extern const Score Logbase;
extern const Score Logzero;  // something close to -infinity for a computer

extern const double Max_deviation;
extern const double dPi;

extern const char* convert_int_to_sequence_type[3];

extern const int Fasta_line_length;

extern char*  executable_path;
extern char*  temporary_path;

#endif
