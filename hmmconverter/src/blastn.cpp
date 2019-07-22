 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: run blatn and blastx, and obtain local alignment for reduced space searching
 */

#include <stdio.h>
#include <math.h>
#include "blastn.h"
#include "multi_array_2.h"
#include "define.h"
#include "sequence.h"

/*
  RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/blastn.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
*/

Tube<int> convert_matches_to_tube(// input
				  const int    length_x,
				  const int    length_y,
				  const int    n_of_matches,
				  const Match* matches, 
				  const int    min_radius   // minimum value of tube radius				  
				  )  
{

    if ((n_of_matches > 0) && (matches != NULL)) {

	int i, j; 
	
	int length_covered_x = 0;
	int jump             = 0;

	int max_gaps = 1;
	int min_overlap = 1;
	
	if (max_gaps == 1) {
	
	    jump = 2;
	    
	    for (i=0; i<n_of_matches; i++) {
		length_covered_x += matches[i].End_x() - matches[i].Start_x() + 1;
	    }
	}
	else {
	
	    jump             = 1;
	    length_covered_x = matches[n_of_matches-1].End_x() - matches[0].Start_x() + 1;
	}
	
	// get lower and upper points
	
	int* upper_x = new int[n_of_matches*2];
	int* lower_x = new int[n_of_matches*2];
	int* upper_y = new int[n_of_matches*2];
	int* lower_y = new int[n_of_matches*2];
	
	for (i=0; i<n_of_matches; i++) {

	    lower_x[2*i]   = matches[i].Start_x();           
	    lower_y[2*i]   = matches[i].Start_y();  
	    lower_x[2*i+1] = matches[i].End_x();             
	    lower_y[2*i+1] = matches[i].End_y();    
	    
	    upper_x[2*i]   = lower_x[2*i];
	    upper_y[2*i]   = lower_y[2*i];
	    upper_x[2*i+1] = lower_x[2*i+1];
	    upper_y[2*i+1] = lower_y[2*i+1];
	}
    
	float real_rise  = 0.0;
	float max_rise   = 0.0;
	int   radius     = 0;
	int   next       = 1;
	
	if (max_gaps == 1) { // if no tube between matches is to be used => consider only matches
	    next = 2;
	}
	
	for (i=0; i<(2*n_of_matches)-1; i+=next) {
	    
	    real_rise = ((float)(lower_y[i+1] - lower_y[i]))/((float)(lower_x[i+1] - lower_x[i]));
	    if (real_rise > max_rise) {max_rise = real_rise;}

	}

	radius = (int)(ceil(max_rise/2.))-1;
	
	if (radius < min_radius) {
	    
	    radius = min_radius;

	}
	if ((min_overlap > 0) && (((2*radius)-1) < min_overlap)) {
	    
	    radius = (int)(ceil(((float)(min_overlap)+1.0)/2.0));

	}
	
	for (i=0; i<n_of_matches; i++) {
	    
	    lower_y[2*i]   -= radius;  
	    lower_y[2*i+1] -= radius;  
	    
	    upper_y[2*i]   += radius;
	    upper_y[2*i+1] += radius;
	          
	}

	// get tube in between lower and upper points
	
	int offset_x       = 0;
	int rise           = 0;
	int lower_offset_y = 0;
	int upper_offset_y = 0;
	
	int covered_x      = 0;
	
	vector<vector<int> > upper(length_covered_x, 2);
	vector<vector<int> > lower(length_covered_x, 2);
	
	int count_x = 0;
	
	int local_x       = 0;
	int local_y_lower = 0;
	int local_y_upper = 0;
	
	float real_rise_l = 0.0; 
	float real_rise_u = 0.0; 

	for (i=0; i<(2*n_of_matches)-1; i+=jump) {
	    
	    offset_x       = lower_x[i]; // = upper_x[i]
	    real_rise_l    = ((float)(lower_y[i+1] - lower_y[i]))/((float)(lower_x[i+1] - lower_x[i]));
	    real_rise_u    = ((float)(upper_y[i+1] - upper_y[i]))/((float)(upper_x[i+1] - upper_x[i]));
	    lower_offset_y = lower_y[i];
	    upper_offset_y = upper_y[i];
	    
	    covered_x      = lower_x[i+1] - lower_x[i];
	    
	    if ((i == (2*n_of_matches)-2) || (jump == 2)) {covered_x += 1;}
	    
	    for (j=0; j<covered_x; j++) {
		
		local_x       = offset_x+j;
		local_y_lower = (int)ceilf((float)(lower_offset_y) + (float)(j) * real_rise_l);
		local_y_upper = (int)ceilf((float)(upper_offset_y) + (float)(j) * real_rise_u);
		
		lower[count_x][0] = local_x;                       // x
		lower[count_x][1] = max(0, local_y_lower);         // y (check boundaries
		
		upper[count_x][0] = local_x;                       // x
		upper[count_x][1] = min(local_y_upper, length_y);  // y (check boundaries)	
		count_x++;
	    }
	}
       	
	// delete memory 
	
	if (upper_x) delete [] upper_x;
	if (upper_y) delete [] upper_y;
	if (lower_x) delete [] lower_x;
	if (lower_y) delete [] lower_y;
	
	upper_x = NULL;
	upper_y = NULL;
	lower_x = NULL;
	lower_y = NULL;
	
	Tube<int> final_tube = convert_boundaries_to_tube(length_x, length_y, upper, lower);
		
	return(final_tube);
    }
    else {
	
	Tube<int> empty_tube;
	return(empty_tube);
    }
}

int get_blastn_results(// input 
		       Sequence* sequence_x,
		       Sequence* sequence_y,
		       // output
		       int*    n_of_matches,
		       Match** matches,
		       model_parameters* const MP)
{
    // initialise output variables
    
    (*n_of_matches) = 0;
    if ((*matches)) delete [] (*matches);
    (*matches) = NULL;

    // call function
    
    int    n_matches         = 0;
    Score* scores_of_matches = NULL;
    int*   x_start_positions = NULL;
    int*   x_end_positions   = NULL;
    int*   y_start_positions = NULL;
    int*   y_end_positions   = NULL;
    
    int check = get_blastn_results(sequence_x, sequence_y,
				   &n_matches,
				   &scores_of_matches,
				   &x_start_positions,
				   &x_end_positions,
				   &y_start_positions,
				   &y_end_positions,
				   MP);
    
    // get matches in desired format
    
    if ((check == 0) && (n_matches > 0)) {
	
	// select matches for which start <= end positions
	
	int i;
	
	int count_matches = 0;
	
	for (i=0; i<n_matches; i++) {
	    if ((x_start_positions[i] <= x_end_positions[i]) &&
		(y_start_positions[i] <= y_end_positions[i])) {
		
		count_matches++;
	    }
	}
	
	(*n_of_matches) = count_matches;
	(*matches)      = new Match[count_matches];
	
	int count = 0;
	
	for (i=0; i<n_matches; i++) {
	    if ((x_start_positions[i] <= x_end_positions[i]) &&
		(y_start_positions[i] <= y_end_positions[i])) {
		
		Match next_match(x_start_positions[i], y_start_positions[i],
				 x_end_positions[i], y_end_positions[i],
				 scores_of_matches[i]);
		(*matches)[count] = next_match;
		
		count++;
	    }
	}
    }
    
    // release memory 
    
    if (scores_of_matches) delete [] scores_of_matches; 
    scores_of_matches = NULL;
    
    if (x_start_positions) delete [] x_start_positions; 
    x_start_positions = NULL;
    
    if (x_end_positions) delete [] x_end_positions; 
    x_end_positions = NULL;
    
    if (y_start_positions) delete [] y_start_positions; 
    y_start_positions = NULL;

    if (y_end_positions) delete [] y_end_positions; 
    y_end_positions = NULL;
    
    return(check);
}

int get_tblastx_results(// input 
		       Sequence* sequence_x,
		       Sequence* sequence_y,
		       // output
		       int*    n_of_matches,
		       Match** matches,
		       model_parameters* const MP)
{
    // initialise output variables

    (*n_of_matches) = 0;
    if ((*matches)) delete [] (*matches);
    (*matches) = NULL;
    
    // call function
    
    int    n_matches         = 0;
    Score* scores_of_matches = NULL;
    int*   x_start_positions = NULL;
    int*   x_end_positions   = NULL;
    int*   y_start_positions = NULL;
    int*   y_end_positions   = NULL;
    
    int check = get_tblastx_results(sequence_x, sequence_y,
				    &n_matches,
				    &scores_of_matches,
				    &x_start_positions,
				    &x_end_positions,
				    &y_start_positions,
				    &y_end_positions,
				    MP);

    // get matches in desired format
    
    if ((check == 0) && (n_matches > 0)) {
	
	// select matches for which start <= end positions
	
	int i;
	
	int count_matches = 0;
	
	for (i=0; i<n_matches; i++) {
	    if ((x_start_positions[i] <= x_end_positions[i]) &&
		(y_start_positions[i] <= y_end_positions[i])) {
		
		count_matches++;
	    }
	}
	
	(*n_of_matches) = count_matches;
	(*matches)      = new Match[count_matches];
	
	int count = 0;
	
	for (i=0; i<n_matches; i++) {
	    if ((x_start_positions[i] <= x_end_positions[i]) &&
		(y_start_positions[i] <= y_end_positions[i])) {
		
		Match next_match(x_start_positions[i], y_start_positions[i],
				 x_end_positions[i], y_end_positions[i],
				 scores_of_matches[i]);
		(*matches)[count] = next_match;
		
		count++;
	    }
	}
    }
    
    // release memory 
    
    if (scores_of_matches) delete [] scores_of_matches; 
    scores_of_matches = NULL;
    
    if (x_start_positions) delete [] x_start_positions; 
    x_start_positions = NULL;
    
    if (x_end_positions) delete [] x_end_positions; 
    x_end_positions = NULL;

    if (y_start_positions) delete [] y_start_positions; 
    y_start_positions = NULL;
    
    if (y_end_positions) delete [] y_end_positions; 
    y_end_positions = NULL;
    
    return(check);
}

int find_maximal_subset_of_matches(//input
				   const int    n_of_matches,
				   const Match* matches,
				   // output
				   int*    new_n_of_matches,
				   Match** new_matches) {

    // NOTE: the scores of matches must be positive for this procedure to work (finds maximum)
    
    
    if ((matches != NULL) && (n_of_matches>0)) {
	
	int from, to, i, max_from;
	double max_from_score, remaining_score;
	
	double* viterbi = new double[n_of_matches+2];
	int*    memory  = new int[n_of_matches+2];
	
	for (i=0; i<n_of_matches+2; i++) {
	    
	    viterbi[i] = 0.0;
	    memory[i]  = 0;
	}
	
	// initialisation (start -> matches)
	
	for (i=1; i<n_of_matches+1; i++) {
      
	    viterbi[i] = matches[i-1].Score();

	}
	
	// recursion (among matches)	
	
	for (to=1; to<n_of_matches+1; to++) {
	    
	    max_from       = 0;
	    max_from_score = viterbi[max_from];
	    
	    for (from=1; from<n_of_matches+1; from++) {		
		
		if ((matches[from-1] < matches[to-1]) &&
		    (viterbi[from] > max_from_score)) {

		    max_from       = from;
		    max_from_score = viterbi[from];
		    

		}
	    }
	    
	    viterbi[to] += max_from_score;
	    memory[to]   = max_from;

	}
	
	// termination (matches -> end)
	
	max_from       = 1;
	max_from_score = viterbi[1];;

	
	for (from=2; from<n_of_matches+1; from++) {
	   
	    
	    if (viterbi[from] > max_from_score) {
		
		max_from       = from;
		max_from_score = viterbi[from];
		
	    }
	}
	
	viterbi[n_of_matches+1] = max_from_score;
	memory[n_of_matches+1]  = max_from;
	
	// traceback
	
	(*new_n_of_matches) = 0;
	if ((*new_matches)) delete [] (*new_matches);
	(*new_matches) = NULL;
	
	Match* new_matches_reverse = new Match[n_of_matches];
	int    new_n_matches       = 0;
	
	remaining_score = viterbi[n_of_matches+1];
	from = memory[n_of_matches+1];
       
	
	while (from > 0) {
	    
	    new_matches_reverse[new_n_matches] = matches[from-1];
	    new_n_matches++;

	    from = memory[from];
	    
	}
	
	(*new_n_of_matches) = new_n_matches;
	(*new_matches)      = new Match[new_n_matches];
		
	for (i=new_n_matches-1; i>=0; i--) {
	    
	    (*new_matches)[(new_n_matches-1)-i] = new_matches_reverse[i];
	    
	}
	
	if (viterbi) delete [] viterbi;
	viterbi = NULL;
	
	if (memory) delete [] memory;
	memory = NULL;
	
	if (new_matches_reverse) delete [] new_matches_reverse;
	new_matches_reverse = NULL;
	
    }
    return(0);
}

int discard_blastn_matches_below_score_threshold(// input
						 const Score blastn_score_threshold,
						 // input and output
						 int* n_of_matches,
						 Score** scores_of_matches,
						 int**   x_start_positions,
						 int**   x_end_positions,
						 int**   y_start_positions,
						 int**   y_end_positions) // checked
{
    // note: - this function will not complain if fed zero matches

    int check = 0;
    
    if (blastn_score_threshold < 0) {
	cout << "ERROR: select_blastn_matches : blastn_score_threshold ("
	     << blastn_score_threshold << ") < 0.\n";
	check++;}
    
    if ((check == 0)                 &&
	((*n_of_matches) > 0)        &&
	(scores_of_matches != NULL)  &&
	(x_start_positions != NULL)  &&
	(x_end_positions   != NULL)  &&      
	(y_start_positions != NULL)  &&
	(y_end_positions   != NULL)  &&
	(blastn_score_threshold > 0)) {
       
	int i, j, restart = 0;

  RESTART_LOOP:


	for (i=restart; i<(*n_of_matches); i++) {
	    
	    if ((*scores_of_matches)[i] < blastn_score_threshold) {
		
		for (j=i+1; j<(*n_of_matches); j++) {
		    
		    (*x_start_positions)[j-1] = (*x_start_positions)[j];
		    (*x_end_positions)[j-1]   = (*x_end_positions)[j];
		    (*y_start_positions)[j-1] = (*y_start_positions)[j];
		    (*y_end_positions)[j-1]   = (*y_end_positions)[j];
		    (*scores_of_matches)[j-1] = (*scores_of_matches)[j];
		}
		
		(*n_of_matches)--;
		restart = i;
		goto RESTART_LOOP;
	    }
	}
	
    }
    return(check);
}

int get_blastn_results(// input 
		       Sequence* sequence_x,
		       Sequence* sequence_y,
		       // output
		       int* n_of_matches,
		       Score** scores_of_matches,
		       int** x_start_positions,
		       int** x_end_positions,
		       int** y_start_positions,
		       int** y_end_positions,
		       model_parameters* const MP)     
{
    // note : - function returns 0 if o.k. 
    //        - pointers of output-values will be NULL if not o.k.
    //        - release memory for the arrays of positions outside this function
    //        - reports if zero matches are found (see "NOTE...")
    
    int check=0;
    
    const double blast_lambda = 1.0/3.615;   // WashU blastn on fork
    const double blast_k      = 0.173;       // WashU blastn on fork
    
    // check if sequence_x and sequence_y exist
    
    if ((sequence_x == NULL) || (sequence_x->length() <1))
    {
	cout << "ERROR: get_blastn_results : sequence_x NULL.\n";
	check+=1;
    }
    if ((sequence_y == NULL) || (sequence_y->length() <1))
    {
	cout << "ERROR: get_blastn_results : sequence_y NULL.\n";
	check+=1;
    }
    
    if (check==0)
    {

	int i = 0;
      
	// fixed parameters
      
	const int   max_number_of_hsp_in_blastn = 10000; // max value for hsp per db in blastn (blastn default is 100)
	char* max_hsp_blastn = new char[Max_word_length];
	sprintf(max_hsp_blastn, "%d", max_number_of_hsp_in_blastn);

	const int   max_number_of_matches=10000; // maximal number of matches in blastn that can be dealt with here
	const int   max_line_length=500;        // maximal length of each line in the blastn results file
	const int   max_word_length=50;         // maximal length of each word in the blastn results file

	// create fasta-files for the two sequences

	char* filename_x          = new char[Max_word_length];
	char* filename_y          = new char[Max_word_length];
	char* filename_blastn     = new char[Max_word_length];
	
	char* command_line        = new char[max_line_length];
	
	// 0.) get temporary file names
	
	// buh : The memory allocated by each tempnam is deleted using unlink.
	// Using purify there is no memory leak, but purify complains
	// ("Freeing mismatched memory") about the 'filename_x delete []' etc. lines below.
	// (see *+*+*)

	strcpy(filename_x         , "XXXXXXXX");
	strcpy(filename_y         , "XXXXXXXX");
	strcpy(filename_blastn    , "XXXXXXXX");
	
	mkstemp(filename_x);
	mkstemp(filename_y);
	mkstemp(filename_blastn);
	
	// 1.) write sequences to file_x and file_y with names filename_x and filename_y
	//     and close files
	
	FILE* file_x = fopen(filename_x, "wt");
	if (sequence_x->get_ac() != NULL)
	{
	    fprintf(file_x, "%s%s%s", ">", sequence_x->get_ac(), "\n");
	}
	else
	{
	    fprintf(file_x, "%s", ">sequence_x\n");
	}
	{
	    int count = 0; 
	    for (i=0; i<sequence_x->length(); i++) {
		if (count == Fasta_line_length) {
		    fputs("\n", file_x);           
		    count = 0;                     
		}                                
		fputc(MP->get_Alphabet_name(sequence_x->letter(i)), file_x);
		count++;                         
	    }
	    fputs("\n", file_x);               
	}
	fclose(file_x);
	
	FILE* file_y = fopen(filename_y, "wt");
	if (sequence_y->get_ac() != NULL)
	{
	    fprintf(file_y, "%s%s%s", ">", sequence_y->get_ac(), "\n");
	}
	else
	{
	    fprintf(file_y, "%s", ">sequence_y\n");
	}
	{
	    int count = 0; 
	    for (i=0; i<sequence_y->length(); i++) {
		if (count == Fasta_line_length) { 
		    fputs("\n", file_y);            
		    count = 0;                      
		}                                 
		fputc(MP->get_Alphabet_name(sequence_y->letter(i)), file_y);
		count++;                          
	    }
	    fputs("\n", file_y);                
	}
	fclose(file_y);
	
	strcpy(command_line, executable_path);
	strcat(command_line, "pressdb ");
	
	strcat(command_line, filename_x);
	strcat(command_line, " >& /dev/null");

	FILE* delete_file_pressdb_x=popen(command_line, "r");
	pclose(delete_file_pressdb_x);
	
	strcpy(command_line, executable_path);
	strcat(command_line, "pressdb ");
	
	strcat(command_line, filename_y);
	strcat(command_line, " >& /dev/null"); 

	FILE* delete_file_pressdb_y=popen(command_line, "r");
	pclose(delete_file_pressdb_y);
	
	// 4.) use fasta files with names filename_x and _y to generate blastn result
	
	strcpy(command_line, executable_path);
	strcat(command_line, "blastn ");
	strcat(command_line, filename_y);
	strcat(command_line, " ");
	strcat(command_line, filename_x);

	strcat(command_line, " -top -hspmax ");
	strcat(command_line, max_hsp_blastn);
	strcat(command_line, " | ");
	strcat(command_line, executable_path);
	strcat(command_line, "MSPcrunch -d -w - | sort -nr ");
	
	strcat(command_line, " | perl -ane 'for $i (0..7) {print $F[$i], \" \";} print \"\\n\";' > ");
	
	// equal to the following command line used in a shell: 
	// perl -ane 'for $i (0..7) {print $F[$i], " ";} print "\n";'
	
	strcat(command_line, filename_blastn);

	FILE* execute_blastn_command=popen(command_line, "r");
	pclose(execute_blastn_command);

	// 5.) delete files with names filename_x and _y, filename_x, filename_x.??? (.csq,
	//     .nhd, .ntb) and the corresponding _y files;
	//     release memory
	
	unlink(filename_x);
	unlink(filename_y);
	
	strcpy(command_line, "rm -f ");
	strcat(command_line, filename_x);
	strcat(command_line, ".???");
	FILE* delete_pressdb_files_x=popen(command_line, "r");
	pclose(delete_pressdb_files_x);
	
	strcpy(command_line, "rm -f ");
	strcat(command_line, filename_y);
	strcat(command_line, ".???");
	FILE* delete_pressdb_files_y=popen(command_line, "r");
	pclose(delete_pressdb_files_y);
	
	// read file with blastn results
  
	FILE* file_with_blastn_results=fopen(filename_blastn, "r"); 
	
	// count number of matches
	
	int number_of_matches=0;
	int check_number_of_matches=0;
	
	strcpy(command_line, "less ");
	strcat(command_line, filename_blastn);
	strcat(command_line, " | wc -l");
	FILE* match_count_command=popen(command_line, "r");
	fscanf(match_count_command, "%i", &number_of_matches);
	pclose(match_count_command);

	// variables that will hold the matches generated by blastn 
	
	int* start_x_positions=NULL;
	int* end_x_positions=NULL;
	int* start_y_positions=NULL;
	int* end_y_positions=NULL;
	Score* scores=NULL;
	
	if (number_of_matches == 0) {
	    
	    cout << "NOTE: get_blastn_results : number_of_matches (" << number_of_matches
		 << ") is zero.\n";
	}
	else if ((number_of_matches > 0) && (number_of_matches < max_number_of_matches)) {
	    
	    // allocate memory for matches
	    
	    start_x_positions = new int[number_of_matches];
	    end_x_positions   = new int[number_of_matches];
	    start_y_positions = new int[number_of_matches];
	    end_y_positions   = new int[number_of_matches];
	    scores            = new Score[number_of_matches];
	    
	    // variables for each line of blastn results
	    
	    int start_x=0;                  // start of match in x sequence 
	    int end_x=0;                    // end of match in x sequence
	    int start_y=0;                  // start of match in y sequence 
	    int end_y=0;                    // end of match in y sequence
	    int score=0;                    // will hold score of match
	    Score dScore = 0.0;         // lalala
	    
	    float number=0.0;               // ? <------------------- buh
	    
	    // variables to process blastn output
	    
	    char* word_x = new char[max_word_length+1];
	    char* word_y = new char[max_word_length+1];
	    char* line   = new char[max_line_length+1];  
	    
	    char* word_format= new char[100];
	    sprintf(word_format, "%s%i%s", "%", max_word_length, "s");
	    char* line_format= new char[100];
	    sprintf(line_format, "%s%i%s%i%s", "%i%f%i%i%", max_word_length, "s%i%i%", max_word_length,"s");
	    
	    // read lines of blastn results file
	
	    int match_number             = 0;
	    double constant              = 0.0;

	    constant = log_Stirling(number_of_matches)/(double)(number_of_matches);

	    while (! feof(file_with_blastn_results)) {

		check_number_of_matches++;
		
		// --------------------------------------------------
		// read values of one line into variables
		
		fscanf(file_with_blastn_results, line_format, 
		       &score, &number, 
		       &start_x, &end_x, word_x, 
		       &start_y, &end_y, word_y);

	  // convert scores into normalized scores (* blast_lambda) and substract log(something) and add approximate constant 
          // in order to have the correct summands for calculating the sum statistics
	  // of the best scoring set of matches (see Blast book chapter 4.6 and Table 4-2 and blast
	  // manual blast.pdf page 21)

		dScore = blast_lambda * (double)(score) - log(blast_k * (double)(sequence_x->length()) * (double)(sequence_y->length()))
		    + constant;

		if (match_number < number_of_matches) {

		    start_x_positions[match_number]=start_x-1; // blast sequence positions start at 1, whereas ours start at 0
		    end_x_positions[match_number]=end_x-1; 
		    start_y_positions[match_number]=start_y-1; 
		    end_y_positions[match_number]=end_y-1; 
		    scores[match_number]=static_cast<Score>(dScore); 
		}
	  
		match_number++;
	  
		// initialise variables
	  
		start_x=0; 
		end_x=0;   
		start_y=0; 
		end_y=0;   
		score=0;   
		number=0.0;

		// --------------------------------------------------
	    } // while loop over all blastn matches
	  
	    if (check_number_of_matches != (number_of_matches+1)) {

		cout << "ERROR: get_blastn_results : number of lines read from file " 
		     << check_number_of_matches-1 << " is not equal to number_of_matches ("
		     << number_of_matches << "). This error is probably due to a change in the blastn format.\n";
		check++;
	    }
	    
	    // create output values
	    
	    *n_of_matches      = number_of_matches;
	    *scores_of_matches = scores;
	    *x_start_positions = start_x_positions;
	    *x_end_positions   = end_x_positions;
	    *y_start_positions = start_y_positions;
	    *y_end_positions   = end_y_positions;
	    
	    
	    if (line) delete [] line;
	    line = NULL;
	    if (word_x) delete [] word_x;
	    word_x = NULL;
	    if (word_y) delete [] word_y;
	    word_y = NULL;
	    if (word_y) delete [] word_y;
	    word_y = NULL;
	    if (word_format) delete [] word_format;
	    word_format = NULL;
	    if (line_format) delete [] line_format;
	    line_format = NULL;
	    
	}
	else { // if the number of matches in blastn is too large
	    
	    cout << "ERROR: get_blastn_results : number_of_matches (" << number_of_matches 
		 << ") larger than max_number_of_matches (" << max_number_of_matches 
		 << "). Increase the value of max_number_of_matches.\n";
	    check+=1;
	}

	// print results

	
	for (i=0; i<number_of_matches; i++) {
	
	    if (scores[i] < 0.0) {
		cout << "ERROR: get_blastn_results : scores[" << i << "] = " << scores[i] << " < 0.\n" << flush;
		check++;
	    }

	}


	// close and delete file with blastn results, release memory

	fclose(file_with_blastn_results);
	unlink(filename_blastn); 
	
	if (command_line)    delete [] command_line;
	command_line = NULL;
	
	if (max_hsp_blastn)       delete [] max_hsp_blastn;
	max_hsp_blastn = NULL;

	// *+*+*
	if (filename_x)           delete [] filename_x;         
	if (filename_y)           delete [] filename_y;         
	if (filename_blastn)      delete [] filename_blastn;
	// *+*+*
	
	filename_x          = NULL; 
	filename_y          = NULL; 
	filename_blastn     = NULL; 
	
	if (check == 0) {
	
	    // orientation of the two sequences is the same when blastn takes place
	    // as sequences are already reverse-complemented
	    // => sequences_in_same_orientation = 1;

	    const int sequences_in_same_orientation = 1;

	    check += select_blastn_matches(// input
		sequences_in_same_orientation,
		// input & output                 
		n_of_matches,
		scores_of_matches,
		x_start_positions,
		x_end_positions,
		y_start_positions,
		y_end_positions);
		      
	    if (check != 0) {
		cout << "ERROR: blastn : get_blastn_results : in function select_blastn_matches\n" << flush;
	    }

	}
    }
    // **********************************************************************

    if (check != 0) // set pointers of output values to NULL if something went wrong
    {
	n_of_matches      = NULL;
	scores_of_matches = NULL;
	x_start_positions = NULL;
	x_end_positions   = NULL;
	y_start_positions = NULL;
	y_end_positions   = NULL;
    }			  
    return(check);
}

int get_tblastx_results(// input 
		       Sequence* sequence_x,
		       Sequence* sequence_y,
		       // output
		       int* n_of_matches,
		       Score** scores_of_matches,
		       int** x_start_positions,
		       int** x_end_positions,
		       int** y_start_positions,
		       int** y_end_positions,
		       model_parameters* const MP)     
{
    // note : - function returns 0 if o.k. 
    //        - pointers of output-values will be NULL if not o.k.
    //        - release memory for the arrays of positions outside this function
    //        - reports if zero matches are found (see "NOTE...")

    int check=0;
  
    // check if sequence_x and sequence_y exist

    if ((sequence_x == NULL) || (sequence_x->length() <1))
    {
	cout << "ERROR: get_tblastx_results : sequence_x NULL.\n";
	check+=1;
    }
    if ((sequence_y == NULL) || (sequence_y->length() <1))
    {
	cout << "ERROR: get_tblastx_results : sequence_y NULL.\n";
	check+=1;
    }
   
    if (check==0)
    {

	int i = 0;
	
	const double blast_lambda = 1.0/2.1829;  // WashU tblastx on fork
	const double blast_k      = 0.135;       // WashU tblastx on fork
	
	// fixed parameters
      
	const int   max_number_of_hsp_in_blastn = 10000; // max value for hsp per db in blastn (blastn default is 100)
	char* max_hsp_blastn = new char[Max_word_length];
	sprintf(max_hsp_blastn, "%d", max_number_of_hsp_in_blastn);

	const int   max_number_of_matches=10000; // maximal number of matches in blastn that can be dealt with here
	const int   max_line_length=500;        // maximal length of each line in the blastn results file
	const int   max_word_length=50;         // maximal length of each word in the blastn results file


	// create fasta-files for the two sequences

	char* filename_x          = new char[Max_word_length];
	char* filename_y          = new char[Max_word_length];
	char* filename_blastn     = new char[Max_word_length];
	
	char* command_line        = new char[max_line_length];
	
	// 0.) get temporary file names
	
	// buh : The memory allocated by each tempnam is deleted using unlink.
	// Using purify there is no memory leak, but purify complains
	// ("Freeing mismatched memory") about the 'filename_x delete []' etc. lines below.
	// (see *+*+*)

	strcpy(filename_x         , "XXXXXXXX");
	strcpy(filename_y         , "XXXXXXXX");
	strcpy(filename_blastn    , "XXXXXXXX");

	mkstemp(filename_x);
	mkstemp(filename_y);
	mkstemp(filename_blastn);
	
	// 1.) write sequences to file_x and file_y with names filename_x and filename_y
	//     and close files
	
	FILE* file_x = fopen(filename_x, "wt");
	if (sequence_x->get_ac() != NULL)
	{
	    fprintf(file_x, "%s%s%s", ">", sequence_x->get_ac(), "\n");
	}
	else
	{
	    fprintf(file_x, "%s", ">sequence_x\n");
	}
	{
	    int count = 0; 
	    for (i=0; i<sequence_x->length(); i++) {
		if (count == Fasta_line_length) {
		    fputs("\n", file_x);           
		    count = 0;                     
		}                                
		fputc(MP->get_Alphabet_name(sequence_x->letter(i)), file_x);
		count++;                         
	    }
	    fputs("\n", file_x);               
	}
	fclose(file_x);
	
	FILE* file_y = fopen(filename_y, "wt");
	if (sequence_y->get_ac() != NULL)
	{
	    fprintf(file_y, "%s%s%s", ">", sequence_y->get_ac(), "\n");
	}
	else
	{
	    fprintf(file_y, "%s", ">sequence_y\n");
	}
	{
	    int count = 0; 
	    for (i=0; i<sequence_y->length(); i++) {
		if (count == Fasta_line_length) { 
		    fputs("\n", file_y);            
		    count = 0;                      
		}                                 
		fputc(MP->get_Alphabet_name(sequence_y->letter(i)), file_y);
		count++;                          
	    }
	    fputs("\n", file_y);                
	}
	fclose(file_y);
	
	strcpy(command_line, executable_path);
	strcat(command_line, "pressdb ");
	
	strcat(command_line, filename_x);
	strcat(command_line, " >& /dev/null");

	FILE* delete_file_pressdb_x=popen(command_line, "r");
	pclose(delete_file_pressdb_x);
	
	strcpy(command_line, executable_path);
	strcat(command_line, "pressdb ");
	
	strcat(command_line, filename_y);
	strcat(command_line, " >& /dev/null"); 
   
	FILE* delete_file_pressdb_y=popen(command_line, "r");
	pclose(delete_file_pressdb_y);
	
	// 4.) use fasta files with names filename_x and _y to generate blastn result
	
	strcpy(command_line, executable_path);
	strcat(command_line, "tblastx ");
	strcat(command_line, filename_y);
	strcat(command_line, " ");
	strcat(command_line, filename_x);
	strcat(command_line, " -matrix ");
	strcat(command_line, executable_path);
	strcat(command_line, "BLOSUM62");
#ifdef _NOGAP
	strcat(command_line, " -nogap ");
#endif
	strcat(command_line, " -top -hspmax ");
	strcat(command_line, max_hsp_blastn);
	strcat(command_line, " | ");
	strcat(command_line, executable_path);
	strcat(command_line, "MSPcrunch -d -w - | sort -nr ");
	
	strcat(command_line, " | perl -ane 'for $i (0..7) {print $F[$i], \" \";} print \"\\n\";' > ");
	
	// equal to the following command line used in a shell: 
	// perl -ane 'for $i (0..7) {print $F[$i], " ";} print "\n";'

	strcat(command_line, filename_blastn);

	FILE* execute_blastn_command=popen(command_line, "r");
	pclose(execute_blastn_command);
	
	// 5.) delete files with names filename_x and _y, filename_x, filename_x.??? (.csq,
	//     .nhd, .ntb) and the corresponding _y files;
	//     release memory
	
	unlink(filename_x);
	unlink(filename_y);
      
	strcpy(command_line, "rm -f ");
	strcat(command_line, filename_x);
	strcat(command_line, ".???");
	FILE* delete_pressdb_files_x=popen(command_line, "r");
	pclose(delete_pressdb_files_x);
      
	strcpy(command_line, "rm -f ");
	strcat(command_line, filename_y);
	strcat(command_line, ".???");
	FILE* delete_pressdb_files_y=popen(command_line, "r");
	pclose(delete_pressdb_files_y);
	
	// read file with blastn results
	
	FILE* file_with_blastn_results=fopen(filename_blastn, "r"); 
	
	// count number of matches
	
	int number_of_matches=0;
	int check_number_of_matches=0;
	
	strcpy(command_line, "less ");
	strcat(command_line, filename_blastn);
	strcat(command_line, " | wc -l");
	FILE* match_count_command=popen(command_line, "r");
	fscanf(match_count_command, "%i", &number_of_matches);
	pclose(match_count_command);

	// variables that will hold the matches generated by blastn 
  
	int* start_x_positions=NULL;
	int* end_x_positions=NULL;
	int* start_y_positions=NULL;
	int* end_y_positions=NULL;
	Score* scores=NULL;
	
	if (number_of_matches == 0) {
	    
	    cout << "NOTE: get_tblastx_results : number_of_matches (" << number_of_matches
		 << ") is zero.\n";
	}
	else if ((number_of_matches > 0) && (number_of_matches < max_number_of_matches)) {
	    
	    // allocate memory for matches
	    
	    start_x_positions = new int[number_of_matches];
	    end_x_positions   = new int[number_of_matches];
	    start_y_positions = new int[number_of_matches];
	    end_y_positions   = new int[number_of_matches];
	    scores            = new Score[number_of_matches];
	    
	    // variables for each line of blastn results
	    
	    int start_x=0;                  // start of match in x sequence 
	    int end_x=0;                    // end of match in x sequence
	    int start_y=0;                  // start of match in y sequence 
	    int end_y=0;                    // end of match in y sequence
	    int score=0;                    // will hold score of match
	    Score dScore = 0.0;         // lalala
	    
	    float number=0.0;               // ? <------------------- buh
	    
	    // variables to process blastn output
	    
	    char* word_x = new char[max_word_length+1];
	    char* word_y = new char[max_word_length+1];
	    char* line   = new char[max_line_length+1];  
	    
	    char* word_format= new char[100];
	    sprintf(word_format, "%s%i%s", "%", max_word_length, "s");
	    char* line_format= new char[100];
	    sprintf(line_format, "%s%i%s%i%s", "%i%f%i%i%", max_word_length, "s%i%i%", max_word_length,"s");

	    // read lines of blastn results file
	    
	    int match_number             = 0;
	    double constant              = 0.0;
	
	    constant = log_Stirling(number_of_matches)/(double)(number_of_matches);

	    while (! feof(file_with_blastn_results)) {

		check_number_of_matches++;

		// --------------------------------------------------
		// read values of one line into variables
		
		fscanf(file_with_blastn_results, line_format, 
		       &score, &number, 
		       &start_x, &end_x, word_x, 
		       &start_y, &end_y, word_y);

		// convert scores into normalized scores (* blast_lambda) and substract log(something) and add approximate constant 
		// in order to have the correct summands for calculating the sum statistics
		// of the best scoring set of matches (see Blast book chapter 4.6 and Table 4-2 and blast
		// manual blast.pdf page 21)
	  
		dScore = blast_lambda * (double)(score) - log(blast_k * (double)(sequence_x->length()) * (double)(sequence_y->length()))
		    + constant;

		if (match_number < number_of_matches) {
		    start_x_positions[match_number]=start_x-1; // blast sequence positions start at 1, whereas ours start at 0
		    end_x_positions[match_number]=end_x-1; 
		    start_y_positions[match_number]=start_y-1; 
		    end_y_positions[match_number]=end_y-1; 
		    scores[match_number]=static_cast<Score>(dScore); 
		}
	  
		match_number++;
		
		// initialise variables
		
		start_x=0; 
		end_x=0;   
		start_y=0; 
		end_y=0;   
		score=0;   
		number=0.0;
		
		// --------------------------------------------------
	    } // while loop over all blastn matches
	    
	    if (check_number_of_matches != (number_of_matches+1)) {
		
		cout << "ERROR: get_tblastx_results : number of lines read from file " 
		     << check_number_of_matches-1 << " is not equal to number_of_matches ("
		     << number_of_matches << "). This error is probably due to a change in the blastn format.\n";
		check++;
	    }
	    
	    // create output values
	    
	    *n_of_matches      = number_of_matches;
	    *scores_of_matches = scores;
	    *x_start_positions = start_x_positions;
	    *x_end_positions   = end_x_positions;
	    *y_start_positions = start_y_positions;
	    *y_end_positions   = end_y_positions;
	

	    if (line) delete [] line;
	    line = NULL;
	    if (word_x) delete [] word_x;
	    word_x = NULL;
	    if (word_y) delete [] word_y;
	    word_y = NULL;
	    if (word_y) delete [] word_y;
	    word_y = NULL;
	    if (word_format) delete [] word_format;
	    word_format = NULL;
	    if (line_format) delete [] line_format;
	    line_format = NULL;
	    
	}
	else { // if the number of matches in blastn is too large
	    
	    cout << "ERROR: get_tblastx_results : number_of_matches (" << number_of_matches 
		 << ") larger than max_number_of_matches (" << max_number_of_matches 
		 << "). Increase the value of max_number_of_matches.\n";
	    check+=1;
	}
	
	// print results


	for (i=0; i<number_of_matches; i++) {
	    
	    if (scores[i] < 0.0) {
		cout << "ERROR: get_tblastx_results : scores[" << i << "] = " << scores[i] << " < 0.\n" << flush;
		check++;
	    }

	}


	// close and delete file with blastn results, release memory

	fclose(file_with_blastn_results);
	unlink(filename_blastn); 
      
	if (command_line)    delete [] command_line;
	command_line = NULL;
	
	if (max_hsp_blastn)       delete [] max_hsp_blastn;
	max_hsp_blastn = NULL;
	
	// *+*+*
	if (filename_x)           delete [] filename_x;         
	if (filename_y)           delete [] filename_y;         
	if (filename_blastn)      delete [] filename_blastn;
	// *+*+*
	
	filename_x          = NULL; 
	filename_y          = NULL; 
	filename_blastn     = NULL; 
	
	if (check == 0) {
	    
	    // orientation of the two sequences is the same when blastn takes place
	    // as sequences are already reverse-complemented
	    // => sequences_in_same_orientation = 1;
	    
	    const int sequences_in_same_orientation = 1;
	
	    check += select_blastn_matches(// input
		sequences_in_same_orientation,
		// input & output                 
		n_of_matches,
		scores_of_matches,
		x_start_positions,
		x_end_positions,
		y_start_positions,
		y_end_positions);
	    
	    if (check != 0) {
		cout << "ERROR: blastn : get_tblastx_results : in function select_blastn_matches\n" << flush;
	    }
	}
    }
    // **********************************************************************
    
    if (check != 0) // set pointers of output values to NULL if something went wrong
    {
	n_of_matches      = NULL;
	scores_of_matches = NULL;
	x_start_positions = NULL;
	x_end_positions   = NULL;
	y_start_positions = NULL;
	y_end_positions   = NULL;
    }			  
    return(check);
}

int select_blastn_matches(// input
			  const int sequences_in_same_orientation, //
			  // 1 = yes, 0 == no
			  // input & output                 
			  int*    n_of_matches,
			  Score** scores_of_matches,
			  int**   x_start_positions,
			  int**   x_end_positions,
			  int**   y_start_positions,
			  int**   y_end_positions) // checked
{
  // note: this function will not complain if fed zero matches

    int check = 0;
    
    if ((sequences_in_same_orientation != 1) && (sequences_in_same_orientation != 0)) {
	cout << "ERROR: select_blastn_matches : sequences_in_same_orientation ("
	     << sequences_in_same_orientation << ") must be either 1 (for yes) or 0 (for no).\n";
	check++;
    }
    if ((check == 0)                &&
	((*n_of_matches) > 0)       &&
	(scores_of_matches != NULL) &&
	(x_start_positions != NULL) &&
	(x_end_positions   != NULL) &&      
	(y_start_positions != NULL) &&
	(y_end_positions   != NULL)) {
	
	int i = 0;
	int j = 0;
	int discard = 0;
	int restart = 0;
    
	int new_n_of_matches = (*n_of_matches);
	
    RESTART_LOOP:
	
	for (i=restart; i<new_n_of_matches; i++) {
	    
	    discard = 0;

	    if ((((*x_end_positions)[i] - (*x_start_positions)[i])) *
		(((*y_end_positions)[i] - (*y_start_positions)[i])) == 0) {
		
		// skip this case (pointlike match)
	    }
	    else if ((sequences_in_same_orientation == 1) &&
		     (((((*x_end_positions)[i] - (*x_start_positions)[i])) *
		       (((*y_end_positions)[i] - (*y_start_positions)[i]))) < 0)) {
		
		discard++;
	    }
	    else if ((sequences_in_same_orientation == 0) &&
		     (((((*x_end_positions)[i] - (*x_start_positions)[i])) *
		       (((*y_end_positions)[i] - (*y_start_positions)[i]))) > 0)) {
		
		discard++;
	    }
	    
	    if (discard != 0) {

		for (j=i+1; j<new_n_of_matches; j++) {
		    
		    (*scores_of_matches)[j-1] = (*scores_of_matches)[j]; 
		    (*x_start_positions)[j-1] = (*x_start_positions)[j]; 
		    (*x_end_positions)[j-1]   = (*x_end_positions)[j];   
		    (*y_start_positions)[j-1] = (*y_start_positions)[j]; 
		    (*y_end_positions)[j-1]   = (*y_end_positions)[j];   
		}
		
		new_n_of_matches--;
		
		restart = i;
		goto RESTART_LOOP;
	    }
	}
	if (check == 0) {
	    (*n_of_matches) = new_n_of_matches;
	}
    }
    if (check != 0) {
	
	if ((*scores_of_matches)) delete [] (*scores_of_matches);
	if ((*x_start_positions)) delete [] (*x_start_positions);
	if ((*x_end_positions))   delete [] (*x_end_positions);  
	if ((*y_start_positions)) delete [] (*y_start_positions);
	if ((*y_end_positions))   delete [] (*y_end_positions);  
	
	(*scores_of_matches) = NULL;
	(*x_start_positions) = NULL;
	(*x_end_positions)   = NULL;
	(*y_start_positions) = NULL;
	(*y_end_positions)   = NULL;
	
	(*n_of_matches) = 0;
	
    }
    return(check);
}



