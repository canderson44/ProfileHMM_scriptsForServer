/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define functions which read in information on genes or emission probs
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/input_from_files.cpp,v 1.3 2008/08/17 07:20:45 natural Exp $
*/

#include <stdio.h>  
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
#include <string.h>    
#include <stdlib.h>    
#include <fstream.h>
#include "input_from_files.h"

int get_next_sequence(FILE* sequence_file,
		      // output for sequence x
		      char** name_of_sequence_x,
		      char** sequence_x,
		      int* length_of_sequence_x,
		      int* start_of_sequence_x,
		      int* end_of_sequence_x,
		      model_parameters* const MP)
{
  // note: - release memory for ac_of_gene and sequence in program which calls this function
  //       - o.k. if return-value 0, else not o.k.
  //       - the input sequence-file is assumed to come in the following format:
  //
  //                         - there are no blank lines (if there are any, they will be ignored)
  //                         - the header line looks like '>seqname start-end (length:'forward' or 'reverse':number)'
  //                         - each sequence directly follows its header line
  //                         - the sequences come in pairs
  //                         - the '7' in the header_format is the length of 'reverse' and 'forward'

    int check=0;

    if (sequence_file == NULL)
    {
	cout << "ERROR: get_next_sequence: input file sequence_file is NULL.\n";
	check++;
    }

    if (!check)
    {

	int i =0;
	// parameters
	
	const int max_word_length = 100;
	const int max_line_length = 1000;
	
	long int position_in_file      = 0;
	long int last_position_in_file = 0;
	
	// allocate memory for output variables
	
	(*name_of_sequence_x) = new char[max_word_length+1];
	
	// variables for checking

	int check_header_x = 0;
	int check_seq_x    = 0;
	
	// allocate memory for internal variables
	
	(*sequence_x) = NULL;
	int check_x_l = 0;
	
	int number_1 = 0;
	int number_2 = 0;
	
	char* word             = new char[max_word_length+1];
	
	char* word_format= new char[10];
	sprintf(word_format, "%s%d%s", "%", max_word_length, "s");

	char* line = new char[max_line_length];
	
	char** Items = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	int NoOfItems = 0;

	char** NumberItems = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    NumberItems[i] = new char[Max_word_length];
	    strcpy(NumberItems[i]," ");
	}
	int NoOfNumberItems = 0;
	
	for (int data_set = 0; data_set < 1; data_set++)
	{      
	    // initialise some variables
	    
	    number_1 = 0;
	    number_2 = 0;
	    
	    // read line
	    
	    position_in_file=ftell(sequence_file);

	    // read one line if either max_line_length-1 characters were read or end of line was encountered
	    
	    fgets(line, max_line_length-1, sequence_file);
	    last_position_in_file=ftell(sequence_file);

	    // get first word in that line

	    sscanf(line, word_format, word);

	    if (strcmp(word, "//")==0)
	    {

		data_set--;
	    }
	    else
	    {

		// read header
	    
		check+=splitstring(line,NoOfItems,&Items,' ');
		if(check)
		{
		    cout<<"Error: get_next_sequence: error occured in splitstring for line:"
			<<line<<endl;
		}
		
		if(NoOfItems<2)
		{
		    cout<<"Error: get_next_sequence: NoOfItems("<<NoOfItems<<")"
			<<"< 2, invalid header line:"
			<<line<<endl;
		    check++;
		}

		if(!check)
		{
		    //get range
		    check+=splitstring(Items[1],NoOfNumberItems,&NumberItems,'-');
		    if(check)
		    {
			cout<<"Error: get_next_sequence: error occured in splitstring for range:"
			    <<Items[1]<<endl;
		    }
		
		    if(NoOfNumberItems!=2)
		    {
			cout<<"Error: get_next_sequence: NoOfNumberItems("<<NoOfNumberItems<<")"
			    <<"!= 2, invalid range:"
			<<Items[1]<<endl;
			check++;
		    }else{
			
			//get range
			
			number_1 = atoi(NumberItems[0]);
			number_2 = atoi(NumberItems[1]);
			
			if((number_1==0)&&(strcmp(NumberItems[0],"0")))
			{
			    cout<<"Error: get_next_sequence: invalid range: "
				<<Items[1]<<endl;
			    check++;
			}
			
			if((number_2==0)&&(strcmp(NumberItems[1],"0")))
			{
			    cout<<"Error: get_next_sequence: invalid range: "
				<<Items[1]<<endl;
			    check++;
			}
		    }
		}
		// check header and save information in ouput variables if header is o.k.
		
		if ((Items[0][0] != '>')       ||      // header should start with '>'
		    (number_1 < 0)        ||      // start position should 0 or larger
		    (number_2 < 0)        ||      // end position should 0 or larger
		    (number_1 > number_2) ||      // start <= end
		    (static_cast<int>(strlen(Items[0])) == 0))  // sequence name should exist
	      // or 'reverse'
		{
		    if (data_set == 0) 
		    {

			check_header_x++;
		    }
		}
		else
		{
		    if (data_set == 0) 
		    {
			for(i=1; i<strlen(Items[0]); i++)
			{
			    (*name_of_sequence_x)[i-1] = Items[0][i];
			}
			(*name_of_sequence_x)[i-1] = '\0';
			
			(*length_of_sequence_x)      = number_2-number_1+1;
			(*start_of_sequence_x)       = number_1;
			(*end_of_sequence_x)         = number_2;
			
			(*sequence_x) = new char[(*length_of_sequence_x)+1];
		    }
		}
		
		// read sequence 
		
		int  c_number = 0;
		char c = '\0';
		
		while ((c != '>') && (! feof(sequence_file)))
		{		   
		    c_number = getc(sequence_file);
		    c = static_cast<unsigned char>(c_number);
		    
		    // ignore end of lines
		    if (c == '\n') 
		    {		
		    }
		    // test for end of file
		    else if (feof(sequence_file))
		    {
		    }
		    // stop before '>'
		    else if (c == '>') 
		    {
			ungetc(c_number, sequence_file);
			break;
		    }
		    // record error if letter does not belong to alphabet
#ifndef _IGNORE_NON_ALPHABET
		    else if (convert_alphabet_to_int(MP,c) < 0) 
		    {
			cout<<"Can not recognize charactor : "<<c<<endl;
			if (data_set == 0) {
			    check_seq_x++;
			}
		    }
#endif
		    // write to strstream if everything o.k.
		    else 
		    {
			if (data_set == 0) {		
			    if (check_x_l <= (*length_of_sequence_x)) {(*sequence_x)[check_x_l] = c;}
			    check_x_l++;
			}

		    }
		}     

	    } // if this is not a comment line
	} // loop over two data sets
	
	// check headers and sequences of both data sets and write information to output variables if everything is o.k.
	
	(*sequence_x)[check_x_l] = 0;
	
	const int length_x = strlen((*sequence_x));

	if ((check_header_x == 0) &&
	    (check_seq_x == 0)    &&
	    (length_x == (*length_of_sequence_x)))
	{

	}
	else // if there is something wrong
	{

	    if ((*name_of_sequence_x)) delete [] (*name_of_sequence_x);
	    (*name_of_sequence_x) = NULL;
	    if ((*sequence_x)) delete [] (*sequence_x);
	    (*sequence_x) = NULL;
	    
	    (*length_of_sequence_x)      = 0;
	    (*start_of_sequence_x)       = 0;
	    (*end_of_sequence_x)         = 0;
	    
	    check++;
	}
	
	// release memory
	if(Items)
	{
	    for(int i=0; i<Max_number_of_items; i++)
	    {
		if(Items[i]) delete [] Items[i];
		Items[i] = NULL;
	    }
	    delete [] Items;
	}
	Items = NULL;
	
	if(NumberItems)
	{
	    for(i=0; i<Max_number_of_items; i++)
	    {
		if(NumberItems[i]) delete [] NumberItems[i];
		NumberItems[i] = NULL;
	    }
	    delete [] NumberItems;
	}
	NumberItems = NULL;

	if (word) delete [] word;
	word = NULL;
	if (word_format) delete [] word_format;
	word_format = NULL;

	if (line) delete [] line;
	line = NULL;
    }
    return(check);
}

int get_single_sequence(FILE* sequence_files, 
			Sequence* const x, 
			model_parameters* const MP)
{
    int check = 0;

    if (MP->get_Number_of_States() < 2)
    {
	cout << "ERROR: get_single_sequence: number_of_states_in_hmm ("
	     << MP->get_Number_of_States() << ") < 2.\n" << flush;
	check++;
    }

    if (sequence_files == NULL) {
	cout << "ERROR: get_single_sequence: sequence_files is NULL.\n" << flush;
	check++;
    }
    if (x == NULL) {
	cout << "ERROR: get_single_sequence: Sequence x is NULL.\n" << flush;
	check++;
    }

    if (check == 0) {
	
	char* name_of_sequence_x      = NULL;
	char* sequence_x              = NULL;
	int length_of_sequence_x      = 0;
	int start_of_sequence_x       = 0;
	int end_of_sequence_x         = 0;
	int orientation_of_sequence_x = 0;
	
	check += get_next_sequence(sequence_files,
				   // info on sequence x
				   &name_of_sequence_x, &sequence_x, &length_of_sequence_x,
				   &start_of_sequence_x, &end_of_sequence_x,MP);
	if (check) {
	    cout << "ERROR: get_single_sequence: error occurred in get_next_sequence.\n" << flush;
	}
	
	if (!check) {
	    
	    if(MP->is_SpecialEmit())
	    {
		Sequence new_x(1,
			       sequence_x,
			       MP);
		(*x) = new_x;
	    }else{

		Sequence new_x(sequence_x,MP);
		(*x) = new_x;
	    }
	    x->set_ac(name_of_sequence_x);
	    x->set_sequence_start_and_end(start_of_sequence_x, end_of_sequence_x);
	    
	}	  
	
	// release memory
	
	if (name_of_sequence_x) delete [] name_of_sequence_x;
	name_of_sequence_x = NULL;
	if (sequence_x) delete [] sequence_x;
	sequence_x = NULL;
    }
    if (check) { // if checks failed, return empty sequences
	
	Sequence new_x;
	(*x) = new_x;
    }
    return(check);
}

int get_next_pair_of_sequences(FILE* sequence_file_with_pairs_of_sequences, 
			       // output for sequence x
			       char** name_of_sequence_x,
			       char** sequence_x,
			       int* length_of_sequence_x,
			       int* start_of_sequence_x,
			       int* end_of_sequence_x,
			       // output for sequence y
			       char** name_of_sequence_y,
			       char** sequence_y,
			       int* length_of_sequence_y,
			       int* start_of_sequence_y,
			       int* end_of_sequence_y)
{
  // note: - release memory for ac_of_gene and sequence in program which calls this function
  //       - o.k. if return-value 0, else not o.k.
  //       - the input sequence-file is assumed to come in the following format:
  //
  //                         - there are no blank lines (if there are any, they will be ignored)
  //                         - the header line looks like '>seqname start-end (length:'forward' or 'reverse':number)'
  //                         - each sequence directly follows its header line
  //                         - the sequences come in pairs
  //                         - the '7' in the header_format is the length of 'reverse' and 'forward'

    int check=0;

    if (sequence_file_with_pairs_of_sequences == NULL)
    {
	cout << "ERROR: get_next_pair_of_sequences: input file sequence_file_with_pairs_of_sequences is NULL.\n";
	check++;
    }
    
    if (!check)
    {
	int i =0;
	// parameters
	
	const int max_word_length = 100;
	const int max_line_length = 1000;
      
	long int position_in_file      = 0;
	long int last_position_in_file = 0;
	
	// allocate memory for output variables
	
	(*name_of_sequence_x) = new char[max_word_length+1];
	(*name_of_sequence_y) = new char[max_word_length+1];
	
	// variables for checking
	
	int check_header_x = 0;
	int check_seq_x    = 0;
	int check_header_y = 0;
	int check_seq_y    = 0;

	// allocate memory for internal variables
	
	(*sequence_x) = NULL;
	(*sequence_y) = NULL;
	int check_l_x = 0;
	int check_l_y = 0;
	
	int number_1 = 0;
	int number_2 = 0;

	char* word             = new char[max_word_length+1];
	char* word_format= new char[10];
	sprintf(word_format, "%s%d%s", "%", max_word_length, "s");
	char* line = new char[max_line_length];

	char** Items = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	int NoOfItems = 0;

	char** NumberItems = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    NumberItems[i] = new char[Max_word_length];
	    strcpy(NumberItems[i]," ");
	}
	int NoOfNumberItems = 0;
	
	for (int data_set = 0; data_set < 2; data_set++)
	{      
	    // initialise some variables
	    
	    number_1 = 0;
	    number_2 = 0;
	  
	    // read line
	    
	    position_in_file=ftell(sequence_file_with_pairs_of_sequences);
	    
	    // read one line if either max_line_length-1 characters were read or end of line was encountered
	    
	    fgets(line, max_line_length-1, sequence_file_with_pairs_of_sequences);
	    last_position_in_file=ftell(sequence_file_with_pairs_of_sequences);
	    
	    // get first word in that line
	    
	    sscanf(line, word_format, word);
	    
	    if (strcmp(word, "//")==0)
	    {
		data_set--;
	    }
	    else
	    {
		// read header
		
		check+=splitstring(line,NoOfItems,&Items,' ');
		if(check)
		{
		    cout<<"Error: get_next_pair_of_sequence: error occured in splitstring for line:"
			<<line<<endl;
		}
		
		if(NoOfItems<2)
		{
		    cout<<"Error: get_next_pair_of_sequence: NoOfItems("<<NoOfItems<<")"
			<<"< 2, invalid header line:"
			<<line<<endl;
		    check++;
		}

		if(!check)
		{
		    //get range
		    check+=splitstring(Items[1],NoOfNumberItems,&NumberItems,'-');
		    if(check)
		    {
			cout<<"Error: get_next_pair_of_sequence: error occured in splitstring for range:"
			    <<Items[1]<<endl;
		    }
		
		    if(NoOfNumberItems!=2)
		    {
			cout<<"Error: get_next_pair_of_sequence: NoOfNumberItems("<<NoOfNumberItems<<")"
			    <<"!= 2, invalid range:"
			<<Items[1]<<endl;
			check++;
		    }else{
			
			//get range
			
			number_1 = atoi(NumberItems[0]);
			number_2 = atoi(NumberItems[1]);
			
			if((number_1==0)&&(strcmp(NumberItems[0],"0")))
			{
			    cout<<"Error: get_next_pair_of_sequence: invalid range: "
				<<Items[1]<<endl;
			    check++;
			}
			
			if((number_2==0)&&(strcmp(NumberItems[1],"0")))
			{
			    cout<<"Error: get_next_pair_of_sequence: invalid range: "
				<<Items[1]<<endl;
			    check++;
			}
		    }
		}

		// check header and save information in ouput variables if header is o.k.
		
		if ((Items[0][0] != '>')       ||                    // header should start with '>'
		    (number_1 < 0)        ||                    // start position should 0 or larger
		    (number_2 < 0)        ||                    // end position should 0 or larger
		    (number_1 > number_2) ||                    // start <= end
		    (static_cast<int>(strlen(Items[0])) == 0))  // sequence name should exist
		   
		{
		    if (data_set == 0) 
		    {

			check_header_x++;
		    }
		    else if (data_set == 1)
		    {

			check_header_y++;
		    }
		}
		else
		{
		    if (data_set == 0) 
		    {
			for(i = 1; i<strlen(Items[0]) ;i++)
			{
			    (*name_of_sequence_x)[i-1] = Items[0][i];
			}
			(*name_of_sequence_x)[i-1] = '\0';
			(*length_of_sequence_x)      = number_2-number_1+1;
			(*start_of_sequence_x)       = number_1;
			(*end_of_sequence_x)         = number_2;
			
			(*sequence_x) = new char[(*length_of_sequence_x)+1];
		    }
		    else if (data_set == 1) 
		    {
			for(i = 1; i<strlen(Items[0]) ;i++)
			{
			    (*name_of_sequence_y)[i-1] = Items[0][i];
			}
			(*name_of_sequence_y)[i-1] = '\0';
			(*length_of_sequence_y)      = number_2-number_1+1;
			(*start_of_sequence_y)       = number_1;
			(*end_of_sequence_y)         = number_2;
			
			(*sequence_y) = new char[(*length_of_sequence_y)+1];
		    }
		}
		
		// read sequence 
	      
		int  c_number = 0;
		char c = '\0';
		
		while ((c != '>') && (! feof(sequence_file_with_pairs_of_sequences)))
		{
		    c_number = getc(sequence_file_with_pairs_of_sequences);
		    c = static_cast<unsigned char>(c_number);
		    
		    // ignore end of lines
		    if (c == '\n') 
		    {}
		    // test for end of file
		    else if (feof(sequence_file_with_pairs_of_sequences))
		    {}
		    // stop before '>'
		    else if (c == '>') 
		    {
			ungetc(c_number, sequence_file_with_pairs_of_sequences);
			break;
		    }
		    // record error if letter does not belong to alphabet
#ifdef _IGNORE_NON_ALPHABET
		    else if (check_if_in_alphabet(alphabet, c) != 0) 
		    {
			if (data_set == 0) {check_seq_x++;}
			else if (data_set == 1) {check_seq_y++;}
		    }
#endif
		    // write to strstream if everything o.k.
		    else 
		    {
			if (data_set == 0) {
			  
			    if (check_l_x <= (*length_of_sequence_x)) {(*sequence_x)[check_l_x] = c;}
			    check_l_x++;
			}
			else if (data_set == 1) {
			    
			    if (check_l_y <= (*length_of_sequence_y)) {(*sequence_y)[check_l_y] = c;}
			    check_l_y++;
			}

		    }
		}      

	    } // if this is not a comment line
	} // loop over two data sets
	
	// check headers and sequences of both data sets and write information to output variables if everything is o.k.
	
	(*sequence_x)[check_l_x] = 0;
	(*sequence_y)[check_l_y] = 0;
	
	const int length_x = strlen((*sequence_x));
	const int length_y = strlen((*sequence_y));

	if ((check_header_x == 0) &&
	    (check_header_y == 0) &&
	    (check_seq_x == 0)    &&
	    (check_seq_y == 0)    &&
	    (length_x == (*length_of_sequence_x)) &&
	    (length_y == (*length_of_sequence_y)))
	{

	}
	else // if there is something wrong with one or both of the data set, reset values of output variables and report error
	{

	    if ((*name_of_sequence_x)) delete [] (*name_of_sequence_x);
	    (*name_of_sequence_x) = NULL;
	    if ((*sequence_x)) delete [] (*sequence_x);
	    (*sequence_x) = NULL;
	    
	    (*length_of_sequence_x)      = 0;
	    (*start_of_sequence_x)       = 0;
	    (*end_of_sequence_x)         = 0;
	    
	    if ((*name_of_sequence_y)) delete [] (*name_of_sequence_y);
	    (*name_of_sequence_y) = NULL;
	    if ((*sequence_y)) delete [] (*sequence_y);
	    (*sequence_y) = NULL;
	    
	    (*length_of_sequence_y)      = 0;
	    (*start_of_sequence_y)       = 0;
	    (*end_of_sequence_y)         = 0;
	    
	    check++;
	}
	
	// release memory

	if(Items)
	{
	    for(int i=0; i<Max_number_of_items; i++)
	    {
		if(Items[i]) delete [] Items[i];
		Items[i] = NULL;
	    }
	    delete [] Items;
	}
	Items = NULL;
	
	if(NumberItems)
	{
	    for(i=0; i<Max_number_of_items; i++)
	    {
		if(NumberItems[i]) delete [] NumberItems[i];
		NumberItems[i] = NULL;
	    }
	    delete [] NumberItems;
	}
	NumberItems = NULL;
	
	if (word) delete [] word;
	word = NULL;
	if (word_format) delete [] word_format;
	word_format = NULL;
	if (line) delete [] line;
	line = NULL;
    }
    
    return(check);
}

int get_sequence_pair(FILE* sequence_file_with_pairs_of_sequences,
		      Sequence* const x,
		      Sequence* const y,
		      model_parameters* MP)
{
    int check = 0;

    if (MP->get_Number_of_States() < 2)
    {
	cout << "ERROR: get_sequence_pair: number_of_states_in_hmm ("
	     << MP->get_Number_of_States() << ") < 2.\n" << flush;
	check++;
    }
    
    if (sequence_file_with_pairs_of_sequences == NULL)
    {
	cout << "ERROR: get_sequence_pair: sequence_file_with_pairs_of_sequences is NULL.\n" << flush;
	check++;
    }
    
    if (x == NULL)
    {
	cout << "ERROR: get_sequence_pair: Sequence x is NULL.\n" << flush;
	check++;
    }
    if (y == NULL)
    {
	cout << "ERROR: get_sequence_pair: Sequence y is NULL.\n" << flush;
	check++;
    }

    if (!check)
    {
	char* name_of_sequence_x      = NULL;
	char* sequence_x              = NULL;
	int length_of_sequence_x      = 0;
	int start_of_sequence_x       = 0;
	int end_of_sequence_x         = 0;
      
	char* name_of_sequence_y      = NULL;
	char* sequence_y              = NULL;
	int length_of_sequence_y      = 0;
	int start_of_sequence_y       = 0;
	int end_of_sequence_y         = 0;
	
	check += get_next_pair_of_sequences(sequence_file_with_pairs_of_sequences,
					    // info on sequence x
					    &name_of_sequence_x,
					    &sequence_x,
					    &length_of_sequence_x,
					    &start_of_sequence_x,
					    &end_of_sequence_x,
					    // info on sequence y
					    &name_of_sequence_y,
					    &sequence_y,
					    &length_of_sequence_y,
					    &start_of_sequence_y,
					    &end_of_sequence_y);
	
	
	if (check)
	{
	    cout << "ERROR: get_sequence_pair: error occurred in get_next_pair_of_sequences.\n" << flush;
	}
	
	if (!check)
	{	    
	    if(MP->is_SpecialEmit())
	    {
		Sequence new_x(1,
			       sequence_x,
			       MP);
		Sequence new_y(2,
			       sequence_y,
			       MP);
		
		(*x) = new_x;
		(*y) = new_y;

	    }else{
		Sequence new_x(sequence_x,MP);
		Sequence new_y(sequence_y,MP);
		
		(*x) = new_x;
		(*y) = new_y;

	    }

	    x->set_ac(name_of_sequence_x);
	    y->set_ac(name_of_sequence_y);

	    x->set_sequence_start_and_end(start_of_sequence_x, end_of_sequence_x);
	    y->set_sequence_start_and_end(start_of_sequence_y, end_of_sequence_y);

	}	
	// release memory
	
	if (name_of_sequence_x) delete [] name_of_sequence_x;
	name_of_sequence_x = NULL;
	if (sequence_x) delete [] sequence_x;
	sequence_x = NULL;
	
	if (name_of_sequence_y) delete [] name_of_sequence_y;
	name_of_sequence_y = NULL;
	if (sequence_y) delete [] sequence_y;
	sequence_y = NULL;
    }
    if (check) // if checks failed, return empty sequences
    {
	Sequence new_x;
	Sequence new_y;
	
	(*x) = new_x;
	(*y) = new_y;
    }
    return(check);
}

int get_sequence_decoding_parameters(char* const input_filename,
				     char** name_input_annotation_file,	  
				     char** name_input_sequence_file,
				     char** name_input_tube_file,				 
				     int& algorithm,
				     unsigned long long int& volume,
				     int& radius,
				     char** name_output_file)				      
{
    int check=0;
    TiXmlDocument doc(input_filename);
    bool loadOkay = doc.LoadFile();
    if(!loadOkay){
	cout << "ERROR: get_sequence_decoding_parameters : cannot open file " << input_filename << ".\n" << flush;
	check++;
    }
    TiXmlNode * Node = 0;
    TiXmlElement * Element = 0;
    bool decoding = true;

    Node = doc.FirstChild("HMMConverter");
    if(!Node)
    {
        check++;
	return check;
    }

    Node = Node->FirstChild("sequence_analysis");
    // no sequence analysis has to be done
    if(!Node)  
    {
	decoding = false;
    }else{
	Node = Node->FirstChild("sequence_decoding");
	if(!Node)
	{
	    decoding = false;
	}
    }

    if(!decoding)
    {
	if(*name_input_annotation_file) delete[] (*name_input_annotation_file);
	(*name_input_annotation_file) = NULL;
	if(*name_input_sequence_file) delete[] (*name_input_sequence_file);
	(*name_input_sequence_file) = NULL;
	algorithm = -1;
	volume = 0;
	radius = 0;
	
	if(*name_output_file) delete [] (*name_output_file);
	(*name_output_file) = NULL;

	if(*name_input_tube_file) delete[] (*name_input_tube_file);
	(*name_input_tube_file) = NULL;

	return check;
    }

    if(!check)
    {    
	// read parameters for sequence_decoding
	
	// read algorithm

	Element = Node->FirstChildElement("algorithm");
	
	if(!Element)
	{
	    algorithm = -1;
	}else{
	
	    // read algorithm
	    if(!Element->Attribute("alg"))
	    {
		cout<<"Error : get_sequence_decoding_parameters : "
		    <<"alg attribute in algorithm tag is NULL "<<endl;
		check++;
	    }else{	    
		algorithm = atoi(Element->Attribute("alg"));
		if((algorithm<0)||(algorithm>3)){
		    cout<<"Error : get_sequence_decoding_parameters : "
			<<"alg attribute in the algorithm tag "
			<<algorithm<<" not in the valid range (0..3) "<<endl;
		    check++;
		}else if((algorithm==0)
			 &&(strcmp(Element->Attribute("alg"),"0")))
		{
		    cout<<"Error : get_sequence_decoding_parameters : "
			<<"alg attribute in the algorithm tag "
			<<Element->Attribute("alg")<<" not in the valid range (0..3) "<<endl;
		    check++;
		}
	    }
	    	    
	    // read volume	    
	    if(!Element->Attribute("MaxVolume"))
	    {
		cout<<"Error : get_sequence_decoding_parameters "
		    <<"MaxVolume attribute in algorithm tag is NULL "<<endl;
		check++;
	    }else{
		volume = atol(Element->Attribute("MaxVolume"));
		
		if(((volume==0) &&(strcmp(Element->Attribute("MaxVolume"),"0")))
		   ||(volume<0))
		{
		    cout<<"Error : get_sequence_decoding_parameters : "
			<<"MaxVolume in the algorithm tag "
			<<Element->Attribute("MaxVolume")<<" not in the valid range "<<endl;
		    check++;
		}
	    }

	    
	    if(((algorithm==1)||(algorithm==2))&&(!check)) //blastn or tblastx
	    {
		// read radius	
	
		if((!Element->Attribute("radius")))
		{
		    cout<<"ERROR : get_sequence_decoding_parameters: "
			<<" radius attribute algorithm tag is NULL while"
			<<" algorithm = "<<algorithm<<"."<<endl;
		    check++;
		}else{
		    radius = atoi(Element->Attribute("radius"));
		    if((radius==0)&&(strcmp(Element->Attribute("radius"),"0")))
		    {	       	    
			cout<<"ERROR : get_sequence_decoding_parameters: "
			    <<" radius attribute in algorithm tag"
			    <<" is invalid("<<Element->Attribute("radius")<<endl;
			check++;
		    }
		}		
	    }
	   	    
	}
    }

    // read input_files
    if(!check)
    {
	Element = Node->FirstChildElement("input_files");
	if(!Element)

	{
	    cout<<"Error : get_sequence_decoding_parameters : "
		<<"input_files tag is NULL "<<endl;
	    check++;		
	    return check;
	}	

	// get input file
	    
	if(!Element->Attribute("SeqFile"))
	{
	    cout<<"ERROR: get_sequence_decoding_parameters: "
		<<" SeqFile attribute is NULL in input_files tag!"<<endl;   
	}else
	{
	    (*name_input_sequence_file) = Nullstrcpy(Element->Attribute("SeqFile"),check);	       
	    if(check)
	    {
		cout<<"ERROR: get_sequence_decoding_parameters "
		    <<" error in Nullstrcpy for SeqFile attribute "
		    <<" in input_files tag!"<<endl;
	    }  		
	    
	}

	// get annotation input file
		
	if(Element->Attribute("AnnFile"))
	{
	    (*name_input_annotation_file) = Nullstrcpy(Element->Attribute("AnnFile"),check);
	    if(check)
	    {
		cout<<"ERROR: get_sequence_decoding_parameters: "
		    <<" error in Nullstrcpy for AnnFile attribute "
		    <<" in input_files tag!"<<endl;
	    }  		
	}

	// read tube file
	if((algorithm==3)&&(!check))
	{
	    if(!Element->Attribute("TubeFile"))
	    {
		cout<<"ERROR : get_sequence_decoding_parameters : "
		    <<" TubeFile attribute in input_files tag is NULL "<<endl;
		check++;
	    }else{
		(*name_input_tube_file) = Nullstrcpy(Element->Attribute("TubeFile"),check);
		if(check)
		{
		    cout<<"ERROR: get_sequence_decoding_parameters "
			<<" error in Nullstrcpy for TubeFile attribute "
			<<" in input_files tag!"<<endl;
		}
	    }	
	}
    }
 
    // read viterbi output
    if(!check)
    {

	Element = Node->FirstChildElement("output_files");
	if(!Element){
	    cout<<"No Output file is required for sequence decoding. "<<endl;	    
	}else{	   

	    if(!Element->Attribute("OutFile")){
		cout<<"Error : get_sequence_decoding_parameters : "
		    <<"OutFile attribute in output_files tag is NULL."<<endl;
		check++;
	    }else{
		(*name_output_file)=Nullstrcpy(Element->Attribute("OutFile"),
					       check);		
		if(check)
		{
		    cout<<"Error : get_sequence_decoding_parameters : "
			<<"error in Nullstrcpy for OutputFile attribute in"
			<<"output_files tag"<<endl;
		}		    
	    }					       
	}
    }
 
    Node = 0;
    Element = 0;
    // release memory if error
    if(check)
    {
	algorithm = -1;

	if(*name_input_sequence_file) delete[] (*name_input_sequence_file);
	(*name_input_sequence_file) = NULL;
       
	if(*name_input_tube_file) delete[] (*name_input_tube_file);
	(*name_input_tube_file) = NULL;

	if(*name_input_annotation_file) delete[] (*name_input_annotation_file);
	(*name_input_annotation_file) = NULL;

	if(*name_output_file) delete [] (*name_output_file);
	(*name_output_file) = NULL;

    }
    return check;
}

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
                                      char** name_XMLfile)
{
    int check=0;
    TiXmlDocument doc(input_filename);
    bool loadOkay = doc.LoadFile();
    if(!loadOkay){
	cout << "ERROR: get_parameter_training_parameters : cannot open file " << input_filename << ".\n" << flush;
	check++;
    }
    TiXmlNode * Node = 0;
    TiXmlElement * Element = 0;
    bool training = true;

    Node = doc.FirstChild("HMMConverter");
    if(!Node)
    {
        check++;
        return check;
    }

    Node = Node->FirstChild("sequence_analysis");
    // no sequence analysis has to be done
    if(!Node)  
    {
	training = false;
    }else{
	Node = Node->FirstChild("parameter_training");
	if(!Node)
	{
	    training = false;
	}
    }

    if(!training)
    {
	if(*name_input_annotation_file) delete[] (*name_input_annotation_file);
	(*name_input_annotation_file) = NULL;
	if(*name_input_sequence_file) delete[] (*name_input_sequence_file);
	(*name_input_sequence_file) = NULL;
	algorithm = -1;
	volume = 0;
	radius = 0;
	
	if(*name_transition_prob_file) delete [] (*name_transition_prob_file);
	(*name_transition_prob_file) = NULL;

	if(*name_emission_prob_file) delete [] (*name_emission_prob_file);
	(*name_emission_prob_file) = NULL;
	
	if(*name_XMLfile) delete [] (*name_XMLfile);
	(*name_XMLfile) = NULL;

	if(*name_input_tube_file) delete[] (*name_input_tube_file);
	(*name_input_tube_file) = NULL;

	return check;
    }

    if(!check)
    {    
	// read parameters for sequence_decoding
	
	// read algorithm

	Element = Node->FirstChildElement("algorithm");
	
	if(!Element)
	{
	    algorithm = -1;
	}else{
	
	    // read algorithm
	    if(!Element->Attribute("alg"))
	    {
		cout<<"Error : get_parameter_training_parameters : "
		    <<"alg attribute in algorithm tag is NULL "<<endl;
		check++;
	    }else{	    
		algorithm = atoi(Element->Attribute("alg"));
		if((algorithm<0)||(algorithm>2)){
		    cout<<"Error : get_parameter_training_parameters : "
			<<"alg attribute in the algorithm tag "
			<<algorithm<<" not in the valid range (0..2) "<<endl;
		    check++;
		}else if((algorithm==0)
			 &&(strcmp(Element->Attribute("alg"),"0")))
		{
		    cout<<"Error : get_parameter_training_parameters : "
			<<"alg attribute in the algorithm tag "
			<<Element->Attribute("alg")<<" not in the valid range (0..1) "<<endl;
		    check++;
		}
	    }
	    	    
	    // read volume	    
	    if(!Element->Attribute("MaxVolume"))
	    {
		cout<<"Error : get_parameter_training_parameters "
		    <<"MaxVolume attribute in algorithm tag is NULL "<<endl;
		check++;
	    }else{
		volume = atol(Element->Attribute("MaxVolume"));
		
		if(((volume==0) &&(strcmp(Element->Attribute("MaxVolume"),"0")))
		   ||(volume<0))
		{
		    cout<<"Error : get_parameter_parameters : "
			<<"MaxVolume in the algorithm tag "
			<<Element->Attribute("MaxVolume")<<" not in the valid range "<<endl;
		    check++;
		}
	    }

	    // read radius
	    if(!Element->Attribute("radius"))
	    {
	        radius = 0;
	    }else{
		radius = atoi(Element->Attribute("radius"));
		
		if(((radius==0) &&(strcmp(Element->Attribute("radius"),"0")))
		   ||(radius<0))
		{
		    cout<<"Error : get_parameter_parameters : "
			<<"radius in the algorithm tag "
			<<Element->Attribute("radius")<<" not in the valid range "<<endl;
		    check++;
		}
	    }

	    // read threshold	
	    if(!Element->Attribute("threshold"))
	    {
		// set default
		threshold = 0;
	    }else {
		threshold = atof(Element->Attribute("threshold"));
		if((threshold==0)&&(strcmp(Element->Attribute("threshold"),"0")))
		{	       	    
		    cout<<"ERROR : get_parameter_training_parameters : "
			<<" threshold attribute in algorithm tag"
			<<" is invalid("<<Element->Attribute("threshold")<<")"<<endl;
		    check++;
		}
	    }

	    // read Maxiter
	    if(!Element->Attribute("Maxiter"))
	    {
		// set default
		Maxiter = 10;
	    }else {
		Maxiter = atoi(Element->Attribute("Maxiter"));
		if((Maxiter==0)&&(strcmp(Element->Attribute("Maxiter"),"0")))
		{	       	    
		    cout<<"ERROR : get_parameter_training_parameters : "
			<<" Maxiter attribute in algrotihm tag"
			<<" is invalid("<<Element->Attribute("Maxiter")<<")"<<endl;
		    check++;
		}
	    }
	   	
	    if(!Element->Attribute("SamplePaths"))
	    {
		// set default
		Samples = 0;
	    }else {
		Samples = atoi(Element->Attribute("SamplePaths"));
		if((Samples==0)&&(strcmp(Element->Attribute("SamplePaths"),"0")))
		{	       	    
		    cout<<"ERROR : get_parameter_training_parameters : "
			<<" SamplePaths attribute in algrotihm tag"
			<<" is invalid("<<Element->Attribute("SamplePaths")<<")"<<endl;
		    check++;
		}
	    }
    
	}
    }

    // read input_files
    if(!check)
    {
	Element = Node->FirstChildElement("input_files");
	if(!Element)

	{
	    cout<<"Error : get_parameter_training_parameters : "
		<<"input_files tag is NULL "<<endl;
	    check++;		
	    return check;
	}	

	// get input file
	    
	if(!Element->Attribute("SeqFile"))
	{
	    cout<<"ERROR: get_parameter_training_parameters: "
		<<" SeqFile attribute is NULL in input_files tag!"<<endl;   
	}else
	{
	    (*name_input_sequence_file) = Nullstrcpy(Element->Attribute("SeqFile"),check);	       
	    if(check)
	    {
		cout<<"ERROR: get_parameter_training_parameters "
		    <<" error in Nullstrcpy for SeqFile attribute "
		    <<" in input_files tag!"<<endl;
	    }  		
	    
	}

	// get annotation input file
		
	if(Element->Attribute("AnnFile"))
	{
	    (*name_input_annotation_file) = Nullstrcpy(Element->Attribute("AnnFile"),check);
	    if(check)
	    {
		cout<<"ERROR: get_parameter_training_parameters: "
		    <<" error in Nullstrcpy for AnnFile attribute "
		    <<" in input_files tag!"<<endl;
	    }  		
	}

	// read tube file
	if(Element->Attribute("TubeFile"))
	{	  
	    (*name_input_tube_file) = Nullstrcpy(Element->Attribute("TubeFile"),check);
	    if(check)
	    {
		cout<<"ERROR: get_parameter_training_parameters "
		    <<" error in Nullstrcpy for TubeFile attribute "
		    <<" in input_files tag!"<<endl;
		
	    }	    
	}
    }
 
    // read output parameters
    if(!check)
    {
	Element = Node->FirstChildElement("output_files");
	if(!Element){
	    cout<<"No Output file is required for parameter training. "<<endl;	    
	}else{	   

	    if(Element->Attribute("XMLFile")){
		(*name_XMLfile)=Nullstrcpy(Element->Attribute("XMLFile"),
					   check);		
		if(check)
		{
		    cout<<"Error : get_parameter_training_parameters : "
			<<"error in Nullstrcpy for XMLFile attribute in"
			<<"output_files tag"<<endl;
		}		    
	    }	

	    if(Element->Attribute("TProbFile")){
		(*name_transition_prob_file)=Nullstrcpy(Element->Attribute("TProbFile"),
							check);		
		if(check)
		{
		    cout<<"Error : get_parameter_training_parameters : "
			<<"error in Nullstrcpy for TProbFile attribute in"
			<<"output_files tag"<<endl;
		}		    
	    }	
			
	    if(Element->Attribute("EProbFile")){
		(*name_emission_prob_file)=Nullstrcpy(Element->Attribute("EProbFile"),
						      check);		
		if(check)
		{
		    cout<<"Error : get_parameter_training_parameters : "
			<<"error in Nullstrcpy for EProbFile attribute in"
			<<"output_files tag"<<endl;
		}		    
	    }	
	}
    }
 
    Node = 0;
    Element = 0;
    // release memory if error
    if(check)
    {
	algorithm = -1;

	if(*name_input_sequence_file) delete[] (*name_input_sequence_file);
	(*name_input_sequence_file) = NULL;

	if(*name_input_tube_file) delete[] (*name_input_tube_file);
	(*name_input_tube_file) = NULL;
       
	if(*name_input_annotation_file) delete[] (*name_input_annotation_file);
	(*name_input_annotation_file) = NULL;

	if(*name_transition_prob_file) delete [] (*name_transition_prob_file);
	(*name_transition_prob_file) = NULL;

	if(*name_emission_prob_file) delete [] (*name_emission_prob_file);
	(*name_emission_prob_file) = NULL;
	
	if(*name_XMLfile) delete [] (*name_XMLfile);
	(*name_XMLfile) = NULL;

    }
    return check;
}

int get_tube_from_file(const char* filename,
		       const char* seq_x_name,
		       const char* seq_y_name,
		       const int seq_x_length,
		       const int seq_y_length,
		       Tube<int>* tTube)
{
    int check = 0;
    // check input
    if(!filename)
    {
	cout<<"ERROR: get_tube_from_file: input filename is NULL."<<endl;
	check++;
    }
    if(!seq_x_name)
    {
	cout<<"ERROR: get_tube_from_file: input seq_x_name is NULL."<<endl;
	check++;
    }
    if(!seq_y_name)
    {
	cout<<"ERROR: get_tube_from_file: input seq_y_name is NULL."<<endl;
	check++;
    }
    if(seq_x_length<0)
    {
	cout<<"ERROR: get_tube_from_file: input seq_x_length("
	    <<seq_x_length<<") out of range."<<endl;		
	check++;
    }
    if(seq_y_length<0)
    {
	cout<<"ERROR: get_tube_from_file: input seq_y_length("
	    <<seq_y_length<<") out of range."<<endl;		
	check++;
    }
    
    if(check)
    {
	return check;
    }

    int i = 0;
    char** Items = new char*[Max_number_of_items];
    for(i=0; i<Max_number_of_items; i++)
    {
	Items[i] = new char[Max_word_length];
	strcpy(Items[i]," ");
    }
    int NoOfItems = 0;

    char** ItemsForLabel = new char*[Max_number_of_items];
    for(i=0; i<Max_number_of_items; i++)
    {
	ItemsForLabel[i] = new char[Max_word_length];
	strcpy(ItemsForLabel[i]," ");
    }
    int NoOfItemsForLabel = 0;

    char* line = new char[Max_line_length];
    
    char* seqpair = new char[Max_word_length];
    // set seqpair
    strcpy(seqpair,">");
    strcat(seqpair,seq_x_name);
    strcat(seqpair,"&");
    strcat(seqpair,seq_y_name);
    bool seqmatch = false;
    bool block = false;

    FILE* tube_file = fopen(filename,"rt");
    if(!tube_file)
    {
	cout<<"ERROR: get_tube_from_file: can not open file :"
	    <<filename<<", no tube is got. "<<endl;
	check++;
    }
    int countline = 0;

    int n_of_matches = 0;
    Match* matches = NULL;
    matches = new Match[Max_number_of_gtf_lines];

    int n_of_matches_new = 0;
    Match* matches_new = NULL;      

    int radius = 0;

    int xstart = -1;
    int ystart = -1;
    int xend = -1;
    int yend = -1;

    int length_covered_x = 0;
    
    while((!feof(tube_file)) && (!check))
    {
	fgets(line, Max_line_length-1, tube_file);
	
	check += splitstring(static_cast<const char*>(line),NoOfItems,&Items,' ');

	if(check)
	{
	    cout<<"ERROR: get_tube_from_file: error in splitstring on line : "
		<<countline<<endl;
	    break;
	}

	if((seqmatch)&&(Items[0][0]!='>')) // not header line
	{
	    if(NoOfItems!=4)
	    {
		cout<<"ERROR: get_tube_from_file: interruption on line("
		    <<countline<<") : "<<endl<<line<<endl;		
		check++;
		break;
	    }else{
		// check range

		xstart = atoi(Items[0]);
		if((xstart==0)&&(strcmp(Items[0],"0")))
		{
		    cout<<"ERROR: get_tube_from_file: interruption on line("
			<<countline<<") : "<<"xstart("<<Items[0]
			<<") out of range."<<endl;
		    check++;
		    break;
		}

		ystart = atoi(Items[1]);
		if((ystart==0)&&(strcmp(Items[1],"0")))
		{
		    cout<<"ERROR: get_tube_from_file: interruption on line("
			<<countline<<") : "<<"ystart("<<Items[1]
			<<") out of range."<<endl;
		    check++;
		    break;
		}
		
		xend = atoi(Items[2]);
		if((xend==0)&&(strcmp(Items[2],"0")))
		{
		    cout<<"ERROR: get_tube_from_file: interruption on line("
			<<countline<<") : "<<"xend("<<Items[2]
			<<") out of range."<<endl;
		    check++;
		    break;
		}
		
		yend = atoi(Items[3]);
		if((yend==0)&&(strcmp(Items[3],"0 ")))
		{
		    cout<<"ERROR: get_tube_from_file: interruption on line("
			<<countline<<") : "<<"yend("<<Items[1]
			<<") out of range."<<endl;
		    check++;
		    break;
		}

		if((xstart<0)||(ystart<0))
		{
		    cout<<"ERROR: get_tube_from_file: interruption on line("
			<<countline<<") : "<<"xstart("<<xstart<<") or "
			<<"ystart ("<<ystart<<") out of range."<<endl;			
		    check++;
		    break;
		}else if((xend>seq_x_length)||(yend>seq_y_length))
		{
		    cout<<"ERROR: get_tube_from_file: interruption on line("
			<<countline<<") : "<<"xend("<<xend<<") or "
			<<"yend ("<<yend<<") out of range."<<endl;			     
		    check++;
		    break;
		}		
		//xstart, ystart, xend, yend       	
		if(block)
		{
		    length_covered_x += xend-xstart+1;		    
		}
		Match next_match(xstart,ystart,xend,yend,0);
		matches[n_of_matches] = next_match;
		n_of_matches++;		
	    }	    
	}	

	if(NoOfItems<=0)
	{	    
	    check++;
	    break;
	}else{

	    if(Items[0][0]=='>') // header line
	    {	
		if(seqmatch)
		{
		    break;
		}

		if(NoOfItems==1)
		{
		    block = true;
		}else if(NoOfItems==2)
		{
		    block = false;
		}else
		{
		    cout<<"ERROR: get_tube_from_file: "
			<<"error in headline format interruption on line("
			<<countline<<") : "<<endl<<line<<endl;			 
		    check++;
		    break;
		}
		if(!strcmp(Items[0],seqpair))
		{	
		    if(!block)
		    {
			// get radius
			check+=splitstring(Items[1],NoOfItemsForLabel,&ItemsForLabel,':');
			
			if(check)
			{
			    cout<<"ERROR: get_tube_from_file: "
				<<"error in splitstring for radius, interruption on line("
				<<countline<<") : "<<endl<<line<<endl;		
			}
			
			radius = atoi(ItemsForLabel[1]);
			
			if((radius<0)||((radius==0)&&strcmp("0",ItemsForLabel[1])))
			{
			    cout<<"ERROR: get_tube_from_file: "
				<<"error in getting radius interruption on line("
				<<countline<<") : "<<endl<<line<<endl;		
			    check++;
			    break;
			}
		    }
	
		    seqmatch = true;
		}	       
	    }
	}		
    }
       
    fclose(tube_file);
    
    if(!check)
    {
	Tube<int> tmp_tube;

	if(seqmatch)
	{		  
	    if(!block)
	    {
		check+= find_maximal_subset_of_matches(n_of_matches, matches, &n_of_matches_new, &matches_new);
		if(check)
		{
		    cout<<"ERROR: get_tube_from_file: error in find_maximal_subset_of_matches.\n"<<flush;   
		}	    
		
		tmp_tube = convert_matches_to_tube(seq_x_length,seq_y_length,n_of_matches_new,
						   matches_new, radius); //1 is MaxGap, will remove
	    }else{
		
		vector<vector<int> > upper(length_covered_x,2);
		vector<vector<int> > lower(length_covered_x,2);
		
		int count_x = 0;
		
		for(i=0; i<n_of_matches; i++)
		{
		    if(check)
		    {
			break;
		    }
		    int count = matches[i].Start_x();		
		    while(count <= matches[i].End_x())
		    {
			if(count_x>=length_covered_x)
			{
			    cout<<"Error: get_tube_from_file: count_x("<<count_x
				<<") >= length_covered_x("<<length_covered_x<<")."<<endl;
			    check++;
			    break;
			}
			lower[count_x][0] = count;
			lower[count_x][1] = matches[i].Start_y();
			upper[count_x][0] = count;
			upper[count_x][1] = matches[i].End_y();
			count++;
			count_x++;		    
		    }
		}		

		tmp_tube = convert_boundaries_to_tube(seq_x_length,seq_y_length,upper,lower);
	    }
	}
	(*tTube) = tmp_tube;
    }

    if(check)
    {
	Tube<int> tmp_tube;
	(*tTube) = tmp_tube;
    }

    if(Items)
    {
	for(i=0; i<Max_number_of_items; i++)
	{
	    if(Items[i]) delete [] Items[i];
	    Items[i] = NULL;
	}
	delete [] Items;
    }
    Items = NULL;
    NoOfItems = 0;

    if(ItemsForLabel)
    {
	for(i=0; i<Max_number_of_items; i++)
	{
	    if(ItemsForLabel[i]) delete [] ItemsForLabel[i];
	    ItemsForLabel[i] = NULL;
	}
	delete [] ItemsForLabel;
    }
    ItemsForLabel = NULL;
    NoOfItemsForLabel = 0;

    if(line) delete [] line;
    line = NULL;

    if(seqpair) delete [] seqpair;
    seqpair = NULL;

    if(matches) delete [] matches;
    matches = NULL;
    
    if(matches_new) delete [] matches_new;
    matches_new = NULL;
    return check;
}



