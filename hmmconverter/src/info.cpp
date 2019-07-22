/* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare info-class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/info.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

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
#include <string.h>
#include <string>
#include "info.h"

/* define constructors */

Info::Info()
{
    info_number_of_sources = 0;
    info_number_of_lines = NULL;
    
    info_number_of_annotation_labels = 0;
    info_number_of_each_annotation_label = NULL;

    info_seq_names       = NULL;              
    info_source_names    = NULL;	    

    start_positions = NULL;        
    end_positions   = NULL;          
    info_scores          = NULL;

    info_annotation_labels = NULL;

    output          = NULL;
  
}

Info::Info(const int n_of_sources,
	   int* const n_of_lines,
	   model_parameters* const MP,
	   int& check
    )
{
    check = 0;
    
    if (n_of_sources <= 0)
    {
	cout << "ERROR: class Info::constructor : n_of_sources (" << n_of_sources << ") is <= 0.\n";
	check++;
    }
    if(!n_of_lines)
    {
	cout << "ERROR: class Info:: constructor : n_of_lines is NULL.\n ";
	check++;
    }
   
    if (check == 0)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	info_number_of_sources = n_of_sources;
	info_number_of_lines = new int[info_number_of_sources];
	for(i=0; i<info_number_of_sources; i++){
	    info_number_of_lines[i] = n_of_lines[i];
	}
	
	info_number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
	
	info_number_of_each_annotation_label = new int [info_number_of_annotation_labels];
	
	for(i=0;i<info_number_of_annotation_labels; i++){
	    info_number_of_each_annotation_label[i] = MP->get_Annotation_Label_size(i);
	}
	
	info_source_names    = new char*     [info_number_of_sources];	
	info_seq_names       = new char**    [info_number_of_sources];
        
	info_annotation_labels      = new bool***[info_number_of_sources];
	info_scores                 = new double*** [info_number_of_sources];
	start_positions = new int*      [info_number_of_sources];  
	end_positions   = new int*      [info_number_of_sources];  	
	
	for (i=0; i<info_number_of_sources; i++)
	{
	  
	    info_source_names[i]      = new char     [Max_word_length];	    
	    info_seq_names[i]         = new char*    [info_number_of_lines[i]];
	  
	    info_annotation_labels[i] = new bool**   [info_number_of_lines[i]];
	    info_scores[i]            = new double** [info_number_of_lines[i]];

	    start_positions[i]   = new int      [info_number_of_lines[i]];  
	    end_positions[i]     = new int      [info_number_of_lines[i]];

	    strcpy(info_source_names[i], " ");
	    
	    for(j=0; j<info_number_of_lines[i]; j++){
	
		info_seq_names[i][j]         = new char    [Max_word_length];
		
		info_annotation_labels[i][j] = new bool*   [info_number_of_annotation_labels];
		info_scores[i][j]            = new double* [info_number_of_annotation_labels];

		strcpy(info_seq_names[i][j], " ");			
		start_positions[i][j] = 0;  
		end_positions[i][j]   = 0;
   
		for(k=0; k<info_number_of_annotation_labels; k++){
		    info_annotation_labels[i][j][k] = new bool   [info_number_of_each_annotation_label[k]];
		    info_scores[i][j][k]            = new double [info_number_of_each_annotation_label[k]];
		    for(l=0; l<info_number_of_each_annotation_label[k]; l++){
			info_annotation_labels[i][j][k][l] = false;
			info_scores[i][j][k][l] = 0;
		    }
		}
	    }
	    
	}
	
	output          = NULL;
    }
}

Info::Info(const Info &s)
{
    *this=s;
}

Info::~Info()
{
    int i = 0;
    int j = 0;
    int k = 0;
    
    if(info_source_names){
	for(i=0; i<info_number_of_sources; i++){
	    if(info_source_names[i]) delete[] info_source_names[i];
	    info_source_names[i] = NULL;
	}
	delete [] info_source_names;
    }
    
    info_source_names = NULL;
    
    
    if(info_seq_names){
	for(i=0; i<info_number_of_sources; i++){
	    if(info_seq_names[i]){
		for(j=0; j<info_number_of_lines[i]; j++){
		    if(info_seq_names[i][j]) delete[] info_seq_names[i][j];
		    info_seq_names[i][j] = NULL;
		}
		delete [] info_seq_names[i];
	    }
	    info_seq_names[i] = NULL;
	}
	delete[] info_seq_names;
    }
    info_seq_names = NULL;

    
    if(info_annotation_labels){
	for(i=0; i<info_number_of_sources; i++){
	    if(info_annotation_labels[i]){
		for(j=0; j<info_number_of_lines[i]; j++){
		    if(info_annotation_labels[i][j]){
			for(k=0; k<info_number_of_annotation_labels; k++){
			    if(info_annotation_labels[i][j][k]) delete[] info_annotation_labels[i][j][k];
			    info_annotation_labels[i][j][k] = NULL;
			}
			delete[] info_annotation_labels[i][j];
		    }
		    info_annotation_labels[i][j] = NULL;
		}
		delete[] info_annotation_labels[i];
	    }
	    info_annotation_labels[i] = NULL;
	}
	delete[] info_annotation_labels;
    }
    info_annotation_labels = NULL;
    
    if(info_scores){
	for(i=0; i<info_number_of_sources; i++){
	    if(info_scores[i]){
		for(j=0; j<info_number_of_lines[i]; j++){
		    if(info_scores[i][j]){
			for(k=0; k<info_number_of_annotation_labels; k++){
			    if(info_scores[i][j][k]) delete[] info_scores[i][j][k];
			    info_scores[i][j][k] = NULL;
			}
			delete[] info_scores[i][j];
		    }
		    info_scores[i][j] = NULL;
		}
		delete[] info_scores[i];
	    }
	    info_scores[i] = NULL;
	}
	delete[] info_scores;
    }
    info_scores = NULL;

    if(start_positions){
	for(i=0; i<info_number_of_sources; i++)
	{
	    if(start_positions[i]) delete[] start_positions[i];
	    start_positions[i] = NULL;
	}
	delete[] start_positions;
    }
    start_positions = NULL;
    
    if(end_positions){
	for(i=0; i<info_number_of_sources; i++){
	    if(end_positions[i]) delete[] end_positions[i];
	    end_positions[i] = NULL;
	}
	delete[] end_positions;
    }
    end_positions = NULL;

    if(info_number_of_lines)  delete [] info_number_of_lines;
    info_number_of_lines = NULL;
    
    if(info_number_of_each_annotation_label) delete [] info_number_of_each_annotation_label;
    info_number_of_each_annotation_label = NULL;

    if (output)          delete [] output;
    output = NULL;

    info_number_of_sources = 0;
    info_number_of_annotation_labels = 0;
    
}

Info & Info::operator = (const Info &s)
{
    if ( this != &s)
    {
	this->~Info();

	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	info_number_of_sources = s.info_number_of_sources;
	info_number_of_lines = new int[info_number_of_sources];
	for(i=0; i<info_number_of_sources; i++){
	    info_number_of_lines[i] = s.info_number_of_lines[i];
	}
	
	info_number_of_annotation_labels = s.info_number_of_annotation_labels;
	info_number_of_each_annotation_label = new int[info_number_of_annotation_labels];
	for(i=0; i<info_number_of_annotation_labels; i++){
	    info_number_of_each_annotation_label[i] = s.info_number_of_each_annotation_label[i];
	}
	
	if (output) delete [] output;
	output = NULL;
	
	if (s.output != NULL)
	{
	    if (static_cast<int>(strlen(s.output)) > 0)
	    {
		const int length = static_cast<int>(strlen(s.output))+1;
		output = new char[length];
		strcpy(output, s.output);
	    }
	}	

	if (info_number_of_sources > 0)
	{
	    info_source_names      = new char*      [info_number_of_sources];
	    info_seq_names         = new char**     [info_number_of_sources];        
	    
	    info_annotation_labels = new bool***    [info_number_of_sources];
	    info_scores            = new double***  [info_number_of_sources];
	  
	    start_positions   = new int*       [info_number_of_sources];  
	    end_positions     = new int*       [info_number_of_sources];

	    for (i=0; i<info_number_of_sources; i++)
	    {
		info_source_names[i]       = new char      [Max_word_length];
		info_seq_names[i]          = new char*     [info_number_of_lines[i]];
		
		info_annotation_labels[i]  = new bool**    [info_number_of_lines[i]];
		info_scores[i]             = new double**  [info_number_of_lines[i]];
		
		start_positions[i]    = new int       [info_number_of_lines[i]];
		end_positions[i]      = new int       [info_number_of_lines[i]];

		strcpy(info_source_names[i], s.info_source_names[i]);
	
		for(j=0; j<info_number_of_lines[i]; j++){
		    info_seq_names[i][j]         = new char     [Max_word_length];
		    
		    info_annotation_labels[i][j] = new bool*    [info_number_of_annotation_labels];
		    info_scores[i][j]            = new double*  [info_number_of_annotation_labels];

		    strcpy(info_seq_names[i][j], s.info_seq_names[i][j]);
		    start_positions[i][j] = s.start_positions[i][j];
		    end_positions[i][j]   = s.end_positions[i][j];
		 
		    for(k=0; k<info_number_of_annotation_labels; k++){
			info_annotation_labels[i][j][k] = new bool   [info_number_of_each_annotation_label[k]];
			info_scores[i][j][k]            = new double [info_number_of_each_annotation_label[k]];
			for(l=0; l<info_number_of_each_annotation_label[k]; l++){
			    info_annotation_labels[i][j][k][l] = s.info_annotation_labels[i][j][k][l];
			    info_scores[i][j][k][l]            = s.info_scores[i][j][k][l];
			}
		    }
		}
	    }
	}
    }
    return(*this);
}

// mutators

int Info::set_info_number_of_sources(const int n_of_sources)
{
    int check = 0;
    if(n_of_sources<0){
	cout << "ERROR : class Info::set_number_of_sources : n_of_sources (" 
	     << n_of_sources << ") out of range.\n";
	check++;
    }else{
	info_number_of_sources = n_of_sources;
    }
    return check;
}

int Info::set_info_number_of_lines(const int source, const int n_of_lines)
{
    int check = 0;
    if((source<0)||(source >= info_number_of_sources)){
	cout << "ERROR : class Info::set_info_number_of_lines : source (" << info_number_of_sources << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }
    if(n_of_lines<0){
	cout << "ERROR : class Info::set_info_number_of_lines : n_of_lines (" 
	     << n_of_lines << ") out of range.\n";
	check++;
    }
    info_number_of_lines[source] = n_of_lines;
    return check;
}

int Info::set_info_seq_names(const int source, const int line, const char* seq_name)
{
    int check = 0;
    
    if((source<0)||(source >= info_number_of_sources)){
	cout << "ERROR : class Info::set_seq_name : source (" << info_number_of_sources << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if ((line < 0) || (line > info_number_of_lines[source]-1))
	{
	    cout << "ERROR : class Info::set_seq_name : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if (seq_name != NULL)
    {
	if (static_cast<int>(strlen(seq_name)) > Max_word_length)
	{
	    cout << "ERROR : class Info::set_seq_name : length of seq_name " << seq_name << " (" 
		 <<  strlen(seq_name) << ") too long ( > Max_word_length = " << Max_word_length << ").\n";
	    check++;
	}
    }
    else
    {
	cout << "ERROR : class Info::set_seq_name : seq_name is NULL.\n";
	check++;
    }
    
    if (check == 0)
    {
	strcpy(info_seq_names[source][line], seq_name);
    }
    return(check);
}

int Info::set_info_source_names(const int index, const char* source_name)
{
    int check = 0;
    
    if ((index < 0) || (index > info_number_of_sources-1))
    {
	cout << "ERROR : class Info::set_source_name : index (" << index << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }
    
    if (source_name != NULL)
    {
	if (static_cast<int>(strlen(source_name)) > Max_word_length)
	{
	    cout << "ERROR : class Info::set_source_name : length of source_name " << source_name << " (" 
		 <<  strlen(source_name) << ") too long ( > Max_word_length = " << Max_word_length << ").\n";
	    check++;
	}
    }
    else
    {
	cout << "ERROR : class Info::set_source_name : source_name is NULL.\n";
	check++;
    }
    
    if (check == 0)
    {
	strcpy(info_source_names[index], source_name);
    }
    return(check);
}

int Info::set_info_number_of_annotation_labels(const int n_of_annotation_labels)
{
    int check = 0;
    if(n_of_annotation_labels<0){
	cout << "ERROR : class Info::set_info_number_of_annotation_labels : n_of_annotation_labels (" 
	     << n_of_annotation_labels << ") out of range.\n";
	check++;
    }else{
	info_number_of_annotation_labels = n_of_annotation_labels;
    }
    return check;
}

int Info::set_info_number_of_each_annotation_label(const int type, const int n_of_label)
{
    int check = 0;
    if ((type < 0) || (type > info_number_of_annotation_labels-1))
    {
	cout << "ERROR : class Info::set_info_number_of_each_annotation_label : type (" << type << ") out of range [0,"
	     << info_number_of_annotation_labels-1 << "].\n";
	check++;
    }
    if(!check){
	info_number_of_each_annotation_label[type] = n_of_label;
    }
    return check;
}
    
int Info::set_info_annotation_labels(const int source, const int line, const int type, const int label, const bool on)
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::set_info_annotation_labels : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::set_info_annotation_labels : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if ((type < 0) || (type > info_number_of_annotation_labels-1))
    {
	cout << "ERROR : class Info::set_info_annotation_labels : type (" << type << ") out of range [0,"
	     << info_number_of_annotation_labels-1 << "].\n";
	check++;
    }else{
	if((label<0)||(label>=info_number_of_each_annotation_label[type])){
	    cout << "ERROR : class Info::set_info_annotation_labels : label (" << label << ") out of range [0,"
		 << info_number_of_each_annotation_label[type]-1 << "].\n";
	    check++;
	}
    }
    
    if (check == 0)
    {
	info_annotation_labels[source][line][type][label] = on;
    }
    return(check);
}

int Info::set_info_scores(const int source, const int line, const int type, const int label, const double s)
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::set_info_scores : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::set_info_scores : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if ((type < 0) || (type > info_number_of_annotation_labels-1))
    {
	cout << "ERROR : class Info::set_info_scores : type (" << type << ") out of range [0,"
	     << info_number_of_annotation_labels-1 << "].\n";
	check++;
    }else{
	if((label<0)||(label>=info_number_of_each_annotation_label[type])){
	    cout << "ERROR : class Info::set_info_scores : label (" << label << ") out of range [0,"
		 << info_number_of_each_annotation_label[type]-1 << "].\n";
	    check++;
	}
    }
    if (check == 0)
    {
	info_scores[source][line][type][label] = s;
    }
    return(check);
}

int Info::set_start_positions(const int source, const int line, const int start_pos)
{
    int check = 0;

    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::set_start_positions : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::set_start_positions : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if (check == 0)
    {
	start_positions[source][line] = start_pos;
    }
    return(check);
}

int Info::set_end_positions(const int source, const int line, const int end_pos)
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::set_end_positions : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::set_end_positions : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if (check == 0)
    {
	end_positions[source][line] = end_pos;
    }
    return(check);
}

int Info::set_output(const char* o)
{
    int check = 0;
    if(!o){
	cout<<"ERROR: class Info:: set_output, input o is NULL.\n";
	check++;
    }
    if(!check){
	if(!output){
	    output = Nullstrcpy(o,check);
	    if(check){
		cout<<"ERROR: class Info:: set_output, can not copy o: "<<o<<"to output "<<endl;
	    }
	}else{
	    strcpy(output,o);
	}
    }
    return check;
}

int Info::get_info_number_of_sources(void) const 
{
    return (info_number_of_sources);
}   

int Info::get_info_number_of_lines(const int source) const 
{
    int check = 0;
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::get_numebr_of_lines : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }
    if(!check){
	return(info_number_of_lines[source]);
    }else{
	return -1;
    }
}

char* Info::get_info_seq_names (const int source, const int line) const
{
    int check = 0;
    if((source<0)||(source >= info_number_of_sources)){
	cout << "ERROR : class Info::get_info_seq_names : source (" << info_number_of_sources << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if ((line < 0) || (line > info_number_of_lines[source]-1))
	{
	    cout << "ERROR : class Info::get_info_seq_names : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if(!info_seq_names){
	cout<<"ERROR : class Info::get_info_seq_names : info_seq_names is NULL.\n";
	check++;
    }
    if(!check){
	return info_seq_names[source][line];
    }else{
	return NULL;
    }
}

char* Info::get_info_source_names (const int source) const
{
    int check = 0;
    if((source<0)||(source >= info_number_of_sources)){
	cout << "ERROR : class Info::get_info_source_names : source (" << info_number_of_sources << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }
    if(!info_source_names){
	cout<<"ERROR : class Info::get_info_source_names : info_source_names is NULL.\n";
	check++;
    }
    if(!check){
	return info_source_names[source];
    }else{
	return NULL;
    } 
}

int Info::get_info_number_of_annotation_labels() const
{
    return (info_number_of_annotation_labels);
}

int Info::get_info_number_of_each_annotation_label(const int type) const
{
    int check = 0;
    if ((type < 0) || (type > info_number_of_annotation_labels-1))
    {
	cout << "ERROR : class Info::get_info_number_of_each_annotation_label : type (" << type << ") out of range [0,"
	     << info_number_of_annotation_labels-1 << "].\n";
	check++;
    }
    if(!check){
	return info_number_of_each_annotation_label[type];
    }else{
	return -1;
    }
}

bool Info::get_info_annotation_labels(const int source, const int line, const int type, const int label) const
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::get_info_annotation_labels : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::get_info_annotation_labels : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if ((type < 0) || (type > info_number_of_annotation_labels-1))
    {
	cout << "ERROR : class Info::get_info_annotation_labels : type (" << type << ") out of range [0,"
	     << info_number_of_annotation_labels-1 << "].\n";
	check++;
    }else{
	if((label<0)||(label>=info_number_of_each_annotation_label[type])){
	    cout << "ERROR : class Info::get_info_annotation_labels : label (" << label << ") out of range [0,"
		 << info_number_of_each_annotation_label[type]-1 << "].\n";
	    check++;
	}
    }
    if(!check){
	return info_annotation_labels[source][line][type][label];
    }else{
	return false;
    }
}

double Info::get_info_scores(const int source, const int line, const int type, const int label) const
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::get_info_scores : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::get_info_scores : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if ((type < 0) || (type > info_number_of_annotation_labels-1))
    {
	cout << "ERROR : class Info::get_info_scores : type (" << type << ") out of range [0,"
	     << info_number_of_annotation_labels-1 << "].\n";
	check++;
    }else{
	if((label<0)||(label>=info_number_of_each_annotation_label[type])){
	    cout << "ERROR : class Info::get_info_scores : label (" << label << ") out of range [0,"
		 << info_number_of_each_annotation_label[type]-1 << "].\n";
	    check++;
	}
    }
    if (check == 0)
    {
	return info_scores[source][line][type][label];
    }else{
	return 0;
    }
}

int Info::get_start_positions (const int source, const int line) const
{
    int check = 0;

    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::get_start_positions : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::get_start_positions : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    if(!check){
	return start_positions[source][line];
    }else{
	return -1;
    }
}

int Info::get_end_positions (const int source, const int line) const
{
     int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout << "ERROR : class Info::get_end_positions : source (" << source << ") out of range [0,"
	     << info_number_of_sources-1 << "].\n";
	check++;
    }else{
	if((line<0)||(line>=info_number_of_lines[source])){
	    cout << "ERROR : class Info::get_end_positions : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;
	}
    }
    
    if(!check){
	return end_positions[source][line];
    }else{
	return -1;
    }
}

char* Info:: get_output() const
{
    return output;
}

int Info::remove_source(const int source)
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout<<" ERROR : class Info::remove_source : source ("<<source<<") out of range [0,"
	    << info_number_of_sources-1 <<"].\n";
	check++;
    }

    if(!check)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	
	for (i = source; i<info_number_of_sources-1; i++)
	{
	    strcpy(info_source_names[i],info_source_names[i+1]);
	    
	    //remove memory first
	    if(info_seq_names[i]){
		for(j=0; j<info_number_of_lines[i]; j++){
		    if(info_seq_names[i][j]) delete [] info_seq_names[i][j];
		    info_seq_names[i][j] = NULL;
		}
		delete[] info_seq_names[i];
	    }
	    info_seq_names[i] = NULL;

	    if(info_annotation_labels[i]){
		for(j=0;j<info_number_of_lines[i];j++){
		    if(info_annotation_labels[i][j]){
			for(k=0;k<info_number_of_annotation_labels;k++){
			    if(info_annotation_labels[i][j][k]) delete [] info_annotation_labels[i][j][k];
			    info_annotation_labels[i][j][k] = NULL;
			}
			delete [] info_annotation_labels[i][j];
		    }
		    info_annotation_labels[i][j] = NULL;
		}
		delete[] info_annotation_labels[i];
	    }
	    info_annotation_labels[i] = NULL;

	    if(info_scores[i]){
		for(j=0;j<info_number_of_lines[i];j++){
		    if(info_scores[i][j]){
			for(k=0;k<info_number_of_annotation_labels;k++){
			    if(info_scores[i][j][k]) delete [] info_scores[i][j][k];
			    info_scores[i][j][k] = NULL;
			}
			delete [] info_scores[i][j];
		    }
		    info_scores[i][j] = NULL;
		}
		delete[] info_scores[i];
	    }
	    info_scores[i] = NULL;

	    if(start_positions[i]) delete[] start_positions[i];
	    start_positions[i] = NULL;
	    if(end_positions[i])  delete [] end_positions[i];
	    end_positions[i] = NULL;	 
	    
	    info_number_of_lines[i] = info_number_of_lines[i+1];

	    // reallocate memory
	    info_seq_names[i] = new char*[info_number_of_lines[i]];
	    for(j=0;j<info_number_of_lines[i];j++){
		info_seq_names[i][j] = new char[Max_word_length];
		strcpy(info_seq_names[i][j]," ");
	    }
	    
	    info_annotation_labels[i] = new bool**[info_number_of_lines[i]];
	    for(j=0; j<info_number_of_lines[i]; j++){
		info_annotation_labels[i][j] = new bool*[info_number_of_annotation_labels];
		for(k=0; k<info_number_of_annotation_labels; k++){
		    info_annotation_labels[i][j][k] = new bool[info_number_of_each_annotation_label[k]];
		    info_annotation_labels[i][j][k] = false;
		}
	    }
	    
	    info_scores[i] = new double**[info_number_of_lines[i]];
	    for(j=0; j<info_number_of_lines[i]; j++)
	    {
		info_scores[i][j] = new double*[info_number_of_annotation_labels];
		for(k=0; k<info_number_of_annotation_labels; k++)
		{
		    info_scores[i][j][k] = new double[info_number_of_each_annotation_label[k]];
		    info_scores[i][j][k] = 0;
		}
	    }
	    
	    start_positions[i] = new int[info_number_of_lines[i]];
	    end_positions[i]   = new int[info_number_of_lines[i]];

	    for(j=0; j<info_number_of_lines[i]; j++)
	    {
		start_positions[i][j] = 0;
		end_positions[i][j]   = 0;
	    }

	    // shift lines
	    for(j=0; j<info_number_of_lines[i]; j++)
	    {
		strcpy(info_seq_names[i][j], info_seq_names[i+1][j]);
		for(k=0;k<info_number_of_annotation_labels; k++)
		{
		    for(l=0;l<info_number_of_each_annotation_label[k]; l++)
		    {
			info_annotation_labels[i][j][k][l] = info_annotation_labels[i+1][j][k][l];
			info_scores[i][j][k][l]            = info_scores[i][j][k][l];
		    }
		}
		start_positions[i][j] = start_positions[i+1][j];
		end_positions[i][j]   = end_positions[i+1][j];
	    }    
	}
	if(info_number_of_sources-1>=0)
	{
	    // remove memory
	    if(info_seq_names[info_number_of_sources-1]){
		for(j=0; j<info_number_of_lines[info_number_of_sources-1]; j++){
		    if(info_seq_names[info_number_of_sources-1][j]) delete [] info_seq_names[info_number_of_sources-1][j];
		    info_seq_names[info_number_of_sources-1][j] = NULL;
		}
		delete[] info_seq_names[info_number_of_sources-1];
	    }
	    info_seq_names[info_number_of_sources-1] = NULL;

	    if(info_source_names[info_number_of_sources-1]) delete[] info_source_names[info_number_of_sources-1];
	    info_source_names[info_number_of_sources-1] = NULL;

	    if(info_annotation_labels[info_number_of_sources-1]){
		for(j=0;j<info_number_of_lines[info_number_of_sources-1];j++){
		    if(info_annotation_labels[info_number_of_sources-1][j]){
			for(k=0;k<info_number_of_annotation_labels;k++){
			    if(info_annotation_labels[info_number_of_sources-1][j][k]) delete [] info_annotation_labels[info_number_of_sources-1][j][k];
			    info_annotation_labels[info_number_of_sources-1][j][k] = NULL;
			}
			delete [] info_annotation_labels[info_number_of_sources-1][j];
		    }
		    info_annotation_labels[info_number_of_sources-1][j] = NULL;
		}
		delete[] info_annotation_labels[info_number_of_sources-1];
	    }
	    info_annotation_labels[info_number_of_sources-1] = NULL;

	    if(info_scores[info_number_of_sources-1]){
		for(j=0;j<info_number_of_lines[info_number_of_sources-1];j++){
		    if(info_scores[info_number_of_sources-1][j]){
			for(k=0;k<info_number_of_annotation_labels;k++){
			    if(info_scores[info_number_of_sources-1][j][k]) delete [] info_scores[info_number_of_sources-1][j][k];
			    info_scores[info_number_of_sources-1][j][k] = NULL;
			}
			delete [] info_scores[info_number_of_sources-1][j];
		    }
		    info_scores[info_number_of_sources-1][j] = NULL;
		}
		delete[] info_scores[info_number_of_sources-1];
	    }
	    info_scores[info_number_of_sources-1] = NULL;

	    if(start_positions[info_number_of_sources-1]) delete[] start_positions[info_number_of_sources-1];
	    start_positions[info_number_of_sources-1] = NULL;
	    if(end_positions[info_number_of_sources-1])  delete [] end_positions[info_number_of_sources-1];
	    end_positions[info_number_of_sources-1] = NULL;
	   
	}
	if(info_number_of_sources-1==0){
	    if(info_source_names) delete [] info_source_names;
	    info_source_names = NULL;

	    if(info_seq_names) delete [] info_seq_names;
	    info_seq_names = NULL;
	    
	    if(info_annotation_labels) delete [] info_annotation_labels;
	    info_annotation_labels = NULL;

	    if(info_scores) delete [] info_scores;
	    info_scores = NULL;
	    if(start_positions) delete [] start_positions;
	    start_positions = NULL;
	    if(end_positions) delete [] end_positions;
	    end_positions = NULL;
	 
	    if(info_number_of_lines) delete [] info_number_of_lines;
	    info_number_of_lines = NULL;
	    if(info_number_of_each_annotation_label) delete [] info_number_of_each_annotation_label;
	    info_number_of_each_annotation_label = NULL;
				
	}
	info_number_of_sources--;
    }
}

int Info::remove_line (const int source, const int line)
{
    int check = 0;
    
    if((source<0)||(source>=info_number_of_sources))
    {
	cout<<" ERROR : class Info::remove_line : source ("<<source<<") out of range [0,"
	    << info_number_of_sources-1 <<"].\n";
	check++;
    }else{	
	if ((line < 0) || (line >= info_number_of_lines[source]))
	{
	    cout << "ERROR : class Info::remove_line : line (" << line << ") out of range [0,"
		 << info_number_of_lines[source]-1 << "].\n";
	    check++;	
	}
    }
    
    if(!check){
	if (info_number_of_lines[source] > 0)
	{
	    int i = 0;
	    int j = 0;
	    int k = 0;

	    for (i = line; i<info_number_of_lines[source]-1; i++)
	    {
		strcpy(info_seq_names[source][i], info_seq_names[source][i+1]);

		for(j=0; j<info_number_of_annotation_labels;j++){
		    for(k=0; k<info_number_of_each_annotation_label[j]; k++){
			info_annotation_labels[source][i][j][k] = info_annotation_labels[source][i+1][j][k];
			info_scores[source][i][j][k]            = info_scores[source][i][j][k];
		    }
		}
		start_positions[source][i] = start_positions[source][i+1];  
		end_positions[source][i]   = end_positions[source][i+1];      		
	
         
	    }
	    
	    if (info_number_of_lines[source]-1 > 0)
	    {
		if (info_seq_names[source][info_number_of_lines[source]-1])  
		    delete [] info_seq_names[source][info_number_of_lines[source]-1];
		info_seq_names[source][info_number_of_lines[source]-1] = NULL;

		if (info_annotation_labels[source][info_number_of_lines[source]-1])
		{
		    for(i=0; i<info_number_of_annotation_labels; i++){
			if(info_annotation_labels[source][info_number_of_lines[source]-1][i])
			    delete [] info_annotation_labels[source][info_number_of_lines[source]-1][i];
			info_annotation_labels[source][info_number_of_lines[source]-1][i] = NULL;
			
		    }		    
		    delete [] info_annotation_labels[source][info_number_of_lines[source]-1];
		}
	    }
	    info_number_of_lines[source]--;
	    
	    if (info_number_of_lines[source]-1 == 0)
	    {
		this->remove_source(source);
	    }		    
	}
    }
    return(check);
}

int Info::preselect_lines(const int type, const int number_of_labels, int* const preselect_labels)
{
    int check = 0;
    
    if (number_of_labels <= 0)
    {
	cout << "ERROR: Info::preselect_lines: number_of_labels (" 
	     << number_of_labels << ").\n" << flush;
	check++;
    }
    if((type<0)||(type>=info_number_of_annotation_labels)){
	cout<< "ERROR: Info::preselect_lines: type : "<<type
	    <<"out of range [0.."<<info_number_of_annotation_labels<<"]\n.";
	check++;
    }
    if (preselect_labels == NULL)
    {
	cout << "ERROR: Info::preselect_lines: array preselect_labels is NULL.\n" << flush;
	check++;
    }

    if (check == 0)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	    
	int      found_label         = 0;
	
	for(i=0; i<info_number_of_sources; i++)
	{
	    if (info_number_of_lines[i] > 0)
	    {
		for (j=0; j<info_number_of_lines[i]; j++)
		{
		    found_label = 0;
	    
		    for (k=0; k<number_of_labels; k++)
		    {			
			if (info_annotation_labels[i][j][type][preselect_labels[k]] == true)
			{
			   
			    found_label++;
			    break;
			}
		    }
		    
		    
		    if (found_label == 0)
		    {
			this->remove_line(i,j);
			j--;
		
		    }
		}
	    }
	}
    }
    return(check);
}

int Info::preselect_lines(const int number_of_labels, char** const preselect_labels, model_parameters* const MP)
{
    int check = 0;
    
    if (number_of_labels <= 0)
    {
	cout << "ERROR: Info::preselect_lines: number_of_labels (" 
	     << number_of_labels << ").\n" << flush;
	check++;
    }

    if (preselect_labels == NULL)
    {
	cout << "ERROR: Info::preselect_lines: array preselect_labels is NULL.\n" << flush;
	check++;
    }

    if (check == 0)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int label_index = -1;
	int label_type = -1;
	int found_label = 0;
	    
	for(i=0; i<info_number_of_sources; i++)
	{
	    if (info_number_of_lines[i] > 0)
	    {
		for (j=0; j<info_number_of_lines[i]; j++)
		{
		    found_label = 0;
		    label_index = -1;
		    label_type = -1;
		    
		    for (k=0; k<number_of_labels; k++)
		    {		
			for(l=0;l<info_number_of_annotation_labels;l++){
			    label_index = convert_labelname_to_int(MP->get_Annotation_Label(),l,preselect_labels[k]);
			    if(label_index!=-1){
				label_type = l;
				break;
			    }
			}
			if(label_index!=-1){
			    if (info_annotation_labels[i][j][label_type][label_index] == true)
			    {
				found_label++;
				break;
			    }
			}
		    }		    
		    
		    if (found_label == 0)
		    {
			this->remove_line(i,j);
			j--;
		    }
		}
	    }
	}
    }
    return(check);
}

int Info::get_boundary_positions(int** const number_of_boundary_positions, 
				 int*** const boundary_positions) const
{
    int check = 0;

    int i = 0;
    int j = 0;
    
    if(*number_of_boundary_positions) delete[] (*number_of_boundary_positions);
    (*number_of_boundary_positions) = NULL;
	
    if(*boundary_positions){
	for(i=0; i<info_number_of_sources; i++){
	    if((*boundary_positions)[i]) delete[] (*boundary_positions)[i];
	    (*boundary_positions)[i] = NULL;
	}	       
	delete [] (*boundary_positions);
    }
    (*boundary_positions) = NULL;
    
    if (check == 0)
    {
	
	if(info_number_of_sources>0){
	    (*number_of_boundary_positions) = new int[info_number_of_sources];
	    (*boundary_positions) = new int*[info_number_of_sources];
	}
	
	for(i=0; i<info_number_of_sources; i++){	
	    (*number_of_boundary_positions) = 0;
	    
	    if (info_number_of_lines[i] > 0)
	    {
		
		int* start_positions_copy = new int[info_number_of_lines[i]];
		for (j=0; j<info_number_of_lines[i]; j++)
		{
		    start_positions_copy[j] = start_positions[i][j];
		}
		
		int* end_positions_plus_one = new int[info_number_of_lines[i]];
		for (j=0; j<info_number_of_lines[i]; j++)
		{
		    end_positions_plus_one[j] = end_positions[i][j]+1;
		}
		
		int  length_array = 0;
		int* array        = NULL;

		check += sort_elements_of_array_by_increasing_order(info_number_of_lines[i],
								    start_positions_copy);
		
		if (check != 0)
		{
		    cout << "ERROR: Info::get_boundary_positions: occurred in function "
			 << "sort_elements_of_array_by_increasing_order for array start_positions_copy.\n";
		}
	    
		if (check == 0)
		{
		    check += sort_elements_of_array_by_increasing_order(info_number_of_lines[i],
									end_positions_plus_one);
		    
		    if (check != 0)
		    {
			cout << "ERROR: Info::get_boundary_positions: occurred in function "
			     << "sort_elements_of_array_by_increasing_order for array end_positions_plus_one.\n";
		    }
		}
		
		if (check == 0)
		{
		    check += sort_elements_of_two_arrays_and_create_new_array(// input
			info_number_of_lines[i],
			start_positions_copy,
			info_number_of_lines[i],
			end_positions_plus_one,
			// output
			&length_array,
			&array);
		    if (check != 0)
		    {
			cout << "ERROR: Info::get_boundary_positions: occurred in function "
			     << "sort_elements_of_two_arrays_and_create_new_array.\n";
		    }
		}
		
		if (check == 0)
		{
		    check += make_elements_of_array_unique(// input and output
			&length_array,
			array);
		
		    if (check != 0)
		    {
			cout << "ERROR: Info::get_boundary_positions: occurred in function "
			     << "make_elements_of_array_unique.\n";
		    }
		    
		}
		
		if (check == 0)
		{
		    (*number_of_boundary_positions)[i] = length_array;
		    
		    if (length_array > 0)
		    {
			(*boundary_positions)[i] = new int[length_array];
			
			for (j=0; j<length_array; j++)
			{
			    (*boundary_positions)[i][j] = array[j];
			}
		    }
		}
	    
		if (start_positions_copy) delete [] start_positions_copy;
		start_positions_copy = NULL;
		
		if (end_positions_plus_one) delete [] end_positions_plus_one;
		end_positions_plus_one = NULL;
		
		if (array) delete [] array;
		array = NULL;
	    }
	}
	if(!check){
	    if(*number_of_boundary_positions) delete [] (*number_of_boundary_positions);
	    (*number_of_boundary_positions) = NULL;
	    if(*boundary_positions){
		for(i=0; i<info_number_of_sources; i++){
		    if((*boundary_positions)[i]) delete [] (*boundary_positions)[i];
		    (*boundary_positions)[i] = NULL;
		}
		delete [] (*boundary_positions);
	    }
	    (*boundary_positions) = NULL;
	}	    
    }
    return(check);
}

int Info::get_details_for_position(// input
    const Sequence* const seq,
    const int             seq_position,
    model_parameters* const  MP,
    // output
    bool****              annotation_label,
    double****            score
    ) const
{
    int check = 0;

    if (seq->get_ac() == NULL)
    {
	cout << "ERROR: Info::get_details_for_position: ac of Sequence seq is NULL.\n" << flush;
	check++;
    }
    
    // initialise the output variables
    
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    bool score_label_set = false;

    bool* set_score = new bool[info_number_of_annotation_labels];

    for(i=0; i<info_number_of_annotation_labels; i++){
	set_score[i] = false;
    }

    if(!(*annotation_label))
    {
	(*annotation_label) = new bool**[info_number_of_sources];
	for(i=0; i<info_number_of_sources; i++)
	{
	    (*annotation_label)[i] = new bool*[info_number_of_annotation_labels];
	    for(j=0; j<info_number_of_annotation_labels; j++){
		(*annotation_label)[i][j] = new bool[info_number_of_each_annotation_label[j]];
		for(k=0;k<info_number_of_each_annotation_label[j];k++){
		    (*annotation_label)[i][j][k] = false;
		}
	    }
	}
    }

    if(!(*score)){
	(*score) = new double **[info_number_of_sources];
	for(i=0;i<info_number_of_sources;i++)
	{
	    (*score)[i] = new double *[info_number_of_annotation_labels];
	    for(j=0; j<info_number_of_annotation_labels; j++)
	    {
		(*score)[i][j] = new double[info_number_of_each_annotation_label[j]];
		for(k=0; k<info_number_of_each_annotation_label[j]; k++)
		{
		    (*score)[i][j][k]  = 0;
		}
	    }
	}
    }

    if ((check == 0) && (info_number_of_sources > 0))
    {

	bool* set_score = new bool[info_number_of_annotation_labels];
	int set_phase = -1;
	int ori_phase = -1;
	
	// select lines which contain the desired position for each source
	for(i=0; i<info_number_of_sources; i++)
	{
	    if(check){
		break;
	    }
	    ori_phase = -1;
	    set_phase = -1;
	    
	    for(j=0; j<info_number_of_annotation_labels; j++){
		set_score[j] = false;
	    }
	    
	    if(info_number_of_lines[i]>0)
	    {
		for (j=0; j<info_number_of_lines[i]; j++)
		{
		    if(check){
			break;
		    }
		    if ((cmp_nocase(info_seq_names[i][j], seq->get_ac()) == 0) &&
			(start_positions[i][j] <= seq_position)           && 
			(end_positions[i][j]   >= seq_position))
		    {
			
			for(k=0; k<info_number_of_annotation_labels; k++)
			{
			    if(check){
				break;
			    }
			    score_label_set = MP->get_score_of_Annotation_Label(k);
			    for(l=0; l<info_number_of_each_annotation_label[k]; l++)
			    {
				if(info_annotation_labels[i][j][k][l]){			
				    if((score_label_set==true) &&
				       ((*annotation_label)[i][k][l]==true)&&
				       ((*score)[i][k][l]!=info_scores[i][j][k][l]))
				    {
					cout<<"ERROR: class Info: get_details_for_position, contradiction label on line : "<<j
					    <<" position : "<<seq_position<<" label("<<MP->get_Annotation_Label_setname(k)<<","
					    <<MP->get_Annotation_Label_name(k,l)<<")"<<endl;
					cout<<"original score : "<<(*score)[i][k][l]<<", new score : "<<info_scores[i][j][k][l]<<endl;
					check++;
					break;
				    }
				    else{					    
					(*annotation_label)[i][k][l] = info_annotation_labels[i][j][k][l];
					(*score)[i][k][l]            = info_scores[i][j][k][l];
					set_score[k] = true;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	if(set_score) delete [] set_score;
	set_score = NULL;
    }
    if (check != 0)
    {
	for(i=0;i<info_number_of_sources;i++){
	    for(j=0; j<info_number_of_annotation_labels; j++){
		for(k=0; k<info_number_of_each_annotation_label[j]; k++){
		    (*annotation_label)[i][j][k] = false;
		    (*score)[i][j][k] = 0;
		}
	    }
	}
    }
    return(check);
}

int Info::get_details_for_next_label(// input
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
    ) const
{
    // note : direction == 1  (use direction == seq->get_orientation()), 
    //        direction == -1 (use direction == - seq->get_orientation()), 

    int check = 0;

    if (seq->get_ac() == NULL)
    {
	cout << "ERROR: Info::get_details_for_next_label: ac of Sequence seq is NULL.\n" << flush;
	check++;
    }
  
    if ((direction != 1) && (direction != -1))
    {
	cout << "ERROR: Info::get_details_for_next_label: direction (" 
	     << direction << ") has to be 1 or -1.\n" << flush;
	check++;
    }

    if(seq_annotation_label==NULL)
    {
	cout << "ERROR: Info::get_details_for_next_label: seq_annotation_label is NULL.\n" << flush;
	check++;
    }

    if(seq_score==NULL)
    {
	cout << "ERROR: Info::get_details_for_next_label: seq_score is NULL.\n" << flush;
	check++;
    }
    // initialise the output variables
  
    int i = 0;
    int j = 0;
    int k = 0;
  
    (*position)       = 0;
   
    if(!(*annotation_label))
    {
	(*annotation_label) = new bool**[info_number_of_sources];
	for(i=0; i<info_number_of_sources; i++)
	{
	    (*annotation_label)[i] = new bool*[info_number_of_annotation_labels];
	    for(j=0; j<info_number_of_annotation_labels; j++){
		(*annotation_label)[i][j] = new bool[info_number_of_each_annotation_label[j]];
		for(k=0;k<info_number_of_each_annotation_label[j];k++){
		    (*annotation_label)[i][j][k] = false;
		}
	    }
	}
    }

    if(!(*score)){
	(*score) = new double **[info_number_of_sources];
	for(i=0;i<info_number_of_sources;i++)
	{
	    (*score)[i] = new double *[info_number_of_annotation_labels];
	    for(j=0; j<info_number_of_annotation_labels; j++)
	    {
		(*score)[i][j] = new double[info_number_of_each_annotation_label[j]];
		for(k=0; k<info_number_of_each_annotation_label[j]; k++)
		{
		    (*score)[i][j][k]  = 0;
		}
	    }
	}
    }
    if ((check == 0) && (info_number_of_sources > 0))
    {
	int orientation = direction;
	int pos         = 0;
	int start_pos   = 0;
	int end_pos     = 0;
	
	if (orientation == 1) {
	    start_pos = seq->get_start_position();
	    end_pos   = seq->get_end_position();
	}
	else if (orientation == -1) {
	    start_pos = seq->get_end_position();
	    end_pos   = seq->get_start_position();
	}

	for (pos = seq_position + orientation; pos != end_pos + orientation; pos += orientation) 
	{

	    for(i=0; i<info_number_of_sources; i++)
	    {
		for(j=0; j<info_number_of_annotation_labels; j++)
		{
		    for(k=0; k<info_number_of_each_annotation_label[j]; k++)
		    {
			(*annotation_label)[i][j][k] = false;
			(*score)[i][j][k] = 0;
		    }
		}
	    }
	    (*position) = pos;
	    
	    check += this->get_details_for_position(// input
		seq,
		pos, //abs coordinates 
		MP,
		// output
		annotation_label,
		score
		);
	  	    
	    bool next_label = false;

	    // if any label from any source is different from previous label,
	    // we consider as different, even same label with different score
	    for(i=0; i<info_number_of_sources; i++)
	    {
		if(next_label)
		{
		    break;
		}
		for(j=0; j<info_number_of_annotation_labels; j++)
		{
		    if(next_label)
		    {
			break;
		    }
		    for(k=0; k<info_number_of_each_annotation_label[j]; k++)
		    {			
			if (((*annotation_label)[i][j][k] != seq_annotation_label[i][j][k])
			    ||((*score)[i][j][k]!=seq_score[i][j][k]))
			{
			    next_label = true;
			    break;
			}
		    }
		}
	    }   	   
	    if(next_label)
	    {
		break;
	    }
	}
	
	// if loop is empty take info from seq_position
	
	if (seq_position == end_pos) {

	    (*position) = seq_position;
	    
	    check += this->get_details_for_position(// input
		seq,
		seq_position, //abs coordinates 
		MP,
		// output
		annotation_label,
		score
		);
	}

    }
    if (check != 0) {

	(*position)       = 0;

	for(i=0;i<info_number_of_sources;i++){
	    for(j=0; j<info_number_of_annotation_labels; j++){
		for(k=0; k<info_number_of_each_annotation_label[j]; k++){
		    (*annotation_label)[i][j][k] = false;
		    (*score)[i][j][k] = 0;
		}
	    }
	}
    }
    return(check);
}

int Info::get_annotation_of_sequence(// input
    const Sequence* const seq,
    model_parameters* const MP,
    // output
    bool*****              seq_annotation_labels,
    double*****            seq_scores
    ) const
{
    // note : the output arrays have an entry for every position in the sequence
    //        and have the orientation of the sequence (i.e. start codon comes before stop codon)
    
    int check = 0;
    
    if (seq->get_ac() == NULL)
    {
	cout << "ERROR: Info::get_annotation_of_sequence: ac of Sequence seq is NULL.\n" << flush;
	check++;
    }
    if (seq->length() == 0)
    {
	cout << "ERROR: Info::get_annotation_of_sequence: length of sequence is 0.\n" << flush;
	check++;
    }
 
    if((*seq_annotation_labels)!=NULL)
    {
	cout << "ERROR: Info::get_annotation_of_sequence: seq_annotation_labels not NULL.\n" << flush;
	check++;
    }
    if((*seq_scores)!=NULL)
    {
	cout << "ERROR: Info::get_annotation_of_sequence: seq_scores not NULL.\n" << flush;
	check++;
    }
    if ((check == 0) && (info_number_of_sources > 0))
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int length = seq->length();
	// initialise the output variables
	
	(*seq_annotation_labels) = new bool   ***[length];
	(*seq_scores)            = new double ***[length];
	for(i=0;i<length;i++)
	{
	    (*seq_annotation_labels)[i] = new bool   **[info_number_of_sources];
	    (*seq_scores)[i]            = new double **[info_number_of_sources];
	    for(j=0; j<info_number_of_sources; j++){
		(*seq_annotation_labels)[i][j] = new bool   *[info_number_of_annotation_labels];
		(*seq_scores)[i][j]                = new double *[info_number_of_annotation_labels];
		for(k=0; k<info_number_of_annotation_labels; k++){
		    (*seq_annotation_labels)[i][j][k] = new bool   [info_number_of_each_annotation_label[k]];
		    (*seq_scores)[i][j][k]            = new double [info_number_of_each_annotation_label[k]];
		    for(l=0; l<info_number_of_each_annotation_label[k]; l++){
			(*seq_annotation_labels)[i][j][k][l] = false;
			(*seq_scores)[i][j][k][l] = 0;
		    }
		}
	    }
	}

	int pos            = 0;
	int same_label_pos = 0;
	int start_pos      = 0;
	int end_pos        = 0;
	int started        = 0;
	int next_start_pos = 0;
	
	bool*** annotation_label = new bool**[info_number_of_sources];
	double*** score = new double ** [info_number_of_sources];
	for(i=0;i<info_number_of_sources;i++){
	    annotation_label[i] = new bool * [info_number_of_annotation_labels];
	    score[i] = new double * [info_number_of_annotation_labels];
	    for(j=0; j<info_number_of_annotation_labels; j++){
		annotation_label[i][j] = new bool[info_number_of_each_annotation_label[j]];
		score[i][j] = new double[info_number_of_each_annotation_label[j]];
		for(k=0;k<info_number_of_each_annotation_label[j]; k++){
		    annotation_label[i][j][k] = false;
		    score[i][j][k] = 0;
		}
	    }
	}
	bool*** tmp_annotation_label = new bool**[info_number_of_sources];
	double*** tmp_score = new double ** [info_number_of_sources];
	for(i=0;i<info_number_of_sources;i++){
	    tmp_annotation_label[i] = new bool * [info_number_of_annotation_labels];
	    tmp_score[i] = new double * [info_number_of_annotation_labels];
	    for(j=0; j<info_number_of_annotation_labels; j++){
		tmp_annotation_label[i][j] = new bool[info_number_of_each_annotation_label[j]];
		tmp_score[i][j] = new double[info_number_of_each_annotation_label[j]];
		for(k=0;k<info_number_of_each_annotation_label[j]; k++){
		    tmp_annotation_label[i][j][k] = false;
		    tmp_score[i][j][k] = 0;
		}
	    }
	}
	
	int      next_position       = 0;
            
	start_pos = seq->get_start_position();
	end_pos   = seq->get_end_position();

	pos = start_pos;
	
	while (pos != end_pos + 1) 
	{      
	    if (check == 0) 
	    {
		for(i=0; i<info_number_of_sources; i++){
		    for(j=0;j<info_number_of_annotation_labels; j++){
			for(k=0; k<info_number_of_each_annotation_label[j]; k++){
			    annotation_label[i][j][k] = false;
			    score[i][j][k] = 0;
			}
		    }
		}
		check += this->get_details_for_position(// input
		    seq,
		    pos, //abs coordinates 
		    MP,
		    // output
		    &annotation_label,
		    &score
		    );
     
		if (check != 0) 
		{
		    cout << "ERROR: Info::get_annotation_of_sequence: error in function "
			 << "Info::get_details_for_position for position (" << pos 
			 << ").\n";
		}
		
		
		if (check == 0) 
		{
		
		    for(i=0; i<info_number_of_sources; i++){
			for(j=0;j<info_number_of_annotation_labels; j++){
			    for(k=0; k<info_number_of_each_annotation_label[j]; k++){
				tmp_annotation_label[i][j][k] = false;
				tmp_score[i][j][k] = 0;
			    }
			}
		    }
		    check += this->get_details_for_next_label(// input
			seq,
			1,   // get next label
			pos,
			annotation_label,
			score,
			MP,
			// output
			&next_position,
			&tmp_annotation_label,
			&tmp_score
			);
		  
		    if (check != 0) {
			cout << "ERROR: Info::get_annotation_of_sequence: error in function "
			     << "Info::get_details_for_next_label for direction 1 and position (" << pos 
			     << ").\n";
		    }
		}//check == 0, get_detail_for_next_label, with orientation == 1
		

		// loop over all positions until next label start and set both label and phase		
		if (next_position == pos) {
		    next_start_pos = pos + 1;
		}
		else {
		    next_start_pos = next_position;
		}
	
		for (same_label_pos = pos; same_label_pos != next_start_pos; same_label_pos += 1) 
		{

		    // get annotation labels
		    for(i=0;i<info_number_of_sources;i++){
			for(j=0;j<info_number_of_annotation_labels;j++){
			    for(k=0;k<info_number_of_each_annotation_label[j];k++){
				(*seq_annotation_labels)[abs(same_label_pos-start_pos)][i][j][k] = annotation_label[i][j][k];
				(*seq_scores)[abs(same_label_pos-start_pos)][i][j][k] = score[i][j][k];
			    }
			}
		    }
		} // after same label pos;
		pos += 1;
	    } // if check ==0, at the beginning of while loos
	    
	    // *********************************************************************
	    
	    if (check != 0) 
	    {
		break;
	    }
	} // for-loop over pos
	
	if (check != 0) 
	{	            
	    if(*seq_annotation_labels){
		for(i=0;i<length;i++)
		{
		    if((*seq_annotation_labels)[i]){
			for(j=0; j<info_number_of_sources; j++){
			    if((*seq_annotation_labels)[i][j]){
				for(k=0; k<info_number_of_annotation_labels; k++){
				    if((*seq_annotation_labels)[i][j][k]) delete[] (*seq_annotation_labels)[i][j][k];
				    (*seq_annotation_labels)[i][j][k] = NULL;
				}
				delete [] (*seq_annotation_labels)[i][j];
			    }
			    (*seq_annotation_labels)[i][j];
			}
			delete [] (*seq_annotation_labels)[i];
		    }
		    (*seq_annotation_labels)[i] = NULL;
		}
		delete[] (*seq_annotation_labels);
	    }
	    (*seq_annotation_labels) = NULL;

	    if(*seq_scores){
		for(i=0;i<length;i++)
		{
		    if((*seq_scores)[i]){
			for(j=0; j<info_number_of_sources; j++){
			    if((*seq_scores)[i][j]){
				for(k=0; k<info_number_of_annotation_labels; k++){
				    if((*seq_scores)[i][j][k]) delete[] (*seq_scores)[i][j][k];
				    (*seq_scores)[i][j][k] = NULL;
				}
				delete [] (*seq_scores)[i][j];
			    }
			    (*seq_scores)[i][j];
			}
			delete [] (*seq_scores)[i];
		    }
		    (*seq_scores)[i] = NULL;
		}
		delete[] (*seq_scores);
	    }
	    (*seq_scores) = NULL;
	}// for check != 0	 
	
	if(annotation_label){
	    for(i=0;i<info_number_of_sources;i++)
	    {
		if((annotation_label)[i]){
		    for(j=0; j<info_number_of_annotation_labels; j++){
			if(annotation_label[i][j]) delete[] annotation_label[i][j];
			annotation_label[i][j] = NULL;
		    }
		    delete [] annotation_label[i];
		}
		annotation_label[i];
	    }
	    delete [] annotation_label;
	}
	annotation_label = NULL;
	
	if(score)
	{
	    for(i=0;i<info_number_of_sources;i++)
	    {
		if(score[i]){
		    for(j=0; j<info_number_of_annotation_labels; j++){
			if(score[i][j]) delete [] score[i][j];
			score[i][j] = NULL;
		    }			    
		    delete [] score[i];
		}
		score[i] = NULL;
	    }
	    delete[] score;
	}
	score = NULL;

	if(tmp_annotation_label){
	    for(i=0;i<info_number_of_sources;i++)
	    {
		if(tmp_annotation_label[i]){
		    for(j=0; j<info_number_of_annotation_labels; j++){
			if(tmp_annotation_label[i][j]) delete[] tmp_annotation_label[i][j];
			tmp_annotation_label[i][j] = NULL;
		    }
		    delete [] tmp_annotation_label[i];
		}
		tmp_annotation_label[i];
	    }
	    delete [] tmp_annotation_label;
	}
	tmp_annotation_label = NULL;
	
	if(tmp_score){
	    for(i=0;i<info_number_of_sources;i++)
	    {
		if(tmp_score[i]){
		    for(j=0; j<info_number_of_annotation_labels; j++){
			if(tmp_score[i][j]) delete [] tmp_score[i][j];
			tmp_score[i][j] = NULL;
		    }			    
		    delete [] tmp_score[i];
		}
		tmp_score[i] = NULL;
	    }
	    delete[] tmp_score;
	}
	tmp_score = NULL;
    }

    return (check);
}
  
int Info::get_length_of_output(void) const
{
    if (output != NULL)
    {
	return(static_cast<int>(strlen(output)));
    }
    else
    {
	return(0);
    }
}

int Info::get_line_of_output(// input
    const int start_of_line,
    const int length_of_line,
    // output
    char* line) const
{
    // line must have memory to store length_of_line+1 letters
    
    int check = 0;
    int length_of_output = this->get_length_of_output();
    
    if ((start_of_line < 0) || (start_of_line > length_of_output))
    {
	cout << "ERROR: Info::get_line_of_output: start_of_line (" << start_of_line 
	     << ") out of range [0, " << length_of_output << "].\n";
	check++;
    }
    if ((length_of_line) < 0)
    {
	cout << "ERROR: Info::get_line_of_output: length_of_line (" << length_of_line << ") <= 0.\n";
	check++;
    }
    if (line == NULL)
    {
	cout << "ERROR: Info::get_line_of_output: line is NULL.\n";
	check++;
    }
    
    if (check == 0)
    {
	check += get_substr(start_of_line, length_of_line, 
			    output, line);
    }
    return(check);
}

int Info::print_source (const int source,
			model_parameters* const MP, 
			std::ostream &o) const
{
    int i = 0;
    int check = 0;

    if ((source < 0) || (source > info_number_of_sources-1))
    {
	cout << "ERROR : Info::print_source : source (" << source 
	     << ") out of range [0,"<<info_number_of_sources-1<<"].\n"<<flush;
	check++;
    }

    for (i=0; i<info_number_of_lines[source]; i++)
    {
	o<<"source : "<<info_source_names[source]<<endl;
	this->print_line(source,i,MP, o);
    }
    return check;
}

int Info::print_line(const int source,
		     const int line, 
		     model_parameters* const MP,
		     std::ostream &o) const
{
    int check = 0;
    int i = 0;
    int j = 0;
  
    if ((source < 0) || (source > info_number_of_sources-1))
    {
	cout << "ERROR : Info::print_line : source (" << source 
	     << ") out of range [0,"<<info_number_of_sources-1<<"].\n"<<flush;
	check++;
    }else{
	if((line<0)||(line > info_number_of_lines[source])){
	    cout <<" ERROR : Info:: print_line : line ("<<line
		 <<") out of range [0,"<<info_number_of_lines[source]<<"].\n"<<flush;
	    check++;
	}
    }
    if (o == NULL)
    {
	cout << "ERROR : Info::print_line : std::ostream o is NULL.\n";
	check++;
    }
    
    if (check == 0)
    {
	o << "[" << line 
	  << "]\t : name " << info_seq_names[source][line] 
	  <<"|\t"<<endl;
	for(i=0; i<info_number_of_annotation_labels; i++)
	{	    
	    o<<MP->get_Annotation_Label_setname(i)<<" : "<<endl;
	    for(j=0; j<info_number_of_each_annotation_label[i]; j++)
	    {
		if(info_annotation_labels[source][line][i][j]==true)
		{
		    o<<MP->get_Annotation_Label_name(i,j)<<" : "<<info_scores[source][line][i][j]<<endl;
		}
	    }
	}
	o<< "\t| " << start_positions[source][line] 
	 << "\t -> | " << end_positions[source][line] 
	 <<"|\t"; 
	o<<"\n";
	
    }
    return check;
}
 
void Info::print (model_parameters* const MP, std::ostream &o) const
{
    int i = 0;

    for (i=0; i<info_number_of_sources; i++)
    {
	this->print_source(i,MP, o);
    }
    return;
}
