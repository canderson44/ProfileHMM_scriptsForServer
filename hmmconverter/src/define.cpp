/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define constants, typedefs and enums used in the pairhmm, sequence classes 
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/define.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "define.h"
#include <ctype.h>
#include <math.h>

using namespace std;
        
int NumOfChar(const char* const instr)
{
    if(!instr){
	cout<<"The input string is a null string! "<<endl;
	return -1;
    }
    int length = strlen(instr);
    if(length>Max_word_length){
	 cout<<"Input String has length exceed "<<Max_word_length<<" return -1"<<endl;
	 return -1;
    }
    return length;
}

char* Nullstrcpy(const char* str2, int& check)
{
    char* str1 = NULL;
    int length = 0;

    if(!str2){
	cout<<"Error: Nullstrcpy::input str2 is NULL !"<<endl;
	check++;
    }
    
    if(check==0){
	length = strlen(str2);
	str1 = new char[length+1];
	strcpy(str1,str2);
	str1[length] = '\0';
    }
    return str1;
}

int break_id_into_tokens(const char* const id, char*** tokens, int& num_of_token){
    int check = 0;
    int length = NumOfChar(id);
    if(length<=0){
	cout<<"ERROR:: Break id into tokens: input id : "<<id<<endl;
	check++;
    }

    if(!check){
	int i = 0;
	int j = 0;
	int count = 0;
	int first = 0;
	int last = 0;
	int* charnum = new int[MaxTokens];
	for(i= 0; i<length; i++){
	    if(id[i]=='.'){
		charnum[count]=last-first;
		count++;
		first = last;
	    }else{
		last++;
	    }
	}	
	charnum[count]=last-first;
	(*tokens) = new char*[count+1];
	num_of_token = count+1;
	int tmp = 0;
	// the actual count is count+1
	for(i=0;i<=count;i++){
	    (*tokens)[i]=new char[charnum[i]+1];
	    for(j=0;j<charnum[i];j++){
		(*tokens)[i][j]=id[tmp];
		tmp++;
	    }
	    if(id[tmp]=='.'){
		tmp++;
	    }
	    (*tokens)[i][j]='\0';
	    if(tmp>length){
		cout<<"ERROR: Break_id_into_tokens : input id: "<<id<<endl;
		check++;
		break;
	    }
	}
	if(charnum) delete[] charnum;
	charnum=NULL;
    }
    return check;
}

void convert_to_base_and_invert(int number, const int base, array<int>* index)
{

    int coefficient=0;
    int value_of_coefficient=0;
    
    // initialise index with zeros
    
    for (int i=0; i<index->GetDimension(0); i++)
    {index->SetElement(i, 0);}
    
    while (number>0) {
	value_of_coefficient=(int)((float)(number)/(powf((float)(base), (float)(coefficient))))%base;
	number-=value_of_coefficient* (int)(powf((float)(base), (float)(coefficient)));
	index->SetElement(coefficient, value_of_coefficient);
	coefficient++;
    }
    return;
}

array<int> get_indices(int linear_index, int* dim, int d, int& check)
{
    int length=dim[0];
    int i = 0;
    for (i=1; i<d; i++)
    {
	length*=dim[i];
    }
   
    if (!dim)
    {
	cout << "ERROR: get_indices : dim array does not exist.\n";
	check++;
    }
    if (d<1)
    {
	cout << "ERROR: get_indices : number of dimensions d " << d << "<1.\n";
	check++;
    }
    for (i=0; i<d; i++)
    { 
	if (!dim[i]) {
	    cout << "ERROR: get_indices : dim i = " << i << " does not exist.\n";
	    check++;
	    break;
	}
	if (dim[i]==0) {
	    cout << "ERROR: get_indices : dim i = " << i << " = 0.\n";
	    check++;
	    break;
	}
    }
    if ((linear_index<0) || (linear_index>(length-1)))
    {
	cout << "ERROR: get_indices : linear_index " << linear_index << " out of range (0.." << length-1 << ").\n";
	check++;
    }
    
    array<int> indices;
    if(!check){
	indices.SetNumberofDimensions(1);
	indices.SetDimension(0,d);
  
	int index;
	int last_linear_index;
	
	for (int j=d-1; j>-1; j--) {
	    last_linear_index=linear_index;
	    
	    linear_index=static_cast<int>(floor( (static_cast<float>(linear_index)/(static_cast<float>(dim[j]))) ) );
	    index=0;
	    index=static_cast<int>(nint(static_cast<float>(dim[j])*
					((static_cast<float>(last_linear_index)/(static_cast<float>(dim[j])))-
					 (static_cast<float>(linear_index))
					    )
				       )
		);
	    indices.SetElement(j,index);
	}
    }
    return(indices);
}

int generate_random_letter(const int alphabet) //modified
{
    double random = rand();
    return (int)floor((double)(alphabet*random));
}

int get_substr(const int start, const int length, 
	       const char* const string, char* const substring)
{
    // substring must have memory to store length+1 letters

    int check = 0;
    
    if (string == NULL){
	cout << "get_substr : string is NULL/\n";
	check++;
    }
    if (substring == NULL){
	cout << "get_substr : substring is NULL/\n";
	check++;
    }
    if (start < 0){
	cout << "get_substr : start (" << start << ") is < 0.\n";
	check++;
    }
    if (length <= 0){
	cout << "get_substr : length (" << length << ") is <= 0.\n";
	check++;
    }

    if (check == 0) {

	int i = 0;
	int length_of_string = static_cast<int>(strlen(string)); 
// length of string excluding '\0' at end
	
	for (i=start; i< start+min(length_of_string-start, length)+1; i++) {
	    substring[i-start] = string[i];
	}
	substring[min(length_of_string-start,length)] = '\0';
    }
    return(check);
}

Tube<int> convert_boundaries_to_tube(const int length_x,
				     const int length_y,
				     vector<vector<int> > & upper_bounds, 
				     vector<vector<int> > & lower_bounds) {
  
    // returns empty tube if upper and lower bounds are empty

    int check = 0;

    Tube<int> return_tube(length_x, 2);
    
    if (length_x <= 0) {
	cout << "ERROR: convert_boundaries_to_tube: length_x (" << length_x << ") <= 0.\n";
	check++;
    }
    if (length_y <= 0) {
	cout << "ERROR: convert_boundaries_to_tube: length_y (" << length_y << ") <= 0.\n";
	check++;
    }
    if ((check == 0) && (! (upper_bounds.empty() && lower_bounds.empty()))) {
	
	const int x = 0;
	const int y = 1;
	
	int i, j;
	
	// check if upper bounds are ordered both in x and y if they are in range [0, length_x-1] or [0, length_y-1]
	
	for (i=1; i<upper_bounds.size(); i++) {
	    
	    if (! ((upper_bounds[i][x] > upper_bounds[i-1][x]) &&
		   (upper_bounds[i][y] >= upper_bounds[i-1][y]))) {
		
		cout << "ERROR: convert_boundaries_to_tube: upper_bounds[" << i << "] ([x],[y]) (" 
		     << upper_bounds[i][x] << ", " << upper_bounds[i][y] << ") not > (x) / >= (y)  [" << i-1 << "] ([x],[y]) (" 
		     << upper_bounds[i-1][x] << ", " << upper_bounds[i-1][y] << ")\n";
		check++;
	    }
	    if (! ((upper_bounds[i][x] >= 0)  && (upper_bounds[i][x] <= length_x))) {
		
		cout << "ERROR: convert_boundaries_to_tube: upper_bounds[" << i << "][x] ("
		     << upper_bounds[i][x] << ") out of range [0, " << length_x << "].\n";
		check++;
	    }
	    if (! ((upper_bounds[i][y] >= 0)  && (upper_bounds[i][y] <= length_y))) {
	
		cout << "ERROR: convert_boundaries_to_tube: upper_bounds[" << i << "][y] ("
		     << upper_bounds[i][y] << ") out of range [0, " << length_y << "].\n";
		check++;
	    }
	}
	
	// check if lower bounds are ordered both in x and y if they are in range
	
	for (i=1; i<lower_bounds.size(); i++) {
	    
	    if (! ((lower_bounds[i][x] > lower_bounds[i-1][x]) &&
		   (lower_bounds[i][y] >= lower_bounds[i-1][y]))) {
		
		cout << "ERROR: convert_boundaries_to_tube: lower_bounds[" << i << "] ([x],[y]) (" 
		     << lower_bounds[i][x] << ", " << lower_bounds[i][y] << ") not > (x) / >= (y) [" << i-1 << "] ([x],[y]) (" 
		     << lower_bounds[i-1][x] << ", " << lower_bounds[i-1][y] << ")\n";
		check++;
	    }
	    if (! ((lower_bounds[i][x] >= 0)  && (lower_bounds[i][x] <= length_x))) {
		
		cout << "ERROR: convert_boundaries_to_tube: lower_bounds[" << i << "][x] ("
		     << lower_bounds[i][x] << ") out of range [0, " << length_x << "].\n";
		check++;
	    }
	    if (! ((lower_bounds[i][y] >= 0)  && (lower_bounds[i][y] <= length_y))) {
		
		cout << "ERROR: convert_boundaries_to_tube: lower_bounds[" << i << "][y] ("
		     << lower_bounds[i][y] << ") out of range [0, " << length_y << "].\n";
		check++;
	    }
	}
    
	// get tube and make sure it anables to go from (0,0) to (length_x, length_y)
	
	if (check == 0) {
	    
	    Tube<int> tube(length_x, 2); 
	    
	    int index_l = 0;
	    int index_u = 0;
	    
	    int lower   = 0;
	    int lower_x = lower_bounds[index_l][x];
	    int upper   = upper_bounds[index_u][y];
	    int upper_x = upper_bounds[index_u][x];
	    
	    for (i=0; i<length_x; i++) {
		
		tube.SetElement(i, 0, lower); // lower bound at position x along the x sequence
		tube.SetElement(i, 1, upper); // upper bound at position x along the x sequence
		
		// check that this interval i is not empty (minimum length 1)
		
		if (lower > upper) {
		    
		    cout << "ERROR: convert_boundaries_to_tube: tube[" << i << "][0] ("
			 << lower << ") >= tube[" << i << "][1] (" << upper << ") => cannot .\n";
		    check++;
		}
		
		// check that this interval i and the previous interval i-1 have at least pointlike overlap
		// and that previous interval has "smaller coordinates" than this one
		
		if (i > 0) {
		    
		    int previous_lower = tube.GetElement(i-1, 0);
		    int previous_upper = tube.GetElement(i-1, 1);
		    
		    if ((previous_upper < lower) || (previous_lower > upper)) {

			cout << "ERROR: convert_boundaries_to_tube: tube[" << i << "] interval (" << lower
			     << ", " << upper << ") does not have correct overlap with previous interval [" << i-1 << "] (" 
			     << previous_lower << ", " << previous_upper << "). Cannot go from (0, 0) to (" << length_x 
			     << ", " << length_y << ").\n";
			check++;
		    }
		    else {
			if (previous_upper > upper) {
			    tube.SetElement(i-1, 1, upper);
			}
			if (lower < previous_lower) {
			    lower = previous_lower;
			    tube.SetElement(i, 0, lower);
			}
		    }
		}
		
		
		if (i < length_x) 
		{
		    if (i == upper_x) {
			
			index_u++;
			
			if (index_u < upper_bounds.size()) {
			    
			    upper   = upper_bounds[index_u][y];
			    upper_x = upper_bounds[index_u][x];
			}
			else {
			    
			    upper   = length_y-1;
			    upper_x = length_x-1;
			}
		    }
		    if (i == (lower_x-1)) {
			
			index_l++;
			
			if (index_l < lower_bounds.size()) {
			    
			    lower   = lower_bounds[index_l-1][y];
			    lower_x = lower_bounds[index_l][x];
			}
			else {
			    
			    if ((index_l-1) < lower_bounds.size()) {
				
				lower   = lower_bounds[index_l-1][y];
				lower_x = length_x-1;
			    }
			}
		    }
		} // if i < length_x		
		tube.SetElement(i, 1, upper); // upper bound at position x along the x sequence
	    } // i loop [0, length_x-1]

	    if (check == 0) {
		
		return_tube = tube;
	    }
	} // if check == 0
    } // if input check == 0 and upper and lower not empty
    
    return return_tube;
}

double log_Stirling(const int n) { // Stirling approximation of log(n!), see Bronstein page 103

    double log_n_fak = 0.0;
    double d_n       = (double)(n);

    if (n > 0) {
	log_n_fak = ((d_n + 0.5) * log(d_n)) - d_n + (0.5 * log(2.0 * dPi));
    }
    else {
	cout << "ERROR: log_Stirling: n (" << n << ") <= 0.\n";
    }
    return(log_n_fak);
}

int bitshift(const int value)
{
    // note: function returns positive value (>0) if successful, else returns negative value (<0)
    
    int check=0;

    if (value < 0) {
	cout << "ERROR: bitshift : cannot calculate bitshift for value (" << value
	     << ") < 0 .\n";
	check--;
    }
    if (check == 0) {

	int base     = 2;
	int exponent = 1;
	int number   = base;
	
	while (value > number) {
	    number *= base;
	    exponent++;
	}
	return(exponent);
    }
    else {
	return(check);
    }
}

double calculate_number_from_score(const Score  score_value)
{
    double return_value = 0.0;

    return_value = pow(Base, score_value);
    
    return(return_value);
}

Score calculate_score(const double number)
{
    Score return_score = Logzero;

    if (number > 0)
    {
	return_score = log(number)/log(Base);
    }
    else // if number == 0
    {
	return_score = Logzero;
    }

    return(return_score);
}

int set_executable_path(const char* const new_executable_path) 
{
    int check = 0;
    
    if (new_executable_path == NULL) {
	cout << "ERROR: set_executable_path : new_executable_path is NULL.\n" << flush;
	check++;
    }
    if (executable_path != NULL) {
	cout << "ERROR: set_executable_path : executable_path has already been set.\n" << flush;;
	check++;
    }
    
    if (check == 0) {
	
	executable_path = new char[Max_line_length];
	strcpy(executable_path, new_executable_path);       
	
    }
    return(check);
}

int set_temporary_path(const char* const new_temporary_path) 
{
    int check = 0;

    if (new_temporary_path == NULL) {
	cout << "ERROR: set_temporary_path : new_temporary_path is NULL.\n" << flush;
	check++;
    }
    if (temporary_path != NULL) {
	cout << "ERROR: set_temporary_path : temporary_path has already been set.\n" << flush;
	check++;
    }
    
    if (check == 0) {
	
	temporary_path = new char[Max_line_length];
	strcpy(temporary_path, new_temporary_path);
	
#ifdef _PRINT
	cout << "set_temporary_path: temporary_path = " << temporary_path << "\n";
#endif
	
    }
    return(check);
}


void convert_to_base(int number, const int base, array<int>* index)
{
    
    int length_of_index=index->GetDimension(0);
    int coefficient=0;
    int value_of_coefficient=0;
    
    // initialise index with zeros
    
    for (int i=0; i<index->GetDimension(0); i++)
    {index->SetElement(i, 0);}
    
    while (number>0) {
	value_of_coefficient=(int)((float)(number)/(powf((float)(base), (float)(coefficient))))%base;
	number-=value_of_coefficient* (int)(powf((float)(base), (float)(coefficient)));
	index->SetElement(length_of_index-1 - coefficient, value_of_coefficient);
	coefficient++;
    }
    return;
}


int compare_state_type(int delta_x, int delta_y, int delta_x2, int delta_y2)
{
    int state_type = 0;
    int state_type2 = 0;
    if((delta_x==0)&&(delta_y==0)){
	state_type = 0;
    }else if((delta_x>0)&&(delta_y==0)){
	state_type = 1;
    }else if((delta_x==0)&&(delta_y>0)){
	state_type = 2;
    }else if((delta_x>0)&&(delta_y>0)){
	state_type = 3;
    }else{
	state_type = -1;
    }

    if((delta_x2==0)&&(delta_y2==0)){
	state_type2 = 0;
    }else if((delta_x2>0)&&(delta_y2==0)){
	state_type2 = 1;
    }else if((delta_x2==0)&&(delta_y2>0)){
	state_type2 = 2;
    }else if((delta_x2>0)&&(delta_y2>0)){
	state_type2 = 3;
    }else{
	state_type2 = -1;
    }

    if(state_type==state_type2){
	return 1;
    }else{
	return 0;
    }
}

int get_n_chars(const int n, const char letter, char** word)
{
    int check = 0;

    if ((*word) != NULL)
    {
	cout << "ERROR: get_n_chars : word has to be NULL.\n";
	check++;
    }
    if (check == 0) {
	if (n > 0) {
	    int i = 0;
	    (*word) = new char[n+1];
	    
	    for (i=0; i<n; i++)
	    {(*word)[i] = letter;}
	    (*word)[n] = '\0';
	}
	else {
	    (*word) = new char[1];
	    (*word)[0] = '\0';
	}
    }
    return(check);
}

int splitstring(const char* instr, int& NoOfItems, char*** Items, char splitchar)
{
   
    int i=0;
    int check = 0;
    NoOfItems = 0;
    char splitchar2;
    if(splitchar==' ')
    {
	splitchar2 = '\t';
    }else{
	splitchar2 = splitchar;
    }

    if(!instr)
    {
	cout<<"Error : splitstring: input instr is NULL!.\n";
	check++;
	return check;
    }
    while(instr[i]!='\0')
    {
	if((instr[i]==splitchar)||(instr[i]=='\n')||(instr[i]==splitchar2)){
	    if((i-1)>=0)
	    {
		if((instr[i-1]!=splitchar)&&(instr[i-1]!='\n')&&(instr[i-1]!=splitchar2))
		{
		    NoOfItems++;
		}
	    }
	}
	i++;
    }
    if(instr[i]=='\0'){
	if((i-1)>=0)
	{
	    if((instr[i-1]!=splitchar)&&(instr[i-1]!='\n')&&(instr[i-1]!=splitchar2))
	    {
		NoOfItems++;
	    }
	}
    }
   
    i= 0 ;
    int j=0;
    int start_index=0;
    while(instr[i]!='\0')
    {
	if((instr[i]==splitchar)||(instr[i]=='\n')||(instr[i]==splitchar2)){
	    if((i-1)>=0)
	    {
		if((instr[i-1]!=splitchar)&&(instr[i-1]!='\n')&&(instr[i-1]!=splitchar2))
		{
		    int length = i-start_index;
		    
		    for(int k=0; k<length; k++)
		    {
			(*Items)[j][k] = instr[k+start_index];
		    }
		    (*Items)[j][length] = '\0';
		    j++;
		    start_index = i+1;
		}else{
		    start_index++;
		}
	    }
	}
	i++;
    }
    if(instr[i]=='\0'){
	if((i-1)>=0)
	{
	    if((instr[i-1]!=splitchar)&&(instr[i-1]!='\n')&&(instr[i-1]!=splitchar2))
	    {
		int length = i-start_index;
	
		for(int k=0; k<length; k++)
		{
		    (*Items)[j][k] = instr[k+start_index];
		}
		(*Items)[j][length] = '\0';
	    }
	}
    }
    return(check);
}

int cmp_nocase(const char* s1, const char* s2)
{
    int check = 0;

    int length_1 = strlen(s1);
    int length_2 = strlen(s2);
    
    if (length_1 == length_2) {
	for (int i=0; i<length_1; i++) {
	    if (toupper(s1[i]) != toupper(s2[i])) {
		check++;
		break;
	    }
	}
    }
    else {
	check++;
    }
    return(check);
}

int sort_elements_of_array_by_increasing_order(const int length_array,
					       int* const array)
{
    int check = 0;
    
    if ((array != NULL) && (length_array > 0)) {
	
	int i = 0;
	int j = 0;
	
	int min       = 0;
	int min_index = 0;
	int copy      = 0;
	
	for (i=0; i<length_array-1; i++) {
	    
	    min       = array[i];
	    min_index = i;
	    
	    for (j=i+1; j<length_array; j++) {
		if (array[j] < min) {
		    min       = array[j];
		    min_index = j;
		}
	    }
	    if (min_index != i) {
		copy             = array[min_index];
		array[min_index] = array[i];
		array[i]         = copy;
	    }
	}
    }
    return(check);
}

int sort_elements_of_two_arrays_and_create_new_array(// input
						     const int        length_array_1,
						     const int* const array_1,
						     const int        length_array_2,
						     const int* const array_2,
						     // output
						     int*             length_array,
						     int**      const array)
{
    // note : - if one elements appears several times it will appear the same number of times in
    //          the output array
    //        - array_1 and array_2 have to be already sorted by increasing numbers and the 
    //          the output array is also sorted by increasing numbers
    
    int check = 0;
    
    if (((array_1 == NULL) || (length_array_1 <= 0)) &&
	((array_2 == NULL) || (length_array_2 <= 0))) {
	cout << "ERROR : sort_elements_of_two_arrays_and_create_new_array : array_1 and array_2 are either NULL "
	     << "or have length <= 0.\n";
	check++;
    }
    else {
	if ((array_1 != NULL) && (length_array_1 > 0)) {
	    int i = 0;
	    for (i=0; i<length_array_1-1; i++) {
		if (array_1[i+1] < array_1[i]) {
		    cout << "ERROR : sort_elements_of_two_arrays_and_create_new_array : array_1 does not "
			 << "contain increasing numbers.\n";
		    check++;
		    break;
		}
	    }
	}
	if ((array_2 != NULL) && (length_array_2 > 0)) {
	    int i = 0;
	    for (i=0; i<length_array_2-1; i++) {
		if (array_2[i+1] < array_2[i]) {
		    cout << "ERROR : sort_elements_of_two_arrays_and_create_new_array : array_2 does not "
			 << "contain increasing numbers.\n";
		    check++;
		    break;
		}
	    }
	}
    }
    if ((*array) != NULL) {
	cout << "ERROR : sort_elements_of_two_arrays_and_create_new_array : array has to be NULL.\n";
	check++;
    }
    if (check == 0) {
	
	// get elements of output array
	
	int i = 0;
	int j = 0;
	
	int pos_1 = 0;
	int pos_2 = 0;
	int pos   = 0;
	
	// allocate memory for output array, initialise things
	
	if (((array_1 != NULL) && (length_array_1 > 0)) ||
	    ((array_2 != NULL) && (length_array_2 > 0))) {
	    
	    (*array) = new int[length_array_1 + length_array_2];
	}
	if (((array_1 != NULL) && (length_array_1 > 0)) &&
	    ((array_2 != NULL) && (length_array_2 > 0))) {
	    
	START_LOOP:
	    
	    for (i=pos_1; i<length_array_1; i++) {
		for (j=pos_2; j<length_array_2; j++) {
		    if (array_1[i] <= array_2[j]) {
			
			(*array)[pos] = array_1[i];
			pos_1 = i+1;
			pos_2 = j;
		    }
		    else {
			(*array)[pos] = array_2[j];
			pos_1 = i;
			pos_2 = j+1;
		    }
		    pos++;
		    goto START_LOOP;
		}
	    }
	    // check if there are any elements left
	  
	    if ((pos_1 < length_array_1) &&
		(pos_2 < length_array_2)) {
		cout << "ERROR : sort_elements_of_two_arrays_and_create_new_array\n";
	    }
	    else {
		if (pos_1 < length_array_1) {
		    for (i=pos_1; i<length_array_1; i++) {
			(*array)[pos] = array_1[i];
			pos++;
		    }
		}
		else if (pos_2 < length_array_2) {
		    for (i=pos_2; i<length_array_2; i++) {
			(*array)[pos] = array_2[i];
			pos++;
		    }
		}
	    }
	    (*length_array) = pos;
	}
	else if (((array_1 != NULL) && (length_array_1 > 0)) &&
		 ((array_2 == NULL) || (length_array_2 <= 0))) {
	    for (i=0; i<length_array_1; i++) {
		(*array)[pos] = array_1[i];
		pos++;
	    }
	    (*length_array) = pos;
	}
	else if (((array_1 == NULL) || (length_array_1 <= 0)) &&
		 ((array_2 != NULL) && (length_array_2 >  0))) {
	    for (i=0; i<length_array_2; i++) {
		(*array)[pos] = array_2[i];
		pos++;
	    }
	    (*length_array) = pos;
	}
    }
    return(check);
}

int make_elements_of_array_unique(int* length_array,
				  int* array)
{
    int check = 0;

    if ((array != NULL) && ((*length_array) > 0)) {
	int i = 0;
	for (i=0; i<(*length_array)-1; i++) {
	    if (array[i+1] < array[i]) {
		cout << "ERROR : make_elements_of_array_unique : array does not \n"
		     << "contain increasing numbers.\n";
		check++;
		break;
	    }
	}
    }
    
    if ((array != NULL) && ((*length_array) > 0) && (check == 0)) {
	
	int i = 0;
	int j = 0;
	int k = 0;
	
	int start_index = 0;
	
	for (i=0; i<(*length_array); i++) {
	    
	    start_index = i+1;
	    
	RESTART:
	    
	    for (j=start_index; j<(*length_array); j++) {
		if (array[j] == array[i]) {
		    for (k=j+1; k<(*length_array); k++)
		    {array[k-1] = array[k];}
		    
		    (*length_array)--;
		    start_index = j;
		    
		    goto RESTART;
		}
	    }
	}
	if ((*length_array) == 0) {
	    if (array) delete [] array;
	    array = NULL;
	}
    }
    return(check);
}

int swap_string(char** a, char** b){
    int check = 0;
    char* tmp_buffer = new char[Max_word_length];
    if((!*a)||(!*b)){
	check++;
	return check;
    }
    strcpy(tmp_buffer,(*a));
    strcpy((*a),(*b));
    strcpy((*b),tmp_buffer);
    if(tmp_buffer) delete[] tmp_buffer;
    tmp_buffer = NULL;
    
    return check;

}

int check_undefine(bool* const labels,
		   int   length)
{
    int i = 0;
    for(i=0; i<length; i++)
    {
	if(labels[i])
	{
	    return false;
	}
    }
    return true;
}

bool get_tag(const char* line,
	     const char* tag)
{
    if(!line)
    {
	cout<<"Error : get_tag : input line is NULL."<<endl;
	return false;
    }
    if(!tag)
    {
	cout<<"Error : get_tag : input tag is NULL."<<endl;
	return false;
    }

    int line_length = strlen(line);
    int tag_length = strlen(tag);

    int i =0;
    while(i<line_length)
    {
	if(line[i]=='<')
	{
	    bool match = true;
	    int j = 0;
	    i++;
	    while((i<line_length)&&(j<tag_length))
	    {
		if(line[i]!=tag[j])
		{
		    match =false;
		    break;
		}
		i++;
		j++;
	    }
	    if(match)
	    {
		return true;
	    }else{
		return false;
	    }
	}
	i++;
    }
    return false;
}
