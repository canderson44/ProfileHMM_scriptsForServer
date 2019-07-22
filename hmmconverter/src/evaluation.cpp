#include <ctype.h>
#include "evaluation.h"

 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
 */

/*
  RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/evaluation.cpp,v 1.3 2008/12/14 10:39:23 natural Exp $
*/

int get_details_for_position(// input
			     const char*  const seq_ac,
			     const int          seq_position,
			     model_parameters* const MP,
			     //			     
			     const int          number_of_sources,      
			     int*         const number_of_lines,
			     char***      const seq_names,
			     char**       const source_names,
			     bool****     const annotation_labels,
			     double****   const scores,			 
			     int**        const start_positions,
			     int**        const end_positions,
			     // output
			     bool****     const annotation_label,
			     double****   const score			  
			     )
{
    int check = 0;
    
    if (seq_ac == NULL)
    {
	cout << "ERROR: get_details_for_position: seq_ac is NULL.\n" << flush;
	check++;
    }
    
    if (number_of_sources <= 0)
    {
	cout << "ERROR: get_details_for_position: number_of_sources (" << number_of_sources
	     << ") must be > 0.\n" << flush;
	check++;
    }
    if (number_of_lines == NULL)
    {
	cout << "ERROR: get_details_for_position: number_of_lines are NULL.\n"<<flush;
	check++;
    }
    if (seq_names == NULL)
    {
	cout << "ERROR: get_details_for_position: seq_names are NULL.\n" << flush;
	check++;
    }
    if (source_names == NULL)
    {
	cout << "ERROR: get_details_for_position: source_names are NULL.\n" << flush;
	check++;
    }
    if (annotation_labels == NULL)
    {
	cout << "ERROR: get_details_for_position: annotation_labels are NULL.\n" << flush;
	
	check++;
    }
    if (scores == NULL)
    {
	cout << "ERROR: get_details_for_position: scores are NULL.\n" << flush;
	
	check++;
    }

    if (start_positions == NULL)
    {
	cout << "ERROR: get_details_for_position: start_positions are NULL.\n" << flush;
	check++;
    }
    if (end_positions == NULL)
    {
	cout << "ERROR: get_details_for_position: end_positions are NULL.\n" << flush;
	check++;
    }

    if((*annotation_label)!=NULL)
    {
	cout<<"ERROR: get_details_for_position: annotaiton_label should be NULL.\n"<<flush;
	check++;
    }
    if((*score)!=NULL)
    {
	cout<<"ERROR: get_details_for_position: score should be NULL.\n"<<flush;
	check++;
    }
    
    // initialise the output variables
    
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    
    int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
    int* number_of_each_annotation_label = new int[number_of_annotation_labels];
    for(i=0; i<number_of_annotation_labels; i++)
    {
	number_of_each_annotation_label[i] = MP->get_Annotation_Label_size(i);
    }

    bool score_label_set = false;
  
    (*annotation_label) = new bool**[number_of_sources];
    for(i=0; i<number_of_sources; i++)
    {
	(*annotation_label)[i] = new bool*[number_of_annotation_labels];
	for(j=0; j<number_of_annotation_labels; j++){
	    (*annotation_label)[i][j] = new bool[MP->get_Annotation_Label_size(j)];
	    for(k=0;k<MP->get_Annotation_Label_size(j);k++){
		(*annotation_label)[i][j][k] = false;
	    }
	}
    }

    if(!(*score)){
	(*score) = new double **[number_of_sources];
	for(i=0;i<number_of_sources;i++)
	{
	    (*score)[i] = new double *[number_of_annotation_labels];
	    for(j=0; j<number_of_annotation_labels; j++)
	    {
		(*score)[i][j] = new double[MP->get_Annotation_Label_size(j)];
		for(k=0; k<MP->get_Annotation_Label_size(j); k++)
		{
		    (*score)[i][j][k]  = 0;
		}
	    }
	}
    }

    if ((check == 0)&&(number_of_sources>0))
    {
	bool* set_score = new bool[number_of_annotation_labels];

	// select lines which contain the desired position
	for(i=0; i<number_of_sources; i++)
	{
	    if(check)
	    {
		break;
	    }

	    if(number_of_lines[i] > 0)
	    {
		for (j=0; j<number_of_lines[i]; j++)
		{
		    if(check){
			break;
		    }

		    if ((cmp_nocase(seq_names[i][j], seq_ac) == 0) &&
			(start_positions[i][j] <= seq_position)    && 
			(end_positions[i][j]   >= seq_position))
		    {
	
			for(k=0; k<number_of_annotation_labels; k++)
			{
			    set_score[k] = false;
			    if(check){
				break;
			    }
			    score_label_set = MP->get_score_of_Annotation_Label(k);
			    for(l=0; l<number_of_each_annotation_label[k]; l++)
			    {
				if(annotation_labels[i][j][k][l]){
				    set_score[k] = true;
				    if((score_label_set == true) &&
				       ((*annotation_label)[i][k][l]==true)&&
				       ((*score)[i][k][l]==scores[i][j][k][l]))
				    {
					cout<<"ERROR: get_details_for_position, contradiction label on line : "<<j
					    <<" position : "<<seq_position<<" label : "<<MP->get_Annotation_Label_name(j,l)<<endl;
					check++;
					break;
				    }else if((score_label_set == false) &&
					     (set_score[k] == true))
				    {
					cout<<"ERROR: get_details_for_position, contradiction label on line : "<<j
					    <<" position : "<<seq_position<<" label : "<<MP->get_Annotation_Label_name(k,l)<<endl;
					check++;
					break;
				    }else{
					(*annotation_label)[i][k][l] = annotation_labels[i][j][k][l];
					(*score)[i][k][l]            = scores[i][j][k][l];
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
	for(i=0;i<number_of_sources;i++){
	    for(j=0; j<number_of_annotation_labels; j++){
		for(k=0; k<number_of_each_annotation_label[j]; k++){
		    (*annotation_label)[i][j][k] = false;
		    (*score)[i][j][k] = 0;
		}
	    }
	}
    }
    return(check);
}

int preselect_lines(// input
		    const int             number_of_labels,
		    const int             label_type,
		    int*  const           preselect_labels,
		    model_parameters* const  MP,
		    // input and output
		    int            number_of_sources,
		    int*     const number_of_lines,
		    char***  const seq_names,
		    char**   const source_names,
		    bool**** const annotation_labels,
		    int**    const start_positions,
		    int**    const end_positions
		    ) // optional
{
    int check = 0;

    if (number_of_labels <= 0)
    {
	cout << "ERROR: preselect_lines: number_of_labels (" 
	     << number_of_labels << ").\n" << flush;
	check++;
    }

    if ((label_type<0)||(label_type>=MP->get_Total_Number_of_Annotation_Labels()))
    {
	cout << "ERROR: preselect_lines: label_type (" 
	     << label_type << ") out of range[0,"
	     << MP->get_Total_Number_of_Annotation_Labels()<<"].\n" << flush;
	check++;
    }

    if (preselect_labels == NULL)
    {
	cout << "ERROR: preselect_lines: labels is NULL.\n" << flush;
	check++;
    }
    
    if (number_of_sources <= 0)
    {
	cout << "ERROR: preselect_lines: number_of_sources (" << number_of_sources
	     << ") must be > 0.\n" << flush;
	check++;
    }

    if(number_of_lines==NULL)
    {
	cout<<"ERROR: preselect_lines: number_of_lines are NULL.\n"<<flush;
	check++;
    }
    
    if (seq_names == NULL)
    {
	cout << "ERROR: preselect_lines: seq_names are NULL.\n" << flush;
	check++;
    }
    if (source_names == NULL)
    {
	cout << "ERROR: preselect_lines: source_names are NULL.\n" << flush;
	check++;
    }
  
    if (annotation_labels == NULL)
    {
	cout << "ERROR: preselect_lines: annotation_labels are NULL.\n" << flush;
	check++;
    }

    if (start_positions == NULL)
    {
	cout << "ERROR: preselect_lines: start_positions are NULL.\n" << flush;
	check++;
    }
    if (end_positions == NULL)
    {
	cout << "ERROR: preselect_lines: end_positions are NULL.\n" << flush;
	check++;
    }

    if (check == 0)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;
	
	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
	int      found_label         = 0;

	for (i=0; i<number_of_sources; i++ )
	{
	    if(number_of_lines[i] > 0)
	    {
		for(j=0; j<number_of_lines[i]; j++)
		{
		    found_label = 0;
		    
		    for (k=0; k<number_of_labels; k++)
		    {
			if (annotation_labels[i][j][label_type][preselect_labels[k]] == true)
			{
			    found_label++;
			    break;
			}
		    }	    	    
		    if (found_label == 0) 
		    {

			for (k=j+1; k< number_of_lines[i]; k++)
			{

			    strcpy(seq_names[i][k-1], seq_names[i][k]);          
			    // no need to move source
			    for(l=0; l<number_of_annotation_labels; l++)
			    {
				for(m=0; m<MP->get_Annotation_Label_size(l); m++)
				{
				    annotation_labels[i][k-1][l][m]= annotation_labels[i][k][l][m];
				}
			    }
			    			    
			    start_positions[i][k-1] = start_positions[i][k];
			    end_positions[i][k-1]   = end_positions[i][k];  
			}			
		    }
		    number_of_lines[i]--;
		    j--;
		}
	    }
	}
    }
    return(check);
}

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
    )
{
    int check = 0;

    if (number_of_labels <= 0)
    {
	cout << "ERROR: preselect_lines: number_of_labels (" 
	     << number_of_labels << ").\n" << flush;
	check++;
    }

    if (preselect_labels == NULL)
    {
	cout << "ERROR: preselect_lines: labels is NULL.\n" << flush;
	check++;
    }
    
    if (number_of_sources <= 0)
    {
	cout << "ERROR: preselect_lines: number_of_sources (" << number_of_sources
	     << ") must be > 0.\n" << flush;
	check++;
    }

    if(number_of_lines==NULL)
    {
	cout<<"ERROR: preselect_lines: number_of_lines are NULL.\n"<<flush;
	check++;
    }
    
    if (seq_names == NULL)
    {
	cout << "ERROR: preselect_lines: seq_names are NULL.\n" << flush;
	check++;
    }
    if (source_names == NULL)
    {
	cout << "ERROR: preselect_lines: source_names are NULL.\n" << flush;
	check++;
    }  
    if (annotation_labels == NULL)
    {
	cout << "ERROR: preselect_lines: annotation_labels are NULL.\n" << flush;
	check++;
    }
    if (start_positions == NULL)
    {
	cout << "ERROR: preselect_lines: start_positions are NULL.\n" << flush;
	check++;
    }
    if (end_positions == NULL)
    {
	cout << "ERROR: preselect_lines: end_positions are NULL.\n" << flush;
	check++;
    }
  
    if (check == 0)
    {
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;
	int label_index = -1;
	int label_type = -1;
	
	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
	int      found_label         = 0;

	for (i=0; i<number_of_sources; i++)
	{
	    if(number_of_lines[i] > 0)
	    {
		for(j=0; j<number_of_lines[i]; j++)
		{
		    found_label = 0;
		    label_index = -1;
		    label_type = -1;
		    
		    for (k=0; k<number_of_labels; k++)
		    {
			for(l=0;l<number_of_annotation_labels;l++){
			    label_index = convert_labelname_to_int(MP->get_Annotation_Label(),l,preselect_labels[k]);
			    if(label_index!=-1)
			    {
				label_type = l;
				break;
			    }
			}      	
			if(label_index!=-1)
			{
			    if (annotation_labels[i][j][label_type][label_index] == true)
			    {
				found_label++;
				break;
			    }
			}
		    }	    
	    
		    if (found_label == 0) 
		    {
			for (k=j+1; k< number_of_lines[i]; k++)			
			{
			    strcpy(seq_names[i][k-1], seq_names[i][k]);          
			    // no need to move source
			    for(l=0; l<number_of_annotation_labels; l++)
			    {
				for(m=0; m<MP->get_Annotation_Label_size(l); m++)
				{
				    annotation_labels[i][k-1][l][m]= annotation_labels[i][k][l][m];
				}
			    }		       

			    start_positions[i][k-1] = start_positions[i][k];
			    end_positions[i][k-1]   = end_positions[i][k]; 
			}
			
		    }
		    number_of_lines[i]--;
		    j--;
		}
	    }
	}
    }
    return(check);
}

int read_special_file_for_sequence(// input
    const char*     const file_name,
    const char*     const seq_ac,
    const int             seq_start_position,
    const int             seq_end_position,
    const int             position_offset,
    model_parameters*  const MP,
    // output
    int*        const number_of_allocated_lines,
    int*        const number_of_sources,
    int**       const number_of_lines,
    char****    const seq_names,
    char***     const source_names,
    bool*****   const annotation_labels,
    double***** const scores,
    int***      const start_positions,
    int***      const end_positions
    )
{
    // note: - variables that will be used to store output information have to be NULL when calling this function
    //       - upon successfull call of this function the memory of the output variables has to be 
    //         released within the calling program
    //       - gtf-comment lines start with '//' or '#'
    //       - it is assumed that there will be at most (length of sequence)/2 gtf lines for a sequence
    
    int check=0;
    
    // check input 
    
    if (file_name == NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: file_name is NULL.\n" << flush;
	check++;
    }
    if (seq_ac == NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: seq_ac is NULL.\n" << flush;
	check++;
    }
    if (seq_start_position > seq_end_position)
    {
	cout << "ERROR: read_special_file_for_sequence: seq_start_position (" << seq_start_position 
	     << ") > seq_end_position (" << seq_end_position << ").\n" << flush;
	check++;
    }
    
    // check that variables for output are NULL
    
    if ((*seq_names) != NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: seq_names should be NULL.\n" << flush;
	check++;
    }
    if ((*source_names) != NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: source_names should be NULL.\n" << flush;
	check++;
    }
    if ((*annotation_labels) != NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: annotation_labels should be NULL.\n" << flush;
	check++;
    }
    if ((*scores) != NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: scores should be NULL.\n" << flush;
	check++;
    }
 
    if ((*start_positions) != NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: start_positions should be NULL.\n" << flush;
	check++;
    }
    if ((*end_positions) != NULL)
    {
	cout << "ERROR: read_special_file_for_sequence: end_positions should be NULL.\n" << flush;
	check++;
    }
    // open file
    
    FILE* special_file = fopen(file_name, "rt");
    
    if(!special_file)
    {
	cout<<"ERROR: read_special_file_for_sequence: cannot read special file : "
	    <<file_name<<".\n"<<flush;
	check++;
    }
        
    if (check == 0)    
    {		
	// variables

	int max_number_of_lines = static_cast<int>(static_cast<float>(seq_end_position-seq_start_position+1)/2.);

	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();


	int i =0;
	int j =0;

	int NoOfItems = 0;
	char** Items = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	int NoOfItemsForLabel = 0;
	char** ItemsForLabel = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    ItemsForLabel[i] = new char[Max_word_length];
	    strcpy(ItemsForLabel[i]," ");
	}

	int NoOfItemsForScore = 0;
	char** ItemsForScore = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    ItemsForScore[i] = new char[Max_word_length];
	    strcpy(ItemsForScore[i]," ");
	}
	
	char** line_seq_names = NULL;
	bool*** line_annotation_labels = NULL;
	double*** line_scores = NULL;

	int*   line_start_positions = NULL;
	int*   line_end_positions = NULL;

	check+= initialize_labels(max_number_of_lines,
				  MP,
				  &line_seq_names,
				  &line_annotation_labels,
				  &line_scores,
				  //&line_other_labels,
				  &line_start_positions,
				  &line_end_positions);

	if(check)
	{
	    cout<<"ERROR: read_special_label_for_sequence, initialize_labels for line "<<endl;
	}

	char** line_source_names       = new char*[max_number_of_lines];
	for(i=0; i<max_number_of_lines; i++)
	{
	    line_source_names[i] = new char[Max_word_length];
	    strcpy(line_source_names[i]," ");
	}
	
	char* line                     = new char[Max_line_length];
	int   feature_start_position   = 0;
	
	// information on this sequence
	
	char*  seq_name_constraint       = new char[Max_word_length];
	strcpy(seq_name_constraint, seq_ac);
	int    start_position_constraint = seq_start_position;
	int    end_position_constraint   = seq_end_position;
	
	int    number_of_line            = 0;
	
	bool valid = false;
	
	// allocate memory for output 
	
	// need at most (length of sequence)/2 gtf lines for sequence
	
	while ((!feof(special_file)) && (check == 0))
	{
	    // initialise variables

	    strcpy(line, " "); 
	    
	    strcpy(line_seq_names[number_of_line], " ");
	    strcpy(line_source_names[number_of_line], " ");
	    line_start_positions[number_of_line] = 0;
	    line_end_positions[number_of_line] = 0;
	    
	    for(i=0;i<number_of_annotation_labels;i++)
	    {
		for(j=0;j<MP->get_Annotation_Label_size(i); j++)
		{
		    line_scores[number_of_line][i][j] = 0;
		}
	    }
	    
	    for(i=0; i<number_of_annotation_labels; i++)
	    {
		for(j=0; j<MP->get_Annotation_Label_size(i); j++)
		{
		    line_annotation_labels[number_of_line][i][j] = false;
		}
	    }       

	    // read one line if either max_line_length-1 characters were read or end of line was encountered	    
	    fgets(line, Max_line_length-1, special_file);
  
	    // decide of which type this line is 1.) comment, 2.) gtf_format line, else (is skipped)
	    check+=  splitstring(static_cast<const char*>(line), NoOfItems, &Items,' ');

	    valid = false;
	    
	    if(check){
		cout<<"ERROR: read_special_file_for_sequence: splitstring for line : "<<line<<"."<<endl;
		break;
	    }
	    
	    if (!((strcmp(Items[0], "//") == 0) || (strcmp(Items[0], "#") == 0)))	    
	    {
		int found_keyword  = 0;
		
		int tmp_pos = 0;
		
		tmp_pos = 4;
		
		int tmp_label_index = -1;
		int feature_label_size = 0;
		// keep reading the annotation label until it ends
	      
		while(1){
		   
		    if((tmp_pos>=NoOfItems)||(found_keyword)){
			break;
		    }

		    // end of annotation labels
		    if(!strcmp(Items[tmp_pos],"|")){
			break;
		    }
		    
		    check+=  splitstring(Items[tmp_pos], NoOfItemsForLabel, &ItemsForLabel,':');
		    
		    if(check){
			cout<<"ERROR: read_special_file_for_sequence: splitstring for Items["<<tmp_pos
			    <<"] : "<<Items[tmp_pos]<<"."<<endl;
			break;
		    }

		    if(NoOfItemsForLabel<2)
		    {
			cout<<"ERROR: read_special_file_for_sequence: line : "<<number_of_line<<", feature label : "<<Items[tmp_pos]<<" not in a right format."<<endl;
			check++;
			break;
		    }
		    
		    tmp_label_index = convert_typename_to_int(MP->get_Annotation_Label(),ItemsForLabel[0],number_of_annotation_labels);
		    
		    if(tmp_label_index==-1){
			tmp_pos++;
			continue;
		    }else{
			feature_label_size =MP->get_Annotation_Label_size(tmp_label_index);
			for(i=0; i<feature_label_size; i++){
			    if(cmp_nocase(ItemsForLabel[1],MP->get_Annotation_Label_name(tmp_label_index,i))==0)
			    {
				found_keyword++;
				break;
			    }			    
			}
			if(found_keyword==0){
			    tmp_pos++;
			}
		    }
		}

		if (found_keyword != 0)
		{	    		   
		    // read information from special_format line
		    
		    // get seq_name
		    strcpy(line_seq_names[number_of_line],Items[0]);
		    
		    // get source_name
		    strcpy(line_source_names[number_of_line],Items[1]);
		    
		    // get start position
		    line_start_positions[number_of_line] = atoi(Items[2]);
		    
		    // get end position
		    line_end_positions[number_of_line] = atoi(Items[3]);

		    // check if this is a line in valid gtf_format
		    
		    if ((line_start_positions[number_of_line] <= line_end_positions[number_of_line])&&(!check))
		    {

			feature_start_position = 0;
		
			feature_start_position = line_start_positions[number_of_line] + position_offset;
		       
			if ((cmp_nocase(line_seq_names[number_of_line], seq_name_constraint) == 0)      &&
			    //
			    // if there is an overlap between the constraint interval and the interval
			    // defined by the line
			    //
			    ( ! ((feature_start_position > end_position_constraint) ||
				 (line_end_positions[number_of_line]   < start_position_constraint))))
			{
			    // prune start or end position to start or end constraints, if necessary
			    
			    if ((feature_start_position <  start_position_constraint) &&
				(line_end_positions[number_of_line]   <= end_position_constraint))
			    {
				line_start_positions[number_of_line] = start_position_constraint - position_offset;

			    }
			    else if ((feature_start_position >=  start_position_constraint) &&
				     (line_end_positions[number_of_line]   >   end_position_constraint))
			    {
				line_end_positions[number_of_line] = end_position_constraint + position_offset;
		      
			    }

			    valid = true;			    
			}
		    }
		    else
		    {
			cout << "   ERROR: read_special_file_for_sequence: special_line is corrupted, stop "
			     << "reading this special_file, return empty arrays.\n" << flush;
			
			cout << "   special_line[" << number_of_line << "] :\n" << flush;
			cout << "--------------------------------------------\n" << flush;
			cout << "   line_seq_names                       = " << line_seq_names[number_of_line]
			     << "\n   source_name                    = " << line_source_names[number_of_line]<<endl;
			for(i=0;i<number_of_annotation_labels; i++){
			    for(j=0;j<MP->get_Annotation_Label_size(i);j++){
				if(line_annotation_labels[number_of_line][i][j]==true){
				    cout<<MP->get_Annotation_Label_setname(i)<<" ="
					<<MP->get_Annotation_Label_name(i,j)
					<<" socre : "<<line_scores[number_of_line][i][j]<<endl;
				}
			    }
			}				

			cout << "\n   start_position                 = " << line_start_positions[number_of_line]
			     << "\n   end_position                   = " << line_end_positions[number_of_line]
			     << flush;
			cout << "\n--------------------------------------------\n" << flush;
			check++;
			break;
		    }

		    if(valid)
		    {
			// get series of annotation_labels

			while(1){

			    if(check){
				break;
			    }
			    
			    if(tmp_pos>=NoOfItems){
				break;
			    }			
			    // end of annotation labels
			    if(!strcmp(Items[tmp_pos],"|")){
				tmp_pos++;
				break;
			    }
			    
			    check+=  splitstring(Items[tmp_pos], NoOfItemsForLabel, &ItemsForLabel,':');
			    if(check){
				cout<<"ERROR: read_special_file_for_sequence: splitstring for Items["<<tmp_pos
				    <<"] : "<<Items[tmp_pos]<<"."<<endl;
				break;
			    }

			    if(NoOfItemsForLabel<2)
			    {
				cout<<"ERROR: read_special_file_for_sequence: feature label : "<<Items[tmp_pos]<<" not in a right format."<<endl;
				check++;
				break;
			    }
			    
			    tmp_label_index = convert_typename_to_int(MP->get_Annotation_Label(),ItemsForLabel[0],number_of_annotation_labels);
			    
			    if(tmp_label_index==-1){
				tmp_pos++;
				continue;
			    }else{

				tmp_pos++;
				check+= splitstring(Items[tmp_pos],NoOfItemsForScore, &ItemsForScore,':');
				if(check)
				{
				    cout<<"ERROR: read_special_file_for_sequence: splitstring for Items["<<tmp_pos
					<<"] : "<<Items[tmp_pos]<<"."<<endl;
				    break;
				}

				if(NoOfItemsForScore!=NoOfItemsForLabel-1)
				{
				    cout<<"ERROR: read_special_file_for_sequence: score label : "<<Items[tmp_pos]
					<<" for feature label : "<<ItemsForLabel[0]<<"not in valid format"<<endl;
				    check++;
				    break;
				}
				
				feature_label_size =MP->get_Annotation_Label_size(tmp_label_index);
				//found_keyword = 0;
						
				for(j=1; j<NoOfItemsForLabel; j++)
				{
				    for(i=0; i<feature_label_size; i++)
				    {
					if(cmp_nocase(ItemsForLabel[j],MP->get_Annotation_Label_name(tmp_label_index,i))==0)
					{
					    if(line_annotation_labels[number_of_line][tmp_label_index][i]==true){
						cout<<"Error:: read_special_file_for_sequence : "
						    <<"same label has been defined at the same by position by the same source!"<<endl;
						check++;
						break;
					    }else{
						line_annotation_labels[number_of_line][tmp_label_index][i]=true;
						//found_keyword ++ ;
					    }
					
					    //tmp_pos++;
					    if(!strcmp(ItemsForScore[j-1],"."))
					    {
						line_scores[number_of_line][tmp_label_index][i] = 0;
						//tmp_pos++;
					    }else if(!MP->get_score_of_Annotation_Label(tmp_label_index)){
						cout<<"Error:: read_special_file_for_sequence : score("<<Items[tmp_pos]<<") is set for non score label("<<MP->get_Annotation_Label_setname(tmp_label_index)<<")"<<endl;
						check++;
					    }else{
						// check if the score is valid
						//double tmp_score = atof(Items[tmp_pos]);
						double tmp_score = atof(ItemsForScore[j-1]);
						
						if((tmp_score<0)||(tmp_score>1))
						{
						    cout<<"Error:: read_special_file_for_sequence : score("<<ItemsForScore[j-1]
							<<") out of range [0..1].\n";
						    cout<<"Interrupted line : \n"
							<<line<<endl;
						    check++;
						    tmp_score = 0;
						}else if((tmp_score==0)&&(strcmp(ItemsForScore[j-1],"0")))
						{
						    cout<<"Error:: read_special_file_for_sequence : score("<<ItemsForScore[j-1]
							<<") not valid, should be a numerical number in between 0 and 1.\n";
						    cout<<"Interrupted line : \n"
							<<line<<endl;
						    check++;
						    tmp_score = 0;
						}
						line_scores[number_of_line][tmp_label_index][i]=tmp_score;
						
						//tmp_pos++;
					    }				    
					    break;
					}				    
				    }
				}				
				/*
				if(found_keyword == 0){
				    tmp_pos++;
				}			
				*/
				tmp_pos++;
			    }
			}
   
			if(!check)
			{
			    number_of_line++;
			
			    if ((number_of_line == max_number_of_lines)&&(!feof(special_file)))
			    {
				cout << "   ERROR: read_special_file_for_sequence: max_number_of_lines ("
				     << max_number_of_lines << ") == number_of_line ("
				     << number_of_line << "). Increase the max value.\n" << flush;
				check++;
				break;
			    }
			}
		    }
		}		
	    }
	}	
	// close file
	fclose(special_file);

	// release memory
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

	if(ItemsForScore)
	{
	    for(i=0; i<Max_number_of_items; i++)
	    {
		if(ItemsForScore[i]) delete [] ItemsForScore[i];
		ItemsForScore[i] = NULL;
	    }
	    delete[] ItemsForScore;
	}
	ItemsForScore = NULL;
	NoOfItemsForScore = 0;

	if(!check)
	{
	    int k = 0;
	    int l =0;
	    int label_size = 0;
	    // initialize output variables	    
	    (*number_of_allocated_lines) = number_of_line;
	    (*number_of_sources) = 0;
	    (*number_of_lines) = NULL;  // counts numbers of gtf_lines stored in the output arrays from each source   

	    check+=sort_lines_by_sources(//input 
		number_of_line,
		MP,
		&line_seq_names,
		&line_source_names,
		&line_annotation_labels,
		&line_scores,
		&line_start_positions,
		&line_end_positions,
		//output
		number_of_sources,
		number_of_lines);
	            
	    if(check){
		cout<<"ERROR : read_special_file_for_sequence : error in sort_lines_by_sources "<<endl;
	    }

	    if(!check){
		
		check+=initialize_labels((*number_of_sources),
					 (*number_of_lines),
					 MP,
					 seq_names,
					 source_names,
					 annotation_labels,
					 scores,
					 start_positions,
					 end_positions);
		if(check){
		    cout<<"ERROR : read_special_label_for_sequence, initialize_labels for source "<<endl;
		}

		//set output variables
		int source_number = 0;
		int count = 0;
	       
		for(i=0; i<number_of_line; i++)
		{
		    if(i!=0){
			if(strcmp(line_source_names[i-1],line_source_names[i])!=0){
			    source_number++;
			    if(source_number>=(*number_of_sources)){
				cout<<"Error:: read_special_file_for_sequence : source_number("<<source_number+1
				    <<") exceed number of source : "<<(*number_of_sources)<<endl;
				check++;
				break;
			    }
			    strcpy((*source_names)[source_number],line_source_names[i]);
			    count = 0;	
			}
		    }else{
			strcpy((*source_names)[source_number],line_source_names[i]);
		    }
		    
		    strcpy((*seq_names)[source_number][count],line_seq_names[i]);
		    
		    for(j=0;j<number_of_annotation_labels;j++){
			label_size = MP->get_Annotation_Label_size(j);
			for(k=0;k<label_size;k++){
			    (*annotation_labels)[source_number][count][j][k] = line_annotation_labels[i][j][k];
			    (*scores)[source_number][count][j][k] = line_scores[i][j][k];
			}
		    }
		    
		    (*start_positions)[source_number][count] = line_start_positions[i];
		    (*end_positions)[source_number][count]   = line_end_positions[i];
	
		    count++;
		    
		    if(count>(*number_of_lines)[source_number]){
			cout<<"Error:: read_special_file_for_sequence: count : "<<count+1 
			    <<" exceed the limit of line("<<source_number<<") numebr_of_lines["<<source_number<<"] : "
			    <<(*number_of_lines)[source_number]<<".\n"<<flush;
			check++;
			break;
		    }
		}

		// if check != 0 delete the output information
		
		if (check)
		{		
		    // release memory
		    
		    check+=release_memory_for_labels(number_of_sources,
						     number_of_lines,
						     MP,
						     seq_names,
						     source_names,
						     annotation_labels,
						     scores,
						     start_positions,
						     end_positions);		    
		  
		    if(check){
			cout<<"ERROR : read_special_label_for_sequence, release_memory_for__labels for source "<<endl;
		    }	

		    (*number_of_allocated_lines) = 0;
		}
	    }
	}
	// delete temporary varibales

	if(line) delete[] line;
	line = NULL;

	if(seq_name_constraint) delete[] seq_name_constraint;
	seq_name_constraint = NULL;

	check += release_memory_for_labels(&max_number_of_lines,
					   MP,
					   &line_seq_names,
					   &line_annotation_labels,
					   &line_scores,
					   &line_start_positions,
					   &line_end_positions);
	
	if(check){
	    cout<<"ERROR : read_special_label_for_sequence, release_memory_for_labels for line "<<endl;
	}	

	// delete source_names
	if(line_source_names){
	    for(i=0; i<max_number_of_lines; i++){
		if(line_source_names[i]) delete[] line_source_names[i];
		line_source_names[i] = NULL;
	    }
	    delete[] line_source_names;
	}
	line_source_names = NULL;
    }else{
	(*number_of_allocated_lines) = 0;
    }
    return(check);
}

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
     )
{
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    if(number_of_sources<0){
	cout<<"ERROR : initialize_labels : number_of_sources("<<number_of_sources<<") < 0 .\n"<<flush;
	check++;
    }else{
	for(i=0; i<number_of_sources; i++)
	{
	    if(number_of_lines[i]<=0){
		cout<<"ERROR: initialize_labels : number_of_lines["<<i<<"]("<<number_of_lines[i]
		    <<") <= 0.\n"<<flush;
		check++;
	    }
	}
    }

    if ((*seq_names) != NULL)
    {
	cout << "ERROR: initialize_labels: seq_names should be NULL.\n" << flush;
	check++;
    }
    if ((*source_names) != NULL)
    {
	cout << "ERROR:: initialize_labels source_names should be NULL.\n" << flush;
	check++;
    }
    if ((*annotation_labels) != NULL)
    {
	cout << "ERROR: initialize_labels: annotation_labels should be NULL.\n" << flush;
	check++;
    }
    if ((*scores) != NULL)
    {
	cout << "ERROR: initialize_labels: scores should be NULL.\n" << flush;
	check++;
    }

    if ((*start_positions) != NULL)
    {
	cout << "ERROR: initialize_labels: start_positions should be NULL.\n" << flush;
	check++;
    }
    if ((*end_positions) != NULL)
    {
	cout << "ERROR: initialize_labels: end_positions should be NULL.\n" << flush;
	check++;
    }

    if(number_of_sources==0)
    {
	return check;
    }

    if(!check){

	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
	//int number_of_other_labels = MP->get_Total_Number_of_Other_Labels();
	int label_size = -1;

	(*seq_names)       = new char**[number_of_sources];
	for(i=0;i<number_of_sources;i++){
	    (*seq_names)[i] = new char*[number_of_lines[i]]; 
	    for(j=0;j<number_of_lines[i];j++){
		(*seq_names)[i][j] = new char[Max_word_length];
		strcpy((*seq_names)[i][j]," ");
	    }
	}
	(*source_names)    = new char*[number_of_sources]; 
	for(i=0;i<number_of_sources;i++){
	    (*source_names)[i] = new char[Max_word_length];
	    strcpy((*source_names)[i]," ");
	}
	
	(*annotation_labels)     = new bool***[number_of_sources]; 
	(*scores)                = new double***[number_of_sources];
	for(i=0;i<number_of_sources ;i++){
	    (*annotation_labels)[i] = new bool**[number_of_lines[i]];
	    (*scores)[i]            = new double**[number_of_lines[i]];
	    for(j=0;j<number_of_lines[i];j++){
		(*annotation_labels)[i][j] = new bool*[number_of_annotation_labels];
		(*scores)[i][j]             = new double*[number_of_annotation_labels];
		for(k=0;k<number_of_annotation_labels;k++){
		    label_size = MP->get_Annotation_Label_size(k);
		    (*annotation_labels)[i][j][k] = new bool[label_size];
		    (*scores)[i][j][k]             = new double[label_size];
		    for(l=0;l<label_size;l++){
			(*annotation_labels)[i][j][k][l] = false;
			(*scores)[i][j][k][l]             = 0;
		    }
		}
	    }
	}
	
	(*start_positions) = new int*[number_of_sources];
	(*end_positions)   = new int*[number_of_sources];
	for(i=0; i<number_of_sources; i++){
	    (*start_positions)[i] = new int[number_of_lines[i]]; 
	    (*end_positions)[i]   = new int[number_of_lines[i]];
	    for(j=0; j<number_of_lines[i]; j++)
	    {
		(*start_positions)[i][j] = 0;
		(*end_positions)[i][j] = 0;
	    }
	}
    }
    return check;
}

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
    )
{
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;

    if(number_of_lines<=0){
	cout<<"ERROR : initialize_labels : number_of_lines("<<number_of_lines<<") <=0 .\n"<<flush;
	check++;
    }

    if ((*seq_names) != NULL)
    {
	cout << "ERROR: initialize_labels: seq_names should be NULL.\n" << flush;
	check++;
    }
   
    if ((*annotation_labels) != NULL)
    {
	cout << "ERROR: initialize_labels: annotation_labels should be NULL.\n" << flush;
	check++;
    }
    if ((*scores) != NULL)
    {
	cout << "ERROR: initialize_labels: scores should be NULL.\n" << flush;
	check++;
    }

    if ((*start_positions) != NULL)
    {
	cout << "ERROR: initialize_labels: start_positions should be NULL.\n" << flush;
	check++;
    }
    if ((*end_positions) != NULL)
    {
	cout << "ERROR: initialize_labels: end_positions should be NULL.\n" << flush;
	check++;
    }

    if(!check)
    {
	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
	int label_size = -1;

	(*seq_names)       = new char*[number_of_lines];
	for(i=0;i<number_of_lines;i++){
	    (*seq_names)[i] = new char[Max_word_length];
	    strcpy((*seq_names)[i]," ");   
	}
	
	(*annotation_labels)     = new bool**[number_of_lines]; 
	(*scores)                = new double**[number_of_lines];
	for(i=0;i<number_of_lines ;i++)
	{	   
	    (*annotation_labels)[i] = new bool*[number_of_annotation_labels];
	    (*scores)[i]             = new double*[number_of_annotation_labels];
	    for(j=0;j<number_of_annotation_labels;j++){
		label_size = MP->get_Annotation_Label_size(j);
		(*annotation_labels)[i][j]  = new bool[label_size];
		(*scores)[i][j]              = new double[label_size];
		for(k=0;k<label_size;k++){
		    (*annotation_labels)[i][j][k] = false;
		    (*scores)[i][j][k]             = 0;
		}					    
	 
	    }
	}
	(*start_positions) = new int[number_of_lines];
	(*end_positions)   = new int[number_of_lines];
	for(i=0; i<number_of_lines; i++){
	    (*start_positions)[i] = 0; 
	    (*end_positions)[i]   = 0;
	}
    }
    return check;    
}

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
    )
{
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    
    if((*number_of_sources)<0)
    {
	cout<<"ERROR: release_memory : number_of_sources("<<number_of_sources<<") < 0.\n"<<flush;
	check++;
    }else{
	for(i=0;i<(*number_of_sources); i++){
	    if((*number_of_lines)[i]<0){
		cout<<"ERROR: release_memory : number_of_lines["<<i<<"]("<<(*number_of_lines)[i]
		    <<") < 0.\n"<<flush;
		check++;
	    }
	}
    }
    if(!check)
    {
	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
	
	// delete seq_names
	if(*seq_names){
	    for(i=0;i<(*number_of_sources);i++){
		if((*seq_names)[i]){
		    for(j=0; j<(*number_of_lines)[i];j++){
			if((*seq_names)[i][j]) delete[] (*seq_names)[i][j];
			(*seq_names)[i][j] = NULL;
		    }
		    delete[] (*seq_names)[i];
		}
		(*seq_names)[i] = NULL;
	    }
	    delete[] (*seq_names);
	}
	(*seq_names) = NULL;

	// delete annotation_labels
	
	if(*annotation_labels){
	    for(i=0; i<(*number_of_sources); i++){
		if((*annotation_labels)[i]){
		    for(j=0;j<(*number_of_lines)[i];j++){
			if((*annotation_labels)[i][j]){
			    for(k=0; k<number_of_annotation_labels; k++){
				if((*annotation_labels)[i][j][k]) delete[] (*annotation_labels)[i][j][k];
				(*annotation_labels)[i][j][k] = NULL;
			    }
			    delete[] (*annotation_labels)[i][j];
			}
			(*annotation_labels)[i][j] = NULL;
		    }
		    delete[] (*annotation_labels)[i];
		}
		(*annotation_labels)[i] = NULL;
	    }
	    delete[] (*annotation_labels);
	}
	(*annotation_labels) = NULL;
			
	// delete scores
	if(*scores){
	    for(i=0; i<(*number_of_sources); i++){
		if((*scores)[i]){
		    for(j=0;j<(*number_of_lines)[i];j++){
			if((*scores)[i][j]){
			    for(k=0; k<number_of_annotation_labels; k++){
				if((*scores)[i][j][k]) delete[] (*scores)[i][j][k];
				(*scores)[i][j][k] = NULL;
			    }
			    delete[] (*scores)[i][j];
			}
			(*scores)[i][j] = NULL;
		    }
		    delete[] (*scores)[i];
		}
		(*scores)[i] = NULL;
	    }
	    delete[] (*scores);
	}
	(*scores) = NULL;
		    
	// delete start_positions
	if (*start_positions){
	    for(i=0; i<(*number_of_sources); i++){
		if((*start_positions)[i]) delete[] (*start_positions)[i];
		(*start_positions)[i] = NULL;
	    }
	    delete[] (*start_positions);
	}
	(*start_positions) = NULL;		    
	// delete end_positions
	if (*end_positions){
	    for(i=0; i<(*number_of_sources); i++){
		if((*end_positions)[i]) delete[] (*end_positions)[i];
		(*end_positions)[i] = NULL;
	    }
	    delete[] (*end_positions);
	}
	(*end_positions) = NULL;

	// delete source_names
	if(*source_names){
	    for(i=0; i<(*number_of_sources); i++){
		if((*source_names)[i]) delete [] (*source_names)[i];
		(*source_names)[i] = NULL;
	    }
	    delete[] (*source_names);
	}
	(*source_names) = NULL;
	
	// delete number_of_lines
	if (*number_of_lines) delete[] (*number_of_lines);
	(*number_of_lines) = NULL;
	(*number_of_sources) = 0;
    }
    return 0;
}

int release_memory_for_labels(
    // input
    int*           const number_of_lines,
    model_parameters* const MP,
    // output
    char***      const seq_names,
    bool****     const annotation_labels,
    double****   const scores,
    int**        const start_positions,
    int**        const end_positions
    )
{
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    
    if((*number_of_lines)<0)
    {
	cout<<"ERROR: release_memory : number_of_lines("<<number_of_lines<<") < 0.\n"<<flush;
	check++;
    }
 
    if(!check)
    {
	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();	

	// delete seq_names
	if(*seq_names){
	    for(i=0;i<(*number_of_lines);i++){
		if((*seq_names)[i]) delete[] (*seq_names)[i];
		(*seq_names)[i] = NULL;
	    }
	    delete[] (*seq_names);
	}
	(*seq_names) = NULL;

	// delete annotation_labels	
	if(*annotation_labels){
	    for(i=0; i<(*number_of_lines); i++){
		if((*annotation_labels)[i]){
		    for(j=0; j<number_of_annotation_labels; j++){
			if((*annotation_labels)[i][j]) delete[] (*annotation_labels)[i][j];
			(*annotation_labels)[i][j] = NULL;
		    }
		    delete[] (*annotation_labels)[i];
		}
		(*annotation_labels)[i] = NULL;
	    }
	    delete[] (*annotation_labels);
	}
	(*annotation_labels) = NULL;
 
	// delete scores
	if(*scores){
	    for(i=0; i<(*number_of_lines); i++){
		if((*scores)[i]){
		    for(j=0; j<number_of_annotation_labels; j++){
			if((*scores)[i][j]) delete[] (*scores)[i][j];
			(*scores)[i][j] = NULL;
		    }              
		    delete[] (*scores)[i];
		}
		(*scores)[i] = NULL;
	    }
	    delete[] (*scores);
	}
	(*scores) = NULL;
	            
	// delete start_positions
	if (*start_positions)    delete[] (*start_positions);
	(*start_positions) = NULL;		    
	// delete end_positions
	if (*end_positions)  delete[] (*end_positions);
	(*end_positions) = NULL;

	(*number_of_lines) = 0;
    }
    return 0;
}

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
    int**                number_of_lines)
{
    int check = 0;
    bool finished = true;
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int label_size = 0;
    int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
    
    // perform bubble sort
    for(i=0; i<number_of_line;i++){
	finished = true;
	for(j=0;j<number_of_line-i-1;j++){
	    if(strcmp((*line_source_names)[j],(*line_source_names)[j+1])>0){
		check+=swap_string(&((*line_source_names)[j]),&((*line_source_names)[j+1]));
		check+=swap_string(&((*line_seq_names)[j]),&((*line_seq_names)[j+1])); // swap seq_name
		for(k=0;k<number_of_annotation_labels;k++){
		    label_size = MP->get_Annotation_Label_size(k);
		    for(l=0;l<label_size;l++){
			swap(&((*line_annotation_labels)[j][k][l]),&((*line_annotation_labels)[j+1][k][l])); // swap word_annotation_label
			swap(&((*line_scores)[j][k][l]),&((*line_scores)[j+1][k][l])); // swap score
		    }
		}

		swap(&((*line_start_positions)[j]),&((*line_start_positions)[j+1]));  // swap start position
		swap(&((*line_end_positions)[j]),&((*line_end_positions)[j+1]));  // swap end position
		finished = false;
	    }
	}
	if(finished){
	    break;
	}
    }
    
    int n_of_source = 0;
    int count = 0;
    
    for(i=0; i<number_of_line; i++){
	if(i!=0){
	    if(strcmp((*line_source_names)[i-1],(*line_source_names)[i])!=0){
		n_of_source++;	       
	    }
	}else{
	    n_of_source = 1;
	}
    }
    (*number_of_sources) = n_of_source;

    if(!(*number_of_lines))
    {
	(*number_of_lines) = new int[(*number_of_sources)];
    }
    int source_count = 0;

    for(i=0; i<number_of_line; i++)
    {
	if(i!=0){
	    if(strcmp((*line_source_names)[i-1],(*line_source_names)[i])!=0){
		(*number_of_lines)[source_count] = count;
		source_count++;
		count = 1;
	    }else{
		count++;
	    }
	}else{
	    count=1;
	}
    }
    (*number_of_lines)[source_count] = count;
    return check;
}

int get_positions_of_feature_start_and_end(// input
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
    int*** flags_end)   // absolute end_positions of feature
{
    // note: - flags_start and _end have to be NULL when entering this function and the memory
    //         allocated for flags_start and _end within this function must be released in the
    //         program which calls this function 
    
    int check = 0;
    
    if (x == NULL)
    {
	cout << "ERROR: get_positions_of_feature_start_and_end : Sequence x is NULL.\n" << flush;
	check++;
    }
    else
    {
	if (x->get_end_position() < x->get_start_position())
	{
	    cout << "ERROR: get_positions_of_feature_start_and_end : start_position ("
		 << x->get_start_position() << ") of Sequence x is > end_position ("
		 << x->get_end_position() << ").\n" << flush;
	    check++;
	}
	if (x->get_ac() == NULL)
	{
	    cout << "ERROR: get_positions_of_feature_start_and_end : ac of Sequence x is NULL.\n" << flush;
	    check++;
	}
    }
    if (file_name == NULL)
    {
	cout << "ERROR: get_positions_of_feature_start_and_end : file_name is NULL.\n" << flush;
	check++;
    }
    
    if ((*flags_start) != NULL)
    {
	cout << "ERROR: get_positions_of_feature_start_and_end : flags_start has to be NULL.\n" << flush;
	check++;
    }
    if ((*flags_end) != NULL)
    {
	cout << "ERROR: get_positions_of_feature_start_and_end : flags_end has to be NULL.\n" << flush;
	check++;
    }
    
    if (check == 0)
    {
	FILE* special_file   = fopen(file_name, "rt");
	
	int i = 0;
	int j = 0;
	const int max_number_of_sources = 10;
	const int max_number_of_lines = 3000;
	
	// info about sequence x
	
	const int length         = x->length();
	const int start_position = x->get_start_position();
	const int end_position   = x->get_end_position();
	const int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
		
	// get gtf-lines which correspond to the sequence x
	
	int        number_of_sources = 0;
	int*       number_of_lines   = NULL;
	char***    seq_names         = NULL;
	char**     source_names      = NULL;
	int**      start_positions   = NULL;
	int**      end_positions     = NULL;
	bool****   annotation_labels = NULL;
	double**** scores            = NULL;
	//char****   other_labels      = NULL;

	check += read_special_file(// input
	    special_file,
	    file_name,
	    max_number_of_lines,
	    MP,
	    // output
	    &number_of_sources,
	    &number_of_lines,
	    &seq_names,
	    &source_names,
	    &annotation_labels,
	    &scores,
	    &start_positions,
	    &end_positions
	    );
	
	if (check != 0)
	{
	    cout << "ERROR: get_positions_of_feature_start_and_end : error occurred in function read_gtf_file.\n" << flush;
	}

	// order those lines and select only those that are within the range defined by the
	// sequence
	
	char*** seq_name     = NULL;
	int number_of_source =1;
	int* number_of_line   = new int[number_of_source];
	number_of_line[0] = 1;
	
	if (check == 0)
	{	    
	    seq_name             = new char**[number_of_source];
	    for(i=0; i<number_of_source; i++)
	    {
		seq_name[i]         = new char*[number_of_line[i]]; 
		for(j=0; j<number_of_line[i]; j++)
		{
		    seq_name[i][j] = new char[Max_word_length];
		}
	    }

	    int** start_position = new int*[number_of_source];
	    for(i=0; i<number_of_source; i++)
	    {
		start_position[i] = new int[number_of_line[i]];
	    }
	    int** end_position   = new int*[number_of_source];
	    for(i=0; i<number_of_source; i++)
	    {
		end_position[i] = new int[number_of_line[i]];
	    }
	    
	    strcpy(seq_names[0][0], x->get_ac()); 		  
	    start_positions[0][0] = x->get_start_position();
	    end_positions[0][0]   = x->get_end_position();

	    check += apply_restriction_set_to_lines(// input
		max_number_of_sources,
		max_number_of_lines,
		MP,
		// restriction set
		&number_of_source,
		&number_of_line,
		&seq_name,        
		&start_position,  
		&end_position,   
		// gtf-set
		&number_of_sources,
		&number_of_lines,
		&seq_names,        
		&source_names,     
		&annotation_labels,
		&scores,
		&start_positions,  
		&end_positions);       
	    
	    
	    if(number_of_line) delete [] number_of_line;
	    number_of_line = NULL;
			    
	    if (start_position){
		for(i=0; i<number_of_source; i++){
		    if(start_position[i]) delete [] start_position[i];
		    start_position[i] = NULL;
		}
		delete [] start_position;
	    }
	    start_position = NULL;
	    
	    if (end_position){
		for(i=0; i<number_of_source; i++){
		    if(end_position[i]) delete [] end_position[i];
		    end_position[i] = NULL;
		}
		delete [] end_position;
	    }
	    end_positions = NULL;
	    
	    if (check != 0)
	    {
		cout << "ERROR: get_positions_of_feature_start_and_end : error occurred in function \n"
		     << "apply_restriction_set_to_lines.\n" << flush;
	    }
	}	
	
	// search for the requested transition along the whole sequence
	
	if (check == 0)
	{
	    // allocate memory for output array and initialise array with 0s
	    
	    (*output_number_of_sources) = number_of_sources;
	    (*number_of_start_flags) = new int[(*output_number_of_sources)];
	    (*number_of_end_flags)   = new int[(*output_number_of_sources)];
	    for(i=0; i<(*output_number_of_sources); i++){
		(*number_of_start_flags)[i] = 0;
		(*number_of_end_flags)[i] = 0;
	    }
	    (*flags_start)           = new int*[(*output_number_of_sources)];
	    (*flags_end)             = new int*[(*output_number_of_sources)];
	    for(i=0; i<(*output_number_of_sources); i++){
		(*flags_start)[i] = new int[length];
		(*flags_end)[i]   = new int[length];
		for(j=0; j<length; j++){
		    (*flags_start)[i][j] = 0;
		    (*flags_end)[i][j]   = 0;
		}
	    }
	    
	    // determine flags

	    // go through all gtf-lines and select the relevant line with the feature_of_interest   
	    for(i=0; i<number_of_sources; i++){		    
		for (j=0; j<number_of_lines[i]; j++)
		{
		    if (annotation_labels[i][j][feature_label_index][feature_of_interest])
		    {			
			(*flags_start)[i][(*number_of_start_flags)[i]] = start_positions[i][j];
			(*flags_end)[i][(*number_of_end_flags)[i]] = end_positions[i][j];
			(*number_of_start_flags)[i]++;
			(*number_of_end_flags)[i]++;
		    }
		}
	    }	    
	} // if check == 0
	
	// delete memory	
	check+= release_memory_for_labels(&number_of_sources,
					  &number_of_lines,
					  MP,
					  &seq_names,
					  &source_names,
					  &annotation_labels,
					  &scores,
					  &start_positions,
					  &end_positions);

	if(seq_name){
	    for(i=0; i<number_of_source; i++){
		if(seq_name[i]){
		    for(j=0;j<number_of_line[i];j++){
			if(seq_name[i][j]) delete[] seq_name[i][j];
			seq_name[i][j] = NULL;
		    }
		    delete [] seq_name[i];
		}
		seq_name[i] = NULL;
	    }
	    delete [] seq_name;
	}
	seq_name = NULL;
					  	
	fclose(special_file);
    } // if initial checks were passed and function could be entered
    
    // delete output, if checks were not passed
    
    if (check != 0)
    {
	if (*number_of_start_flags) delete [] (*number_of_start_flags);
	(*number_of_start_flags) = NULL;
	if(*number_of_end_flags) delete [] (*number_of_end_flags);
	(*number_of_end_flags) = NULL;
	
	if ((*flags_start)){
	    for(int i=0; i<(*output_number_of_sources); i++){
		if((*flags_start)[i]) delete [] (*flags_start)[i];
		(*flags_start)[i] = NULL;
	    }
	    delete [] (*flags_start);
	}
	(*flags_start) = NULL;
	
	if ((*flags_end)){
	    for(int i=0; i<(*output_number_of_sources); i++){
		if((*flags_end)[i]) delete [] (*flags_end)[i];
		(*flags_end)[i] = NULL;
	    }
	    delete [] (*flags_end);
	}
	(*flags_end) = NULL;

	(*output_number_of_sources) = 0;
    }
    
    return(check);
}

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
    // sources
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
    double****   const score)
{
    // note: this function assumes that gtf-lines are ordered 
    //       (ordered mean ordered according to function order_gtf_lines)
    
    int check = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    // input interval

    if ((*seq_name) == NULL)
    {
	cout << "ERROR: get_feature_label : seq_name is NULL.\n" << flush;
	check++;
    }
    
    if ((*end_position) < (*start_position))
    {
	cout << "ERROR: get_feature_label : end_position (" << (*end_position) 
	     << ") < start_position (" << (*start_position) << ")\n" << flush;
	check++;
    }
    
    
    // input prediction set

    if ((*number_of_sources) < 1)
    {
	cout << "ERROR: get_feature_label : number_of_sources (" << (*number_of_sources) << ") < 1.\n" << flush;
	check++;
    }

    for(i=0; i<(*number_of_sources); i++)
    {
	if ((*start_line)[i] < 0)
	{
	    cout << "ERROR: get_feature_label : start_line["<<i<<"] (" << (*start_line)[i] 
		 << ") < 0.\n" << flush;
	    check++;
	}
	
	if ((*end_line)[i] < (*start_line)[i])
	{
	    cout << "ERROR: get_feature_label : end_line["<<i<<"] (" << (*end_line)[i] 
		 << ") < gtf_start_line["<<i<<"] (" << (*start_line)[i] << ").\n" << flush;
	    check++;
	}
    }
    
    // gtf-lines
    
    if ((*seq_names) == NULL)
    {
	cout << "ERROR: get_feature_label : seq_names is NULL.\n" << flush;
	check++;
    }
    if ((*source_names) == NULL)
    {
	cout << "ERROR: get_feature_label : source_names is NULL.\n" << flush;
	check++;
    }

    if ((*annotation_labels) == NULL)
    {
	cout << "ERROR: get_feature_label : annotation_labels is NULL.\n" << flush;
	check++;
    }
    if ((*scores) == NULL)
    {
	cout << "ERROR: get_feature_label : scores is NULL.\n" << flush;
	check++;
    }

    if ((*start_positions) == NULL)
    {
	cout << "ERROR: get_feature_label : start_positions is NULL.\n" << flush;
	check++;
    }
    if ((*end_positions) == NULL)
    {
	cout << "ERROR: get_feature_label : end_positions is NULL.\n" << flush;
	check++;
    }
   
    // check some consistencies

    if ((position < (*start_position)) || (position > (*end_position)))
    {
	cout << "ERROR: get_feature_label : either position (" << position << ") < start_position ("
	     << (*start_position) << ") or > end_position (" << (*end_position) << ")\n" << flush;
	check++;
    }

    for(i=0; i<(*number_of_sources); i++){
    
	if ((*number_of_lines)[i] > (*end_line)[i])
	{
	    if (min((*start_positions)[i][(*start_line)[i]], (*start_positions)[i][(*end_line)[i]]) 
		< (*start_position))
	    {
		cout << "ERROR: get_feature_label : min(start_positions["<<i
		     <<"][start_line["<<i<<"]] = start_positions["<<i<<"]["
		     << (*start_line)[i] << "] (" << (*start_positions)[i][(*start_line)[i]]
		     << "), start_positions["<<i<<"][gtf_end_line["<<i<<"]] = start_positions["<<i<<"]["
		     << (*end_line)[i] << "] (" << (*start_positions)[i][(*end_line)[i]]
		     << ")) < start_position (" << (*start_position) << ").\n" << flush;
		check++;
	    }
	    if (max((*end_positions)[i][(*start_line)[i]], (*end_positions)[i][(*end_line)[i]]) 
		> (*end_position))
	    {
		cout << "ERROR: get_feature_label : max(end_positions["<<i<<"][gtf_start_line["
		     <<i<<"] = end_positions["<<i<<"]["
		     << (*start_line)[i] << "] (" << (*end_positions)[i][(*start_line)[i]]
		     << "), end_positions["<<i<<"][end_line["<<i<<"]] = end_positions["
		     << (*end_line)[i] << "] (" << (*end_positions)[i][(*end_line)[i]]
		     << ")) > end_position (" << (*end_position) << ").\n" << flush;
		check++;
	    }
	    
	    if ((strcmp((*seq_name), (*seq_names)[i][(*start_line)[i]]) != 0) ||
		(strcmp((*seq_name), (*seq_names)[i][(*end_line)[i]]) != 0))
	    {
		cout << "ERROR: get_feature_label : either seq_name (" << (*seq_name) 
		     << ") != seq_names["<<i<<"][start_line["<<i<<"]] = seq_names[" 
		     <<i<<"]["<< (*start_line)[i]
		     << "] (" << (*seq_names)[i][(*start_line)[i]] 
		     << ") or seq_name (" << (*seq_name) 
		     << ") != seq_names["<<i<<"][end_line["<<i<<"]] = seq_names["
		     <<i<<"]["<< (*end_line)[i]
		     << "] (" << (*seq_names)[i][(*end_line)[i]] << ").\n" << flush;
		check++;
	    }
	}
	else
	{
	    cout << "ERROR: get_feature_label : number_of_lines["<<i<<"] (" << (*number_of_lines)[i]
		 << ") <= end_line["<<i<<"] (" << (*end_line)[i] << ")\n" << flush;
	    check++;
	}
    }
    
    if (check == 0)
    {
	// max_number_of_labels = maximum number of labels one piece of sequence can have
	
	// buh const int max_number_of_labels = 3;
	
	const int max_number_of_labels = 10;
	
	int i = 0;
	int j = 0;

	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();

	if(!(*feature)){
	    (*feature) = new bool**[(*number_of_sources)];
	    for(i=0; i<(*number_of_sources); i++){
		(*feature)[i] = new bool*[number_of_annotation_labels];
		for(j=0;j<number_of_annotation_labels; j++){
		    (*feature)[i][j] = new bool[MP->get_Annotation_Label_size(j)];
		    for(k=0; k<MP->get_Annotation_Label_size(j); k++){
			(*feature)[i][j][k] = false;
		    }
		}
	    }
	}

	if(!(*score)){
	    (*score) = new double**[(*number_of_sources)];
	    for(i=0; i<(*number_of_sources); i++){
		(*score)[i] = new double*[number_of_annotation_labels];
		for(j=0; j<number_of_annotation_labels; j++){
		    (*score)[i][j] = new double[MP->get_Annotation_Label_size(j)];
		    for(k=0; k<MP->get_Annotation_Label_size(j); k++){
			(*score)[i][j][k] = 0;
		    }
		}
	    }
	}

	bool* set_score = new bool[number_of_annotation_labels];
	
	for (i=0; i<(*number_of_sources); i++)
	{
	    if(check){
		break;
	    }

	    for(j=0; j<number_of_annotation_labels; j++){
		set_score[j] = false;
	    }
	 
	    for (int line = (*start_line)[i]; line<(*end_line)[i]+1; line++)
	    {
		if(check){
		    break;
		}
		if ((position >= (*start_positions)[i][line]) &&
		    (position <= (*end_positions)[i][line]))
		{		   
		    // assign values
		    for(j=0; j<number_of_annotation_labels; j++){
			for(k=0; k<MP->get_Annotation_Label_size(j); k++){
			    if((*annotation_labels)[i][line][j][k])
			    {
				if(MP->get_score_of_Annotation_Label(j))
				{
				    if(((*feature)[i][j][k]==true)&&
				       ((*score)[i][j][k]!=(*scores)[i][line][j][k]))
				    {
					cout<<"ERROR: get_feature_label : same feature_label was assigned by the same source with same score." <<endl;
					check++;
					break;
				    }
				    (*feature)[i][j][k] = (*annotation_labels)[i][line][j][k];
				    (*score)[i][j][k] = (*scores)[i][line][j][k];
				    set_score[j] = true;
				}else{
				    if(set_score[j]){
					if(!(*feature)[i][j][k])
					{
					    cout<<"ERROR: get_feature_label : different label was assigned with the same source while score is not allowed for this feature label set."<<endl;
					}
				    }else{
					(*feature)[i][j][k] = (*annotation_labels)[i][line][j][k];
					set_score[j] = true;
				    }
				}
			    }		    
			}			 
		    }
		}
		else
		{
		    cout << "ERROR: get_feature_label : position (" << position 
			 << ") has already more than max_number_of_labels (" << max_number_of_labels
			 << ") labels associated with it.\n" << flush;
		    check++;
		    break;
		}
	    }
	}
	
	if (check != 0)
	{
	    for(i=0; i<(*number_of_sources); i++){
		for(j=0; j<number_of_annotation_labels; j++){
		    for(k=0; k<MP->get_Annotation_Label_size(j); k++)
		    {
			(*feature)[i][j][k] = false;
			(*scores)[i][j][k] = 0;
		    }
		}
	    }
	}      
	if(set_score) delete[] set_score;
	set_score = NULL;
    }   
    return(check);
}

int read_special_file(// input
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
    )
{
    // note: - variables that will be used to store output information have to be NULL when calling this function
    //       - upon successfull call of this function the memory of the output variables has to be released within the
    //         calling program
    //       - gtf-comment lines start with '//' or '#'
    
    int check=0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    
    int tmp_number_of_lines             = 0;
    
    if (file_ptr == NULL)
    {
	cout << "ERROR: read_special_file: file_ptr is NULL.\n" << flush;
	check++;
    }
    if (file_name == NULL)
    {
	cout << "ERROR: read_special_file: file_name is NULL.\n" << flush;
	check++;
    }
    if (max_number_of_lines < 1)
    {
	cout << "ERROR: read_special_file: max_number_of_lines (" << max_number_of_lines << ") < 1.\n" << flush;
	check++;
    }
    

    // check that the number of lines in gtf_file does not exceed the max_number_of_gtf_lines

    if ((file_ptr != 0) && (file_name != NULL))
    {
	char* command_line = new char[Max_word_length];
	
	strcpy(command_line, "less ");
	strcat(command_line, file_name);
	strcat(command_line, " | wc -l");
	FILE* line_count_command=popen(command_line, "r");
	fscanf(line_count_command, "%i", &tmp_number_of_lines);
	pclose(line_count_command);
	
	if (command_line) delete [] command_line;
	command_line = NULL;

	if (tmp_number_of_lines > max_number_of_lines)
	{
	    cout << "ERROR: read_special_file: number_of_lines in file (" << tmp_number_of_lines
		 << ") > max_number_of_lines = " << max_number_of_lines << ". Increase this max value.\n" << flush;
	    check++;
	}
    }
    if ((*seq_names) != NULL)
    {
	cout << "ERROR: read_special_file: seq_names should be NULL.\n" << flush;
	check++;
    }
    if ((*source_names) != NULL)
    {
	cout << "ERROR: read_special_file: source_names should be NULL.\n" << flush;
	check++;
    }
    
    if((*annotation_labels)!=NULL){
	cout << "ERRORL read_special_file: annotation_labels should be NULL.\n" <<flush;
	check++;
    } 
    if((*scores)!=NULL){
	cout << "ERRORL read_special_file: scores should be NULL.\n" <<flush;
	check++;
    }

    if ((*start_positions) != NULL)
    {
	cout << "ERROR: read_special_file: start_positions should be NULL.\n" << flush;
	check++;
    }
    if ((*end_positions) != NULL)
    {
	cout << "ERROR: read_special_file: end_positions should be NULL.\n" << flush;
	check++;
    }
   
    if (check == 0)
    {
	
	int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();

	int NoOfItems = 0;	
	char** Items = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    Items[i] = new char[Max_word_length];
	    strcpy(Items[i]," ");
	}
	
	int NoOfItemsForLabel = 0;
	char** ItemsForLabel = new char*[Max_number_of_items];
	for(i=0; i<Max_number_of_items; i++)
	{
	    ItemsForLabel[i] = new char[Max_word_length];
	    strcpy(ItemsForLabel[i]," ");
	}

	char** line_seq_names;
	bool*** line_annotation_labels;	
	double*** line_scores; 
	int*   line_start_positions;
	int*   line_end_positions;
	check+= initialize_labels(max_number_of_lines,
				  MP,
				  &line_seq_names,
				  &line_annotation_labels,
				  &line_scores,
				  &line_start_positions,
				  &line_end_positions);
	if(check){
	    cout<<"ERROR : read_special_labels, initialize_labels for line "<<endl;
	}	
	
	char** line_source_names  = new char*[max_number_of_lines];
	for(i=0; i<max_number_of_lines; i++){
	    line_source_names[i] = new char[Max_word_length];
	    strcpy(line_source_names[i]," ");
	}
	
	char* line                     = new char[Max_line_length];

	int number_of_line = 0;

	bool valid = false;
		
	while ((! feof(file_ptr)) && (check == 0))
	{
	    // initialise variables
	           
	    strcpy(line, " ");        
	    strcpy(line_seq_names[number_of_line], " ");        
	    strcpy(line_source_names[number_of_line], " ");        
	    line_start_positions[number_of_line] = 0;
	    line_end_positions[number_of_line] = 0;

	    for(i=0;i<number_of_annotation_labels;i++)
	    {
		for(j=0;j<MP->get_Annotation_Label_size(i); j++)
		{
		    line_scores[number_of_line][i][j] = 0;
		}
	    }
	    
	    for(i=0; i<number_of_annotation_labels; i++)
	    {
	        for(j=0; j<MP->get_Annotation_Label_size(i); j++)
		{
		    line_annotation_labels[number_of_line][i][j] = false;
		}
	    } 

	    // read one line if either max_line_length-1 characters were read or end of line was encountered
	    
	    fgets(line, Max_line_length-1, file_ptr);
	    
	    // decide of which type this line is 1.) comment, 2.) empty line or 3.) gtf_format line

	    check+=  splitstring(line, NoOfItems, &Items,' ');

	    if(check){
		cout<<"ERROR: read_special_file: splitstring : "<<line<<".\n";
		break;
	    }

	    valid = false;
	    
	    if (!((strcmp(Items[0], "//") == 0) || (strcmp(Items[0], "#") == 0)))
	    {
		
		int found_keyword  = 0;

		int tmp_pos = 0;
		
		tmp_pos = 4;

		int tmp_label_index = -1;
		int feature_label_size = 0;
		
		// keep reading the annotation label until it ends
		while(1){
		    
		    if((tmp_pos>=NoOfItems)||(found_keyword)){
			break;
		    }
		    
		    // end of annotation labels
		    if(!strcmp(Items[tmp_pos],"|")){
			break;
		    }
		    
		    check+=  splitstring(Items[tmp_pos], NoOfItemsForLabel, &ItemsForLabel,':');
		    
		    if(check){
			cout<<"ERROR: read_special_file: splitstring for Items["<<tmp_pos
			    <<"] : "<<Items[tmp_pos]<<"."<<endl;
			break;
		    }
		    
		    tmp_label_index = convert_typename_to_int(MP->get_Annotation_Label(),ItemsForLabel[0],number_of_annotation_labels);
		    if(tmp_label_index==-1){
			tmp_pos++;
			continue;
		    }else{
			feature_label_size =MP->get_Annotation_Label_size(tmp_label_index);
			for(i=0; i<feature_label_size; i++){
			    if(cmp_nocase(ItemsForLabel[1],MP->get_Annotation_Label_name(tmp_label_index,i))==0)
			    {
				found_keyword++;
				break;
			    }	        	    
			}
			if(found_keyword==0){
			    tmp_pos++;
			}
		    }
		}

		if (found_keyword != 0)
		{
		    // read information from special_format line

		    // get seq_name
		    strcpy(line_seq_names[number_of_line],Items[0]);
		    
		    // get source_name
		    strcpy(line_source_names[number_of_line],Items[1]);
		    
		    // get start position
		    line_start_positions[number_of_line] = atoi(Items[2]);
		    
		    // get end position
		    line_end_positions[number_of_line] = atoi(Items[3]);

		    if ((line_start_positions[number_of_line] <= line_end_positions[number_of_line]) && (!check))                             
		    {		       		     
			
			valid = true;

			if ((number_of_line == max_number_of_lines)&&(!feof(file_ptr)))
			{
			    cout << "   ERROR: read_special_file: max_number_of_lines ("
				 << max_number_of_lines << ") == number_of_line ("
				 << number_of_line << "). Increase the max value.\n" << flush;
			    check++;
			    break;
			}
			
		    }
		    else
		    {
			cout << "ERROR: read_special_file: special_line is corrupted, stop "
			     << "reading this special_file, return empty arrays.\n" << flush;
			
			cout << "   special_line[" << number_of_line << "] :\n" << flush;
			cout << "--------------------------------------------\n" << flush;
			cout << "   line_seq_names                       = " << line_seq_names[number_of_line]
			     << "\n   source_name                    = " << line_source_names[number_of_line]<<endl;
			for(i=0;i<number_of_annotation_labels; i++){
			    for(j=0;j<MP->get_Annotation_Label_size(i);j++){
				if(line_annotation_labels[number_of_line][i][j]==true){
				    cout<<MP->get_Annotation_Label_setname(i)<<" ="
					<<MP->get_Annotation_Label_name(i,j)
					<<" score : "<<line_scores[number_of_line][i][j]<<endl;
				}
			    }
			}				

			cout << "\n   start_position                 = " << line_start_positions[number_of_line]
			     << "\n   end_position                   = " << line_end_positions[number_of_line]
			     << flush;
			cout << "\n--------------------------------------------\n" << flush;
			check++;
			break;
		    }	

		    if(valid)
		    {
			// get annotation label
			
			// get series of annotation_labels
			while(1){
			    
			    if(check){
				break;
			    }
			    
			    if(tmp_pos>=NoOfItems){
				break;
			    }			
			    // end of annotation labels
			    if(!strcmp(Items[tmp_pos],"|")){
				tmp_pos++;
				break;
			    }
			    
			    check+=  splitstring(Items[tmp_pos], NoOfItemsForLabel, &ItemsForLabel,':');
			    if(check){
				cout<<"ERROR: read_special_file: splitstring for Items["<<tmp_pos
				    <<"] : "<<Items[tmp_pos]<<"."<<endl;
				break;
			    }
			    
			    tmp_label_index = convert_typename_to_int(MP->get_Annotation_Label(),ItemsForLabel[0],number_of_annotation_labels);
			    if(tmp_label_index==-1){
				tmp_pos++;
				continue;
			    }else{
				feature_label_size =MP->get_Annotation_Label_size(tmp_label_index);
				found_keyword = 0;
				for(i=0; i<feature_label_size; i++){
				    if(cmp_nocase(ItemsForLabel[1],MP->get_Annotation_Label_name(tmp_label_index,i))==0){
					if(line_annotation_labels[number_of_line][tmp_label_index][i]==true){
					    cout<<"Error:: read_special : "
						<<"same label has been defined at the same by position by the same source!"<<endl;
					    check++;
					    break;
					}else{
					    line_annotation_labels[number_of_line][tmp_label_index][i]=true;
					    found_keyword ++ ;
					}
					
					tmp_pos++;
					if(!strcmp(Items[tmp_pos],".")){
					    scores[number_of_line][tmp_label_index][i] = 0;
					    tmp_pos++;
					}else if(!MP->get_score_of_Annotation_Label(tmp_label_index)){
					    cout<<"Error:: read_special_file : score("<<Items[tmp_pos]<<") is set for non score label("<<MP->get_Annotation_Label_setname(tmp_label_index)<<")"<<endl;
					    check++;
					}else{
					    // check if the score is valid
					    double tmp_score = atof(Items[tmp_pos]);
					    
					    if((tmp_score<0)||(tmp_score>1))
					    {
						cout<<"Error:: read_special_file : score("<<Items[tmp_pos]
						    <<") out of range [0..1].\n";
						cout<<"Interrupted line : \n"
						    <<line<<endl;
						check++;
						tmp_score = 0;
					    }else if((tmp_score==0)&&(strcmp(Items[tmp_pos],"0")))
					    {
						cout<<"Error:: read_special_file : score("<<Items[tmp_pos]
						    <<") not valid, should be a numerical number in between 0 and 1.\n";
						cout<<"Interrupted line : \n"
						    <<line<<endl;
						check++;
						tmp_score = 0;
					    }
					    line_scores[number_of_line][tmp_label_index][i]=tmp_score;					     		  
					    tmp_pos++;
					}	        		    
					break;
				    }
				}
				if(found_keyword == 0){
				    tmp_pos++;
				}			
			    }
			}	
					
			number_of_line++;
		    }
		}
	    }
	}
	// close file
	fclose(file_ptr);

	// release memory
	if(Items){
	    for(i=0; i<Max_number_of_items; i++){
		if(Items[i]) delete [] Items[i];
		Items[i] = NULL;
	    }
	    delete[] Items;
	}
	Items = NULL;
	NoOfItems = 0;

	if(ItemsForLabel){
	    for(i=0; i<Max_number_of_items; i++)
	    {
		if(ItemsForLabel[i]) delete [] ItemsForLabel[i];
		ItemsForLabel[i] = NULL;
	    }
	    delete [] ItemsForLabel;
	}
	ItemsForLabel = NULL;
	NoOfItemsForLabel = 0;

	if(!check)
	{
	    int label_size = 0;
	    // initialize output variables	    
	    (*number_of_sources) = 0;
	    (*number_of_lines) = NULL;  // counts numbers of gtf_lines stored in the output arrays from each source   	    
	    check+=sort_lines_by_sources(//input 
		number_of_line,
		MP,
		&line_seq_names,
		&line_source_names,
		&line_annotation_labels,
		&line_scores,
		&line_start_positions,
		&line_end_positions,
		//output
		number_of_sources,
		number_of_lines);
	    
	    if(check){
		cout<<"ERROR : read_special_file : error in sort_lines_by_sources "<<endl;
	    }

	    if(!check)
	    {
		check+= initialize_labels((*number_of_sources),
					  (*number_of_lines),
					  MP,
					  seq_names,
					  source_names,
					  annotation_labels,
					  scores,
					  start_positions,
					  end_positions);
		if(check){
		    cout<<"ERROR : read_special_labels, initialize_labels for source "<<endl;
		}
		
		//set output variables
		int source_number = 0;
		int count = 0;
	    
		for(i=0; i<number_of_line; i++){
		    if(i!=0){
			if(strcmp(line_source_names[i-1],line_source_names[i])!=0){
			    source_number++;
			    if(source_number>=(*number_of_sources)){
				cout<<"Error:: read_special_file_for_sequence : source_number("<<source_number+1
				    <<") exist number of source : "<<(*number_of_sources)<<endl;
				strcpy((*source_names)[source_number],line_source_names[i]);
				count = 0;		     
			    }
			}
		    }else{
			strcpy((*source_names)[source_number],line_source_names[i]);
		    }
		    
		    strcpy((*seq_names)[source_number][count],line_seq_names[i]);
		    
		    for(j=0;j<number_of_annotation_labels;j++){
			label_size = MP->get_Annotation_Label_size(j);
			for(k=0;k<label_size;k++){
			    (*annotation_labels)[source_number][count][j][k] = line_annotation_labels[i][j][k];
			    (*scores)[source_number][count][j][k] = line_scores[i][j][k];
			}
		    }
	
		    (*start_positions)[source_number][count] = line_start_positions[i];
		    (*end_positions)[source_number][count]   = line_end_positions[i];		
		    count++;
		    
		    if(count>(*number_of_lines)[i]){
			cout<<"Error:: read_special_file_for_sequence: count : "<<count+1 
			    <<" exceed the limit of line("<<i<<") numebr_of_lines["<<i<<"] : "
			    <<(*number_of_lines)[i]<<".\n"<<flush;
			check++;
			break;
		    }
		}

		// if check != 0 delete the output information	
	
		if (check != 0)
		{
		    // release memory
		    check+=release_memory_for_labels(number_of_sources,
						     number_of_lines,
						     MP,
						     seq_names,
						     source_names,
						     annotation_labels,
						     scores,
						     start_positions,
						     end_positions);

		    if(check){
			cout<<"ERROR : read_special_labels, release_memory_for_labels for source "<<endl;
		    }	
		}
	    }
	}
    
	// delete temporary varibales

	if(line) delete[] line;
	line = NULL;

	check += release_memory_for_labels(&max_number_of_lines,
					   MP,
					   &line_seq_names,
					   &line_annotation_labels,
					   &line_scores,
					   &line_start_positions,
					   &line_end_positions);
	if(check){
	    cout<<"ERROR : read_special_labels, release_memory_for_labels for line "<<endl;
	}	

	// delete line_source_names
	if(line_source_names){
	    for(i=0; i<max_number_of_lines; i++){
		if(line_source_names[i]) delete[] line_source_names[i];
		line_source_names[i] = NULL;
	    }
	    delete[] line_source_names;
	}
	line_source_names = NULL;
    }else
    {
	// if error occurred during initial checks, set number_of_lines = 0
	// note that program which calls this function is responsible for deleting any memory
	// that was allocated for the output variables before calling this function	
	(*number_of_sources) = 0;
    }   
    return(check);
}

int apply_restriction_set_to_lines(// input
    const int max_number_of_sources,
    const int max_number_of_lines,
    model_parameters* const MP,
    // variables that define the subset to be looked at
    int*      const res_number_of_sources,
    int**     const res_number_of_lines,
    char****  const res_seq_names,        
    int***    const res_start_positions,  
    int***    const res_end_positions,    
    // variables that shall be restricted to above subset
    // variables are used as input, are modified and then used as output
    int*        const number_of_sources,
    int**       const number_of_lines,
    char****    const seq_names,
    char***     const source_names,
    bool*****   const annotation_labels,
    double***** const scores,
    int***      const start_positions,
    int***      const end_positions
    )
{
    int check = 0;
    if (max_number_of_sources < 1)
    {
	cout << "ERROR: apply_restriction_set_to_lines : max_number_of_sources (" 
	     << max_number_of_sources << ") < 1.\n" << flush;
	check++;
    }
    if (max_number_of_lines < 1)
    {
	cout << "ERROR: apply_restriction_set_to_lines : max_number_of_lines (" 
	     << max_number_of_lines << ") < 1.\n" << flush;
	check++;
    }
    if((*res_number_of_sources)<1)
    {
	cout << "ERROR: apply_restriction_set_to_lines : res_number_of_sources (" 
	     << (*res_number_of_sources) << ") < 1.\n" << flush;
	check++;
    }
    if ((*res_number_of_lines) ==NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : res_number_of_lines is NULL.\n"<<flush;
	check++;
    }
    if ((*res_seq_names) == NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : res_seq_names is NULL.\n" << flush;
	check++;
    }
    if ((*res_start_positions) == NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : res_start_positions is NULL.\n" << flush;
	check++;
    }
    if ((*res_end_positions) == NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : res_end_positions is NULL.\n" << flush;
	check++;
    }
    
    if ((*number_of_sources) < 1)
    {
	cout << "ERROR: apply_restriction_set_to_lines : number_of_sources (" 
	     << (*number_of_sources) << ") < 1.\n" << flush;
	check++;
    }
    if ((*number_of_lines) ==NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : number_of_lines is NULL.\n"<<flush;
	check++;
    }
    if ((*seq_names) == NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : seq_names is NULL.\n" << flush;
	check++;
    }
    if ((*source_names) == NULL)
    {
      cout << "ERROR: apply_restriction_set_to_lines : source_names is NULL.\n" << flush;
      check++;
    }
   
    if((*annotation_labels)==NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : annotation_labels is NULL.\n" << flush;
	check++;
    }

    if ((*start_positions) == NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : start_positions is NULL.\n" << flush;
	check++;
    }
    if ((*end_positions) == NULL)
    {
	cout << "ERROR: apply_restriction_set_to_lines : end_positions is NULL.\n" << flush;
	check++;
    }

    int i = 0;    
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    
    int number_of_annotation_labels = MP->get_Total_Number_of_Annotation_Labels();
    if (check == 0)
    {
	// check that lines for each source in restriction set have no mutual overlap 	
	for(i=0; i< (*res_number_of_sources); i++)
	{
	    for (j=0; j<(*res_number_of_lines)[i]; j++)
	    {
		for (k=0; k<(*res_number_of_lines)[i]; k++)
		{
		    if (j != k)
		    {
			// [start[i],end[i]] and [start[j],end[j]] have an overlap
			// i.e. if they do not have no overlap
			
			if ((strcmp((*res_seq_names)[i][j], (*seq_names)[i][k]) == 0)  &&
			    (! (((*res_start_positions)[i][j] > (*end_positions)[i][k]) ||
				((*end_positions)[i][j] < (*start_positions)[i][k]))))
			{
			    cout << "ERROR: apply_restriction_set_to_gtf_lines: there exists an overlap between "
				 << "line j (" << j << ") = [" << (*start_positions)[i][j]
				 << ", " << (*end_positions)[i][k] << "] and line k (" << k
				 << ") = " << (*start_positions)[i][k]
				 << ", " << (*end_positions)[i][k] << "] of the restriction set in source i("<<i<<").\n" << flush;
			    
			    check++;
			    break;
			}
			
		    }
		}
	    }
	}
	
	if (check == 0)
	{
	    bool same = true;
	    // check if the restriction set is the same as the original set
	    if((*number_of_sources)!=(*res_number_of_sources)){
		same = false;
	    }else{
		for(i=0; i<(*number_of_sources); i++){
		    if(!same){
			break;
		    }
		    if((*number_of_lines)[i]!=(*res_number_of_lines)[i]){
			same = false;
			break;
		    }else{
			for(j=0;j<(*number_of_lines)[i];j++){
			    if((*start_positions)[i]!=(*res_start_positions)[i]){
				same = false;
				break;
			    }
			}
		    }
		}
	    }
	    
	    if (!same)
	    {
		// allocate memory for the variables which will hold the new lines

		int* new_number_of_lines = new int[max_number_of_sources];
		for(i=0; i<max_number_of_sources; i++){
		    new_number_of_lines[i] = max_number_of_lines;
		}

		char*** new_seq_names;
		char** new_source_names;	    
		bool**** new_annotation_labels;    
		double**** new_scores;
		int** new_start_positions;
		int** new_end_positions;
		
		check+= initialize_labels(max_number_of_sources,
					  new_number_of_lines,
					  MP,
					  &new_seq_names,
					  &new_source_names,
					  &new_annotation_labels,
					  &new_scores,
					  &new_start_positions,
					  &new_end_positions);
		
		int  new_number_of_sources = 0;

		// determine set of gtf_lines that are restricted to the area defined
		// by the intervals of the restriction set
		
		bool set_line = false; // count if the is line in the restriction set from a source
		int scount              = 0; // counts lines in new set of gtf_lines
		int lcount                = 0;

		for(i=0; i<(*res_number_of_sources); i++)
		{
		    lcount = 0;
		    set_line = false;
		    
		    for (j=0; j<(*number_of_lines)[i]; j++)
		    {		    
			for (k=0; k<(*res_number_of_lines)[i]; k++)
			{
			    // notation:
			    // [j] := [j--j] := [start[j],end[j]]
			    // [k] := [k--k] := [start[k],end[k]]
			
			    // [j] and [k] have an overlap
			    // i.e. if they do not have no overlap
			
			    if ((strcmp((*seq_names)[i][j], (*res_seq_names)[i][k]) == 0) &&
				(! (((*start_positions)[i][j] > (*res_end_positions)[i][k]) ||
				    ((*end_positions)[i][j] < (*res_start_positions)[i][k]))))
			    {
				if (((*start_positions)[i][j] > (*res_start_positions)[i][k]) &&
				    ((*end_positions)[i][j]   < (*res_end_positions)[i][k]))
				{
				    // if
				    // 1.) case: [k---[j---j]---k] => [j] does not overlap with any other [k]
				    //           because the [k]'s don't overlap
				    
				    // [j] can remain as it is
				
				    strcpy(new_seq_names[scount][lcount],(*seq_names)[i][j]);
				    strcpy(new_source_names[scount],(*source_names)[i]);
				    
				    for(l=0; l<number_of_annotation_labels; l++){
					for(m=0;m<MP->get_Annotation_Label_size(m);m++){
					    new_annotation_labels[scount][lcount][l][m]= (*annotation_labels)[i][j][l][m];
					    new_scores[scount][lcount][l][m] = (*scores)[i][j][l][m];
					}
				    }
				    
				    new_start_positions[scount][lcount] = (*start_positions)[i][j];
				    new_end_positions[scount][lcount]   = (*end_positions)[i][j];   
				    set_line = true;
				    lcount++;
				    
				}
				else
				{
				    // if
				    // 2.) case: [j---[k---j]---k] ||
				    // 3.) case: [k---[j---k]---j] ||
				    // 4.) case: [j---[k---k]---j]
				    
				    // adjust start and/or end of [j] by setting
				    // start[j] = max[start[j],start[k]] and
				    // end[j]   = min[end[j], end[k]]
				    
				    strcpy(new_seq_names[scount][lcount],(*seq_names)[i][j]);
				    strcpy(new_source_names[scount],(*source_names)[i]);
				    for(l=0; l<number_of_annotation_labels; l++){
					for(m=0; m<MP->get_Annotation_Label_size(l); m++){
					    new_annotation_labels[scount][lcount][l][m], (*annotation_labels)[i][j][l][m];
					}
				    }
				   
				    new_start_positions[scount][lcount] = max((*start_positions)[i][j],(*res_start_positions)[i][k]);
				    new_end_positions[scount][lcount]   = min((*end_positions)[i][j],(*res_end_positions)[i][k]);
				    
				    set_line = true;
				    lcount++;				
				}
			    } // [j] and [k] overlap
			} // loop over k lines of restriction set
		    } // loop over j lines of original set
		    if(set_line){
			scount++;
		    }
		    new_number_of_lines[i] = lcount;
		} // loop over i sources
		new_number_of_sources = scount;
		
		if(new_number_of_sources!=(*number_of_sources))
		{
		    check+= release_memory_for_labels(number_of_sources,
						      number_of_lines,
						      MP,
						      seq_names,
						      source_names,
						      annotation_labels,
						      scores,
						      start_positions,
						      end_positions);
						      				
		    (*number_of_sources) = new_number_of_sources;
		    for(i=0; i<(*number_of_sources);i++){
			(*number_of_lines)[i] = new_number_of_lines[i];
		    }
		
		    //allocate memory
		    check+= initialize_labels((*number_of_sources),
					      (*number_of_lines),
					      MP,
					      seq_names,
					      source_names,
					      annotation_labels,
					      scores,
					      start_positions,
					      end_positions);
		}else{
		    for(i=0; i<(*number_of_sources); i++)
		    {
			if(new_number_of_lines[i]!=(*number_of_lines)[i])
			{
			    check+= release_memory_for_labels(&((*number_of_lines)[i]),
							      MP,
							      &((*seq_names)[i]),
							      &((*annotation_labels)[i]),
							      &((*scores)[i]),
							      &((*start_positions)[i]),
							      &((*end_positions)[i]));

			    (*number_of_lines)[i] = new_number_of_lines[i];
			    
			    check+= initialize_labels((*number_of_lines)[i],
						      MP,
						      &((*seq_names)[i]),
						      &((*annotation_labels)[i]),
						      &((*scores)[i]),
						      &((*start_positions)[i]),
						      &((*end_positions)[i]));
			    
			}
		    }
		}
		
                // copy new lines into original lines	
		for(i=0; i<(*number_of_sources); i++)
		{
		    strcpy((*source_names)[i], new_source_names[i]);
		    for (j=0; j<(*number_of_lines)[i]; j++)
		    {
			strcpy((*seq_names)[i][j],new_seq_names[i][j]);		   
			for(k=0; k<number_of_annotation_labels; k++){
			    for(l=0; l<MP->get_Annotation_Label_size(k); l++){
				(*annotation_labels)[i][j][k][l] = new_annotation_labels[i][j][k][l];
				(*scores)[i][j][k][l]            = new_scores[i][j][k][l];
			    }
			}
		
			(*start_positions)[i][j]    = new_start_positions[i][j];
			(*end_positions)[i][j]      = new_end_positions[i][j]; 
		    }
		}
			
		// delete memory of new lines
		
		check+= release_memory_for_labels(&new_number_of_sources,
						  &new_number_of_lines,
						  MP,
						  &new_seq_names,
						  &new_source_names,
						  &new_annotation_labels,
						  &new_scores,
						  &new_start_positions,
						  &new_end_positions);	
		
	    } // if restriction set and gtf set are not identical
	} // if check == 0 (if intervals of restriction set do not overlap)
    } // if initial checks on input and output variables o.k.
    
    // delete gtf-lines if function was not successfull
    
    if (check != 0)
    {
	check+= release_memory_for_labels(number_of_sources,
					  number_of_lines,
					  MP,
					  seq_names,
					  source_names,
					  annotation_labels,
					  scores,
					  start_positions,
					  end_positions);
    }    
    return(check);
}
