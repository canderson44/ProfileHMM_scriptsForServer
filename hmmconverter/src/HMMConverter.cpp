/*
  Author: Philip Lam
  Date: 16/06/2007
  Purpose: the main program for HMMConverter

  RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/HMMConverter.cpp,v 1.3 2008/12/14 10:36:54 natural Exp $
*/

#include <fstream>
#include <iostream>
#include <ctype.h>

#include "tinyxml.h"
#include "define.h"
#include "models.h"
#include "hmm_state.h"
#include "sequence.h"
#include "input_from_files.h"
#include "blastn.h"
#include "tube.h"
#include "scoring_functions.h"
#include "evaluation.h"
#include "match.h"
#include "transitionprob.h"
#include "emissionprob.h"
#include "model_parameters.h"

using namespace std;

int main(int argc, char** argv){

    // argv[1] Input directory
    // argv[2] Output directory
    // argv[3] mass xml file name

    int check = 0;
    // read from argv
    
    char cIn_dir[Max_line_length];   
    char cOut_dir[Max_line_length];
    char cIn_xml[Max_line_length];
    
    strcpy(cIn_dir, argv[1]);

    strcpy(cOut_dir, argv[2]);
    
    strcpy(cIn_xml, cIn_dir);
    strcat(cIn_xml, argv[3]);

    check+= set_executable_path(cIn_dir);
    if(check){
	cout<<"ERROR: set_executable_path("<<cIn_dir<<"). Abort.\n"<<flush;
	return check;
    }

    check+= set_temporary_path(cIn_dir);
    if(check){
	cout<<"ERROR: set_temporary_path("<<cIn_dir<<"). Abort.\n"<<flush;
	return check;
    }

    //char* cIn_xml = "pairhmm_full.xml";

    //cIn_xml = "pairhmm_full.xml";

    model_parameters MP(cIn_xml,check);
 

    //if(!check)
	//MP.print(cout);
    TransitionProb TranProb(cIn_dir,cIn_xml,check);

    //if(!check)
//	TranProb.print(cout);
    EmissionProb EmitProb(cIn_dir,cIn_xml,&MP,check);

    //  if(!check)
//	EmitProb.print(cout);

    Hmm real(&MP);
    check+=get_hmm(cIn_xml,&real,&MP,&TranProb,&EmitProb);
 
    if(check){
	cout<<"Error ! get_hmm. Abort.\n"<<flush;
	return check;
    }

    Hmm mirror_of_real(real,-1,0);

    check += mirror_of_real.check_consistency();
    if(check){
	cout<<"Error! mirror_of_real did not pass the consistency check. Abort.\n" <<flush;
	return check;
    }

    check+=mirror_of_real.calculate_scores_from_probs();
    if(check){
	cout <<"Error : mirror_of_model: function pairhmm:: calculate_scores_from_probs failed. Abort.\n"<<flush;
	return check;
    }

    /*
    check+=TranProb.get_parameters_for_training(cIn_xml,&MP,&real);
    if(!check)
    {
	TranProb.print(cout);
    }
    */

    /*
    check+=EmitProb.get_parameters_for_training(cIn_xml,&MP,&real);
    if(!check)
    {
	EmitProb.print(cout);
    } 
    */
    //real.print(cout);    

    // Analyse Sequence
    
    // ####################################################################

    ofstream Ofile;

    int NoOfItems = 0;
    char** Items = new char*[Max_number_of_items];
    for(int i =0; i<Max_number_of_items; i++)
    {
	Items[i] = new char[Max_word_length];
	strcpy(Items[i]," ");
    }

    /*
    char* TmpcIn_annotation=NULL;
    char* TmpcIn_sequence = NULL;
    char* TmpcIn_output = NULL;
    char* TmpcIn_transition_prob = NULL;
    char* TmpcIn_emission_prob = NULL;
    int number_of_output = 0;
    char* TmpcIn_tube = NULL;
    unsigned long long int LILimit_volume = 0;
    int Algorithm = 0;
    int Training_alg = 0;
    int IRadius = 0;          // tube radius in bp

    int IOutputFormat = 0;
    int Ispecial = 0;
    */

    // parameters for sequence decoding
    char* TmpcIn_decoding_annotation_file = NULL;
    char* TmpcIn_decoding_sequence_file = NULL;
    char* TmpcIn_decoding_tube_file = NULL;
    char* TmpcIn_decoding_output_file = NULL;
    int number_of_output = 0;
    int decoding_algorithm = -1;
    int radius = 0;
    int Ispecial = 0;
    unsigned long long int max_decoding_volume =0;

    check+= get_sequence_decoding_parameters(cIn_xml,
					     &TmpcIn_decoding_annotation_file,
					     &TmpcIn_decoding_sequence_file,
					     &TmpcIn_decoding_tube_file,
					     decoding_algorithm,
					     max_decoding_volume,
					     radius,
					     &TmpcIn_decoding_output_file);
	         
    if(check)
    {
	cout<<"Error: get_sequence_decoding_parameters : abort.\n";
	return check;
    }
    
    if((!MP.is_PairHMM())
       &&((decoding_algorithm==1)||(decoding_algorithm==2)||(decoding_algorithm==3)))
    {
	cout<<"Error : the model is not a PairHMM, can not use algorithm("
	    <<decoding_algorithm<<").\n";
	check++;
	return check;
    }

    // check if the model has the special emission structure 
    // and the prior information file is defined
    if((MP.is_SpecialEmit())&&(TmpcIn_decoding_annotation_file!=NULL))
    {
	Ispecial = 1;
    }else{
	Ispecial = 0;
    }

    //cout<<"\n-------------------------------------------------------------------\n\n";

    // parameters for parameter training
    char* TmpcIn_training_annotation_file = NULL;
    char* TmpcIn_training_sequence_file = NULL;
    char* TmpcIn_training_tube_file = NULL;
    int training_algorithm = -1;
    unsigned long long int max_training_volume = 0;
    double threshold = 0;
    int Maxiter = 0;
    int Tradius = 0;
    int SamplePaths = 0;
    char* TmpcIn_transition_prob_file = NULL;
    char* TmpcIn_emission_prob_file = NULL;
    char* TmpcIn_XMLfile = NULL;
    int trainingcheck = 0;
 
    check+= get_parameter_training_parameters(cIn_xml,
					      &TmpcIn_training_annotation_file,
					      &TmpcIn_training_sequence_file,
					      &TmpcIn_training_tube_file,
					      training_algorithm,
					      max_training_volume,
					      Tradius,
					      threshold,
					      Maxiter,					      
					      SamplePaths,
					      &TmpcIn_transition_prob_file,
					      &TmpcIn_emission_prob_file,
					      &TmpcIn_XMLfile);

    if(check)
    {
	cout<<"Error: get_parameter_training_parameters : abort.\n";
	return check;
    }
    
    if(training_algorithm<0)
    {
	cout<<"No training algorithm needed to be run"<<endl;
    }else 
    {
	// exact path for input files
	
	char cIn_training_sequence_file[Max_line_length];
	if(TmpcIn_training_sequence_file)
	{
	    strcpy(cIn_training_sequence_file,cIn_dir);
	    strcat(cIn_training_sequence_file,TmpcIn_training_sequence_file);
	    delete[] TmpcIn_training_sequence_file;
	    TmpcIn_training_sequence_file = NULL;
	}else{
	    strcpy(cIn_training_sequence_file," ");
	}
	
	char cIn_training_tube_file[Max_line_length];
	if(TmpcIn_training_tube_file)
	{
		if(!((!cmp_nocase(TmpcIn_training_tube_file,"TBlastx"))||
			(!cmp_nocase(TmpcIn_training_tube_file,"Blastn"))))
		{
			strcpy(cIn_training_tube_file,cIn_dir);
			strcat(cIn_training_tube_file,TmpcIn_training_tube_file);
			delete[] TmpcIn_training_tube_file;
			TmpcIn_training_tube_file = NULL;
		}else{
			strcpy(cIn_training_tube_file,TmpcIn_training_tube_file);
			delete[] TmpcIn_training_tube_file;
			TmpcIn_training_tube_file = NULL;
		}
	}else{
	    strcpy(cIn_training_tube_file, " ");
	}
	
	char cIn_training_annotation_file[Max_line_length];
	if(TmpcIn_training_annotation_file)
	{
	    strcpy(cIn_training_annotation_file,cIn_dir);
	    strcat(cIn_training_annotation_file,TmpcIn_training_annotation_file);
	    delete[] TmpcIn_training_annotation_file;
	    TmpcIn_training_annotation_file = NULL;
	}else{
	    strcpy(cIn_training_annotation_file," ");
	}
	
	if(training_algorithm==0)
	{	    	    	    
	    trainingcheck+= real.BaumWelchTraining(cIn_xml,				       
						   &TranProb,
						   &EmitProb,
						   &MP,
						   cIn_training_sequence_file,
						   cIn_training_annotation_file,
						   cIn_training_tube_file,
						   max_training_volume,
						   Tradius,
						   Maxiter,
						   threshold);	

	    if(trainingcheck)
	    {
		cout<<"Error occur in function BaumWelchTraining, training abort. "<<endl;
	    }	    
	}else if(training_algorithm==1)
	{
	    // Posterior training
	    trainingcheck+= real.PosteriorTraining(cIn_xml,				       
						   &TranProb,
						   &EmitProb,
						   &MP,
						   cIn_training_sequence_file,
						   cIn_training_annotation_file,
						   cIn_training_tube_file,
						   max_training_volume,
						   Tradius,
						   Maxiter,
						   SamplePaths);	

	    if(trainingcheck)
	    {
		cout<<"Error occur in function PosteriorTraining, training abort. "<<endl;
	    }
	}else if(training_algorithm==2){
	    
	    // Viterbi training
	    trainingcheck+= real.ViterbiTraining(cIn_xml,				       
						 &TranProb,
						 &EmitProb,
						 &MP,
						 cIn_training_sequence_file,
						 cIn_training_annotation_file,
						 cIn_training_tube_file,
						 max_training_volume,
						 Tradius,
						 Maxiter);	

	    if(trainingcheck)
	    {
		cout<<"Error occur in function ViteriTraining, training abort. "<<endl;
	    }
	}else{
	    cout<<"Error: training_algorithm("<<training_algorithm<<") out of range. "
		<<"No training algorithm is run"<<endl;
	    trainingcheck++;
	}
	
	// print result
	if(!trainingcheck)
	{
	    // print result for FTP
	    if(TmpcIn_transition_prob_file)
	    {
		char cIn_transition_prob_file[Max_line_length];
		strcpy(cIn_transition_prob_file,cOut_dir);
		strcat(cIn_transition_prob_file,TmpcIn_transition_prob_file);	      
		if(TmpcIn_transition_prob_file) delete [] TmpcIn_transition_prob_file;
		TmpcIn_transition_prob_file = NULL;
		ofstream Ofile;
		Ofile.open(cIn_transition_prob_file);
		real.transition_prob_output(Ofile,&TranProb);		
	    }
	    
	    // print result for emission prob
	    if(TmpcIn_emission_prob_file)
	    {
		char cIn_emission_prob_file[Max_line_length];
		strcpy(cIn_emission_prob_file,cOut_dir);
		strcat(cIn_emission_prob_file,TmpcIn_emission_prob_file);	      
		if(TmpcIn_emission_prob_file) delete [] TmpcIn_emission_prob_file;
		TmpcIn_emission_prob_file = NULL;
		ofstream Ofile;
		Ofile.open(cIn_emission_prob_file);
		real.emission_prob_output(Ofile,&MP,&EmitProb);		
	    }
	    
	    // print result for transition prob
	    if(TmpcIn_XMLfile)
	    {
		char cIn_XMLfile[Max_line_length];
		strcpy(cIn_XMLfile,cOut_dir);
		strcat(cIn_XMLfile,TmpcIn_XMLfile);	      
		if(TmpcIn_XMLfile) delete [] TmpcIn_XMLfile;
		TmpcIn_XMLfile = NULL;
		ofstream Ofile;
		Ofile.open(cIn_XMLfile);
		real.XML_output(cIn_xml,Ofile,&MP,&TranProb);		
	    }		
	}
    }
    
    if(decoding_algorithm<0)
    {
	
	cout<<"No decoding algorithm needed to be run"<<endl;
	return check; // no sequence analysis is needed
    }  

    /*
    strcpy(cSource,MP.get_Model_Name());
    strcat(cSource,"Alg");
    char cAlgorithm[2];
    cAlgorithm[0]= Algorithm+'0'-0;
    cAlgorithm[1]='\0';
    strcat(cSource,cAlgorithm);
    */
    
    char cIn_decoding_sequence_file[Max_line_length];
    if(TmpcIn_decoding_sequence_file)
    {
	strcpy(cIn_decoding_sequence_file,cIn_dir);
	strcat(cIn_decoding_sequence_file,TmpcIn_decoding_sequence_file);
	delete[] TmpcIn_decoding_sequence_file;
	TmpcIn_decoding_sequence_file = NULL;
    }

    char cIn_decoding_tube_file[Max_line_length];
    if(TmpcIn_decoding_tube_file)
    {
	strcpy(cIn_decoding_tube_file,cIn_dir);
	strcat(cIn_decoding_tube_file,TmpcIn_decoding_tube_file);
	delete[] TmpcIn_decoding_tube_file;
	TmpcIn_decoding_tube_file = NULL;
    }
    
    char cIn_decoding_annotation_file[Max_line_length];
    if(TmpcIn_decoding_annotation_file)
    {
	strcpy(cIn_decoding_annotation_file,cIn_dir);
	strcat(cIn_decoding_annotation_file,TmpcIn_decoding_annotation_file);
	delete[] TmpcIn_decoding_annotation_file;
	TmpcIn_decoding_annotation_file = NULL;
    }

    // determine the number of output 
    if(MP.get_Total_Number_of_Annotation_Labels()>0)
    {
	number_of_output = 2;
    }else{
	number_of_output = 1;
    }

    char** cIn_decoding_output_files = new char*[number_of_output];
    for(int i=0; i<number_of_output; i++){
	cIn_decoding_output_files[i]=new char[Max_line_length];
    }
    if(TmpcIn_decoding_output_file){
	check+=splitstring(TmpcIn_decoding_output_file,NoOfItems,&Items,'.');
	if(check)
	{
	    cout<<"Output file not declared in right format : "
		<<TmpcIn_decoding_output_file<<endl;
	}	
	if(NoOfItems!=2)
	{
	    cout<<"Output file not declared in right format : "
		<<TmpcIn_decoding_output_file<<endl;
	    check++;
	}
	if(!check)
	{
	    for(int i=0; i<number_of_output; i++){
		strcpy(cIn_decoding_output_files[i],cOut_dir);
		strcat(cIn_decoding_output_files[i],Items[0]);	    

	    }		
	    if(number_of_output>1)
	    {
		strcat(cIn_decoding_output_files[0],"_1");
		strcat(cIn_decoding_output_files[0],".");
		strcat(cIn_decoding_output_files[0],Items[1]);
		if(number_of_output>1)
		{
		    strcat(cIn_decoding_output_files[1],"_2");
		    strcat(cIn_decoding_output_files[1],".");
		    strcat(cIn_decoding_output_files[1],Items[1]);
		}
	    }else{
		strcat(cIn_decoding_output_files[0],".");
		strcat(cIn_decoding_output_files[0],Items[1]);
	    }
	}
	delete[] TmpcIn_decoding_output_file;
	TmpcIn_decoding_output_file = NULL;
    }	
   
    FILE* fIn_sequence = fopen(cIn_decoding_sequence_file,"rt");

    if(!fIn_sequence){
	cout<<"ERROR: can not open file : "<<cIn_decoding_sequence_file<<endl;
	check++;
	return check;
    }

    int iSPCount = 0;

    while(!feof(fIn_sequence))
    {
    
	Sequence sX, sY;	

	if(MP.is_PairHMM())
	{
	    check+= get_sequence_pair(fIn_sequence, &sX, &sY, &MP);
	    if(check){
		cout<<"Error: function sequence:: get_sequence_pair failed. Abort.\n"<<flush;
		return check;
	    }
	
	}else{
	    check+= get_single_sequence(fIn_sequence,&sX,&MP);
	    if(check){
		cout<<"Error: function sequence:: get_single_sequence failed. Abort.\n"<<flush;
		return check;
	    }
	}
		
	if(Ispecial)
	{	
	    check += sX.set_annotation(cIn_decoding_annotation_file, &MP);
	    
	    if (check != 0) 
	    {
		cout << "ERROR: Sequence::set_annotation for x sequence of sequence pair " 
		     << iSPCount << ", "<< " and annotation file " 
		     << cIn_decoding_annotation_file << ".\n" << flush;
		break;
	    }
	    
	    //sX.print_annotation_linewise(cout,&MP);
	    
	    check += set_scores_for_annotated_sequence(&sX,
						       0, // (0 = x, 1 = y)
						       &real);						 
	    //sX.print_emission_scores_with_sequence_info(cout,&MP);
	    
	    if (check != 0) {
		cout << "ERROR: set_scores_for_annotated_sequence for x sequence of sequence pair " 
		     << iSPCount << " and annotation file " 
		     << cIn_decoding_annotation_file << ".\n" << flush;
		break;
	    }
	    
	    
	   
	    if(MP.is_PairHMM())
	    {
		//check += sY.set_constraint(1);

		//if(check)
		//{
		//    cout<<"ERROR: sY set_constraint"<<endl;
		//    break;
		//}

		check += sY.set_annotation(cIn_decoding_annotation_file, &MP);
		
		if (check != 0) 
		{
		    cout << "ERROR: Sequence::set_annotation for x sequence of sequence pair " 
			 << iSPCount << ", "<< " and annotation file " 
			 << cIn_decoding_annotation_file << ".\n" << flush;
		    break;
		}
		
		//sY.print_annotation_linewise(cout,&MP);
		
		check += set_scores_for_annotated_sequence(&sY,
							   1, // (0 = x, 1 = y)
							   &real);				      
 	
		//sY.print_emission_scores_with_sequence_info(cout,&MP);
	  	
		if (check != 0) {
		    cout << "ERROR: set_scores_for_annotated_sequence for x sequence of sequence pair " 
			 << iSPCount << " and annotation file " 
			 << cIn_decoding_annotation_file << ".\n" << flush;
		    break;
		}
	    }
	}
	// sX.print();
	//sY.print();	

	Tube<int> tTube;

	
	if((decoding_algorithm==1)||(decoding_algorithm==2)){
	    
	    Tube<int> tBlast_tube;
	    int iCount_matches, iCount_matches_new, iI;
	    unsigned long long int dTube_volume, dViterbi_volume;
	    Match* mMatches = NULL;
	    Match* mMatches_new = NULL;
	    
	    if(training_algorithm ==1 ){
		// get blastn matches

		cout<<"get blastn matches\n"<<flush;
		check += get_blastn_results(&sX,&sY, &iCount_matches, &mMatches, &MP);
		
		if(check){
		    cout<<"ERROR: get_blastn_results: abort.\n"<<flush;
		    break;
		}
	    }else{

		// get tblastx matches

		cout<<"get tblastx matches \n"<<flush;
		check+= get_tblastx_results(&sX,&sY,&iCount_matches,&mMatches,&MP);	
		if(check){
		    cout<<"ERROR: get_tblastx_results: abort.\n"<<flush;
		    break;
		}
	    }

	    check+= find_maximal_subset_of_matches(iCount_matches, mMatches, &iCount_matches_new, &mMatches_new);
	
	    if(check){
		cout<<"ERROR: find_maximal_subset_of_matches: abort.\n"<<flush;
		break;
	    }	    
	      
	    tBlast_tube = convert_matches_to_tube(sX.length(),sY.length(),iCount_matches_new,mMatches_new, radius);
	    dTube_volume= 0;

	    for(iI=0;iI<sX.length();iI++){
		dTube_volume+=(tBlast_tube.GetElement(iI,1)-tBlast_tube.GetElement(iI,0)+1);
	    }

	    dTube_volume += real.get_number_of_states();

	    dViterbi_volume = sX.length()*sY.length();
	    dViterbi_volume *= real.get_number_of_states();

	    if(mMatches_new) delete [] mMatches_new;
	    mMatches_new=NULL;

	    if(mMatches) delete []mMatches;
	    mMatches = NULL;

	    tTube = tBlast_tube;

	}

	if(decoding_algorithm==3){
	    Tube<int> tUser_tube;
	    int iI;
	    unsigned long long int dTube_volume, dViterbi_volume;
	    
	    cout<<"get tube \n"<<flush;

	    check+= get_tube_from_file(cIn_decoding_tube_file,sX.get_ac(),sY.get_ac(),sX.length(),sY.length(),&tUser_tube);

	    if(tUser_tube.Empty())
	    {
		cout<<"No tube is defined for the sequence pair"<<endl;
	    }
	    
	    /*
	    dTube_volume= 0;
	    	    
	    for(iI=0;iI<sX.length();iI++){
		dTube_volume+=(tUser_tube.GetElement(iI,1)-tUser_tube.GetElement(iI,0)+1);
	    }
	    
	    dTube_volume += real.get_number_of_states();

	    dViterbi_volume = sX.length()*sY.length();
	    dViterbi_volume *= real.get_number_of_states();
	    */

	    tTube = tUser_tube;

	}
	    
	// normal viterbi and Hirschberg

	int tried_hirschberg_empty = 0;

	while(1)
	{
	    if(!check){
		if(training_algorithm==0){
		    tried_hirschberg_empty =1;
		}
 	    
		check += real.Hirschberg_tube(mirror_of_real, sX, sY, max_decoding_volume, tTube); 
		cout <<"after real.Hirschberg_tube: check = "<<check<<"\n"<<flush;
	
		if (!check) {		      
		    cout << "log_odds score    = " << real.get_score() << "\n" << flush;
		}
		break;
	    }else{
		if(tried_hirschberg_empty==0){
		
		    cout<<"ERROR: Hirschberg_tube failed.\n";
		    cout<<"NOTE: starting Hirschberg algorithm with empty tube.\n";

		    tTube.Reset();
		    tried_hirschberg_empty=1;
		    check=0;
		}else{
		    cout<<"ERROR: Hirschberg_tube failed with empty tube.\n";
		    break;
		}
	    }
	}
	
	// print viterbi result
	for(int i=0; i<number_of_output; i++)
	{
	
	    ofstream Ofile;

	    if(iSPCount==0)
	    {
	
		Ofile.open(cIn_decoding_output_files[i]);
		  
	    }else{
		
		Ofile.open(cIn_decoding_output_files[i],ios::app);
		
	    }
	    if(i==0)
	    {
		real.output1(Ofile,&MP,sX.get_ac(),sY.get_ac());
	    }else{
		real.output2(Ofile,&MP,sX.get_ac(),sX.length(),sY.get_ac(),sY.length());
	    }	    
	}
	Ofile.close();
	
	iSPCount++;
    } // while end of file with sequences not reached

    // free memory
       
    if(TmpcIn_decoding_annotation_file) delete [] TmpcIn_decoding_annotation_file;
    TmpcIn_decoding_annotation_file = NULL;

    if(TmpcIn_decoding_sequence_file) delete [] TmpcIn_decoding_sequence_file;
    TmpcIn_decoding_sequence_file = NULL;

    if(TmpcIn_decoding_tube_file) delete[] TmpcIn_decoding_tube_file;
    TmpcIn_decoding_tube_file = NULL;
       
    if(cIn_decoding_output_files)
    {
	for(int i =0; i<number_of_output; i++){
	    if(cIn_decoding_output_files[i]) delete[] cIn_decoding_output_files[i];
	    cIn_decoding_output_files[i] = NULL;
	}
	delete[] cIn_decoding_output_files;
    }
    cIn_decoding_output_files = NULL;
    
    if(TmpcIn_decoding_output_file) delete[] TmpcIn_decoding_output_file;
    TmpcIn_decoding_output_file = NULL;

    if(TmpcIn_training_annotation_file) delete [] TmpcIn_training_annotation_file;
    TmpcIn_training_annotation_file = NULL;

    if(TmpcIn_training_sequence_file) delete [] TmpcIn_training_sequence_file;
    TmpcIn_training_sequence_file = NULL;

    if(TmpcIn_training_tube_file) delete[] TmpcIn_training_tube_file;
    TmpcIn_training_tube_file = NULL;

    if(TmpcIn_transition_prob_file) delete[] TmpcIn_transition_prob_file;
    TmpcIn_transition_prob_file = NULL;

    if(TmpcIn_emission_prob_file) delete[] TmpcIn_emission_prob_file;
    TmpcIn_emission_prob_file = NULL;
    
    if(TmpcIn_XMLfile) delete[] TmpcIn_XMLfile;
    TmpcIn_XMLfile = NULL;

    number_of_output = 0;
    
    if(Items)
    {
	for(int i=0; i<Max_number_of_items; i++)
	{
	    if(Items[i]) delete [] Items[i];
	    Items[i] = NULL;
	}
	if(Items) delete [] Items;
	Items = NULL;
    }

    NoOfItems = 0;
    if(executable_path) delete [] executable_path;
    executable_path = NULL;
    if(temporary_path) delete[] temporary_path;
    temporary_path = NULL;

    return check;
}
