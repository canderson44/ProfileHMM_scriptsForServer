/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define the emissionprob class

   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/emissionprob.cpp,v 1.4 2008/12/14 10:39:23 natural Exp $
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "emissionprob.h"

using namespace std;

//For emission_prob class

int EmissionProb:: set_FEP_from_file(const char* In_dir,
				     const char* xmlfile,
				     model_parameters* const MP)
{
    ifstream instream;
    int check = 0;
    int TNumOfEProb = 0;
    int size = 0;
    int dim = 0;
    int LCount = 0;
    int alphint = 0;
    	    
    int NoOfItems = 0;
    char** Items = new char*[Max_number_of_items];
    for(int i=0; i<Max_number_of_items; i++)
    {
	Items[i] = new char[Max_word_length];
	strcpy(Items[i]," ");
    }
    
    Prob prob = 0;
    Prob Tprob = 0;
    double pseudoprob = 0;
    char* line = new char[Max_line_length];
    char* filename = new char[Max_word_length];
    TiXmlNode* Node = 0;
    TiXmlNode* ENode = 0;
    TiXmlElement* EElement = 0;
    TiXmlElement* LElement = 0;
    int NumOfTokens=0;
    char** tokens=NULL;

    TiXmlDocument doc(xmlfile);
    bool loadOkay = doc.LoadFile();
    
    if(!loadOkay){
	cout<<"Error: set_emissiont_prob_from_file : cannot open file "<<xmlfile<<".\n"<<flush;
	check++;
	return check;
    }
    
    Node = doc.FirstChild("HMMConverter");
    assert(Node);
    Node = Node->FirstChild("model");
    assert(Node);
    ENode = Node->FirstChild("Emission_Probs");
    assert(ENode);
    EElement = Node->FirstChildElement("Emission_Probs");
    assert(EElement);

    // set tname of EP
    FEP_tname = Nullstrcpy(EElement->Attribute("id"),check);
    if(check){
	cout<<"Error:: set_emission_prob_from_file : "
	    <<"error in Nullstrcpy, id Attribute inf Emission_Probs tag : "
	    <<EElement->Attribute("id")<<" is invalid"<<endl;
    }

	
    //set size of EP
    EElement->Attribute("size",&size);
    if(size<=0){
	cout<<"Error: set_emission_prob_from_file :: size : "<<size<<endl;
	check++;
    }
    
    if(!check)
    {
	FEPsize = size;
	FEP = new FEProb[FEPsize];
	
	for(int i=0; i<FEPsize; i++)
	{
	    FEP[i].name = NULL;
	    FEP[i].dim = 0; 
	    FEP[i].train = false;
	    FEP[i].exp = NULL;
	}

	//get the file for this set of emission prob

	// executable_path is from user
	if(!EElement->Attribute("file"))
	{
	    cout<<"Error: set_FEP_from_file : "
		<<"file attribute in Emission_Probs tag is NULL"<<endl;
	    check++;
	}else{
	    strcpy(filename,In_dir);
	    strcat(filename,EElement->Attribute("file"));
	}

	//LElement=ENode->FirstChildElement("label");
	LCount = 0;
    }

    LCount = 0;

    if(!check)
    {
	//open file for setting the emission probs

	instream.open(filename);
	if(!instream)
	{
	    cout<<"Error: set_emission_prob_from_file, can not read file : "<<filename<<endl;
	    check++;
	}
	
	while(!instream.eof())
	{    
	    
	    // check consistency of EP[ProbCount].id

	    instream.getline(line,Max_line_length-1);
	    //instream>>buffer;
	    check+= splitstring(line,NoOfItems,&Items,' ');	    

	    if(NoOfItems==2) // no name and no train
	    {
		if(atoi(Items[1])==0)
		{
		    cout<<"Error:: set_FEP_from_xml : "
			<<"FEP["<<LCount<<"].dim : "<<Items[1]
			<<" should be an integer >0."<<endl;
		    check++;
		    break;
		}
		FEP[LCount].name=Nullstrcpy("Not Defined",check);
		if(check)
		{
		    cout<<"Error:: set_FEB_from_xml : Nullstrcpy, "
			<<"FEP["<<LCount<<"].name : Not Defined"<<endl;
		    break;
		}		
		FEP[LCount].dim = atoi(Items[1]);
		FEP[LCount].train = false;

	    }else if(NoOfItems == 3) // with name or train
	    {
		if(!strcmp(Items[2],"train"))
		{
		    if(atoi(Items[1])==0)
		    {
			cout<<"Error:: set_FEP_from_xml : "
			    <<"FEP["<<LCount<<"].dim : "<<Items[1]
			    <<" should be an integer >0."<<endl;
			check++;
			break;
		    }
		    FEP[LCount].dim = atoi(Items[1]);		    
		    FEP[LCount].train = true;
		}else{
		    if(atoi(Items[2])==0)
		    {
			cout<<"Error:: set_FEP_from_xml : "
			    <<"FEP["<<LCount<<"].dim : "<<Items[2]
			    <<" should be an integer >0."<<endl;
			check++;
			break;
		    }
		    FEP[LCount].name=Nullstrcpy(Items[1],check);
		    if(check)
		    {
			cout<<"Error:: set_FEP_from_xml : Nullstrcpy,"
			    <<"FEP["<<LCount<<"].name : "<<Items[1]<<endl;
			break;
		    }
		    FEP[LCount].dim = atoi(Items[2]);
		}
	    }else if(NoOfItems==4){
		if(atoi(Items[2])==0)
		{
		    cout<<"Error:: set_FEP_from_xml : "
			<<"FEP["<<LCount<<"].dim : "<<Items[2]
			<<" should be an integer >0."<<endl;
		    check++;
		    break;
		}
		FEP[LCount].name=Nullstrcpy(Items[1],check);
		if(check)
		{
		    cout<<"Error:: set_FEP_from_xml : Nullstrcpy,"
			<<"FEP["<<LCount<<"].name : "<<Items[1]<<endl;
		    break;
		}
		FEP[LCount].dim = atoi(Items[2]);
		if(!strcmp(Items[3],"train"))
		{
		    FEP[LCount].train = true;
		}else{
		    cout<<"Error:: set_FEP_from_xml : "
			<<"FEP["<<LCount<<"] : "<<Items[3]
			<<" : unrecognized label."<<endl;
		    check++;
		    break;
		}
	    }else{
		cout<<"Error:: set_FEP_from_xml : line : "
		    <<line<<" not in specified format"<<endl;
		check++;
	    }	   
	    
	    check+=break_id_into_tokens(Items[0],&tokens,NumOfTokens);
	    
	    if((NumOfTokens>2)||(NumOfTokens<1))
	    {
		check++;
	    }
	    if(check)
	    {
		cout<<"Error: set_FEP_from_xml, break_id_into_tokens" <<endl;
		break;
	    }
		
	    if(strcmp(tokens[0],FEP_tname))
	    {
		cout<<"Error: set_FEP_from_xml, id in .txt file : "<<tokens[0]
		    <<" does not match with id in .xml file : "<<FEP_tname<<endl;
		check++;
		break;
	    }
	    
	    if(atoi(tokens[1])!=LCount)
	    {
		cout<<"Error: set_FEP_from_xml:: emission prob id  from file "<<tokens[1]
		    <<" not match with LCount : "<<LCount<<endl;
		check++;
		break;
	    }
	    
	    // Get Emission Prob

	    array<int> TmpIndex(1);
	    TmpIndex.SetDimension(0,FEP[LCount].dim);
	    FEP[LCount].prob.SetNumberofDimensions(FEP[LCount].dim);
	    FEP[LCount].pseudoprob.SetNumberofDimensions(FEP[LCount].dim);
	    
	    size = MP->get_Alphabet_size();
	    for(int i =0; i<FEP[LCount].dim; i++){
		FEP[LCount].prob.SetDimension(i,size,0);
		FEP[LCount].pseudoprob.SetDimension(i,size,0);
	    }
	    Tprob = 0; // calculate the sum of the prob declare in the file, should be sum up to 1

	    while(1)
	    {		
		if(instream.eof()){
		    /*
		    cout<<"Error:: set_emission_prob_from_file : input file : "
			<<filename<<" invalid "<<endl;
		    
		    check++;
		    */
		    break;
		}
		
		//instream>>buffer;		
		//cout<<buffer<<endl;

		instream.getline(line,Max_line_length-1);
		
		check+= splitstring(line,NoOfItems,&Items,' ');

		if(NoOfItems==0)
		{
		//if((!strcmp(Items[0],"::"))||(!strcmp(Items[0],"--")))
		    break;
		}else{
		    if((NoOfItems<2)||(NoOfItems>3))
		    {
			cout<<"ERROR: set_FEP_from_file:: invalid line format:" <<endl
			    <<line<<endl;
		    }else
		    {
			// check consistency of EP[ProbCount].dim
			if(FEP[LCount].dim!=strlen(Items[0])){
			    cout<<"Error: set_FEP_from_file:: dimension of buffer("<<Items[0]<<") : "
				<<strlen(Items[0])<<" not match with FEP["<<LCount<<"].dim : "
				<<FEP[LCount].dim<<endl;
			    check++;
			    break;
			}
			//instream>>prob;
			prob = atof(Items[1]);
			if(NoOfItems==3)
			{
			    pseudoprob = atof(Items[2]);
			}else{
			    pseudoprob = 0;
			}
		    }
		}
		
		if((prob==0)&&(pseudoprob==0))
		{
		    continue;
		}
		
		for(int j=0; j<FEP[LCount].dim; j++)
		{
		    char alph = Items[0][j];
		    alphint = convert_alphabet_to_int(MP,alph);
		    TmpIndex.SetElement(j,alphint);
		}
		FEP[LCount].prob.SetElement(TmpIndex,static_cast<Prob>(prob));
		if(pseudoprob!=0)
		{
		    FEP[LCount].pseudoprob.SetElement(TmpIndex,pseudoprob);
		}
		Tprob+=prob;
	    }
	    LCount++;
	    if(abs(Tprob-1.0)>Max_deviation){
		cout<<"Error: set_FEP_from_file,prob sum up not = 1"<<endl;
		cout<<"Emission_prob["<<LCount<<"]"<<endl;
		check++;
		break;
	    }
	    NoOfItems = 0;
	}
	instream.close();
    }
  
    if(line) delete[] line;
    line = NULL;

    if(filename) delete[] filename;
    filename = NULL;
    Node = 0;
    ENode = 0;
    EElement = 0;
    LElement = 0;
    if(tokens)
    {
	for(int i = 0; i<NumOfTokens; i++){
	    if(tokens[i]) delete [] tokens[i];
	    tokens[i] = NULL;
	}
	delete [] tokens;
    }
    tokens = NULL;

    if(Items){
	for(int i=0; i<Max_number_of_items; i++)
	{
	    if(Items[i]) delete[] Items[i];
	    Items[i] = NULL;
	}
	delete [] Items;
    }
    Items = NULL;
    NoOfItems = 0;
    
    if(check)
    {
	if(FEP)
	{
	    for(int i =0; i<FEPsize; i++)
	    {
		if(i<=LCount)
		{
		    if(FEP[i].name) delete[] FEP[i].name;
		    FEP[i].name = NULL;
		    if(FEP[i].exp) delete[] FEP[i].exp;
		    FEP[i].exp = NULL;
		}
	    }
	    delete[] FEP;
	}
	FEP = NULL;
	if(FEP_tname) delete [] FEP_tname;
	FEP_tname = NULL;
    }
    return(check);
}

int EmissionProb::get_FEP_parameters_for_training(const char* filename,
						  int& NumOfAutoEP,
						  int** AutoEP)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: get_FTP_parameters_for_training : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * SequenceNode = 0;
    TiXmlNode * TmpNode = 0; 
    TiXmlElement* Element = 0;
    int count = 0;
    int index = -1;
    bool trainall = false;
    NumOfAutoEP = 0;
    const int MaxListLength = 200;

    if(*AutoEP) delete (*AutoEP);
    (*AutoEP) = NULL;
    (*AutoEP) = new int [MaxListLength];
    for(int i =0; i<MaxListLength; i++)
    {
	(*AutoEP)[i] = -1;
    }
  
    // get train and exp for FTP
    SequenceNode = doc.FirstChild("HMMConverter");
    if(!SequenceNode)
    {
        check++;
        return check;
    }

    SequenceNode = SequenceNode->FirstChild("sequence_analysis");
    if(!SequenceNode)
    {
	return check;
    }
    SequenceNode = SequenceNode->FirstChild("parameter_training");
    if(!SequenceNode)
    {
	return check;
    }
    SequenceNode = SequenceNode->FirstChild("train_free_parameters");
    if(!SequenceNode)
    {	    
	NumOfAutoEP = FEPsize;
	for(int i=0; i<FEPsize; i++)
	{
	    //FEP[i].train = true;
	    if(FEP[i].exp) delete[] FEP[i].exp; // get from state x+1;
	    FEP[i].exp = NULL;
	    (*AutoEP)[i] = i;
	}	
	return check;
    }
    
    TmpNode = SequenceNode->FirstChild("FreeEmissionParameters");
    Element = SequenceNode->FirstChildElement("FreeEmissionParameters");
    if(TmpNode)
    {
	Element = TmpNode->FirstChildElement(FEP_tname);
	if(!Element)
	{
	    cout<<"ERROR: EmissionProb:: get_FEP_parameters_for_training: No <"
		<<FEP_tname<<"> tag under <FreeEmissionParameters> tag."<<endl;
	    check++;
	    return check;
	}
	count = 0;
	
	while((Element)&&(!check))
	{
	    // get idref and set the train field
	    if(!Element->Attribute("idref"))
	    {
		cout<<"ERROR: EmissionProb:: get_FEP_parameters_for_training: No "
		    <<"idref attribute in the "<<count<<"-th <"<<FEP_tname
		    <<"> tag!"<<endl;
		check++;
	    }
	    
	    // set train
	    if(!check)
	    {
		index = convert_FEP_id_to_int(this,Element->Attribute("idref"));
		if((index<0)||(index>=FEPsize))
		{
		    cout<<"ERROR: EmissionProb:: get_FEP_parameters_for_training: "
			<<"id: "<<Element->Attribute("idref")
			<<"out of range[0.."<<FEPsize-1<<"] in the "
			<<count<<"-th tag!"<<endl;
		    check++;
		}else{
		    FEP[index].train = true;
		}		    
	    }
	    // set exp
	    if(Element->Attribute("exp"))
	    {
		if(FEP[index].exp) delete[] FEP[index].exp;
		FEP[index].exp = NULL;
		FEP[index].exp = Nullstrcpy(Element->Attribute("exp"),check);
		if(check)
		{
		    cout<<"ERROR: EmissionProb:: get_FEP_parameters_for_training: "
			<<" error in Nullstrcpy for "
			<<Element->Attribute("exp")<<" in "
			<<count<<"-th "<<"<"<<FEP_tname<<"> tag!"<<endl;
		    check++;
		}
		
	    }else{
		// auto detect
		(*AutoEP)[NumOfAutoEP] = index;
		NumOfAutoEP++;
		if(FEP[index].exp) delete[] FEP[index].exp;
		FEP[index].exp = NULL;
		
	    }
	    
	    Element = Element->NextSiblingElement();
	    count++;
	}
	    	
    }else{
	// no need to train
	/*
	for(int i=0; i<FEPsize; i++)
	{
	    FEP[i].train = false;
	}
	*/
    }

    SequenceNode = 0;
    TmpNode = 0; 
    Element = 0;
    return check;
}

int EmissionProb::get_GEP_parameters_for_training(const char* filename,
						  model_parameters* const MP,
						  int NumOfAutoEP)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: get_GEP_parameters_for_training : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * ModelNode = 0;
    TiXmlNode * SequenceNode = 0;
    TiXmlNode * TmpNode = 0;
    TiXmlNode * Node = 0;
    TiXmlNode * fromNode = 0;
    TiXmlElement* Element = 0;
    TiXmlElement* fromElement = 0;
    int count = 0;
    int emitcount = 0;
    int index = -1;
    int fromindex = -1;
    
    SequenceNode = doc.FirstChild("HMMConverter");
    if(!SequenceNode)
    {
        check++;
        return check;
    }

    SequenceNode = SequenceNode->FirstChild("sequence_analysis");
    if(!SequenceNode)
    {
	return check;
    }
    SequenceNode = SequenceNode->FirstChild("parameter_training");
    if(!SequenceNode)
    {
	return check;
    }
    SequenceNode = SequenceNode->FirstChild("train_free_parameters");
    if(!SequenceNode)
    {
	return check;
    }
    TmpNode = SequenceNode->FirstChild("GroupEmissions");
    if(!TmpNode)
    {
	return check;
    }
    
    Element = SequenceNode->FirstChildElement("GroupEmissions");
    
    if(!Element->Attribute("id"))
    {
	cout<<"ERROR: EmissionProb:: get_GEP_parameters_for_training: No id attribute in "
	    <<"the <GroupEmissions> tag!"<<endl;
	check++;
    }else{
	if(GEP_tname) delete[] GEP_tname;
	GEP_tname = NULL;
	GEP_tname = Nullstrcpy(Element->Attribute("id"),check);
	if(check)
	{
	    cout<<"ERROR: EmissionProb:: get_GEP_parameters_for_training: "
		<<" error in Nullstrcpy for : "<<Element->Attribute("id")
		<<" in <GroupEmissions> tag."<<endl;			
	}
    }

    if(!check)
    {

	Element = TmpNode->FirstChildElement(GEP_tname);
	if(!Element)
	{
	    cout<<"ERROR: EmissionProb:: get_GRP_parameters_for_training: No <"
		<<GEP_tname<<"> tag under <GroupEmissions> tag."<<endl;
	    check++;
	}

	count = 0;
	// count the number of GEP
	while((Element)&&(!check))
	{
	    Element = Element->NextSiblingElement();
	    count++;		
	}
	GEPsize = count+NumOfAutoEP;
	GEP = new GEProb[GEPsize+NumOfAutoEP];
	Element = TmpNode->FirstChildElement(GEP_tname);
	Node = TmpNode->FirstChild(GEP_tname);
	count = 0;
	while((Element)&&(!check))
	{
	    if(count>=GEPsize)
	    {
		cout<<"ERROR: EmissionProb:: get_GEP_parameters_for_training: "
		    <<"number of <"<<GEP_tname<<"> tags count("<<count<<") out of range[0.."
		    <<GEPsize-1<<"]"<<endl;
		check++;
		break;
	    }
	    // get idref and set the train field
	    if(!Element->Attribute("id"))
	    {
		cout<<"ERROR: EmissionProb:: get_GEP_parameters_for_training: No id "
		    <<"attribute in the "<<count<<"-th <"<<GEP_tname<<"> tag!"<<endl;
		check++;
	    }
	    if(!check)
	    {
		index = convert_GEP_id_to_int(this,Element->Attribute("id"));
		if((index<0)||(index>=GEPsize))
		{
		    cout<<"ERROR: EmissionProb:: get_GEP_parameters_for_training: "
			<<"id: "<<Element->Attribute("id")<<"out of range[0.."
			<<GEPsize-1<<"] in the "<<count<<"-th tag!"<<endl;
		    check++;
		}
	    }
	    // get from
	    
	    // get NumOffrom
	    emitcount = 0;
	    fromNode = Node->FirstChild("from");
	    fromElement = Node->FirstChildElement("from");
	    while((fromNode)&&(!check))
	    {
		// add count from the number of transition from corresponding state

		fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));

		if((fromindex<0)||(fromindex>=MP->get_Number_of_States()))
		{
		    cout<<"ERROR:: get_GEP_parameters_for_training : "
			<<" idref attribute at the "<<count<<"-th "<<GEP_tname
			<<" tag("<<fromElement->Attribute("idref")<<"out of "
			<<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
		    check++;
		}else{
		    emitcount++;
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
		
	    }

	    GEP[index].NumOffrom = emitcount;

	    GEP[index].from = new int[GEP[index].NumOffrom];
	  
	    // set the from and to arrays 
	    emitcount = 0;
	    fromNode = Node->FirstChild("from");
	    fromElement = Node->FirstChildElement("from");
	    while((fromNode)&&(!check))
	    {

		fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));

		if((fromindex<0)||(fromindex>=MP->get_Number_of_States()))
		{
		    cout<<"ERROR:: get_GEP_parameters_for_training : "
			<<" idref attribute at the "<<count<<"-th "<<GEP_tname
			<<" tag("<<fromElement->Attribute("idref")<<"out of "
			<<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
		    check++;
		}else{
		    GEP[index].from[emitcount] = fromindex;
		    emitcount++;
		}
			
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
		
	    }
	    Element = Element->NextSiblingElement();
	    count++;
	}
    }
    if(check)
    {
	if(GEP)
	{
	    for(int i=0; i<GEPsize; i++)
	    {
		if(GEP[i].from) delete[] GEP[i].from;
		GEP[i].from = NULL;
	    }
	    delete [] GEP;
	}	
	GEP = NULL;
    }
    return check;    
}

int EmissionProb::auto_get_parameters_for_training(Hmm* const p,
						   const int NumOfAutoEP,
						   int* AutoEP)
{
    int check = 0;
    int start_index = -1;
    int count = 0;
    int NumOfStates = p->get_number_of_states();

    if(NumOfAutoEP==0)
    {
	return check;
    }
    char* buffer = new char[Max_word_length];
    char* buffer2 = new char[Max_word_length];
    if(GEPsize==0)
    {
	if(GEP) delete[] GEP;
	GEP = NULL;
	GEP = new GEProb[NumOfAutoEP];
	GEPsize = NumOfAutoEP;
	if(GEP_tname)  delete[] GEP_tname;
	GEP_tname = NULL;
	GEP_tname = Nullstrcpy("GEP",check);
	start_index = 0;
    }else{
	start_index = GEPsize - NumOfAutoEP;
    }

    for(int i = start_index; i<GEPsize; i++)
    {
	// set FEP[AutoEP[i-start_index]].exp
	int curindex = i-start_index;

	if(FEP[AutoEP[curindex]].exp) delete[] FEP[AutoEP[curindex]].exp;
	FEP[AutoEP[curindex]].exp = NULL;
	
	sprintf(buffer,"%d",i);
	FEP[AutoEP[curindex]].exp = Nullstrcpy(GEP_tname,check);
	
	if(check)
	{
	    cout<<"ERROR: EmissionProb:: auto_get_parameters_for_training : "
		<<"error in Nullstrcpy of FEP["<<AutoEP[i-start_index]<<"].exp("
		<<strcat(GEP_tname,buffer)<<")" <<endl;
	    break;
	}
	strcat(FEP[AutoEP[curindex]].exp,".");
	strcat(FEP[AutoEP[curindex]].exp,buffer);

	// get NumOffrom of GEP[i];
	count = 0;
	sprintf(buffer,"%d",AutoEP[curindex]);
	strcpy(buffer2,FEP_tname);	
	strcat(buffer2,".");
	strcat(buffer2,buffer);

	
	for(int j =0; j<NumOfStates; j++)
	{
	    if((*p)[j]->get_emission_probs_expression())
	    {
		if(!cmp_nocase((*p)[j]->get_emission_probs_expression(),buffer2))
		{
		    count++;
		}
	    }
	}

	GEP[i].NumOffrom = count;
	if(count==0)
	{
	    cout<<"ERROR: EmissionProb:: auto_get_parameters_for_training "
		<<"GEP["<<i<<"].NumOffrom is 0!"<<endl;
	    check++;
	    break;
	}else{
	    GEP[i].from = new int[GEP[i].NumOffrom];
	}
	
	// get from 

	count = 0;
	for(int j =0; j<NumOfStates; j++)
	{
	    if((*p)[j]->get_emission_probs_expression())
	    {
		if(!cmp_nocase((*p)[j]->get_emission_probs_expression(),buffer2))
		{
		    GEP[i].from[count] = j;
		    count++;
		}
	    }
	}

    }
    if(buffer) delete[] buffer;
    buffer = NULL;
    if(buffer2) delete [] buffer2;
    buffer2 = NULL;
    return check;
}

int EmissionProb::check_consistency_of_GEP(Hmm* const p)
{
    int check = 0;
    int index = -1;
    for(int i =0; i<FEPsize; i++)
    {
	if(FEP[i].train)
	{
	    index = convert_GEP_id_to_int(this,FEP[i].exp);
	    if(index<0)
	    {
		cout<<"ERROR: EmissionProb:: check_consistency_of_GEP "
		    <<"input GEP ";
		if(!FEP[i].exp)
		{
		    cout<<" is NULL!. "<<endl;
		}else{
		    cout<<" : "<<FEP[i].exp<<" is invalid."<<endl;
		}
		check++;
		break;
	    }
	    
	    for(int j=0; j<GEP[index].NumOffrom; j++)
	    {
		// check dimension
		
		if(FEP[i].dim != (*p)[GEP[index].from[j]]->get_number_of_dimensions_of_emission_probs())
		{
		    cout<<"ERROR: EmissionProb:: check_consistency_of_GEP "
			<<"dim of FEP["<<i<<"]("<<FEP[i].dim<<") ! = "
			<<"dim of State["<<GEP[index].from[j]<<"]("
			<<(*p)[GEP[index].from[j]]->get_number_of_dimensions_of_emission_probs()
			<<")."<<endl;
		    check++;
		}
	    }
	}
    }
    return check;
}

EmissionProb :: EmissionProb()
{
    FEPsize = 0;
    FEP_tname = 0;
    FEP = NULL;
    GEPsize = 0;
    GEP_tname = NULL;
    GEP = NULL;
}

EmissionProb:: EmissionProb(const char* In_dir,
			    const char* xmlfile, 
			    model_parameters* MP, 
			    int& check){
    if(check)
    {
	cout<<"ERROR:: EmissionProb constructor : check("<<check<<") > 0"<<endl;
    }else
    {
	GEPsize = 0;
	GEP_tname = NULL;
	GEP = NULL;
	check+=set_FEP_from_file(In_dir,xmlfile,MP);
    }
    if(check){
	cout<<"Error:: EmissionProb constructor : set_emission_prob_from_file "<<endl;
	this->~EmissionProb();
    }
}
    
EmissionProb::~EmissionProb()
{    
    if(FEP)
    {
	for(int i=0; i<FEPsize; i++)
	{
	    if(FEP[i].name) delete[] FEP[i].name;
	    FEP[i].name = NULL;
	    if(FEP[i].exp) delete[] FEP[i].exp;
	    FEP[i].exp = NULL;
	}
	delete [] FEP;
    }
    FEP = NULL;
    if(FEP_tname) delete[] FEP_tname;
    FEP_tname = NULL;
    FEPsize = 0;

    if(GEP)
    {
	for(int i =0; i<GEPsize; i++)
	{
	    if(GEP[i].from) delete[] GEP[i].from;
	    GEP[i].from = NULL;
	}
	delete [] GEP;
    }
    GEP = NULL;
    if(GEP_tname) delete[] GEP_tname;
    GEP_tname = NULL;
    GEPsize = 0;

}

//accessors
 
int EmissionProb::get_FEPsize() const
{
    return FEPsize;
}


char* EmissionProb::get_FEP_tname() const
{
    return FEP_tname;
}

char* EmissionProb::get_FEP_name(const int i) const 
{
    if((i<0)||(i>=FEPsize))
	return NULL;
    return FEP[i].name;
}

int EmissionProb::get_FEP_dim(const int i) const
{
    if((i<0)||(i>=FEPsize))
	return -1;
    return FEP[i].dim;
}

Prob EmissionProb:: get_FEP_prob(const int i, const array<int> &indices) const 
{
    if((i<0)||(i>=FEPsize))
	return -1;
    return FEP[i].prob.GetElement(indices);
}

Prob EmissionProb:: get_FEP_prob(const int i, const int linear_index) const 
{
    if((i<0)||(i>=FEPsize))
	return -1;
    return FEP[i].prob.GetElement(linear_index);
}

array<Prob> EmissionProb:: get_FEP_prob(const int i) const 
{
    if((i<0)||(i>=FEPsize)){
	array<Prob> temp(0);
	return temp;
    }
    return FEP[i].prob;
}

Prob EmissionProb:: get_FEP_pseudoprob(const int i, const array<int> &indices) const 
{
    if((i<0)||(i>=FEPsize))
	return -1;
    return FEP[i].pseudoprob.GetElement(indices);
}

Prob EmissionProb:: get_FEP_pseudoprob(const int i, const int linear_index) const 
{
    if((i<0)||(i>=FEPsize))
	return -1;
    return FEP[i].pseudoprob.GetElement(linear_index);
}

array<Prob> EmissionProb:: get_FEP_pseudoprob(const int i) const 
{
    if((i<0)||(i>=FEPsize)){
	array<Prob> temp(0);
	return temp;
    }
    return FEP[i].pseudoprob;
}

FEProb EmissionProb:: get_FEP(const int i) const
{
    return FEP[i];
}

bool EmissionProb:: is_FEP_train(const int i) const
{
    if((i<0)||(i>=FEPsize))
    {
	cout<<"EmissionProb:: is_FEP_train: input index("<<i
	    <<") out of range[0.."<<FEPsize-1<<"]."<<endl;
	return false;
    }
    return FEP[i].train;
}

char* EmissionProb:: get_FEP_exp(const int i) const
{
    if((i<0)||(i>=FEPsize))
    {
	cout<<"EmissionPRob:: get_FEP_exp: input index("<<i
	    <<") out of range[0.."<<FEPsize-1<<"]."<<endl;
	return NULL;
    }
    return FEP[i].exp;
}

// for GEP

int EmissionProb:: get_GEPsize() const
{
    return GEPsize;
}

char* EmissionProb:: get_GEP_tname() const
{
    return GEP_tname;
}

int EmissionProb:: get_GEP_NumOffrom(const int i) const
{
    if((i<0)||(i>=GEPsize))
    {
	cout<<"EmissionProb:: get_GEP_NumOffrom: input index("<<i
	    <<") out of range[0.."<<GEPsize-1<<"]."<<endl;
	return -1;
    }
    return GEP[i].NumOffrom;
}

int EmissionProb::get_GEP_from(const int i, const int j) const
{
    if((i<0)||(i>=GEPsize))
    {
	cout<<"EmissionProb:: get_GEP_from: input index("<<i
	    <<") out of range[0.."<<GEPsize-1<<"]."<<endl;
	return -1;
    }
    if((j<0)||(j>=GEP[i].NumOffrom))
    {
	cout<<"EmissionProb:: get_GEP_from: input index("<<j
	    <<") out of range[0.."<<GEP[i].NumOffrom-1<<"]."<<endl;
	return -1;
    }
    return GEP[i].from[j];
}

// mutator

int EmissionProb::set_FEPsize(const int s)
{
    if(s<0){
	cout<<"Error! set_FEPSize:: "
	    <<"s "<<s<<" out of range! "<<endl;
	return 1;
    }
    FEPsize=s;
    return 0;
}

int EmissionProb::set_FEP_tname(char* const tn)
{
    int check = 0;

    if(FEP_tname) delete[] FEP_tname;
    FEP_tname=NULL;
    FEP_tname = Nullstrcpy(tn,check);
    if(check){
	cout<<"Error! set_FEP_tname:: "<<endl;
	cout<<"tn : "<<tn<<endl;
    }
    return check;
}

int EmissionProb::set_FEP_name(const int i, char* const n)
{
    int check = 0;
    
    if((i>=FEPsize)||(i<0)){
	cout<<"Error! set_FEP_name:: "<<endl;
	cout<<"i "<<i<<" exceed FEPsize "<<FEPsize<<endl;
	check++;
	return check;
    }

    if(FEP[i].name) delete[] FEP[i].name;
    FEP[i].name=NULL;
    FEP[i].name = Nullstrcpy(n,check);
    if(check){
	cout<<"Error! set_Emission_Prob_name:: "<<endl;
	cout<<"Input FEP["<<i<<"].name : "<<n<<endl;
    }
    return check;
}

int EmissionProb::set_FEP_dim(const int i, const int d)
{
    int check = 0;
    
    if((i>=FEPsize)||(i<0)){
	cout<<"Error! set_FEP_dim:: "<<endl;
	cout<<"i "<<i<<" exceed FEPsize "<<FEPsize<<endl;
	check++;
    }

    if(check!=0){
	cout<<"Error! set_Emission_Prob_dim:: "<<endl;
	cout<<"Input FEP["<<i<<"].dim : "<<d<<endl;
    }else{
	FEP[i].dim = d;
    }
    return check;
}

int EmissionProb::set_FEP_prob(const int i, const array<int> &indices, const Prob p)
{
    
    if((i>=FEPsize)||(i<0)){
	cout<<"Error! set_FEP_prob:: "<<endl;
	cout<<"i "<<i<<" exceed FEPsize "<<FEPsize<<endl;
	return 1;
    }

    FEP[i].prob.SetElement(indices,p);
    return 0;
}

int EmissionProb::set_FEP_prob(const int i, int linear_index, const Prob p)
{

    if((i>=FEPsize)||(i<0)){
	cout<<"Error:: set_FEP_prob: "<<endl;
	cout<<"i "<<i<<" exceed FEPsize "<<FEPsize<<endl;
	return 1;
    }

    FEP[i].prob.SetElement(linear_index,p);
    return 0;
}

int EmissionProb::set_FEP_pseudoprob(const int i, const array<int> &indices, const Prob p)
{
    
    if((i>=FEPsize)||(i<0)){
	cout<<"Error! set_FEP_pseudocout:: "<<endl;
	cout<<"i "<<i<<" exceed FEPsize "<<FEPsize<<endl;
	return 1;
    }

    FEP[i].pseudoprob.SetElement(indices,p);
    return 0;
}

int EmissionProb::set_FEP_pseudoprob(const int i, int linear_index, const Prob p)
{

    if((i>=FEPsize)||(i<0)){
	cout<<"Error:: set_FEP_pseudoprob: "<<endl;
	cout<<"i "<<i<<" exceed FEPsize "<<FEPsize<<endl;
	return 1;
    }

    FEP[i].pseudoprob.SetElement(linear_index,p);
    return 0;
}

int EmissionProb::set_FEP_train(const int i, const bool t)
{
    if((i<0)||(i>=FEPsize))
    {
	cout<<"Error:: set_FEP_train: input index("<<i
	    <<") out of range[0.."<<FEPsize-1<<"]."<<endl;
	return 1;
    }

    FEP[i].train = t;
    return 0;
}

int EmissionProb:: set_FEP_exp(const int i, char* const e)
{
    int check = 0;
    if((i<0)||(i>=FEPsize))
    {
	cout<<"Error:: set_FEP_exp: input index("<<i
	    <<") out of range[0.."<<FEPsize-1<<"]."<<endl;
	check++;
    }

    if(!check)
    {
	if(FEP[i].exp) delete[] FEP[i].exp;
	FEP[i].exp=NULL;
	FEP[i].exp = Nullstrcpy(e,check);
	if(check){
	    cout<<"Error! set_FEP_tname:: "<<endl;
	    cout<<"e : "<<e<<endl;
	}
    }
    
    return check;
}

int EmissionProb::set_GEPsize(const int s)
{
    if(s<0){
	cout<<"Error! set_GEPSize:: "
	    <<"s "<<s<<" out of range! "<<endl;
	return 1;
    }
    GEPsize=s;
    return 0;
}

int EmissionProb::set_GEP_tname(char* const tn)
{
    int check = 0;

    if(GEP_tname) delete[] GEP_tname;
    GEP_tname=NULL;
    GEP_tname = Nullstrcpy(tn,check);
    if(check){
	cout<<"Error! set_GEP_tname:: "<<endl;
	cout<<"tn : "<<tn<<endl;
    }
    return check;
}

int EmissionProb::set_GEP_NumOffrom(const int i, const int num)
{
    if((i<0)||(i>=GEPsize))
    {
	cout<<"Error:: set_GEP_NumOffrom: input index("<<i
	    <<") out of range[0.."<<GEPsize-1<<"]."<<endl;
	return 1;	
    }
    GEP[i].NumOffrom = num;
    return 0;
}

int EmissionProb::set_GEP_from(const int i, const int j, const int f)
{
    if((i<0)||(i>=GEPsize))
    {
	cout<<"Error:: set_GRP_from: input index("<<i
	    <<") out of range[0.."<<GEPsize-1<<"]."<<endl;
	return 1;	
    }
    if((j<0)||(j>=GEP[i].NumOffrom))
    {
	cout<<"Error:: set_GEP_from: input index("<<j
	    <<") out of range[0.."<<GEP[i].NumOffrom-1<<"]."<<endl;
	return 1;	
    }
    GEP[i].from[j] = f;
    return 0;
}

// operators
EmissionProb & EmissionProb::operator = (const EmissionProb &ep)
{
    int check = 0;
    if(this != &ep)
    {
	if(ep.FEP)
	{
	    if(FEP)
	    {
		for(int i =0; i<FEPsize;i++)
		{
		    if(FEP[i].name) delete[] FEP[i].name;
		    FEP[i].name = NULL;
		    if(FEP[i].exp) delete[] FEP[i].exp;
		    FEP[i].exp = NULL;
		}
		delete[] FEP;
	    }
	    FEP = NULL;
	    FEPsize = ep.FEPsize;
	    if(FEP_tname) delete[] FEP_tname;
	    FEP_tname = NULL;
	    FEP_tname = Nullstrcpy(ep.FEP_tname,check);
	    FEP = new FEProb[FEPsize];
	    for(int i =0; i<FEPsize; i++)
	    {
		FEP[i].name = Nullstrcpy(ep.FEP[i].name,check);
		FEP[i].dim = ep.FEP[i].dim;
		FEP[i].train = ep.FEP[i].train;
		FEP[i].prob = ep.FEP[i].prob;
		FEP[i].exp = Nullstrcpy(ep.FEP[i].exp,check);
	    }
	}else{
	    FEPsize = 0;
	    if(FEP_tname) delete[] FEP_tname;
	    FEP_tname = NULL;
	    if(FEP)
	    {
		for(int i =0; i<FEPsize;i++)
		{
		    if(FEP[i].name) delete[] FEP[i].name;
		    FEP[i].name = NULL;
		    if(FEP[i].exp) delete[] FEP[i].exp;
		    FEP[i].exp = NULL;
		}
		delete[] FEP;
	    }
	    FEP = NULL;
	}
	
	if(ep.GEP)
	{
	    if(GEP)
	    {
		for(int i =0; i<GEPsize;i++)
		{
		    if(GEP[i].from) delete[] GEP[i].from;
		    GEP[i].from = NULL;
		}
		delete[] GEP;
	    }
	    GEP = NULL;
	    GEPsize = ep.GEPsize;
	    if(GEP_tname) delete[] GEP_tname;
	    GEP_tname = NULL;
	    GEP_tname = Nullstrcpy(ep.GEP_tname,check);
	    GEP = new GEProb[GEPsize];
	    for(int i =0; i<GEPsize; i++)
	    {
		GEP[i].NumOffrom = ep.GEP[i].NumOffrom;
		for(int j=0; j<GEP[i].NumOffrom; j++)
		{
		    GEP[i].from[j] = ep.GEP[i].from[j];
		}
	    }
	}else{
	    GEPsize = 0;
	    if(GEP_tname) delete[] GEP_tname;
	    GEP_tname = NULL;
	    if(GEP)
	    {
		for(int i =0; i<GEPsize;i++)
		{
		    if(GEP[i].from) delete[] GEP[i].from;
		    GEP[i].from = NULL;
		}
		delete[] GEP;
	    }
	    GEP = NULL;
	}	
    }
    return (*this);
}

//Other functions

int EmissionProb::get_parameters_for_training(const char* filename,
					      model_parameters* const MP,
					      Hmm* const p)
{
    int check = 0;
    if(!filename)
    {
	cout<<"ERROR: EmissionProb:: get_parameters_for_training: "
	    <<"input filename is NULL"<<endl;
	check++;
	return check;
    }

    int NumOfAutoEP = 0;
    int* AutoEP = NULL;
    check += get_FEP_parameters_for_training(filename,NumOfAutoEP,&AutoEP);
    check += get_GEP_parameters_for_training(filename,MP,NumOfAutoEP);
    check += auto_get_parameters_for_training(p,NumOfAutoEP,AutoEP);

    check += check_consistency_of_GEP(p);

  
    if(AutoEP) delete [] AutoEP;
    AutoEP = NULL;

    return check;
}

void EmissionProb::print(std:: ostream &o) const
{
    o<<"Emission Probs"<<endl;
    o<<"------------------------------------------------------"<<endl<<endl;
    o<<"Free_Emission_Prob_Name : "<<FEP_tname<<endl;
    o<<"-----------------------------------------------------"<<endl;
    for(int i=0;i<FEPsize;i++)
    {
	o<<"FEP["<<i<<"].name : "<<FEP[i].name<<endl;
	o<<"FEP["<<i<<"].dim : "<<FEP[i].dim<<endl;
	o<<"FEP["<<i<<"].train : "<<FEP[i].train<<endl;
	if(FEP[i].exp)
	{
	    o<<"FEP["<<i<<"].exp : "<<FEP[i].exp<<endl;
	}
	o<<"dimensions(FEP["<<i<<"].prob) : ";
	FEP[i].prob.PrintDimensions(o);
	o<<"FEP["<<i<<"].prob : "<<endl;
	FEP[i].prob.PrintonlyNonZerowithIndices(o);
    }
    o<<"---------------------------------------------------"<<endl;
    o<<"Free Emission Parameters to be trained "<<endl;
    o<<"Number of Group_Emission_Prob : "<<GEPsize<<endl;
    for(int i = 0; i<GEPsize; i++)
    {
	o<<"Group("<<i<<") : "<<endl;
	o<<"from: ";
	for(int j=0; j<GEP[i].NumOffrom; j++)
	{
	    o<<"S."<<GEP[i].from[j];
	    if(j!=GEP[i].NumOffrom-1)
	    {
		o<<" + ";
	    }
	}   
	o<<endl;
    }
    o<<endl;
    return;
}

int convert_FEP_id_to_int(EmissionProb* const EP, const char* buffer)
{
    char** tokens;
    int num_of_tokens;
    int check = 0;
    check+=break_id_into_tokens(buffer,&tokens,num_of_tokens);
    if(check)
    {
	//cout<<"ERROR: convert_FEP_id_to_number : buffer : "<<buffer<<endl;
	return -1;
    }
    if((num_of_tokens>2)||(num_of_tokens<=1)){
	/*
	cout<<"Error: convert_FEP_id_to_number : num_of_tokens != 2 "<<endl;
	cout<<"buffer : "<<buffer<<endl;
	check++;
	*/
	return -1;
    }
  
    if(strcmp(tokens[0],EP->get_FEP_tname())){
	/*
	cout<<"Error: convert_FEP_id_to_number : tokens[0] doesn't match";
	cout<<"with FEP_tname"<<endl;
	cout<<"FEP_tname : "<<EP->get_FEP_tname()<<endl;
	cout<<"tokens[0] : "<<tokens[0]<<endl;
	*/
	check++;
    }

    int FEP_id = atoi(tokens[1]);
    if((FEP_id<0)||(FEP_id>EP->get_FEPsize())){
	/*
	cout<<"Error: convert_FEP_id_to_int : FEP_id not in a valid format "<<endl;
	cout<<"FEP_id : "<<tokens[1]<<endl;
	*/
	check++;
    }
    if(!check){
	return FEP_id;
    }else{
	return -1;
    }
}

int convert_GEP_id_to_int(EmissionProb* const EP, const char* buffer)
{
    if(!buffer)
    {
	return -1;
    }

    char** tokens;
    int num_of_tokens;
    int check = 0;
    check+=break_id_into_tokens(buffer,&tokens,num_of_tokens);
    if(check)
    {
	//cout<<"ERROR: convert_GEP_id_to_number : buffer : "<<buffer<<endl;
	return -1;
    }
    if((num_of_tokens>2)||(num_of_tokens<=1)){
	/*
	cout<<"Error: convert_GEP_id_to_number : num_of_tokens != 2 "<<endl;
	cout<<"buffer : "<<buffer<<endl;
	check++;
	*/
	return -1;
    }
  
    if(strcmp(tokens[0],EP->get_GEP_tname())){
	/*
	cout<<"Error: convert_GEP_id_to_number : tokens[0] doesn't match";
	cout<<"with GEP_tname"<<endl;
	cout<<"GEP_tname : "<<EP->get_GEP_tname()<<endl;
	cout<<"tokens[0] : "<<tokens[0]<<endl;
	*/
	check++;
    }

    int GEP_id = atoi(tokens[1]);
    if((GEP_id<0)||(GEP_id>EP->get_GEPsize())){
	/*
	cout<<"Error: convert_GEP_id_to_int : GEP_id not in a valid format "<<endl;
	cout<<"GEP_id : "<<tokens[1]<<endl;
	*/
	check++;
    }
    if(!check){
	return GEP_id;
    }else{
	return -1;
    }
}
