 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/transitionprob.cpp,v 1.4 2008/12/14 10:39:23 natural Exp $
 */


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "transitionprob.h"
#include "Stack.h"

using namespace std;

//For transition class

// private function

int TransitionProb::set_FTP_from_xml(const char* In_dir,
				     const char* filename)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: set_FTP_from_xml : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * RootNode = 0;
    TiXmlNode * TmpNode = 0; // Node of Transition_Probs
    TiXmlElement* Element = 0; // Element of Transition_Probs

    int length=0;
    int size = 0;
    int count = -1;
    double initprob = 0;

    char* line = new char[Max_line_length];
    
    int NoOfItems = 0;
    char** Items = new char*[Max_number_of_items];
    for(int i=0; i<Max_number_of_items; i++)
    {
	Items[i] = new char[Max_word_length];
	strcpy(Items[i]," ");
    }

    int NumOfTokens = 0;
    char** tokens = NULL;
    
    RootNode = doc.FirstChild("HMMConverter");
    assert(RootNode);
    RootNode = RootNode->FirstChild("model");
    assert(RootNode);

    if(!RootNode->FirstChild("Transition_Probs"))
    {
	cout<<"No Transition Probabilities are defined "<<endl;
	/*
	TP.tname = NULL;
	TP.size = 0;
	TP.name = NULL;
	*/
	FTPsize = 0;
	FTP = NULL;
	FTP_tname = NULL;
    }else{
	TmpNode = RootNode->FirstChild("Transition_Probs");
	assert(TmpNode);
	
	Element = RootNode->FirstChildElement("Transition_Probs");
	assert(Element);
	
	if(!Element->Attribute("id")){
	    cout<<"Error : set_FTP_from_xml : Element->Attribute(id) is NULL "<<endl;
	    check++;
	}else{
	    FTP_tname=Nullstrcpy(Element->Attribute("id"),check);
	    if(check){
		cout<<"Error : set_FTP_from_xml : id attribute "
		    <<"in the TransitionProb tag : "<<Element->Attribute("id")
		    <<"is not a valid"<<endl;
	    }
	}
	
	if(!check){
	    Element->Attribute("size",&size);
	    FTPsize=size;
	    FTP = new FTProb[FTPsize];
	    for(int i=0; i<FTPsize; i++)
	    {
		FTP[i].name = NULL;
		FTP[i].prob = 0;
		FTP[i].train = false;
		FTP[i].exp = NULL;
	    }
	    
	    if(Element->Attribute("file"))
	    {
		// read from .txt file
		char* tran_filename= new char[Max_word_length];
		
		// executable_path is from user
		strcpy(tran_filename,In_dir);
		strcat(tran_filename,Element->Attribute("file"));

		//open file for setting the transition probs
		ifstream instream;
		instream.open(tran_filename);
		if(!instream){
		    cout<<"Error: set_FTP_from_file, can not read file : "<<filename<<endl;
		    check++;
		}
		
		count=0;
		while(!instream.eof())
		{
		    // get transition probs
		    if(check){
			break;
		    }

		    instream.getline(line,Max_line_length-1);
		    if(!strcmp(line,"")) //empty line
		    {
			break;
		    }
		    
		    if(count>=FTPsize){
			cout<<"Error : set_FTP_from_xml : count : "<<count
			    <<" >= FTPsize : "<<FTPsize<<endl;
			check++;
			break;
		    }

		    check+= splitstring(line,NoOfItems,&Items,' ');	    

		    if(NoOfItems==2) // no name
		    {
			FTP[count].name=Nullstrcpy("Not Defined",check);
			if(check)
			{
			    cout<<"Error:: set_FTP_from_file : Nullstrcpy, "
				<<"FTP["<<count<<"].name : Not Defined"<<endl;
			    break;
			}
			double tempparam = atof(Items[1]);
			if((tempparam==0)&&(strcmp(Items[1],"0")))
			{
			    cout<<"Error:: set_transition_prob_from_file : "
				<<"invalid line :"<<line<<endl;
			}else{
			    FTP[count].prob = tempparam;
			}
			FTP[count].pseudocount = 0;
		    }else if(NoOfItems == 3) // with name or pseudocount
		    {
			if((atof(Items[1])!=0) 
			   ||(!strcmp(Items[1],"0"))) // pseudocount
			{
			    double tempparam = atof(Items[1]);
			    if((tempparam==0)&&(strcmp(Items[1],"0")))
			    {
				cout<<"Error:: set_transition_prob_from_file : "
				    <<"invalid line :"<<line<<endl;
				check++;
				FTP[count].prob = 0;
				break;
			    }else{
				FTP[count].prob = tempparam;
			    }	
			    double temppseudoparam = atof(Items[2]);
			    if((temppseudoparam==0)&&(strcmp(Items[2],"0")))
			    {
				cout<<"Error:: set_transition_prob_from_file : "
				    <<"invalid line :"<<line<<endl;
				check++;
				FTP[count].pseudocount = 0;
				break;
			    }else{
				FTP[count].pseudocount = temppseudoparam;
			    }			    
			}else{ //name
			    FTP[count].name=Nullstrcpy(Items[1],check);
			    if(check)
			    {
				cout<<"Error:: set_transition_prob_from_file : Nullstrcpy,"
				    <<"TP.name["<<count<<"] : "<<Items[1]<<endl;
				break;
			    }
			    double tempparam = atof(Items[2]);
			    if((tempparam==0)&&(strcmp(Items[2],"0")))
			    {
				cout<<"Error:: set_transition_prob_from_file : "
				    <<"invalid line :"<<line<<endl;
				check++;
				FTP[count].prob = 0;
				break;
			    }else{
				FTP[count].prob = tempparam;
			    }			    
			    FTP[count].pseudocount = 0;
			}
		    }else if(NoOfItems == 4) 
		    {
			FTP[count].name=Nullstrcpy(Items[1],check);
			if(check)
			{
			    cout<<"Error:: set_transition_prob_from_file : Nullstrcpy,"
				<<"TP.name["<<count<<"] : "<<Items[1]<<endl;
			    break;
			}
			double tempparam = atof(Items[2]);
			if((tempparam==0)&&(strcmp(Items[2],"0")))
			{
			    cout<<"Error:: set_transition_prob_from_file : "
				<<"invalid line :"<<line<<endl;
			    check++;
			    FTP[count].prob = 0;
			    break;
			}else{
			    FTP[count].prob = tempparam;
			}
			double temppseudoparam = atof(Items[3]);
			if((temppseudoparam==0)&&(strcmp(Items[3],"0")))
			{
			    cout<<"Error:: set_transition_prob_from_file : "
				<<"invalid line :"<<line<<endl;
			    check++;
			    FTP[count].pseudocount = 0;
			    break;
			}else{
			    FTP[count].pseudocount = temppseudoparam;
			}	
		    }
		    else{
			cout<<"Error:: set_transition_prob_from_file : line : "
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
			cout<<"Error: set_FTP_from_xml, break_id_into_tokens" <<endl;
			break;
		    }
		    
		    if(strcmp(tokens[0],FTP_tname))
		    {
			cout<<"Error: set_FTP_from_xml, id in .txt file : "<<tokens[0]
			    <<" does not match with id in .xml file : "<<FTP_tname<<endl;
			check++;
			break;
		    }
		    
		    if(atoi(tokens[1])!=count)
		    {
			cout<<"Error: set_FTP_from_file:: TP.id  from file "<<tokens[1]
			    <<" not match with LCount : "<<count<<endl;
			check++;
			break;
		    }
		    
		    count++;
		    NoOfItems = 0;
		}
		instream.close();
		if(tran_filename) delete[] tran_filename;
		tran_filename = NULL;
		if(count!=FTPsize){
		    cout<<"Error : set_FTP_from_xml : count : "<<count
			<<" != FTPsize : "<<FTPsize<<endl;
		    check++;
		}
	    }else{
		
		Element = TmpNode->FirstChildElement("label");
		count=-1;
		while(Element)
		{
		    assert(Element);
		    
		    if (count >= FTPsize ){
			cout<<"Error : set_FTP_from_xml : number of elements : "<<count;
			cout<<"FTPsize : "<<FTPsize<<endl;
			check++;
			break;
		    }
		
		    // read name
		    if (!Element->Attribute("name")){
			FTP[count].name=Nullstrcpy("Not Defined",check);
			if(check){
			    cout<<"Error:: set_FTP_from_file : Nullstrcpy, "
				<<"FTP["<<count<<"].name : Not Defined"<<endl;
			    break;
			}
		    }else if(strlen(Element->Attribute("name"))>Max_word_length){
			cout<<"Error : set_FTP_from_xml : length of Element->Attribute(name) : "<<strlen(Element->Attribute("name"));
			cout<<" exceed Max_word_length : "<<Max_word_length<<endl;
			check++;
			break;
		    }else{		
			count++;
			FTP[count].name=Nullstrcpy(Element->Attribute("name"),check);
			if(check){
			    cout<<"Error : set_FTP_from_xml : "
				<<"name attribute in "<<FTP_tname<<" tag : "
				<<Element->Attribute("name")
				<<"is not a valid string."<<endl;
			    break;
			}			
		    }
		    // read init prob
		    if(!Element->Attribute("initprob")){
			FTP[count].prob=0.0;
		    }else{
			Element->Attribute("initprob",&initprob);
			FTP[count].prob=initprob;
		    }
		    Element = Element->NextSiblingElement();
		}
	    }
	}
    }
    if(line) delete[] line;
    line=NULL;
    
    if(tokens)
    {
	for(int i = 0; i<NumOfTokens; i++){
	    if(tokens[i]) delete [] tokens[i];
	    tokens[i] = NULL;
	}
	delete [] tokens;
    }
    tokens = NULL;
    NumOfTokens = 0;

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

    RootNode = 0;
    TmpNode = 0;
    Element = 0;
    if(check)
    {
	for(int i =0; i<FTPsize; i++)
	{	    
	    if(i<=count)
	    {
		if(FTP[i].name) delete[] FTP[i].name;
		FTP[i].name = NULL;
		if(FTP[i].exp) delete[] FTP[i].exp;
		FTP[i].exp = NULL;
	    }
	}
	if(FTP) delete[] FTP;
	FTP = NULL;
	FTPsize = 0;
    }

    return check;
}

int TransitionProb::get_FTP_parameters_for_training(const char* filename,
						    model_parameters* const MP)
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
	return check;
    }
    TmpNode = SequenceNode->FirstChild("FreeTransitionParameters");
    if(TmpNode)
    {
	Element = TmpNode->FirstChildElement(FTP_tname);
	if(!Element)
	{
	    cout<<"ERROR: TransitionProb:: get_FTP_parameters_for_training: No <"<<FTP_tname<<"> tag "
		<<"under <FreeTransitionParameters> tag."<<endl;
	    check++;
	    return check;
	}
	count = 0;

	while((Element)&&(!check))
	{
	    // get idref and set the train field
	    if(!Element->Attribute("idref"))
	    {
		cout<<"ERROR: TransitionProb:: get_FTP_parameters_for_training: No idref attribute in "
		    <<"the "<<count<<"-th <"<<FTP_tname<<"> tag!"<<endl;
		check++;
	    }else{	    
		index = convert_FTP_id_to_int(this,Element->Attribute("idref"));
		if((index<0)||(index>=FTPsize))
		{
		    cout<<"ERROR: TransitionProb:: get_FTP_parameters_for_training: "
			<<"id: "<<Element->Attribute("idref")
			<<"out of range[0.."<<FTPsize-1<<"] in the "
			<<count<<"-th tag!"<<endl;
		    check++;
		}
	    }
	    
	    if(!Element->Attribute("exp"))
	    {
		cout<<"ERROR: TransitionProb:: get_FTP_parameters_for_training: No exp attribute in "
		    <<"the "<<count<<"-th <"<<FTP_tname<<"> tag!"<<endl;
		check++;
	    }	    // set train and exp field
	    if(!check)
	    {
		FTP[index].train = true;
		if(FTP[index].exp) delete [] FTP[index].exp;
		FTP[index].exp = NULL;

		FTP[index].exp = postfix(Element->Attribute("exp"),check);
		if(check)
		{
		    cout<<"ERROR: TransitionProb:: get_FTP_parameters_for_training: "
			<<" error in postfix for "<<Element->Attribute("exp")<<" in "
			<<count<<"-th "<<"<"<<FTP_tname<<"> tag!"<<endl;
		    check++;		    
		}
		index = -1;
	    }
	    Element = Element->NextSiblingElement();
	    count++;
	}
   
    }

    SequenceNode = 0;
    TmpNode = 0; 
    Element = 0;
    return check;
}

int TransitionProb::get_GTP_parameters_for_training(const char* filename,
						    model_parameters* const MP,
						    Hmm* p)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: get_GTP_parameters_for_training : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * ModelNode = 0;
    TiXmlNode * SequenceNode = 0;
    TiXmlNode * TmpNode = 0;
    TiXmlNode * Node = 0;
    TiXmlNode * fromNode = 0;
    TiXmlNode * toNode = 0;
    TiXmlElement* Element = 0;
    TiXmlElement* fromElement = 0;
    TiXmlElement* toElement = 0;
    int count = 0;
    int trancount = 0;
    int index = -1;
    int fromindex = -1;
    int toindex = -1;

    bool trainFTP = false;

    for(int i=0; i<FTPsize; i++)
    {
	if(this->is_FTP_train(i))
	{
	    trainFTP = true;
	    break;
	}
    }
    
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

    TmpNode = SequenceNode->FirstChild("GroupTransitions");
    if(!TmpNode)
    {
	if((SequenceNode->FirstChild("FreeTransitionParameters"))&&(trainFTP))
	{
	    cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: No <GroupTransitions>"
		<<" tag while there exist free transition parameters to be trained"<<endl;
	    check++;
	}

	if(GTP)
	{
	    for(int i=0; i<GTPsize; i++)
	    {
		if(GTP[i].from) delete[] GTP[i].from;
		GTP[i].from = NULL;
		if(GTP[i].to) delete[] GTP[i].to;
		GTP[i].to = NULL;
		if(GTP[i].Overfrom) delete[] GTP[i].Overfrom;
		GTP[i].Overfrom = NULL;
		if(GTP[i].Overto) delete[] GTP[i].Overto;
		GTP[i].Overto = NULL;
	    }
	    delete[] GTP;
	}
	GTP = NULL;
	GTPsize = 0;
	
	return check;
    }
    
    char* buffer = new char[Max_word_length];

    Element = SequenceNode->FirstChildElement("GroupTransitions");
    
    if(!Element->Attribute("id"))
    {
	cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: No id attribute in "
	    <<"the <GroupTransitions> tag!"<<endl;
	check++;
    }else{

	if(GTP_tname) delete[] GTP_tname;
	GTP_tname = NULL;
	GTP_tname = Nullstrcpy(Element->Attribute("id"),check);
	if(check)
	{
	    cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
		<<" error in Nullstrcpy for : "<<Element->Attribute("id")
		<<" in <GroupTransitions> tag."<<endl;			
	}
    }

    if(!check)
    {

	Element = TmpNode->FirstChildElement(GTP_tname);
	if(!Element)
	{
	    cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: No <"<<GTP_tname<<"> tag"
		<<" under <GroupTransitions> tag."<<endl;
	    check++;
	}

	count = 0;
	// count the number of GTP
	while((Element)&&(!check))
	{
	    Element = Element->NextSiblingElement();
	    count++;		
	}
	GTPsize = count;
	GTP = new GTProb[GTPsize];
	Element = TmpNode->FirstChildElement(GTP_tname);
	Node = TmpNode->FirstChild(GTP_tname);
	count = 0;
	while((Element)&&(!check))
	{
	    if(count>=GTPsize)
	    {
		cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
		    <<"number of <"<<GTP_tname<<"> tags count("<<count<<") out of range[0.."
		    <<GTPsize-1<<"]"<<endl;
		check++;
		break;
	    }

	    // get idref and set the train field
	    if(!Element->Attribute("id"))
	    {
		cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: No id attribute in "
		    <<"the "<<count<<"-th <"<<GTP_tname<<"> tag!"<<endl;
		check++;
	    }
	    if(!check)
	    {
		index = convert_GTP_id_to_int(this,Element->Attribute("id"));
		if((index<0)||(index>=GTPsize))
		{
		    cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
			<<"id: "<<Element->Attribute("id")
			<<"out of range[0.."<<GTPsize-1<<"] in the "
			<<count<<"-th tag!"<<endl;
		    check++;
		}
	    }
	    // get from and to
	    
	    // get NumOffromto
	    trancount = 0;
	    fromNode = Node->FirstChild("from");
	    fromElement = Node->FirstChildElement("from");
	    while((fromNode)&&(!check))
	    {
		toElement = fromNode->FirstChildElement("to");
		while((toElement)&&(!check))
		{
		    if(toElement->Attribute("idref"))
		    {
			if(!cmp_nocase(toElement->Attribute("idref"),"All"))
			{
			    // add count from the number of transition from corresponding state
			    fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
			    if((fromindex<0)||(fromindex>=MP->get_Number_of_States()))
			    {
				cout<<"ERROR:: get_GTP_parameters_for_training : "
				    <<" idref attribute at the "<<count<<"-th "<<GTP_tname
				    <<" tag("<<fromElement->Attribute("idref")<<"out of "
				    <<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
				check++;
			    }else{				
				// count number of out transition from this state except to the end state
				for(int i=0; i<(*p)[fromindex]->get_number_of_next_states(); i++)
				{
				    if((*p)[fromindex]->get_number_of_next_state(i)!=MP->get_Number_of_States()-1)
				    {
					trancount++;
				    }
				}
			    }
			}else{
			    trancount++;
			}

		    }else{
			cout<<"ERROR:: get_GTP_parameters_for_training : "
			    <<" no idref attribute at the "<<count<<"-th "<<GTP_tname
			    <<" tag!"<<endl;
			check++;
		    }
		    toElement = toElement->NextSiblingElement();
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();		
	    }
	    GTP[index].NumOffromto = trancount;

	    GTP[index].from = new int[GTP[index].NumOffromto];
	    GTP[index].to = new int[GTP[index].NumOffromto];
	    // set the from and to arrays 
	    trancount = 0;
	    fromNode = Node->FirstChild("from");
	    fromElement = Node->FirstChildElement("from");
	    while((fromNode)&&(!check))
	    {
  
		toElement = fromNode->FirstChildElement("to");
		fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
		if((fromindex<0)||(fromindex>=MP->get_Number_of_States()))
		{
		    cout<<"ERROR:: get_GTP_parameters_for_training : "
			<<" idref attribute at the "<<count<<"-th "<<GTP_tname
			<<" tag("<<fromElement->Attribute("idref")<<"out of "
			<<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
		    check++;
		    break;
		}
		while((toElement)&&(!check))
		{
		    if(trancount>=GTP[index].NumOffromto)
		    {
			cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
			    <<"number of transitions for group("<<index<<") out of range[0.."
			    <<GTP[index].NumOffromto-1<<"]"<<endl;
			check++; 
			break;
		    }
		    
		    if(toElement->Attribute("idref"))
		    {
			if(!cmp_nocase(toElement->Attribute("idref"),"All"))
			{
			    // add count from the number of transition from corresponding state		
			    
			    for(int i=0; i<(*p)[fromindex]->get_number_of_next_states(); i++)
			    {
				if((*p)[fromindex]->get_number_of_next_state(i)==MP->get_Number_of_States()-1)
				{
				    continue;
				}

				if((trancount>=GTP[index].NumOffromto))
				{
				    cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
					<<"number of transitions for group("<<index<<") out of range[0.."
					<<GTP[index].NumOffromto-1<<"]"<<endl;
				    check++; 
				    break;
				}

				GTP[index].from[trancount] = fromindex;
				GTP[index].to[trancount] = (*p)[fromindex]->get_number_of_next_state(i);
				trancount++;
			       							
			    }
			}else{
			    toindex = convert_State_id_to_int(MP,toElement->Attribute("idref"));
			    if((toindex<0)||(toindex>=MP->get_Number_of_States()))
			    {
				cout<<"ERROR:: get_GTP_parameters_for_training : "
				    <<" idref attribute at the "<<count<<"-th "<<GTP_tname
				    <<" tag("<<toElement->Attribute("idref")<<"out of "
				    <<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
				check++;
				break;
			    }else{
				GTP[index].from[trancount] = fromindex;
				GTP[index].to[trancount] = toindex;
				trancount++;
			    }
			}
		    }else{
			cout<<"ERROR:: get_GTP_parameters_for_training : "
			    <<" no idref attribute at the "<<count<<"-th "<<GTP_tname
			    <<" tag!"<<endl;
			check++;
		    }
		    toElement = toElement->NextSiblingElement();
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
		
	    }
	   
	    // get overfrom and overto
	    // get NumOfOverfromto

	    trancount = 0;
	    fromNode = Node->FirstChild("Overfrom");
	    fromElement = Node->FirstChildElement("Overfrom");
	    while((fromNode)&&(!check))
	    {
		toElement = fromNode->FirstChildElement("Overto");
		while((toElement)&&(!check))
		{
		    if(toElement->Attribute("idref"))
		    {
			if(!cmp_nocase(toElement->Attribute("idref"),"All"))
			{
			    // add count from the number of transition from corresponding state
			    fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
			    if((fromindex<0)||(fromindex>=MP->get_Number_of_States()))
			    {
				cout<<"ERROR:: get_GTP_parameters_for_training : "
				    <<" idref attribute at the "<<count<<"-th "<<GTP_tname
				    <<" tag("<<fromElement->Attribute("idref")<<"out of "
				    <<"range[0.."<<MP->get_Number_of_States()<<"]."<<endl;
				check++;
			    }else{		
				// count number of out transition from this state except to the end state
				for(int i=0; i<(*p)[fromindex]->get_number_of_next_states(); i++)
				{
				    if((*p)[fromindex]->get_number_of_next_state(i)!=MP->get_Number_of_States()-1)
				    {
					trancount++;
				    }
				}
			    }
			}else{
			    trancount++;
			}
		    }else{
			cout<<"ERROR:: get_GTP_parameters_for_training : "
			    <<" no idref attribute at the "<<count<<"-th "<<GTP_tname
			    <<" tag!"<<endl;
			check++;
		    }
		    toElement = toElement->NextSiblingElement();
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
		
	    }
	    
	    GTP[index].NumOfOverfromto = trancount;
	    GTP[index].Overfrom = new int[GTP[index].NumOfOverfromto];
	    GTP[index].Overto = new int[GTP[index].NumOfOverfromto];

	    // set the overfrom and overto arrays
	    trancount = 0;
	    fromNode = Node->FirstChild("Overfrom");
	    fromElement = Node->FirstChildElement("Overfrom");
	    while((fromNode)&&(!check))
	    {
  
		toElement = fromNode->FirstChildElement("Overto");
		fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
		if((fromindex<0)||(fromindex>=MP->get_Number_of_States()))
		{
		    cout<<"ERROR:: get_GTP_parameters_for_training : "
			<<" idref attribute at the "<<count<<"-th "<<GTP_tname
			<<" tag("<<fromElement->Attribute("idref")<<"out of "
			<<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
		    check++;
		    break;
		}
		while((toElement)&&(!check))
		{
		    if(trancount>=GTP[index].NumOfOverfromto)
		    {
			cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
			    <<"number of transitions for group("<<index<<") out of range[0.."
			    <<GTP[index].NumOffromto-1<<"]"<<endl;
			check++; 
			break;
		    }
		    
		    if(toElement->Attribute("idref"))
		    {
			if(!cmp_nocase(toElement->Attribute("idref"),"All"))
			{
			    // add count from the number of transition from corresponding state		
			    
			    for(int i=0; i<(*p)[fromindex]->get_number_of_next_states(); i++)
			    {
				if((*p)[fromindex]->get_number_of_next_state(i)==MP->get_Number_of_States()-1)
				{
				    continue;
				}
				if(trancount>=GTP[index].NumOfOverfromto)
				{
				    cout<<"ERROR: TransitionProb:: get_GTP_parameters_for_training: "
					<<"number of transitions for group("<<index<<") out of range[0.."
					<<GTP[index].NumOfOverfromto-1<<"]"<<endl;
				    check++; 
				    break;				
				}
								
				GTP[index].Overfrom[trancount] = fromindex;
				GTP[index].Overto[trancount] = (*p)[fromindex]->get_number_of_next_state(i);
				trancount++; // skip the transition to end state
				
			    }
			}else{
			    toindex = convert_State_id_to_int(MP,toElement->Attribute("idref"));
			    if((toindex<0)||(toindex>=MP->get_Number_of_States()))
			    {
				cout<<"ERROR:: get_GTP_parameters_for_training : "
				    <<" idref attribute at the "<<count<<"-th "<<GTP_tname
				    <<" tag("<<toElement->Attribute("idref")<<"out of "
				    <<"range[0.."<<MP->get_Number_of_States()-1<<"]."<<endl;
				check++;
				break;
			    }else{
				GTP[index].Overfrom[trancount] = fromindex;
				GTP[index].Overto[trancount] = toindex;
				trancount++;
			    }
			}
		    }else{
			cout<<"ERROR:: get_GTP_parameters_for_training : "
			    <<" no idref attribute at the "<<count<<"-th "<<GTP_tname
			    <<" tag!"<<endl;
			check++;
		    }
		    toElement = toElement->NextSiblingElement();
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
		
	    }
	    
	    count++;
	    Element = Element->NextSiblingElement();
	    Node = Node->NextSibling();
	}
    }

    if(check)
    {
	if(GTP)
	{
	    for(int i=0; i<GTPsize; i++)
	    {
		if(GTP[i].from) delete[] GTP[i].from;
		GTP[i].from = NULL;
		if(GTP[i].to) delete[] GTP[i].to;
		GTP[i].to = NULL;
		if(GTP[i].Overfrom) delete[] GTP[i].Overfrom;
		GTP[i].Overfrom = NULL;
		if(GTP[i].Overto) delete[] GTP[i].Overto;
		GTP[i].Overto = NULL;
	    }
	    delete[] GTP;
	}	
	GTP = NULL;
    }

    return check;
}

int TransitionProb::get_TTP_parameters_for_training(const char* filename,
						    model_parameters* const MP,
						    Hmm* const p)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    int count = 0;

    int fromindex = -1;
    int toindex = -1;
    int statecount = 0;
    int NumOfStates = MP->get_Number_of_States();
    if(!loadOkay){
	cout << "ERROR: get_TTP_parameters_for_training : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }

    TiXmlNode * ModelNode = 0;
    TiXmlNode * TmpNode = 0; 
    TiXmlNode * fromNode = 0;
    TiXmlElement* Element = 0;
    TiXmlElement* fromElement = 0;
    TiXmlElement* toElement = 0;
    
    ModelNode = doc.FirstChild("HMMConverter");
    assert(ModelNode);
    ModelNode = ModelNode->FirstChild("model");
    assert(ModelNode);
    
    TmpNode = ModelNode->FirstChild("Transitions");
    Element = ModelNode->FirstChildElement("Transitions");
    if(!Element)
    {
	cout<<"ERROR: get_TTP_parameters_for_training : no <Transitions> tag in XML file"<<endl;
	check++;
	return check;
    }

    if(Element->Attribute("train"))
    {
	if(!cmp_nocase(Element->Attribute("train"),"All"))
	{	    
	    // train all transitions

	    fromElement = TmpNode->FirstChildElement("from");
	    fromNode = TmpNode->FirstChild("from");
	    
	    // count TTPsize
	    count = 0;
	    statecount = 0;
	 
	    while((fromNode)&&(!check))
	    {
		if(statecount==0)
		{
		    fromNode = fromNode->NextSibling(); // skip begin state
		    fromElement = fromElement->NextSiblingElement();
		    statecount++;
		    continue;
		}
		if(statecount==NumOfStates-1)
		{
		    break;
		}
		
		// train all transitions from this state
		fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
		// assume index should be correct
		// have been checked when building the model
		// calcualte count
		toElement = fromNode->FirstChildElement("to");
		while(toElement)
		{
		    toindex = convert_State_id_to_int(MP,toElement->Attribute("idref"));
		    if(toindex==NumOfStates-1)
		    {
			toElement = toElement->NextSiblingElement();
			continue;
		    }else{
			count++;
			toElement = toElement->NextSiblingElement();
		    }
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();	
		statecount++;
	    }
	    
	    TTPsize = count;
	    TTP = new TTProb[TTPsize];
	    
	    fromNode = TmpNode->FirstChild("from");
	    fromElement = TmpNode->FirstChildElement("from");
	    
	    count = 0;	    
	    // get from and to 
	    statecount = 0;	 
	    while((fromNode)&&(!check))
	    {
		if(statecount==0)
		{
		    fromNode = fromNode->NextSibling(); // skip begin state
		    fromElement = fromElement->NextSiblingElement();
		    statecount++;
		    continue;
		}
		if(statecount==NumOfStates-1)
		{
		    break;
		}			
		// train all transitions from this state
		fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
		// assume index should be correct
		// have been checked when building the model
		
		toElement = fromNode->FirstChildElement("to");
		while(toElement)
		{
		    toindex = convert_State_id_to_int(MP,toElement->Attribute("idref"));
		    if(toindex==NumOfStates-1)
		    {
			toElement = toElement->NextSiblingElement();
			continue;
		    }
		    TTP[count].from = fromindex;
		    TTP[count].to = toindex;
		    TTP[count].score = (*p)[TTP[count].to]->get_transition_score(fromindex);
		    if(toElement->Attribute("pseudoprob"))
		    {				 
			double temppseudoprob = 0;
			temppseudoprob = atof(toElement->Attribute("pseudoprob"));
			if((temppseudoprob<0)||(temppseudoprob>1))
			{
			    cout<<"ERROR: get_TTP_parameters_for_training, pseudoprob for transition "
				<<fromindex<<"->"<<toindex<<" : "<<temppseudoprob
				<<" is not valid."<<endl;
			    check++;
			    break;
			}else if((temppseudoprob==0)
				 &&(strcmp(toElement->Attribute("pseudoprob"),"0")))
			{
			    cout<<"ERROR: get_TTP_parameters_for_training, pseudoprob for transition "
				<<fromindex<<"->"<<toindex<<" : "<<toElement->Attribute("pseudoprob")
				<<" is not valid."<<endl;
			    check++;
			    break;
			}else{
			    TTP[count].pseudoprob = temppseudoprob;
			}
		    }else{
			TTP[count].pseudoprob = 0.0;
		    } 
		    count++;
		    toElement = toElement->NextSiblingElement();
		}

		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
		statecount++;
	    }	      
	}else if(!cmp_nocase(Element->Attribute("train"),"1"))
	{
	    // train some of the transitions
	    fromElement = TmpNode->FirstChildElement("from");
	    fromNode = TmpNode->FirstChild("from");
	    
	    // count TTPsize
	    count = 0;
	    statecount = 0;
	 
	    while((fromNode)&&(!check))
	    {
		if(statecount==0)
		{
		    fromNode = fromNode->NextSibling(); // skip begin state
		    fromElement = fromElement->NextSiblingElement();
		    statecount++;
		    continue;
		}
		if(statecount==NumOfStates-1)
		{
		    break;
		}
		if(fromElement->Attribute("train"))
		{		    
		    if(!cmp_nocase(fromElement->Attribute("train"),"1"))
		    {
			// train all transitions from this state
			fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
			// assume index should be correct
			// have been checked when building the model
			// calculate count
			toElement = fromNode->FirstChildElement("to");
			while(toElement)
			{
			    toindex = convert_State_id_to_int(MP,toElement->Attribute("idref"));
			    if(toindex==NumOfStates-1)
			    {
				toElement = toElement->NextSiblingElement();
				continue;
			    }else{
				count++;
				toElement = toElement->NextSiblingElement();
			    }
			}			
		    }else if(cmp_nocase(fromElement->Attribute("train"),"0"))
		    {
			cout<<"ERROR: get_TTP_parameters_for_training, train attribute at state("
			    <<convert_State_id_to_int(MP,fromElement->Attribute("idref"))
			    <<" should be either y or n but not "<<fromElement->Attribute("train")
			    <<endl;
			check++;  
			break;
		    }
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();	
		statecount++;
	    }
	    
	    TTPsize = count;
	    TTP = new TTProb[TTPsize];
	    
	    fromNode = TmpNode->FirstChild("from");
	    fromElement = TmpNode->FirstChildElement("from");
	    
	    count = 0;
	    
	    // get from and to 

	    statecount = 0;
	 
	    while((fromNode)&&(!check))
	    {
		if(statecount==0)
		{
		    fromNode = fromNode->NextSibling(); // skip begin state
		    fromElement = fromElement->NextSiblingElement();
		    statecount++;
		    continue;
		}
		if(statecount==NumOfStates-1)
		{
		    break;
		}
		
		if(fromElement->Attribute("train"))
		{
		    if(!cmp_nocase(fromElement->Attribute("train"),"1"))
		    {
			// train all transitions from this state
			fromindex = convert_State_id_to_int(MP,fromElement->Attribute("idref"));
			// assume index should be correct
			// have been checked when building the model
			
			toElement = fromNode->FirstChildElement("to");
			while(toElement)
			{
			    if(!strcmp(toElement->Attribute("idref"),"All"))			    
			    {
				
			    }else{
				toindex = convert_State_id_to_int(MP,toElement->Attribute("idref"));
				if(toindex == NumOfStates-1)
				{
				    toElement = toElement->NextSiblingElement();
				    continue;
				}
				TTP[count].from = fromindex;
				TTP[count].to = toindex;
				TTP[count].score = (*p)[TTP[count].to]->get_transition_score(fromindex);
				if(toElement->Attribute("pseudoprob"))
				{				 
				    double temppseudoprob = 0;
				    temppseudoprob = atof(toElement->Attribute("pseudoprob"));
				    if((temppseudoprob<0)||(temppseudoprob>1))
				    {
					cout<<"ERROR: get_TTP_parameters_for_training, pseudoprob for transition "
					    <<fromindex<<"->"<<toindex<<" : "<<temppseudoprob
					    <<" is not valid."<<endl;
					check++;
					break;
				    }else if((temppseudoprob==0)
					     &&(strcmp(toElement->Attribute("pseudoprob"),"0")))
				    {
					cout<<"ERROR: get_TTP_parameters_for_training, pseudoprob for transition "
					    <<fromindex<<"->"<<toindex<<" : "<<toElement->Attribute("pseudoprob")
					    <<" is not valid."<<endl;
					check++;
					break;
				    }else{
					TTP[count].pseudoprob = temppseudoprob;
				    }
				}else{
				    TTP[count].pseudoprob = 0.0;
				} 
			    }
			    toElement = toElement->NextSiblingElement();
			}

			for(int i=0; i<(*p)[fromindex]->get_number_of_next_states(); i++)
			{
			    TTP[count].from = fromindex;
			    TTP[count].to = (*p)[fromindex]->get_number_of_next_state(i);
			    TTP[count].score = (*p)[TTP[count].to]->get_transition_score(fromindex);
			    count++;			    
			}
		    }else if(cmp_nocase(fromElement->Attribute("train"),"0"))
		    {
			cout<<"ERROR: get_TTP_parameters_for_training, train attribute at state("
			    <<convert_State_id_to_int(MP,fromElement->Attribute("idref"))
			    <<" should be either y or n but not "<<fromElement->Attribute("train")
			    <<endl;
			check++;  
			break;
		    }
		}
		fromNode = fromNode->NextSibling();
		fromElement = fromElement->NextSiblingElement();
	    }	    
	}else if(cmp_nocase(Element->Attribute("train"),"0")) // train != y or n
	{ 
	    cout<<"ERROR: get_TTP_parameters_for_training, train attribute at <Transitions> "
		<<" should be either y or n but not "<<Element->Attribute("train")<<endl;
	    check++;
	}
    }else{
        // not training any transition
	TTPsize = 0;
	if(TTP) delete [] TTP;
	TTP = NULL;
    }
 
    if(check)
    {
	TTPsize = 0;
	if(TTP) delete[] TTP;
	TTP = NULL;
    }

    ModelNode = 0;
    TmpNode = 0; 
    fromNode = 0;
    Element = 0;
    fromElement = 0;
    toElement = 0;
    return check;
}

int TransitionProb::check_consistency_of_GTP(Hmm* const p)
{
    int check = 0;
    for(int i =0; i<GTPsize; i++)
    {
	for(int j =0; j<GTP[i].NumOffromto; j++)
	{
	    if((*p)[GTP[i].to[j]]->get_transition_prob(GTP[i].from[j])<0)
	    {
		cout<<"ERROR:: TransitionProb: check_consistency_of_GTP "
		    <<"GTP["<<i<<"].from["<<j<<"] to "<<"GTP["<<i<<"].to["<<j<<"]"
		    <<endl
		    <<"S."<<GTP[i].from[j]<<" to "
		    <<"S."<<GTP[i].to[j]<<" doesn't exist!"<<endl;
		check++;		
	    }
	}	
	for(int j =0; j<GTP[i].NumOfOverfromto; j++)
	{
	    if((*p)[GTP[i].Overto[j]]->get_transition_prob(GTP[i].Overfrom[j])<0)
	    {
		cout<<"ERROR:: TransitionProb: check_consistency_of_GTP "
		    <<"GTP["<<i<<"].Overfrom["<<j<<"]("<<GTP[i].Overfrom[j]<<") to "
		    <<"GTP["<<i<<"].Overto["<<j<<"]("<<GTP[i].Overto[j]<<") doesn't exist!"<<endl;
		check++;		
	    }
	}
    }
    return check;
}

//constructor
TransitionProb::TransitionProb(){
    FTPsize = 0;
    FTP_tname = NULL;
    FTP = NULL;
    GTPsize = 0;
    GTP_tname = NULL;
    GTP = NULL;
    TTPsize = 0;
    TTP = NULL;
}

TransitionProb:: TransitionProb(const char* In_dir,
				const char* filename, 
				int& check)
{
    if(check){
	cout<<"Error: TransitionProb constructor, input variable check : "<<check<<endl;
    }else{
	GTPsize = 0;
	GTP = NULL;
	GTP_tname = NULL;
	TTPsize = 0;
	TTP = NULL;
	check += set_FTP_from_xml(In_dir,filename);
    }
    if(check){
	this->~TransitionProb();
    }
}

TransitionProb::~TransitionProb()
{        
    if(FTP)
    {
	for(int i = 0; i<FTPsize; i++)
	{
	    if(FTP[i].name) delete[] FTP[i].name;
	    FTP[i].name = NULL;
	    if(FTP[i].exp) delete [] FTP[i].exp;
	    FTP[i].exp = NULL;
	}
	delete [] FTP;
    }
    FTP = NULL;
    FTPsize = 0;
    if(FTP_tname) delete[] FTP_tname;
    FTP_tname = NULL;
 
    if(GTP)
    {
	for(int i=0; i<GTPsize; i++)
	{
	    GTP[i].NumOffromto = 0;
	    GTP[i].NumOfOverfromto = 0;
	    if(GTP[i].from) delete[] GTP[i].from;
	    GTP[i].from = NULL;
	    if(GTP[i].to) delete[] GTP[i].to;
	    GTP[i].to = NULL;
	    if(GTP[i].Overfrom) delete[] GTP[i].Overfrom;
	    GTP[i].Overfrom = NULL;
	    if(GTP[i].Overto) delete[] GTP[i].Overto;
	    GTP[i].Overto = NULL;
	}
	delete[] GTP;
    }
    GTP = NULL;
    GTPsize = 0;
    if(GTP_tname) delete [] GTP_tname;
    GTP_tname = NULL;
    
    if(TTP) delete[] TTP;
    TTP = NULL;
    TTPsize = 0;
    
}

//accessors

// for FTP 
int TransitionProb::get_FTPsize() const
{
    return FTPsize;
}

char* TransitionProb::get_FTP_tname() const
{
    return FTP_tname;
}

char* TransitionProb::get_FTP_name(const int i) const 
{
    if((i<0)||(i>=FTPsize))
	return NULL;
    return FTP[i].name;
}

Prob TransitionProb:: get_FTP_prob(const int i) const 
{
    if((i<0)||(i>=FTPsize))
	return -1;
    return FTP[i].prob;
}

FTProb TransitionProb:: get_FTP(const int i) const
{
    if((i<0)||(i>=FTPsize))
    {
	FTProb temp;
	temp.name = NULL;
	temp.prob = 0.0;
	temp.train = false;
	temp.exp = NULL;
    }
    return FTP[i];
}

bool TransitionProb:: is_FTP_train(const int i) const
{
    if((i<0)||(i>=FTPsize))
    {
	return false;
    }
    return FTP[i].train;
}

char* TransitionProb:: get_FTP_exp(const int i) const
{
    if((i<0)||(i>=FTPsize))
    {
	return NULL;
    }
    return FTP[i].exp;
}

Prob TransitionProb:: get_FTP_pseudocount(const int i) const
{
    if((i<0)||(i>=FTPsize))
    {
	return 0;
    }
    return FTP[i].pseudocount;
}

// For GTP
int TransitionProb::get_GTPsize() const
{
    return GTPsize;
}

char* TransitionProb::get_GTP_tname() const
{
    return GTP_tname;
}

int TransitionProb:: get_GTP_NumOffromto(const int i) const
{
    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }
    return GTP[i].NumOffromto;
}

int TransitionProb:: get_GTP_NumOfOverfromto(const int i) const
{
    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }
    return GTP[i].NumOfOverfromto;
}

int TransitionProb:: get_GTP_from(const int i, const int j) const
{
    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }
    if((j<0)||(j>=GTP[i].NumOffromto))
    {
	return -1;
    }
    return GTP[i].from[j];
}

int TransitionProb:: get_GTP_to(const int i, const int j) const
{
    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }
    if((j<0)||(j>=GTP[i].NumOffromto))
    {
	return -1;
    }
    return GTP[i].to[j];
}

int TransitionProb:: get_GTP_Overfrom(const int i, const int j) const
{
    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }
    if((j<0)||(j>=GTP[i].NumOfOverfromto))
    {
	return -1;
    }
    return GTP[i].Overfrom[j];
}

int TransitionProb:: get_GTP_Overto(const int i, const int j) const
{
    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }
    if((j<0)||(j>=GTP[i].NumOfOverfromto))
    {
	return -1;
    }
    return GTP[i].Overto[j];
}

bool TransitionProb:: get_GTP_trained(const int i) const
{

    if((i<0)||(i>=GTPsize))
    {
	return false;
    }

    return GTP[i].trained;
}

double TransitionProb:: get_GTP_score(const int i) const
{

    if((i<0)||(i>=GTPsize))
    {
	return -1;
    }

    return GTP[i].score;
}

// for TTP
int TransitionProb:: get_TTPsize() const
{
    return TTPsize;
}

int TransitionProb:: get_TTP_from(const int i) const
{
    if((i<0)||(i>=TTPsize))
    {
	return -1;
    }
    return TTP[i].from;
}

int TransitionProb:: get_TTP_to(const int i) const
{
    if((i<0)||(i>=TTPsize))
    {
	return -1;
    }
    return TTP[i].to;
}

double TransitionProb:: get_TTP_score(const int i) const
{
    if((i<0)||(i>=TTPsize))
    {
	return 1;
    }
    return TTP[i].score;
}

bool TransitionProb:: get_TTP_trained(const int i) const
{
    if((i<0)||(i>=TTPsize))
    {
	return false;
    }
    return TTP[i].trained;
}

double TransitionProb:: get_TTP_pseudoprob(const int i) const
{
    if((i<0)||(i>=TTPsize))
    {
	return 1;
    }
    return TTP[i].pseudoprob;
}

//mutators , return 0 if no errors

int TransitionProb::set_FTPsize(const int s)
{
    if(s<0){
	cout<<"Error:: TransitionProb:: set_FTPsize, s : "<<s
	    <<" out of range! "<<endl;
	return 1;
    }
    FTPsize=s;
    return 0;
}

int TransitionProb::set_FTP_tname(char* const tn)
{
    int check = 0;

    if(FTP_tname) delete[] FTP_tname;
    FTP_tname=NULL;
    FTP_tname = Nullstrcpy(tn,check);
    if(check!=0){
	cout<<"Error! set_Transition_Prob:: set_FTP_tname"<<endl;
	cout<<"tn : "<<tn<<endl;
    }
    return check;
}

int TransitionProb::set_FTP_name(const int i, char* const n)
{
    int check = 0;
    
    if((i>=FTPsize)||(i<0)){
	cout<<"Error! set_FTP_name:: "<<endl;
	cout<<"i "<<i<<" exceed FTPsize "<<FTPsize<<endl;
	check++;
	return check;
    }

    if(FTP[i].name) delete[] FTP[i].name;
    FTP[i].name=NULL;
    FTP[i].name = Nullstrcpy(n,check);
    if(check!=0){
	cout<<"Error! set_FTP_name:: "<<endl;
	cout<<"Input FTP["<<i<<"].name : "<<n<<endl;
    }
    return check;
}

int TransitionProb::set_FTP_prob(const int i, const Prob p)
{
    if((i>=FTPsize)||(i<0)){
	cout<<"Error! set_FTP_prob:: "<<endl;
	cout<<"i "<<i<<" exceed FTPsize "<<FTPsize<<endl;
	return 1;
    }

    FTP[i].prob = p;
    return 0;
}

int TransitionProb::set_FTP_train(const int i, const bool t)
{
    if((i<0)||(i>=FTPsize))
    {
	cout<<"ERROR: set_FTP_train:: input("<<i
	    <<") out of range[0.."<<FTPsize-1<<"]."<<endl;
	return 1;
    }
    FTP[i].train = t;
    return 0;
}

int TransitionProb::set_FTP_exp(const int i, char* const e)
{
    int check = 0;
    
    if((i>=FTPsize)||(i<0)){
	cout<<"Error! set_FTP_exp:: "<<endl;
	cout<<"i "<<i<<" exceed FTPsize "<<FTPsize<<endl;
	check++;
	return check;
    }

    if(FTP[i].exp) delete[] FTP[i].exp;
    FTP[i].exp=NULL;
    FTP[i].exp = Nullstrcpy(e,check);
    if(check!=0){
	cout<<"Error! set_FTP_exp:: "<<endl;
	cout<<"Input FTP["<<i<<"].exp : "<<e<<endl;
    }
    return check;
}

int TransitionProb::set_FTP_pseudocount(const int i, const Prob p)
{
    int check = 0;
    
    if((i>=FTPsize)||(i<0)){
	cout<<"Error! set_FTP_pseudocount:: "<<endl;
	cout<<"i "<<i<<" exceed FTPsize "<<FTPsize<<endl;
	check++;
	return check;
    }
    FTP[i].pseudocount = p;
    return check;
}

// for GTP

int TransitionProb::set_GTPsize(const int s)
{
    if(s<0)
    {
	cout<<"ERROR: TransitionProb:: set_GTPsize: input("<<s
	    <<") out of range."<<endl;
	return 1;
    }
    GTPsize = s;
    return 0;
}

int TransitionProb::set_GTP_tname(char* const tn)
{
    int check = 0;

    if(GTP_tname) delete[] GTP_tname;
    GTP_tname=NULL;
    GTP_tname = Nullstrcpy(tn,check);
    if(check!=0){
	cout<<"Error! set_Transition_Prob:: set_GTP_tname"<<endl;
	cout<<"tn : "<<tn<<endl;
    }
    return check;
}

int TransitionProb::set_GTP_NumOffromto(const int i, const int num)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_NumOffromto:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
    GTP[i].NumOffromto = num;
    return 0;
}

int TransitionProb::set_GTP_NumOfOverfromto(const int i, const int num)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_NumOfOverfromto:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
    GTP[i].NumOfOverfromto = num;
    return 0;
}

int TransitionProb::set_GTP_from(const int i, const int j, const int f)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_from:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
    if((j<0)||(j>=GTP[i].NumOffromto))
    {
	cout<<"ERROR: set_GTP_from:: input("<<j
	    <<") out of range[0.."<<GTP[i].NumOffromto<<"]."<<endl;
	return 1;
    }
    GTP[i].from[j] = f;
    return 0;
}

int TransitionProb::set_GTP_to(const int i, const int j, const int t)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_to:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
    if((j<0)||(j>=GTP[i].NumOffromto))
    {
	cout<<"ERROR: set_GTP_to:: input("<<j
	    <<") out of range[0.."<<GTP[i].NumOffromto<<"]."<<endl;
	return 1;
    }
    GTP[i].to[j] = t;
    return 0;
}

int TransitionProb::set_GTP_Overfrom(const int i, const int j, const int Of)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_Overfrom:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
    if((j<0)||(j>=GTP[i].NumOfOverfromto))
    {
	cout<<"ERROR: set_GTP_Overfrom:: input("<<j
	    <<") out of range[0.."<<GTP[i].NumOfOverfromto<<"]."<<endl;
	return 1;
    }
    GTP[i].Overfrom[j] = Of;
    return 0;
}

int TransitionProb::set_GTP_Overto(const int i, const int j, const int Ot)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_Overto:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
    if((j<0)||(j>=GTP[i].NumOfOverfromto))
    {
	cout<<"ERROR: set_GTP_Overto:: input("<<j
	    <<") out of range[0.."<<GTP[i].NumOfOverfromto<<"]."<<endl;
	return 1;
    }
    GTP[i].Overto[j] = Ot;
    return 0;
}

int TransitionProb::set_GTP_trained(const int i, const bool t)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_trained:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
 
    GTP[i].trained = t;
    return 0;
}

int TransitionProb::set_GTP_score(const int i, const double s)
{
    if((i<0)||(i>=GTPsize))
    {
	cout<<"ERROR: set_GTP_score:: input("<<i
	    <<") out of range[0.."<<GTPsize-1<<"]."<<endl;
	return 1;
    }
 
    GTP[i].score = s;
    return 0;
}

// for TTP
int TransitionProb:: set_TTPsize(const int s)
{
    if(s<0)
    {
	cout<<"ERROR: TransitionProb:: set_TTPsize: input("<<s
	    <<") out of range."<<endl;
	return 1;
    }
    TTPsize = s;
    return 0;
}
int TransitionProb:: set_TTP_from(const int i, const int f)
{
    if((i<0)||(i>=TTPsize))
    {
	cout<<"ERROR: TransitionProb:: set_TTP_from: input("<<i
	    <<") out of range[0.."<<TTPsize<<"]."<<endl;
	return 1;
    }
    TTP[i].from = f;
    return 0;
}

int TransitionProb:: set_TTP_to(const int i, const int t)
{
    if((i<0)||(i>=TTPsize))
    {
	cout<<"ERROR: TransitionProb:: set_TTP_to: input("<<i
	    <<") out of range[0.."<<TTPsize<<"]."<<endl;
	return 1;
    }
    TTP[i].to = t;
    return 0;
}

int TransitionProb:: set_TTP_score(const int i, const double s)
{
    if((i<0)||(i>=TTPsize))
    {
	cout<<"ERROR: TransitionProb:: set_TTP_score: input("<<i
	    <<") out of range[0.."<<TTPsize<<"]."<<endl;
	return 1;
    }
    TTP[i].score = s;
    return 0;
}
 
int TransitionProb:: set_TTP_trained(const int i, const bool t)
{
    if((i<0)||(i>=TTPsize))
    {
	cout<<"ERROR: TransitionProb:: set_TTP_trained: input("<<i
	    <<") out of range[0.."<<TTPsize<<"]."<<endl;
	return 1;
    }
    TTP[i].trained = t;
    return 0;
}

int TransitionProb:: set_TTP_pseudoprob(const int i, const double p)
{
    if((i<0)||(i>=TTPsize))
    {
	cout<<"ERROR: TransitionProb:: set_TTP_pseudoprob: input("<<i
	    <<") out of range[0.."<<TTPsize<<"]."<<endl;
	return 1;
    }
    TTP[i].pseudoprob = p;
    return 0;
}
   
// operator
TransitionProb & TransitionProb:: operator = (const TransitionProb &tp)
{
    int check = 0;
    if(this!=&tp)
    {    
	if(tp.FTP)
	{
	    if(FTP)
	    {
		for(int i =0; i<FTPsize; i++)
		{
		    if(FTP[i].name) delete [] FTP[i].name;
		    FTP[i].name = NULL;
		    if(FTP[i].exp) delete [] FTP[i].exp;
		    FTP[i].exp = NULL;
		}
		delete [] FTP;
	    }
	    FTP = NULL;
	    FTPsize = tp.FTPsize;
	    if(FTP_tname) delete [] FTP_tname;
	    FTP_tname = NULL;
	    FTP_tname=Nullstrcpy(tp.FTP_tname,check);
	    FTP = new FTProb[FTPsize];
	    for(int i=0 ;i<FTPsize; i++)
	    {
		FTP[i].name = Nullstrcpy(tp.FTP[i].name,check);
		FTP[i].prob = tp.FTP[i].prob;
		FTP[i].train = tp.FTP[i].train;
		FTP[i].exp = Nullstrcpy(tp.FTP[i].exp,check);
		FTP[i].pseudocount = tp.FTP[i].pseudocount;
	    }
	}else{
	    FTPsize = 0;
	    if(FTP_tname) delete[] FTP_tname;
	    FTP_tname = NULL;
	    if(FTP)
	    {
		for(int i =0; i<FTPsize; i++)
		{
		    if(FTP[i].name) delete [] FTP[i].name;
		    FTP[i].name = NULL;
		    if(FTP[i].exp) delete [] FTP[i].exp;
		    FTP[i].exp = NULL;
		}
		delete [] FTP;
	    }
	    FTP = NULL;
	}

	if(tp.GTP)
	{
	    if(GTP)
	    {
		for(int i =0; i<GTPsize; i++)
		{
		    if(GTP[i].from) delete [] GTP[i].from;
		    GTP[i].from = NULL;
		    if(GTP[i].to) delete [] GTP[i].to;
		    GTP[i].to = NULL;
		    if(GTP[i].Overfrom) delete [] GTP[i].Overfrom;
		    GTP[i].Overfrom = NULL;
		    if(GTP[i].Overto) delete [] GTP[i].Overto;
		    GTP[i].Overto = NULL;		    
		}
		delete [] GTP;
	    }
	    GTP = NULL;
	    GTPsize = tp.GTPsize;
	    if(GTP_tname) delete [] GTP_tname;
	    GTP_tname = NULL;
	    GTP_tname=Nullstrcpy(tp.GTP_tname,check);
	    GTP = new GTProb[GTPsize];
	    for(int i=0 ;i<GTPsize; i++)
	    {
		GTP[i].NumOffromto = tp.GTP[i].NumOffromto;
		GTP[i].NumOfOverfromto = tp.GTP[i].NumOfOverfromto;
		for(int j=0; j<GTP[i].NumOffromto; j++)
		{
		    GTP[i].from[j] = tp.GTP[i].from[j];
		    GTP[i].to[j] = tp.GTP[i].to[j];
		}
		for(int j=0; j<GTP[i].NumOfOverfromto; j++)
		{
		    GTP[i].Overfrom[j] = tp.GTP[i].Overfrom[j];
		    GTP[i].Overto[j] = tp.GTP[i].Overto[j];
		}
		GTP[i].trained = tp.GTP[i].trained;
		GTP[i].score = tp.GTP[i].score;
	    }
	}else{
	    GTPsize = 0;
	    if(GTP_tname) delete[] GTP_tname;
	    GTP_tname = NULL;
	    if(GTP)
	    {		
		for(int i =0; i<GTPsize; i++)
		{
		    if(GTP[i].from) delete [] GTP[i].from;
		    GTP[i].from = NULL;
		    if(GTP[i].to) delete [] GTP[i].to;
		    GTP[i].to = NULL;
		    if(GTP[i].Overfrom) delete [] GTP[i].Overfrom;
		    GTP[i].Overfrom = NULL;
		    if(GTP[i].Overto) delete [] GTP[i].Overto;
		    GTP[i].Overto = NULL;		    
		}
		delete [] GTP;
	    }
	    GTP = NULL;
	}

	if(tp.TTP)
	{
	    if(TTP) delete [] TTP;
	    TTP = NULL;
	    TTPsize = tp.TTPsize;
	    TTP = new TTProb[TTPsize];
	    for(int i =0; i<TTPsize; i++)
	    {
		TTP[i].from = tp.TTP[i].from;
		TTP[i].to = tp.TTP[i].to;
		TTP[i].score = tp.TTP[i].score;
		TTP[i].trained = tp.TTP[i].trained;
		TTP[i].pseudoprob = tp.TTP[i].pseudoprob;
	    }	    
	}else{
	    TTPsize = 0;
	    if(TTP) delete [] TTP;
	    TTP = NULL;
	}
    }
    return (*this);
}

// other functions

int TransitionProb::get_parameters_for_training(const char* filename,
						model_parameters* const MP,
						Hmm* const p)
{
    int check = 0;
    if(!filename)
    {
	cout<<"ERROR: TransitionProb:: get_parameters_for_training: "
	    <<"input filename is NULL"<<endl;
	check++;
	return check;
    }
    
    check += get_FTP_parameters_for_training(filename,MP);
    check += get_GTP_parameters_for_training(filename,MP,p);  
    check += get_TTP_parameters_for_training(filename,MP,p);
    
    check += check_consistency_of_GTP(p);
    
    return check;
}

int TransitionProb::derive_FTPs_from_GTPs() 
{
    int check = 0;
    int i =0;

    for(i=0; i<FTPsize; i++)
    {
	if(FTP[i].train)
	{
	    FTP[i].prob = evaluate_FTP_expression(FTP[i].exp,check);
	}
    }
    
    if(check)
    {
	cout<<"ERROR:: TransitionProb class:: derive_FTPs_from_GTPs, "
	    <<"error in evluate_FTP_expression : "<<FTP[i].exp<<endl;
    }
    
    return(check);
}

double TransitionProb::evaluate_FTP_expression(const char* post, int check)
{
    char token;
    double a,b,result;
    Stack<double> opStack;
    int j=0;
    
    if(!post){
	cout<<"Error: evaluate_FTP_expression, input postfix is NULL. "<<endl;
	check++;
	return -1;
    }
    
    char* postfixExp = new char[Max_word_length]; 
    char** posttokens = NULL;
    int NumPosttoken = 0;
    double GTP_prob = -1;
    int postlength = strlen(post);
    int pos=0;

    for(int i=0; i<postlength; i++)
    {
	if(check){
	    break;
	}
	token=post[i];
	switch(token)
	{
	    case ' ': break; // cut space
		
	    case '+' : case '-' : 
	    case '*' : case '/' : 
		if(opStack.empty()){
		    check++;
		    cout<<"Error: evaluate_transition_expression, no operand in the stack! "<<endl;
		    break;
		}else{
		    a=opStack.top();
		    opStack.pop();
		    if(opStack.empty()){
			check++;
			cout<<"Error: evaluate_FTP_exprssion, no operand in the stack! "<<endl;
			break;
		    }
		    else{
			b=opStack.top();
			opStack.pop();
		    }
		}
		if(token=='+')
		    opStack.push(b+a);
		else if(token=='-')
		    opStack.push(b-a);
		else if(token=='*')
		    opStack.push(b*a);
		else if(token=='/')
		    opStack.push(b/a);
		else{
		    check++;
		    cout<<"Error: evaluate_FTP_expression, Unrecognized operator"<<endl;
		}
		break;
		
	    default:
		pos = 0;
		postfixExp[pos] = token;
		pos++;
		for(;;)
		{
		    if((post[i+1]==' ')||(i>=postlength))
			break;
		    i++;
		    token = post[i];
		    postfixExp[pos]=token;
		    pos++;
		}
		postfixExp[pos] = '\0';
		if(isdigit(postfixExp[0])){ // A number
		    opStack.push(atof(postfixExp));
		}else{ // An expression of transition prob, i.e. TP.1
		    check+=break_id_into_tokens(postfixExp,&posttokens,NumPosttoken);
		    if((NumPosttoken>2)||(NumPosttoken<1)){
			check++;
		    }
		    if(check){
			cout<<"Error: evaluate_FTP_expression, break_id_into_tokens" <<endl;
			break;
		    }
		    if(strcmp(GTP_tname,posttokens[0])){
			cout<<"Error: evaluate_FTP_expression, break_id_into_tokens, can not match transition type"<<endl;
			cout<<"TP->tname : "<<GTP_tname<<" posttokens[0] : "<<posttokens[0]<<endl;
			check++;
			break;
		    }
		    GTP_prob = GTP[atoi(posttokens[1])].score;
		    opStack.push(GTP_prob);
		    // free memory
		    for(j=0;j<NumPosttoken;j++){
			if(posttokens[j]) delete[] posttokens[j];
			posttokens[j] = NULL;
		    }
		    if(posttokens) delete[] posttokens;
		    posttokens = NULL;
		} 		
	}
    }
    result = opStack.top();
    opStack.pop();
    if(!opStack.empty()){
	check++;
	cout<<"Error: evaluate_FTP_expression, Unexpected element in the stack"<<endl;
	cout<<opStack.top()<<endl;
    }
    // free memory
    if(postfixExp) delete[] postfixExp;
    postfixExp = NULL;
    return result;
}

void TransitionProb::print(std:: ostream &o) const
{
    o<<"Transition Probs"<<endl;
    o<<"------------------------------------------------------"<<endl<<endl;
    if(FTP)
    {
	o<<"Free_Transition_Prob_name : "<<FTP_tname<<endl;
	o<<"-----------------------------------------------------"<<endl;
	for(int i=0;i<FTPsize;i++)
	{
	    o<<FTP[i].name<<" "<<FTP[i].prob;
	    if(FTP[i].train)
	    {
		o<<" train";
	    }
	    o<<endl;
	}
	o<<endl;
    }
    if(GTP)
    {
	o<<"Group_Transition_Prob_name : "<<GTP_tname<<endl;
	o<<"-----------------------------------------------------"<<endl;
	for(int i=0; i<GTPsize;i++){
	    o<<"Group "<<i<<": "<<endl;
	    for(int j=0; j<GTP[i].NumOffromto; j++)
	    {
		if(j!=0)
		{
		    o<<" + ";
		}
		o<< GTP[i].from[j] <<"->"<<GTP[i].to[j];
	    }
	    o<<endl;
	    o<<"Over"<<endl;	    
	    for(int j=0; j<GTP[i].NumOfOverfromto; j++)
	    {
		if(j!=0)
		{
		    o<<" + ";
		}
		o<<GTP[i].Overfrom[j] <<"->"<<GTP[i].Overto[j];
	    }
	    o<<endl;
	}
	o<<"----------------------------------------------------"<<endl;
    }
    if(TTP)
    {
	o<<"Transitions to be trained"<<endl;
	o<<"----------------------------------------------------"<<endl;
	o<<"Number of transitions to be trained : "<<TTPsize<<endl;
	for(int i =0; i<TTPsize; i++)
	{
	    o<<"transition "<<i<<": "<<TTP[i].from<<" -> "<<TTP[i].to
	     <<" : "<<TTP[i].score<<endl;
	}
    }
    return;
}

int convert_FTP_id_to_int(TransitionProb* const TP, const char* buffer)
{
    char** tokens;
    int num_of_tokens;
    int check = 0;
    check+=break_id_into_tokens(buffer,&tokens,num_of_tokens);
    if(check)
    {
	return -1;
    }
    if((num_of_tokens>2)||(num_of_tokens<=1)){
	return -1;
    }

    if(strcmp(tokens[0],TP->get_FTP_tname())){
	check++;
    }

    int FTP_id = atoi(tokens[1]);
    if((FTP_id<0)||(FTP_id>=TP->get_FTPsize())){
	check++;
    }
    if(!check){
	return FTP_id;
    }else{
	return -1;
    }
}

int convert_GTP_id_to_int(TransitionProb* const TP, const char* buffer)
{
    char** tokens;
    int num_of_tokens;
    int check = 0;
    check+=break_id_into_tokens(buffer,&tokens,num_of_tokens);
    if(check)
    {
	//cout<<"ERROR: convert_GTP_id_to_int : buffer : "<<buffer<<endl;
	return -1;
    }
    if((num_of_tokens>2)||(num_of_tokens<=1)){
	return -1;
    }

    if(strcmp(tokens[0],TP->get_GTP_tname())){
	check++;
    }

    int GTP_id = atoi(tokens[1]);
    if((GTP_id<0)||(GTP_id>=TP->get_GTPsize())){
	check++;
    }
    if(!check){
	return GTP_id;
    }else{
	return -1;
    }
}
