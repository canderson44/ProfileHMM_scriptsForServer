/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define different pairhmm models
*/

#include "models.h"
#include "input_from_files.h"
#include "sequence.h"
#include "hmm_state.h"
#include "Stack.h"
using namespace std;

int get_previous_states(// input
    const int state_number,
    int**  const connection_number,
    int*** const connection_map,
    // output
    int*  number_of_previous_states,
    array<int>* previous_states,
    model_parameters* const MP)
{
    int check = 0;
    int num_of_states = MP->get_Number_of_States();

    if ((state_number < 0) || (state_number > (num_of_states-1)))
    {
	cout << "ERROR: get_previous_states: state_number (" << state_number 
	     << ") out of range (0..." << num_of_states-1 << ").\n" << flush;
	check++;
    }
    if ((*connection_number) == NULL)
    {
	cout << "ERROR: get_previous_states: connection_number is NULL.\n" << flush;
	check++;
    }
    if ((*connection_map) == NULL)
    {
	cout << "ERROR: get_previous_states: connection_map is NULL.\n" << flush;
	check++;
    }
    if (previous_states == NULL)
    {
	cout << "ERROR: get_previous_states: array previous_states is NULL.\n" << flush;
	check++;
    }
    else
    {
	if (previous_states->GetNumberofDimensions() != 1)
	{
	    cout << "ERROR: get_previous_states: array previous_states must be 1-dimensional ("
		 << previous_states->GetNumberofDimensions() << ").\n" << flush;
	    check++;
	}
    }
    if (check == 0)
    {
	int i,k;
	
	// count number of states which connect to state state_number
	
	int number_of_states_connected_to_state = 0;
      
	for (i=0; i< num_of_states; i++)
	{
	    for (k=0; k<(*connection_number)[i]; k++)
	    {
		if ((*connection_map)[i][k] == state_number) {number_of_states_connected_to_state++;}
	    }
	}
	
	// copy number of states which are connected to state state_number into array previous_states
	
	previous_states->SetDimension(0, number_of_states_connected_to_state);
	(*number_of_previous_states) = number_of_states_connected_to_state;

	int count = 0;

	for (i=0; i<num_of_states; i++)
	{
	    for (k=0; k<(*connection_number)[i]; k++)
	    {
		if ((*connection_map)[i][k] == state_number) 
		{
		    previous_states->SetElement(count, i);
		    count++;
		}
	    }
	}
    }
    return(check);
}

int get_next_states(// input
    const int          state_number,
    int**  const connection_number,
    int*** const connection_map,
    // output
    int*        number_of_next_states,
    array<int>* next_states,
    model_parameters* const MP)
{
    int check = 0;
    int num_of_states = MP->get_Number_of_States();

    if ((state_number < 0) || (state_number > (num_of_states-1)))
    {
	cout << "ERROR: get_next_states: state_number (" << state_number 
	     << ") out of range (0..." << num_of_states-1 << ").\n" << flush;
	check++;
    }
    if ((*connection_number) == NULL)
    {
	cout << "ERROR: get_next_states: connection_number is NULL.\n" << flush;
	check++;
    }
    if ((*connection_map) == NULL)
    {
	cout << "ERROR: get_next_states: connection_map is NULL.\n" << flush;
	check++;
    }
    if (next_states == NULL)
    {
	cout << "ERROR: get_next_states: array next_states is NULL.\n" << flush;
	check++;
    }
    else
    {
	if (next_states->GetNumberofDimensions() != 1)
	{
	    cout << "ERROR: get_next_states: array next_states must be 1-dimensional ("
		 << next_states->GetNumberofDimensions() << ").\n" << flush;
	    check++;
	}
    }
    if (check == 0)
    {
	// copy number of states to which state is connected into array next_states
	
	next_states->SetDimension(0, (*connection_number)[state_number]);
	(*number_of_next_states) = (*connection_number)[state_number];
	
	for (int k=0; k<(*connection_number)[state_number]; k++)
	{
	    next_states->SetElement(k, (*connection_map)[state_number][k]);
	}
    }
    return(check);
}

int get_connection_map(const char* file_name,
		       int**  const array_connection_number,
		       int*** const array_connection_map,
		       model_parameters* const MP)
{
    int check = 0;
    int i, j, k = 0;
    int NumNextState=0;
    int state_count=0;
    int NumState = MP->get_Number_of_States();
    
    TiXmlDocument doc(file_name);
    bool loadOkay = doc.LoadFile();
    TiXmlNode * node=0; 
    TiXmlNode * TransitionNode=0;
    TiXmlNode * FromNode=0;
    TiXmlElement * FromElement=0;
    TiXmlElement * ToElement = 0;
    TiXmlElement * CountElement = 0;
    
    if(!loadOkay){
	cout << "ERROR: models.cpp : cannot open file " << file_name << ".\n" << flush;
	check++;
    }else{
	if ((*array_connection_number) != NULL) {
	    cout << "ERROR:get_connection_map: array connection_number != NULL.\n" << flush;
	    check++;
	}
	if ((*array_connection_map) != NULL) {
	    cout << "ERROR:get_connection_map: array connection_map != NULL.\n" << flush;
	    check++;
	}
	if (check == 0) {

	    (*array_connection_number) = new int[NumState];
	    (*array_connection_map)    = new int*[NumState];
		
	    for (i=0; i<NumState; i++) {
		(*array_connection_map)[i] = new int[NumState];
	    }
	    
	    for (i=0; i<NumState; i++) {
		(*array_connection_number)[i] = 0;
		
		for (j=0; j<NumState; j++) {
		    (*array_connection_map)[i][j] = 0;
		}
	    }
   
	    node = doc.FirstChild("HMMConverter");
	    assert(node);
	    node = node->FirstChild("model");
	    assert(node);
	    TransitionNode=node->FirstChild("Transitions");
	    assert(TransitionNode);
	    FromNode=TransitionNode->FirstChild("from");
	    assert(FromNode);
	    FromElement = TransitionNode->FirstChildElement("from");
	    state_count=0;
	    i=0;
	    while(FromElement){
		if(check){
		    break;
		}
		assert(FromElement);
	
		i=convert_State_id_to_int(MP,FromElement->Attribute("idref"));
		if(i<0){
		    cout<<"Error : get_connection_map : "
			<<"idref attribute in the "<<state_count
			<<"-th"<<"transition tag : "
			<<FromElement->Attribute("idref")
			<<" is invalid"<<endl;
		    check++;
		    break;
		}
		
		// get NumNextState
		CountElement = FromNode->FirstChildElement("to");
		NumNextState = 0;
		while(CountElement)
		{
		    NumNextState++;
		    if((CountElement->Attribute("idref"))&&(NumNextState==1))
		    {
			if(!cmp_nocase(CountElement->Attribute("idref"),"ALL"))
			{
			    NumNextState = MP->get_Number_of_States()-2;
			    break;
			}
		    }
		    CountElement = CountElement->NextSiblingElement();
		}
	
		(*array_connection_number)[i]=NumNextState;
		ToElement=FromNode->FirstChildElement("to");
		k=0;
		while(ToElement){
		    assert(ToElement);
		    if(!strcmp(ToElement->Attribute("idref"),"All")){
			for(j=0;j<NumState-1;j++){
			    if(j!=0){
				(*array_connection_map)[i][k]=j;
				k++;
			    }
			}
		    }else{
			j=convert_State_id_to_int(MP,ToElement->Attribute("idref"));
			(*array_connection_map)[i][k]=j;
			k++;
		    }
		    ToElement = ToElement-> NextSiblingElement();
		}
	
		FromNode = FromNode -> NextSibling();
		FromElement = FromElement -> NextSiblingElement();
		state_count++;
	    }
	    if(state_count>NumState){
		cout<<"Error:: get_connection_map : Number of states : "<<NumState;
		cout<<" less than number of states defined for the transitions : "<<state_count<<endl;
		check++;
	    }
	}
    }
    node=0; 
    TransitionNode=0;
    FromNode=0;
    FromElement=0;
    ToElement = 0;

    return(check);
}

int set_up_topology_of_model(const char* filename,
			     Hmm* const real,
			     model_parameters* const MP)
{
    int check = 0;
    
    // general variables
    
    // set up states
    
    int*  connection_number  = NULL;
    int** connection_map     = NULL;
    
    
    check+=get_connection_map(filename,&connection_number,
			      &connection_map,MP);
    
    int n_of_previous_states = 0;
    array<int> previous_states(1); 
  
    int NumStateLabel = MP->get_Total_Number_of_Annotation_Labels();
    int NumState = MP->get_Number_of_States();

    TiXmlNode * RootNode=0; 
    TiXmlNode * StateNode=0;
    TiXmlNode * StateLabelNode = 0; 
    TiXmlNode * SpecialNode = 0;
    TiXmlElement * StateElement=0;
    TiXmlElement * StateLabelElement=0;
    TiXmlElement * SpecialElement = 0;
    
    char** StateLabelType = new char*[NumStateLabel];
    array<int> *state_label= new array<int>[NumStateLabel];

    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    
    if(!loadOkay){
	cout << "ERROR: set_up_topology_of_model : cannot open file " << filename << ".\n" << flush;
	check++;
    }
    
    // set AnnotationLabelType
    
    for(int i=0;i<NumStateLabel; i++){
	StateLabelType[i] = Nullstrcpy(MP->get_Annotation_Label_setname(i),check);
	if(check){
	    cout<<"ERROR: set_up_topology_of_model, StateLabelType["<<i<<"] : "<<MP->get_Annotation_Label_setname(i)<<endl;
	    break;
	}
    }
    
    // to get the details for all the states	
    
    if(check==0){
	
	int i = 0;
	int xdim = 0;
	int ydim = 0;
	int StateCount = 0;
	int LabelDim = 0;
	int LabelCount = 0;
	int LabelId = 0;
	bool special = false;
	int n_of_child_state = 0;
	int n_of_virtural_state = MP->get_Number_of_States()+1;
	int StateId = -1;
	
	RootNode = doc.FirstChild("HMMConverter");
	assert(RootNode);
	RootNode = RootNode->FirstChild("model");
	assert(RootNode);
	StateNode = RootNode->FirstChild("States");
	assert(StateNode);
      
	StateElement = StateNode->FirstChildElement("State");
	StateNode = StateNode->FirstChild("State");
	
	while(StateElement)
	{	    
	    if(check){
		break;
	    }
	    assert(StateNode);
	    assert(StateElement);
	    // get the dimensions of the state
	    if(StateCount>=NumState){
		cout<<"Error : set_up_topology_of_models, number of state : "<<StateCount;
		cout<<" exceed Number_of_States : "<<NumState<<endl;
		check++;
		break;
	    }
	    
	    if(!StateElement->Attribute("xdim"))
	    {	
		xdim = 0;
	    }else{
		StateElement->Attribute("xdim",&xdim);
	    }
	    
	    if(!StateElement->Attribute("ydim"))
	    {
		ydim = 0;
	    }else if((!MP->is_PairHMM())&&(StateElement->Attribute("ydim")))
	    {
		int tempdim;
		StateElement->Attribute("yhim",&tempdim);
		if(tempdim>0)
		{
		    cout<<"Error : set_up_topology_of_models, number of state : "<<StateCount;
		    cout<<" ydim("<<tempdim<<") > 0 while the model is not PairHMM."<<endl;
		    check++;
		    break;
		}else{
		    ydim = 0;
		}
	    }else{
		StateElement->Attribute("ydim",&ydim);
	    }

	    if(!StateElement->Attribute("id"))
	    {
		cout<<"Error : set_up_topology_to_model, id attribute in "
		    <<"State tag is NULL "<<endl;
		check++;
	    }else{
		StateId = convert_State_id_to_int(MP,StateElement->Attribute("id"));
		if(StateId<0)
		{
		    cout<<"Error : set_up_topology_to_model, id attribute "
			<<StateElement->Attribute("id")
			<<"in State tag is invalid "<<endl;
		    check++;
		}
	    }

	    // get special
	    
	    if(!StateElement->Attribute("special"))
	    {
		special = false;
	    }else{
		if(!cmp_nocase(StateElement->Attribute("special"),"1"))
		{
		    special = true;
		}else if(!cmp_nocase(StateElement->Attribute("special"),"0"))
		{
		    special = false;
		}else{
		    cout<<"Error : set_up_topology_to_model,state("<<StateId
			<<") special field must be y or n "<<endl;
		    check++;
		    break;
		}		
	    }

	    // get annotation labels
	    for(i=0;i<NumStateLabel;i++){
		if(check){
		    break;
		}
		StateLabelNode = StateNode->FirstChild(MP->get_Annotation_Label_setname(i));
		
		if(!StateLabelNode)
		{
		    state_label[i].SetNumberofDimensions(1);
		    state_label[i].SetDimension(0,0);
		    continue;
		}
	
		LabelDim = xdim+ydim;

		StateLabelElement = StateLabelNode->FirstChildElement("label");
		assert(StateLabelElement);
		if(LabelDim<0){
		    cout<<"Error: set_up_topology_of_model, LabelDim of : "<<MP->get_Annotation_Label_setname(i)<<" < 0"<<endl;
		    check++;
		    break;
		}else if(LabelDim==0)
		{
		    LabelDim = 1;   // allow setting labels for silent state 
		}
		state_label[i].SetNumberofDimensions(1);
		state_label[i].SetDimension(0,LabelDim);
		LabelCount = 0;
		while(StateLabelElement){
		    if(check){
			break;
		    }
		    if(LabelCount >= LabelDim){
			cout<<"Error: set_up_topology_of_model, LabelCount : "<<LabelCount<< " >= LabelDim : "<<LabelDim<<endl;
		      check++;
		      break;
		    }
	
		    LabelId =convert_Annotation_Label_id_to_int(MP,
							   StateLabelElement->Attribute("idref"),
							   i);
		    if(LabelId<0){
			cout<<"Error: set_up_topology_of_model : "
			    <<"idref attribute of Annotation_Label tag "
			    <<StateLabelElement->Attribute("idref")
			    <<"of"<<StateId<<"=th state is invalid"<<endl;
			check++;
			break;
		    }
		    state_label[i].SetElement(LabelCount,LabelId);
		    LabelCount++;
		    StateLabelElement = StateLabelElement->NextSiblingElement();
		}
		// report error as label is not set for all the dimension
		if(!((LabelCount==LabelDim)||(LabelCount==1)))
		{
		    cout<<"Error: set_up_topologu_of_model, number of state : "<<StateCount;
		    cout<<" LabelCount("<<LabelCount<<") not equal to 1 or LabelDim("<<LabelDim<<")."<<endl;
		    check++;
		    break;
		}else if((LabelCount!=LabelDim)&&(LabelCount==1))
		{
		    while(LabelCount<LabelDim)
		    {
			state_label[i].SetElement(LabelCount,LabelId);
			LabelCount++;
		    }
		}
	    }
	    
	    check+=get_previous_states(StateId,
				       &connection_number,
				       &connection_map,
				       &n_of_previous_states,
				       &previous_states,
				       MP);
	 	
	  
	    Hmm_State temp(xdim,
			   ydim,
			   MP->get_Alphabet_size(),
			   StateId,
			   NumState,
			   n_of_previous_states,
			   previous_states,			
			   NumStateLabel,
			   state_label,
			   StateLabelType,
			   MP);
				   
	    check+=real->implement_state(StateCount, &temp);		    
	    
	    if(special)
	    {
		if((xdim>0)&&(ydim>0)){
		    (*real)[StateCount]->set_special_emissions(n_of_virtural_state,n_of_virtural_state+1);
		    n_of_virtural_state+=2;
		    n_of_child_state+=2;
		}else if(xdim>0){
		    (*real)[StateCount]->set_special_emissions(StateCount);
		    n_of_child_state++;
		}else if(ydim>0){
		    (*real)[StateCount]->set_special_emissions(StateCount);
		    n_of_child_state++;
		}
	    }
	    
	    MP->set_Number_of_Special_Emissions(n_of_child_state);
	    MP->set_Max_n_of_Child_State(n_of_virtural_state);
	    
	    if(check){
		cout<<"Error: set_up_topology_of_model : problems implementing state : "<<StateCount<<endl;
		cout<<"Stop implementing the remaining states. "<<endl;
		break;
	    }
	    StateNode = StateNode->NextSibling();
	    StateElement = StateElement->NextSiblingElement();
	    StateCount++;
	}
	if(StateCount!=NumState){
	    cout<<"Error: set_up_topology_of_model, inconsistency in number of states."<<endl;
	    cout<<"StateCount : "<<StateCount<<" NumState : "<<NumState<<endl;
	    check++;
	}
    }

    if(!check){
	check+= real->build_connected_pairhmm();
	if(check){
	    cout<<"Error:: models.cpp: errors in build_connected_pairhmm() "<<endl;
	}
    }
    
    //free memory
    
    RootNode=0; 
    StateNode=0;
    StateLabelNode = 0;
    SpecialNode = 0;
    StateElement=0;
    StateLabelElement=0;
    SpecialElement = 0;
    
    if (connection_map) { 
	for (int i=0; i<connection_number[i]; i++) {
	    if (connection_map[i]) delete [] connection_map[i];
	  connection_map[i] = NULL;
	}
	delete [] connection_map;
	connection_map = NULL;
    }
    
    if (connection_number) delete [] connection_number;
    connection_number = NULL;
       
    if(state_label) delete[] state_label;
    state_label = NULL;
    
    for(int i=0; i<NumStateLabel; i++)
    {
	if(StateLabelType[i]) delete [] StateLabelType[i];
	StateLabelType[i] = NULL;
    }
    if(StateLabelType) delete [] StateLabelType;
    StateLabelType = NULL;
    
    return(check);
}

double evaluate_transition_expression(const char* post, TransitionProb* const TP, int& check)
{
    char token;
    double a,b,result;
    Stack<double> opStack;
    int j=0;
    
    if(!post){
	cout<<"Error: evaluate_transition_expression, input postfix is NULL. "<<endl;
	check++;
	return -1;
    }
    
    char* postfixExp = new char[Max_word_length]; 
    char** posttokens = NULL;
    int NumPosttoken = 0;
    double TP_prob = -1;
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
			cout<<"Error: evaluate_transition_exprssion, no operand in the stack! "<<endl;
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
		    cout<<"Error: evaluate_transition_expression, Unrecognized operator"<<endl;
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
			cout<<"Error: evaluate_transition_expression, break_id_into_tokens" <<endl;
			break;
		    }
		    if(strcmp(TP->get_FTP_tname(),posttokens[0])){
			cout<<"Error: evaluate_transition_expression, break_id_into_tokens, can not match transition type"<<endl;
			cout<<"TP->tname : "<<TP->get_FTP_tname()<<" posttokens[0] : "<<posttokens[0]<<endl;
			check++;
			break;
		    }
		    TP_prob = TP->get_FTP_prob(atoi(posttokens[1]));
		    opStack.push(TP_prob);
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
	cout<<"Error: evaluate_transition_expression, Unexpected element in the stack"<<endl;
	cout<<opStack.top()<<endl;
    }
    // free memory
    if(postfixExp) delete[] postfixExp;
    postfixExp = NULL;
    return result;
}

int set_up_transition_probs_of_model(Hmm* const real,
				     const char* filename,
				     model_parameters* const MP)

  
{
    int check = 0;
    
    
    if (real == NULL){
	cout << "ERROR: set_up_transition_probs_of_model: Hmm real is NULL.\n" << flush;
	check++;
	return check;
    }
    // set transition probs 
    // get transition_probs_expression
    
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    
    if(!loadOkay){
	cout << "ERROR: set_up_transitions_of_model : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    int i, j, k = 0;
    int NumNextState=0;
    int StateCount=0;
    int NumState = MP->get_Number_of_States();
    
    if (check == 0) {

	TiXmlNode * RootNode=0; 
	TiXmlNode * TransitionNode=0;
	TiXmlNode * FromNode=0;
	TiXmlElement * FromElement=0;
	TiXmlElement * ToElement = 0;
	RootNode = doc.FirstChild("HMMConverter");
	assert(RootNode);
	RootNode = RootNode->FirstChild("model");
	assert(RootNode);
	TransitionNode=RootNode->FirstChild("Transitions");
	assert(TransitionNode);
	FromNode=TransitionNode->FirstChild("from");
	assert(FromNode);
	FromElement = TransitionNode->FirstChildElement("from");
	StateCount=0;
	i=0;
	while(FromElement){
	    if(check){
		break;
	    }
	    assert(FromElement);

	    i=convert_State_id_to_int(MP,FromElement->Attribute("idref"));
	    if(i<0){
		cout<<"Error: set_up_transition_prob_of_model, convert_State_id_to_int, state "<<StateCount<<endl;
		check++;
		break;
	    }
	  
	    ToElement=FromNode->FirstChildElement("to");
	    k=0;
	    while(ToElement){
		assert(ToElement);

		if(!strcmp(ToElement->Attribute("idref"),"All")){
		    for(j=0;j<NumState-1;j++){
			if(j!=0){
		 
			    (*real)[i]->set_transition_probs_expression(j,postfix(ToElement->Attribute("exp"),check));
 
			    if(check){
				cout<<"Error: set_up_transition_prob_of_model: postfix in "
				    <<"the "<<StateCount<<"-th state"<<endl;
				break;
			    }
			   
			    k++;
			}
		    }
		}else{
		    
		    j=convert_State_id_to_int(MP,ToElement->Attribute("idref"));
		    if(j<0){
			cout<<"Error: set_up_transition_prob_of_model : "
			    <<"idref attribute in To tag of "<<StateCount
			    <<"-th state : "<<ToElement->Attribute("idref")
			    <<"is invalid."<<endl;
			check++;
			break;
		    }

		    (*real)[i]->set_transition_probs_expression(j,postfix(ToElement->Attribute("exp"),check));
		    
		    if(check){
			cout<<"Error: set_up_transition_prob_of_model: postfix in "
			    <<"the "<<StateCount<<"-th state"<<endl;
			break;
		    }

		    k++;
		}
		ToElement = ToElement-> NextSiblingElement();
	    }

	    FromNode = FromNode->NextSibling();
	    FromElement = FromElement->NextSiblingElement();
	    StateCount++;
	}
	if(StateCount > NumState){
	    cout<<"Error: set_up_transitions_of_model, number of states defined for transitions: "<<StateCount;
	    cout<<" >  number of states defined in the model : "<<NumState<<endl;
	    check++;
	}

	RootNode=0; 
	TransitionNode=0;
	FromNode=0;
	FromElement=0;
	ToElement = 0;
    }
    return(check);
}

Prob sum_over(Hmm* const real,
	      EmissionProb* const EP,
	      const int From,
	      const int CurState,
	      const int index,
	      const bool state,
	      int& check
	      )
{
    Prob probability=0;
    
    if(real==NULL){
	cout<<"ERROR: Sum_over: Hmm real is NULL.\n"<<flush;
	check++;
    }
    if(!check){
	// initialize a vector to store different combination of charactors
	// with the same dimension of the emission_prob vector of from state
	int SumOverDim = (*real)[CurState]->get_num_sum_over();
	int FixedDim = (*real)[CurState]->get_sum_over_FromPos_dim(0);
	long alphabet = (*real)[CurState]->get_alphabet();
	long NumSumOver = static_cast<long>(pow(static_cast<float>(alphabet),static_cast<float>(SumOverDim)));
	long NumFixed = static_cast<long>(pow(static_cast<float>(alphabet),static_cast<float>(FixedDim)));
	if(state)
	{
	    if((SumOverDim+FixedDim)!=(*real)[From]->get_number_of_dimensions_of_emission_probs()){
		cout<<"ERROR! In sum_over, dimension of SumOverDim + FixedDim != dimension of emission prob in the from state.\n";
		cout<<"SumOverDim = "<<SumOverDim<<" FixedDim = "<<FixedDim<<endl;
		cout<<"Dimension of emission prob in the from_state = "<<(*real)[From]->get_number_of_dimensions_of_emission_probs()<<endl;
		check++;
	    }
	}else{
	    if((SumOverDim+FixedDim)!=EP->get_FEP_dim(From))
	    {
		cout<<"ERROR! In sum_over, dimension of SumOverDim + FixedDim != dimension of emission prob in the from state.\n";
		cout<<"SumOverDim = "<<SumOverDim<<" FixedDim = "<<FixedDim<<endl;
		cout<<"Dimension of the free emission parameter = "<<EP->get_FEP_dim(From)<<endl;
		check++;
	    };
	}
	// set corr_pos
	if(!check){
	    int i;
	    
	    int temp;
	    array<int> FromIndex;
	    FromIndex.SetNumberofDimensions(1);
	    FromIndex.SetDimension(0,SumOverDim+FixedDim);
	    
	    int j;
	    probability =0;
	    temp = index;
	    // set fixed value
	    for(i=0;i<FixedDim;i++){
		FromIndex.SetElement((*real)[CurState]->get_sum_over_FromPos(i),temp%alphabet);
		temp=temp/alphabet;
	    }
	    if(state)
	    {
		for(i=0; i<NumSumOver;i++){
		    // set sum_over_value
		    temp = i;
		    for(j=0; j<SumOverDim;j++){
			FromIndex.SetElement((*real)[CurState]->get_sum_over_pos(j),temp%alphabet);
			temp = temp/alphabet;
		    }
		    probability+=(*real)[From]->get_emission_prob(FromIndex);
		}
	    }else{
		for(i=0; i<NumSumOver;i++){
		    // set sum_over_value
		    temp = i;
		    for(j=0; j<SumOverDim;j++){
			FromIndex.SetElement((*real)[CurState]->get_sum_over_pos(j),temp%alphabet);
			temp = temp/alphabet;
		    }
		    probability += EP->get_FEP_prob(From,FromIndex);
		}
	    }
	}
    }
    return (probability);
}

int set_up_emission_probs_of_model(Hmm* const real,
				   const char* filename,
				   model_parameters* const MP,
				   EmissionProb* const EP)

{

    int check = 0;
  
    // set emission probs 
    // get emission_probs_expression
    if(real==NULL){
	cout<<"ERROR: set_up_emission_probs_of_model: Hmm is NULL.\n"<<flush;
	check++;
	return check;
    }
	
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
	
    if(!loadOkay){
	cout << "ERROR: pairhmm.xml : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }

    int i, j, k = 0;
    int NumNextState=0;
    int NumState = 0;
    int StateCount=0;
    int NumSum=0;
    int index = 0;
    char* StrExp;
    
    char* SumPos = new char[Max_word_length];
    char* ExpBuffer = new char[Max_word_length];

    char* str_id = new char[Max_word_length];
    TiXmlNode * node=0;
    TiXmlNode * EmissionsNode=0;
    TiXmlNode * StateNode=0;
    TiXmlNode * EmissionProbNode=0;
    TiXmlElement * EmissionProbElement=0;
    TiXmlElement * StateElement=0;
    TiXmlElement * SumOverElement=0;
    TiXmlElement * ProductElement=0;

    node = doc.FirstChild("HMMConverter");
    assert(node);
    node = node->FirstChild("model");
    assert(node);
    StateNode=node->FirstChild("States");
    assert(StateNode);
    StateElement = node->FirstChildElement("States");
    assert(StateElement);
    StateElement->Attribute("size",&NumState);
    StateElement=StateNode->FirstChildElement("State");
    assert(StateElement);
    StateNode=StateNode->FirstChild("State");
    assert(StateNode);
    
    int StateId;
    i=0; // checking
    while(StateElement){

	assert(StateElement);

	StateId = convert_State_id_to_int(MP,StateElement->Attribute("id"));
	if(StateId<0){
	    cout<<"Error: set_up_emission_probs_of_model, convert_State_id_to_int : "<<endl;
	    cout<<"StateCount : "<<StateCount<<endl;
	    break;
	}
	EmissionProbNode = StateNode->FirstChild("State_Emission_Probs");
	if(!EmissionProbNode) 
	{
	    if((*real)[StateId]->get_letters_to_read()==0) // silent state
	    {
		(*real)[StateId]->set_emission_probs_expression("none");
		(*real)[StateId]->set_emission_sum_over(false);
		(*real)[StateId]->set_emission_product(false);	    
	    }else{
		if(StateId<=EP->get_FEPsize())
		{
		 
		    strcpy(ExpBuffer,EP->get_FEP_tname());
		    strcat(ExpBuffer,".");
		    sprintf(str_id,"%d",StateId-1);
		    strcat(ExpBuffer,str_id);
		    
		    (*real)[StateId]->set_emission_probs_expression(ExpBuffer);
		    (*real)[StateId]->set_emission_sum_over(false);
		    (*real)[StateId]->set_emission_product(false);
		}else{
		    cout<<"Error : models.cpp set_up_emission_probs_of_model : state("
			<<StateId<<") has no emission prob declared, while it reads "
			<<(*real)[StateId]->get_letters_to_read()<<" letters, "
			<<"and the StateId exceed the number of declared emission probs : "
			<<EP->get_FEPsize()<<endl;
		    check++;
		    break;
		}
	    }
	}else{
	    EmissionProbElement = StateNode->FirstChildElement("State_Emission_Probs");
	    assert(EmissionProbElement);
	    
	    // read transition expressions
	    if((!EmissionProbElement->Attribute("GetFrom"))&&
	       (!EmissionProbElement->Attribute("SumOver"))&&
	       (!EmissionProbElement->Attribute("Product")))
	    {
		cout<<"Error: set_up_emission_probs_of_model, there is no valid attribute set for Emission_Probs of state("<<StateId<<")."<<endl;
		check++;
	    }else{ 
		if(EmissionProbElement->Attribute("GetFrom"))
		{	
		    (*real)[StateId]->set_emission_probs_expression(EmissionProbElement->Attribute("GetFrom"));
		}else{
		    (*real)[StateId]->set_emission_probs_expression("none");
		}
	    
		if(EmissionProbElement->Attribute("SumOver")){
		    // read the sum over field
	
		    if(!cmp_nocase(EmissionProbElement->Attribute("SumOver"),"1"))
		    {
			(*real)[StateId]->set_emission_probs_expression("SumOver");
			(*real)[StateId]->set_emission_sum_over(true);
			// set variables related to sum_over
			SumOverElement = EmissionProbNode->FirstChildElement("SumOver");
	
			index = convert_State_id_to_int(MP,SumOverElement->Attribute("From"));
			if(index<0)
			{
			    index = convert_FEP_id_to_int(EP,SumOverElement->Attribute("From"));
			    if(index<0)
			    {
				cout<<"Error:: set_emission_prob_of_model : "
				    <<"From attribute in "<<StateId<<"-th Sumover tag : "
				    <<SumOverElement->Attribute("From")
				    <<" is invalid "<<endl;
				check++;
			    }else{ // get from a parameter
				(*real)[StateId]->set_sum_over_FromParam(index);
			    }
			}else{ // get from a state
			    (*real)[StateId]->set_sum_over_FromState(index);
			}

			check+=(*real)[StateId]->set_sum_over_ThisPos(SumOverElement->Attribute("ThisPos"));
			
			int from = (*real)[StateId]->get_sum_over_FromState();
			if(from>0)
			{
			    NumSum = (*real)[from]->get_letters_to_read()
				- (*real)[StateId]->get_sum_over_ThisPos_dim(0);
			}else{
			    from = (*real)[StateId]->get_sum_over_FromParam();
			    NumSum = EP->get_FEP_dim(from)
				-(*real)[StateId]->get_sum_over_ThisPos_dim(0);
			}
			
			(*real)[StateId]->set_num_sum_over(NumSum);
			
			(*real)[StateId]->set_sum_over_FromPos(SumOverElement->Attribute("FromPos"));
			
			// compute the vector sum_over_pos
			k=0;
			j=0;
			
			int dim = (*real)[StateId]->get_sum_over_FromPos_dim(0);
			
			for(i=0;i<dim+NumSum;i++){
			    if(j<dim){
				if((*real)[StateId]->get_sum_over_FromPos(j)==i){
				    j++;
				}else{
				    sprintf(&SumPos[k],"%d",i);
				    k++;
				    if(k<2*(NumSum-1)){
					SumPos[k]=' ';
					k++;
				    }
				}
			    }else{
				sprintf(&SumPos[k],"%d",i);
				k++;
				if(k<2*(NumSum-1)){
				    SumPos[k]=' ';
				    k++;
				}
			    }
			}
			SumPos[k+1]='\0';
			(*real)[StateId]->set_sum_over_pos(SumPos);
		    
		    }else if(!cmp_nocase(EmissionProbElement->Attribute("SumOver"),"0")){
			(*real)[StateId]->set_emission_sum_over(false);
		    }else{
			cout<<"Error! new_models.cpp set_up_emission_prob_to_models: "<<endl;
			cout<<"combination of sum_buffer unregconized !"<<endl;
			cout<<"SumOver : "<<EmissionProbElement->Attribute("SumOver");
			check++;
		    }
		}else{		    
		    (*real)[StateId]->set_emission_sum_over(false);	    		
		}
		
		if(EmissionProbElement->Attribute("Product"))
		{
		    // read the product field
		    
		    if(!cmp_nocase(EmissionProbElement->Attribute("Product"),"1"))
		    {
			if((*real)[StateId]->get_emission_sum_over())
			{
			    (*real)[StateId]->set_emission_probs_expression("SumOver&Product");
			}else{
			    (*real)[StateId]->set_emission_probs_expression("Product");
			}
			(*real)[StateId]->set_emission_product(true);
			// set variables related to product site
			
			ProductElement = EmissionProbNode->FirstChildElement("Product");
			TiXmlElement* tempProductElement = ProductElement;

			int countproduct = 0;
			while(tempProductElement)
			{
			    countproduct++;
			    tempProductElement = tempProductElement->NextSiblingElement();
			}

			check+= (*real)[StateId]->set_number_of_products(countproduct);
		
			int product_index = 0;
			while((ProductElement)&&(!check))
			{
		
			    index = convert_State_id_to_int(MP,ProductElement->Attribute("From"));
			    if(index<0)
			    {
				index = convert_FEP_id_to_int(EP,ProductElement->Attribute("From"));
				if(index<0)
				{
				    cout<<"Error:: set_emission_prob_of_model : "
					<<"State attribute of the "<<StateId<<"-th"
					<<" Product tag : "
					<<ProductElement->Attribute("From")<<"is invalid"<<endl;
				    check++;
				}else{
				    (*real)[StateId]->set_product_FromParam(product_index,index);
				}
			    }else{				
				(*real)[StateId]->set_product_FromState(product_index,index);
			    }

			    (*real)[StateId]->set_product_ThisPos(product_index,
								  ProductElement->Attribute("ThisPos"));
			    
			    (*real)[StateId]->set_product_FromPos(product_index,
								  ProductElement->Attribute("FromPos"));

			    ProductElement = ProductElement->NextSiblingElement();
			    product_index++;
			}
			
		    }else if(!cmp_nocase(EmissionProbElement->Attribute("Product"),"0")){
			(*real)[StateId]->set_emission_product(false);
		    }else{
			cout<<"Error! new_models.cpp set_up_emission_prob_to_models: "<<endl;
			cout<<"combination of sum_buffer and product_buffer unregconized !"<<endl;
			cout<<"Product : "<<EmissionProbElement->Attribute("Product")<<endl;
			check++;
		    }
		}else{
		    (*real)[StateId]->set_emission_product(false);
		}
	    }
	}
	StateCount++;
	StateNode = StateNode->NextSibling();
	StateElement = StateElement->NextSiblingElement();
    }
    
    node=0;
    EmissionsNode=0;
    StateNode=0;
    EmissionProbNode=0;
    EmissionProbElement=0;
    StateElement=0;
    SumOverElement=0;
    ProductElement=0;

    if(ExpBuffer) delete [] ExpBuffer;
    ExpBuffer=NULL;
    if(SumPos) delete [] SumPos;
    SumPos = NULL;

    if(str_id) delete [] str_id;
    str_id = NULL;
    
    return(check);
}

int derive_emission_probs(Hmm* const real, 
			  model_parameters* const MP,
			  EmissionProb* const EP){
    int check=0;
    if ( real == NULL ){
	cout<<"ERROR: Derive_emission_probs: Paihmm real is NULL.\n"<<flush;
	check++;
	return check;
    }

    if( EP == NULL){
	cout<<"ERROR: Derive_emission_probs: EmissionProb EP is NULL.\n"<<flush;
	check++;
	return check;
    }

    if( MP == NULL){
	cout<<"ERROR: Derive_emission_probs: model_parameters MP is NULL.\n"<<flush;
	check++;
	return check;
    }
    
    int NumState,alphabet;
    int i,j,k,l,m,temp;
    int index1, index2;
    bool SumOver , Product;
    char* exp = new char[Max_word_length];

    NumState=MP->get_Number_of_States();
    Prob Tprob;
    Prob SumProb, ProductProb;
    for(i=0; i<NumState;i++)
    {
	if(check){
	    break;
	}
	strcpy(exp,(*real)[i]->get_emission_probs_expression());
	SumOver = (*real)[i]->get_emission_sum_over();
	Product = (*real)[i]->get_emission_product();
	alphabet = (*real)[i]->get_alphabet();
	if(!cmp_nocase(exp,"none"))
	{
	    // nothing has to be done
	    if((SumOver)||(Product)){
		cout<<"Error:: derive emission prob : state("<<i<<") "
		    <<"exp("<<exp<<") SumOver("<<SumOver<<") Product("<<Product<<")"<<endl;
		check++;
	    }
	    
	    if((*real)[i]->get_letters_to_read()!= 0)
	    {
		cout<<"Error:: derive_emission_prob : state("<<i<<") "
		    <<"exp("<<exp<<") while the state is supposed to read("
		    <<(*real)[i]->get_letters_to_read()<<") characters"<<endl;
		check++;
	    }
	    
	}else if(!cmp_nocase(exp,"SumOver")){
	    // SumOver
	    // only need to sum_over
	    Tprob = 0;
	    int SumOverFrom = (*real)[i]->get_sum_over_FromState();
	    bool state = true;
	    if(SumOverFrom<0)
	    {
		SumOverFrom = (*real)[i]->get_sum_over_FromParam();
		state = false;
	    }
	    int CurFixedDim = (*real)[i]->get_sum_over_ThisPos_dim(0);
	    int CurStateDim = (*real)[i]->get_number_of_dimensions_of_emission_probs();
	    
	    array<int> CurIndex;
	    CurIndex.SetNumberofDimensions(1);
	    CurIndex.SetDimension(0,CurStateDim);
	    
	    long NumFixed = static_cast<long>(pow(static_cast<float>(alphabet),static_cast<float>(CurFixedDim)));
	    for(j=0;j<NumFixed;j++){
		SumProb = 0;
		temp = j;
		// convert the index j to the "alphabet_based" representation
		for(k=0;k<CurFixedDim;k++){
		    CurIndex.SetElement((*real)[i]->get_sum_over_ThisPos(k),temp%alphabet);
		    temp = temp / alphabet;
		}
		SumProb=static_cast<Prob>(sum_over(real,EP,SumOverFrom,i,j,state,check));
		if(check){
		    cout<<"Error:: derive_emission_prob, sum_over : "<<endl;
		    cout<<"SumOverFrom : "<<SumOverFrom<<" CurState : "<<i<<" index : "<<j<<endl;
		    break;
		}
		(*real)[i]->set_emission_prob(CurIndex,SumProb);
		Tprob+=SumProb;
	    }
	    //check if the sum of those sum_over probs consistent
	    
	    if(abs(Tprob-1.)>Max_deviation){
		cout<<"Tprob : "<<Tprob<<endl;
		cout<<"ERROR: emission of state("<<i<<") != 1 (|sum_of_probs-1| = "
		    <<abs(Tprob-1.)<<"). Resetting all emission probs.\n";
		(*real)[i]->reset_emission_probs();
		check++;
		break;
	    }
	}else if(!cmp_nocase(exp,"SumOver&Product")){
	    // SumOver * Product
	    // need both sum over and multiplication
	    Tprob = 0;
	    
	    bool SumOverstate = true;
	    int SumOverFrom = (*real)[i]->get_sum_over_FromState();
	    if(SumOverFrom<0)
	    {
		SumOverFrom = (*real)[i]->get_sum_over_FromParam();
		SumOverstate = false;
	    }
	    int CurFixedDim = (*real)[i]->get_sum_over_FromPos_dim(0);
	    int CurStateDim = (*real)[i]->get_number_of_dimensions_of_emission_probs();
	  
	    array<int> CurIndex;
	    CurIndex.SetNumberofDimensions(1);
	    CurIndex.SetDimension(0,CurStateDim);
	       
	    long NumCur = static_cast<long>(pow(static_cast<float>(alphabet),static_cast<float>(CurStateDim)));

	    int number_of_products = (*real)[i]->get_number_of_products();
	   		
   	    for(j=0; j<NumCur; j++)
	    {
		// Set the index of product to the current emission vector, and the corr product
		temp = j;
		for(k=0; k<CurStateDim; k++)
		{
		    CurIndex.SetElement(k,temp%alphabet);
		    temp = temp/alphabet;
		    
		}
		ProductProb = 1.0;
		
		for(k=0; k<number_of_products; k++)
		{
		    bool Productstate = true;
		    int ProductFrom = (*real)[i]->get_product_FromState(k);
		    if(ProductFrom<0)
		    {
			ProductFrom = (*real)[i]->get_product_FromParam(k);
			Productstate = false;
		    }
		    int ProductDim = (*real)[i]->get_product_FromPos_dim(k,0);
		    		        
		    array<int> ProductCorrIndex;
		    ProductCorrIndex.SetNumberofDimensions(1);
		    ProductCorrIndex.SetDimension(0,ProductDim);
		    		
		    for(l=0; l<ProductDim; l++)
		    {
			int alph = CurIndex.GetElement((*real)[i]->get_product_ThisPos(k,l));
			ProductCorrIndex.SetElement((*real)[i]->get_product_FromPos(k,l),alph);
		    }
		    if(Productstate)
		    {
			ProductProb *= static_cast<Prob>((*real)[ProductFrom]->get_emission_prob(ProductCorrIndex));
		    }else{
			ProductProb *= static_cast<Prob>(EP->get_FEP_prob(ProductFrom,ProductCorrIndex));
		    }			    
		}
		
		// calculate the index;
		k = 0;
		for(l=0; l<CurFixedDim; l++)
		{
		    int alph = CurIndex.GetElement((*real)[i]->get_sum_over_ThisPos(l));
		    k = k + alph*static_cast<long>(pow(4,l));
		 
		}
		SumProb=static_cast<Prob>(sum_over(real,EP,SumOverFrom,i,k,SumOverstate,check));
		if(check){
		    cout<<"Error:: derive_emission_prob, sum_over : "<<endl;
		    cout<<"SumOverFrom : "<<SumOverFrom<<" CurState : "<<i<<" index : "<<j<<endl;
		    break;
		}
		// finished calculating the sumover prob, now the product prob
		
		// the combination of the two
		SumProb = static_cast<Prob>(SumProb*ProductProb);
	
		(*real)[i]->set_emission_prob(CurIndex,SumProb);
		Tprob+=SumProb;
	    }			
		    
	    if(abs(Tprob-1.)>Max_deviation){
		cout<<"Tprob : "<<Tprob<<endl;
		cout<<"ERROR: emission of state("<<i<<") != 1 (|sum_of_probs-1| = "
		    <<abs(Tprob-1.)<<"). Resetting all emission probs.\n";
		check++;
		(*real)[i]->reset_emission_probs();
		break;
	    }    
	    
	}else if(!cmp_nocase(exp,"Product"))
	{
	    // only product

	    Tprob = 0;

	    int CurStateDim = (*real)[i]->get_number_of_dimensions_of_emission_probs();
	    int number_of_products = (*real)[i]->get_number_of_products();

	    // check dimension
	    int product_dimension = 0;
	    for(j=0; j<number_of_products; j++)
	    {
		product_dimension+= (*real)[i]->get_product_FromPos_dim(j,0);
	    }
	    if(product_dimension!=CurStateDim)
	    {
		cout<<"ERROR: derive_emission_prob: state("<<i<<") product dimension doesn't match with state dimesion."<<endl;
		check++;
	    }

	    if(!check)
	    {
		
		long NumCurState  = static_cast<long>(pow(static_cast<float>(alphabet),static_cast<float>(CurStateDim)));

		array<int> CurIndex;
		CurIndex.SetNumberofDimensions(1);
		CurIndex.SetDimension(0,CurStateDim);	   

		for(j=0; j<NumCurState; j++)
		{
		    // get a particular curstate index
		    temp = j;
		    for(k=0; k<CurStateDim; k++)
		    {
			CurIndex.SetElement(k,temp%alphabet);
			temp = temp/alphabet;
			
		    }
		    // get the prob for the particular emission from the product states
		    ProductProb = 1.0;
			
		    for(k=0; k<number_of_products; k++)
		    {
			bool Productstate = true;
			int ProductFrom = (*real)[i]->get_product_FromState(k);
			if(ProductFrom<0)
			{
			    ProductFrom = (*real)[i]->get_product_FromParam(k);
			    Productstate = false;
			}
			int ProductDim = (*real)[i]->get_product_FromPos_dim(k,0);

			array<int> ProductCorrIndex;
			ProductCorrIndex.SetNumberofDimensions(1);
			ProductCorrIndex.SetDimension(0,ProductDim);
			
			for(l=0; l<ProductDim; l++)
			{
			    int alph = CurIndex.GetElement((*real)[i]->get_product_ThisPos(k,l));
			    ProductCorrIndex.SetElement((*real)[i]->get_product_FromPos(k,l),alph);
			}
			if(Productstate)
			{
			    ProductProb *= static_cast<Prob>((*real)[ProductFrom]->get_emission_prob(ProductCorrIndex));
			}else{
			    ProductProb *= static_cast<Prob>(EP->get_FEP_prob(ProductFrom,ProductCorrIndex));
			}			    
		    }
			
		    (*real)[i]->set_emission_prob(CurIndex,ProductProb);
		    
		    Tprob += ProductProb;
		    
		}

		if(abs(Tprob-1.)>Max_deviation){
		    cout<<"Tprob : "<<Tprob<<endl;
		    cout<<"ERROR: emission of state("<<i<<") : != 1 (|sum_of_probs-1| = "
			<<abs(Tprob-1.)<<"). Resetting all emission probs.\n";
		    check++;
		    (*real)[i]->reset_emission_probs();
		    break;
		}    

	    }	   

	}else if(convert_State_id_to_int(MP,exp)>=0){
	    // Get emission from other state
	    index1 = convert_State_id_to_int(MP,exp);
	    if(index1 < 0){
		cout<< "Error:: derive_emission_probs, in the case exp[0] == S "
		    <<"exp : "<<exp<<endl;
		cout<<"index1 : "<<index1<<endl;
		check++;
		break;
	    }

	    (*real)[i]->set_emission_prob((*real)[index1]->get_emission_prob());
 
	}else if(convert_FEP_id_to_int(EP,exp)>=0){
	    index1=convert_FEP_id_to_int(EP,exp);
	    if(index1 < 0){
		cout<<"Error:: derive_emission_probs, in the case exp[0] == E "
		    <<"exp : "<<exp<<endl;
		cout<<"index1 : "<<index1<<endl;
		check++;
		break;
	    }
	    // Get emission from EP
	    (*real)[i]->set_emission_prob(EP->get_FEP_prob(index1));
	}else{
	    cout<<"ERROR:: derive_emission_prob : exp invalid : "<<exp<<endl;
	    check++;
	    break;
	}
	
    }
    if(exp) delete[] exp;
    exp = NULL;

    return(check);
}

int derive_transition_probs_from_transition_parameters(Hmm* const real, 
						       model_parameters* const MP,
						       TransitionProb* const TP) {

    int check = 0;

    if (real == NULL) {
	cout << "ERROR: derive_transition_probs_from_transition_parameters: Hmm real is NULL.\n" << flush;
	check++;
	return check;
    }
    int NumState;
    int NumNextStates;
    int NextState;
    
    NumState = MP->get_Number_of_States();
    for (int i = 0; i<NumState; i++){
	NumNextStates = (*real)[i]->get_number_of_next_states();
	for (int j = 0; j<NumNextStates; j++){
	    NextState=(*real)[i]->get_number_of_next_state(j);
	    (*real)[i]->set_transition_prob(NextState, 
					    evaluate_transition_expression((*real)[i]->get_transition_probs_expressions(NextState),TP,check));
	}
    }
    return(check);
}

int get_hmm(const char* filename,
		Hmm* const real,
		model_parameters* const MP,
		TransitionProb* const TP,
		EmissionProb * const EP)
{
	
    int check=0;

    int i = 0;
    int NumState = MP->get_Number_of_States();
    if (check==0) {
       
	// initialise pairhmm
	
	for (i=0; i<NumState; i++) {
	    (*real)[i]->completely_reset_transition_probs();
	    (*real)[i]->completely_reset_emission_probs();
	    (*real)[i]->completely_reset_transition_probs_expressions();
	    (*real)[i]->completely_reset_emission_probs_expression();
	    (*real)[i]->completely_reset_transition_scores();
	    (*real)[i]->completely_reset_emission_scores();
	    (*real)[i]->reset_viterbi_strip();
	    (*real)[i]->reset_viterbi_rectangle();
	    (*real)[i]->reset_viterbi_scores();
	}
	
	// set pair

	check += real->set_pair(MP->is_PairHMM()); 

	// set up topology of model and define states
	
	check += set_up_topology_of_model(filename,real, MP);
	
	if (check) {
	    cout << "ERROR: get_hmm: set_up_topology_of_model\n" << flush;
	}
	
	// read in transition expression (and priors)		
	if (!check){
	    check+=set_up_transition_probs_of_model(real,filename,MP); //read in the expression
	    if (check) {
		cout << "ERROR: get_hmm: set_up_transition_probs_of_model.\n" << flush;
	    }
	}

	if(TP){
	    // derive transition probs from transition parameters and set them in pairhmm
	    if (!check) {			
		check += derive_transition_probs_from_transition_parameters(real,MP,TP); 
		if (check) {
		    cout << "ERROR: get_hmm: derive_transition_probs_from_transition_parameters.\n" << flush;
		}
	    }
	}

	// read in emission probs and set them in pairhmm
	if (!check ) {
	    check += set_up_emission_probs_of_model(real,filename,MP,EP);
	    if (check) {
		cout << "ERROR: get_hmm: set_up_emission_probs_of_model.\n" << flush;
	    }
	}
	
	// derive emission probs 
	if(!check){
	    check+= derive_emission_probs(real, MP, EP);
	    if(check){
		cout<<"ERROR: get_hmm: derive_emission_probs.\n"<<flush;
	    }
	}

	// derive scores from probs
	
	
	if(!check){
	    check+= real->calculate_scores_from_probs();
	    if(check){
		cout<<"ERROR: get_hmm: Hmm::calculate_scores_from_probs.\n"<<flush;
	    }
	}
	
	// check consistency of pairhmm
	
	if (!check) {

	    check += real->check_consistency_of_probs();
	    
	    if (check) {
		cout << "ERROR: get_hmm: Hmm::check_consistency_of_probs().\n" << flush;
	    }
	}
	    
    }
    
    // initialise pairhmm, if checks failed
    
    if (check) {
	for (i=0; i<NumState; i++) {
	    (*real)[i]->completely_reset_transition_probs();
	    (*real)[i]->completely_reset_emission_probs();
	    (*real)[i]->completely_reset_transition_probs_expressions();
	    (*real)[i]->completely_reset_emission_probs_expression();
	    (*real)[i]->completely_reset_transition_scores();
	    (*real)[i]->completely_reset_emission_scores();
	    (*real)[i]->reset_viterbi_strip();
	    (*real)[i]->reset_viterbi_rectangle();
	    (*real)[i]->reset_viterbi_scores();	}
    }
    return(check);
}

char* postfix(const char* exp, int& check)
{
    char token, topToken;
    Stack<char> opStack;
    const char BLANK = ' ';
    if(!exp){
	cout<<"Error : postfix, input string exp is NULL "<<endl;
	check++;
	return NULL;
    }
    int length = strlen(exp);
    char* postfixExp = new char[Max_word_length];
    int pos = 0;
    postfixExp[pos] = '\0';
    for(int i=0; i<length; i++)
    {
	token=exp[i];
	switch(token)
	{
	    case ' ': break; //skip blank 
	    case '(': opStack.push(token);
		break;
	    case ')':
		for(;;){
		    assert(!opStack.empty());
		    topToken = opStack.top();
		    opStack.pop();
		    if(topToken == '(') 
			break;
		    postfixExp[pos]=BLANK;
		    pos++;
		    postfixExp[pos]=topToken;
		    pos++;
		}
		break;
	    case '+' : case '-' : 
	    case '*' : case '/' : case '%':
		for (;;){
		    if(opStack.empty()||opStack.top()=='('||
		       (token == '*'|| token=='/' || token=='%')&&
		       (opStack.top()=='+'|| opStack.top()=='-'))
		    {
			opStack.push(token);
			break;
		    }
		    else{
			topToken = opStack.top();
			opStack.pop();
			postfixExp[pos]=BLANK;
			pos++;
			postfixExp[pos]=topToken;
			pos++;
		    }
		}
		break;
	    default:
		postfixExp[pos]=BLANK;
		pos++;
		postfixExp[pos]=token;
		pos++;
		for(;;){
		    if((exp[i+1]=='(')||(exp[i+1]==')')||(exp[i+1]=='+')
		       ||(exp[i+1]=='-')||(exp[i+1]=='*')||(exp[i+1]=='/')
		       ||(i>=length-1)){
			break;
		    }
		    i++;
		    token = exp[i];
		    postfixExp[pos]=token;
		    pos++;
		}
		break;
	}
    }
    //pop remaining operators on the stack
    for(;;){
	if(opStack.empty()) 
	    break;
	topToken = opStack.top();
	opStack.pop();
	if(topToken != '(')
	{
	    postfixExp[pos] = BLANK;
	    pos++;
	    postfixExp[pos] = topToken;
	    pos++;
	}else{
	    check++;
	    cout<<"Error: postfix, in infix expression"<<endl;
	    break;
	}
    }

    postfixExp[pos]='\0';
    length = strlen(postfixExp);
    char* RpostfixExp = new char[length+1];
    strcpy(RpostfixExp,postfixExp);
    // free memory
    if(postfixExp) delete[] postfixExp;
    postfixExp = NULL;

    return RpostfixExp;
}
