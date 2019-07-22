 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#include<iostream>
#include<stdlib.h>
#include"model_parameters.h"

using namespace std;

//For model_parameters class

//Constructors
model_parameters :: model_parameters()
{
    pair = false;
    Model_Name = NULL;
    SpecialEmit = false;

    Number_of_Special_Emissions = 0;

    Total_Number_of_Annotation_Labels = 0;
    Annotation_Label = NULL;

    Alphabet.size = 0;
    Alphabet.name = NULL;

    Number_of_States = 0 ;
    State_Name = NULL;

}

model_parameters :: model_parameters(char* const m_type,
			       const bool p,
			       const bool se,
			       const int n_of_special_emissions,
			       const int total_n_of_a_labels,
			       Label* const a_label,
			       const int total_n_of_m_elements,
			       char** const m_element,
			       const Letter alph,
			       const int n_of_states,
			       char** const s_name)
{
    int i = 0;
    int j = 0;
    int check = 0;

    Model_Name = Nullstrcpy(m_type,check);
    if(check!=0){
	cout<<"ERROR: model_parameters Constructor : m_type : "<<m_type<<endl;
	check++;
    }

    pair = p;
    SpecialEmit = se;

    Number_of_Special_Emissions = n_of_special_emissions;

    Total_Number_of_Annotation_Labels = total_n_of_a_labels;
    Annotation_Label = new Label[Total_Number_of_Annotation_Labels];
    for(i=0;i<Total_Number_of_Annotation_Labels;i++){
	if(a_label[i].size<=0){
	    cout<<"ERROR: model_parameters Constructor : a_label["<<i<<"].size : "<<a_label[i].size<<endl;
	    check++;
	    break;
	}
	Annotation_Label[i].size = a_label[i].size;
	Annotation_Label[i].setname = Nullstrcpy(a_label[i].setname,check);
	if(check!=0){
	    cout<<"ERROR: model_parameters Constructor : a_label["<<i<<"].setname : "<<a_label[i].setname<<endl;
	    check++;
	    break;
	}
	for(j=0;j<Annotation_Label[i].size;j++){
	    Annotation_Label[i].name[j]=Nullstrcpy(a_label[i].name[j],check);
	    if(check!=0){
		cout<<"ERROR: model_parameters Constructor : a_label["<<i<<"].name["<<j<<"] : "<<a_label[i].name[j]<<endl;
		check++;
		break;
	    }
	}
    }
    
    if(alph.size<=0){
	cout<<"ERROR: model_parameters Constructor : alph.size : "<<alph.size<<endl;
	check++;
    }
    if(!check){
	Alphabet.size = alph.size;
    }

    for(i=0;i<Alphabet.size;i++){
	Alphabet.name[i]=alph.name[i];
    }
    
    Number_of_States = n_of_states;
    State_Name = new char*[Number_of_States];
    for(i=0;i<Number_of_States;i++){
	State_Name[i]=Nullstrcpy(s_name[i],check);
	if(check!=0){
	    cout<<"ERROR: model_parameters Constructor: s_name[i] : "<<s_name[i]<<endl;
	    break;
	}
    }
    if(check){
	this->~model_parameters();
    }
}

model_parameters :: model_parameters(char* filename, int& check)
{
    
    //get model elements
    check+= set_model_element_from_xml(filename);

    //get alphabet
    check+= set_alphabet_from_xml(filename);

    //get annotation label
    check+= set_annotation_labels_from_xml(filename);

    check+= set_state_name_from_xml(filename);

    if(check){
	this->~model_parameters();
    }
}

model_parameters:: ~model_parameters()
{
    int i = 0;
    int j = 0;
    
    if(Model_Name) delete[] Model_Name;
    Model_Name;
    pair = false;
    SpecialEmit = false;

    Number_of_Special_Emissions = 0;
 
    for(i=0; i<Total_Number_of_Annotation_Labels; i++)
    {
	if(Annotation_Label[i].setname) delete[] Annotation_Label[i].setname;
	Annotation_Label[i].setname = NULL;
	for(j=0; j<Annotation_Label[i].size; j++){
	    if(Annotation_Label[i].name[j]) delete[] Annotation_Label[i].name[j];
	    Annotation_Label[i].name[j] = NULL;
	}
	Annotation_Label[i].size = 0;
    }
    if(Annotation_Label) delete[] Annotation_Label;
    Annotation_Label = NULL;

    Total_Number_of_Annotation_Labels = 0;
    
    if(Alphabet.name) delete[] Alphabet.name;
    Alphabet.name = NULL;
    Alphabet.size = 0;

    for(i = 0; i < Number_of_States ; i++){
	if(State_Name[i]) delete[] State_Name[i];
	State_Name[i] = NULL;
    }
    if(State_Name) delete[] State_Name;
    State_Name = NULL;

    Number_of_States = 0;

}

// Private functions

int model_parameters::set_model_element_from_xml(char* const filename)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: set_model_element_from_xml : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * RootNode = 0;
    TiXmlNode * TmpNode = 0;
    TiXmlElement* Element = 0;
    int count = -1;
    // create memory for Model_Element

    RootNode = doc.FirstChild("HMMConverter");
    assert(RootNode);
    RootNode = RootNode->FirstChild("model");
    assert(RootNode);

    // get Model Type
    Element = RootNode->FirstChildElement("Model_Type");
    assert(Element);
    
    if(!Element->Attribute("name")){
	cout<<"Error : set_model_element_from_xml : "
	    <<"type attribute in Model_Type tag is NULL"<<endl;
	check++;
    }else{
	Model_Name=Nullstrcpy(Element->Attribute("name"),check);
	if(check){
	    cout<<"Error : set_model_element_from_xml : "
		<<"error in Nullstrcpy for type attribute in "
		<<"in Model_Type tag."<<endl;
	}		
    }

    if(!Element->Attribute("pair")){
	pair = false;                                 // default as HMM
    }else{
	if(!cmp_nocase(Element->Attribute("pair"),"1"))
	{
	    pair = true;
	}else if(!cmp_nocase(Element->Attribute("pair"),"0"))
	{
	    pair = false;
	}else
	{
	    cout<<"ERROR: set_model_element_from_xml : "
		<<"pair attribute in Model_Type tag "
		<<Element->Attribute("pair")
		<<" not equal to either  y or n "<<endl;
	    check++;
	}
    }

    if(!Element->Attribute("SpecialEmission")){
	SpecialEmit = false;
	Number_of_Special_Emissions = 0;
    }else{	    
	if(!cmp_nocase(Element->Attribute("SpecialEmission"),"1"))
	{
	    SpecialEmit = true;
	}else if(!cmp_nocase(Element->Attribute("SpecialEmission"),"0"))
	{
	    SpecialEmit = false;
	    Number_of_Special_Emissions = 0;
	}else{
	    cout<<"ERROR: set_model_element_from_xml : "
		<<"SpecialEmission attribute in Model_Type tag "
		<<Element->Attribute("SpecialEmission")
		<<" not equal to either y or n "<<endl;
	    check++;	   
	}	
    }

    RootNode = 0;
    TmpNode = 0;
    Element = 0;
  
    return check;
}

int model_parameters::set_alphabet_from_xml(char* const filename)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: set_alphabet_from_xml : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * RootNode = 0;
    TiXmlNode * ANode = 0; // Node of Alphabets
    TiXmlElement* AElement = 0; // Element of Alphabets

    int length=0;
    int Asize = 0;
    int Acount = -1;
    int Acases = -1;
    
    RootNode = doc.FirstChild("HMMConverter");
    assert(RootNode);
    RootNode = RootNode->FirstChild("model");
    assert(RootNode);

    ANode = RootNode->FirstChild("Alphabets");
    assert(ANode);
    
    AElement = RootNode->FirstChildElement("Alphabets");
    assert(AElement);   

    if(AElement->Attribute("cases"))
    {
	AElement->Attribute("cases",&Acases);
	if((Acases<0)||(Acases>1)){
	    cout<<"Error: set_alphabets_from_xml, Alphabet.cases : "<<Acases
		<<" out of range [0,1] "<<endl;
	    check++;
	}
    }else{
	Acases = 0;
    }
    Alphabet.cases = Acases;
    
    if(!AElement->Attribute("set"))
    {
	cout<<"Error : set_alphabet_from_xml : no valid alphabet is set "<<endl;
	check++;
    }

    if(!check)
    {

	if(strlen(AElement->Attribute("set"))<=0)
	{
	    cout<<"Error : set_alphabet_from_xml : no valid alphabet is set : "
		<<AElement->Attribute("set")<<endl;
	    check++;
	}else{
		Alphabet.size = strlen(AElement->Attribute("set"));
	}
			
	if(!check)
	{
	    Alphabet.name = new char[Alphabet.size];
	    for(int i =0; i<Alphabet.size; i++)
	    {
		Alphabet.name[i] = AElement->Attribute("set")[i];
	    }
	}
    }

    RootNode = 0;
    ANode = 0;
    AElement = 0;
    
    if(check)
    {

	if(Alphabet.name) delete[] Alphabet.name;
	Alphabet.name = NULL;
	Alphabet.size = 0;

    }
    return check;
}

int model_parameters::set_annotation_labels_from_xml(char* const filename)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: set_annotation_labels_from_xml : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * RootNode = 0;
    TiXmlNode * TmpNode1 = 0; 
    TiXmlNode * TmpNode2 = 0; 
    TiXmlNode * TmpNode3 = 0; 
    TiXmlElement* Element1 = 0;
    TiXmlElement* Element2 = 0;
    int length=0;
    int size = 0;
    int count1 = -1;
    int count2 = -1;
  
    RootNode = doc.FirstChild("HMMConverter");
    assert(RootNode);
    RootNode = RootNode->FirstChild("model");
    assert(RootNode);

    TmpNode1 = RootNode->FirstChild("Annotation_Labels");
    if(!TmpNode1)
    {
	Total_Number_of_Annotation_Labels = 0;
	Annotation_Label = NULL;
	return check;
    }

    // get Total_Number_of_Annotation_Labels
    TiXmlElement* CountElement = 0;
    size = 0;
    CountElement = TmpNode1->FirstChildElement("Annotation_Label");
    while(CountElement)
    {
	size++;
	CountElement = CountElement->NextSiblingElement();
    }

    if(size<=0)
    {
	cout<<"Error: set_annotation_labels_from_xml, Total_Number_of_Annotation_Labels : "<<size<<endl;
	check++;
	return check;
    }


    Total_Number_of_Annotation_Labels = size;

    Annotation_Label = new Label[Total_Number_of_Annotation_Labels];    

    TmpNode2 = TmpNode1->FirstChild("Annotation_Label");
    assert(TmpNode2);

    Element1 = TmpNode1->FirstChildElement("Annotation_Label");

    while((Element1)&&(!check))
    {
	assert(Element1);

	if(count1 >= Total_Number_of_Annotation_Labels ){
	    cout<<"Error : set_annotation_labels_from_xml : number of elements : "<<count1;
	    cout<<"Total_Number_of_Annotation_Labels : "<<Total_Number_of_Annotation_Labels<<endl;
	    check++;
	    break;
	}

	// get setname
	if(!Element1->Attribute("name")){
	    cout<<"Error : set_annotation_labels_from_xml : "
		<<"name attribute in Annotation_Label tag is NULL"<<endl;
	    check++;
	    break;
	}else{
	    count1++;
	    Annotation_Label[count1].setname=Nullstrcpy(Element1->Attribute("name"),check);
	    if(check){
		cout<<"Error : set_annotation_labels_from_xml : "
		    <<"error in Nullstrcpy for name attribute in Annotation_Label tag"<<endl;
		break;
	    }	    
	}

	// get score
	if(Element1->Attribute("score"))
	{
	    if(!cmp_nocase(Element1->Attribute("score"),"1"))
	    {
		Annotation_Label[count1].score = true;
	    }else{
		Annotation_Label[count1].score = false;
	    }
	}else{
	    Annotation_Label[count1].score = false;
	}

	//get size
	CountElement = TmpNode2->FirstChildElement("label");
	size = 0;
	while(CountElement)
	{
	    size++;
	    CountElement = CountElement->NextSiblingElement();
	}

	Annotation_Label[count1].size=size;
	Annotation_Label[count1].name = new char* [size];

	Element2=TmpNode2->FirstChildElement("label");

	count2=-1;
	// get annotation_label
	while((Element2)&&(!check))
	{
	    assert(Element2);
	    if (count2 >= Annotation_Label[count1].size ){
		cout<<"Error : set_annotation_labels_from_xml : number of elements("<<count2
		    <<") exceed Annotation_Label["<<count1<<"].size("<<Annotation_Label[count1].size
		    <<")."<<endl;
		check++;
		break;
	    }
	    if (!Element2->Attribute("name")){
		cout<<"Error : set_annotation_labels_from_xml : Element2->Attribute is NULL "<<endl;
		check++;
		break;
	    }else{		
		count2++;
		Annotation_Label[count1].name[count2]=Nullstrcpy(Element2->Attribute("name"),check);
		if(check){
		    cout<<"Error : set_state_label_from_xml "
			<<"error in Nullstrcpu for name attribute in Annotation_Label tag "<<endl;
		    break;		    
		}
	    }
	    Element2 = Element2->NextSiblingElement();
	}
	TmpNode2 = TmpNode2->NextSibling();
	Element1 = Element1->NextSiblingElement();
    }

    RootNode = 0;
    TmpNode1 = 0;
    TmpNode2 = 0;
    TmpNode3 = 0;
    Element1 = 0;
    Element2 = 0;
    if(check){
	for(int i =0; i<=count1;i++){
	    if(Annotation_Label[i].setname) delete[] Annotation_Label[i].setname;
	    Annotation_Label[i].setname = NULL;
	    for(int j =0; j<Annotation_Label[i].size; j++){
		if(i!=count1){
		    if(Annotation_Label[i].name[j]) delete[] Annotation_Label[i].name[j];
		    Annotation_Label[i].name[j] = NULL;
		}else{
		    if(j<=count2){
			if(Annotation_Label[i].name[j]) delete[] Annotation_Label[i].name[j];
			Annotation_Label[i].name[j] = NULL;
		    }
		}
	    }
	    Annotation_Label[i].size = 0;
	}
	if(Annotation_Label) delete[] Annotation_Label;
	Annotation_Label = NULL;
    }
    return check;
}

int model_parameters::set_state_name_from_xml(char* const filename)
{
    TiXmlDocument doc(filename);
    bool loadOkay = doc.LoadFile();
    int check = 0;
    if(!loadOkay){
	cout << "ERROR: set_state_name_from_xml : cannot open file " << filename << ".\n" << flush;
	check++;
	return check;
    }
    TiXmlNode * RootNode = 0;
    TiXmlNode * TmpNode = 0;
    TiXmlElement* Element = 0;
    TiXmlElement* CountElement = 0;
    int count = -1;
 
    // create memory for State_Name

    RootNode = doc.FirstChild("HMMConverter");
    assert(RootNode);
    RootNode = RootNode->FirstChild("model");
    assert(RootNode);

    TmpNode = RootNode->FirstChild("States");
    assert(TmpNode);

    Element = RootNode->FirstChildElement("States");
    assert(Element);
  
    //get Number_of_States
    CountElement = TmpNode->FirstChildElement("State");
    Number_of_States = 0;
    while(CountElement)
    {
	Number_of_States++;
	CountElement = CountElement->NextSiblingElement();
    }
    
    if(Number_of_States<=0){
	cout<<"Error : set_state_name_from_xml, Number_of_States : "<<Number_of_States<<endl;
	check++;
	return check;
    }
    
    State_Name = new char*[Number_of_States];

    Element = TmpNode->FirstChildElement("State");
    
    TmpNode = TmpNode->FirstChild("State");
    while(Element){
	assert(Element);
	if(count >= Number_of_States ){
	    cout<<"Error : set_state_name_from_xml : number of states : "<<count;
	    cout<<" Number_of_States : "<<Number_of_States<<endl;
	    check++;
	    break;
	}
	if(!Element->Attribute("name")){
	    cout<<"Error : set_state_name_from_xml : Element->Attribute(name) is NULL"<<endl;
	    check++;
	    break;
	}else{
	    count++;
	    State_Name[count]=Nullstrcpy(Element->Attribute("name"),check);
	    if(check){
		cout<<"Error : set_state_name_from_xml : "
		    <<"error in Nullstrcpy for name attribute in State tag."<<endl;
		break;	    
     	    }
	   
	}
	Element = Element->NextSiblingElement();
    }

    RootNode = 0;
    TmpNode = 0;
    Element = 0;
    if(check){
	for(int i = 0; i <=count; i++){
	    if(State_Name[i]) delete[] State_Name[i];
	    State_Name[i] = NULL;
	}
	if(State_Name) delete[] State_Name;
	State_Name = NULL;
    }
    return check;  
}

// accessors

char* model_parameters::get_Model_Name() const
{
    return Model_Name;
}

bool model_parameters::is_PairHMM() const
{
    return pair;
}

bool model_parameters::is_SpecialEmit() const
{
    return SpecialEmit;
}

int model_parameters::get_Number_of_Special_Emissions() const
{
    return Number_of_Special_Emissions;
}

int model_parameters::get_Max_n_of_Child_State() const
{
    return Max_n_of_Child_State;
}

int model_parameters::get_Total_Number_of_Annotation_Labels() const
{
    return Total_Number_of_Annotation_Labels;
}

int model_parameters:: get_Annotation_Label_size(int i) const
{
    if((i<0)||(i>=Total_Number_of_Annotation_Labels))
	return -1;
    return Annotation_Label[i].size;
}

char* model_parameters:: get_Annotation_Label_setname(int i) const
{
    if((i<0)||(i>=Total_Number_of_Annotation_Labels))
	return NULL;
    return Annotation_Label[i].setname;
}

char* model_parameters:: get_Annotation_Label_name(int i, int j) const
{
    if((i<0)||(i>=Total_Number_of_Annotation_Labels))
	return NULL;
    if((j<0)||(j>=Annotation_Label[i].size))
	return NULL;
    return Annotation_Label[i].name[j];
}

char** model_parameters:: get_Annotation_Label_name(int i) const
{
    if((i<0)||(i>=Total_Number_of_Annotation_Labels))
	return NULL;
    return Annotation_Label[i].name;
}

Label* model_parameters:: get_Annotation_Label() const
{
	return Annotation_Label;
}

bool model_parameters::get_score_of_Annotation_Label(const int i) const
{
    if((i<0)||(i>=Total_Number_of_Annotation_Labels))
	return false;
    return Annotation_Label[i].score;
}

bool model_parameters::get_Alphabet_cases() const
{
    return Alphabet.cases;
}

int model_parameters::get_Alphabet_size() const
{
    return Alphabet.size;
}

char model_parameters::get_Alphabet_name(const int i) const
{
    if((i<0)|(i>=Alphabet.size))
	return static_cast<char>(0);
    return Alphabet.name[i];
}

char* model_parameters:: get_Alphabet_name() const
{
    return Alphabet.name;
}

Letter model_parameters:: get_Alphabet() const
{
    return Alphabet;
}

int model_parameters::get_Number_of_States() const
{
    return Number_of_States;
}

char* model_parameters:: get_State_Name(const int i) const
{
    return State_Name[i];
}

char** model_parameters:: get_State_Name() const
{
    return State_Name;
}
// mutator

int model_parameters::set_Model_Name(char* const name)
{
    int check = 0;
    Model_Name = Nullstrcpy(name,check);
    if(check!=0){
	cout<<"ERROR: model_parameters::set_Model_Name, can not copy input type("
	    <<name<<").\n"<<flush;
    }
    return check;
}

int model_parameters::set_pair(const bool p)
{
    pair = p;
    return 0;
}

int model_parameters::set_SpecialEmit(const bool se)
{
    SpecialEmit = se;
    return 0;
}

int model_parameters::set_Number_of_Special_Emissions(const int n_of_special_emissions)
{
    int check=0;
    if(n_of_special_emissions<0)
    {
	check++;
    }
    Number_of_Special_Emissions = n_of_special_emissions;
    return check;
}

int model_parameters::set_Max_n_of_Child_State(const int max_n_of_child_state)
{
    int check=0;
    if(max_n_of_child_state<0)
    {
	check++;
    }
    Max_n_of_Child_State = max_n_of_child_state;
    return check;
}

int model_parameters::set_Total_Number_of_Annotation_Labels(const int number){
    Total_Number_of_Annotation_Labels = number; 
    return 0;
}

int model_parameters::set_Annotation_Label_size(const int i, const int s){
    if((i>=Total_Number_of_Annotation_Labels)||(i<0)){
	cout<<"Error! set_Annotation_Label_size:: "<<endl;
	cout<<"i "<<i<<" exceed Total_Number_of_Annotation_Labels : "<<Total_Number_of_Annotation_Labels<<endl;
	return 1;
    }
    Annotation_Label[i].size=s;
    return 0;
}

int model_parameters:: set_Annotation_Label_setname(const int i, char* const tn){

    int check = 0;

    if((i>=Total_Number_of_Annotation_Labels)||(i<0)){
	cout<<"Error! set_Annotation_Label_setname:: "<<endl;
	cout<<"i "<<i<<" exceed Total_Number_of_Annotation_Labels : "<<Total_Number_of_Annotation_Labels<<endl;
	check++;
    }
    
    Annotation_Label[i].setname= Nullstrcpy(tn,check);

    if(check){
	cout<<"Error! set_Annotation_Label_setname:: "<<endl;
	cout<<"Input string not valid : "<<tn<<endl;
    }
 
    return check;
}

int model_parameters:: set_Annotation_Label_name(const int i, const int j, char* const n){
    int check = 0;

    if((i>=Total_Number_of_Annotation_Labels)||(i<0)){
	cout<<"Error! set_Annotation_Label_name:: "<<endl;
	cout<<"i "<<i<<" exceed Total_Number_of_Annotation_Labels : "<<Total_Number_of_Annotation_Labels<<endl;
	check++;
    }
    
    if((j>=Annotation_Label[i].size)||(j<0)){
	cout<<"Error! set_Annotation_Label_name:: "<<endl;
	cout<<"j "<<j<<" exceed Annotation_Label[i].size "<<Annotation_Label[i].size<<endl;
	check++;
    }

    Annotation_Label[i].name[j]=Nullstrcpy(n,check);

    if(check){
	cout<<"Error! set_Annotation_Label_name:: "<<endl;
	cout<<"Input String invalid : "<<n<<endl;
    }

    return check;
}

int model_parameters::set_score_of_Annotation_Label(const int i, const bool s)
{
    int check = 0;

    if((i>=Total_Number_of_Annotation_Labels)||(i<0)){
	cout<<"Error! set_score_of_Annotation_Label:: "<<endl;
	cout<<"i "<<i<<" exceed Total_Number_of_Annotation_Labels : "<<Total_Number_of_Annotation_Labels<<endl;
	check++;
    }

    if(!check)
    {
	Annotation_Label[i].score = s;
    }

    return check;
}

int model_parameters:: set_Alphabet_cases(const bool c){
    Alphabet.cases = c;
    return 0;
}

int model_parameters:: set_Alphabet_size(const int s){
    Alphabet.size = s;
    return 0;
}

int model_parameters:: set_Alphabet_name(const int i, const char n){
    int check = 0;

    if((i>=Alphabet.size)||(i<0)){
	cout<<"Error! set_Alphabet_name:: "<<endl;
	cout<<"i "<<i<<" exceed Alphabet.size "<<Alphabet.size<<endl;
	check++;
    }
    
    if(!check){
	Alphabet.name[i]=n;
    }
           
    return check;

}

int model_parameters:: set_Number_of_States(const int number){
    Number_of_States = number;
    return 0;
}

int model_parameters:: set_State_Name(const int i, char* const n)
{
    int check = 0;
    if ((i>=Number_of_States)||(i<0)){
	cout<<"Error! set_State_Name:: "<<endl;
	cout<<"i "<<i<<" exceed Number_of_States : "<<Number_of_States<<endl;
	return 1;
    }

    State_Name[i]=Nullstrcpy(n,check);
    
    if(check){
       cout<<"Error! set_State_Name:: "<<endl;
       cout<<"Input string not valid : "<<n<<endl;
    }

    return check;
  
}

void model_parameters::print(std::ostream &o) const
{
    int i = 0;
    int j = 0;

    o<<"Model Name : "<<Model_Name<<endl;
    if(pair){
	o<<"This is a PairHMM "<<endl;
    }
    o<<"Special Emission : ";
    if(SpecialEmit){
	o<<"on";
    }else{
	o<<"off";
    }
    o<<endl;
    o<<"Number of Special Emission : "<<Number_of_Special_Emissions<<endl;

    o<<"Number of Annotation Labels : "<<Total_Number_of_Annotation_Labels<<endl;
    for(i=0;i<Total_Number_of_Annotation_Labels;i++){
	o<<"Annotation_Label["<<i<<"] : "<<Annotation_Label[i].setname<<endl;
	for(j=0;j<Annotation_Label[i].size;j++){
	    o<<"Annotation_Label["<<i<<"].name["<<j<<"] : "<<Annotation_Label[i].name[j]<<endl;
	}
	o<<"-----------------------------------------------------"<<endl;
    }
    o<<endl;
    o<<"-----------------------------------------------------"<<endl;
    o<<endl;

    o<<"Alphabet.size : "<<Alphabet.size<<endl;
    o<<"Alphabet.cases : "<<Alphabet.cases<<endl;
    for(i=0;i<Alphabet.size;i++){
	o<<"Alphabet.name["<<i<<"] : "<<Alphabet.name[i]<<endl;
    }
    o<<"-----------------------------------------------------"<<endl;
    o<<endl;
    
    o<<"Number of States : "<<Number_of_States<<endl;
    for(i=0;i<Number_of_States;i++){
	o<<"State["<<i<<"] : "<<State_Name[i]<<endl;
    }
    o<<endl;

    return;
}

int convert_state_name_to_int(model_parameters* const MP, char* const buffer)
{
    int size = MP->get_Number_of_States();
    if(!buffer){
	return -1;
    }
    for(int i = 0; i<size; i++){
	if(!strcmp(buffer, MP->get_State_Name(i))){
	    return i;
	}
    }
    return -1;
}

int convert_alphabet_to_int(model_parameters* const MP, const char buffer)
{
    int size = MP->get_Alphabet_size();
    bool cases = MP->get_Alphabet_cases();
    if(!buffer){
	return -1;
    }
    if(!cases){
	for(int i = 0; i<size; i++){
	    if(toupper(buffer)==toupper( MP->get_Alphabet_name(i))){
		return i;
	    }
	}
    }else{
	for(int i = 0; i<size; i++){
	    if(buffer ==  MP->get_Alphabet_name(i)){
		return i;
	    }
	}
    }
    return -1;
}

int convert_State_id_to_int(model_parameters* const MP, const char*  buffer)
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
    if(strcmp(tokens[0],"S")){
	check++;
    }
    int state_num = atoi(tokens[1]);
 
    if(!check){
	return state_num;
    }else{
	return -1;
    }
}
 
int convert_Annotation_Label_id_to_int(model_parameters* const MP, const char* buffer, int type)
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
    int num_of_state_label = MP->get_Total_Number_of_Annotation_Labels();
    if((type<0)||(type>=num_of_state_label)){
	return -1;
    }
    if(strcmp(tokens[0],MP->get_Annotation_Label_setname(type))){

	check++;	
    }

    int state_label_id = atoi(tokens[1]);
    if((state_label_id<0)||(state_label_id>MP->get_Annotation_Label_size(type)))
    {
	check++;
    }
    if(!check){
	return state_label_id;
    }else{
	return -1;
    }
} 

int convert_Annotation_Label_set_id_to_int(model_parameters* const MP, const char* buffer)
{
    char** tokens;
    int num_of_tokens;
    int check = 0;
    check+=break_id_into_tokens(buffer,&tokens,num_of_tokens);
    if(check)
    {
	return -1;
    }
    if((num_of_tokens>2)||(num_of_tokens<=1))
    {
	check++;
	return -1;
    }
    if(strcmp(tokens[0],"SL"))
    {
	check++;
    }
    int state_label_num = atoi(tokens[1]);
    if((state_label_num<0)||(state_label_num > MP->get_Total_Number_of_Annotation_Labels()))
    {

	check++;
    }
    if(!check){
	return state_label_num;
    }else{
	return -1;
    }
}
