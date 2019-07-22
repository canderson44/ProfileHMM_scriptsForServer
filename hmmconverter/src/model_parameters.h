 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
 */

#ifndef _model_parameters_h
#define _model_parameters_h

#include <iostream>
#include "tinyxml.h"
#include "define.h"

struct Label{
    int size;
    char* setname;
    char** name;
    bool score;
};

struct Letter{
    int size;
    bool cases;
    char* name;
};


class model_parameters{
 private:

    // Model Type

    char* Model_Name;
    bool pair;
    bool SpecialEmit;

    int Number_of_Special_Emissions;
    int Max_n_of_Child_State;

    // Annotation Labels

    int Total_Number_of_Annotation_Labels;
    Label* Annotation_Label;

    // Alphabet

    Letter Alphabet;

    // States

    int Number_of_States;
    char** State_Name;

    // Private functions
    
    int set_model_element_from_xml(char* const filename);
    int set_alphabet_from_xml(char* const filename);
    int set_annotation_labels_from_xml(char* const filename);
    int set_state_name_from_xml(char* const filename);

 public:

    // Constructors

    model_parameters();

    model_parameters(char* const m_type,
		  const bool p,
		  const bool se,
		  const int n_of_special_emissions,	
		  const int total_n_of_a_labels,
		  Label* const a_label,	
		  const int total_n_of_m_elements,
		  char** const m_elements,
		  const Letter alph,
		  const int n_of_states,
		  char** const s_name);
 
    model_parameters(char* file_name, int& check);

    // Destructor
    ~model_parameters();

    // accessors 

    char* get_Model_Name() const;
    bool is_PairHMM() const;
    bool is_SpecialEmit() const;

    int get_Number_of_Special_Emissions() const;
    int get_Max_n_of_Child_State() const;

    int get_Total_Number_of_Annotation_Labels() const;
    int get_Annotation_Label_size(const int i) const;
    char* get_Annotation_Label_setname(const int i) const;
    char* get_Annotation_Label_name(const int i, const int j) const;
    char** get_Annotation_Label_name(const int i) const;
    Label* get_Annotation_Label() const;
    bool get_score_of_Annotation_Label(const int i) const;

    bool get_Alphabet_cases() const;
    int get_Alphabet_size()const;
    char get_Alphabet_name(const int i) const;
    char* get_Alphabet_name() const;
    Letter get_Alphabet() const;

    int get_Number_of_States() const;
    char* get_State_Name(const int i) const;
    char** get_State_Name() const;

    // mutators

    int set_Model_Name(char* const name);
    int set_pair(const bool p);
    int set_SpecialEmit(const bool se);

    int set_Number_of_Special_Emissions(const int n_of_special_emissions);
    int set_Max_n_of_Child_State(const int max_n_of_child_state);
    
    int set_Total_Number_of_Annotation_Labels(const int number);
    int set_Annotation_Label_size(const int i, const int s);
    int set_Annotation_Label_setname(const int i, char * const sn);
    int set_Annotation_Label_name(const int i, const int j, char* const n);
    int set_score_of_Annotation_Label(const int i, const bool s);

    int set_Alphabet_cases(const bool c);
    int set_Alphabet_size(const int s);
    int set_Alphabet_name(const int i, const char n);   

    int set_Number_of_States(const int number);
    int set_State_Name(const int i, char* const name);
    
    // Other Functions
    
    void print(std::ostream &o) const;

};

int convert_state_name_to_int(model_parameters* const MP, char* const buffer);
int convert_alphabet_to_int(model_parameters* const MP, const char buffer);
int convert_State_id_to_int(model_parameters* const MP, const char* buffer);
int convert_Annotation_Label_id_to_int(model_parameters* const MP, const char* buffer, int type);
int convert_Annotation_Label_set_id_to_int(model_parameters* const MP, const char* buffer );

#endif
