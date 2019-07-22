/*
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define the emissionprob class

   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/emissionprob.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _emissionprob_h
#define _emissionprob_h

#include <iostream>
#include "define.h"
#include "model_parameters.h"
#include "tinyxml.h"
#include "hmm.h"

struct FEProb{    
    int dim;    
    char* name;
    array<Prob> prob;
    array<Prob> pseudoprob;
    bool train;
    char* exp;
};

struct GEProb{
    int NumOffrom;
    int* from;
};

class EmissionProb{

 private:
    int FEPsize;
    char* FEP_tname;
    FEProb* FEP;
    int GEPsize;
    char* GEP_tname;
    GEProb* GEP;

    int set_FEP_from_file(const char* In_dir,
			  const char* xmlfile,
			  model_parameters* const MP);

        
    int get_FEP_parameters_for_training(const char* filename,
					int& NumOfAutoEP,
					int** AutoEP);

    int get_GEP_parameters_for_training(const char* filename,
					model_parameters* const MP,
					int NumOfAutoEP);

    int auto_get_parameters_for_training(Hmm* const p,
					 const int NumOfAutoEP,
					 int* AutoEP);

    int check_consistency_of_GEP(Hmm* const p);

    
 public:
    //constructor
    EmissionProb();
    EmissionProb(const char* In_dir,
		 const char* xmlfile, 
		 model_parameters* const MP, 
		 int& check);
    ~EmissionProb();

    //accessors

    //For FEP
    int get_FEPsize() const;
    char* get_FEP_tname() const;
    char* get_FEP_name(const int i) const ;
    int get_FEP_dim(const int i) const;
    Prob get_FEP_prob(const int i, const array<int> &indices) const;
    Prob get_FEP_prob(const int i, const int linear_index) const;
    array<Prob> get_FEP_prob(const int i)const;  
    
    Prob get_FEP_pseudoprob(const int i, const array<int> &indices) const;
    Prob get_FEP_pseudoprob(const int i, const int linear_index) const;
    array<Prob> get_FEP_pseudoprob(const int i) const;
    
    FEProb get_FEP(const int i) const;
    bool is_FEP_train(const int i) const;
    char* get_FEP_exp(const int i) const;
    
    //For GEP
    int get_GEPsize() const;
    char* get_GEP_tname() const;
    int get_GEP_NumOffrom(const int i) const;
    int get_GEP_from(const int i, const int j) const;

    //mutators

    // For FEP
    int set_FEPsize(const int s);
    int set_FEP_tname(char* const tn);
    int set_FEP_name(const int i, char* const n);
    int set_FEP_dim(const int i, const int d);
    int set_FEP_prob(const int i, const array<int> &indices, const Prob p);
    int set_FEP_prob(const int i, const int linear_index, const Prob p);
    int set_FEP_pseudoprob(const int i, const array<int> &indices, const Prob p);
    int set_FEP_pseudoprob(const int i, const int linear_index, const Prob p);
    
    int set_FEP_train(const int i,const bool t);
    int set_FEP_exp(const int i, char* const e);

    // For GEP
    int set_GEPsize(const int s);
    int set_GEP_tname(char* const tn);
    int set_GEP_NumOffrom(const int i, const int num);
    int set_GEP_from(const int i, const int j, const int f);

    //operators
    EmissionProb & operator = (const EmissionProb &ep);
    
    // other functions

    int get_parameters_for_training(const char* filename,
				    model_parameters* const MP,
				    Hmm* const p);

    

    void print(std:: ostream &o) const;
 
};

int convert_FEP_id_to_int(EmissionProb* const EP, const char* buffer);

int convert_GEP_id_to_int(EmissionProb* const EP, const char* buffer);

#endif
