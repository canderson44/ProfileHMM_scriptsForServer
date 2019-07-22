 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/transitionprob.h,v 1.3 2008/12/14 10:39:23 natural Exp $
 */

#ifndef _transitionprob_h
#define _transitionprob_h

#include <iostream>
#include "define.h"
#include "model_parameters.h"
#include "tinyxml.h"
#include "hmm.h"
#include "models.h"

struct FTProb{
    char* name;
    Prob prob;
    bool train;
    char* exp;
    Prob pseudocount;
};

struct GTProb{
    int NumOffromto;
    int* from;
    int* to;
    int NumOfOverfromto;
    int* Overfrom;
    int* Overto;
    bool trained;
    double score;
};

struct TTProb{
    int from;
    int to;
    double score;
    bool trained;
    double pseudoprob;
};

class TransitionProb{
 private:
    int FTPsize;
    int GTPsize;
    int TTPsize;
    char* FTP_tname;
    char* GTP_tname;
    FTProb* FTP;
    GTProb* GTP;
    TTProb* TTP;
    int set_FTP_from_xml(const char* In_dir,
			 const char* filename);   

    int get_FTP_parameters_for_training(const char*  filename,
					model_parameters* const MP);
    int get_GTP_parameters_for_training(const char*  filename,
					model_parameters* const MP,
					Hmm* const p);
    int get_TTP_parameters_for_training(const char*  filename,
					model_parameters* const MP,
					Hmm* const p);   

    int check_consistency_of_GTP(Hmm* const p);
    

 public:
    //constructor
    TransitionProb();

    TransitionProb(const char* In_dir,
		   const char* filename, 
		   int& check);
    ~TransitionProb();

    //accessors
    // for FTP
    int get_FTPsize() const;
    char* get_FTP_tname() const;
    char* get_FTP_name(const int i) const ;
    Prob get_FTP_prob(const int i) const ;
    FTProb get_FTP(const int i) const;
    bool is_FTP_train(const int i) const;
    char* get_FTP_exp(const int i) const;
    Prob get_FTP_pseudocount(const int i) const;

    // for GTP
    int get_GTPsize() const;
    char* get_GTP_tname() const;
    int get_GTP_NumOffromto(const int i) const;
    int get_GTP_NumOfOverfromto(const int i) const;
    int get_GTP_from(const int i, const int j) const;
    int get_GTP_to(const int i, const int j) const;
    int get_GTP_Overfrom(const int i, const int j) const;
    int get_GTP_Overto(const int i, const int j) const;
    bool get_GTP_trained(const int i) const;
    double get_GTP_score(const int i) const;
    
    // for TTP
    int get_TTPsize() const;
    int get_TTP_from(const int i) const;
    int get_TTP_to(const int i) const;
    double get_TTP_score(const int i) const;
    bool get_TTP_trained(const int i) const;
    double get_TTP_pseudoprob(const int i) const;
    
    //mutators
    // for FTP
    int set_FTPsize(const int s);
    int set_FTP_tname(char* const tn);
    int set_FTP_name(const int i, char* const n);
    int set_FTP_prob(const int i, const Prob p);
    int set_FTP_train(const int i, const bool t);
    int set_FTP_exp(const int i, char* const e);
    int set_FTP_pseudocount(const int i, const Prob p);

    // for GTP
    int set_GTPsize(const int s);
    int set_GTP_tname(char* const tn);
    int set_GTP_NumOffromto(const int i, const int num);
    int set_GTP_NumOfOverfromto(const int i, const int num);
    int set_GTP_from(const int i, const int j, const int f);
    int set_GTP_to(const int i, const int j, const int t);
    int set_GTP_Overfrom(const int i, const int j, const int Of);
    int set_GTP_Overto(const int i, const int j, const int Ot);    
    int set_GTP_trained(const int i, const bool t);
    int set_GTP_score(const int i, const double s);

    // for TTP
    int set_TTPsize(const int s);
    int set_TTP_from(const int i, const int f);
    int set_TTP_to(const int i, const int t);
    int set_TTP_score(const int i, const double s);
    int set_TTP_trained(const int i,const bool t);
    int set_TTP_pseudoprob(const int i, const double p);
    //operators
    TransitionProb & operator = (const TransitionProb &tp);

    // other functions
    
    int get_parameters_for_training(const char* filename,
				    model_parameters* const MP,
				    Hmm* const p);

    int derive_FTPs_from_GTPs();
	
    double evaluate_FTP_expression(const char* post, int check);

    void print(std:: ostream &o) const;
};


int convert_FTP_id_to_int(TransitionProb* const TP, const char* buffer);

int convert_GTP_id_to_int(TransitionProb* const TP, const char* buffer);

#endif
