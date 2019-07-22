 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare match-class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/match.h,v 1.1.1.1 2008/07/13 09:39:09 natural Exp $
 */

#ifndef _match_h
#define _match_h

#include <stdio.h>  

#define min(a,b)  (((a) <= (b)) ? (a) : (b))
#define max(a,b)  (((a) >= (b)) ? (a) : (b))


class Match
{
 private:

    int    x_start;
    int    y_start;
    int    x_end;
    int    y_end;
    double score;
    int    number;
    
 public:
    
    // constructors and destructor
    
    Match();                               
    Match(const int start_x, const int start_y, const int end_x, const int end_y, const double match_score);
    Match(const Match &m);
    ~Match();

    // functions

    double Score(void)   const;
    int    Start_x(void) const;
    int    Start_y(void) const;
    int    End_x(void)   const;
    int    End_y(void)   const;
    void   print(std::ostream &o) const;

    // operators
    
    Match & operator = (const Match &m);
    int operator == (const Match &m) const;
    int operator <  (const Match &m) const;
    int operator >  (const Match &m) const;
};

#endif
