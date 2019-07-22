/* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: declare match-class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/match.cpp,v 1.1.1.1 2008/07/13 09:39:09 natural Exp $

 */

#include <iostream>
#include "match.h"

/* define constructors */

Match::Match() {

  x_start = 0;
  y_start = 0;
  x_end   = 0;
  y_end   = 0;
  score   = 0.0;
}

Match::Match(const int start_x, const int start_y, const int end_x, const int end_y, const double match_score) {

  x_start = start_x;
  y_start = start_y;
  x_end   = end_x;
  y_end   = end_y;
  score   = match_score;
}

Match & Match::operator = (const Match &m) {

  if (this != &m) {

    x_start = m.x_start;
    y_start = m.y_start;
    x_end   = m.x_end;  
    y_end   = m.y_end;  
    score   = m.score;  
  }
  return(*this);
}

int Match::operator == (const Match &m) const {

    // two matches are equal if the same of identical

    if (this != &m) {

	if ((x_start == m.x_start) &&
	    (y_start == m.y_start) &&
	    (x_end   == m.x_end)   &&
	    (y_end   == m.y_end)   &&
	    (score   == m.score)) {
	    
	    return(1);
	}
	else {
	    return(0);
	}
    }
    return(1);
}

int Match::operator < (const Match &m) const {

    if (this != &m) {
	
	if ((max(x_start, x_end) < min(m.x_start, m.x_end)) && 
	    (max(y_start, y_end) < min(m.y_start, m.y_end)))
	{
	    return(1);
	}
	else 
	{
	    return(0);
	}
    }
    return(0);
}

int Match::operator > (const Match &m) const {

  if (this != &m) {

    if ((max(m.x_start, m.x_end) < min(x_start, x_end)) && 
	(max(m.y_start, m.y_end) < min(y_start, y_end)))
      {return(1);}
    else {return(0);}
  }
  return(0);
}

Match::Match(const Match &m) {

  *this = m;
}

Match::~Match() {}

double Match::Score(void)   const {return(score);}
int    Match::Start_x(void) const {return(x_start);}
int    Match::Start_y(void) const {return(y_start);}
int    Match::End_x(void)   const {return(x_end);}
int    Match::End_y(void)   const {return(y_end);}

void Match::print(std::ostream &o) const {

  o << "(" << x_start << ", " << y_start << ", " << x_end << ", " << y_end << ", " << score << ")\n";

  return;
}

