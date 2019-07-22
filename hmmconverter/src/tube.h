/* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define two-dimensional tube
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/blastn.h,v 1.3 2008/12/14 10:39:23 natural Exp $
*/

#ifndef _tube_h
#define _tube_h

#include <iostream>

template<class T>
class Tube
{
  public:

  Tube<T>       ( void );                                   
  Tube<T>       ( const int length_x, const int length_y); 
  Tube<T>       ( const Tube<T> &v );                       
 ~Tube<T>      ( void );                                   

  int             Lx ( void ) const;          
  int             Ly ( void ) const;          
  int             Empty ( void ) const;       
  int             Reset ( void );       

  void            SetElement       ( const int index_x, const int index_y, const T &v );
  T               GetElement       ( const int index_x, const int index_y) const;       

  Tube<T>   &     operator=        ( const Tube<T> &v);
  int             Equal            ( const Tube<T> &v) const;
  
 private:

  T**    data;
  int    L_x;
  int    L_y;
};

template<class T> Tube<T>::Tube ( void ) {

  data = NULL;
  L_x  = 0;
  L_y  = 0;
}

template<class T> Tube<T>::Tube ( const int length_x, const int length_y) {

  int j, k;

  data = NULL;
  L_x  = 0;
  L_y  = 0;

  if ((length_x > 0) && (length_y > 0)) {

    L_x = length_x;
    L_y = length_y;

    data = new T*[L_x];
    for (j = 0; j < L_x; j++) {
      data[j] = new T[L_y];
      for (k = 0; k < L_y; k++) {
	data[j][k] = static_cast<T>(0);
      }
    }
  }
  else {
    std::cout << "ERROR: Tube<T>:Tube(const int length_x, const int length_y): length_x (" << length_x
	 << ") or length_y (" << length_y << ") are 0.\n" << std::flush;
  }
}

template<class T> Tube<T>::Tube ( const Tube<T> &v ) {

  data = NULL;
  L_x  = 0;
  L_y  = 0;

  *this=v;
}

template<class T> Tube<T>::~Tube ( void ) {

  int i, j;

  if (data) {
    for (i=0; i < L_x; i++) {
      if (data[i]) {
	delete [] data[i]; 
	data[i] = NULL;
      }
    }
    if (data) delete [] data;
    data = NULL;
  }
  L_x = 0;
  L_y = 0;
}

template<class T> int Tube<T>::Lx (void) const {
  return(L_x);
}

template<class T> int Tube<T>::Ly (void) const {
  return(L_y);
}

template<class T> int Tube<T>::Empty (void) const {
  if (data) {
    return(0);
  }
  else {
    return(1);
  }
}

template<class T> int Tube<T>::Reset( void ) {

  int i, j;

  if (data) {
    for (i=0; i < L_x; i++) {
      if (data[i]) {
	delete [] data[i]; 
	data[i] = NULL;
      }
    }
    if (data) delete [] data;
    data = NULL;
  }
  L_x = 0;
  L_y = 0;

  return(0);
}


template<class T> T Tube<T>::GetElement( const int index_x, const int index_y ) const {

  if (data && (0 <= index_x) && (index_x < L_x) && (0 <= index_y) && (index_y < L_y)) {

    return(data[index_x][index_y]);
  }
  else {
    std::cout << "ERROR: Tube<T>:GetElement: index_x (" << index_x << ") out of range [0, " << L_x-1
      << "] or index_y (" << index_y << ") out of range [0, " << L_y-1 << "] or data (" << data
	 << ") is NULL.\n";
    return(0);
  }
}

template<class T> void Tube<T>::SetElement( const int index_x, const int index_y, const T &v) {

  if (data && (0 <= index_x) && (index_x < L_x) && (0 <= index_y) && (index_y < L_y)) {

    data[index_x][index_y] = v;
  }
  else {
    std::cout << "ERROR: Tube<T>:SetElement: index_x (" << index_x << ") out of range [0, " << L_x-1
      << "] or index_y (" << index_y << ") out of range [0, " << L_y-1 << "] or data (" << data
	 << ") is NULL.\n";
  }
  return;
}

template<class T> int Tube<T>::Equal( const Tube<T> &v) const {

  int i, j;
  int equal = 1;
  
  if (this == &v) {
    equal = 1;
  }
  else {
    if (data && v.data) {
      if ((L_x == v.L_x) && (L_y == v.L_y)) {
	
	for (i=0; i<L_x; i++) {
	  for (j=0; j<L_y; j++) {
	    if (data[i][j] != v.data[i][j]) {
	      equal = 0;
	      break;
	    }
	  }
	  if (equal == 0) {break;}
	}
      }
      else {
	equal = 0;
      }
    }
    else {
      equal = 0;
    }
  }
  return(equal);
}

template<class T> Tube<T> & Tube<T>::operator= ( const Tube<T> &v)
{
  int i, j, k;

  if (this != &v) {

    if (data) {
      for (i=0; i < L_x; i++) {
	if (data[i]) {
	  delete [] data[i]; 
	  data[i] = NULL;
	}
      }
      if (data) delete [] data;
      data = NULL;
    }

    L_x = v.L_x;
    L_y = v.L_y;
    
    if(v.data)
    {
	data = new T*[L_x];
	for (j = 0; j < L_x; j++) {
	    data[j] = new T[L_y];
	    for (k = 0; k < L_y; k++) {
		data[j][k] = v.data[j][k];
	    }
	}
    }
  }
  return(*this);
}

#endif
