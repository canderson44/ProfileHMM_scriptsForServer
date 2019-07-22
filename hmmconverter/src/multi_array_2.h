/*
   Authors: Irmtraud M Meyer
   Copyright: Irmtraud M Meyer (1999-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: define multi-dimensional arrays using recursive templates
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/multi_array_2.h,v 1.1.1.1 2008/07/13 09:39:09 natural Exp $
*/

#ifndef _multi_array_2_h
#define _multi_array_2_h

#include <iostream>
#include <math.h>
#include "parameters.h"

template<class T>
class array
{
 public:

  array<T>       ( void );
  array<T>       ( const int d );        
  array<T>       ( const array<T> &v );
  ~array<T>      ( void );

  void             SetDimension     ( int i, int n );
  void             SetDimension     ( int i, int n, T initial_value );
  void             SetNumberofDimensions (int i);

  int              GetDimension     ( int i ) const;
  int              GetNumberofDimensions (void) const;
  int              Empty(void) const;

  void             SetElement       ( const array<int> &indices, const T &v );
  void             SetElement       ( const int linear_index, const T &v );

  T                GetElement       ( const array<int> &indices) const;
  T                GetElement       ( const int linear_index ) const;

  array<T>   &     operator=        ( const array<T> &v);
  int              Equal            ( const array<T> &v) const;
  array<T>         operator*        ( const array<T> &b) const;
  array<T>         SpecialMultiply  ( const array<T> &b) const; 				      
  
  void             Print            ( std::ostream &o ) const;
  void             PrintwithoutNewLine ( std::ostream &o ) const;
  void             PrintwithIndices ( std::ostream &o ) const;
  void             PrintwithIndices ( std::ostream &o, const char* int_to_char_mapping) const;
  void             PrintwithIndicesRestrict ( std::ostream &o, const T dontprint) const;
  void             PrintonlyNonZerowithIndices ( std::ostream &o ) const;
  void             PrintonlyNonLogzerowithIndices ( std::ostream &o ) const;
  void             PrintDimensions  ( std::ostream &o ) const;
  void             ResetData        (void);
  void             ResetData        (T initial_value);
  void             ResetDataandDimensions (void);
  void             Reset            (void);

 private:

  T*     data;
  int*   dim;
  int    number_of_dimensions; 
};

void convert_to_base(int number, const int base, array<int>* index);

template<class T> array<T>::array ( void )
{
  data = NULL;
  dim = NULL;
  number_of_dimensions=0;
}

template<class T> array<T>::array ( const int d )
{
#ifdef _DEBUG
  if (d<1)
    {std::cout << "<ERROR: array::Constructor: number of dimensions cannot be " << d << " < 1\n";}
#endif
  data = NULL;
  dim = NULL;
  number_of_dimensions=d;

  if (d>0)
    {
      dim = new int[d];  
      for (int j=0; j<number_of_dimensions; j++)
	{dim[j]=0;}
    }
}

template<class T> array<T>::array ( const array<T> &v )
{
  dim=NULL;
  data=NULL;
  number_of_dimensions=0;

  *this=v;
}

template<class T> array<T>::~array ( void )
{
  if (data) delete [] data;
  data=NULL;
  if (dim) delete [] dim;
  dim=NULL;
  number_of_dimensions=0;
}

template<class T> void array<T>::SetNumberofDimensions (int i)
{
#ifdef _DEBUG
  if (i<0)
    {std::cout << "<ERROR:array::SetNumberofDimensions: number of dimensions cannot be " << i << " < 0\n";}
  if (number_of_dimensions>0)
    {std::cout << "<ERROR:array::SetNumberofDimensions: number of dimensions has already been set to " << number_of_dimensions << "\n";}
#endif
  if (dim) delete [] dim;
  dim=NULL;
  if (data) delete [] data;
  data=NULL;

  number_of_dimensions=i;
  if (number_of_dimensions>0) 
    {
      dim= new int[number_of_dimensions];
      for (int j=0; j<number_of_dimensions; j++)
	{dim[j]=0;}
    }
  return;
}

template<class T> void array<T>::SetDimension ( int i, int n )
{
#ifdef _DEBUG
  if (i<0)
    {std::cout << "<ERROR:array::SetDimension: number of dimension cannot be " << i << " < 0\n";}
  if (n<1 && number_of_dimensions>0)
    {std::cout << "<ERROR:array::SetDimension: length of dimension i = " << i << " cannot be " << n << " < 1\n";}
  if ( (!dim) || (i>number_of_dimensions-1) )
    {
      if (!dim)
	{std::cout << "<ERROR:array::SetDimension: array dim has not been filled \n";}
      if (i>number_of_dimensions-1)      
	{std::cout << "<ERROR:array::SetDimension: dimensions i = " << i << " out of range. Has to be "
	      << " < " << number_of_dimensions-1 << " (number_of_dimensions-1)\n";}
    }
  else
    {
      int check=0;
      for (int j=0; j<i; j++)
	{if (dim[j]<1) check+=1;}
      if (check!=0)
	{std::cout << "<ERROR:array::SetDimension: length of a previous dimensions i < " << i << " has not been set\n";}
    }
#endif
  dim[i]=n;

  if (data) delete [] data;
  data = NULL;

  int length=dim[0];
  for (int j=1;j<number_of_dimensions; j++)
    {length*=dim[j];}

  if (length>0)
    {
      data= new T [length];
      for (int j=0;j<length;j++)
	{data[j]=static_cast<T>(0);}
    }
  return;
}

template<class T> void array<T>::SetDimension ( int i, int n, T initial_value )
{
#ifdef _DEBUG
  if (i<0)
    {std::cout << "<ERROR:array::SetDimension: number of dimension cannot be " << i << " < 0\n";}
  if (n<1 && number_of_dimensions>0)
    {std::cout << "<ERROR:array::SetDimension: length of dimension i = " << i << " cannot be " << n << " < 1\n";}
  if ( (!dim) || (i>number_of_dimensions-1) )
    {
      if (!dim)
	{std::cout << "<ERROR:array::SetDimension: array dim has not been filled \n";}
      if (i>number_of_dimensions-1)      
	{std::cout << "<ERROR:array::SetDimension: dimensions i = " << i << " out of range. Has to be "
	      << " < " << number_of_dimensions-1 << " (number_of_dimensions-1)\n";}
    }
  else
    {
      int check=0;
      for (int j=0; j<i; j++)
	{if (dim[j]<1) check+=1;}
      if (check!=0)
	{std::cout << "<ERROR:array::SetDimension: length of a previous dimensions i < " << i << " has not been set\n";}
    }
#endif
  dim[i]=n;

  if (data) delete [] data;
  data = NULL;

  int length=dim[0];
  for (int j=1;j<number_of_dimensions; j++)
    {length*=dim[j];}

  if (length>0)
    {
      data= new T [length];
      for (int j=0;j<length;j++)
	{data[j]=initial_value;}
    }
  return;
}

template<class T> int array<T>::GetDimension ( int i ) const
{
#ifdef _DEBUG
  if (i<0)
    {std::cout << "<ERROR:array::GetDimension: number of dimension cannot be " << i << " < 0\n";}
  if (i>number_of_dimensions-1)      
    {std::cout << "<ERROR:array::GetDimension: dimensions i = " << i << " out of range. Has to be "
	  << " < " << number_of_dimensions-1 << " (number_of_dimensions-1)\n";}
#endif
  return(dim[i]);
}

template<class T> int array<T>::Empty ( void ) const
{
  if (data == NULL) {
    return(1);
  }
  else {
    return(0);
  }
}

template<class T> T array<T>::GetElement( const int i ) const
{
#ifdef _DEBUG
  if (i<0)
    {std::cout << "<ERROR:array::GetElement: there is no element with linear index i = " << i << " < 0\n";}

  int max_index=dim[0];
  for (int j=1;j<number_of_dimensions; j++)
    {max_index*=dim[j];}
  max_index-=1;

  if (i>max_index)
    {std::cout << "<ERROR:array::GetElement: there is no element with linear index i = " << i 
	    << " > max_index = " << max_index << "\n";}
#endif
  return(data[i]);
}

template<class T> T array<T>::GetElement ( const array<int> &indices) const
{
#ifdef _DEBUG
  if ( (indices.GetNumberofDimensions()!=1) || (indices.GetDimension(0)!=number_of_dimensions) )
    {
      if (indices.GetNumberofDimensions()!=1)
	{std::cout << "<ERROR:array::GetElement: array of indices is no vector, but has number of dimensions "
	      << indices.GetNumberofDimensions() << "\n";}
      if (indices.GetDimension(0)!=number_of_dimensions)
	{std::cout << "<ERROR:array::GetElement: number of indices is not " 
	      << number_of_dimensions << " (number_of_dimensions), but "
	      << indices.GetDimension(0) << "\n";}
    }
  else
    {
      int check=0;
      for (int j=0; j<number_of_dimensions; j++)
	{
	  if ( (indices.GetElement(j)<0) || (indices.GetElement(j)>(dim[j]-1)) )
	    {check+=1;}
	}
      if (check!=0)
	{std::cout << "<ERROR:array::GetElement: at least one of the indices of the index array is out of range\n";}
    }
  int max_index=dim[0];
  for (int j=1;j<number_of_dimensions; j++)
    {max_index*=dim[j];}
  max_index-=1;
#endif

  int linear_index=0;
  if (indices.GetDimension(0)==1) 
    {linear_index=indices.GetElement(0);}
  else
    {
      linear_index=dim[1]*indices.GetElement(0)+indices.GetElement(1);
      for (int j=2;j<number_of_dimensions;j++)
	{linear_index=dim[j]*linear_index+indices.GetElement(j);}
    }

#ifdef _DEBUG
  if (linear_index>max_index)
    {std::cout << "ERROR:array:GetElement: linear_index = " << linear_index << " > max_index = " << max_index << "\n";}
#endif

  return(data[linear_index]);
}

template<class T> void array<T>::SetElement ( const array<int> &indices, const T &v )
{
#ifdef _DEBUG
  if ( (indices.GetNumberofDimensions()!=1) || (indices.GetDimension(0)!=number_of_dimensions) )
  {
      if (indices.GetNumberofDimensions()!=1)
      {std::cout << "<ERROR:array::SetElement: array of indices is no vector, but has number of dimensions "
		 << indices.GetNumberofDimensions() << "\n";}
      if (indices.GetDimension(0)!=number_of_dimensions)
      {std::cout << "<ERROR:array::SetElement: number of indices is not " 
		 << number_of_dimensions << " (number_of_dimensions), but "
		 << indices.GetDimension(0) << "\n";}
  }
  else
  {
      int check=0;
      for (int j=0; j<number_of_dimensions; j++)
      {
	  if ( (indices.GetElement(j)<0) || (indices.GetElement(j)>(dim[j]-1)) )
	  {
	      check+=1;
	  }
      }
      if (check!=0)
      {
	  std::cout << "<ERROR:array::SetElement: at least one of the indices of the index array is out of range\n";
      }
  }
  int max_index=dim[0];
  for (int j=1;j<number_of_dimensions; j++)
  {
      max_index*=dim[j];
  }
  max_index-=1;
#endif

  int linear_index=0;
  if (indices.GetDimension(0)==1) 
  {
      linear_index=indices.GetElement(0);
  }
  else
  {
      linear_index=dim[1]*indices.GetElement(0)+indices.GetElement(1);
      for (int j=2;j<number_of_dimensions;j++)
      {
	  linear_index=dim[j]*linear_index+indices.GetElement(j);}
  }

#ifdef _DEBUG
  if (linear_index>max_index)
  {
      std::cout << "ERROR:array:SetElement: linear_index = " << linear_index << " > max_index = " << max_index << "\n";
  }
#endif

  data[linear_index]=v;
  return;
}

template<class T> void array<T>::SetElement( const int i, const T &v )
{
#ifdef _DEBUG
  if (i<0)
    {std::cout << "<ERROR:array::SetElement: there is no element with linear index i = " << i << " < 0\n";}

  int max_index=dim[0];
  for (int j=1;j<number_of_dimensions; j++)
    {max_index*=dim[j];}
  max_index-=1;

  if (i>max_index)
    {std::cout << "<ERROR:array::SetElement: there is no element with linear index i = " << i 
	    << " > max_index = " << max_index << "\n";}
#endif
  data[i]=v;
  return;
}

template<class T> int array<T>::GetNumberofDimensions (void) const
{
  return(number_of_dimensions);
}

template<class T> array<T> array<T>::operator* ( const array<T> &b) const
{
  // a == this 
  // a = (n_a x m_a) matrix
  // b = (n_b x m_b) matrix
  // need m_a == n_b in order to multiply them
  //
  // note: if multiplication cannot be done, return array with 0 dimensions

  if ((number_of_dimensions==2) && (b.number_of_dimensions==2) &&
      (dim[1]==b.dim[0])) 
    {
      // create new matrix for product of multiplication

      array<T> product(2);
      product.SetDimension(0, dim[0]);
      product.SetDimension(1, b.dim[1]);
      
      // fill elements of new matrix: product(i,j) = sum_k (k=1..m_a) [a(i,k)*b(k,j)]

      array<int> index(1);
      index.SetDimension(0, 2);

      array<int> index_a(1);
      index_a.SetDimension(0, 2);

      array<int> index_b(1);
      index_b.SetDimension(0, 2);

      T element=static_cast<T>(0);

      for (int i=0; i<dim[0]; i++)
	{
	  for (int j=0; j<b.dim[1]; j++)
	    {
	      element=static_cast<T>(0);
	      
	      for (int k=0; k<dim[1]; k++)
		{
		  index_a.SetElement(0, i);
		  index_a.SetElement(1, k);

		  index_b.SetElement(0, k);
		  index_b.SetElement(1, j);

		  element+=this->GetElement(index_a) * b.GetElement(index_b);
		}

	      index.SetElement(0, i);
	      index.SetElement(1, j);

	      product.SetElement(index, element);
	    }
	}
      return(product);
    }
  else // if multiplication cannot be done, return array with 0 dimensions
    {
      array<T> product;
      return(product);
    }
}

template<class T> array<T> array<T>::SpecialMultiply( const array<T> &b) const
{
  // input :
  //
  // a == this 
  // a =  matrix with 2d dimensions, length in each dimension is l,
  // sum of all elements of a == 1
  // b =  matrix with 2d dimensions, length in each dimension is l
  // sum of all elements of b == 1
  //
  // notation i := (i1, i2, ..., i_d)
  //          j := (j1, j2, ..., j_d)

  // function fills elements of new matrix (number of dimensions = 2d, length = l):
  // new_element(i1,i2,..id,j1,j2,..jd) = sum_c1=1^l sum_c2=1^l
  // ... sum_cd=1^l this(i1,i2,...id, c1,c2,...cd) * 
  // b(c1, c2,...cd, j1,j2,...jd)
  //
  // note: if multiplication cannot be done, function returns array with 0 dimensions

  int check=0;

  int length       = dim[0];
  int dimensions   = number_of_dimensions;
  int dimensions_2 = static_cast<int>(static_cast<float>(number_of_dimensions)/2.);

  array<T> product(dimensions);
  for (int i=0; i<dimensions; i++)
    {product.SetDimension(i, length);}

  if ((number_of_dimensions > 0)    &&
      (number_of_dimensions%2 == 0) && 
      (number_of_dimensions == b.number_of_dimensions))
    {
      // more checks

      for (int i=0; i<number_of_dimensions; i++)
	{
	  if (dim[i] != b.dim[i]) 
	    {
	      check++;
	      break;}
	  else 
	    {
	      if ((dim[i] != length) ||
		  (b.dim[i] != length))
		{
		  check++;
		  break;}
	    }
	}
      if (check != 0)
	{
	  std::cout << "ERROR: class array<T>:: SpecialMultiply: arrays must "
	       << "have same length in each dimension and "
	       << "these lengths must be the same for both arrays.\n" << std::flush;
	}
    }
  else
    {
      std::cout << "ERROR: class array<T>:: SpecialMultiply: arrays must "
	   << "have same number of dimensions and the number of "
	   << "dimensions must be an even number.\n" << std::flush;
    }
  
  // check that both matrices are positive definite and their
  // elements add up to one

  if (check == 0)
    {
      int max = static_cast<int>(pow( static_cast<float>(length), static_cast<float>(dimensions)));
      T element             = static_cast<T>(0);
      T sum_of_all_elements = static_cast<T>(0);
      
      // check this matrix

      {      
	for (int l_index=0; l_index<max; l_index++) // loop over all possible indices
	  {
	    element = static_cast<T>(0);
	    element = this->GetElement(l_index);
	    
	    if (element < 0) 
	      {
		std::cout << "ERROR: array<T>::SpecialMultiply: elements of this matrix not >= 0\n" << std::flush;
		check++;
		break;
	      }
	    else
	      {
		sum_of_all_elements += element;
	      }
	  }
      }
      
      if (fabs(sum_of_all_elements - static_cast<T>(1.)) > Max_deviation)
	{
	  std::cout << "ERROR: array<T>::SpecialMultiply: this matrix is not normalised "
	       << "i.e. the sum of all its elements is not 1\n" << std::flush;
	  check++;
	}

      // check matrix b

      sum_of_all_elements = static_cast<T>(0);
      
      {
	for (int l_index=0; l_index<max; l_index++) // loop over all possible indices
	  {
	    element = static_cast<T>(0);
	    element = this->GetElement(l_index);
	    
	    if (element < 0) 
	      {
		std::cout << "ERROR: array<T>::SpecialMultiply: elements of matrix b not >= 0\n" << std::flush;
		check++;
		break;
	      }
	    else
	      {
		sum_of_all_elements += element;
	      }
	  }
      }

      if (fabs(sum_of_all_elements - static_cast<T>(1.)) > Max_deviation)
	{
	  std::cout << "ERROR: array<T>::SpecialMultiply: matrix b is not normalised "
	       << "i.e. the sum of all its elements is not 1\n" << std::flush;
	  check++;
	}
    }
  
  if (check == 0)
    {
#ifdef _PRINT
      std::cout << "array<T>::SpecialMultiply:\n"
	   << "----------------------------------------------------------------------\n" << std::flush;
#endif
      // variables
      
      // s_max = number of different d/2-tuples
      
      int s_max = static_cast<int>(pow( static_cast<float>(length), static_cast<float>(dimensions_2)));	  
      
      array<int> index(1);    // for indices of a d-tuple
      index.SetDimension(0, dimensions);
      
      array<int> index_1(1);  // for indices of a d-tuple
      index_1.SetDimension(0, dimensions);
      
      array<int> index_2(1);  // for indices of a d-tuple
      index_2.SetDimension(0, dimensions);

      array<int> s_index(1);  // for indices of a d/2-tuple
      s_index.SetDimension(0, dimensions_2);

      array<int> s_1_index(1);  // for indices of a d/2-tuple
      s_1_index.SetDimension(0, dimensions_2);

      array<int> s_2_index(1);  // for indices of a d/2-tuple
      s_2_index.SetDimension(0, dimensions_2);

      array<int> s_3_index(1);  // for indices of a d/2-tuple
      s_3_index.SetDimension(0, dimensions_2);

      T sum_of_row                = static_cast<T>(0);
      T sum_of_all_elements       = static_cast<T>(0);
      T element                   = static_cast<T>(0);
      T product_element           = static_cast<T>(0);

      // copy this matrix and matrix b into new matrices this_t
      // and b_t

      array<T> this_t = *this;
      array<T> b_t    = b;

      // convert this_t and b_t into transition matrices using:
      //  
      //       this_t(i, j) = this(i, j) / sum_j' this(i, j')
      //       b_t(i,j)     = b(i,j) / sum_j' b(i, j')
      
      // 1.) calculate this_t:
      
      // loop over all d/2-tuples i

      array<T>* matrix_pointer[2];
      matrix_pointer[0] = &this_t;
      matrix_pointer[1] = &b_t;	  
	  
      for (int matrix=0; matrix<2; matrix++)
	{
#ifdef _PRINT
	  if (matrix == 0)
	    {
	      std::cout << "transforming this emission prob matrix into "
		   << "transition prob matrix this_t\n" << std::flush;
	    }
	  else if (matrix == 1)
	    {
	      std::cout << "transforming emission prob matrix b into "
		   << "transition prob matrix b_t\n" << std::flush;
	    }
#endif
	  for (int s_1_l_index=0; s_1_l_index<s_max; s_1_l_index++)
	    {
	      convert_to_base(s_1_l_index, length, &s_1_index);

	      // a) calculate element = sum_j' this(i,j')

#ifdef _PRINT
	      std::cout << "   calculate element = sum_j' matrix(";
	      s_1_index.PrintwithoutNewLine(std::cout);		    
	      std::cout << ", j')\n" << std::flush;
#endif
	      element=static_cast<T>(0);

	      {	  
		for (int s_2_l_index=0; s_2_l_index<s_max; s_2_l_index++) 
		  //
		  // loop over all d/2-tuples j'
		  //
		  {
		    convert_to_base(s_2_l_index, length, &s_2_index);
		      
		    for (int i=0; i<dimensions_2; i++) 
		      {
			index.SetElement(i,              s_1_index.GetElement(i));
			index.SetElement(i+dimensions_2, s_2_index.GetElement(i));
		      }
		    // element += this(i, j')

		    element+=matrix_pointer[matrix]->GetElement(index);

#ifdef _PRINT
		    std::cout << "      (i, j') = (";
		    s_1_index.PrintwithoutNewLine(std::cout);		    
		    std::cout << ", ";
		    s_2_index.PrintwithoutNewLine(std::cout);
		    std::cout << ") = " 
			 << matrix_pointer[matrix]->GetElement(index)
			 << " += element = "
			 << element << "\n" << std::flush;
#endif
		  }
	      }
#ifdef _PRINT
	      std::cout << "   i              = ";
	      s_1_index.Print(std::cout);
	      std::cout << "   sum_j' (i, j') = " << element << "\n" << std::flush;
#endif

	      sum_of_row = static_cast<T>(0);

	      // b.)
	      // loop over all j and set 
	      // this_t(i,j) = this(i, j) / sum_j' this(i, j')
	      //             = this(i, j) / element
	      
	      {
		for (int s_2_l_index=0; s_2_l_index<s_max; s_2_l_index++) 
		  {
		    convert_to_base(s_2_l_index, length, &s_2_index);
		  
		    for (int i=0; i<dimensions_2; i++) 
		      {
			index.SetElement(i,              s_1_index.GetElement(i));
			index.SetElement(i+dimensions_2, s_2_index.GetElement(i));
		      }
#ifdef _PRINT
		    if (matrix == 0)
		      {
			std::cout << "      this_t(i, j) = (";
		      }
		    else if (matrix == 1)
		      {
			std::cout << "      b_t(i, j) = (";
		      }
		  
		    index.PrintwithoutNewLine(std::cout);		    
		    std::cout << ") = " 
			 << matrix_pointer[matrix]->GetElement(index)
			 << " / element (" 
			 << element << ") = ";
#endif
		    if (element != 0)
		      {
			matrix_pointer[matrix]->SetElement(index, matrix_pointer[matrix]->GetElement(index)/element);
			sum_of_row += matrix_pointer[matrix]->GetElement(index);
		      }
		    else if ((element == 0) &&
			     (matrix_pointer[matrix]->GetElement(index) == 0))
		      {sum_of_row += matrix_pointer[matrix]->GetElement(index);} 
		    else
		      {
			if (matrix == 0)
			  {
			    std::cout << "ERROR: array<T>::SpecialMultiply: sum_j' this(i, j') for i = ";
			    s_1_index.PrintwithoutNewLine(std::cout);
			    std::cout << " = 0 but this(i, j) (";
			    index.PrintwithoutNewLine(std::cout);
			    std::cout << ") = " 
				 << matrix_pointer[matrix]->GetElement(index) 
				 << " != 0\n" << std::flush;
			  }
			else if (matrix == 1)
			  {
			    std::cout << "ERROR: array<T>::SpecialMultiply: sum_j' b(i, j') for i = ";
			    s_1_index.PrintwithoutNewLine(std::cout);
			    std::cout << " = 0 but b(i, j) (";
			    index.PrintwithoutNewLine(std::cout);
			    std::cout << ") = " 
				 << matrix_pointer[matrix]->GetElement(index) 
				 << " != 0\n" << std::flush;
			  }
			check++;
		      }
		  
#ifdef _PRINT
		    std::cout << matrix_pointer[matrix]->GetElement(index) << "\n" << std::flush;
#endif
		  }
	      }
	    }
	}

      // now that we have calculated this_t and b_t, calculate
      // final matrix which will be saved in matrix product
      
      if (check == 0)
	{	
#ifdef _PRINT      
	  std::cout << "calculate final matrix:\n" << std::flush;
	  std::cout << "calculate final product(i, j) = sum_c [ this_t(i, c)* b_t(c, j) ] * element \n" << std::flush;
	  std::cout << "          final product(i, j) = sum_c [ product_element ]         * element \n" << std::flush;
#endif
	  sum_of_all_elements = static_cast<T>(0);

	  // loop over all first d indices i

	  for (int s_1_l_index=0; s_1_l_index<s_max; s_1_l_index++)
	    {
	      convert_to_base(s_1_l_index, length, &s_1_index);

	      // a) calculate element = this(i) = sum_j' this(i,j')

	      element=static_cast<T>(0);

#ifdef _PRINT
	      std::cout << "   calculate this(";
	      s_1_index.PrintwithoutNewLine(std::cout);	      
	      std::cout << ") = sum_j' this(";
	      s_1_index.PrintwithoutNewLine(std::cout);	      
	      std::cout << ", j')\n" << std::flush;
#endif
	      {
		for (int s_2_l_index=0; s_2_l_index<s_max; s_2_l_index++) 
		  {
		    convert_to_base(s_2_l_index, length, &s_2_index);
		      
		    for (int i=0; i<dimensions_2; i++) 
		      {
			index.SetElement(i,              s_1_index.GetElement(i));
			index.SetElement(i+dimensions_2, s_2_index.GetElement(i));
		      }
		    element+=this->GetElement(index);
#ifdef _PRINT
		    std::cout << "      element += this(";
		    index.PrintwithoutNewLine(std::cout);		    
		    std::cout << ") (" 
			 << this->GetElement(index)
			 << ") = " 
			 << element << "\n" << std::flush;
#endif
		  }
	      }
#ifdef _PRINT
	      std::cout << "   element = this(";
	      s_1_index.PrintwithoutNewLine(std::cout);	      
	      std::cout << ") = sum_j' this(";
	      s_1_index.PrintwithoutNewLine(std::cout);	      
	      std::cout << ", j') = " << element << "\n" << std::flush;
#endif
	      // b.) calculate product(i, j) = sum_c this_t(i, c)
	      // * b_t(c, j) * element
	      // loop over all last d indices j
		
	      {  
		for (int s_2_l_index=0; s_2_l_index<s_max; s_2_l_index++)  
		  {
		    convert_to_base(s_2_l_index, length, &s_2_index);
		      
		    product_element = static_cast<T>(0);

		    for (int i=0; i<dimensions_2; i++) 
		      {
			// index = (i,j)
			index.SetElement(i,              s_1_index.GetElement(i));
			index.SetElement(i+dimensions_2, s_2_index.GetElement(i));
		      }
		     
		    // loop over d indices c

		    for (int s_3_l_index=0; s_3_l_index<s_max; s_3_l_index++) 
		      {
			convert_to_base(s_3_l_index, length, &s_3_index);

			for (int i=0; i<dimensions_2; i++) 
			  {
			    // index_1 = (i,c)
			    index_1.SetElement(i,              s_1_index.GetElement(i));
			    index_1.SetElement(i+dimensions_2, s_3_index.GetElement(i));
			    // index_2 = (c,j)
			    index_2.SetElement(i,              s_3_index.GetElement(i));
			    index_2.SetElement(i+dimensions_2, s_2_index.GetElement(i));
			  }
			  
			product_element += this_t.GetElement(index_1) * b_t.GetElement(index_2);
#ifdef _PRINT		      
			std::cout << "      product_element += this_t(";
			index_1.PrintwithoutNewLine(std::cout);		    
			std::cout << ") (" 
			     << this_t.GetElement(index_1)
			     << ") * b_t(" ;
			index_2.PrintwithoutNewLine(std::cout);		    
			std::cout << ") = " 
			     << product_element << "\n" << std::flush;
#endif
		      }
#ifdef _PRINT
		    std::cout << "   product(";
		    index.PrintwithoutNewLine(std::cout);		    		  
		    std::cout << ") = product_element (" 
			 << product_element
			 << ") * element ("
			 << element 
			 << ") = " << std::flush;
#endif
		    product_element *= element;
		    sum_of_all_elements += product_element;
		    product.SetElement(index, product_element);
#ifdef _PRINT
		    std::cout << product.GetElement(index) << "\n" << std::flush;
#endif
		  }
	      }
	    }

#ifdef _PRINT
	  std::cout << "check sum_of_all_elements for final matrix = " 
	       << sum_of_all_elements << "\n";
#endif
	  if (fabs(sum_of_all_elements - static_cast<T>(1.)) > Max_deviation)
	    {
	      std::cout << "ERROR: array<T>::SpecialMultiply: the sum of element for the final matrix is " 
		   << sum_of_all_elements 
		   << "and |sum_of_all_element - 1.| = "
		   << fabs(sum_of_all_elements -static_cast<T>(1.))
		   << " > Max_deviation = " << Max_deviation
		   << "\n" << std::flush;
	      check++;
	    }
	}
    }
  
  if (check == 0)
    {
      return(product);
    }
  else // return empty matrix if checks failed
    {
      product.~array<T>();
      array<T> product;
      return(product);
    }
}

template<class T> int array<T>::Equal( const array<T> &v) const
{
  int arrays_are_equal = 1;

  if (this == &v) {
    arrays_are_equal = 1;
  }
  else {

    // check that number of dimensions is the same in the two arrays
    
    if (this->GetNumberofDimensions() != v.GetNumberofDimensions()) {

      arrays_are_equal = 0;
#ifdef _PRINT
      std::cout << "arrays do not have the same number of dimensions => two arrays are not the same.\n";
#endif
    }
    else {

      int d = this->GetNumberofDimensions();
      int i = 0;
      int m = 1;

#ifdef _PRINT
      std::cout << "d = " << d << "\n";
#endif

      // check that length of every dimension is the same in the
      // two arrays

      for (i=0; i<d; i++) {
	if (this->GetDimension(i) != v.GetDimension(i)) {
	  arrays_are_equal = 0;
#ifdef _PRINT
	  std::cout << "dimension " << i << " dimension in this array (" 
	       << this->GetDimension(i) << ") != dimension in input array (" 
	       << v.GetDimension(i) << ") => two arrays are not the same.\n";
#endif
	  break;
	}
	else {
	  m *= this->GetDimension(i);
	}
      }

#ifdef _PRINT
      std::cout << "m = " << m << "\n";
#endif
      // check that every entry in the two arrays has the same value

      if (arrays_are_equal == 1) {
#ifdef _PRINT
	std::cout << "number of elements in array = " << m << "\n";
#endif
	for (i=0; i<m; i++) {
	  if (this->GetElement(i) != v.GetElement(i)) {
#ifdef _PRINT
	    std::cout << "element with linear index " << i 
		 << " = " << this->GetElement(i) 
		 << " in this array, but = " << v.GetElement(i)
		 << ") in input array => two arrays are not the same.\n";
#endif
	    arrays_are_equal = 0;
	    break;
	  }
#ifdef _PRINT
	  else {
	    std::cout << i << "\t" << this->GetElement(i) << " = "
		 << v.GetElement(i) << "\n";
	  }
#endif
	}
      }
    }

  }
  return(arrays_are_equal);
}


template<class T> array<T> & array<T>::operator= ( const array<T> &v)
{
  if ( this != &v )
    {
      if (data) delete [] data;
      data=NULL;
      if (dim)  delete [] dim;
      dim=NULL;

      number_of_dimensions=v.number_of_dimensions;

      // create dim and copy its elements

      if (number_of_dimensions>0)
	{
	  dim= new int[number_of_dimensions];

	  for (int i=0; i<number_of_dimensions; i++)
	    {
	      dim[i]=v.dim[i];
	    }

	  // create data and copy its elements
	  
	  int length=dim[0];
	  for (int j=1;j<number_of_dimensions; j++)
	    {length*=dim[j];}

	  if (length>0)
	    {
	      data= new T [length];	 
	      for (int i=0; i<length; i++)
		{
		  data[i]=v.data[i];
		}
	    }
	}
    }
  return(*this);
}


template<class T> void array<T>::Print ( std::ostream &o ) const
{
  PrintwithoutNewLine(o);
  o << "\n";
  return;
}                                  

template<class T> void array<T>::PrintwithoutNewLine ( std::ostream &o ) const
{
  if ((dim) && (data) && (number_of_dimensions!=0))
    {
      int d=number_of_dimensions;
      int i,j,n;
      int delta_bracket=0;
      
      int *intern_indices= new int[d+1];
      for (j=0; j<d+1;j++)
	{intern_indices[j]=0;}
      
      int length=dim[0];
      for (j=1;j<d; j++)
	{length*=dim[j];}
      
      
      for (i=0;i<d;i++)
	{o << "[";} 
      
      for (i=0;i<length;i++)
	{
	  // opening brackets [:
	  
	  for (n=0;n<delta_bracket;n++)
	    {o << "[";} 
	  
	  o << data[i];
	  
	  // closing brackets ]:
	  
	  intern_indices[d]+=1;
	  delta_bracket=0;

	  for (n=d; n>0; n--)
	    {
	      if (intern_indices[n]==dim[n-1])
		{
		  intern_indices[n-1]+=1;
		  delta_bracket+=1;
		  o << "]";
		  if (n==d-1) {o << "\n";}
		  for (j=n;j<d+1;j++)
		    {intern_indices[j]=0;}
		}
	    }
	  if (delta_bracket==0) {o << ", ";}
	}
      if (intern_indices) delete [] intern_indices;
      intern_indices = NULL;
    }
  return;
}                                  
 
template<class T> void array<T>::PrintDimensions ( std::ostream &o ) const
{
  if (dim && (number_of_dimensions>0))
    {
      for (int i=0; i<number_of_dimensions-1; i++)
	{
	  o << dim[i] << " x ";
	}
      o << dim[number_of_dimensions-1] << "\n";
    }
  else
    {
      o << "--\n";
    }
  return;
}

template<class T> void array<T>::PrintwithIndices ( std::ostream &o, const char* int_to_char_mapping) const 
{
  if ((dim) && (data) && (number_of_dimensions!=0))
    {
      int d=number_of_dimensions;
      int i,j,n;

      array<int> indices(1);
      indices.SetDimension(0,d);
      
      int *intern_indices= new int[d+1];
      for (j=0; j<d+1;j++)
	{intern_indices[j]=0;}
      
      int length=dim[0];
      for (j=1;j<d; j++)
	{length*=dim[j];}
      
      for (i=0;i<length;i++)
	{
	  for (n=0; n<d; n++)
	    {indices.SetElement(n,intern_indices[n+1]);}

	  for (n=0; n<d; n++) 
	    {o << int_to_char_mapping[indices.GetElement(n)];}
	  
	  o << " " << data[i] << "\n"; 
	  
	  // closing brackets ]:
	  
	  intern_indices[d]+=1;
	  
	  for (n=d; n>0; n--)
	    {
	      if (intern_indices[n]==dim[n-1])
		{
		  intern_indices[n-1]+=1;
		  for (j=n;j<d+1;j++)
		    {intern_indices[j]=0;}
		}
	    }
	}
      if (intern_indices) delete [] intern_indices;
      intern_indices = NULL;
    }
  return;
}

template<class T> void array<T>::PrintwithIndices ( std::ostream &o ) const
{
  if ((dim) && (data) && (number_of_dimensions!=0))
    {
      int d=number_of_dimensions;
      int i,j,n;

      array<int> indices(1);
      indices.SetDimension(0,d);
      
      int *intern_indices= new int[d+1];
      for (j=0; j<d+1;j++)
	{intern_indices[j]=0;}
      
      int length=dim[0];
      for (j=1;j<d; j++)
	{length*=dim[j];}
      
      for (i=0;i<length;i++)
	{
	  for (n=0; n<d; n++)
	    {indices.SetElement(n,intern_indices[n+1]);}
	  
	  indices.PrintwithoutNewLine(o); 
	  o << " " << data[i] << "\n";
	  
	  // closing brackets ]:
	  
	  intern_indices[d]+=1;
	  
	  for (n=d; n>0; n--)
	    {
	      if (intern_indices[n]==dim[n-1])
		{
		  intern_indices[n-1]+=1;
		  for (j=n;j<d+1;j++)
		    {intern_indices[j]=0;}
		}
	    }
	}
      if (intern_indices) delete [] intern_indices;
      intern_indices = NULL;
    }
  return;
}                                  

template<class T> void array<T>::PrintwithIndicesRestrict ( std::ostream &o, const T dontprint) const
{
  if ((dim) && (data) && (number_of_dimensions!=0))
    {
      int d=number_of_dimensions;
      int i,j,n;

      array<int> indices(1);
      indices.SetDimension(0,d);
      
      int *intern_indices= new int[d+1];
      for (j=0; j<d+1;j++)
	{intern_indices[j]=0;}
      
      int length=dim[0];
      for (j=1;j<d; j++)
	{length*=dim[j];}
      
      for (i=0;i<length;i++)
	{
	  for (n=0; n<d; n++)
	    {indices.SetElement(n,intern_indices[n+1]);}
	  
	  if (data[i] != dontprint) {

	    indices.PrintwithoutNewLine(o); 
	    o << " " << data[i] << "\n";
	  }	  

	  // closing brackets ]:
	  
	  intern_indices[d]+=1;
	  
	  for (n=d; n>0; n--)
	    {
	      if (intern_indices[n]==dim[n-1])
		{
		  intern_indices[n-1]+=1;
		  for (j=n;j<d+1;j++)
		    {intern_indices[j]=0;}
		}
	    }
	}
      if (intern_indices) delete [] intern_indices;
      intern_indices = NULL;
    }
  return;
}                                  

template<class T> void array<T>::PrintonlyNonZerowithIndices ( std::ostream &o ) const
{
  if ((dim) && (data) && (number_of_dimensions!=0))
    {
      int d=number_of_dimensions;
      int i,j,n;

      array<int> indices(1);
      indices.SetDimension(0,d);
      
      int *intern_indices= new int[d+1];
      for (j=0; j<d+1;j++)
	{intern_indices[j]=0;}
      
      int length=dim[0];
      for (j=1;j<d; j++)
	{length*=dim[j];}
      
      for (i=0;i<length;i++)
	{
	  for (n=0; n<d; n++)
	    {indices.SetElement(n,intern_indices[n+1]);}
	  
	  if (data[i] != static_cast<T>(0))
	    {
	      indices.PrintwithoutNewLine(o); 
	      o << " " << data[i] << "\n";
	    }

	  // closing brackets ]:
	  
	  intern_indices[d]+=1;
	  
	  for (n=d; n>0; n--)
	    {
	      if (intern_indices[n]==dim[n-1])
		{
		  intern_indices[n-1]+=1;
		  for (j=n;j<d+1;j++)
		    {intern_indices[j]=0;}
		}
	    }
	}
      if (intern_indices) delete [] intern_indices;
      intern_indices = NULL;
    }
  return;
}                                  

template<class T> void array<T>::PrintonlyNonLogzerowithIndices ( std::ostream &o ) const
{
  if ((dim) && (data) && (number_of_dimensions!=0))
    {
      int d=number_of_dimensions;
      int i,j,n;

      array<int> indices(1);
      indices.SetDimension(0,d);
      
      int *intern_indices= new int[d+1];
      for (j=0; j<d+1;j++)
	{intern_indices[j]=0;}
      
      int length=dim[0];
      for (j=1;j<d; j++)
	{length*=dim[j];}
      
      for (i=0;i<length;i++)
	{
	  for (n=0; n<d; n++)
	    {indices.SetElement(n,intern_indices[n+1]);}
	  
	  if (data[i] != static_cast<T>(Logzero))
	    {
	      indices.PrintwithoutNewLine(o); 
	      o << " " << data[i] << "\n";
	    }

	  // closing brackets ]:
	  
	  intern_indices[d]+=1;
	  
	  for (n=d; n>0; n--)
	    {
	      if (intern_indices[n]==dim[n-1])
		{
		  intern_indices[n-1]+=1;
		  for (j=n;j<d+1;j++)
		    {intern_indices[j]=0;}
		}
	    }
	}
      if (intern_indices) delete [] intern_indices;
      intern_indices = NULL;
    }
  return;
}                                  

template<class T> void array<T>::ResetData(void)
{
  if (data) delete [] data;
  data=NULL;

  if (dim)
    {
      int length=dim[0];
      for (int j=1;j<number_of_dimensions; j++)
	{length*=dim[j];}
      
      if (length>0)
	{
	  data= new T [length];      
	  for (int j=0;j<length;j++)
	    {data[j]=0;}
	}
    }
  return;
}

template<class T> void array<T>::ResetData(T initial_value)
{
  if (data) delete [] data;
  data=NULL;

  if (dim)
    {
      int length=dim[0];
      for (int j=1;j<number_of_dimensions; j++)
	{length*=dim[j];}
      
      if (length>0)
	{
	  data= new T [length];      
	  for (int j=0;j<length;j++)
	    {data[j]=initial_value;}
	}
    }
  return;
}


template<class T> void array<T>::ResetDataandDimensions (void)
{
  if (data) delete [] data;
  data=NULL;

  if (dim)
    {
      for (int i=0; i<number_of_dimensions; i++)
	{dim[i]=0;}
    }
  return;
}


template<class T> void array<T>::Reset(void)
{
  // has the same action as the destructor

  if (data) delete [] data;
  data=NULL;
  if (dim)  delete [] dim;
  dim=NULL;

  number_of_dimensions=0;
  return;
}

#endif
