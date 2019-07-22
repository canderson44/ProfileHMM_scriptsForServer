 /* 
   Authors: Irmtraud M Meyer and Philip Lam
   Copyright: Irmtraud M Meyer (1999-2009) and Philip Lam (2007-2009)
   License: licensed under the GNU General Public License version 3 (GPLv3)
   Purpose: the stack class
   RCS-Info: $Header: /ubc/cs/home/n/natural/cvs/HMMConverter/Stack.h,v 1.1.1.1 2008/07/13 09:39:09 natural Exp $
 */

#include <iostream>
#include <cassert>

#ifndef STACK
#define STACK

template <typename StackElement>
class Stack
{
 public:
    
    Stack(int numElements = 128);
    
    Stack(const Stack<StackElement> & original);
         
    ~Stack();
         
    const Stack<StackElement> & operator = (const Stack<StackElement> & original);
         
    bool empty() const;
    
    void push(const StackElement & value);
         
    void display(ostream & out) const;
         
    StackElement top() const;
    
    void pop();
             
 private:

    int myCapacity,
	myTop;
    StackElement * myArray;
};

#include <new>

template <typename StackElement>
Stack<StackElement>::Stack(int numElements)
{
     assert(numElements > 0);
     myCapacity = numElements;
     
     myArray = new(nothrow) StackElement[myCapacity];
     if(myArray != 0)
	 myTop = -1;
     else
     {
         cerr << "Inadequate memory to allocate stack \n"
	     " -- terminating execution\n";
         exit(1);
     }
}

template <typename StackElement>
Stack<StackElement>:: Stack(const Stack<StackElement> & original)
 : myCapacity(original.myCapacity), myTop(original.myTop)
{
  myArray = new(nothrow) StackElement[myCapacity];
  if(myArray != 0)
    for(int pos = 0; pos <= myTop; pos++)
      myArray[pos] = original.myArray[pos];
  else{
    cerr<< "*Inadequate memory to allocate stack ***\n";
    exit(1);
  }
}

template<typename StackElement>
inline Stack<StackElement>::~Stack()
{
	delete[] myArray;
}

template<typename StackElement>
const Stack<StackElement> & Stack<StackElement>::operator=(const Stack<StackElement> & rightHandSide)
{
    if(this != &rightHandSide)
    {
	if(myCapacity != rightHandSide.myCapacity)
	{
	    delete[] myArray;
	    
	    myCapacity = rightHandSide.myCapacity;
	    myArray = new StackElement[myCapacity];
	    if(myArray == 0)
	    {
		cerr << "*** Inadequate memory ***\n";
		exit(1);
	    }
	}
	
	myTop = rightHandSide.myTop;
	for(int pos = 0; pos <= myTop; pos++)
	    myArray[pos] = rightHandSide.myArray[pos];
    }
    return *this;
}

template<typename StackElement>
inline bool Stack<StackElement>::empty() const
{
    return(myTop==-1);
}

template<typename StackElement>
inline void Stack<StackElement>::push(const StackElement & value)
{
    if(myTop < myCapacity -1)
    {
	++myTop;
	myArray[myTop] = value;
    }
    else
    {
	cerr<< "*** Stack full -- can't add new value ***\n"
	    "Must increase value of STACK_CAPACITY in DoubleStack.h\n";
	exit(1);
    }
}

template<typename StackElement>
inline void Stack<StackElement>::display(ostream & out) const
{
    for(int i = myTop; i >= 0; i--)
	out<< myArray[i] << endl;
}

template<typename StackElement>
inline StackElement Stack<StackElement>::top() const
{
    if(!empty())
	return (myArray[myTop]);
    else
    {
	cerr << "*** Stack is empty -- returning garbage value ***\n";
	StackElement garbage=0; 
	return garbage;
    }
}

template<typename StackElement>
inline void Stack<StackElement>:: pop()
{
    if(!empty())
        myTop--;
    else
	cerr <<"*** Stack is empty - can't remove a value ***\n";
}
 
#endif
