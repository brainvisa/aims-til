#ifndef TIL_PTR_H
#define TIL_PTR_H

// Includes from TIL library
#include "til/til_common.h"
#include "til/labelClasses.h"
#include "til/templateTools.h"

// Ptr and ConstPtr replace standard C pointers in order to have
// automatic deallocation (garbage collection).

// Ptr and ConstPtr are used to point at objects deriving from SmartObject.
// When pointing to a new SmartObject (i.e. via a constructor, or an operator=)
// they register to this object.
// When stopping to point to a SmartObject (i.e. via a destructor or an operator=)
// they unregister to this object. If they are the last to unregister, it means no
// pointer is pointing to this object anymore, and the object is destroyed.

// Ptr allows the user to modify the object pointed at, ConstPtr only allows
// const-compatible function and members call. Ptr should not point to a const
// object: ConstPtr should be used in that case.

// A Ptr can be used whenever a ConstPtr is needed.

//#define SMARTOBJ_DEBUG

// Namespace 

namespace til {


/// Smart pointer for constant objects deriving from SmartObject.
template < class T >
class ConstPtr : public SmartPtr_label
{

public: // typedefs

	typedef ConstPtr<T> Self;
	typedef T DataType;

public: // constructors & destructor

	/// Default constructor, pointing to NULL.
	// TODO: Is it better to initialize with default constructor of T
	// instead? It seems at least safer, because a call to a member of T
	// via operator-> will be always safe.
	ConstPtr() { m_pointer=0; }

	ConstPtr(const ConstPtr<T> &ptr) : m_pointer(ptr.m_pointer)
	{
#ifdef SMARTOBJ_DEBUG
		std::cout << "ConstPtr constr(ptr)" <<std::endl;
#endif
		m_pointer->subscribe();
	}


	// Deals with ConstPtr<T> ptr = p
	ConstPtr(const T *p) : m_pointer(const_cast<T*>(p))
	{
#ifdef SMARTOBJ_DEBUG
		std::cout << "ConstPtr constr(p)" << std::endl;
#endif
		m_pointer->subscribe();
	}

	template < class U >
	ConstPtr(const U *p) : m_pointer(const_cast<U*>(p))
	{
#ifdef SMARTOBJ_DEBUG
		std::cout << "ConstPtr constr(p)" << std::endl;
#endif
		m_pointer->subscribe();
	}

	// Destructor
	// NB: the single following 'virtual' is the only reason why sizeof(Ptr) is
	// 8 and not 4. And thus we have to use &Ptr everywhere, instead of
	// simple Ptr. Is that what we want? Is it bad to remove the virtual
	// here if no data is added to inheriting classes?
	virtual ~ConstPtr()
	{ 
		if (m_pointer)
		{
			m_pointer->unsubscribe();
			m_pointer = 0;
		}
	}

	

public: // functions

	/// Convert a ConstPrt into a const T*
	operator const T*(void) const { return m_pointer; }

	// We can write *ptr to get the object
	// Add an exception here if p==0;
	const T& operator*(void) const { return *m_pointer; }


	// We can access members via ptr->
	const T* operator->(void) const { return m_pointer; }

	// We can write ptr1 = ptr2
	ConstPtr& operator=(const ConstPtr<T> &p)
	{
		return operator=((const T*)(p));
	}

	/// Assignement to a pointer.
	/// We unsubscribe to the previous object we pointed to, if any.
	/// Note that this object will destroy itself if we were the last pointer
	/// to point to it (and its lock mode is off).
	/// Note that this can therefore be used to free an object before the end of
	/// the current block, by setting the pointer to NULL.
	ConstPtr& operator=(const T* p)
	{
		// First check that the new pointer is different from the old.
		if (m_pointer != p)
		{
			// If we were pointing to an object, unsubscribe from it.
			if (m_pointer) 
			{
				m_pointer->unsubscribe();
			}
			
			// Now point to the new object
			m_pointer = const_cast<T*>(p);
			
			// If the new object is not NULL, subscribe to it
			if (m_pointer)
			{
				m_pointer->subscribe();
			}
		}
		return *this;
	}


	// NB: *apparently* the following operators are not needed and even
	// worse should not be there. When we compare a Ptr with a T* the compiler
	// apparently cannot decide which operator use, between those and the
	// built-in one between two T*. So, without these, the Ptr are cast
	// into T* and compared as such. Fine with me, but might not always
	// be good. Dunno how to change that though.

	
	/// Checks whether two pointers point to the same object
	bool operator==(const T* p) const			{ return m_pointer == p;}
	bool operator==(const ConstPtr<T> &p) const	{ return this->operator==((const T*)p);}

	/// Checks whether two pointers do not point to the same object
	bool operator!=(const T* p) const			{ return !this->operator==(p); }
	bool operator!=(const ConstPtr<T> &p) const	{ return !this->operator==(p); }
	

	// Checks whether the pointer points to something
	// depreciated: just type p != 0, just like ordinary pointers
	// plus, isALlocated is a bad terminology: Ptr could point to
	// an image object, that is not allocated...
	// bool isAllocated() const { return m_pointer != 0; }


protected: // data

	T * m_pointer;
};



/// Smart pointer for objects that can be modified.
/// Simply overloads some of ConstPtr's operators to let us do non-const operations
/// on the object.
template < class T >
class Ptr : public ConstPtr<T>
{
public: // typedefs

	typedef Ptr<T> Self;

public: // constructors & destructor

	Ptr() : ConstPtr<T>() {};
	Ptr(T* p) : ConstPtr<T>(p) {};
	template < class U >
	Ptr(U* p) : ConstPtr<T>(p) {};

public: // operators

	// Basically the same operators than ConstPtr's, except that they allow non-const
	// operations on the data.

	operator T*(void) const { return this->m_pointer; }
	T& operator*(void) const { return *(this->m_pointer); }
	T* operator->(void) const { return this->m_pointer; }
	Ptr& operator=(T* p)
	{
		((ConstPtr<T>*)(this))->operator=(p);
		return *this;
	}
};


} // namespace


#endif

