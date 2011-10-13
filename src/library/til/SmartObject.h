#ifndef TIL_SMARTOBJECT_H
#define TIL_SMARTOBJECT_H

#include <iostream>

#include "til/til_common.h"


// SmartObject is a base class used for garbage collection.
// A class should derive from SmartObject if garbage collection is needed
// for this class.
// SmartObject contains nothing but a number that counts the number of 
// pointers pointing to this object. It provided also method for the pointers
// to increase and decrease this number.

// Ptr and ConstPtr classes implements this functionality.


// Namespace 

namespace til {

/// Base class for all classes needing reference counting based garbage collection.
/// Use in conjunction with Ptr and ConstPtr.
class TIL_API SmartObject
{

public: // Constructors and destructor

	/// Default constructor.
	/// By default, reference count is zero and lock is off
	SmartObject(void)
	{
#ifdef SMARTOBJ_DEBUG
		std::cout << "SmartObj constr" << std::endl;
#endif
		m_refCount = 0;
		m_lock = false;
	}

	/// Destructor
	/// Do nothing special. I mean, nothing at all.
	virtual ~SmartObject()
	{
#ifdef SMARTOBJ_DEBUG
		std::cout << "SmartObj destr" << std::endl;
#endif
	}

public: // functions

	/// Register to this object.
	/// (unfortunately, register is a reserved word :(
	/// Basically increase the reference count
	void subscribe(void)
	{
		++m_refCount;
#ifdef SMARTOBJ_DEBUG
		std::cout << "Increase to " << m_refCount << std::endl;
#endif
	}

	/// Unregister to this object.
	/// Decrease the reference count, and delete this object if it falls
	/// to zero and if the lock is off
	void unsubscribe(void)
	{
#ifdef SMARTOBJ_DEBUG
		std::cout << "Decrease to " << m_refCount - 1 << std::endl;
#endif
		if ((--m_refCount == 0) && (!m_lock)) delete this;
	}

	/// Get the number of pointers that have registered to this object
	int getReferenceCount(void)
	{
		return m_refCount;
	}

	/// Disable garbage collection.
	/// Prevents the object from being automatically deleted. Reference counting
	/// is still active though, in case one will unlock the object
	void lock() { m_lock = true; }

	/// Re-enable garbage collection.
	/// Note that if the reference count is null when re-enabling garbage collection,
	/// the object is immediately destroyed
	void unlock()
	{
		m_lock = false;
		if (m_refCount==0) delete this;
	}

private: // data
	
	int m_refCount;
	bool m_lock;
};


} // namespace


#endif

