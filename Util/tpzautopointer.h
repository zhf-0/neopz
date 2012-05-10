/**
 * @file
 * @brief Contains declaration of the TPZAutoPointer class which has Increment and Decrement actions are mutexed by this mutex.
 */
//
// C++ Interface: tpzautopointer
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZAUTOPOINTER_H
#define TPZAUTOPOINTER_H

#include "pzp_thread.h"

/**
 * \addtogroup util
 * @{
 */


/**
 * @brief Increment and Decrement actions are mutexed by this mutex
 */
extern pthread_mutex_t gAutoPointerMutex;

#define AP_MUTEX_ARRAY_SZ 512

#define PROFILE_AP_MUTEXES

#ifdef PROFILE_AP_MUTEXES
  extern unsigned long long ap_mutex_accesses[];
#endif

#define AP_MUTEX_HASH_1         \
  addr = (addr >> 32) ^ addr;   \
  addr = (addr >> 16) ^ addr;   \
  addr = (addr >> 8)  ^ addr;   \
  addr = (addr >> 4)  ^ addr;   \
  i = (unsigned) (addr % AP_MUTEX_ARRAY_SZ)

#define AP_MUTEX_HASH_2         \
  addr = (addr >> 8)  ^ addr;   \
  addr = (addr >> 4)  ^ addr;   \
  i = (unsigned) (addr % AP_MUTEX_ARRAY_SZ)

extern pthread_mutex_t gAutoPointerMutexArray[];
inline pthread_mutex_t* get_ap_mutex(void* obj)
{
  //return &gAutoPointerMutex; // Single mutex Approach.
  unsigned i;
  unsigned long long addr = (unsigned long long) obj;
  //  AP_MUTEX_HASH_1;
  AP_MUTEX_HASH_2;
#ifdef PROFILE_AP_MUTEXES
  ap_mutex_accesses[i]++;
#endif
  return &(gAutoPointerMutexArray[i]);
}


/**
 @brief This class implements a reference counter mechanism to administer a dynamically allocated object. \ref util "Utility"
 @author Philippe R. B. Devloo
 */
template<class T>
class TPZAutoPointer{
	
	/** @brief Counter struct */
	template<class T2>
	struct TPZReference
	{
		
		T2 *fPointer;
		int fCounter;
		
		TPZReference()
		{
			fPointer = 0;
			fCounter = 1;
		}
		
		TPZReference(T2 *pointer)
		{
			fPointer = pointer;
			fCounter = 1;
		}
		
		~TPZReference()
		{
			if(fPointer) delete fPointer;
			fPointer = 0;
		}
		
		/** @brief Increment the counter */
		void Increment()
		{
		  PZP_THREAD_MUTEX_LOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__);
		  fCounter++;
		  PZP_THREAD_MUTEX_UNLOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__);
		}
		/** @brief Decrease the counter. If the counter is zero, delete myself */
		void Decrease()
		{
			int should_delete = 0;
		        PZP_THREAD_MUTEX_LOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__);
			fCounter--;
			if(fCounter <= 0) should_delete = 1;
		        PZP_THREAD_MUTEX_UNLOCK(get_ap_mutex((void*) this), __PRETTY_FUNCTION__);
			if(should_delete) 
			{
				delete this;
			}
		}
		
	};
	
	/** @brief The object which contains the pointer and the reference count */
	TPZReference<T> *fRef;
	
public:
	/** @brief Creates an reference counted null pointer */
	TPZAutoPointer()
	{
		fRef = new TPZReference<T>();
	}
	
	/** @brief The destructor will delete the administered pointer if its reference count is zero */
	~TPZAutoPointer()
	{
		fRef->Decrease();
	}
	
	/** @brief This method will create an object which will administer the area pointed to by obj */
	TPZAutoPointer(T *obj)
	{
		fRef = new TPZReference<T>(obj);
	}
	
	/** @brief Share the pointer of the copy */
	TPZAutoPointer(const TPZAutoPointer<T> &copy)
	{
		fRef = copy.fRef;
		fRef->Increment();
	}
	/** @brief Assignment operator */
	TPZAutoPointer &operator=(const TPZAutoPointer<T> &copy)
	{
		if(copy.fRef == fRef) return *this;
		copy.fRef->Increment();
		fRef->Decrease();
		fRef = copy.fRef;
		return *this;
	}
	
	operator bool() const{
		return (this->fRef->fPointer != 0);
	}
	
	operator T&() 
	{
		return *(fRef->fPointer);
	}
	
	T *operator->()
	{
		return fRef->fPointer;
	}
	
	T *operator->() const
	{
		return fRef->fPointer;
	}
	
	operator bool() {
		return fRef->fPointer != 0;
	}
	/** @brief Returns the counter value */
	int Count()
	{
		return fRef->fCounter;
	}
	
	int Count() const
	{
		return fRef->fCounter;
	}
	
};

/** @} */

#endif
