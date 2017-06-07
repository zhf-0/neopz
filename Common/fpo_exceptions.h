/**
 * @file
 * @brief Defines functions to help with invalid arithmetic operations.
 *
 *   This file implements functions previously implemented in http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c in order to improve portability of NeoPZ code and allow feenableexcept function to be called in a macOS environment. These functions were slightly modified.
 */
#ifndef __FPO_EXCEPTIONS_H
#define __FPO_EXCEPTIONS_H

#ifdef WIN32

#include <float.h>

#ifndef _EM_OVERFLOW
#define _EM_OVERFLOW EM_OVERFLOW
#endif//_EM_OVERFLOW

#ifndef _EM_UNDERFLOW
#define _EM_UNDERFLOW EM_UNDERFLOW
#endif//_EM_UNDERFLOW

#ifndef _EM_INVALID
#define _EM_INVALID EM_INVALID
#endif//_EM_INVALID

#ifndef _EM_ZERODIVIDE
#define _EM_ZERODIVIDE EM_ZERODIVIDE
#endif//_EM_ZERODIVIDE

#else 

#include <fenv.h>

#ifdef MACOSX

static int
feenableexcept (unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
    old_excepts;  // previous masks
    
    if ( fegetenv (&fenv) ) return -1;
    old_excepts = fenv.__control & FE_ALL_EXCEPT;
    
    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);
    
    return ( fesetenv (&fenv) ? -1 : old_excepts );
}

#endif //MACOSX
//#endif //DEBUG

#endif //WIN32


struct TExceptionManager {
	private:
#ifdef WIN32
	unsigned int fPrevConfig;
#else
    fenv_t fPrevConfig;
#endif //WIN32
public:
    TExceptionManager(){
#ifdef WIN32
        _controlfp_s(&fPrevConfig, 0, 0);//saves current state of fpu
        _controlfp(1, _EM_OVERFLOW);
        _controlfp(1, _EM_INVALID);
        _controlfp(1, _EM_ZERODIVIDE);
#else
        fegetenv(&fPrevConfig);//saves current state of fpu
        feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
#endif //WIN32
    }
    
    ~TExceptionManager(){
#ifdef WIN32
        unsigned int temp;
        _controlfp_s(&temp, fPrevConfig, _MCW_PC);//restores previous sates of fpu
#else
        fesetenv(&fPrevConfig);//restores previous sates of fpu
#endif //WIN32
    }
};


#endif //__FPO_EXCEPTIONS_H
