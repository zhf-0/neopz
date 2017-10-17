/**
 * @file
 * @brief Contains the TPZSavable methods.
 */

#include "TPZSavable.h"
#include <iostream>                        // for operator<<, basic_ostream
#include <vector>                          // for allocator
#include "TPZStream.h"                     // for TPZStream
#include "pzlog.h"                         // for glogmutex, LOGPZ_ERROR


#include "pzlog.h"
#include "TPZPersistenceManager.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.saveable"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

using namespace std;

int TPZSavable::ClassId() const{
    return -1;
}

std::pair<std::string, long unsigned int> TPZSavable::Version() const {
    return std::make_pair("NeoPZ", 1);
}

//void TPZSavable::Write(TPZStream &buf, int withclassid) 
//{
//	if(withclassid) { 
//		int var = ClassId();
//		if(var == -1)
//		{
//			cout << "TPZSavable::Write with classid -1 expect trouble\n";
//		}
//		buf.Write(&var,1);
//	}
//}

void TPZSavable::Write(TPZStream &buf, int withclassid) const
{
	if(withclassid) { 
		int var = ClassId();
		if(var == -1)
		{
			cout << "TPZSavable::Write const with classid -1 expect trouble\n";
		}
		//buf.Write(&var,1);
	}
}


void TPZSavable::Read(TPZStream &buf, void *context)
{}

void TPZSavable::Register(int classid, TPZRestore_t fun) 
{
#ifndef ELLIPS
	map<int,TPZRestore_t>::iterator it;
	it = Map().find(classid);
	if(it != Map().end()) 
	{
		cout << "TPZSavable::Register duplicate classid " << classid << endl;
		return;
	}
	Map()[classid] = fun;
#endif
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZSavable::Compare(TPZSavable *copy, bool override)
{
	std::stringstream sout;
	sout << "Class id " << ClassId() << " Compare needs to be implemented";
	LOGPZ_ERROR(loggerCheck,sout.str())
	return false;
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZSavable::Compare(TPZSavable *copy, bool override) const
{
	std::stringstream sout;
	sout << "Class id " << ClassId() << " Compare needs to be implemented";
	LOGPZ_ERROR(loggerCheck,sout.str())
	return false;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZSavable>;
#endif


TPZSavable *TPZSavable::CreateInstance(const int &classId) {
    
#ifndef ELLIPS
    map<int,TPZRestore_t>::const_iterator it;
    it = Map().find(classId);
    if(it == Map().end()) {
        std::cout << "TPZSavable trying to restore unknown object " << classId << std::endl;
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " trying to restore unknown object " << classId;
#ifdef LOG4CXX
            LOGPZ_ERROR(logger,sout.str().c_str());
#endif
        }
        return nullptr;
    }
    
    TPZRestore_t fun= it->second;
    return (*fun)();
#else
    return nullptr;
#endif
}
