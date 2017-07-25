#include "pzerror.h"

void DebugStopImpl(const char *fileName, const std::size_t lineN) {
#ifdef WIN32
    //ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif

#if defined(__GNUC__) || defined(__GNUG__)  // GCC
    TPZErrorHandler::printStackTrace();
#endif

    std::cout << "Your chance to put a breakpoint at " << fileName << ":" << lineN << "\n";
    std::bad_exception myex;
    throw myex;

}

bool TPZErrorHandler::is_initialized = false;
std::string TPZErrorHandler::baseProgram;
