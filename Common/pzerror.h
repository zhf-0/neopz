/**
 * @file
 * @brief Defines PZError
 */
#ifndef PZERRORH
#define PZERRORH

#include <iostream>
#include <stdlib.h>

#if defined(__GNUC__) || defined(__GNUG__)  // GCC
#include <execinfo.h>
#include <unistd.h>
#include <cstring>
#endif

/**
 * @ingroup common
 * @brief Defines the output device to error messages and the DebugStop() function.
 */
#define PZError std::cout
/**
 * @ingroup common
 * @brief Returns a message suggesting the user to put a breakpoint in
 */
#define DebugStop() DebugStopImpl(__FILE__, __LINE__)

void DebugStopImpl(const char *fileName, const std::size_t lineN);

class TPZErrorHandler {
public:

    static void init(int argc, char **argv) {
        baseProgram = std::string(argv[0]);
        is_initialized = true;
    }

#if defined(__GNUC__) || defined(__GNUG__)  // GCC

    static void printStackTrace() {
        if (is_initialized) {
            void *array[200];
            size_t size;

            // get void*'s for all entries on the stack
            size = backtrace(array, 200);

            std::string command("addr2line -e ");
            command.append(baseProgram);
            {
                char *temp = (char*) malloc(256);
                for (unsigned int i = 1; i < size; ++i) {
                    sprintf(temp, " %p", (char*) array[i]);
                    command.append(temp);
                }
                free(temp);
            }

            // print out all the frames to stderr
            FILE* pipe = popen(command.c_str(), "r");
            char *line = (char*) malloc(100 * sizeof (char));
            size_t count = 100;
            PZError << std::endl << "Stacktrace:" << std::endl;
            while (getline(&line, &count, pipe) > 5) {
                PZError << line;
            }
            PZError << std::endl;
            free(line);
            pclose(pipe);
        }
    }
#endif
private:
    static bool is_initialized;
    static std::string baseProgram;
};

#endif
