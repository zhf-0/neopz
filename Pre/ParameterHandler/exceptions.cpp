// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "exceptions.h" //deal_ii_include
#include "logstream.h" //deal_ii_include
#include "utilities.h" //deal_ii_include
#include "mpi.h" //deal_ii_include

#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>


namespace deal_II_exceptions
{

  std::string additional_assert_output;

  void set_additional_assert_output (const char *const p)
  {
    additional_assert_output = p;
  }

  bool show_stacktrace = true;

  void suppress_stacktrace_in_exceptions ()
  {
    show_stacktrace = false;
  }

  bool abort_on_exception = true;

  void disable_abort_on_exception ()
  {
    abort_on_exception = false;
  }

}



ExceptionBase::ExceptionBase ()
  :
  file(""),
  line(0),
  function(""),
  cond(""),
  exc(""),
  stacktrace (nullptr),
  n_stacktrace_frames (0),
  what_str("")
{

}



ExceptionBase::ExceptionBase (const ExceptionBase &exc)
  :
  file(exc.file),
  line(exc.line),
  function(exc.function),
  cond(exc.cond),
  exc(exc.exc),
  stacktrace (nullptr), // don't copy stacktrace to avoid double de-allocation problem
  n_stacktrace_frames (0),
  what_str("") // don't copy the error message, it gets generated dynamically by what()
{

}



ExceptionBase::~ExceptionBase () noexcept
{
  free (stacktrace); // free(NULL) is allowed
  stacktrace = nullptr;
}



void ExceptionBase::set_fields (const char *f,
                                const int  l,
                                const char *func,
                                const char *c,
                                const char *e)
{
  file = f;
  line = l;
  function = func;
  cond = c;
  exc  = e;
}

const char *ExceptionBase::what() const noexcept
{
  // If no error c_string was generated so far, do it now:
  if (what_str == "")
    {
      generate_message();
    }

  return what_str.c_str();
}


const char *ExceptionBase::get_exc_name () const
{
  return exc;
}



void ExceptionBase::print_exc_data (std::ostream &out) const
{
  // print a header for the exception
  out << "An error occurred in line <" << line
      << "> of file <" << file
      << "> in function" << std::endl
      << "    " << function << std::endl
      << "The violated condition was: "<< std::endl
      << "    " << cond << std::endl;

  // print the way the additional information message was generated.
  // this is useful if the names of local variables appear in the
  // generation of the error message, because it allows the identification
  // of parts of the error text with what variables may have cause this
  //
  // On the other hand, this is almost never the case for ExcMessage
  // exceptions which would simply print the same text twice: once for
  // the way the message was composed, and once for the additional
  // information. Furthermore, the former of these two is often spread
  // between numerous "..."-enclosed strings that the preprocessor
  // collates into a single string, making it awkward to read. Consequently,
  // elide this text if the message was generated via an ExcMessage object
  if (std::strstr(cond, "ExcMessage") != nullptr)
    out << "The name and call sequence of the exception was:" << std::endl
        << "    " << exc  << std::endl;

  // finally print the additional information the exception provides:
  out << "Additional information: " << std::endl;
}



void ExceptionBase::print_info (std::ostream &out) const
{
  out << "    (none)" << std::endl;
}



void ExceptionBase::print_stack_trace (std::ostream &out) const
{
  if (n_stacktrace_frames == 0)
    return;

  if (deal_II_exceptions::show_stacktrace == false)
    return;

  // if there is a stackframe stored, print it
  out << std::endl;
  out << "Stacktrace:" << std::endl
      << "-----------" << std::endl;

  // print the stacktrace. first omit all those frames that have
  // ExceptionBase or deal_II_exceptions in their names, as these
  // correspond to the exception raising mechanism themselves, rather than
  // the place where the exception was triggered
  int frame = 0;
  while ((frame < n_stacktrace_frames)
         &&
         ((std::string(stacktrace[frame]).find ("ExceptionBase") != std::string::npos)
          ||
          (std::string(stacktrace[frame]).find ("deal_II_exceptions") != std::string::npos)))
    ++frame;

  // output the rest
  const unsigned int first_significant_frame = frame;
  for (; frame < n_stacktrace_frames; ++frame)
    {
      out << '#' << frame - first_significant_frame
          << "  ";

      // the stacktrace frame is actually of the format
      // "filename(functionname+offset) [address]". let's try to get the
      // mangled functionname out:
      std::string stacktrace_entry (stacktrace[frame]);
      const unsigned int pos_start = stacktrace_entry.find('('),
                         pos_end   = stacktrace_entry.find('+');
      std::string functionname = stacktrace_entry.substr (pos_start+1,
                                                          pos_end-pos_start-1);

      // demangle, and if successful replace old mangled string by
      // unmangled one (skipping address and offset). treat "main"
      // differently, since it is apparently demangled as "unsigned int"
      // for unknown reasons :-) if we can, demangle the function name
      stacktrace_entry = stacktrace_entry.substr(0, pos_start)
                         +
                         ": "
                         +
                         functionname;

      // then output what we have
      out << stacktrace_entry
          << std::endl;

      // stop if we're in main()
      if (functionname == "main")
        break;
    }
}



void ExceptionBase::generate_message () const
{
  // build up a c_string with the error message.
  // Guard this procedure with a try block, we shall not throw at this
  // place...
  try
    {
      std::ostringstream converter;

      converter << std::endl
                << "--------------------------------------------------------"
                << std::endl;

      // print out general data
      print_exc_data (converter);
      // print out exception specific data
      print_info (converter);
      print_stack_trace (converter);

      if (!deal_II_exceptions::additional_assert_output.empty())
        {
          converter << "--------------------------------------------------------"
                    << std::endl
                    << deal_II_exceptions::additional_assert_output
                    << std::endl;
        }

      converter << "--------------------------------------------------------"
                << std::endl;

      what_str = converter.str();
    }
  catch (...)
    {
      // On error, resume next. There is nothing better we can do...
      what_str = "ExceptionBase::generate_message () failed";
    }
}


namespace
{
  void internal_abort (const ExceptionBase &exc) noexcept
  {
    // first print the error
    std::cerr << exc.what() << std::endl;

    // then bail out. if in MPI mode, bring down the entire
    // house by asking the MPI system to do a best-effort
    // operation at also terminating all of the other MPI
    // processes. this is useful because if only one process
    // runs into an assertion, then that may lead to deadlocks
    // if the others don't recognize this, or at the very least
    // delay their termination until they realize that their
    // communication with the job that died times out.
    //
    // Unlike std::abort(), MPI_Abort() unfortunately doesn't break when
    // running inside a debugger like GDB, so only use this strategy if
    // absolutely necessary and inform the user how to use a debugger.
    std::abort();
  }
}

namespace deal_II_exceptions
{
  namespace internals
  {

    void issue_error_nothrow (ExceptionHandling,
                              const char       *file,
                              int               line,
                              const char       *function,
                              const char       *cond,
                              const char       *exc_name,
                              ExceptionBase     e) noexcept
    {
      // Fill the fields of the exception object
      e.set_fields (file, line, function, cond, exc_name);
      if (deal_II_exceptions::abort_on_exception)
        internal_abort(e);
      else
        {
          // We are not allowed to throw, and not allowed to abort.
          // Just print the exception name to deallog and continue normally:
          deallog << "Exception: " << e.get_exc_name() << std::endl;
          deallog << e.what() << std::endl;
        }
    }



    void abort (const ExceptionBase &exc)
    {
      if (deal_II_exceptions::abort_on_exception)
        internal_abort(exc);
      else
        {
          // We are not allowed to abort, so just throw the error:
          throw exc;
        }
    }

  } /*namespace internals*/
} /*namespace deal_II_exceptions*/
