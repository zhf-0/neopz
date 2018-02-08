// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
#ifndef dealii_cxx14_algorithm_h
#define dealii_cxx14_algorithm_h



namespace std_cxx14
{

  template<class T>
  constexpr const T &max(const T &a, const T &b)
  {
    return (a < b) ? b : a;
  }

  template<class T, class Compare>
  constexpr const T &max(const T &a, const T &b, Compare comp)
  {
    return (comp(a, b)) ? b : a;
  }

  template<class T>
  constexpr const T &min(const T &a, const T &b)
  {
    return (b < a) ? b : a;
  }

  template<class T, class Compare>
  constexpr const T &min(const T &a, const T &b, Compare comp)
  {
    return (comp(b, a)) ? b : a;
  }

}

#endif // dealii_cxx14_algorithm_h
