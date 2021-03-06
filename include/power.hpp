// power.hpp: calculate x^2, x^3, etc
//
// Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef POWER_HPP
#define POWER_HPP

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// function power2, power3, power4, power6
//

#ifndef BOOST_NO_SFINAE

template<typename T>
T power2(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(std::complex<T> const& t)
{ return power2(real(t)) + power2(imag(t)); }

template<typename T>
T power3(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power3(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power4(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(std::complex<T> const& t)
{ return power2(power2(t)); }

template<typename T>
T power6(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return power3(power2(t)); }

template<typename T>
T power6(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return power3(power2(t)); }

template<typename T>
T power6(std::complex<T> const& t)
{ return power3(power2(t)); }

#else

template<typename T>
T power2(T const& t) { return t * t; }

template<typename T>
T power3(T const& t) { return t * t * t; }

template<typename T>
T power4(T const& t) { return power2(power2(t)); }

template<typename T>
T power6(T const& t) { return power3(power2(t)); }

#endif

#endif // POWER_HPP
