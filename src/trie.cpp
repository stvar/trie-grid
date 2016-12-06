// Copyright (C) 2016  Stefan Vargyas
// 
// This file is part of Trie-Grid.
// 
// Trie-Grid is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Trie-Grid is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Trie-Grid.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __GNUC__
#error we need a GCC compiler
#endif

#include <execinfo.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cerrno>
#include <cctype>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <getopt.h>
#include <unistd.h>

#include <sys/times.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <iterator>
#include <stdexcept>

#define STRINGIFY0(S) #S
#define STRINGIFY(S)  STRINGIFY0(S)

const char program[] = STRINGIFY(PROGRAM);
const char verdate[] = "v0.1 2015-09-07"; // $ date +%F

const char license[] =
"Copyright (C) 2016  Stefan Vargyas.\n"
"License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n"
"This is free software: you are free to change and redistribute it.\n"
"There is NO WARRANTY, to the extent permitted by law.\n";

struct global_options_t
{
    global_options_t() :
        dump_backtrace(false),
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        insertion_threshold(10),
#endif
        utime_prec(3)
    {}

    unsigned dump_backtrace : 1;
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
    size_t insertion_threshold;
#endif
    size_t utime_prec;
};

global_options_t globals;

#define PRINTF_FMT(F) \
    __attribute__ ((format(printf, F, F + 1)))
#define NORETURN \
    __attribute__ ((noreturn))
#define UNUSED \
    __attribute__ ((unused))

#define SYS_UNEXPECT_ERR(M, ...) \
    do { \
        Sys::unexpect_error(__FILE__, __LINE__, M, ## __VA_ARGS__); \
    } \
    while (0)

#ifdef DEBUG
# define SYS_ASSERT(E) \
    do { \
        if (!(E)) \
            Sys::assert_failed(__FILE__, __LINE__, #E); \
    } \
    while (0)
#else
# define SYS_ASSERT(E) \
    do {} while (0)
#endif

namespace Sys {

template<bool b>
struct cxx_assert_t;

template<>
struct cxx_assert_t<true>
{ enum { value }; };

#define CXX_ASSERT(E) \
    do { \
        (void) Sys::cxx_assert_t<(E)>::value; \
    } \
    while (0)

void die(const char *msg, ...)
    PRINTF_FMT(1)
    NORETURN;

void die(const char *msg, ...)
{
    char buf[256];

    va_list args;
    va_start(args, msg);
    vsnprintf(buf, sizeof buf - 1, msg, args);
    va_end(args);
    buf[255] = 0;

    using namespace std;
    cerr << program << ": fatal error: " << buf << endl;

    if (globals.dump_backtrace) {
        const size_t n = 128;
        void* b[n];
        size_t s = backtrace(b, n);
        for (void
            **p = b,
            **e = p + (s % (n + 1));
            p < e; ++ p)
            cerr << *p << endl;
    }

    exit(127);
}

void assert_failed(const char *file, int line, const char *expr)
    NORETURN;

void assert_failed(const char *file, int line, const char *expr)
{
    char buf[256];

    snprintf(buf, sizeof buf - 1, "%s:%d: %s", file, line, expr);
    buf[255] = 0;

    die("assertion failed: %s", buf);
}

void unexpect_error(const char *file, int line, const char *msg, ...)
    PRINTF_FMT(3)
    NORETURN;

void unexpect_error(const char *file, int line, const char *msg, ...)
{
    char buf[256];

    va_list args;
    va_start(args, msg);
    vsnprintf(buf, sizeof buf - 1, msg, args);
    va_end(args);
    buf[255] = 0;

    die("unexpected error:%s:%d: %s", file, line, buf);
}

template<typename T>
void verror(const char* msg, va_list args)
{
    char buf[256];
    vsnprintf(buf, sizeof buf - 1, msg, args);
    buf[255] = 0;

    throw T(buf);
}

template<typename T>
void error(const char* msg, ...)
    PRINTF_FMT(1);

template<typename T>
void error(const char* msg, ...)
{
    va_list args;
    va_start(args, msg);
    verror<T>(msg, args);
    //!!!VA_END va_end(args);
}

void error(const char* msg, ...)
    PRINTF_FMT(1);

void error(const char* msg, ...)
{
    va_list args;
    va_start(args, msg);
    verror<std::runtime_error>(msg, args);
    //!!!VA_END va_end(args);
}

void line_error(size_t ln, const char* msg, ...)
    PRINTF_FMT(2);

void line_error(size_t ln, const char* msg, ...)
{
    char buf[256];

    va_list args;
    va_start(args, msg);
    vsnprintf(buf, sizeof buf - 1, msg, args);
    va_end(args);
    buf[255] = 0;

    error<std::runtime_error>("line #%zu: %s", ln, buf);
}

template<typename T, size_t N>
class array_t
{
public:
    typedef T elem_t;
    typedef T type_t[N];

    array_t(type_t& _a) : a(_a) {}

    elem_t operator[](size_t i) const
    { SYS_ASSERT(i < N); return a[i]; }

private:
    elem_t* a;
};

template<typename T, size_t N, size_t M>
class array_t<T[M], N>
{
public:
    typedef array_t<T, M> elem_t;
    typedef T type_t[N][M];

    array_t(type_t& _a) : a(_a) {}

    elem_t operator[](size_t i) const
    { SYS_ASSERT(i < N); return elem_t(a[i]); }

private:
    typedef T inner_t[M];
    inner_t* a;
};

template<typename T, size_t N>
inline array_t<T, N> array(T (&v)[N])
{ return array_t<T, N>(v); }

template<typename T, size_t N>
inline size_t array_size(T (&v)[N] UNUSED)
{ return N; }

template<typename T, typename A>
inline T* data(std::vector<T, A>& obj)
{ return obj.size() ? &obj[0] : 0; }

template<typename T, typename A>
inline const T* data(const std::vector<T, A>& obj)
{ return obj.size() ? &obj[0] : 0; }

template <typename T>
class Repr;

template <typename T>
inline std::ostream& operator<<(
    std::ostream& ost, const Repr<T>& repr)
{ repr.print(ost); return ost; }

template<typename T>
class ReprPrint
{
public:
    explicit ReprPrint(std::ostream& _ost) :
        ost(_ost)
    {}

    void operator()(const T& val)
    { ost << val; }

private:
    std::ostream& ost;
};

struct repr_type_t
{
    enum char_t {
        quoted_char,
        plain_char
    };

    enum str_t {
        string_str,
        plain_str
    };
};

template<>
class ReprPrint<char>
{
public:
    explicit ReprPrint(std::ostream& _ost) :
        ost(_ost)
    {}

    size_t operator()(char ch, repr_type_t::str_t type)
    {
        const bool str = type == repr_type_t::string_str;

        switch (ch) {
        case '\\':
            return print_out_quoted('\\');
        case '\f':
            return print_out_quoted('f');
        case '\n':
            return print_out_quoted('n');
        case '\r':
            return print_out_quoted('r');
        case '\t':
            return print_out_quoted('t');
        case '\v':
            return print_out_quoted('v');
        case '"':
            if (str)
                return print_out_quoted('"');
            else
                return print_out(ch);
        case '\'':
            if (!str)
                return print_out_quoted('\'');
            else
                return print_out(ch);
        default:
            if (isascii(ch) && !iscntrl(ch))
                return print_out(ch);
            else
                return print_out_hex(ch, str);
        }
    }

    size_t operator()(char ch, repr_type_t::char_t type)
    {
        const bool q = type == repr_type_t::quoted_char;
        size_t n = 0;

        if (q)
            n += print_out('\'');

        n += operator()(ch, q
            ? repr_type_t::plain_str
            : repr_type_t::string_str);

        if (q)
            n += print_out('\'');

        return n;
    }

private:
    size_t print_out(char ch)
    { ost << ch; return 1; }

    size_t print_out_quoted(char ch)
    { ost << '\\' << ch; return 2; }

    size_t print_out_hex(char ch, bool str)
    {
        using namespace std;
        char f = ost.fill();
        streamsize w = ost.width();
        ios_base::fmtflags m = ost.flags();
        unsigned char c = static_cast<unsigned char>(ch);
        ost << "\\x"
            << hex
            << right
            << setfill('0')
            << setw(1 + str)
            << static_cast<unsigned>(c)
            << setfill(f)
            << setw(w);
        ost.flags(m);
        return 3 + (str || c >= '\x10');
    }

    std::ostream& ost;
};

template<>
struct ReprPrint<const char> : public ReprPrint<char>
{
    typedef ReprPrint<char> base_t;

    explicit ReprPrint(std::ostream& _ost) :
        base_t(_ost)
    {}
};

template<typename T>
class ReprPrint<T*>
{
public:
    explicit ReprPrint(std::ostream& _ost) :
        ost(_ost)
    {}

    size_t operator()(T* beg, T* end, bool str)
    {
        if (str)
            ost << '"';
        printer_t r = std::for_each(
            beg, end, printer_t(ost, str));
        if (str)
            ost << '"';
        return str ? r.len + 2 : r.len;
    }

private:
    struct printer_t
    {
        explicit printer_t(std::ostream& _ost, bool _str) :
            repr(_ost), len(0), str(_str)
        {}

        void operator()(char ch)
        { len += repr(ch, str
                    ? repr_type_t::string_str
                    : repr_type_t::plain_str); }

        ReprPrint<T> repr;
        size_t len;
        bool str;
    };

    std::ostream& ost;
};

template<>
class Repr<char>
{
public:
    Repr(char _ch, repr_type_t::char_t _type) :
        ch(_ch), type(_type)
    {}

    void print(std::ostream& ost) const
    { ReprPrint<char> r(ost); r(ch, type); }

private:
    char                ch;
    repr_type_t::char_t type;
};

template<typename T>
class Repr<T*>
{
public:
    Repr(const T* _beg, const T* _end, bool _str) :
        beg(_beg), end(_end), str(_str)
    {}

    Repr(const T* _beg, size_t _size, bool _str) :
        beg(_beg), end(_beg + _size), str(_str)
    {}

    void print(std::ostream& ost) const
    { ReprPrint<const T*> r(ost); r(beg, end, str); }

private:
    const T *beg;
    const T *end;
    bool     str;
};

template<typename T>
inline Repr<T*> repr(T* beg, T* end, bool str)
{ return Repr<T*>(beg, end, str); }

template<typename T>
inline Repr<T*> repr(T* ptr, size_t size, bool str)
{ return Repr<T*>(ptr, size, str); }

inline Repr<char> repr(char ch, bool quote = true)
{ return Repr<char>(
    ch, quote
        ? repr_type_t::quoted_char
        : repr_type_t::plain_char); }

inline Repr<const char*> repr(const char* ptr, bool str = true)
{ return Repr<const char*>(ptr, ptr ? strlen(ptr) : 0, str); }

inline Repr<const char*> repr(
    const std::vector<char>& obj, bool str = true)
{ return Repr<const char*>(data(obj), data(obj) + obj.size(), str); }

inline std::string repr_str(
    const char* beg, const char* end, bool str = true)
{
    std::ostringstream ost;
    ost << repr(beg, end, str);
    return ost.str();
}

inline std::string repr_str(
    const std::vector<char>& obj, bool str = true)
{ return repr_str(data(obj), data(obj) + obj.size(), str); }

struct eq : std::binary_function<const char*, const char*, bool>
{
    bool operator()(const char* key1, const char* key2) const
    { return strcmp(key1, key2) == 0; }
};

} // namespace Sys

namespace Ext {

template<bool, typename>
struct enable_if_t;

template<typename T>
struct enable_if_t<true, T>
{ typedef T type_t; };

template<bool, typename, typename>
struct select_t;

template<typename T, typename U>
struct select_t<false, T, U>
{ typedef T type_t; };

template<typename T, typename U>
struct select_t<true, T, U>
{ typedef U type_t; };

template<typename T>
struct num_traits_t :
    private std::numeric_limits<T>
{
    typedef std::numeric_limits<T> base_t;

    using base_t::is_integer;
    using base_t::is_signed;
    using base_t::digits10;
    using base_t::digits;
};

static const size_t size_max = SIZE_MAX;

// stev: partial implementation of size_cast: only for
// integer types of narrower or equal width than size_t
template<typename T>
inline typename enable_if_t<
    num_traits_t<T>::is_integer &&
    num_traits_t<T>::is_signed,
size_t>::type_t
    size_cast(T v)
{
    CXX_ASSERT(
        num_traits_t<T>::digits <=
        num_traits_t<size_t>::digits &&
        !num_traits_t<size_t>::is_signed);
    SYS_ASSERT(v >= 0);
    return static_cast<size_t>(v);
}

template<typename T>
inline typename enable_if_t<
    num_traits_t<T>::is_integer &&
    !num_traits_t<T>::is_signed,
size_t>::type_t
    size_cast(T v)
{
    CXX_ASSERT(
        num_traits_t<T>::digits <=
        num_traits_t<size_t>::digits &&
        !num_traits_t<size_t>::is_signed);
    return static_cast<size_t>(v);
}

template<typename T>
inline typename enable_if_t<
    num_traits_t<T>::is_integer &&
    !num_traits_t<T>::is_signed,
size_t>::type_t
    digits10(T v)
{
    size_t r = 0;
    do {
        v /= 10;
        r ++;
    }
    while (v > 0);
    return r;
}

inline size_t pow10(size_t n)
{
    SYS_ASSERT(
        size_cast(num_traits_t<size_t>::digits10) >=
        n);
    size_t r = 1;
    while (n --) r *= 10;
    return r;
}

bool parse_size_num(const char* str, size_t& res)
{
    if (!*str || !isdigit(*str))
        return false;
    errno = 0;
    char* e = 0;
    unsigned long long
        r = strtoull(str, &e, 10);
    if (errno || *e || r > size_max)
        return false;
    res = static_cast<size_t>(r);
    return true;
}

template<typename T>
inline uintptr_t ptr_to_int(const T* v)
{
    // stev: converting a pointer to an integer or vice versa:
    // see http://gcc.gnu.org/onlinedocs/gcc-4.0.0/gcc/
    // Arrays-and-pointers-implementation.html
    CXX_ASSERT(
        num_traits_t<uintptr_t>::is_integer &&
        !num_traits_t<uintptr_t>::is_signed &&
        sizeof(T*) == sizeof(uintptr_t));
    return reinterpret_cast<uintptr_t>(v);
}

// stev: straightforward replacement of std::transform which
// guarantees in-order application of the unary operation
template<
    typename InputIter, typename OutputIter, typename UnaryOp>
OutputIter transform(
    InputIter first, InputIter last, OutputIter out, UnaryOp op)
{
    while (first != last)
        *out ++ = op(*first ++);
    return out;
}

template<typename ForwardIter>
ForwardIter is_sorted_until(ForwardIter first, ForwardIter last)
{
    if (first != last) {
        ForwardIter next = first;
        while (++ next != last) {
            if (*next < *first)
                return next;
            first = next;
        }
    }
    return last;
}

using Sys::data;

template<
    typename Char, typename Alloc>
inline size_t offset(
    const std::vector<Char, Alloc>& obj, const Char* ptr)
{ return size_cast(ptr - data(obj)); }

template<
    typename Char, typename Traits, typename Alloc>
std::basic_istream<Char, Traits>& getline(
    std::basic_istream<Char, Traits>& ist,
    std::vector<Char, Alloc>& str,
    Char delim)
{
    typedef Char char_t;

    size_t sz = str.capacity();
    SYS_ASSERT(sz <= str.max_size());
    // stev: enlarge 'str' to full capacity
    str.resize(sz);

#if defined(DEBUG) && \
    defined(DEBUG_GETLINE)
    std::cerr
        << "!!! > !!! sz=" << sz
        << " cap=" << str.capacity()
        << " [" << str.size()
        << "] " << Sys::repr(str)
        << std::endl;
#endif

    char_t *ptr = data(str);
    while (true) {
        // stev: TODO: what if we just doubled 'str's size
        // and the immediately following call to 'getline'
        // returns EOF?
        // Solution to this issue: we should allow calling
        // 'getline' with its 'count' arg set to '0' and if
        // not reached EOF then call 'getline' again with
        // 'count' > 0!
#if defined(DEBUG) && \
    defined(DEBUG_GETLINE)
        std::cerr
            << "!!! g !!! sz=" << sz
            << " ptr=" << offset(str, ptr)
            << " [" << str.size()
            << "] " << Sys::repr(str)
            << std::endl;
#endif
        ist.getline(ptr, sz, delim);
        size_t n = size_cast(ist.gcount());
        // (1) eof:  partial final line
        // (2) fail: partial long line
        // (3) else: complete line ending in '\n'
        const bool eof = ist.eof();
        const bool fail = ist.fail();
#if defined(DEBUG) && \
    defined(DEBUG_GETLINE)
        std::cerr
            << "!!! r !!! n=" << n
            << " eof=" << eof
            << " fail=" << fail
            << " [" << str.size()
            << "] " << Sys::repr(str)
            << std::endl;
#endif
        if (!eof && !fail) {
            // stev: account for '\n'
            SYS_ASSERT(n > 0);
            n --;
        }
        SYS_ASSERT(n <= sz);
        ptr += n;
        sz -= n;
        if (eof || !fail)
            break;
        // stev: !eof && fail
        ist.clear(ist.rdstate() & ~std::ios::failbit);
#if defined(DEBUG) && \
    defined(DEBUG_GETLINE)
        std::cerr
            << "!!! t !!! n=" << n
            << " sz=" << sz
            << " ptr=" << offset(str, ptr)
            << " [" << str.size()
            << "] " << Sys::repr(str)
            << std::endl;
#endif
        if (sz > 1)
            continue;
        sz = str.size();
        if (sz == 0)
            sz = 1;
        SYS_ASSERT(
            sz <= size_max / 2);
        SYS_ASSERT(
            sz <= str.max_size() / 2);
        sz *= 2;
        size_t k = offset(str, ptr);
        str.resize(sz);
        ptr = data(str) + k;
        sz -= k;
    }
    sz = offset(str, ptr);
    SYS_ASSERT(sz <= str.size());
    str.resize(sz);

#if defined(DEBUG) && \
    defined(DEBUG_GETLINE)
    std::cerr
        << "!!! < !!! sz=" << sz
        << " cap=" << str.capacity()
        << " [" << str.size()
        << "] " << Sys::repr(str)
        << std::endl;
#endif

    return ist;
}

template<
    typename Char, typename Traits, typename Alloc>
inline std::basic_istream<Char, Traits>& getline(
    std::basic_istream<Char, Traits>& ist,
    std::vector<Char, Alloc>& str)
{ return getline(ist, str, ist.widen('\n')); }

bool getline(FILE* file, std::vector<char>& str)
{
    struct ptr_t
    {
        ptr_t(size_t s) :
            ptr(s ? static_cast<char*>(malloc(s)) : 0)
        {}
        ~ptr_t() { if (ptr) free(ptr); }

        char*  operator()()         { return ptr; }
        char** operator&()          { return &ptr; }
        char   operator[](size_t k) { return ptr[k]; }

        char* ptr;
    };

    size_t n = 0;
    ptr_t ptr(str.capacity());
    ssize_t r = getline(&ptr, &n, file);
    if (r == -1)
        return false;

    size_t k = size_cast(r);
    SYS_ASSERT(k <= n);
    if (k > 0 && ptr[k - 1] == '\n')
        k --;
    str.assign(ptr(), ptr() + k);

    return true;
}

bool cxx_iostreams_getline(std::vector<char>& str)
{ return getline(std::cin, str); }

bool c_stdio_getline(std::vector<char>& str)
{ return getline(stdin, str); }

} // namespace Ext

namespace Sys {

struct clocks_t
{
public:
    typedef std::pair<size_t, size_t> secs_t;

    clocks_t(const char* _name = 0, size_t _prec = 0) :
        real(std::make_pair(0, 0)),
        user(std::make_pair(0, 0)),
        sys (std::make_pair(0, 0)),
        name(_name),
        prec(_prec)
    {}
    clocks_t(
        secs_t _real,
        secs_t _user,
        secs_t _sys,
        const char* _name = 0,
        size_t _prec = 0) :
        real(_real),
        user(_user),
        sys(_sys),
        name(_name),
        prec(_prec)
    {}

    void print(std::ostream& ost) const
    {
        print(ost, "real", real);
        print(ost, "user", user);
        print(ost, "sys",  sys);
        ost << std::flush;
    }

    secs_t real;
    secs_t user;
    secs_t sys;

    const char* const name;
    const size_t prec;

private:
    void print(
        std::ostream& ost, const char* n, secs_t s) const
    {
        using namespace std;
        char c = ost.fill();
        ost << (name ? name : "")
            << (name ? ": " : "")
            << n
            << setw(6 - strlen(n)) << left << ": "
            << s.first << '.'
            << setw(prec) << setfill('0') << right
            << s.second << setfill(c) << '\n';
    }
};

class ctime_t
{
public:
    ctime_t()
    {}

    clocks_t operator()(const char* n = 0) const;
    clocks_t operator()(const ctime_t& t, const char* n = 0) const;

private:
    struct tms_t : ::tms
    {
        tms_t();
        clock_t tms_rtime;
    };

    struct clocks_impl_t : clocks_t
    {
        clocks_impl_t(
            const tms_t&, const tms_t&,
            const char* = 0);
        static secs_t secs(clock_t);
        static size_t clk_tck();
    };

    tms_t tms;
};

inline clocks_t ctime_t::operator()(const char* n) const
{ return clocks_impl_t(tms, tms_t(), n); }

inline clocks_t ctime_t::operator()(const ctime_t& t, const char* n) const
{ return clocks_impl_t(tms, t.tms, n); }

class utime_t
{
public:
    utime_t() :
        prec(globals.utime_prec)
    { SYS_ASSERT(prec >= 3 && prec <= 6); }

    clocks_t operator()(const char* n = 0) const;
    clocks_t operator()(const utime_t& t, const char* n = 0) const;

private:
    struct timeval_t : timeval
    {
        timeval_t();
        timeval_t(const timeval&);
        timeval_t(time_t, suseconds_t);
        void operator=(const timeval&);
        timeval_t operator-(const timeval_t&) const;
        timeval_t operator+(const timeval_t&) const;
    };
    struct tms_t
    {
        tms_t();
        timeval_t tms_rtime;
        timeval_t tms_utime;
        timeval_t tms_stime;
        timeval_t tms_cutime;
        timeval_t tms_cstime;
    };

    struct clocks_impl_t : clocks_t
    {
        clocks_impl_t(
            const tms_t&, const tms_t&,
            const char* n = 0,
            size_t p = 0);
        static secs_t secs(const timeval_t&, size_t);
    };

    size_t prec;
    tms_t  tms;
};

inline clocks_t utime_t::operator()(const char* n) const
{ return clocks_impl_t(tms, tms_t(), n, prec); }

inline clocks_t utime_t::operator()(const utime_t& t, const char* n) const
{ return clocks_impl_t(tms, t.tms, n, prec); }

inline std::ostream& operator<<(
    std::ostream& ost, const clocks_t& obj)
{
    obj.print(ost);
    return ost;
}

// stev: below we need to make computations of the memory alignment --
// we rely on '__alignof__' -- a built-in non-ISO extension of 'gcc'
// which let inquire for the minimum alignment needed by a type (refer
// to http://gcc.gnu.org/onlinedocs/gcc/Alignment.html).

class time_t
{
public:
    enum type_t {
        none,
        ctime,
        utime,
    };

    time_t(type_t type) : obj(type) {}

    operator bool() const
    { return obj != none; }

    clocks_t operator()(const char* n = 0) const
    { return obj->clocks(n); }

    type_t type() const
    { return obj; }

    const char* name() const
    { return array(types)[obj]; }

    static const char* name(type_t type)
    { return array(types)[type]; }

    static bool type(const char* name, type_t& result);

private:
    time_t(const time_t&);
    time_t& operator=(const time_t&);

    struct base_t
    {
        virtual ~base_t() {}
        virtual clocks_t clocks(const char*) const = 0;
    };
    struct none_proxy_t : base_t
    {
        clocks_t clocks(const char* n) const
        { return clocks_t(n); }
    };
    template<typename T>
    struct proxy_t : base_t, T
    {
        typedef T type_t;

        clocks_t clocks(const char* n) const
        { return (*this)(n); }
    };
    typedef proxy_t<ctime_t> ctime_proxy_t;
    typedef proxy_t<utime_t> utime_proxy_t;

    template<size_t V0, size_t V1>
    struct max_t
    { enum { value = V0 > V1 ? V0 : V1 }; };

    template<size_t V0, size_t V1, size_t V2>
    struct max3_t
    {
        enum {
            value = max_t<
                max_t<V0, V1>::value, V2>::value
        };
    };
    template<typename T0, typename T1, typename T2>
    struct max_alignof_t
    {
        enum {
            value = max3_t<
                __alignof__(T0),
                __alignof__(T1),
                __alignof__(T2)>::value
        };
    };
    template<typename T0, typename T1, typename T2>
    struct max_sizeof_t
    {
        enum {
            value = max3_t<
                sizeof(T0),
                sizeof(T1),
                sizeof(T2)>::value
        };
    };

    class obj_t
    {
    public:
        obj_t(type_t _type) throw();
        ~obj_t();

        operator type_t()    const { return type; }
        base_t* operator->() const { return ptr; }

    private:
        enum {
            buf_size =
                max_alignof_t<
                    none_proxy_t,
                    ctime_proxy_t,
                    utime_proxy_t>::value +
                max_sizeof_t<
                    none_proxy_t,
                    ctime_proxy_t,
                    utime_proxy_t>::value
        };
        typedef char buf_t[buf_size];

        template<typename T>
        base_t* make() throw();

        type_t  type;
        base_t *ptr;
        buf_t   buf;
    };

    static const char* const types[3];

    obj_t obj;
};

ctime_t::tms_t::tms_t()
{
    tms_rtime = times(this);
}

ctime_t::clocks_impl_t::clocks_impl_t(
    const tms_t& b, const tms_t& a, const char* n) :
    clocks_t(
        secs((a.tms_rtime  - b.tms_rtime)),
        secs((a.tms_utime  - b.tms_utime) +
             (a.tms_cutime - b.tms_cutime)),
        secs((a.tms_stime  - b.tms_stime) +
             (a.tms_cstime - b.tms_cstime)),
        n, 3)
{}

inline clocks_t::secs_t ctime_t::clocks_impl_t::secs(clock_t t)
{
    // stev: bash-3.2/lib/sh/clock.c:
    secs_t r;
    size_t c = clk_tck();
    r.first = Ext::size_cast(t / c);
    r.second = Ext::size_cast(t % c);
    r.second = (r.second * 1000) / c;
    if (r.second >= 1000) {
        r.second -= 1000;
        r.first ++;
    }
    return r;
}

size_t ctime_t::clocks_impl_t::clk_tck()
{
    static size_t val = 0;
    if (val == 0)
        val = Ext::size_cast(
            sysconf(_SC_CLK_TCK));
    return val;
}

inline utime_t::timeval_t::timeval_t()
{}

inline utime_t::timeval_t::timeval_t(const timeval& t)
{
    tv_sec = t.tv_sec;
    tv_usec = t.tv_usec;
}

inline utime_t::timeval_t::timeval_t(::time_t sec, suseconds_t usec)
{
    tv_sec = sec;
    tv_usec = usec;
}

inline void utime_t::timeval_t::operator=(const timeval& t)
{
    tv_sec = t.tv_sec;
    tv_usec = t.tv_usec;
}

inline utime_t::timeval_t utime_t::timeval_t::operator-(
    const timeval_t& t) const
{
    // stev: bash-3.2/lib/sh/timeval.c:
    timeval_t r(
        tv_sec  - t.tv_sec,
        tv_usec - t.tv_usec);
    if (r.tv_usec < 0) {
        r.tv_usec += 1000000;
        r.tv_sec --;
        if (r.tv_sec < 0) {
            r.tv_sec = 0;
            r.tv_usec = 0;
        }
    }
    return r;
}

inline utime_t::timeval_t utime_t::timeval_t::operator+(
    const timeval_t& t) const
{
    // stev: bash-3.2/lib/sh/timeval.c:
    timeval_t r(
        tv_sec  + t.tv_sec,
        tv_usec + t.tv_usec);
    if (r.tv_usec >= 1000000) {
        r.tv_usec -= 1000000;
        r.tv_sec ++;
    }
    return r;
}

utime_t::tms_t::tms_t()
{
    gettimeofday(&tms_rtime, 0);
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    tms_utime = r.ru_utime;
    tms_stime = r.ru_stime;
    getrusage(RUSAGE_CHILDREN, &r);
    tms_cutime = r.ru_utime;
    tms_cstime = r.ru_stime;
}

utime_t::clocks_impl_t::clocks_impl_t(
    const tms_t& b, const tms_t& a, const char* n, size_t p) :
    clocks_t(
        secs((a.tms_rtime  - b.tms_rtime), p),
        secs((a.tms_utime  - b.tms_utime) +
             (a.tms_cutime - b.tms_cutime), p),
        secs((a.tms_stime  - b.tms_stime) +
             (a.tms_cstime - b.tms_cstime), p),
        n, p)
{}

inline clocks_t::secs_t utime_t::clocks_impl_t::secs(
    const timeval_t& t, size_t prec)
{
    // stev: bash-3.2/lib/sh/timeval.c:
    secs_t r;
    r.first = Ext::size_cast(t.tv_sec);
    r.second = Ext::size_cast(t.tv_usec % 1000000);
    SYS_ASSERT(prec <= 6);
    size_t p = Ext::pow10(prec);
    size_t s = r.second % p;
    // r.second < 1000000 && p <= 1000000
    // => r.second * p < 1000000 * 1000000
    CXX_ASSERT(
        1000000U < Ext::size_max / 1000000U);
    // => r.second * p < SIZE_MAX
    r.second = (r.second * p) / 1000000;
    if (s >= p / 2)
        r.second ++;
    if (r.second >= p) {
        r.second -= p;
        r.first ++;
    }
    return r;
}

const char* const time_t::types[3] = {
    "none",  // type_t::none
    "ctime", // type_t::ctime
    "utime", // type_t::utime
};

// stev: take note to the fact that the constructor
// of time_t::obj_t *must* have an empty throw spec
// because the class is actually a resource holder:
// it must call the destructor of the object pointed
// to by 'ptr' -- object which it itself created by
// the function 'make' below. The same restriction
// applies to 'make' function too.

time_t::obj_t::obj_t(type_t _type) throw() :
    type(_type), ptr(0)
{
    if (type == none)
        ptr = make<none_proxy_t>();
    else
    if (type == ctime)
        ptr = make<ctime_proxy_t>();
    else
    if (type == utime)
        ptr = make<utime_proxy_t>();
    else
        SYS_UNEXPECT_ERR("unknown type='%d'", type);
}

time_t::obj_t::~obj_t()
{
    if (ptr) ptr->~base_t();
}

bool time_t::type(const char* name, type_t& result)
{
    using namespace std;
    char const* const *b = types;
    char const* const *e = b + array_size(types);
    char const* const *p = find_if(b, e, bind1st(eq(), name));
    SYS_ASSERT(p >= b);
    SYS_ASSERT(p <= e);
    if (p < e) {
        result = static_cast<type_t>(
            Ext::size_cast(p - b));
        return true;
    }
    return false;
}

template<typename T>
time_t::base_t* time_t::obj_t::make() throw()
{
    const size_t a = __alignof__(T);
    const size_t s = sizeof(T);
    CXX_ASSERT(a + s <= buf_size);
    char *p = buf;
    uintptr_t r = Ext::ptr_to_int(p) % a;
    if (r) p += a - r;
    return new (p) T();
}

} // namespace Sys

namespace Trie {

template<typename T>
struct print_t;

template<typename T>
inline print_t<T> print(const T& v, bool quote)
{ return print_t<T>(v, quote); }

template<typename T>
inline std::ostream& operator<<(
    std::ostream& ost, const print_t<T>& obj)
{ obj.print(ost); return ost; }

struct range_t
{
    range_t(size_t _beg, size_t _end) :
        beg(_beg), end(_end)
    { SYS_ASSERT(beg <= end); }

    range_t(size_t _beg) :
        beg(_beg), end(_beg)
    {}

    range_t() :
        beg(0), end(0)
    {}

    range_t& operator+=(size_t sz)
    {
        SYS_ASSERT(beg <= end);
        SYS_ASSERT(beg < Ext::size_max - sz);
        SYS_ASSERT(beg + sz <= end);
        beg += sz;
        return *this;
    }

    range_t operator+(size_t sz) const
    { range_t r = *this; r += sz; return r; }

    range_t intersect(const range_t& rng) const
    { return range_t(
        rng.beg >= beg ? rng.beg : beg,
        rng.end <= end ? rng.end : end); }

    size_t size() const
    { SYS_ASSERT(beg <= end); return end - beg; }

    void print(std::ostream& ost) const
    { ost << '[' << beg << ',' << end << ']'; }

    size_t beg;
    size_t end;
};

inline std::ostream& operator<<(
    std::ostream& ost, const range_t& obj)
{ obj.print(ost); return ost; }

using Sys::data;

typedef std::vector<char> str_t;

inline bool operator<(const str_t& a, const str_t& b)
{
    const size_t n = a.size(), m = b.size();
    if (int r = memcmp(data(a), data(b), n < m ? n : m))
        return r < 0;
    else
        return n < m;
}

struct str_cmp_t
{
    bool operator()(const str_t& a, const str_t& b) const
    { return a < b; }
};

template<>
struct print_t<str_t>
{
    print_t(const str_t& _str, bool _quote) :
        str(_str), quote(_quote)
    {}

    void print(std::ostream& ost) const
    {
        if (quote)
            ost << Sys::repr(str, false);
        else
        if (const char* ptr = data(str))
            ost.write(ptr, str.size());
    }

    const str_t& str;
    bool quote;
};

class repr_t
{
public:
    repr_t(
        const char* _beg, const char* _end, range_t _rng,
        bool _ch = false) :
        beg(_beg), end(_end), rng(_rng), ch(_ch)
    {}

    repr_t intersect(const range_t& rng, bool ch = false) const
    { return repr_t(
        this->beg,
        this->end,
        this->rng.intersect(rng),
        this->ch || ch); }

    range_t range() const
    { return rng; }

    void print(std::ostream& ost) const
    {
        SYS_ASSERT(beg <= end);

        const char* b = beg + rng.beg;
        const char* e = beg + rng.end;

        if (e > end)
            e = end;
        if (b == e - 1 && ch)
            ost << Sys::repr(*b);
        else
        if (b < e)
            ost << Sys::repr(b, e, true);
    }

private:
    const char* const beg;
    const char* const end;
    const range_t rng;
    const bool ch;
};

inline std::ostream& operator<<(
    std::ostream& ost, const repr_t& obj)
{ obj.print(ost); return ost; }

inline repr_t repr(const str_t& str, range_t rng)
{ return repr_t(data(str), data(str) + str.size(), rng); }

inline repr_t repr(const str_t& str)
{ return repr(str, range_t(0, str.size())); }

struct indent_t
{
    indent_t(size_t _level, bool _dots) :
        level(_level), func(_dots ? dots : spaces)
    {}

    void print(std::ostream& ost) const
    { func(ost, level); }

    static void dots(std::ostream& ost, size_t level)
    {
        while (level --)
            ost << std::setw(1) << '.'
                << std::setw(3) << ' ';
    }

    static void spaces(std::ostream& ost, size_t level)
    {
        if (level)
            ost << std::setw(4 * level) << ' ';
    }

    typedef void (*func_t)(std::ostream&, size_t);

    size_t level;
    func_t func;
};

inline std::ostream& operator<<(
    std::ostream& ost, const indent_t& obj)
{ obj.print(ost); return ost; }

typedef std::vector<str_t> grid_t;

inline str_t& grow(grid_t& grid)
{
    using namespace Ext;
    size_t s = size_cast(grid.size());
#if defined(CONFIG_GRID_GROW_RESERVE)
    size_t c = size_cast(grid.capacity());
    SYS_ASSERT(s <= c);
    if (s == c) {
        if (c == 0)
            c = 1;
        SYS_ASSERT(c <= size_max / 2);
        SYS_ASSERT(c <= grid.max_size() / 2);
        grid.reserve(2 * c);
    }
#endif
    SYS_ASSERT(s <= size_max - 1);
    grid.resize(s + 1);
    return grid[s];
}

inline void shrink(grid_t& grid)
{
    size_t s = Ext::size_cast(grid.size());
    SYS_ASSERT(s >= 1);
    grid.resize(s - 1);
}

#if defined(CONFIG_GRID_APPEND_PUSH_BACK) && \
    defined(CONFIG_GRID_APPEND_EMPLACE_BACK)
#error parameters \
    CONFIG_GRID_APPEND_PUSH_BACK and \
    CONFIG_GRID_APPEND_EMPLACE_BACK are exclusive
#elif defined(CONFIG_GRID_APPEND_EMPLACE_BACK)
#define GRID_APPEND(grid, str) \
    grid.emplace_back(str)
#else
#define GRID_APPEND(grid, str) \
    grid.push_back(str)
#endif

typedef bool (*getline_func_t)(str_t&);

template<
    getline_func_t getline>
size_t input_append(
    grid_t& grid, size_t buf_size)
{
    str_t  str;
    SYS_ASSERT(
        buf_size <=
        str.max_size());
    str.reserve(buf_size);
    while (getline(str)) {
        if (str.size() == 0)
            return grid.size() + 1;
        GRID_APPEND(grid, str);
    }
    return 0;
}

template<
    getline_func_t getline>
size_t input_grow(
    grid_t& grid, size_t buf_size)
{
    while (true) {
        str_t& str = grow(grid);
        SYS_ASSERT(
            buf_size <=
            str.max_size());
        str.reserve(buf_size);
        if (!getline(str))
            break;
        if (str.size() == 0)
            return grid.size();
    }
    shrink(grid);
    return 0;
}

enum input_algo_t {
    append_grid,
    grow_grid
};

enum input_lib_t {
    cxx_iostreams,
    c_stdio
};

inline size_t input(
    grid_t& grid, input_algo_t algo, input_lib_t lib, size_t buf_size)
{
    static size_t (*const algos[][2])(grid_t&, size_t) = {{
    //  input_algo_t::append_grid
        input_append<Ext::cxx_iostreams_getline>, // input_lib_t::cxx_iostreams
        input_append<Ext::c_stdio_getline>,       // input_lib_t::c_stdio
    },{
    //  input_algo_t::grow_grid
        input_grow<Ext::cxx_iostreams_getline>,   // input_lib_t::cxx_iostreams
        input_grow<Ext::c_stdio_getline>,         // input_lib_t::c_stdio
    }};
    return Sys::array(algos)[algo][lib](grid, buf_size);
}

template<>
struct print_t<grid_t>
{
    print_t(const grid_t& _grid, bool _quote) :
        grid(_grid), quote(_quote)
    {}

    void print(std::ostream& ost) const
    {
        grid_t::const_iterator
            ptr = grid.begin(),
            end = grid.end();
        while (ptr < end)
            ost << Trie::print(*ptr ++, quote)
                << '\n';
    }

    const grid_t& grid;
    bool quote;
};

struct grid_result_t
{
    grid_result_t(const str_t* _str, size_t _line) :
        res(true), str(_str), line(_line)
    {}
    grid_result_t() :
        res(false), str(0), line(0)
    {}

    operator bool() const
    { return res; }

    bool res;
    const str_t* str;
    size_t line;
};

grid_result_t is_not_sorted(const grid_t& grid)
{
    grid_t::const_iterator
        beg = grid.begin(),
        end = grid.end();
    grid_t::const_iterator
        last = Ext::is_sorted_until(beg, end);
    SYS_ASSERT(last >= beg);
    SYS_ASSERT(last <= end);
    if (last != end)
        return grid_result_t(
            &*last, Ext::size_cast(last - beg) + 1);
    else
        return grid_result_t();
}

grid_result_t is_not_unique(const grid_t& grid)
{
    grid_t::const_iterator
        beg = grid.begin(),
        end = grid.end();
    grid_t::const_iterator
        last = adjacent_find(beg, end);
    SYS_ASSERT(last >= beg);
    SYS_ASSERT(last <= end);
    if (last != end)
        return grid_result_t(
            &*last, Ext::size_cast(last - beg) + 2);
    else
        return grid_result_t();
}

#if defined(CONFIG_GRID_SORT_QUICK3WAY2) && __cplusplus < 201103L
#error the 'quick3way2' algorithm needs a C++11 compiler
#endif

class quick3way_t
{
public:
    static void sort(grid_t& a)
    { if (a.size()) sort(data(a), data(a) + a.size() - 1, 0); }

#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
    static void sort2(grid_t& a)
    { if (a.size()) sort2(data(a), data(a) + a.size() - 1, 0); }
#endif

private:
    static char at(const str_t& s, size_t k)
    { return k < s.size() ? s[k] : 0; }

    typedef unsigned char uchar_t;

    static bool ult(char a, char b)
    { return
        static_cast<uchar_t>(a) <
        static_cast<uchar_t>(b); }

    static bool ugt(char a, char b)
    { return
        static_cast<uchar_t>(a) >
        static_cast<uchar_t>(b); }

    static void sort(str_t* lo, str_t* hi, size_t d);

#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
    static void insertion(str_t* lo, str_t* hi);

    static void sort2(str_t* lo, str_t* hi, size_t d);
#endif
};

// Sedgewick & Wayne: Algorithms, 4th edition
// Algorithm 5.3, Three-way string quicksort, p. 720

void quick3way_t::sort(str_t* lo, str_t* hi, size_t d)
{
    if (hi <= lo)
        return;

    str_t *lt = lo, *gt = hi;
    char v = at(*lo, d);
    str_t *i = lo + 1;
    while (i <= gt) {
        char t = at(*i, d);
        if (ult(t, v))
            (lt ++)->swap(*i ++);
        else
        if (ugt(t, v))
            (gt --)->swap(*i);
        else
            i ++;
    }

    // *lo ... *(lt - 1) < v == *lt ... *gt < *(gt + 1) ... *hi
    sort(lo, lt - 1, d);
    if (v) sort(lt, gt, d + 1);
    sort(gt + 1, hi, d);
}

#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
void quick3way_t::insertion(str_t* lo, str_t* hi)
{
    for (str_t* i = lo + 1; i <= hi; i ++) {
        // *lo <= ... <= *(i - 1)
        str_t* j = i - 1;
        str_t v = std::move(*i);
        for (; j >= lo && v < *j; j --)
            *(j + 1) = std::move(*j);
        // C := number of iterations
        // case C == 0:
        //      lo <= j == i - 1 && *lo <= ... <= *(i - 1) <= v
        // case C == i - j - 1:
        //      lo <= j < i - 1  && *j <= v < *(j + 1) <= ... <= *(i - 1)
        // case C == i - lo:
        //      j == lo - 1 && v < *lo <= ... <= *(i - 1)
        // note: i > lo        
        *(j + 1) = std::move(v);
    } 
}

// Sedgewick: Algorithms in C, 3rd edition
// Program 10.3, Three-way radix quicksort, p. 422

void quick3way_t::sort2(str_t* lo, str_t* hi, size_t d)
{
    if (Ext::size_cast(hi - lo) <=
        globals.insertion_threshold) {
        insertion(lo, hi);
        return;
    }
    char v = at(*hi, d);
    str_t* i = lo - 1;
    str_t* j = hi;
    str_t* p = lo - 1;
    str_t* q = hi;
    while (i < j) {
        while (ult(at(*++ i, d), v))
            ;
        while (ugt(at(*-- j, d), v) && j > lo)
            ;
        if (i > j)
            break;
        i->swap(*j);
        if (at(*i, d) == v)
            i->swap(*++ p);
        if (at(*j, d) == v)
            j->swap(*-- q);
    }
    if (p == q) {
        if (v) sort2(lo, hi, d + 1);
        return;
    }
    if (ult(at(*i, d), v))
        i ++;
    for (str_t* k = lo; k <= p; k ++, j --)
        k->swap(*j);
    for (str_t* k = hi; k >= q; k --, i ++)
        k->swap(*i);
    if (lo <= j)
        sort2(lo, j, d); //!!! without 'if'
    if (i == hi && at(*i, d) == v)
        i ++;
    if (v)
        sort2(j + 1, i - 1, d + 1);
    if (i <= hi)
        sort2(i, hi, d); //!!! without 'if'
}
#endif

enum sort_algo_t {
    quick3way,
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
    quick3way2,
#endif
    stdcxx
};

inline void sort(grid_t& grid, sort_algo_t algo)
{
    if (algo == quick3way)
        quick3way_t::sort(grid);
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
    else
    if (algo == quick3way2)
        quick3way_t::sort2(grid);
#endif
    else
    if (algo == stdcxx)
        std::sort(grid.begin(), grid.end(), str_cmp_t());
    else
        SYS_UNEXPECT_ERR("unexpected algo='%d'", algo);
}

inline void unique(grid_t& grid)
{
    grid_t::iterator
        last = std::unique(grid.begin(), grid.end());
    grid.erase(last, grid.end());
}

class trie_t
{
public:
    trie_t(
#ifdef DEBUG
        bool _debug,
#endif
        bool _dots,
        const grid_t& _grid,
        std::ostream& _ost) :
#ifdef DEBUG
        debug(_debug),
#endif
        dots(_dots),
        beg(data(_grid)),
        end(data(_grid) + _grid.size()),
        max(beg < end
            ? std::max_element(
                beg, end, size_cmp_t())->size()
            : 0),
        ost(_ost)
    {}

    enum gen_type_t {
        compact,
        wide
    };

    void gen_trie(gen_type_t type) const;

private:
    struct size_cmp_t
    {
        bool operator()(const str_t& a, const str_t& b) const
        { return a.size() < b.size(); }
    };

    indent_t indent(size_t level) const
    { return indent_t(level, dots); }

    static char at(const str_t& str, size_t k)
    { return k < str.size() ? str[k] : 0; }

#ifdef DEBUG
    void print_debug(
        const str_t* beg, const str_t* end,
        size_t pos, size_t level) const;
#endif

    const str_t* split(
        const str_t* beg, const str_t* end,
        size_t& pos) const;

    void gen_compact_trie(
        const str_t* beg, const str_t* end,
        size_t pos, size_t level) const;

    void gen_wide_trie(
        const str_t* beg, const str_t* end,
        size_t pos, size_t level) const;

#ifdef DEBUG
    const unsigned     debug : 1;
#endif
    const unsigned     dots : 1;

    const str_t* const beg;
    const str_t* const end;
    const size_t       max;
    std::ostream&      ost;
};

#ifdef DEBUG
void trie_t::print_debug(
    const str_t* beg, const str_t* end,
    size_t pos, size_t level) const
{
    ost << indent(level)
        << "// " << pos << ' ';
    Ext::transform(beg, end - 1,
        std::ostream_iterator<repr_t>(ost, ","),
        static_cast<repr_t(*)(const str_t&)>(&repr));
    ost << repr(*(end - 1)) << std::endl;
}
#endif

inline const str_t* trie_t::split(
    const str_t* beg, const str_t* end,
    size_t& pos) const
{
    for (; pos < max; ++ pos) {
        const str_t* ptr = beg;
        const char c = at(*ptr, pos);
        while (++ ptr < end) {
            if (at(*ptr, pos) != c)
                return ptr;
        }
    }
    return end;
}

void trie_t::gen_compact_trie(
    const str_t* beg, const str_t* end,
    size_t pos, size_t level) const
{
    using namespace std;

    SYS_ASSERT(beg <= end);
    if (beg == end)
        return;

#ifdef DEBUG
    if (debug)
        print_debug(beg, end, pos, level);
#endif

    size_t k = pos;
    const str_t* ptr = split(beg, end, k);

    const bool t = ptr < end;
    const repr_t p = repr(*beg);
    const range_t r(pos, k);

    const bool o = r.intersect(p.range()).size();
    const bool s = beg->size() <= k;

    if (o)
        ost << indent(level);
#ifdef DEBUG
    if (o && debug)
        ost << '[' << r << ',' << p << "] ";
#endif
    if (o)
        ost << p.intersect(r);

    if (o && s)
        ost << ": " << p;

    if (o && t)
        ost << " {\n";

    if (t) {
        gen_compact_trie(beg, ptr, k, level + o);
        gen_compact_trie(ptr, end, k, level + o);
    }

    if (o && t)
        ost << indent(level) << '}';

    if (o)
        ost << '\n';
}

void trie_t::gen_wide_trie(
    const str_t* beg, const str_t* end,
    size_t pos, size_t level) const
{
    using namespace std;

    SYS_ASSERT(beg <= end);
    if (beg == end)
        return;

#ifdef DEBUG
    if (debug)
        print_debug(beg, end, pos, level);
#endif

    size_t k = pos;
    const str_t* ptr = split(beg, end, k);

    const bool t = ptr < end;
    const repr_t p = repr(*beg);
    const range_t r(pos, k);

    const size_t o = r.intersect(p.range()).size();
    const bool s = beg->size() <= k;

    for (size_t i = o, j = pos; i; -- i, ++ j) {
        const range_t r(j, j + 1);

        if (1)
            ost << indent(level + (o - i));
#ifdef DEBUG
        if (debug)
            ost << '[' << r << ',' << p << "] ";
#endif
        if (1)
            ost << p.intersect(r, true);

        if (i == 1 && s)
            ost << ": " << p;

        if (i > 1 || t)
            ost << " {\n";
    }

    if (t) {
        gen_wide_trie(beg, ptr, k, level + o);
        gen_wide_trie(ptr, end, k, level + o);
    }

    for (size_t i = o; i; -- i) {
        if (i < o || t)
            ost << indent(level + i - 1) << '}';

        if (1)
            ost << '\n';
    }
}

void trie_t::gen_trie(gen_type_t type) const
{
    using namespace std;

    SYS_ASSERT(beg <= end);
    if (beg == end)
        return;

    if (type == compact ||
        type == wide)
        ost << "{\n";

    if (type == compact)
        gen_compact_trie(beg, end, 0, 1);
    else
    if (type == wide)
        gen_wide_trie(beg, end, 0, 1);
    else
        SYS_UNEXPECT_ERR("gen-type='%d'", type);

    if (type == compact ||
        type == wide)
        ost << "}\n";
}

} // namespace Trie

class options_t :
    private global_options_t
{
public:
    struct Error : public std::runtime_error
    {
        Error(const char* what) : std::runtime_error(what) {}
    };

    typedef Sys::time_t::type_t time_type_t;

    enum action_t {
        load_only,
        sort_only,
        gen_trie
    };

    enum input_lib_t {
        cxx_iostreams,
        c_stdio,
    };

    enum input_algo_t {
        append_grid,
        grow_grid
    };

    enum sort_algo_t {
        quick3way,
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        quick3way2,
#endif
        stdcxx
    };

    enum sort_type_t {
        unsorted,
        sorted,
        unique
    };

    enum gen_type_t {
        compact,
        wide
    };

    options_t() :
        action(gen_trie),
#ifdef DEBUG
        debug(false),
#endif
        check_sorted(false),
        sync_with_stdio(false),
        quote(false),
        dots(false),

        timings(Sys::time_t::none),
        input_lib(cxx_iostreams),
        input_algo(append_grid),
        sort_algo(quick3way),
        sort_type(unsorted),
        gen_type(wide),
        grid_size(1024),
        buf_size(1024)
    {
        dump_backtrace = false;
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        insertion_threshold = 10;
#endif
        utime_prec = 3;
    }

    void getenv();
    void parse(int argc, char* const argv[]);

    action_t      action;

#ifdef DEBUG
    unsigned      debug : 1;
#endif
    unsigned      check_sorted : 1;
    unsigned      sync_with_stdio : 1;
    unsigned      quote : 1;
    unsigned      dots : 1;

    time_type_t   timings;
    input_lib_t   input_lib;
    input_algo_t  input_algo;
    sort_algo_t   sort_algo;
    sort_type_t   sort_type;
    gen_type_t    gen_type;
    size_t        grid_size;
    size_t        buf_size;

    size_t        argc;
    char *const  *argv;

private:
    static void error(const char* msg, ...) PRINTF_FMT(1);
    static void invalid_env_var(const char* opt_name, const char* opt_arg = 0);
    static void invalid_opt_arg(const char* opt_name, const char* opt_arg);
    static void missing_opt_arg(const char* opt_name);
    static void missing_opt_arg(char opt_name);
    static void invalid_opt(const char* opt_name);
    static void invalid_opt(char opt_name);

    typedef void (*invalid_opt_arg_func_t)(
        const char* opt_name, const char* opt_arg);

    template<invalid_opt_arg_func_t>
    static size_t parse_args_optarg(
        const char* opt_name, const char* opt_arg,
        const char* const* args, size_t n_args);

    template<invalid_opt_arg_func_t>
    static size_t parse_size_num_optarg(
        const char* opt_name, const char* opt_arg);

    template<invalid_opt_arg_func_t>
    static time_type_t parse_timings_optarg(
        const char* opt_name, const char* opt_arg);

    template<invalid_opt_arg_func_t>
    static size_t parse_utime_prec_optarg(
        const char* opt_name, const char* opt_arg);

    void print(std::ostream&) const;

    static void dumpenv();
    static void version();
    static void usage();

    void dump() const;

    class base_var_t;
    class env_var_t;

    class base_var_ref_t;
    class env_var_ref_t;

    static const char* noyes[];
    static const char* actions[];
    static const char* input_libs[];
    static const char* input_algos[];
    static const char* sort_algos[];
    static const char* sort_types[];
    static const char* gen_types[];
};

const char* options_t::noyes[] = {
    "no",
    "yes"
};
const char* options_t::actions[] = {
    "load-only",    // action_t::load_only
    "sort-only",    // action_t::sort_only
    "gen-trie"      // action_t::gen_trie
};
const char* options_t::input_libs[] = {
    "c++iostreams", // input_lib_t::cxx_iostreams
    "cstdio",       // input_lib_t::c_stdio
};
const char* options_t::input_algos[] = {
    "append-grid",  // input_algo_t::append_grid
    "grow-grid",    // input_algo_t::grow_grid
};
const char* options_t::sort_algos[] = {
    "quick3way",    // sort_algo_t::quick3way
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
    "quick3way2",   // sort_algo_t::quick3way2
#endif
    "stdc++",       // sort_type_t::stdcxx
};
const char* options_t::sort_types[] = {
    "unsorted",     // sort_type_t::unsorted
    "sorted",       // sort_type_t::sorted
    "unique"        // sort_type_t::unique
};
const char* options_t::gen_types[] = {
    "compact",      // gen_type_t::compact
    "wide",         // gen_type_t::wide
};

class options_t::base_var_t
{
public:
    enum types_t {
        boolean,
        size,
        action,
        timings,
        utime_prec,
        input_lib,
        input_algo,
        sort_algo,
        sort_type,
        gen_type
    };

    template<types_t>
    struct var_t;

protected:
    base_var_t(const char* _name, types_t _type, bool _quiet) :
        name(_name), type(_type), quiet(_quiet)
    {}

    base_var_t(const std::string& _name, types_t _type, bool _quiet) :
        name(_name), type(_type), quiet(_quiet)
    {}

#if __cplusplus >= 201103L
    base_var_t(std::string&& _name, types_t _type, bool _quiet) :
        name(std::move(_name)), type(_type), quiet(_quiet)
    {}
#endif

    void error(const Error&) const;

    void print_name(std::ostream&) const;

    template<types_t t>
    typename var_t<t>::type_t parse(const char*) const;

    template<types_t t>
    static void print(std::ostream& ost, typename var_t<t>::type_t val)
    { ost << val; }

    void check_type(types_t t) const
    { SYS_ASSERT(type == t); }

    const char* getenv() const
    { return ::getenv(name.c_str()); }

    std::string const name;
    types_t     const type;

private:
    size_t parse_size(const char* val) const
    { check_type(size); return parse_size_num_optarg<invalid_env_var>(
        name.c_str(), val); }

    time_type_t parse_timings(const char* val) const
    { check_type(timings); return parse_timings_optarg<invalid_env_var>(
        name.c_str(), val); }

    size_t parse_utime_prec(const char* val) const
    { check_type(utime_prec); return parse_utime_prec_optarg<invalid_env_var>(
        name.c_str(), val); }

    template<types_t t>
    typename var_t<t>::type_t parse_args(
        const char* val, const char* const* args, size_t n_args) const
    {
        check_type(t); 
        typedef typename var_t<t>::type_t type_t;
        return static_cast<type_t>(parse_args_optarg<invalid_env_var>(
            name.c_str(), val, args, n_args));
    }

    const bool quiet;
};

class options_t::env_var_t :
    protected options_t::base_var_t
{
public:
    typedef options_t::base_var_t base_t;

    template<types_t t>
    static env_var_t make(const char* n)
    { return env_var_t(n, t, true); }

    void print(std::ostream& ost) const;

    static std::string env_name(const char* s)
    { std::string r(s); env_name0(r); return r; }

    static std::string env_name(const std::string& s)
    { std::string r(s); env_name0(r); return r; }

#if __cplusplus >= 201103L
    static std::string env_name(std::string&& s)
    { env_name0(s); return std::move(s); }
#endif

private:
    env_var_t(const char* _name, types_t _type, bool _quiet) :
        base_t(env_name(_name), _type, _quiet)
    {}

    static void env_name0(std::string& s)
    {
        s.insert(0, "trie_");
        Ext::transform(s.begin(), s.end(), s.begin(), toupper);
    }

    template<types_t t>
    void print(std::ostream& ost) const
    {
        if (const char* s = base_t::getenv())
            base_t::print<t>(ost, parse<t>(s));
    }
};

class options_t::base_var_ref_t :
    protected options_t::base_var_t
{
public:
    typedef options_t::base_var_t base_t;

    template<types_t t>
    static base_var_ref_t
        make(const char* n, const typename base_t::var_t<t>::type_t& v)
    { return base_var_ref_t(n, t, &v, true, '-'); }

    void print(std::ostream&) const;

protected:
    base_var_ref_t(
        const char* _name, types_t _type, const void* _ptr,
        bool _quiet, char _sep) :
        base_t(norm(_name, _sep), _type, _quiet),
        ptr(_ptr)
    {}

    base_var_ref_t(
        const std::string& _name, types_t _type, const void* _ptr,
        bool _quiet, char _sep) :
        base_t(norm(_name, _sep), _type, _quiet),
        ptr(_ptr)
    {}

#if __cplusplus >= 201103L
    base_var_ref_t(
        std::string&& _name, types_t _type, const void* _ptr,
        bool _quiet, char _sep) :
        base_t(norm(std::move(_name), _sep), _type, _quiet),
        ptr(_ptr)
    {}
#endif

    template<types_t t>
    const typename base_t::var_t<t>::type_t& ref() const
    {
        typedef 
            typename base_t::var_t<t>::type_t
            type_t;
        check_type(t);
        return *static_cast<const type_t*>(ptr);
    }

    template<types_t t>
    void print(std::ostream& ost) const
    { base_t::print<t>(ost, ref<t>()); }

private:
    static void norm0(std::string& s, char c)
    { if (c != '_') std::replace(s.begin(), s.end(), '_', c); }

    static std::string norm(const char* s, char c)
    { std::string r(s); norm0(r, c); return r; }

    static std::string norm(const std::string& s, char c)
    { std::string r(s); norm0(r, c); return r; }

#if __cplusplus >= 201103L
    static std::string norm(std::string&& s, char c)
    { norm0(s, c); return std::move(s); }
#endif

    const void* ptr;
};

class options_t::env_var_ref_t :
    protected options_t::base_var_ref_t
{
public:
    typedef options_t::base_var_ref_t base_t;

    template<types_t t>
    static env_var_ref_t
        make(const char* n, typename base_t::var_t<t>::type_t& v)
    { return env_var_ref_t(n, t, &v, false, '_'); }

    void getenv() const;

    void print(std::ostream& ost) const
    { base_t::print(ost); }

private:
    env_var_ref_t(
        const char* _name, types_t _type, void* _ptr,
        bool _quiet, char _sep) :
        base_t(env_var_t::env_name(_name), _type, _ptr,
            _quiet, _sep)
    {}

    template<types_t t>
    typename base_t::var_t<t>::type_t& ref() const
    {
        typedef 
            typename base_t::var_t<t>::type_t
            type_t;
        return const_cast<type_t&>(base_t::ref<t>());
    }

    template<types_t t>
    void assign() const
    {
        if (const char* s = base_t::getenv())
            ref<t>() = parse<t>(s);
    }
};

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::boolean>
{ typedef bool type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::size>
{ typedef size_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::action>
{ typedef action_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::timings>
{ typedef time_type_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::utime_prec>
{ typedef size_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::input_lib>
{ typedef input_lib_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::input_algo>
{ typedef input_algo_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::sort_algo>
{ typedef sort_algo_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::sort_type>
{ typedef sort_type_t type_t; };

template<> struct
    options_t::base_var_t::var_t<
    options_t::base_var_t::gen_type>
{ typedef gen_type_t type_t; };

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::boolean>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::boolean>(const char* val) const
{ return parse_args<boolean>(
    val, noyes, Sys::array_size(noyes)); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::action>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::action>(const char* val) const
{ return parse_args<action>(
    val, actions, Sys::array_size(actions)); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::size>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::size>(const char* val) const
{ return parse_size(val); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::timings>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::timings>(const char* val) const
{ return parse_timings(val); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::utime_prec>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::utime_prec>(const char* val) const
{ return parse_utime_prec(val); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::input_lib>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::input_lib>(const char* val) const
{ return parse_args<input_lib>(
    val, input_libs, Sys::array_size(input_libs)); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::input_algo>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::input_algo>(const char* val) const
{ return parse_args<input_algo>(
    val, input_algos, Sys::array_size(input_algos)); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::sort_algo>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::sort_algo>(const char* val) const
{ return parse_args<sort_algo>(
    val, sort_algos, Sys::array_size(sort_algos)); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::sort_type>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::sort_type>(const char* val) const
{ return parse_args<sort_type>(
    val, sort_types, Sys::array_size(sort_types)); }

template<> inline
options_t::base_var_t::var_t<
options_t::base_var_t::gen_type>::type_t
options_t::base_var_t::parse<
options_t::base_var_t::gen_type>(const char* val) const
{ return parse_args<gen_type>(
    val, gen_types, Sys::array_size(gen_types)); }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::boolean>(std::ostream& ost, bool val)
{ ost << Sys::array(noyes)[val]; }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::action>(std::ostream& ost, action_t val)
{ ost << Sys::array(actions)[val]; }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::timings>(std::ostream& ost, time_type_t val)
{ ost << Sys::time_t::name(val); }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::input_lib>(std::ostream& ost, input_lib_t val)
{ ost << Sys::array(input_libs)[val]; }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::input_algo>(std::ostream& ost, input_algo_t val)
{ ost << Sys::array(input_algos)[val]; }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::sort_algo>(std::ostream& ost, sort_algo_t val)
{ ost << Sys::array(sort_algos)[val]; }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::sort_type>(std::ostream& ost, sort_type_t val)
{ ost << Sys::array(sort_types)[val]; }

template<>
inline void 
options_t::base_var_t::print<
options_t::base_var_t::gen_type>(std::ostream& ost, gen_type_t val)
{ ost << Sys::array(gen_types)[val]; }

void options_t::base_var_t::error(const Error& e) const
{
    if (!quiet)
        std::cerr
            << program
            << ": warning: "
            << e.what()
            << std::endl;
}

void options_t::base_var_t::print_name(std::ostream& ost) const
{
    using namespace std;

    const size_t w = 17;
    size_t n = name.size();
    ost << left
        << name
        << setw(w > n ? w - n : 0)
        << ":";
    if (n >= w) ost
        << '\n'
        << setw(w)
        << ' ';
}

void options_t::env_var_t::print(std::ostream& ost) const
{
    static void (env_var_t::* const prints[])(std::ostream&) const = {
        &env_var_t::print<boolean>,    // types_t::boolean
        &env_var_t::print<size>,       // types_t::size
        &env_var_t::print<action>,     // types_t::action
        &env_var_t::print<timings>,    // types_t::timings
        &env_var_t::print<utime_prec>, // types_t::utime_prec
        &env_var_t::print<input_lib>,  // types_t::input_lib
        &env_var_t::print<input_algo>, // types_t::input_algo
        &env_var_t::print<sort_algo>,  // types_t::sort_algo
        &env_var_t::print<sort_type>,  // types_t::sort_type
        &env_var_t::print<gen_type>,   // types_t::gen_type
    };

    bool r = false;
    std::ostringstream o;

    try {
        (this->*Sys::array(prints)[type])(o);
    }
    catch (const Error& e) {
        error(e);
        r = true;
    }

    print_name(ost);
    std::string s = !r ? o.str() : "?";
    ost << (s.size() ? s.c_str() : "-")
        << std::endl;
}

void options_t::dumpenv()
{
#define ENV_VAR_TYPE_D(NAME, TYPE) \
        env_var_t::make<base_var_t::TYPE>(#NAME)
#define ENV_VAR_NAME_D(NAME) \
        ENV_VAR_TYPE_D(NAME, NAME)

    const env_var_t vars[] = {
        ENV_VAR_NAME_D(input_lib),
        ENV_VAR_NAME_D(input_algo),
        ENV_VAR_NAME_D(sort_algo),
        ENV_VAR_NAME_D(timings),
        ENV_VAR_NAME_D(utime_prec),
        ENV_VAR_TYPE_D(grid_size, size),
        ENV_VAR_TYPE_D(buf_size, size),
        ENV_VAR_TYPE_D(dump_backtrace, boolean),
    };

#undef ENV_VAR_TYPE_D
#undef ENV_VAR_NAME_D

    for (const env_var_t*
        p = vars;
        p < vars + Sys::array_size(vars);
        p ++)
        p->print(std::cout);
}

void options_t::base_var_ref_t::print(std::ostream& ost) const
{
    typedef options_t::base_var_t base_var_t;

    static void (base_var_ref_t::* const prints[])(std::ostream&) const = {
        &base_var_ref_t::print<base_var_t::boolean>,    // types_t::boolean
        &base_var_ref_t::print<base_var_t::size>,       // types_t::size
        &base_var_ref_t::print<base_var_t::action>,     // types_t::action
        &base_var_ref_t::print<base_var_t::timings>,    // types_t::timings
        &base_var_ref_t::print<base_var_t::utime_prec>, // types_t::utime_prec
        &base_var_ref_t::print<base_var_t::input_lib>,  // types_t::input_lib
        &base_var_ref_t::print<base_var_t::input_algo>, // types_t::input_algo
        &base_var_ref_t::print<base_var_t::sort_algo>,  // types_t::sort_algo
        &base_var_ref_t::print<base_var_t::sort_type>,  // types_t::sort_type
        &base_var_ref_t::print<base_var_t::gen_type>,   // types_t::gen_type
    };

    print_name(ost);
    (this->*Sys::array(prints)[type])(ost);
    ost << '\n';
}

void options_t::print(std::ostream& ost) const
{
#define REF_VAR_INIT_F(NAME) \
        const bool NAME = this->NAME
#define REF_VAR_TYPE_N(NAME, TYPE, LABEL) \
        base_var_ref_t::make<base_var_t::TYPE>(LABEL, NAME)
#define REF_VAR_TYPE_D(NAME, TYPE) \
        REF_VAR_TYPE_N(NAME, TYPE, #NAME)
#define REF_VAR_NAME_N(NAME, LABEL) \
        REF_VAR_TYPE_N(NAME, NAME, LABEL)
#define REF_VAR_NAME_D(NAME) \
        REF_VAR_TYPE_D(NAME, NAME)

    REF_VAR_INIT_F(sync_with_stdio);
    REF_VAR_INIT_F(check_sorted);
    REF_VAR_INIT_F(quote);
    REF_VAR_INIT_F(dots);
#ifdef DEBUG
    REF_VAR_INIT_F(debug);
#endif
    REF_VAR_INIT_F(dump_backtrace);

    const base_var_ref_t vars[] = {
        REF_VAR_NAME_D(action),
        REF_VAR_NAME_D(input_lib),
        REF_VAR_NAME_N(input_algo, "input-algorithm"),
        REF_VAR_TYPE_D(sync_with_stdio, boolean),
        REF_VAR_NAME_N(sort_algo, "sort-algorithm"),
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        REF_VAR_TYPE_N(insertion_threshold, size, "ins-threshold"),
#endif
        REF_VAR_TYPE_D(check_sorted, boolean),
        REF_VAR_NAME_D(sort_type),
        REF_VAR_NAME_D(gen_type),
        REF_VAR_TYPE_D(grid_size, size),
        REF_VAR_TYPE_D(buf_size, size),
        REF_VAR_NAME_D(timings),
        REF_VAR_NAME_D(utime_prec),
        REF_VAR_TYPE_D(quote, boolean),
        REF_VAR_TYPE_D(dots, boolean),
#ifdef DEBUG
        REF_VAR_TYPE_D(debug, boolean),
#endif
        REF_VAR_TYPE_D(dump_backtrace, boolean),
        REF_VAR_TYPE_D(argc, size),
    };

#undef REF_VAR_INIT_F
#undef REF_VAR_TYPE_N
#undef REF_VAR_TYPE_D
#undef REF_VAR_NAME_N
#undef REF_VAR_NAME_D

    for (const base_var_ref_t*
        p = vars;
        p < vars + Sys::array_size(vars);
        p ++)
        p->print(ost);

    using namespace std;

    const size_t w = 12;

    for (char
        *const *p = argv,
        *const *e = argv + argc;
        p < e; p ++) {
        size_t i =
            Ext::size_cast(p - argv);
        size_t n =
            Ext::digits10(i);
        ost << "argv["
            << i
            << left
            << setw(w > n ? w - n : 0)
            << "]:";
        if (*p) ost
            << Sys::repr(*p, false);
        ost << '\n';
    }
}

void options_t::env_var_ref_t::getenv() const
{
    typedef options_t::base_var_t base_var_t;

    static void (env_var_ref_t::* const assigns[])() const = {
        &env_var_ref_t::assign<base_var_t::boolean>,    // types_t::boolean
        &env_var_ref_t::assign<base_var_t::size>,       // types_t::size
        &env_var_ref_t::assign<base_var_t::action>,     // types_t::action
        &env_var_ref_t::assign<base_var_t::timings>,    // types_t::timings
        &env_var_ref_t::assign<base_var_t::utime_prec>, // types_t::utime_prec
        &env_var_ref_t::assign<base_var_t::input_lib>,  // types_t::input_lib
        &env_var_ref_t::assign<base_var_t::input_algo>, // types_t::input_algo
        &env_var_ref_t::assign<base_var_t::sort_algo>,  // types_t::sort_algo
        &env_var_ref_t::assign<base_var_t::sort_type>,  // types_t::sort_type
        &env_var_ref_t::assign<base_var_t::gen_type>,   // types_t::gen_type
    };

    try {
        (this->*Sys::array(assigns)[type])();
    }
    catch (const Error& e) {
        error(e);
    }
}

void options_t::getenv()
{
#define ENV_VAR_INIT_F(NAME) \
        bool NAME = this->NAME
#define ENV_VAR_DONE_F(NAME) \
        this->NAME = NAME
#define ENV_VAR_TYPE_D(NAME, TYPE) \
        env_var_ref_t::make<base_var_t::TYPE>(#NAME, NAME)
#define ENV_VAR_NAME_D(NAME) \
        ENV_VAR_TYPE_D(NAME, NAME)

    ENV_VAR_INIT_F(dump_backtrace);

    const env_var_ref_t vars[] = {
        ENV_VAR_NAME_D(input_lib),
        ENV_VAR_NAME_D(input_algo),
        ENV_VAR_NAME_D(sort_algo),
        ENV_VAR_NAME_D(timings),
        ENV_VAR_NAME_D(utime_prec),
        ENV_VAR_TYPE_D(grid_size, size),
        ENV_VAR_TYPE_D(buf_size, size),
        ENV_VAR_TYPE_D(dump_backtrace, boolean),
    };

    for (const env_var_ref_t*
        p = vars;
        p < vars + Sys::array_size(vars);
        p ++)
        p->getenv();

    ENV_VAR_DONE_F(dump_backtrace);

#undef ENV_VAR_INIT_F
#undef ENV_VAR_DONE_F
#undef ENV_VAR_TYPE_D
#undef ENV_VAR_NAME_D
}

void options_t::version()
{
    std::cout << program << ": version " << verdate << "\n\n" << license;
}

void options_t::usage()
{
    std::cout <<
        "usage: " << program << " [OPTION]...\n"
        "where the options are:\n"
        "  -L|--load-only       action: only load input and nothing more\n"
        "  -S|--sort-only       action: only sort input and print it out\n"
        "  -G|--gen-trie        action: fully process input (load it, sort it if needed),\n"
        "                         than generate the respective trie on stdout (default)\n"
        "  -C|--check-sorted    check that input is sorted properly\n" 
        "  -s|--sorted          sort type: input is sorted but not uniquely\n"
        "  -u|--unique          sort type: input is sorted uniquely\n"
        "  -U|--unsorted        sort type: input is not sorted (default)\n"
        "     --quick3way       sort algorithm: use 'quick3way' algorithm (default)\n"
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        "     --quick3way2      sort algorithm: use 'quick3way2' algorithm\n"
#endif
        "     --stdc++          sort algorithm: use standard C++ library algorithm\n"
        "  -c|--compact         gen type: generate compact trie on output\n"
        "  -w|--wide            gen type: generate wide trie on output (default)\n"
        "     --c++iostreams    input lib: use C++ I/O Streams library (default)\n"
        "     --cstdio          input lib: use C Standard I/O library\n"
        "     --append-grid     input algorithm: use 'append-grid' algorithm (default)\n"
        "     --grow-grid       input algorithm: use 'grow-grid' algorithm\n"
        "     --[no-]sync-with-stdio\n"
        "                       keep C++ streams in sync with C streams or\n"
        "                         otherwise do not (default not)\n"
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        "     --ins[ertion]-threshold NUM\n"
        "                       the threshold point (a small value) upon which\n"
        "                         the 'quick2way2' algorithm cuts off recursion by\n"
        "                         falling down to sorting by 'insertion' (default: 10)\n"
#endif
        "     --timings [TYPE]  print out on stderr the real/user/system times\n"
        "                         spent executing the main inner tasks if TYPE is\n"
        "                         given as 'ctime' or 'utime'; 'ctime' uses 'times'\n"
        "                         library call, while 'utime' uses 'gettimeofday'\n"
        "                         and 'getrusage' tandem; when the arg is ommited\n"
        "                         use 'ctime'; default is 'none'\n"
        "     --utime-prec NUM  when printing out timings of type 'utime', compute\n"
        "                         the fractional parts with the given precision\n"
        "                         (valid range: 3-6; default: 3)\n"
        "     --grid-size NUM   use NUM as the initial size of the grid buffer;\n"
        "                         the default value is 1024\n"
        "     --buf-size NUM    use NUM as the initial size of the input buffer;\n"
        "                         the default value is 1024\n"
        "  -q|--[no-]quote      when sorting only, print out quoted strings\n"
        "                         or otherwise do not (default not)\n"
        "  -D|--[no-]dots       put indenting dots in structure print outs\n"
        "                         or otherwise do not (default not)\n"
#ifdef DEBUG
        "  -d|--[no-]debug      print some debugging output or otherwise\n"
        "                         do not print debugging output at all (default)\n"
#endif
        "     --[no-]dump-backtrace\n"
        "                       print out the backtrace of the program on fatal error\n"
        "                         or otherwise do not (default not)\n"
        "     --dump-env[iron]-vars\n"
        "                       print environment vars and exit\n"
        "     --dump-options    print options and exit\n"
        "  -v|--version         print version numbers and exit\n"
        "  -?|--help            display this help info and exit\n";
}

inline void options_t::dump() const
{ print(std::cout); }

void options_t::error(const char* msg, ...)
{
    va_list args;
    va_start(args, msg);
    Sys::verror<Error>(msg, args);
    //!!!VA_END va_end(args);
}

void options_t::invalid_env_var(const char* opt_name, const char* opt_arg)
{
    std::ostringstream ost;

    if (opt_arg)
        ost << Sys::repr(opt_arg);
    else
        ost << "null";

    error("ignoring invalid value of environment var %s: %s",
        opt_name, ost.str().c_str());
}

void options_t::invalid_opt_arg(const char* opt_name, const char* opt_arg)
{
    error("invalid argument for '%s' option: '%s'", opt_name, opt_arg);
}

void options_t::missing_opt_arg(const char* opt_name)
{
    error("argument for option '%s' not found", opt_name);
}

void options_t::missing_opt_arg(char opt_name)
{
    error("argument for option '-%c' not found", opt_name);
}

void options_t::invalid_opt(const char* opt_name)
{
    error("invalid command line option '%s'", opt_name);
}

void options_t::invalid_opt(char opt_name)
{
    error("invalid command line option '-%c'", opt_name);
}

template<
    options_t::invalid_opt_arg_func_t
    invalid_opt_arg_func>
size_t options_t::parse_args_optarg(
    const char* opt_name, const char* opt_arg,
    const char* const* args, size_t n_args)
{
    using namespace std;
    char const* const *b = args;
    char const* const *e = b + n_args;
    char const* const *p = find_if(
        b, e, bind1st(Sys::eq(), opt_arg));
    if (p >= e)
        invalid_opt_arg_func(opt_name, opt_arg);
    return Ext::size_cast(p - b);
}

template<
    options_t::invalid_opt_arg_func_t
    invalid_opt_arg_func>
size_t options_t::parse_size_num_optarg(
    const char* opt_name, const char* opt_arg)
{
    size_t r = 0;
    if (!Ext::parse_size_num(opt_arg, r))
        invalid_opt_arg_func(opt_name, opt_arg);
    return r;
}

template<
    options_t::invalid_opt_arg_func_t
    invalid_opt_arg_func>
options_t::time_type_t options_t::parse_timings_optarg(
    const char* opt_name, const char* opt_arg)
{
    time_type_t r = Sys::time_t::ctime;
    if (opt_arg && !Sys::time_t::type(opt_arg, r))
        invalid_opt_arg_func(opt_name, opt_arg);
    return r;
}

template<
    options_t::invalid_opt_arg_func_t
    invalid_opt_arg_func>
size_t options_t::parse_utime_prec_optarg(
    const char* opt_name, const char* opt_arg)
{
    size_t r = 0;
    if (!Ext::parse_size_num(opt_arg, r) || r < 3 || r > 6)
        invalid_opt_arg_func(opt_name, opt_arg);
    return r;
}

void options_t::parse(int argc, char* const argv[])
{
    typedef int opt_t;
    struct opt_type_t
    {
        enum {
            dots         = 'D',
            quote        = 'q',
            debug        = 'd',
            load_only    = 'L',
            sort_only    = 'S',
            gen_trie     = 'G',
            check_sorted = 'C',
            sorted       = 's',
            unique       = 'u',
            unsorted     = 'U',
            compact      = 'c',
            wide         = 'w',
            version      = 'v',
            help         = '?',
            no_quote     = 256,
            quick3way,
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
            quick3way2,
#endif
            stdcxx,
            append_grid,
            grow_grid,
            cxx_iostreams,
            c_stdio,
            timings,
            utime_prec,
            grid_size,
            buf_size,
            sync_with_stdio,
            no_sync_with_stdio,
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
            insertion_threshold,
#endif
            no_dots,
            no_debug,
            dump_backtrace,
            no_dump_backtrace,
            dump_env_vars,
            dump_options
        };
    };

    static const char shorts[] = 
        ":"
#ifdef DEBUG
        "d"
#endif
        "cCDGLqsSuUvw";
    static const struct option longs[] = {
        { "load-only",           0,       0, opt_type_t::load_only },
        { "sort-only",           0,       0, opt_type_t::sort_only },
        { "gen-trie",            0,       0, opt_type_t::gen_trie },
        { "check-sorted",        0,       0, opt_type_t::check_sorted },
        { "sorted",              0,       0, opt_type_t::sorted },
        { "unique",              0,       0, opt_type_t::unique },
        { "unsorted",            0,       0, opt_type_t::unsorted },
        { "quick3way",           0,       0, opt_type_t::quick3way },
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        { "quick3way2",          0,       0, opt_type_t::quick3way2 },
#endif
        { "stdc++",              0,       0, opt_type_t::stdcxx },
        { "compact",             0,       0, opt_type_t::compact },
        { "wide",                0,       0, opt_type_t::wide },
        { "c++iostreams",        0,       0, opt_type_t::cxx_iostreams },
        { "cstdio",              0,       0, opt_type_t::c_stdio },
        { "append-grid",         0,       0, opt_type_t::append_grid },
        { "grow-grid",           0,       0, opt_type_t::grow_grid },
        { "timings",             2,       0, opt_type_t::timings },
        { "utime-prec",          1,       0, opt_type_t::utime_prec },
        { "grid-size",           1,       0, opt_type_t::grid_size },
        { "buf-size",            1,       0, opt_type_t::buf_size },
        { "sync-with-stdio",     0,       0, opt_type_t::sync_with_stdio },
        { "no-sync-with-stdio",  0,       0, opt_type_t::no_sync_with_stdio },
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        { "insertion-threshold", 1,       0, opt_type_t::insertion_threshold },
        { "ins-threshold",       1,       0, opt_type_t::insertion_threshold },
#endif
        { "quote",               0,       0, opt_type_t::quote },
        { "no-quote",            0,       0, opt_type_t::no_quote },
        { "dots",                0,       0, opt_type_t::dots },
        { "no-dots",             0,       0, opt_type_t::no_dots },
#ifdef DEBUG
        { "debug",               0,       0, opt_type_t::debug },
        { "no-debug",            0,       0, opt_type_t::no_debug },
#endif
        { "dump-env",            0,       0, opt_type_t::dump_env_vars },
        { "dump-env-vars",       0,       0, opt_type_t::dump_env_vars },
        { "dump-environ-vars",   0,       0, opt_type_t::dump_env_vars },
        { "dump-options",        0,       0, opt_type_t::dump_options },
        { "dump-backtrace",      0,       0, opt_type_t::dump_backtrace },
        { "no-dump-backtrace",   0,       0, opt_type_t::no_dump_backtrace },
        { "version",             0,       0, opt_type_t::version },
        { "help",                0, &optopt, opt_type_t::help },
        { 0,                     0,       0, 0 },
    };

    struct bits_opts {
        unsigned env: 1;
        unsigned dump: 1;
        unsigned usage: 1;
        unsigned version: 1;
    };
    struct bits_opts bits = {
        env:     0,
        dump:    0,
        usage:   0,
        version: 0,
    };

    opt_t opt;
    opterr = 0;
    optind = 1;
    while ((opt = getopt_long(
        argc, argv, &shorts[0], &longs[0], 0)) != EOF) {
        switch (opt) {
        case opt_type_t::load_only:
            action = load_only;
            break;
        case opt_type_t::sort_only:
            action = sort_only;
            break;
        case opt_type_t::gen_trie:
            action = gen_trie;
            break;
#ifdef DEBUG
        case opt_type_t::debug:
            debug = true;
            break;
        case opt_type_t::no_debug:
            debug = false;
            break;
#endif
        case opt_type_t::dots:
            dots = true;
            break;
        case opt_type_t::no_dots:
            dots = false;
            break;
        case opt_type_t::quote:
            quote = true;
            break;
        case opt_type_t::no_quote:
            quote = false;
            break;
        case opt_type_t::check_sorted:
            check_sorted = true;
            break;
        case opt_type_t::sorted:
            sort_type = sorted;
            break;
        case opt_type_t::unsorted:
            sort_type = unsorted;
            break;
        case opt_type_t::unique:
            sort_type = unique;
            break;
        case opt_type_t::quick3way:
            sort_algo = quick3way;
            break;
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        case opt_type_t::quick3way2:
            sort_algo = quick3way2;
            break;
#endif
        case opt_type_t::stdcxx:
            sort_algo = stdcxx;
            break;
        case opt_type_t::cxx_iostreams:
            input_lib = cxx_iostreams;
            break;
        case opt_type_t::c_stdio:
            input_lib = c_stdio;
            break;
        case opt_type_t::append_grid:
            input_algo = append_grid;
            break;
        case opt_type_t::grow_grid:
            input_algo = grow_grid;
            break;
        case opt_type_t::compact:
            gen_type = compact;
            break;
        case opt_type_t::wide:
            gen_type = wide;
            break;
        case opt_type_t::timings:
            timings =
                parse_timings_optarg<invalid_opt_arg>(
                    "timings", optarg);
            break;
        case opt_type_t::utime_prec:
            utime_prec =
                parse_utime_prec_optarg<invalid_opt_arg>(
                    "utime-prec", optarg);
            break;
        case opt_type_t::grid_size:
            grid_size =
                parse_size_num_optarg<invalid_opt_arg>(
                    "grid-size", optarg);
            break;
        case opt_type_t::buf_size:
            buf_size =
                parse_size_num_optarg<invalid_opt_arg>(
                    "buf-size", optarg);
            break;
        case opt_type_t::sync_with_stdio:
            sync_with_stdio = true;
            break;
        case opt_type_t::no_sync_with_stdio:
            sync_with_stdio = false;
            break;
        case opt_type_t::dump_backtrace:
            dump_backtrace = true;
            break;
        case opt_type_t::no_dump_backtrace:
            dump_backtrace = false;
            break;
        case opt_type_t::dump_env_vars:
            bits.env = 1;
            break;
        case opt_type_t::dump_options:
            bits.dump = 1;
            break;
#if defined(CONFIG_GRID_SORT_QUICK3WAY2)
        case opt_type_t::insertion_threshold:
            insertion_threshold =
                parse_size_num_optarg<invalid_opt_arg>(
                    "insertion-threshold", optarg);
            break;
#endif
        case opt_type_t::version:
            bits.version = 1;
            break;
        case ':': {
            const char* opt = argv[optind - 1];
            if (opt[0] == '-' && opt[1] == '-')
                missing_opt_arg(opt);
            else
                missing_opt_arg(optopt);
            break;
        }
        case 0:
            bits.usage = 1;
            break;
        case '?':
        default:
            if (optopt == 0)
                invalid_opt(argv[optind - 1]);
            else
            if (optopt != '?')
                invalid_opt(optopt);
            else
                bits.usage = 1;
            break;
        }
    }

    argv += optind;
    argc -= optind;

    this->argc = argc;
    this->argv = argv;

    if (bits.version)
        version();
    if (bits.env)
        dumpenv();
    if (bits.dump)
        dump();
    if (bits.usage)
        usage();

    if (bits.env ||
        bits.dump ||
        bits.version ||
        bits.usage)
        exit(0);

    globals = *this;
}

inline const options_t options(int argc, char* const argv[])
{
    options_t opt;
    opt.getenv();
    opt.parse(argc, argv);
    return opt;
}

int main(int argc, char* const argv[])
try
{
    const options_t opt = options(argc, argv);

    using namespace Ext;
    using namespace Sys;
    using namespace Trie;

    grid_t grid;
    SYS_ASSERT(
        opt.grid_size <=
        grid.max_size());
    grid.reserve(opt.grid_size);

    std::iostream::sync_with_stdio(
        opt.sync_with_stdio);

    Sys::time_t time(opt.timings);

    if (size_t r = input(
            grid,
            static_cast<input_algo_t>(opt.input_algo),
            static_cast<input_lib_t>(opt.input_lib),
            opt.buf_size))
        line_error(r, "tries cannot contain empty keys");

    if (time)
        std::cerr << time("reading");

    if (opt.action == options_t::load_only)
        return 0;

    if (opt.check_sorted &&
        opt.action != options_t::sort_only &&
        opt.sort_type != options_t::unsorted) {
        if (grid_result_t r = is_not_sorted(grid))
            line_error(
                r.line, "input is not sorted: %s",
                repr_str(*r.str).c_str());
    }
    if (opt.check_sorted &&
        opt.action != options_t::sort_only &&
        opt.sort_type == options_t::unique) {
        if (grid_result_t r = is_not_unique(grid))
            line_error(
                r.line, "input is not uniquely sorted: %s",
                repr_str(*r.str).c_str());
    }

    if (opt.action == options_t::sort_only ||
        opt.sort_type == options_t::unsorted) {
        Sys::time_t time(opt.timings);

        sort(grid, static_cast<sort_algo_t>(opt.sort_algo));

        if (time)
            std::cerr << '\n' << time("sorting");
    }

    if ((opt.action != options_t::sort_only &&
         opt.sort_type != options_t::unique) ||
        (opt.action == options_t::sort_only &&
         opt.sort_type == options_t::unique)) {
        Sys::time_t time(opt.timings);

        unique(grid);

        if (time)
            std::cerr << '\n' << time("uniqueing");
    }

    if (opt.action == options_t::sort_only) {
        Sys::time_t time(opt.timings);

        std::cout << print(grid, opt.quote);

        if (time)
            std::cerr << '\n' << time("writing");
    }
    else
    if (opt.action == options_t::gen_trie) {
        Sys::time_t time(opt.timings);

        trie_t trie(
#ifdef DEBUG
            opt.debug,
#endif
            opt.dots, grid, std::cout);

        trie.gen_trie(
            static_cast<trie_t::gen_type_t>(opt.gen_type));

        if (time)
            std::cerr << '\n' << time("generating");
    }
    else
        SYS_UNEXPECT_ERR("action='%d'", opt.action);

    if (time)
        std::cerr << '\n' << time("total");

    return 0;
}
catch (const std::exception& err) {
    std::cerr
        << program << ": error: " << err.what() << std::endl;
    return 1;
}


