#ifndef HAVE_NARRAY_H
#define HAVE_NARRAY_H

#include <array>
#include <memory>

#include <assert.h>

#include "meta.h"

enum narray_memorder_t {
  NARRAY_MEMORDER_ROW = 0,   // <- row-major, i.e., last index is fast-running (C/C++)
  NARRAY_MEMORDER_COL = 1,   // <- column-major, i.e., first index is fast-running (Fortran)
};

#ifndef NDEBUG
#define ENABLE_DEBUG 1
#endif // NDEBUG

template<
  int N,
  typename SizeT = size_t,
  narray_memorder_t MEMORDER = NARRAY_MEMORDER_ROW,
  int... Bases>
class NArrayShape {

public:
  using size_type = SizeT;
  constexpr static const int ndims = N;

  using self_t = NArrayShape<N, SizeT, MEMORDER, Bases...>;

  // default ctor
  NArrayShape()
  { }

  template<typename... Args>
  constexpr
  NArrayShape(Args... args) : _size(total_size(args...))
  {
    static_assert(sizeof...(args) == N, "Underspecified extents detected!");
    std::array<size_type, N> extents{ static_cast<size_type>(args)... };

    if (NARRAY_MEMORDER_ROW == MEMORDER) {
      _offsets[N-1] = 1;
      for (int i = N-2; i >= 0; --i) {
        _offsets[i] = _offsets[i+1] * extents[i+1];
      }
    } else if (NARRAY_MEMORDER_COL == MEMORDER){
      _offsets[0] = 1;
      for (int i = 1; i < N; ++i) {
        _offsets[i] = _offsets[i-1] * extents[i-1];
      }
    }
  }

  NArrayShape(self_t&&) = default;

  self_t&
  operator=(self_t&&) = default;

  NArrayShape(const self_t&) = default;

  self_t&
  operator=(const self_t&) = default;

  template<typename... Args>
  constexpr inline
  size_type
  idx(Args... args) const {
    static_assert(sizeof...(args) == N, "Underspecified indices detected!");
    size_type idx = get_idx(args...);
#ifdef ENABLE_DEBUG
    assert(idx >= 0 && idx < _size);
#endif // ENABLE_DEBUG
    return idx;
  }

  constexpr
  size_type size() const {
    return _size;
  }


private:

  static constexpr std::array<SizeT, N> bases = { Bases... };

  template<typename... Args>
  constexpr inline
  size_type get_idx(size_type idx, Args... args) const
  {
    if constexpr (N == 1) {
      // simplest case: only one dimension -> no offset needed
      return (idx-bases[0]);
    } else if constexpr (MEMORDER == NARRAY_MEMORDER_COL) {
      // column-major: first dimension does not need an offset
      return (idx-bases[0]) + get_idx_impl(args...);
    } else {
      constexpr int entry = N - sizeof...(Args) - 1;
      return _offsets[entry]*(idx-bases[entry]) + get_idx_impl(args...);
    }
  }


  template<typename... Args>
  constexpr inline
  size_type get_idx_impl(size_type idx, Args... args) const
  {
    constexpr int entry = N - sizeof...(Args) - 1;
    return _offsets[entry]*(idx-bases[entry]) + get_idx_impl(args...);
  }

  constexpr inline
  size_type get_idx_impl(size_type idx) const {
    // end condition
    if constexpr (MEMORDER == NARRAY_MEMORDER_ROW) {
      // for row-major, the lowest index has no offset
      return idx-bases[N-1];
    }
    return _offsets[N-1]*(idx-bases[N-1]);
  }

  template<typename... Args>
  constexpr static size_type total_size(Args... args) {
    size_t total_size = 1;
    std::array<size_type, N> extents{ static_cast<size_type>(args)... };
    for (int i = N-1; i >= 0; --i) {
      total_size *= extents[i];
    }
    return total_size;
  }

  size_type _size = 0;
  std::array<size_type, N> _offsets = {};
};



/**
 * Variant for which the extends are static, i.e., known at compile-time.
 * The extents are contained in the first \c N elements of the parameter pack
 * \c Values and the specification of the Index-Base in the second part.
 */
template<
  int N,
  typename SizeT,
  narray_memorder_t MEMORDER,
  SizeT... Values>
class StaticNArrayShape {

public:
  using size_type = SizeT;
  constexpr static const int ndims = N;

  using self_t = StaticNArrayShape<N, SizeT, MEMORDER, Values...>;

  constexpr
  StaticNArrayShape()
  { }

  StaticNArrayShape(self_t&&) = default;

  self_t&
  operator=(self_t&&) = default;

  template<typename... Args>
  constexpr inline
  size_type
  idx(Args... args) const {
    static_assert(sizeof...(args) == N, "Underspecified indices detected!");
    size_type idx = get_idx(args...);
#ifdef ENABLE_DEBUG
    assert(idx >= 0 && idx < _size);
#endif // ENABLE_DEBUG
    return idx;
  }

  static constexpr
  size_type size() {
    return _size;
  }

private:

  static constexpr
  std::array<SizeT, N> Extents =
          meta::take<N, SizeT, sizeof...(Values)>(
              std::array<SizeT, sizeof...(Values)>({Values...}));

  static constexpr
  std::array<SizeT, N> Bases =
          meta::drop<sizeof...(Values) - N, SizeT, sizeof...(Values)>(
              std::array<SizeT, sizeof...(Values)>({Values...}));

  static constexpr std::array<size_type, N> compute_offsets()
  {
    std::array<size_type, N> offsets{};
    if (NARRAY_MEMORDER_ROW == MEMORDER) {
      offsets[N-1] = 1;
      for (int i = N-2; i >= 0; --i) {
        offsets[i] = offsets[i+1] * Extents[i+1];
      }
    } else if (NARRAY_MEMORDER_COL == MEMORDER){
      offsets[0] = 1;
      for (int i = 1; i < N; ++i) {
        offsets[i] = offsets[i-1] * Extents[i-1];
      }
    }
    return offsets;
  }

  static constexpr std::array<size_type, N> offsets = compute_offsets();

  template<typename... Args>
  constexpr inline
  size_type get_idx(size_type idx, Args... args) const
  {
    if constexpr (N == 1) {
      // simplest case: only one dimension -> no offset needed
      return (idx-Bases[0]);
    } else if constexpr (MEMORDER == NARRAY_MEMORDER_COL) {
      // column-major: first dimension does not need an offset
      return (idx-Bases[0]) + get_idx_impl(args...);
    } else {
      constexpr int entry = N - sizeof...(Args) - 1;
      return offsets[entry]*(idx-Bases[entry]) + get_idx_impl(args...);
    }
  }


  template<typename... Args>
  constexpr inline
  size_type get_idx_impl(size_type idx, Args... args) const
  {
    constexpr int entry = N - sizeof...(Args) - 1;
    return offsets[entry]*(idx-Bases[entry]) + get_idx_impl(args...);
  }

  constexpr inline
  size_type get_idx_impl(size_type idx) const {
    // end condition
    if constexpr (MEMORDER == NARRAY_MEMORDER_ROW) {
      // for row-major, the lowest index has no offset
      return idx-Bases[N-1];
    }
    return offsets[N-1]*(idx-Bases[N-1]);
  }

  template<typename... Args>
  constexpr static size_type total_size() {
    size_t total_size = 1;
    for (int i = N-1; i >= 0; --i) {
      total_size *= Extents[i];
    }
    return total_size;
  }

  static constexpr size_type _size = total_size();

};


template<typename ValueT, typename ShapeT>
class NArrayMemoryHeap
{
public:

  using size_type  = typename ShapeT::size_type;
  using value_type = ValueT;
  using shape_type = ShapeT;
  constexpr static const int ndims = ShapeT::ndims;

  using self_t = NArrayMemoryHeap<value_type, shape_type>;

  NArrayMemoryHeap(void)
  { }

  template<typename... Args>
  constexpr
  NArrayMemoryHeap(const ShapeT& shape)
  : _mem(std::shared_ptr<value_type[]>(new value_type[shape.size()]))
  { }

  NArrayMemoryHeap(self_t&&) = default;

  self_t&
  operator=(self_t&&) = default;

  constexpr
  value_type* get()
  {
    return _mem.get();
  }

  constexpr
  const value_type* get() const
  {
    return _mem.get();
  }

  void
  allocate(const shape_type& shape){
    _mem = std::shared_ptr<value_type[]>(new value_type[shape.size()]);
  }


private:
  std::shared_ptr<value_type[]> _mem = {};
};

template<typename ValueT, typename ShapeT>
class NArrayMemoryStack
{

};

template<
  typename ValueT,
  int N,
  typename SizeT,
  narray_memorder_t MEMORDER,
  SizeT... Values>
class NArrayMemoryStack<ValueT, StaticNArrayShape<N, SizeT, MEMORDER, Values...>>
{
public:

  using shape_type = StaticNArrayShape<N, SizeT, MEMORDER, Values...>;
  using size_type  = typename shape_type::size_type;
  using value_type = ValueT;
  constexpr static const int ndims = shape_type::ndims;

  using self_t = NArrayMemoryStack<value_type, shape_type>;

  NArrayMemoryStack(void)
  { }

  constexpr
  NArrayMemoryStack(const shape_type&)
  { }

  // move construction is deleted
  NArrayMemoryStack(self_t&&) = default;

  self_t&
  operator=(self_t&&) = default;

  constexpr
  value_type* get()
  {
    return _mem.data();
  }

  constexpr
  const value_type* get() const
  {
    return _mem.data();
  }

  /*
   * Not provided by stack allocation
   */
  void allocate(const shape_type& shape) = delete;

private:
  std::array<value_type, shape_type::size()> _mem;
};


template<
  typename ValueT,
  int N,
  typename ShapeT = NArrayShape<N, size_t, NARRAY_MEMORDER_ROW>,
  typename MemoryT = NArrayMemoryHeap<ValueT, ShapeT>>
class NArray {
public:
  using size_type   = typename ShapeT::size_type;
  using value_type  = ValueT;
  using shape_type  = ShapeT;
  using memory_type = MemoryT;
  constexpr static const int ndims = N;

  using iterator = ValueT*;
  using const_iterator = const ValueT*;

  using self_t = NArray<value_type, N, shape_type, MemoryT>;

  NArray(void)
  { }

  template<typename... Args>
  constexpr
  NArray(Args... args)
  : _shape(args...),
    _mem(_shape)
  {
    static_assert(sizeof...(args) == N, "Underspecified extents detected!");
  }

  NArray(const NArray& narray) : _mem(narray._mem), _shape(narray._shape)
  { }

  NArray(NArray&&) = default;

  self_t&
  operator=(self_t&&) = default;

  self_t&
  operator=(self_t& narray) {
    if (&narray != this) {
      _mem   = narray._mem;
      _shape = narray._shape;
    }
    return *this;
  }

  template<typename... Args>
  void
  allocate(Args... args) {
    _shape = shape_type(args...);
    _mem.allocate(_shape);
  }


  template<typename... Args>
  constexpr inline
  value_type&
  operator()(Args... args) {
    size_type idx = _shape.idx(args...);
    return _mem.get()[idx];
  }


  template<typename... Args>
  constexpr inline
  const value_type&
  operator()(Args... args) const {
    size_type idx = _shape.idx(args...);
    return _mem.get()[idx];
  }

  constexpr inline
  const_iterator
  begin() const {
    return _mem.get();
  }

  constexpr inline
  iterator
  begin() {
    return _mem.get();
  }

  constexpr inline
  const_iterator
  end() const {
    return _mem.get() + _shape.size();
  }

  constexpr inline
  iterator
  end() {
    return _mem.get() + _shape.size();
  }

  const shape_type&
  shape() const {
    return _shape;
  }

private:

  template<typename... Args>
  constexpr static size_type total_size(Args... args) {
    size_t total_size = 1;
    std::array<size_type, N> extents{ static_cast<size_type>(args)... };
    for (int i = N-1; i >= 0; --i) {
      total_size *= extents[i];
    }
    return total_size;
  }

  shape_type       _shape = {};
  memory_type      _mem   = {};
};


template<
  typename ValueT,
  int N,
  typename ShapeT = NArrayShape<N, size_t, NARRAY_MEMORDER_ROW>>
class NArrayView {

public:

  static_assert(N == ShapeT::ndims, "Unequal dimensions in NArrayView and NArrayShape!");

  using size_type  = typename ShapeT::size_type;
  using value_type = ValueT;
  using shape_type = ShapeT;
  constexpr static const int ndims = N;

  using self_t     = NArrayView<ValueT, N, ShapeT>;

  template<typename... Args>
  NArrayView(value_type* mem, Args ...args) : _mem(mem), _shape(args...)
  { }

  template<
    typename ArrayShapeT,
    typename MemoryT = NArrayMemoryHeap<ValueT, ArrayShapeT>>
  NArrayView(NArray<ValueT, N, ArrayShapeT, MemoryT>& narray)
  : _mem(narray.begin()), _shape(narray.shape())
  { }

  NArrayView(self_t&&) = default;

  NArrayView(self_t&)  = default;

  self_t&
  operator=(self_t&&)  = default;

  self_t&
  operator=(self_t&)   = default;

  template<typename... Args>
  constexpr inline
  ValueT&
  operator()(Args... args) {
    size_type idx = _shape.idx(args...);
    return _mem[idx];
  }


  template<typename... Args>
  constexpr inline
  const ValueT&
  operator()(Args... args) const {
    size_type idx = _shape.idx(args...);
    return _mem[idx];
  }



private:

  value_type *_mem   = nullptr;
  const shape_type  _shape;

};


/************************************************************************
 * Some types that are commonly used throughout the code.               *
 ************************************************************************/

/* 1D column-major Array, 1-index-based */
template<typename ValueT, typename SizeT>
using ArrayViewT1 = NArrayView<ValueT, 1, NArrayShape<1, SizeT, NARRAY_MEMORDER_COL, 1>>;
/* 3D column-major Array, 0-index-based */
template<typename ValueT, typename SizeT>
using ArrayViewT3 = NArrayView<ValueT, 3, NArrayShape<3, SizeT, NARRAY_MEMORDER_COL, 0, 0, 0>>;
/* 4D column-major Array, 0-index-based except on the fastest dimension */
template<typename ValueT, typename SizeT>
using ArrayViewT4 = NArrayView<ValueT, 4, NArrayShape<4, SizeT, NARRAY_MEMORDER_COL, 1, 0, 0, 0>>;


/* 1D static array allocated on the stack */
template<
  typename ValueT,
  typename SizeT,
  SizeT    Base0 = 1,
  typename ShapeT = NArrayShape<1, SizeT, NARRAY_MEMORDER_COL, Base0>>
using ArrayT1 = NArray<ValueT, ShapeT::ndims, ShapeT, NArrayMemoryHeap<ValueT, ShapeT>>;


/* 1D static array allocated on the stack */
template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Base0 = 1,
  typename ShapeT = StaticNArrayShape<1, SizeT, NARRAY_MEMORDER_COL, Dim0, Base0>>
using StaticArrayT1 = NArray<ValueT, ShapeT::ndims, ShapeT, NArrayMemoryStack<ValueT, ShapeT>>;

template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Dim1,
  SizeT    Base0 = 1,
  SizeT    Base1 = 1,
  typename ShapeT = StaticNArrayShape<2, SizeT, NARRAY_MEMORDER_COL, Dim0, Dim1, Base0, Base1>>
using StaticArrayT2 = NArray<ValueT, ShapeT::ndims, ShapeT, NArrayMemoryStack<ValueT, ShapeT>>;

template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Dim1,
  SizeT    Dim2,
  SizeT    Base0 = 1,
  SizeT    Base1 = 1,
  SizeT    Base2 = 0,
  typename ShapeT = StaticNArrayShape<3, SizeT, NARRAY_MEMORDER_COL, Dim0, Dim1, Dim2, Base0, Base1, Base2>>
using StaticArrayT3 = NArray<ValueT, ShapeT::ndims, ShapeT, NArrayMemoryStack<ValueT, ShapeT>>;


template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Dim1,
  SizeT    Dim2,
  SizeT    Dim3,
  SizeT    Base0 = 1,
  SizeT    Base1 = 1,
  SizeT    Base2 = 1,
  SizeT    Base3 = 0,
  typename ShapeT = StaticNArrayShape<4, SizeT, NARRAY_MEMORDER_COL, Dim0, Dim1, Dim2, Dim3, Base0, Base1, Base2, Base3>>
using StaticArrayT4 = NArray<ValueT, ShapeT::ndims, ShapeT, NArrayMemoryStack<ValueT, ShapeT>>;


template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Base0 = 1,
  typename ShapeT = StaticNArrayShape<1, SizeT, NARRAY_MEMORDER_COL, Dim0, Base0>>
using StaticArrayViewT1 = NArrayView<ValueT, ShapeT::ndims, ShapeT>;

template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Dim1,
  SizeT    Base0 = 1,
  SizeT    Base1 = 1,
  typename ShapeT = StaticNArrayShape<2, SizeT, NARRAY_MEMORDER_COL, Dim0, Dim1, Base0, Base1>>
using StaticArrayViewT2 = NArrayView<ValueT, ShapeT::ndims, ShapeT>;


template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Dim1,
  SizeT    Dim2,
  SizeT    Base0 = 1,
  SizeT    Base1 = 1,
  SizeT    Base2 = 0,
  typename ShapeT = StaticNArrayShape<3, SizeT, NARRAY_MEMORDER_COL, Dim0, Dim1, Dim2, Base0, Base1, Base2>>
using StaticArrayViewT3 = NArrayView<ValueT, ShapeT::ndims, ShapeT>;


template<
  typename ValueT,
  typename SizeT,
  SizeT    Dim0,
  SizeT    Dim1,
  SizeT    Dim2,
  SizeT    Dim3,
  SizeT    Base0 = 1,
  SizeT    Base1 = 1,
  SizeT    Base2 = 1,
  SizeT    Base3 = 0,
  typename ShapeT = StaticNArrayShape<4, SizeT, NARRAY_MEMORDER_COL, Dim0, Dim1, Dim2, Dim3, Base0, Base1, Base2, Base3>>
using StaticArrayViewT4 = NArrayView<ValueT, ShapeT::ndims, ShapeT>;


#endif // HAVE_NARRAY_H
