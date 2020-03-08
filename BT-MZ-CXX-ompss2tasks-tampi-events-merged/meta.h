#ifndef HAVE_META_H
#define HAVE_META_H

#include <array>
#include <initializer_list>

namespace meta {

// -------------------------------------------------------------------------
// drop
// -------------------------------------------------------------------------

namespace detail {
  template <
    std::size_t    NDrop,
    typename       ValueT,
    std::size_t    NElem,
    std::size_t... Is >
  constexpr std::array<ValueT, (NDrop > NElem) ? 0 : NElem - NDrop >
  drop_impl(
    const std::array<ValueT, NElem> & values,
    std::index_sequence<Is...>) {
    return {{ std::get<NDrop + Is>(values)... }};
  }
} // namespace detail

/**
 * Drops \c d elements from a given sequence of \c N elements with indices
 * \c (0..d..N).
 */
template <
  std::size_t   NDrop,
  typename      ValueT,
  std::size_t   NElem >
constexpr std::array<ValueT, (NDrop > NElem) ? 0 : NElem - NDrop >
drop(
  const std::array<ValueT, NElem> & values) {
  return detail::drop_impl<NDrop, ValueT, NElem>(
           values,
           std::make_index_sequence<(NDrop > NElem) ? 0 : NElem - NDrop>()
          );
}


// -------------------------------------------------------------------------
// take
// -------------------------------------------------------------------------

namespace detail {
  template <
    std::size_t    NTake,
    typename       ValueT,
    std::size_t    NElem,
    std::size_t... Is >
  constexpr std::array<ValueT, (NTake > NElem) ? NElem : NTake >
  take_impl(
    const std::array<ValueT, NElem> & values,
    std::index_sequence<Is...>) {
    return {{ std::get<Is>(values)... }};
  }
} // namespace detail

/**
 * Returns sequence of first \c t elements from a given sequence of size
 * \c N with indices \c (0..t..N).
 */
template <
  std::size_t   NTake,
  typename      ValueT,
  std::size_t   NElem >
constexpr std::array<ValueT, (NTake > NElem) ? NElem : NTake >
take(
  const std::array<ValueT, NElem> & values) {
  return detail::take_impl<NTake, ValueT, NElem>(
           values,
           std::make_index_sequence<
             (NTake > NElem) ? NElem : NTake
           >());
}


// -------------------------------------------------------------------------
// split
// -------------------------------------------------------------------------

template <
  class         ValueT,
  std::size_t   NElemLeft,
  std::size_t   NElemRight >
class split
{
  typedef meta::split<ValueT, NElemLeft, NElemRight> self_t;

  constexpr static std::size_t NElem = NElemLeft + NElemRight;

  // Caveat: copies array values in non-constexpr use cases?
  const std::array<ValueT, NElem> _values;

public:
  constexpr split(
    const std::initializer_list<ValueT> & values)
    : _values(values)
  { }

  constexpr split(
    const std::array<ValueT, NElem> & values)
    : _values(values)
  { }

  constexpr std::array<ValueT, NElemLeft> left() const {
    return take<NElemLeft, ValueT, NElem>(_values);
  }

  constexpr std::array<ValueT, NElemRight> right() const {
    return drop<NElemLeft, ValueT, NElem>(_values);
  }
};

} // namespace meta

#endif // HAVE_META_H
