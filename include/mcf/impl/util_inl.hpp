// vim: expandtab:ts=2:sw=2
#ifndef MCF_UTIL_INL_HPP
#define MCF_UTIL_INL_HPP

#include <mcf/logging.hpp>

namespace mcf {

namespace internal {

template <typename F, typename Compare>
std::result_of_t<F(int)> BinarySearch(F f, int min_bound, int max_bound,
                                      Compare comparator) {
  ASSERT_TRUE(max_bound - min_bound >= 0,
              "Invalid range given (min_bound > max_bound)");
  using ResultType = std::result_of_t<F(int)>;
  print("Running binary search in [", min_bound, ", ", max_bound, "]");

  // Initialize function evalation values.
  ResultType result_at_min = f(min_bound);
  ResultType result_at_max = f(max_bound);

  // Run search until minimum has been found.
  while (min_bound < max_bound) {
    const int center = (min_bound + max_bound) / 2;
    const int next = center + 1;

    ResultType result_at_center = f(center);
    ResultType result_at_next = f(next);

    if (comparator(result_at_next, result_at_center)) {
      // Derivative is smaller than 0, minimum is in range [next, max_bound].
      min_bound = next;
      result_at_min = result_at_next;
      print("New bounds [", min_bound, ", ", max_bound, "]");
    } else {
      // Derivative is greater or equal to 0, minimum is in
      // range [min_bound, center]. This also catches infeasible solutions
      // with infinity cost.
      max_bound = center;
      result_at_max = result_at_center;
      print("New bounds [", min_bound, ", ", max_bound, "]");
    }
  }

  // Bounds converged to a single value, return result.
  return result_at_min;
}

}  // namespace internal

}  // namespace mcf

#endif
