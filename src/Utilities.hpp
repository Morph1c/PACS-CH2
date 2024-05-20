#ifndef UTILITY_HPP
#define UTILITY_HPP
// clang-format off
#include <iostream>
#include <random>
#include <vector>
#include <concepts>
#include <type_traits>


namespace algebra {
/**
 * @brief Enum class to define the storage order of the matrix, the type of norm to be computed and a concepts for accepting
 * only numeric types.
 * 
 */
enum StorageOrder { row, col };

enum NormOrder { frob, one, max };

template <typename T>
concept Numeric = std::is_same_v<T, float> || std::is_same_v<T, double>;;

/**
 * @brief Introduce an ordering relation for an arrays of two entries.
 *        This can be only a partial ordering relation
 * @tparam T is the type of the entries of the array
 */
template <Numeric T>
struct ColOrderComparator {
  bool operator()(const std::array<std::size_t, 2>& x,
                  const std::array<std::size_t, 2>& y) const {

    return (x[1] < y[1]) || ((x[1] == y[1]) && (x[0] < y[0]));
  }
};

/**
 * @brief Default ordering relation for std::map but reimplementing for sake of semplicity
 * 
 * @tparam T Type of the entries.
 */
template <Numeric T>
struct RowOrderComparator {
  bool operator()(const std::array<std::size_t, 2>& x,
                  const std::array<std::size_t, 2>& y) const {
    
    return (x[0] < y[0]) || ((x[0] == y[0]) && (x[1] < y[1]));
  }
};

/**
 * @brief Utility function to generate a random vector of a fixed size
 * for the benchmark test.
 * 
 * @tparam T Type of the entries.
 * @param size Size of the vector.
 * @param lower_bound Lower bounds to draw the values from
 * @param upper_bound Upper bounds to draw the values from
 * @return std::vector<T> Random vector.
 */
template <Numeric T>
std::vector<T> _generate_random_vector(std::size_t size,
                                       double lower_bound = -10,
                                       double upper_bound = 10) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<T> dis(lower_bound, upper_bound);

  std::vector<T> random_vector(size);
//@note avoid annoying warning due to signed/unsigned comparison:
// for (int i = 0u; i < size; ++i) {
// solves the problem
  for (int i = 0; i < size; ++i) {
    random_vector[i] = dis(gen);
  }
  return random_vector;
}
}  // namespace algebra
#endif