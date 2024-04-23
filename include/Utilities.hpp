#ifndef STORAGE_ORDER_HPP
#define STORAGE_ORDER_HPP
#include <iostream>
#include <random>
#include <vector>
namespace algebra {
enum StorageOrder { row, col };
enum NormOrder { frob, one, max };

/**
 * @brief Compare two arrays of size 2 based on their second argument.
 * 
 * @tparam T Type of the entries.
 */
template <typename T>
struct ColOrderComparator {
  bool operator()(const std::array<std::size_t, 2>& a,
                  const std::array<std::size_t, 2>& b) const {
    // #ifdef DEBUG
    //     std::cout << "Col order comparator\n";
    // #endif
    return (a[1] < b[1]) || ((a[1] == b[1]) && (a[0] < b[0]));
  }
};

/**
 * @brief Default, reimplementing lexicographical sorting, just for the sake
 * of completeness.
 * 
 * @tparam T Type of the entries.
 */
template <typename T>
struct RowOrderComparator {
  bool operator()(const std::array<std::size_t, 2>& a,
                  const std::array<std::size_t, 2>& b) const {
    // #ifdef DEBUG
    //     std::cout << "Row compator is used\n";
    // #endif
    return (a[0] < b[0]) || ((a[0] == b[0]) && (a[1] < b[1]));
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
template <class T>
std::vector<T> _generate_random_vector(std::size_t size,
                                       double lower_bound = -10,
                                       double upper_bound = 10) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<T> dis(lower_bound, upper_bound);

  std::vector<T> random_vector(size);
  for (int i = 0; i < size; ++i) {
    random_vector[i] = dis(gen);
  }
  return random_vector;
}
}  // namespace algebra
#endif