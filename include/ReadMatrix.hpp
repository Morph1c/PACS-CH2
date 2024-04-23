#ifndef READ_MATRIX_HPP
#define READ_MATRIX_HPP

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "Utilities.hpp"

namespace algebra {

/**
 * @brief Read a matrix in the matrix-market format.
 *
 * @tparam T Type of the matrix entries.
 * @tparam Store StorageOrder for the matrix, deciding the ordering of the
 * mapping.
 * @param file_name Path to the matrix-market file.
 * @return std::map<std::array<std::size_t, 2>, T,
 * std::conditional_t<Store == StorageOrder::row, RowOrderComparator<T>,
 * ColOrderComparator<T>>> Mapping "(row, col) -> value" which can be directly
 * passed into the constructor.
 */
template <class T, StorageOrder Store>
std::map<std::array<std::size_t, 2>, T,
         std::conditional_t<Store == StorageOrder::row, RowOrderComparator<T>,
                            ColOrderComparator<T>>>
read_matrix(const std::string& file_name) {
  // Define the type of the map based on storage order, i.e. use different
  // comparison operators
  using MapType = std::map<
      std::array<std::size_t, 2>, T,
      std::conditional_t<Store == StorageOrder::row, RowOrderComparator<T>,
                         ColOrderComparator<T>>>;

  std::ifstream file(file_name);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + file_name);
  }

  std::string line;

  while (file.peek() == '%') {
    file.ignore(2048, '\n');
  }
  // std::istringstream iss(line);
  std::size_t num_rows, num_cols, num_elements;

  file >> num_rows >> num_cols >> num_elements;
#ifdef DEBUG
  std::cout << "num_elements: " << num_elements << "\n";
  std::cout << "num_rows: " << num_rows << "\n";
  std::cout << "num_cols: " << num_cols << "\n";
#endif

  MapType entry_value_map;

  for (std::size_t i = 0; i < num_elements; ++i) {
    std::size_t row, col;
    T value;
    file >> row >> col >> value;

#ifdef DEBUG
    std::cout << row << ", " << col << ": " << value << "\n";
#endif
    // we always use the format (row, col) -> value
    // only the comparison operator is different
    entry_value_map[{row - 1, col - 1}] = value;
  }
  return entry_value_map;
}
}  // namespace algebra

#endif