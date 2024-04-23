#ifndef MATRIX_ROW_IMPLEMENTATION_HPP
#define MATRIX_ROW_IMPLEMENTATION_HPP

#include "Matrix.hpp"

/**
 * @brief Helper function to compress a matrix to the CSR-format. NOTE:
 * contrastingly to the suggested implementation this method does not rely on
 * the lower_bound() and upper_bound() methods. Instead, we only use one
 * for-loop and no conditional jumps by exploiting the internal ordering of the
 * dict. The row-counter is always correctly updated since the mapping is
 * ordered based on the rows.
 *
 * @tparam T Type of the matrix entries.
 * @tparam Store StorageOrder, either row or col.
 */
template <class T, StorageOrder Store>
void Matrix<T, Store>::_compress_row() {
#ifdef DEBUG
  std::cout << "Using ROW-MAJOR compression to CSR.\n";
#endif
  // vec1 of length #rows + 1 -> row indices
  // vec2 of length #non-zero-elements -> column index
  // _values: length #non-zero-elements -> actual values

  // #rows + 1 = highest row-number + 2
  std::size_t num_rows = _entry_value_map.rbegin()->first[0] + 2;
  _vec1.resize(num_rows, 0);

  // number of non-zeros are simply the number of map entries
  std::size_t num_non_zeros = _entry_value_map.size();
  _vec2.resize(num_non_zeros);
  _values.resize(num_non_zeros);

#ifdef DEBUG
  std::cout << "num_rows = " << num_rows << "\n";
  std::cout << "vec2.size() = " << _vec2.size() << "\n";
  std::cout << "values.size() = " << _values.size() << "\n";
#endif

  std::size_t num_non_zero = 0;
  // idea: not use conditional jumps
  for (const auto& [k, v] : _entry_value_map) {
    _vec2[num_non_zero] = k[1];  // add the column index
    _values[num_non_zero] = v;   // add the value
    // we just update the count of non-zeros at the curr. row-idx
    // note that the row-idx is automatically incremented
    _vec1[k[0] + 1] = ++num_non_zero;
  }

  // save memory and set flags
  _is_compressed = true;
  _entry_value_map.clear();
}

/**
 * @brief Map the entries from the vectors back to the mapping.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 */
template <class T, StorageOrder Storage>
void Matrix<T, Storage>::_uncompress_row() {
#ifdef DEBUG
  std::cout << "Using ROW-MAJOR uncompression\n";
#endif
  // vec1 of length #rows + 1 -> row indices
  // vec2 of length #non-zero-elements -> column index
  // _values: length #non-zero-elements -> actual values
  std::size_t num_rows = _vec1.size() - 1;
  for (std::size_t row_idx = 0; row_idx < num_rows; ++row_idx) {
    for (std::size_t col_idx = _vec1[row_idx]; col_idx < _vec1[row_idx + 1];
         ++col_idx) {
      // we get the col number and the value accordingly
      _entry_value_map[std::array<std::size_t, 2>{row_idx, _vec2[col_idx]}] =
          _values[col_idx];
    }
  }
  // save memory and set flags
  _is_compressed = false;
  _vec1.clear();
  _vec2.clear();
  _values.clear();
}

/**
 * @brief Helper method to find a compressed element in the case of
 * row-compression.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @param row Row index.
 * @param col Column index.
 * @return const T Value of the element
 */
template <class T, StorageOrder Store>
const T Matrix<T, Store>::_find_compressed_element_row(std::size_t row,
                                                       std::size_t col) const {
#ifdef DEBUG
  std::cout << "Using ROW-MAJOR _find_compressed_element() const version.\n ";
#endif
  for (std::size_t col_idx = _vec1[row]; col_idx < _vec1[row + 1]; ++col_idx) {
    if (_vec2[col_idx] == col) {
#ifdef DEBUG
      std::cout << "Found element: " << row << ", " << col << ".\n";
#endif
      return _values[col_idx];
    }
  }
  return 0;
}

/**
 * @brief Helper method to implement the non-const version of the setter/getter
 * method in the case of row-compression.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @param row Row index.
 * @param col Column index.
 * @return &T Reference to the element.
 */
template <class T, StorageOrder Store>
T& Matrix<T, Store>::_find_compressed_element_row(std::size_t row,
                                                  std::size_t col) {
#ifdef DEBUG
  std::cout
      << "Using ROW-MAJOR _find_compressed_element() non-const version.\n";
#endif
  for (std::size_t col_idx = _vec1[row]; col_idx < _vec1[row + 1]; ++col_idx) {
    if (_vec2[col_idx] == col) {
#ifdef DEBUG
      std::cout << "Found element: " << row << ", " << col << ".\n";
#endif
      return _values[col_idx];
    }
  }
  throw std::invalid_argument(
      "Trying to modify a zero-element in compressed format. Uncompress first");
}

/**
 * @brief Matrix-vector multiplication for the row-compression case.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @param vec Vector x to compute A*x.
 * @return std::vector<T> Result of the Matrix-vector multiplication.
 */
template <class T, StorageOrder Store>
std::vector<T> Matrix<T, Store>::_matrix_vector_row(std::vector<T> vec) const {
  // iterate through the rows, then the elements
  std::vector<T> res;
  res.resize(_vec1.size() - 1, 0);

  for (std::size_t row_idx = 0; row_idx < _vec1.size() - 1; ++row_idx) {
    // get the columns, according to this row
    for (std::size_t col_idx = _vec1[row_idx]; col_idx < _vec1[row_idx + 1];
         ++col_idx) {
      res[row_idx] += (vec[_vec2[col_idx]] * _values[col_idx]);
    }
  }
  return res;
}

/**
 * @brief Compute the max-norm for the row-compression case.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @return T Norm of the matrix.
 */
template <class T, StorageOrder Store>
T Matrix<T, Store>::_max_norm_compressed_row() const {
#ifdef DEBUG
  std::cout << "Max-Norm compressed-ROW.\n";
#endif
  T res = 0;
  for (std::size_t row_idx = 0; row_idx < _vec1.size() - 1; ++row_idx) {
    // get the columns, according to this row
    T norm_curr_row = 0;
    for (std::size_t col_idx = _vec1[row_idx]; col_idx < _vec1[row_idx + 1];
         ++col_idx) {
      norm_curr_row += std::abs(_values[col_idx]);
    }
    res = std::max(res, norm_curr_row);
  }
  return res;
}

/**
 * @brief Compute the one-norm for the col-compression case.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @return T Norm of the matrix.
 */
template <class T, StorageOrder Store>
T Matrix<T, Store>::_one_norm_compressed_row() const {
#ifdef DEBUG
  std::cout << "One-Norm compressed-ROW.\n";
#endif
  std::size_t num_cols = *max_element(std::begin(_vec2), std::end(_vec2));
  std::vector<T> sum_abs_per_col(num_cols, 0);
  for (std::size_t col_idx = 0; col_idx < _vec2.size(); ++col_idx) {
    sum_abs_per_col[_vec2[col_idx]] += std::abs(_values[col_idx]);
  }
  return *max_element(std::begin(sum_abs_per_col), std::end(sum_abs_per_col));
}
#endif