
#ifndef MATRIX_COL_IMPLEMENTATION_HPP
#define MATRIX_COL_IMPLEMENTATION_HPP
#include "Matrix.hpp"

/**
 * @brief Compress a matrix to CSC. NOTE: contrastingly to the suggested
 * implementation this method does not rely on the lower_bound() and
 * upper_bound() methods. Instead, we only use one for-loop and no conditional
 * jumps by exploiting the internal ordering of the dict. The row-counter is
 * always correctly updated since the mapping is ordered based on the rows.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 */
template <class T, StorageOrder Store>
void Matrix<T, Store>::_compress_col() {
#ifdef DEBUG
  std::cout << "Using COL-MAJOR compression to CSC.\n";
#endif
  // NOTE(Assumption): we assume the dict is orderd in (second,first) way
  // so (1,2) is after (2,1) in the ordering

  // vec1 of length #cols + 1 -> col indices
  // vec2 of length #non-zero-elements -> row index
  // _values: length #non-zero-elements -> actual values

  // #cols = highest col-number + 2
  std::size_t num_cols = _entry_value_map.rbegin()->first[1] + 2;
  _vec1.resize(num_cols, 0);

  // number of non-zeros are simply the number of map entries
  std::size_t num_non_zeros = _entry_value_map.size();
  _vec2.resize(num_non_zeros);
  _values.resize(num_non_zeros);

#ifdef DEBUG
  std::cout << "num_rows = " << num_cols << "\n";
  std::cout << "vec2.size() = " << _vec2.size() << "\n";
  std::cout << "values.size() = " << _values.size() << "\n";
#endif

  std::size_t num_non_zero = 0;
  // idea: not use conditional jumps
  for (const auto& [k, v] : _entry_value_map) {
    _vec2[num_non_zero] = k[0];  // add the column index
    _values[num_non_zero] = v;   // add the value
    // we just update the count of non-zeros at the curr. col-idx
    // note that the col-idx is automatically incremented
    _vec1[k[1] + 1] = ++num_non_zero;
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
template <class T, StorageOrder Store>
void Matrix<T, Store>::_uncompress_col() {
#ifdef DEBUG
  std::cout << "Using COL-MAJOR uncompression.\n";
#endif
  // vec1 of length #cols + 1 -> col indices
  // vec2 of length #non-zero-elements -> row index
  // _values: length #non-zero-elements -> actual values
  std::size_t num_cols = _vec1.size() - 1;
  for (std::size_t col_idx = 0; col_idx < num_cols; ++col_idx) {
    for (std::size_t row_idx = _vec1[col_idx]; row_idx < _vec1[col_idx + 1];
         ++row_idx) {
      // we get the col number and the value accordingly
      _entry_value_map[std::array<std::size_t, 2>{_vec2[row_idx], col_idx}] =
          _values[row_idx];
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
 * column-compression.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @param row Row index.
 * @param col Column index.
 * @return const T Value of the element
 */
template <class T, StorageOrder Store>
const T Matrix<T, Store>::_find_compressed_element_col(std::size_t row,
                                                       std::size_t col) const {
#ifdef DEBUG
  std::cout << "Using COL-MAJOR _find_compressed_element() const version.\n";
#endif
  for (std::size_t row_idx = _vec1[col]; row_idx < _vec1[col + 1]; ++row_idx) {
    if (_vec2[row_idx] == row) {
#ifdef DEBUG
      std::cout << "Found element.";
#endif
      return _values[row_idx];
    }
  }
  return 0;
}

/**
 * @brief Helper method to implement the non-const version of the setter/getter
 * method in the case of column-compression.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @param row Row index.
 * @param col Column index.
 * @return &T Reference to the element.
 */
template <class T, StorageOrder Store>
T& Matrix<T, Store>::_find_compressed_element_col(std::size_t row,
                                                  std::size_t col) {
#ifdef DEBUG
  std::cout
      << "Using COL-MAJOR _find_compressed_element() non-const version.\n";
#endif
  for (std::size_t row_idx = _vec1[col]; row_idx < _vec1[col + 1]; ++row_idx) {
    if (_vec2[row_idx] == row) {
#ifdef DEBUG
      std::cout << "Found element.";
#endif
      return _values[row_idx];
    }
  }
  throw std::invalid_argument(
      "Trying to modify a zero-element in compressed format. Uncompress "
      "first");
}

/**
 * @brief Matrix-vector multiplication for the col-compression case.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @param vec Vector x to compute A*x.
 * @return std::vector<T> Result of the Matrix-vector multiplication.
 */
template <class T, StorageOrder Store>
std::vector<T> Matrix<T, Store>::_matrix_vector_col(std::vector<T> vec) const {
  std::vector<T> res;
  // #rows = max value in the row-index vector
  std::size_t num_rows = *max_element(_vec2.begin(), _vec2.end());
  res.resize(num_rows + 1, 0);
  // iterate through the colums
  for (int col_idx = 0; col_idx < _vec1.size() - 1; ++col_idx) {
    for (std::size_t row_idx = _vec1[col_idx]; row_idx < _vec1[col_idx + 1];
         ++row_idx) {
      res[_vec2[row_idx]] += (vec[col_idx] * _values[row_idx]);
    }
  }
  return res;
}

/**
 * @brief Compute the max-norm for the column-compression case.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @return T Norm of the matrix.
 */
template <class T, StorageOrder Store>
T Matrix<T, Store>::_max_norm_compressed_col() const {
#ifdef DEBUG
  std::cout << "Max-Norm compressed-COL.\n";
#endif
  std::size_t num_rows = *max_element(std::begin(_vec2), std::end(_vec2));
  std::vector<T> sum_abs_per_col(num_rows, 0);
  for (std::size_t row_idx = 0; row_idx < _vec2.size(); ++row_idx) {
    sum_abs_per_col[_vec2[row_idx]] += std::abs(_values[row_idx]);
  }
  return *max_element(std::begin(sum_abs_per_col), std::end(sum_abs_per_col));
}

/**
 * @brief Compute the one-norm for the column-compression case.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 * @return T Norm of the matrix.
 */
template <class T, StorageOrder Store>
T Matrix<T, Store>::_one_norm_compressed_col() const {
#ifdef DEBUG
  std::cout << "One-Norm compressed-COLROW.\n";
#endif
  T res = 0;
  for (std::size_t col_idx = 0; col_idx < _vec1.size() - 1; ++col_idx) {
    // get the row, according to this col
    T norm_curr_col = 0;
    for (std::size_t row_idx = _vec1[col_idx]; row_idx < _vec1[col_idx + 1];
         ++row_idx) {
      norm_curr_col += std::abs(_values[row_idx]);
    }
    res = std::max(res, norm_curr_col);
  }
  return res;
}
#endif