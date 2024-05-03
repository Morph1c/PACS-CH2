#ifndef MATRIX_ROW_IMPLEMENTATION_HPP
#define MATRIX_ROW_IMPLEMENTATION_HPP

#include "Matrix.hpp"

/**
 *  Since we are in row-order first we have element with first row and so on
 * The first (the inner indexes), of length the number of rows plus one, contains the starting index (in the values array) for the elements of each row
 * The second vector of indexes (the outer indexes), of length the number of non-zeroes, contains the corresponding column index.
 * @tparam T Type of the matrix entries.
 * @tparam Store StorageOrder, either row or col.
 */
template <Numeric T, StorageOrder Store>
void Matrix<T, Store>::_compress_row() {

  // vec1 of length #rows + 1 -> row indices
  // vec2 of length #non-zero-elements -> column index
  // _values: length #non-zero-elements -> actual values

  // since we are ordering by rows then the last element of the entry_value_map is the highest row number
  // so we have to add 2 since it start with 0 to get the number of rows
  std::size_t num_rows = _entry_value_map.rbegin()->first[0] + 2;
  _inner.resize(num_rows, 0);

  // number of non-zeros are simply the number of map entries
  std::size_t num_non_zeros = _entry_value_map.size();
  _outer.resize(num_non_zeros);
  _values.resize(num_non_zeros);

  std::size_t num_non_zero = 0;

  // implement using lower and upper bound
  for (std::size_t row = 0; row < num_rows - 1; ++row){
    auto low = _entry_value_map.lower_bound({row, std::numeric_limits<std::size_t>::min()});
    auto up = _entry_value_map.upper_bound({row, std::numeric_limits<std::size_t>::max()}); 

    for (auto it = low; it != up; ++it) {
      _outer[num_non_zero] = it->first[1];  // add the column index
      _values[num_non_zero] = it->second;   // add the value
     ++num_non_zero;
    }
    _inner[row + 1] = num_non_zero; // since the next row the non-zero element start from this pos.
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
template <Numeric T, StorageOrder Storage>
void Matrix<T, Storage>::_uncompress_row() {

  // vec1 of length #rows + 1 -> row indices
  // vec2 of length #non-zero-elements -> column index
  // _values: length #non-zero-elements -> actual values
  std::size_t num_rows = _inner.size() - 1;
  for (std::size_t row_idx = 0; row_idx < num_rows; ++row_idx) {
    for (std::size_t col_idx = _inner[row_idx]; col_idx < _inner[row_idx + 1];
         ++col_idx) {
      // we get the col number and the value accordingly
      _entry_value_map[std::array<std::size_t, 2>{row_idx, _outer[col_idx]}] =
          _values[col_idx];
    }
  }
  // save memory and set flags
  _is_compressed = false;
  _inner.clear();
  _outer.clear();
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
template <Numeric T, StorageOrder Store>
const T Matrix<T, Store>::_find_compressed_element_row(std::size_t row,
                                                       std::size_t col) const {

  for (std::size_t col_idx = _inner[row]; col_idx < _inner[row + 1]; ++col_idx) {
    if (_outer[col_idx] == col) {

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
template <Numeric T, StorageOrder Store>
T& Matrix<T, Store>::_find_compressed_element_row(std::size_t row,
                                                  std::size_t col) {

  for (std::size_t col_idx = _inner[row]; col_idx < _inner[row + 1]; ++col_idx) {
    if (_outer[col_idx] == col) {

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
template <Numeric T, StorageOrder Store>
std::vector<T> Matrix<T, Store>::_matrix_vector_row(std::vector<T> vec) const {
  // iterate through the rows, then the elements
  std::vector<T> res;
  res.resize(_inner.size() - 1, 0);

  for (std::size_t row_idx = 0; row_idx < _inner.size() - 1; ++row_idx) {
    // get the columns, according to this row
    for (std::size_t col_idx = _inner[row_idx]; col_idx < _inner[row_idx + 1];
         ++col_idx) {
      res[row_idx] += (vec[_outer[col_idx]] * _values[col_idx]);
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
template <Numeric T, StorageOrder Store>
T Matrix<T, Store>::_max_norm_compressed_row() const {

  T res = 0;
  for (std::size_t row_idx = 0; row_idx < _inner.size() - 1; ++row_idx) {
    // get the columns, according to this row
    T norm_curr_row = 0;
    for (std::size_t col_idx = _inner[row_idx]; col_idx < _inner[row_idx + 1];
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
template <Numeric T, StorageOrder Store>
T Matrix<T, Store>::_one_norm_compressed_row() const {

  std::size_t num_cols = *max_element(std::begin(_outer), std::end(_outer));
  std::vector<T> sum_abs_per_col(num_cols, 0);
  for (std::size_t col_idx = 0; col_idx < _outer.size(); ++col_idx) {
    sum_abs_per_col[_outer[col_idx]] += std::abs(_values[col_idx]);
  }
  return *max_element(std::begin(sum_abs_per_col), std::end(sum_abs_per_col));
}
#endif