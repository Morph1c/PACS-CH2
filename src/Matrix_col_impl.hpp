
#ifndef MATRIX_COL_IMPLEMENTATION_HPP
#define MATRIX_COL_IMPLEMENTATION_HPP
#include "Matrix.hpp"

/**
 * Since we are in col-order first we have element with first col and so on
 * The first (the inner indexes), of length the number of rows plus one, contains the starting index (in the values array) for the elements of each col
 * The second vector of indexes (the outer indexes), of length the number of non-zeroes, contains the corresponding row index.
 *
 * @tparam T Type of the entries.
 * @tparam Store Storage order.
 */
template <class T, StorageOrder Store>
void Matrix<T, Store>::_compress_col() {

  // vec1 of length #cols + 1 -> col indices
  // vec2 of length #non-zero-elements -> row index
  // _values: length #non-zero-elements -> actual values

  // #cols = highest col-number + 2
  std::size_t num_cols = _entry_value_map.rbegin()->first[1] + 2;
  _inner.resize(num_cols, 0);

  // number of non-zeros are simply the number of map entries
  std::size_t num_non_zeros = _entry_value_map.size();
  _outer.resize(num_non_zeros);
  _values.resize(num_non_zeros);



  std::size_t num_non_zero = 0;
  // idea: not use conditional jumps
  //for (const auto& [k, v] : _entry_value_map) {
  //  _outer[num_non_zero] = k[0];  // add the column index
  //  _values[num_non_zero] = v;   // add the value
    // we just update the count of non-zeros at the curr. col-idx
    // note that the col-idx is automatically incremented
  //  _inner[k[1] + 1] = ++num_non_zero;
  //}

  for (std::size_t col = 0; col < num_cols; ++col){
    auto low = _entry_value_map.lower_bound({std::numeric_limits<std::size_t>::min(), col});
    auto up = _entry_value_map.upper_bound({std::numeric_limits<std::size_t>::max(), col}); 
    for (auto it = low; it != up; ++it) {
      _outer[num_non_zero] = it->first[0];  // add the column index
      _values[num_non_zero] = it->second;   // add the value
     ++num_non_zero;
    }
    _inner[col + 1] = num_non_zero; // since the next row the non-zero element start from this pos.
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

  // vec1 of length #cols + 1 -> col indices
  // vec2 of length #non-zero-elements -> row index
  // _values: length #non-zero-elements -> actual values
  std::size_t num_cols = _inner.size() - 1;
  for (std::size_t col_idx = 0; col_idx < num_cols; ++col_idx) {
    for (std::size_t row_idx = _inner[col_idx]; row_idx < _inner[col_idx + 1];
         ++row_idx) {
      // we get the col number and the value accordingly
      _entry_value_map[std::array<std::size_t, 2>{_outer[row_idx], col_idx}] =
          _values[row_idx];
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

  for (std::size_t row_idx = _inner[col]; row_idx < _inner[col + 1]; ++row_idx) {
    if (_outer[row_idx] == row) {

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

  for (std::size_t row_idx = _inner[col]; row_idx < _outer[col + 1]; ++row_idx) {
    if (_outer[row_idx] == row) {

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
  std::size_t num_rows = *max_element(_outer.begin(), _outer.end());
  res.resize(num_rows + 1, 0);
  // iterate through the colums
  for (int col_idx = 0; col_idx < _inner.size() - 1; ++col_idx) {
    for (std::size_t row_idx = _inner[col_idx]; row_idx < _inner[col_idx + 1];
         ++row_idx) {
      res[_outer[row_idx]] += (vec[col_idx] * _values[row_idx]);
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

  std::size_t num_rows = *max_element(std::begin(_outer), std::end(_outer));
  std::vector<T> sum_abs_per_col(num_rows, 0);
  for (std::size_t row_idx = 0; row_idx < _outer.size(); ++row_idx) {
    sum_abs_per_col[_outer[row_idx]] += std::abs(_values[row_idx]);
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

  T res = 0;
  for (std::size_t col_idx = 0; col_idx < _inner.size() - 1; ++col_idx) {
    // get the row, according to this col
    T norm_curr_col = 0;
    for (std::size_t row_idx = _inner[col_idx]; row_idx < _inner[col_idx + 1];
         ++row_idx) {
      norm_curr_col += std::abs(_values[row_idx]);
    }
    res = std::max(res, norm_curr_col);
  }
  return res;
}
#endif