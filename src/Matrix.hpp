#ifndef MATRIX_SPARSE_HPP
#define MATRIX_SPARSE_HPP
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include "Utilities.hpp"

namespace algebra {

  
/**
 * @brief Class representing a sparse matrix, which can be stored in row or
 * column major format. The matrix can be compressed into a compressed sparse
 * row (CSR) or compressed sparse column (CSC) format.
 *
 * @tparam T Type of the entries of the matrix.
 * @tparam Store Storage order of the matrix, either row or column major.
 */
template <class T, StorageOrder Store = StorageOrder::row>
class Matrix {

public:
  // define the type of the matrix to be used as data structure for the matrix rappresentation
  using matrix_type = std::map<
      std::array<std::size_t, 2>, T,
      std::conditional_t<Store == StorageOrder::row, RowOrderComparator<T>,
                         ColOrderComparator<T>>>;
  
private:
  template <NormOrder Norm>
  /**
   * @brief Internal method to wrapt the norm computation in the uncompressed
   * case.
   *
   * @return T Norm of the matrix.
   */
  T _norm_uncompressed() const {
    // ONE NORM
    if constexpr (Norm == NormOrder::one) {
      return _one_norm_uncompressed();
    }
    return _max_norm_uncompressed();
  };

  /**
   * @brief Frobenius norm, i.e. \sqrt{sum_{i, j} \abs{a_{ij}^ 2}}
   *
   * @return T Norm.
   */
  T _frob_norm_uncompressed() const {
    T res = 0;
    for (const auto& [k, v] : _entry_value_map) res += std::norm(v);
    return std::sqrt(res);
  };

  /**
   * @brief Frobenius norm compressed version, i.e. \sqrt{sum_{i, j}
   * \abs{a_{ij}^ 2}}.
   *
   * @return T Norm.
   */
  T _frob_norm_compressed() const {

    T res = 0;
    for (const auto& val : _values) res += std::norm(val);
    return std::sqrt(res);
  }

  /**
   * @brief One norm, i.e. max(sum(abs(x), axis=0)).
   *
   * @return T Norm of the matrix.
   */
  T _one_norm_uncompressed() const {

    std::size_t num_cols = 0;
    if constexpr (Store == StorageOrder::row) {
      for (const auto& [k, v] : _entry_value_map) {
        num_cols = std::max(num_cols, k[1] + 1);
      }
    } else {
      num_cols = _entry_value_map.rbegin()->first[1] + 1;
    }

    std::vector<T> sum_abs_per_col(num_cols, 0.0);
    for (const auto& [k, v] : _entry_value_map) {
      sum_abs_per_col[k[1]] += std::abs(v);
    }
    return *max_element(std::begin(sum_abs_per_col), std::end(sum_abs_per_col));
  };

  /**
   * @brief Max norm, i.e. max(sum(abs(x), axis=1)).
   *
   * @return T Norm of the matrix.
   */
  T _max_norm_uncompressed() const {

    std::size_t num_rows = 0;
    if constexpr (Store == StorageOrder::row) {
      num_rows = _entry_value_map.rbegin()->first[0] + 1;
    } else {
      for (const auto& [k, v] : _entry_value_map) {
        num_rows = std::max(num_rows, k[0] + 1);
      }
    }
    std::vector<T> sum_abs_per_row(num_rows, 0.0);
    for (const auto& [k, v] : _entry_value_map) {
      sum_abs_per_row[k[0]] += std::abs(v);
    }
    return *max_element(std::begin(sum_abs_per_row), std::end(sum_abs_per_row));
  };

  /**
   * @brief Find an element in the uncompressed state, used for the
   * setter/getter methods. Returns the value, if the element is found, else 0.
   *
   * @param row Row index.
   * @param col Column index.
   * @return T
   */
  T _find_uncompressed_element(std::size_t row, std::size_t col) const {
    std::array<std::size_t, 2> to_find{row, col};
    if (auto search = _entry_value_map.find(to_find);
        search != _entry_value_map.end()) {

      return search->second;
    }
    return 0;
  }

  /**
   * @brief Matrix-vector product in the uncompressed state. This type of multiplication
   *  is not the most efficent, maybe better first to compress and then use the compressed multiplication
   *
   * @param vec Vector x to multiply on the right side.
   * @return std::vector<T> Vector y = A*x.
   */
  std::vector<T> _uncompressed_mult(std::vector<T> vect) const {
    std::size_t num_rows, num_cols = 0;
    if constexpr (Store == StorageOrder::row) {

      num_rows = _entry_value_map.rbegin()->first[0];
      num_cols = _entry_value_map.rbegin()->first[1];
    } else {
      for (const auto& [k, v] : _entry_value_map){
        num_rows = std::max(num_rows, k[0]);
        num_cols = std::max(num_cols, k[1]);
      }

    }
    
    //std::cout << "num cols = " << num_cols << "\n";
    std::vector<T> res(num_rows + 1, 0); // should be +1 since the rows starts from 0
  
    for (const auto& [k, v] : _entry_value_map) {

      res[k[0]] += (vect[k[1]] * v);
    }
    return res;
  }

  // specialization to decide via const-expr
  // define the specialization inside different files
  void _compress_row();
  void _uncompress_row();
  const T _find_compressed_element_row(std::size_t row, std::size_t col) const;
  T& _find_compressed_element_row(std::size_t row, std::size_t col);
  std::vector<T> _matrix_vector_row(std::vector<T>) const;
  T _one_norm_compressed_row() const;
  T _max_norm_compressed_row() const;

  void _compress_col();
  void _uncompress_col();
  const T _find_compressed_element_col(std::size_t row, std::size_t col) const;
  T& _find_compressed_element_col(std::size_t row, std::size_t col);
  std::vector<T> _matrix_vector_col(std::vector<T>) const;
  T _one_norm_compressed_col() const;
  T _max_norm_compressed_col() const;
  

  // class attributes
  bool _is_compressed;
  matrix_type& _entry_value_map;

  // internal representations of the values for the compressed formats
  std::vector<std::size_t> _inner;
  std::vector<std::size_t> _outer;
  std::vector<T> _values;

 public:
  
  /**
   * @brief Construct a new Matrix object
   *
   * @param entry_value_map Mapping with (row, col) -> value, the comparison
   * operator and thus the type of the mapping has to be based on the storage
   * order of the matrix. Simply use Matrix::matrix_type type to construct
   * the mapping corectly.
   */
  Matrix(matrix_type& value_map)
      : _is_compressed(false),
        _entry_value_map(value_map),
        _inner(),
        _outer(),
        _values(){};

  /**
   * @brief Construct a new Matrix object, based on a compressed format.
   * See
   *
   * @param vec1 If row-compression, vec1 contains the row-indices, i.e.
   * ROW_INDEX in the wikipedia article in the CSR-format. If col-compression
   * vec1 contains the column indices.
   * @param vec2 If row-compression, vec2 contains the col-indices, i.e.
   * COL_INDEX in the wiki article. If col-compression, vec2 contains the
   * row-indices.
   * @param values Always the same, i.e. the V vector in the wiki article.
   */
  Matrix(std::vector<std::size_t> vec1, std::vector<std::size_t> vec2,
         std::vector<T> values)
      : _is_compressed(true),
        _entry_value_map(),
        _inner(vec1),
        _outer(vec2),
        _values(values){};

  /**
   * @brief Compute the norm of the matrix
   *
   * @tparam Norm options are NormOrder::frob, NormOrder::one, NormOrder::max
   * @return T norm of the matrix
   */
  template <NormOrder Norm>
  T norm() const {

    // FROB norm is the easiest case
    if constexpr (Norm == NormOrder::frob) {
      if (!_is_compressed)
        return _frob_norm_uncompressed();
      else
        return _frob_norm_compressed();
    }

    if (!_is_compressed) return _norm_uncompressed<Norm>();
    // col compression
    if constexpr (Store == StorageOrder::col) {
      if constexpr (Norm == NormOrder::one) {
        return _one_norm_compressed_col();
      }
      return _max_norm_compressed_col();
    }
    // row compression
    if constexpr (Norm == NormOrder::one) {
      return _one_norm_compressed_row();
    }
    return _max_norm_compressed_row();
  };

  /**
   * @brief Non-const getter and setter, in the compressed case only non-zero
   * elements can be changed, else an exception is thrown.
   *  not check if the row, col are out-of-bound, up to the user t check it correctly
   *
   * @param row Index of the row.
   * @param col Index of the col.
   * @return T& Entry of the matrix.
   */
  T& operator()(std::size_t row, std::size_t col) {

    if (!_is_compressed) { // so is the dynamic storage case 
      std::array<std::size_t, 2> find{row, col};
      // either add or override, both is fine
      return _entry_value_map[find];
    }
    if constexpr (Store == StorageOrder::row) {
      // only existing values can be added
      return this->_find_compressed_element_row(row, col);
    } else {
      // only existing values can be added
      return this->_find_compressed_element_col(row, col);
    }
  };

  /**
   * @brief Const getter, returning the value of a matrix entry.
   * i am not checking if the row and col exist effectively, it is up to the user
   *
   * @param row Row index.
   * @param col Column index.
   * @return T Value of the matrix entry.
   */
  T operator()(std::size_t row, std::size_t col) const {

    if (!_is_compressed) {
      return this->_find_uncompressed_element(row, col);
    }
    // compressed case, check for row/col format
    if constexpr (Store == StorageOrder::row) {
      return this->_find_compressed_element_row(row, col);
    }
    return this->_find_compressed_element_col(row, col);
  };

  /**
   * @brief Compute the matrix-vector-product.
   *
   * @param vec Vector x to multiply from the right-hand side.
   * @return std::vector<T> Output vector y, i.e. y = Ax.
   */
  friend std::vector<T> operator*(const Matrix<T, Store> matrix, std::vector<T> vec){
    if (!matrix._is_compressed) {
      return matrix._uncompressed_mult(vec);
    }
    if constexpr (Store == StorageOrder::row) {

      return matrix._matrix_vector_row(vec);
    }

    return matrix._matrix_vector_col(vec);
  };

  
  /**
   * @brief Overload the output operator to print the matrix.
   *
   * @param os Output stream.
   * @param matrix Matrix to print.
   * @return std::ostream& Output stream.
   */

  friend std::ostream& operator<<(std::ostream& os,
                                  const Matrix<T, Store> matrix) {
    if (!matrix.is_compressed()) {
      for (const auto& [k, v] : matrix._entry_value_map) {
        os << "[" << k[0] << ", " << k[1] << "] = " << v << "\n";
      }
      return os;
    }
    os << "The compression format: " << Store << "\n\n";
    os << "inner = \n";
    for (const auto& el : matrix._inner) {
      os << el << ", ";
    }
    os << "\nouter = \n";
    for (const auto& el : matrix._outer) {
      os << el << ", ";
    }
    os << "\nvalues = \n";
    for (const auto& el : matrix._values) {
      os << el << ", ";
    }
    os << "\n";
    return os;
  }

  /**
   * @brief Compress the matrix into row/column sparse format, i.e switch
   * from the internal mapping to a three-vector representation.
   */
  void compress() {
    if constexpr (Store == StorageOrder::row) {
      _compress_row();
    } else {
      _compress_col();
    };
  }
  /**
   * @brief Decompress the matrix from the three-vector format back to the
   * internal mapping format.
   */
  void uncompress() {
    if constexpr (Store == StorageOrder::row) {
      _uncompress_row();
    } else {
      _uncompress_col();
    };
  }

  bool is_compressed() const { return _is_compressed; };
};



// ======= ROW MAJOR OPERATIONS =======
// specialization for row-major, i.e. the default case
// convert internal mapping to compressed sparse row (CSR) format
#include "Matrix_row_impl.hpp"
// ======= COL MAJOR OPERATIONS =======
#include "Matrix_col_impl.hpp"

}  // namespace algebra

#endif
