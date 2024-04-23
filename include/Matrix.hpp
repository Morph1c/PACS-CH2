#ifndef MY_SPARSE_MATRIX_HPP
#define MY_SPARSE_MATRIX_HPP
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
template <class T, StorageOrder Store = StorageOrder::row>
class Matrix {
 public:
  using EntryValuesMap = std::map<
      std::array<std::size_t, 2>, T,
      std::conditional_t<Store == StorageOrder::row, RowOrderComparator<T>,
                         ColOrderComparator<T>>>;

  /**
   * @brief Construct a new Matrix object
   *
   * @param entry_value_map Mapping with (row, col) -> value, the comparison
   * operator and thus the type of the mapping has to be based on the storage
   * order of the matrix. Simply use Matrix::EntryValuesMap type to construct
   * the mapping corectly.
   */
  Matrix(EntryValuesMap& entry_value_map)
      : _is_compressed(false),
        _entry_value_map(entry_value_map),
        _vec1(),
        _vec2(),
        _values(){};

  /**
   * @brief Construct a new Matrix object, based on a compressed format.
   * See
   * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
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
        _vec1(vec1),
        _vec2(vec2),
        _values(values){};

  /**
   * @brief Compute the norm of the matrix
   *
   * @tparam Norm options are NormOrder::frob, NormOrder::one, NormOrder::max
   * @return T norm of the matrix
   */
  template <NormOrder Norm>
  T norm() const {
#ifdef DEBUG
    std::cout << "Calling Norm with " << Norm << "\n";
#endif
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
   * NOTE: This method does for the sake of performance not check any bounds, it
   * is up to the user to call it correctly.
   *
   * @param row Index of the row.
   * @param col Index of the col.
   * @return T& Entry of the matrix.
   */
  T& operator()(std::size_t row, std::size_t col) {
#ifdef DEBUG
    std::cout << "Non-const operator() is called.\n";
#endif
    if (!_is_compressed) {
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
   * NOTE: This method does for the sake of performance not check any bounds, it
   * is up to the user to call it correctly.
   *
   * @param row Row index.
   * @param col Column index.
   * @return T Value of the matrix entry.
   */
  T operator()(std::size_t row, std::size_t col) const {
#ifdef DEBUG
    std::cout << "Const operator() is called.\n";
#endif
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
   * @brief If compressed, print the three vectors individually, if
   * decompressed, print the indices and the corresponding matrix entries.
   *
   * @param os Outputstream.
   * @param matrix Matrix.
   * @return std::ostream& Outputstream, allowing concatination.
   */
  friend std::ostream& operator<<(std::ostream& os,
                                  const Matrix<T, Store> matrix) {
    if (!matrix.is_compressed()) {
      for (const auto& [k, v] : matrix._entry_value_map) {
        os << "[" << k[0] << ", " << k[1] << "] = " << v << "\n";
      }
      return os;
    }
    os << "Compression format: " << Store << "\n\n";
    os << "vec1 = \n";
    for (const auto& el : matrix._vec1) {
      os << el << ", ";
    }
    os << "\nvec2 = \n";
    for (const auto& el : matrix._vec2) {
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
   * @brief Compute the matrix-vector-product.
   *
   * @param vec Vector x to multiply from the right-hand side.
   * @return std::vector<T> Output vector y, i.e. y = Ax.
   */
  std::vector<T> operator*(std::vector<T> vec) const {
    if (!_is_compressed) {
#ifdef DEBUG
      std::cout << "Calling matrix*vector without compression\n";
#endif
      return _matrix_vector_uncompressed(vec);
    }
    if constexpr (Store == StorageOrder::row) {
#ifdef DEBUG
      std::cout << "Calling matrix*vector for row compression\n";
#endif
      return _matrix_vector_row(vec);
    }
#ifdef DEBUG
    std::cout << "Calling matrix*vector for column compression\n";
#endif
    return _matrix_vector_col(vec);
  };

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
#ifdef DEBUG
    std::cout << "Frobenius Norm uncompressed.\n";
#endif
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
#ifdef DEBUG
    std::cout << "Frobenius Norm compressed.\n";
#endif
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
#ifdef DEBUG
    std::cout << "One-Norm uncompressed.\n";
#endif
    std::size_t num_cols = 0;
    if constexpr (Store == StorageOrder::row) {
      for (const auto& [k, v] : _entry_value_map) {
        num_cols = std::max(num_cols, k[1] + 1);
      }
    } else {
      num_cols = _entry_value_map.rbegin()->first[1] + 1;
    }
#ifdef DEBUG
    std::cout << "num_cols = " << num_cols << "\n";
#endif
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
#ifdef DEBUG
    std::cout << "One-Norm compressed.\n";
#endif
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
   * NOTE: This method does for the sake of performance not check any bounds, it
   * is up to the user to call it correctly.
   *
   * @param row Row index.
   * @param col Column index.
   * @return T
   */
  T _find_uncompressed_element(std::size_t row, std::size_t col) const {
    std::array<std::size_t, 2> to_find{row, col};
    if (auto search = _entry_value_map.find(to_find);
        search != _entry_value_map.end()) {
#ifdef DEBUG
      std::cout << "Found the element: " << search->first[0] << ", "
                << search->first[1] << "\n";
#endif
      return search->second;
    }
    return 0;
  }

  /**
   * @brief Matrix-vector product in the uncompressed state. Not the most
   * efficient way, it is recommended to operate on compressed matrices.
   *
   * @param vec Vector x to multiply on the right side.
   * @return std::vector<T> Vector y = A*x.
   */
  std::vector<T> _matrix_vector_uncompressed(std::vector<T> vec) const {
    std::size_t num_rows = 0;
    if constexpr (Store == StorageOrder::row) {
#ifdef DEBUG
      std::cout << "row\n";
#endif
      num_rows = _entry_value_map.rbegin()->first[0];
    } else {
      for (const auto& [k, v] : _entry_value_map)
        num_rows = std::max(num_rows, k[0]);
    }
#ifdef DEBUG
    std::cout << "vec.size()" << vec.size() << "\n";
    // std::cout << "res.size()" << res.size() << "\n";
    std::cout << "num_rows = " << num_rows << "\n";
    // std::cout << _entry_value_map.rbegin()->first[0] + 1 << "\n";
#endif

    std::vector<T> res(num_rows, 0);
    // res.resize(num_rows, 0);
    for (const auto& [k, v] : _entry_value_map) {
#ifdef DEBUG
      std::cout << "k[0] " << k[0] << ", k[1] " << k[1] << ": " << v << "\n";
      std::cout << "vec[k[1]] = " << vec[k[1]] << "\n";
      std::cout << "res[k[0]] = " << res[k[0]] << "\n";
#endif
      res[k[0]] += (vec[k[1]] * v);
    }
    return res;
  }

  // specialization to decide via const-expr
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
  EntryValuesMap& _entry_value_map;

  // internal representations of the values for the compressed formats
  std::vector<std::size_t> _vec1;
  std::vector<std::size_t> _vec2;
  std::vector<T> _values;
};

// ======= ROW MAJOR OPERATIONS =======
// specialization for row-major, i.e. the default case
// convert internal mapping to compressed sparse row (CSR) format
#include "Matrix_row_impl.hpp"
// ======= COL MAJOR OPERATIONS =======
#include "Matrix_col_impl.hpp"
}  // namespace algebra
#endif
