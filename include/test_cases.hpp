
#ifndef TEST_CASES_MATRIX_HPP
#define TEST_CASES_MATRIX_HPP
#include <iostream>
#include <string>

#include "Matrix.hpp"
#include "ReadMatrix.hpp"
#include "Utilities.hpp"
#include "chrono.hpp"

using namespace algebra;

// small utility function
void _print_test_case(StorageOrder Order) {
  std::cout << "\n\n===========================\n";
  std::cout << "Test case for " << Order << "\n";
}

// Test: read a matrix as a matrix-market file and print it.
template <class T, StorageOrder Store>
void test_file_reader(const std::string& file_name) {
  try {
    // Example usage to read matrix in Matrix Market format
    auto matrix_mapping = read_matrix<T, Store>(file_name);
    std::cout << "reading row worked\n";
    for (const auto& [k, v] : matrix_mapping) {
      std::cout << "[" << k[0] << ", " << k[1] << "] = " << v << std::endl;
    }

    auto colMatrix = read_matrix<T, Store>(file_name);
    std::cout << "reading col worked\n";
    for (const auto& [k, v] : colMatrix) {
      std::cout << "[" << k[0] << ", " << k[1] << "] = " << v << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}

// Test: read the matrix, get and set some values, execute compression
template <class T, StorageOrder Order>
void test_basic_operations(const std::string& file_name) {
  _print_test_case(Order);

  auto matrix_mapping = read_matrix<T, Order>(file_name);
  auto matrix = Matrix<T, Order>(matrix_mapping);

  std::cout << "Printing the matrix\n";
  std::cout << matrix;
  // should yield 0 and 4.0
  std::cout << "Calling the setter/getter methods: \n";
  T val0 = matrix(0, 3);  // not present

  T val4 = matrix(3, 3);  // present
  std::cout << "val0 = " << val0 << ", val4 = " << val4 << "\n";

  // Test compression + call operator
  matrix.compress();
  std::cout << "COMPRESSION WORKED\n";
  std::cout << matrix;

  std::cout << "Calling the setter/getter methods: \n";
  T val0_comp = matrix(0, 3);  // not present
  // throws an exception -> works correctly
  // T val11 = 11.0;
  // matrix(0, 4) = val11;

  T val4_comp = matrix(3, 3);  // present
  std::cout << "val0 = " << val0_comp << ", val4 = " << val4_comp << "\n";

  // Test uncompress + call operator
  matrix.uncompress();
  std::cout << "UNCOMPRESSION WORKED\n";
  std::cout << matrix;
}

// Test: matrix-vector multiplication for the small case
template <class T, StorageOrder Order>
void test_multiplication_correctness(const std::string& file_name) {
  _print_test_case(Order);
  auto matrix_mapping = read_matrix<T, Order>(file_name);
  auto matrix = Matrix<T, Order>(matrix_mapping);
  std::vector<T> to_multiply = {42, -17, 100, 0, 73, -5, 21, 8, -33, 55};

  auto res_uncomp = matrix * to_multiply;
  std::cout << "\nUncompressed Multiplication worked\n";
  for (auto i : res_uncomp) std::cout << i << ' ';
  matrix.compress();
  std::cout << matrix << "\n";
  auto res_comp = matrix * to_multiply;
  std::cout << "\nCompressed Multiplication worked\n";

  for (auto i : res_comp) std::cout << i << ' ';
}

// Test: compute all norms
template <class T, StorageOrder Order>
void test_norm(const std::string& file_name) {
  _print_test_case(Order);
  auto matrix_mapping = read_matrix<T, Order>(file_name);
  auto matrix = Matrix<T, Order>(matrix_mapping);

  // uncompressed
  auto frob_uncomp = matrix.template norm<NormOrder::frob>();
  auto one_uncomp = matrix.template norm<NormOrder::one>();
  auto max_uncomp = matrix.template norm<NormOrder::max>();

  std::cout << "Uncompressed Norm: "
            << " frob = " << frob_uncomp << ", one = " << one_uncomp
            << ", max = " << max_uncomp << "\n";

  // compressed version
  matrix.compress();
  auto frob_comp = matrix.template norm<NormOrder::frob>();
  auto one_comp = matrix.template norm<NormOrder::one>();
  auto max_comp = matrix.template norm<NormOrder::max>();
  std::cout << "Compressed Norm: "
            << " frob = " << frob_comp << ", one = " << one_comp
            << ", max = " << max_comp << "\n";
}

// Test: Stop the time for a single multiplication, as well as for the
// compression operations
template <class T, StorageOrder Order>
void benchmark_multiplication() {
  _print_test_case(Order);
  std::string file_name = "./lnsp_131.mtx";
  Timings::Chrono timer;
  // std::vector<T> to_multiply = _generate_random_vector<T>(131);

  std::cout << Order << "-MAJOR UNCOMPRESSED Multiplication took:\n";
  auto matrix_mapping = read_matrix<T, Order>(file_name);
  auto matrix = Matrix<T, Order>(matrix_mapping);

  std::vector<T> to_multiply = _generate_random_vector<T>(131);
  timer.start();
  auto res_raw = matrix * to_multiply;
  timer.stop();
  std::cout << timer;

  std::cout << Order << "-MAJOR COMPRESSED Compression took:\n";
  timer.start();
  matrix.compress();
  timer.stop();
  std::cout << timer;

  std::cout << Order << "-MAJOR COMPRESSED Multiplication took:\n";
  timer.start();
  auto res_compressed = matrix * to_multiply;
  timer.stop();
  std::cout << timer;
}

// Test: large scale benchmark test to average the computation time.
template <class T, StorageOrder Order>
void large_benchmark_multiplication(std::size_t num_runs) {
  std::string file_name = "./lnsp_131.mtx";
  Timings::Chrono timer;
  auto matrix_mapping_raw = read_matrix<T, Order>(file_name);
  auto matrix_mapping_comp = read_matrix<T, Order>(file_name);

  auto matrix_raw = Matrix<T, Order>(matrix_mapping_raw);
  auto matrix_compressed = Matrix<T, Order>(matrix_mapping_comp);
  matrix_compressed.compress();

  // save the times
  double total_time_raw = 0.0;
  double total_time_compressed = 0.0;

  for (std::size_t i = 0; i < num_runs; ++i) {
    std::vector<T> to_multiply = _generate_random_vector<T>(131);

    timer.start();
    auto res_raw = matrix_raw * to_multiply;
    timer.stop();
    total_time_raw += timer.wallTime();

    timer.start();
    auto res_compressed = matrix_compressed * to_multiply;
    timer.stop();
    total_time_compressed += timer.wallTime();
  }

  // Calculate average times
  double avg_time_raw = total_time_raw / num_runs;
  double avg_time_compressed = total_time_compressed / num_runs;

  // Report average times
  std::cout << "Large Benchmark Test for " << Order << "\n";
  std::cout << "Average time for UNCOMPRESSED Multiplication: " << avg_time_raw
            << " micro-seconds\n";
  std::cout << "Average time for COMPRESSED Multiplication: "
            << avg_time_compressed << " micro-seconds\n";
}
#endif