#include <iostream>
#include <string>

#include "Matrix.hpp"
#include "chrono.hpp"
#include "test_cases.hpp"

using namespace algebra;

int main(int argc, char* argv[]) {
  using type_format = double;
  std::string file_name = "./lnsp_131.mtx";
  std::string file_name_small = "./small_example.mtx";

  // Check if an argument was provided
  if (argc > 1) {
    // Use the provided filename
    file_name = argv[1];
  }
  // Test if the filereader works
  // test_file_reader<type_format, StorageOrder::row>(file_name);
  // test_file_reader<type_format, StorageOrder::col>(file_name);

  // Test if constructor, compression and call operator work
  // test_basic_operations<type_format, StorageOrder::row>(file_name_small);
  // test_basic_operations<type_format, StorageOrder::col>(file_name_small);

  // test the norm computation
  // test_norm<type_format, StorageOrder::row>(file_name);
  // test_norm<type_format, StorageOrder::col>(file_name);

  // Test if multiplication works
  // test_multiplication_correctness<type_format, StorageOrder::row>(file_name);
  // test_multiplication_correctness<type_format, StorageOrder::col>(file_name);

  // Benchmark test for the matrix-vector multiplication
  // benchmark_multiplication<type_format, StorageOrder::row>();
  // benchmark_multiplication<type_format, StorageOrder::col>();

  // Large benchmark test with lots of runs
  large_benchmark_multiplication<type_format, StorageOrder::row>(1);
  large_benchmark_multiplication<type_format, StorageOrder::col>(1);
  return 0;
}