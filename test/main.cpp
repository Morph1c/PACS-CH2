#include <iostream>
#include <string>

#include "Matrix.hpp"
#include "chrono.hpp"
#include "Benchmark.hpp"

using namespace algebra;

int main(int argc, char* argv[]) {
  using type_format = double;
  std::string file_name = "./lnsp_131.mtx";
  std::string file_name_small = "./small_example.mtx";
  std::string file_name_very_small = "./very_small_example.mtx";

  std::string complex_file_name = "./lnsp_511.mtx";


  // Check if an argument was provided
  if (argc > 1) {
    // Use the provided filename
    file_name = argv[1];
  }

  Benchmark<type_format, StorageOrder::col> bench;

  // basic operations test
  //bench.test_basic_operations(complex_file_name);
  bench.test_multiplication_correctness(file_name_very_small);
  //bench.test_norm(file_name_very_small);
  // Large benchmark test with lots of runs
  //bench.test_multiplication_correctness(file_name_small);
  //bench.medium_benchmark_multiplication(1);
  //bench.large_benchmark_multiplication(1);

  return 0;
}