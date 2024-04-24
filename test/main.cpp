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

  // Check if an argument was provided
  if (argc > 1) {
    // Use the provided filename
    file_name = argv[1];
  }
  // basic operations test
  // test_basic_operations<type_format, StorageOrder::row>(file_name_small);
  //test_multiplication_correctness<type_format, StorageOrder::row>(file_name);

  // Large benchmark test with lots of runs
  Benchmark<type_format, StorageOrder::row> bench;
  bench.test_multiplication_correctness(file_name_small);
  bench.large_benchmark_multiplication(1);
  bench.large_benchmark_multiplication(1);

  return 0;
}