# Getting started
You can clone the code into your local repo as:
```shell
git clone git@github.com:Morph1c/PACS-CH2.git 
```
Then for compile the code, simply run:
```shell
make 
```
Afterwards you can run the code using the matrix-market file

```sh
./main
```

or using a smaller example.
```sh
./main small_example.mtx
```

To see debugging information, call
```sh
make debug
```

## Provided test cases
It's provided a benchmark template class where you can test:

- ``test_file_reader``: check whether the file is correctly read and the mapping is correctly extracted.
- ``test_basic_operation``: Read the matrix, test setter and getter methods,
test compression and uncompress. Both are checked with the python utility script.
- ``test_multiplication_correctness``: Reads a small test matrix and executes the matrix-vector 
multiplication. Result is checked for correctness using the helper python script. **Only works for the
small test case.**
- ``test_norm``: Compute all of the three provided norms.
- ``benchmark_multiplication``: Stop the time for the compression and the multiplication for both status.
**Only works with the matrix-market matrix.**
- ``large_benchmark_multiplication``: Execute the matrix-vector multiplication between the matrix-market
matrix and a right-hand side for ``num_runs=1000`` times, compute the average afterwards.
**Only works with the matrix-market matrix.** The results are the following:

| | Row major | Col major |
| ---  | --------- | --- |
| Uncompressed  | 6.333 \mu s | 8.083 \mu s  |
| Compressed  | **1.666 \mu s** | 2.792 \mu s  |
| | | |
| Speedup  | **3.815** | 2.895 |


## Design Decisions
To test for maximum speed, I did not implement any bound checks in the getter and
setter methods.

## Possible Extensions/Improvements
- Implement safe getter and setter methods with a bound check.
- Save the number of rows and columns as a class attribute and update them in the
modifying methods. This would simplify some function implementations and ensure
even further speed-ups.
- Make the large benchmark test more generic.
- Extension: work with ``mutable``.


