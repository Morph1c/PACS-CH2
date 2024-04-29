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

-The provided test are 3 ranging from size and are derived from Fluid flow modeling:
- ``lns``: simpler case with a 10x10 matrix
- ``lns_131``: medium complexity case with a 131x131 matrix
- ``lns_511``: high complexity with a 511x511 matrix

## Medium case
| | Row major | Col major |
| ---  | --------- | --- |
| Uncompressed  | 14.831 \mu s | 8.083 \mu s  |
| Compressed  | **5.104 \mu s** | 2.792 \mu s  |
| | | |

## Large case
| | Row major | Col major |
| ---  | --------- | --- |
| Uncompressed  | 42.601 \mu s | 11428.2 \mu s  |
| Compressed  | **6.218 \mu s** | 13.799 \mu s  |
| | | |




