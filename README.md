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

- ``test_file_reader``: check if the file is correctly read and if the entries of the matrix are correct
- ``test_basic_operation``: Read the matrix, test the compress and uncompress method.
- ``test_multiplication_correctness``: Execute matrix-vector, may work only on small size matrix
- ``test_norm``: compute the 3 norms of the matrix
**Only works with the matrix-market matrix.**
- ``-_benchmark_multiplication``: Execute the matrix-vector multiplication between the matrix-market
matrix and a right-hand side (- stands for three different tests)
**Only works with the matrix-market matrix.** 
- ``The provided benchmark_multiplication`` are three ranging from size and are derived from Fluid flow modeling:
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

## Implementation details
- We use a concept called `Numeric` so that you can instatiate only float/double values, this very restrictive due to the fact that our test take as input large matrices of float/double and for avoiding possible conversion erro we decide to 
restrict to this case.
- We implement compression algorithms for row/col ordering using the built-in method upper/lower bound for accessing to the value list in a specified order


