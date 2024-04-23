# Compile and run the code
```shell
make 
```

will compile the code. Afterwards you can run the code using the matrix-market file

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
Each test case is a template function and thus applied to both storage orders.

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


## Note on the compression algorithm
Contrastingly to the suggested solution of using ``lower_bound()`` and ``upper_bound()``
the here presented solution is not using these two commands. Instead, we rely on the ordering
of the mapping and only use a single loop without conditional jumps. By uncommenting on the
debugging flags in the comparison operator, one can see that for computing the ``<>_bound()``
values, internally a linear search is used. Since our map is ordered, we can hide this logic
directly in the arithmetic by calling

```cpp
for (const auto& [k, v] : _entry_value_map) {
_vec2[num_non_zero] = k[0];  // add the column index
_values[num_non_zero] = v;   // add the value
// we just update the count of non-zeros at the curr. col-idx
// note that the col-idx is automatically incremented
_vec1[k[1] + 1] = ++num_non_zero;
}
```
for the column-compression. Since the mapping is sorted based on ``k[1]`` we ensure
that the value in ``_vec1`` is always the current number of non-zero elements and
thus a correct representation.

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

# Questions to the tutor/professor
- Internally, I was working with ``constexpr`` all the time. Alternatively, one could
have worked with two template-classes each specializing one of the storage options. Which
solution is here preferred and why?

# Personal Learnings
- Class member functions can only be template specialized,
if the class itself is already specialized.
- Despite some methods only use ``O(1)`` operations, contiguous memory yields
significant speedups.