# Benchmarks

| Test                                     | Time Taken | NumPy     |
| ---------------------------------------- | ---------- | --------- |
| 10x10 multiply, 1,000,000 iterations     | 0.403451s  | 0.574449s |
| 100x100 multiply, 1,000 iterations       | 0.527973s  | 0.042026s |
| 1000x1000 multiply, 1 iteration          | 0.664094s  | 0.035474s |
| 1000x1000 strassen multiply, 1 iteration | 0.392258s  | 0.035474s |
| 100x100 determinant, 1,000 iterations    | 0.137948s  | 0.046177s |
| 1000x1000 determinant, 1 iteration       | 0.117739s  | 0.014904s |
| 100x100 inverse, 1,000 iterations        | 0.313774s  | 0.131105s |
| 1000x1000 inverse, 1 iteration           | 0.322295s  | 0.051880s |
