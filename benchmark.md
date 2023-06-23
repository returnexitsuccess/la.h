# Benchmarks

| Test                                     | Time Taken | NumPy     |
| ---------------------------------------- | ---------- | --------- |
| 10x10 multiply, 1,000,000 iterations     | 0.375082s  | 0.574449s |
| 100x100 multiply, 1,000 iterations       | 0.495085s  | 0.042026s |
| 1000x1000 multiply, 1 iteration          | 0.642704s  | 0.035474s |
| 1000x1000 strassen multiply, 1 iteration | 0.378002s  | 0.035474s |
| 100x100 determinant, 1,000 iterations    | 0.078891s  | 0.046177s |
| 1000x1000 determinant, 1 iteration       | 0.075856s  | 0.014904s |
| 100x100 inverse, 1,000 iterations        | 0.474230s  | 0.131105s |
| 1000x1000 inverse, 1 iteration           | 0.448225s  | 0.051880s |
