# la.h - Linear Algebra Library in C

This library is far from finished and I am writing it for my own purposes/interests.
There are definitely far better libraries to use, but feel free to use this if you'd like.


## Current/Planned Features

- [x] Matrix Addition
- [x] Scalar Multiplication
- [x] Matrix Multiplication
- [x] Determinant (n! & n^3)
- [x] Row Reduction
- [x] Matrix Inverse
- [x] Transpose
- [x] QR Decomposition
- [ ] Solving Linear Systems
- [ ] Eigenvalues
- [ ] Eigenvectors
- [ ] Diagonalization
- [ ] Alternate types (currently only double)


## Output of test.c

```
a =
0.680375 -0.211234 0.566198 0.596880 0.823295
-0.604897 -0.329554 0.536459 -0.444451 0.107940
-0.045206 0.257742 -0.270431 0.026802 0.904459
0.832390 0.271423 0.434594 -0.716795 0.213938
-0.967399 -0.514226 -0.725537 0.608354 -0.686642

b =
-0.198111 -0.740419 -0.782382 0.997849 -0.563486
0.025865 0.678224 0.225280 -0.407937 0.275105
0.048574 -0.012834 0.945550 -0.414966 0.542715
0.053490 0.539828 -0.199543 0.783059 -0.433371
-0.295083 0.615449 0.838053 -0.860489 0.898654

a + b =
0.482264 -0.951653 -0.216184 1.594729 0.259809
-0.579032 0.348670 0.761739 -0.852387 0.383044
0.003368 0.244908 0.675119 -0.388165 1.447175
0.885880 0.811251 0.235051 0.066264 -0.219433
-1.262482 0.101223 0.112516 -0.252136 0.212012

a * b =
-0.323764 0.174615 0.526330 0.289085 0.346978
0.081746 0.043985 1.085414 -1.132683 0.830946
-0.262971 0.782866 0.590363 -0.795322 0.750793
-0.238246 -0.693087 0.143151 -0.205850 0.344384
0.378267 0.282644 -0.741835 0.612752 -0.870806

|a| = 0.384294

|b| = 0.021203

Echelon form of a:
0.680375 -0.211234 0.566198 0.596880 0.823295
0.000000 -0.517355 1.039846 0.086214 0.839902
-0.000000 0.000000 0.257022 0.107072 1.354808
0.000000 0.000000 0.000000 -1.694865 -4.186203
0.000000 0.000000 0.000000 0.000000 2.506233

Inverse of a:
0.475320 -0.514907 -0.046170 1.463812 0.884239
-1.068258 -0.623726 -0.459988 -2.997072 -2.918616
0.034383 0.448218 -0.962537 -1.721875 -1.692675
0.266967 -0.472822 -0.524217 -1.735642 -0.985516
0.330546 0.300033 0.962145 0.463828 0.399005

Transpose of a:
0.680375 -0.604897 -0.045206 0.832390 -0.967399
-0.211234 -0.329554 0.257742 0.271423 -0.514226
0.566198 0.536459 -0.270431 0.434594 -0.725537
0.596880 -0.444451 0.026802 -0.716795 0.608354
0.823295 0.107940 0.904459 0.213938 -0.686642

QR Decomposition:
Q =
-0.433828 -0.750602 0.063934 0.414041 -0.269948
0.385701 -0.249605 -0.850289 -0.076826 -0.245030
0.028825 0.481808 0.063601 0.381533 -0.785760
-0.530757 0.020782 -0.069512 -0.754683 -0.378797
0.616843 -0.376461 0.513853 -0.327955 -0.325858

R =
-1.568308 -0.489297 -0.724721 0.326107 -0.826564
-0.000000 0.564220 -0.407020 -0.568087 0.053806
0.000000 0.000000 -0.840174 0.780207 -0.349323
-0.000000 -0.000000 0.000000 0.632944 0.741398
-0.000000 0.000000 -0.000000 0.000000 -0.816675
```