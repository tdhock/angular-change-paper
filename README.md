# Figures for paper

* time complexity of angular changepoint
  packages. [PNG](figure-compare-time-sim-nopred-nomem.png),
  [timings code](figure-compare-time-sim-data.R), 
  [figure code](figure-compare-time-sim.R).
* time complexity of FPOP with different loss
  functions. [PNG](figure-geodesicFPOP-robseg-simple.png),
  [timings code](figure-geodesicFPOP-robseg-data.R),
  [figure code](figure-geodesicFPOP-robseg.R).
  
# 20 Feb 2025

Some example loss functions that are created using only the addition
operation, using the L1 geodesic loss (equivalent to FPOP with
infinite penalty). 

[figure-pfpop-add-operations.R](figure-pfpop-add-operations.R) makes

Below we see the best case time complexity, when all data are in a
single tight cluster, there are only two clusters (independent of
number of data points), so each update is a constant time operation,
log-linear time overall in the number of data. Log factor comes from
inserting a new breakpoint into the C++ STL map, which is a log time
operation. Each "cluster" of positive/negative changes in slope is
visually represented by a grey and violet rectangle, and each
iteration of the algorithm is linear in the number of such
clusters. (positive changes for locally convex regions, negative
changes for locally concave)
  
![screenshot](figure-pfpop-add-operations-best.png)

Below shows a data set that explains worst case time complexity:
number of clusters is linear in the number of data points, so
quadratic time overall.

![screenshot](figure-pfpop-add-operations-worst.png)

Below shows a typical data set, with three tight clusters of data
values, resulting in six clusters of positive/negative breakpoints in
loss slope. We see how the global min starts on the left, moves to the
middle, and then the right.

![screenshot](figure-pfpop-add-operations-typical.png)

Below we see an example of clusters of breakpoints which are not part
of a local optimum, at the last iteration. Locally concave or convex
regions can be found in areas of the function that are strictly
increasing or decreasing.

![screenshot](figure-pfpop-add-operations-critx.png)

# 3 Feb 2025

[figure-pfpop-worst-case.R](figure-pfpop-worst-case.R) makes

![screenshot](figure-pfpop-worst-case-typical.png)

The figure above shows what would be probably a "typical" run of the
FPOP algorithm with infinite penalty, using the map data structure, and a constant number of
pointer moves per iteration. We see a single minimum between 0 and 90,
which the constant number of pointer moves is capable of tracking.

![screenshot](figure-pfpop-worst-case.png)

The figure above shows a pathological data set for the FPOP algorithm
with infinite penalty, using the map data structure, and a constant
number of pointer moves per iteration. We see that the min/max is not
correct for iterations 2 and 4, because the constant number of pointer
moves keeps the pointer in a local min/max, but there is a global
min/max which would require more pointer moves.

# 31 Jan 2025

[figure-pfpop-atime.R](figure-pfpop-atime.R) makes

![screenshot](figure-pfpop-atime.png)

The figure above shows asymptotic properties of FPOP using infinite
penalty and an alternative data structure (map), with O(1) pointer
moves per iteration, versus the classic data structure (list), with
O(I) operations per iteration, where I is the number of intervals used
to represent the cost function. Both have O(N) space but with O(N)
total moves, we see that the time complexity of map is O(N log N),
versus O(N^2) for list. (in the best case, worst case both are quadratic)

![screenshot](figure-pfpop-atime-compare.png)

The figure above shows that

* cost: list has same cost values as map in all cases.
* kilobytes: map takes a constant factor more memory than list.
* moves: linear in best case, quadratic in worst case.
* pointers: constant in best case, linear in worst case.
* seconds: map best case is log-linear, others are quadratic.

# 15 Jan 2025

[figure-2d-hmm-sim-noise-data.R](figure-2d-hmm-sim-noise-many-data.R) and
[figure-2d-hmm-sim-noise.R](figure-2d-hmm-sim-noise-many.R) 
makes

![screenshot](figure-2d-hmm-sim-noise-many.png)

The figure shows the advantage of the proposed APART algo, with respect to two baselines:

* SinCosL2BinSeg is the classic binary segmentation heuristic on the
  4d data matrix computed by doing the sin/cos basis expansion on both
  features. The parameter of the binary segmentation algorithm is the
  number of changes, which we selected to agree with the true number
  of changes that were used to generate the data (10, vertical black
  lines). The figure shows that even with little noise, binary
  segmentation results in false negatives (two real change-points in
  black which were not detected), and therefore also false positives
  (change-points detected where there should be none).
* VonMisesHMM is a hidden Markov model, with the number of states
  selected to agree with the true number of states that were used to
  generate the data (3). While the HMM gives reasonable results when
  there is low noise (kappa=5), it gives many false positive
  change-point detections when there is high noise (kappa=1).
* APART is our proposed Approximate Partitioning algorithm, using
  K-means as a pre-processing technique (with K=3), and with the
  penalty chosen in order to recover the correct number of
  change-points. It is clear that APART recovers a model that is
  consistent with the true change-points (vertical black lines), for
  both low and high noise settings.
  
# 3 oct 2024

[msmbuilder.hmm.VonMisesHMM](http://msmbuilder.org/3.4.0/_hmm/msmbuilder.hmm.VonMisesHMM.html)

```shell
conda install -c omnia msmbuilder
```

[source is C++](https://github.com/msmbuilder/msmbuilder/blob/master/msmbuilder/hmm/src/VonMisesHMMFitter.cpp)

[figure_msmbuilder_data.py](figure_msmbuilder_data.py) and 
[figure-2d-hmm-sim-noise-data.R](figure-2d-hmm-sim-noise-data.R) and
[figure-2d-hmm-sim-noise.R](figure-2d-hmm-sim-noise.R) 
makes

![screenshot](figure-2d-hmm-sim-noise.png)

# 19 Sept 2024

[figure-2d-hmm-sim.R](figure-2d-hmm-sim.R) makes

![sim figure](figure-2d-hmm-sim.png)
  
# 3 July 2024

[interesting_signals.py](interesting_signals.py) has indices of real
data with segments near the boundary.

[figure-2d-hmm-real-interactive.R](figure-2d-hmm-real-interactive.R)
makes 
[interactive HMM data viz](https://tdhock.github.io/2024-07-04-HMM-angular-data/)

![screenshot](figure-2d-hmm-real-interactive.png)
  
# 28 June 2024

[figure-2d-hmm.R](figure-2d-hmm.R) makes

![pointer moves figure](figure-2d-hmm.png) 

# 18 Apr 2024

[figure-pointer-moves.R](figure-pointer-moves.R) makes

![pointer moves figure](figure-pointer-moves.png) 

# 16 Apr 2024

[figure-geodesicFPOP-robseg.R](figure-geodesicFPOP-robseg.R) makes

![asymptotic two penalties several losses](figure-geodesicFPOP-robseg-simple.png)

![asymptotic two penalties several losses](figure-geodesicFPOP-robseg.png)

[figure-geodesic-penalties.R](figure-geodesic-penalties.R) makes

![asymptotic four penalties for geodesic loss](figure-geodesic-penalties.png)

# 10 Apr 2024

[figure-approx-algo-1d.R](figure-approx-algo-1d.R) makes

![approx algo penalty 1](figure-approx-algo-1d-1.png)

![approx algo penalty 100](figure-approx-algo-1d-100.png)

![approx algo penalty 10000](figure-approx-algo-1d-10000.png)

![approx algo three penalty values](figure-approx-algo-1d.png)

# 2 Apr 2024

[figure-compare-time-data.R](figure-compare-time-data.R) makes

![real data timings](figure-compare-time.png)

[figure-compare-time-sim-data.R](figure-compare-time-sim-data.R) makes

![simulated data timings](figure-compare-time-sim-nopred-nomem.png)

![simulated data timings](figure-compare-time-sim.png)

![simulated data timings with pred](figure-compare-time-sim-pred.png)

https://cloud.r-project.org/web/packages/moveHMM

https://cloud.r-project.org/web/packages/circular

# 3 Mar 2024

Modified robseg [figure-robseg.R](figure-robseg.R)
https://github.com/tdhock/robust-fpop/tree/interval-count shows that
number of intervals is indeed linear if penalty is large.

```r
      lambda lthreshold lslope intervals  path
       <num>      <num>  <num>     <int> <int>
   1:    100          1      0         3    -1
   2:    100          1      0         5    -1
   3:    100          1      0         7    -1
   4:    100          1      0         9    -1
   5:    100          1      0        11    -1
  ---                                         
 996:    100          1      0       724    -1
 997:    100          1      0       722    -1
 998:    100          1      0       724    -1
 999:    100          1      0       723    -1
1000:    100          1      0       726    -1
      lambda lthreshold lslope intervals  path
       <num>      <num>  <num>     <int> <int>
   1:    0.1          1      0         3    -1
   2:    0.1          1      0         5     1
   3:    0.1          1      0         5     2
   4:    0.1          1      0         5     3
   5:    0.1          1      0         5     4
  ---                                         
 996:    0.1          1      0         5   994
 997:    0.1          1      0         6   996
 998:    0.1          1      0         6   995
 999:    0.1          1      0         7   998
1000:    0.1          1      0         5   999
## Above biweight, below L1.
      lambda lthreshold lslope intervals  path
       <num>      <num>  <num>     <int> <int>
   1:    100          0     -1         3    -1
   2:    100          0     -1         5    -1
   3:    100          0     -1         7    -1
   4:    100          0     -1         9    -1
   5:    100          0     -1        11    -1
  ---                                         
 996:    100          0     -1       787    -1
 997:    100          0     -1       786    -1
 998:    100          0     -1       788    -1
 999:    100          0     -1       787    -1
1000:    100          0     -1       788    -1
      lambda lthreshold lslope intervals  path
       <num>      <num>  <num>     <int> <int>
   1:    0.1          0     -1         3    -1
   2:    0.1          0     -1         7     1
   3:    0.1          0     -1         7     2
   4:    0.1          0     -1         7     3
   5:    0.1          0     -1         7     4
  ---                                         
 996:    0.1          0     -1         7   995
 997:    0.1          0     -1         7   996
 998:    0.1          0     -1         7   997
 999:    0.1          0     -1         7   998
1000:    0.1          0     -1         7   999
> for(lambda in c(100, 0.1)){
+ fit <- PeakSegOptimal::PeakSegFPOP(count.vec, penalty=lambda)
+ fit.dt <- with(fit, data.table(
+ penalty, path=rev(ends.vec), intervals=t(intervals.mat)))
+ print(fit.dt)
+ }
      penalty  path intervals.V1 intervals.V2
        <num> <int>        <int>        <int>
   1:     100    -1            0            1
   2:     100    -1            2            1
   3:     100    -1            2            1
   4:     100    -1            2            1
   5:     100    -1            2            1
  ---                                        
 996:     100    -1            9            7
 997:     100    -1           10            7
 998:     100    -1           11            7
 999:     100    -1           10            7
1000:     100     0           10            8
      penalty  path intervals.V1 intervals.V2
        <num> <int>        <int>        <int>
   1:     0.1    -1            0            1
   2:     0.1    -1            2            1
   3:     0.1    -1            2            3
   4:     0.1    -1            4            4
   5:     0.1    -1            4            5
  ---                                        
 996:     0.1   994            2            4
 997:     0.1   996            4            4
 998:     0.1   997            5            2
 999:     0.1   998            2            3
1000:     0.1   999            3            3
```

Below we show that an ever increasing synthetic data set, with large
penalty, gives linear number of intervals,

```r
> N=10
> count.vec <- as.integer(2^seq(1,N))
> pen.vec <- c(1000000, 0.1)
> for(lambda in pen.vec){
+ fit <- PeakSegOptimal::PeakSegFPOP(count.vec, penalty=lambda)
+ fit.dt <- with(fit, data.table(
+ lambda, path=rev(ends.vec), intervals=t(intervals.mat)))
+ print(fit.dt)
+ rob.fit <- robseg::Rob_seg(count.vec, lambda=lambda, lthreshold=1)
+ print(with(rob.fit, data.table(lambda, t.est, intervals)))
+ }
    lambda  path intervals.V1 intervals.V2
     <num> <int>        <int>        <int>
 1:  1e+06    -1            0            1
 2:  1e+06    -1            1            1
 3:  1e+06    -1            3            1
 4:  1e+06    -1            4            1
 5:  1e+06    -1            5            1
 6:  1e+06    -1            6            1
 7:  1e+06    -1            7            1
 8:  1e+06    -1            8            1
 9:  1e+06    -1            9            1
10:  1e+06     0            9            1
    lambda t.est intervals
     <num> <int>     <int>
 1:  1e+06    10         2
 2:  1e+06    10         3
 3:  1e+06    10         5
 4:  1e+06    10         7
 5:  1e+06    10         9
 6:  1e+06    10        11
 7:  1e+06    10        13
 8:  1e+06    10        15
 9:  1e+06    10        17
10:  1e+06    10        18
    lambda  path intervals.V1 intervals.V2
     <num> <int>        <int>        <int>
 1:    0.1    -1            0            1
 2:    0.1     0            1            1
 3:    0.1     2            2            2
 4:    0.1     3            4            2
 5:    0.1     4            4            2
 6:    0.1     5            4            2
 7:    0.1     6            4            2
 8:    0.1     7            4            2
 9:    0.1     8            4            2
10:    0.1     9            4            2
    lambda t.est intervals
     <num> <int>     <int>
 1:    0.1     1         2
 2:    0.1     2         4
 3:    0.1     3         5
 4:    0.1     4         5
 5:    0.1     5         5
 6:    0.1     6         5
 7:    0.1     7         5
 8:    0.1     8         5
 9:    0.1     9         5
10:    0.1    10         4
```

# 29 Feb 2024

[figure-min-l1-2d.R](figure-min-l1-2d.R) makes

![2D L1 loss for two data points using sf](figure-min-l1-2d.png)

# 23 Feb 2024

[avocados.py](avocados.py) reads avocados/data_1.pkl etc and saves
avocados/trajectory1458/25.csv.gz etc.

# 7 Feb 2024

For functional pruning L2 loss in sin/cos space, we can use this
formula for finding the intersection of two loss functions
https://www.wolframalpha.com/input?i=a%2Bb*sin%28x%29%2Bc*cos%28x%29%3D0

[figure-loss.R](figure-loss.R) makes

![L1 L2 angular vs sin cos loss comparison](figure-loss.png)

[figure-roots-sin-cos.R](figure-roots-sin-cos.R) makes

![roots of sin cos difference computed via atan](figure-roots-sin-cos.png)

[figure-roots-l1-2d.R](figure-roots-l1-2d.R) makes

![2d l1 angular loss for 2 data points](figure-roots-l1-2d-loss.png)

![2d l1 angular loss difference](figure-roots-l1-2d.png)

