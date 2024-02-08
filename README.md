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

