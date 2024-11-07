from msmbuilder.hmm import VonMisesHMM
import sklearn.grid_search
import sys
try:
    cmd, in_csv = sys.argv
except ValueError:
    in_csv = "figure-2d-hmm-sim-data.csv"
import numpy as np
# conda install scikit-learn=0.19.0 #oldest https://anaconda.org/anaconda/scikit-learn/files?page=17
# then conda install omnia::msmbuilder
import pandas as pd
data_df = pd.read_csv(in_csv)
# what is support?
# https://en.wikipedia.org/wiki/Von_Mises_distribution says any interval of length 2*pi
##from msmbuilder.example_datasets import AlanineDipeptide
hmm = VonMisesHMM(n_states=3)
data_array = data_df.to_numpy()
hmm.fit([data_array])
dir(hmm)
hmm.means_ #looks like between -pi and pi
# array([[ 2.0838134,  0.9998006],
#        [-2.7759614, -0.4395614],
#        [ 1.4357036, -0.9330031]], dtype=float32)

# > dcast(true.mean.dt[, angle := ifelse(
# + mean.angle<pi, mean.angle, mean.angle-2*pi
# + )], seg.i ~ dim.i, value.var="angle")
# Key: <seg.i>
#    seg.i         1          2
#    <int>     <num>      <num>
# 1:     1  1.668240 -0.5767475
# 2:     2  2.338123  1.2672049
# 3:     3 -2.683841 -0.6384364
out_prefix = in_csv.replace(".csv", "")
loglik, ( state_array, ) = hmm.predict([data_array])

out_dict = {
    "means":hmm.means_,
    "states":state_array
}
for data_type, out_array in out_dict.items():
    out_df = pd.DataFrame(out_array)
    out_csv = f"{out_prefix}_{data_type}.csv"
    out_df.to_csv(out_csv, index=False)
