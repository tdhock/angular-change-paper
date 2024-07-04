## index 629 is plot
from einops import rearrange
import pandas as pd
import numpy as np
data_path = "/home/truong/code/geodesic-change/data/data_1.pkl.ikwt0u"
df = pd.read_pickle(data_path)
# all_signals, shape (n_signals, n_samples, 2)
all_signals = rearrange(
    df.to_numpy(),
    "n_samples (n_signals n_dims) -> n_signals n_samples n_dims",
    n_dims=2,
)
all_signals = np.deg2rad(all_signals)
all_signals[all_signals < 0] += 2 * np.pi
# a list of problematic signals
interesting_signals = [473, 732, 623, 1092, 240, 146, 590, 250, 378, 629]
