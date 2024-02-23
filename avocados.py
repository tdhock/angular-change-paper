import pandas
import re
import os
import math
import numpy as np
traj_dict = {
    1:1458,
    2:1184
    }
for suffix, traj_id in traj_dict.items():
    data_pkl = f"avocados/data_{suffix}.pkl"
    print(f"unpickling {data_pkl}")
    avocados_df = pandas.read_pickle(data_pkl)
    avocados_df.columns
    # >>> min(avocados_df[0])
    # -179.997620876124
    # >>> max(avocados_df[0])
    # 179.98964711962103
    data_dir = f"avocados/trajectory{traj_id}"
    os.system(f"mkdir -p {data_dir}")
    for col_i in avocados_df.columns:
        one_vec = avocados_df[col_i]
        one_radians = (one_vec+180)*2*math.pi/360
        one_radians.reset_index()
        out_csv = f"{data_dir}/{col_i}.csv.gz"
        print(out_csv)
        one_radians.to_csv(out_csv, header=False, index=False)
# Zoran Stefanic <zoran.stefanic@irb.hr> writes:
# As agreed we are shating with you two dataframes of data (each as a pickled Pandas dataframe, of size 2792x400000 points, which amounts to 8.3 GB of data per dataframe). 
# data_1.pkl represents trajectory: https://alokomp.irb.hr/md/trajectory/1458
# data_2.pkl represents trajectory: https://alokomp.irb.hr/md/trajectory/1184
# Each dataframe represents 400k points in time of 2792 angles. The data is raw, that means it has not been "despiked" using numpy "unwrap"  function for example. The data are in degrees, so you might want to convert it to radians if applying trig functions or similar. For the reference and a bit of a longer story behind it, I attach the html version of Jupyter notebook I have sent you for the first meeting.
# The data is the property of Croatian Science Foundation which finances my project so please use it having that in mind.
