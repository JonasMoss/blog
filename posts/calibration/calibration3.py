import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

chana = pd.DataFrame({
  "prediction" : [0.08000933629,0.08511377598,0.05627307204,0.05656443153,0.1215221966,0.158050244,0.25,0.1691452262,0.025,0.2230206613,0.3395864945,0.3380333811,0.75,0.3250713304,0.4101744615,0.4794329163,0.4479620126,0.5023613803,0.64609314,0.6872966505,0.8074552219,0.9677798693,0.3,0.4,0.5,0.6,0.7,np.nan,0.8,0.9],
  "outcome" : [1,1,1,np.nan,0,1,1,1,0,1,1,np.nan,1,1,np.nan,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1]}) \
  .dropna() \
  .sort_values(by=['prediction']) \
  .reset_index(drop=True)
  
def reg_hist(data, m, x = "prediction", y = "outcome"):
  """
  Make a regression histogram from `data` with `m` evenly spaced bins. Returns
  a tuple of `counts`, `means`, and `bins`.
  - x: Name of the variable on the x-axis.
  - y: Name of the variable on the y-axis.
  """
  bins = np.linspace(0, 1, m + 1)
  grouped = data[y].groupby(pd.cut(data[x], bins=bins))
  counts = grouped.count().to_numpy()
  means = grouped.count().to_numpy()
  return (counts, means, bins)


def plot_reg_hist(counts, means, bins, ax, errors = True):
  """
  Plot a regression histogram with or without error bars.
  """
  means_ = np.concatenate([[0], means])
  ax.step(bins, means_)
  if errors:
    counts_ = np.concatenate([[1], counts])
    sd = np.sqrt(means_ * (1 - means_) / counts_) * (counts_ > 0)
    ax.fill_between(
      bins, 
      np.maximum(means_ + sd * 1.96, 0), 
      np.minimum(means_ - sd * 1.96, 1), 
      step="pre", 
      alpha=0.4)


fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(5, 3))
for (i,j),m in zip(product(range(2), range(2)), [3, 5, 7, 10]):
  counts, means, bins = reg_hist(chana, m = m)
  plot_reg_hist(counts, means, bins, ax = axes[i,j])

fig.tight_layout()    
plt.show()

