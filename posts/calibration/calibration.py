import pandas as pd
import numpy as np
chana = pd.DataFrame({
  "prediction" : [0.08000933629,0.08511377598,0.05627307204,0.05656443153,0.1215221966,0.158050244,0.25,0.1691452262,0.025,0.2230206613,0.3395864945,0.3380333811,0.75,0.3250713304,0.4101744615,0.4794329163,0.4479620126,0.5023613803,0.64609314,0.6872966505,0.8074552219,0.9677798693,0.3,0.4,0.5,0.6,0.7,np.nan,0.8,0.9],
  "outcome" : [1,1,1,np.nan,0,1,1,1,0,1,1,np.nan,1,1,np.nan,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1]}).dropna().sort_values(by=['prediction'])
chana.head()

import statsmodels.formula.api as smf
from scipy.stats import norm
formulas = [
  "~norm.ppf(prediction)",
  "~prediction",
  "~np.log(prediction)",
  "~1"
]

models = [smf.probit("outcome"+formula, data = chana).fit() for formula in formulas]
aics = [model.aic for model in models]

predictions = chana.reset_index(drop=True).join(
  pd.DataFrame({formula: model.predict() for formula, model in zip(formulas, models)}))


import seaborn as sns
import matplotlib.pylab as plt

fig, ax = plt.subplots()
for formula in formulas:
  ax.plot(predictions["prediction"], predictions[formula],label=formula)
ax.set_xlim(left=0, right=1)
ax.set_ylim(bottom=0, top=1)
ax.legend()
          
plt.show()
plt.clf()

from statsmodels.tools import add_constant
import statsmodels.formula.api as smf
from scipy.special import expit, logit


def logistic(data, x = "prediction", y = "outcome", simple = False):
  """
  Fit and plot a logistic regression for calibration data with pointwise
  95% confidence bounds.
  - x: Name of the variable on the x-axis.
  - y: Name of the variable on the y-axis.
  - simple: If True, uses `x` as covariate. If False, uses `logit(x)`.
  """
  formula = f'{y} ~ {x}' if simple else f'{y} ~ logit({x})'
  model = smf.logit(formula, data = data).fit()
  X = add_constant(pd.DataFrame({x: np.linspace(0, 1, 300)}))
  se = np.sqrt(np.array([y @ model.cov_params() @ y for y in np.array(X)]))
  predictions = model.predict(X) 
  df = pd.DataFrame({
    'x': X[x], 
    'y': predictions, 
    'ymin': expit(logit(predictions) - 1.96*se), 
    'ymax': expit(logit(predictions) + 1.96*se)})
  ax = df.plot(x = 'x', y = 'y')
  ax.fill_between(x = df['x'], y1 = df['ymax'], y2 = df['ymin'], alpha=.2)
  ax.set_xlim(left=0, right=1)
  ax.set_ylim(bottom=0, top=1)
  ax.get_legend().remove()
  
rng = np.random.default_rng(seed = 313)
x = np.linspace(0.01, 0.99, 300)
y = rng.binomial(1, x, size = 300)
data = {"outcome":y, "prediction":x}

logistic(data, simple = True)
#logistic(chana)
plt.show()
plt.clf()



fig, ax = plt.subplots()
formula = "~norm.ppf(prediction)"
ax.plot(predictions["prediction"], predictions[formula],label=formula)
ax.set_xlim(left=0, right=1)
ax.set_ylim(bottom=0, top=1)
ax.legend()


bins = np.linspace(0, 1, 9)
in_bin = pd.cut(chana["prediction"], bins=bins)
counts = chana["outcome"].groupby(in_bin).mean().to_numpy()

means = chana["outcome"].groupby(in_bin).mean().to_numpy()
plt.stairs(means, bins)
plt.plot(chana["prediction"], chana["outcome"], 'o', color='black');
plt.show()


from scipy.special import xlogy

def score(m, x):
  bins = np.linspace(0, 1, m+1)
  in_bin = pd.cut(x, bins=bins)
  counts = pd.Series(x).groupby(in_bin).count().to_numpy()
  n = sum(counts)
  j = np.sum(counts > 0)
  return j / n - np.sum(xlogy(counts, counts)) / n - np.log(m)

