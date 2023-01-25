import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit, logit
from statsmodels.tools import add_constant
import statsmodels.api as sm
import pandas as pd

x = np.linspace(-2,2)
X = add_constant(x)
eta = X@np.array([-2.0,0.8])
p = expit(eta)
y = np.random.binomial(n=1,p=p)

model = sm.Logit(y, X).fit()
pred = model.predict(X)

se = np.sqrt(np.array([xx@model.cov_params()@xx for xx in X]))

df = pd.DataFrame({'x':x, 'pred':pred, 'ymin': expit(logit(pred) - 1.96*se), 'ymax': expit(logit(pred) + 1.96*se)})

ax = df.plot(x = 'x', y = 'pred')
df.plot(x = 'x', y = 'ymin', linestyle = 'dashed', ax = ax, color = 'C0')
df.plot(x = 'x', y = 'ymax', linestyle = 'dashed', ax = ax, color = 'C0')
