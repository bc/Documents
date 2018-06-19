import rpy2
print(rpy2.__version__)
from rpy2.robjects.packages import importr

# import rpy2's package module
import rpy2.robjects.packages as rpackages

# import R's utility package
utils = rpackages.importr('hitandrun')

# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# R package names
packnames = ('hitandrun', 'ggplot2')

# R vector of strings
from rpy2.robjects.vectors import StrVector

# Selectively install what needs to be install.
# We are fancy, just because we can.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

# import R's "base" package
base = importr('hitandrun')

# import R's "utils" package
utils = importr('utils')

robjects.r('''
        source("R/feasibility_theory.r")

hit_and_run_delta_constrained(max_delta_per_timestep = 0.2,
							  coefficient1 = 379377442.8395462,
							  coefficient2 = -578605325.96008301,
							  constraint1  = -179578150.00039524,
							  n_samples = 1000,
							  prior_point = c(0.50704678,0.64580454),
							  bounds_tuple_of_numeric = list(c(0,1),c(0,1)),
							  plot=FALSE)
        f(3)
        ''')

robjects.IntVector([1, 2, 3])
robjects.FloatVector([1.1, 2.2, 3.3])

res = robjects.StrVector(['abc', 'def'])
>>> print(res.r_repr())


http://rpy2.readthedocs.io/en/version_2.8.x/introduction.html#getting-started