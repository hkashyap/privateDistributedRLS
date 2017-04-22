# privateDistributedRLS
Distributed linear regression under privacy constraints

Regression analysis is used to estimate the unknown eect of changing
one variable over another. While running a linear regression, it is assumed
that there is a linear relationship between regressand and regressor variables
and this relationship is additive. Various regression algorithms are used for
modeling a data set using linear functions as well as estimating the unknown
model parameters from data. The Recursive Least Square(RLS) algorithm
is one of the most popular regression analysis method with extremely fast
convergence rate. But a regression analysis method does not perform well if
the analysis is performed based on a set of observations representing only a
small portion of independent input variables of the system. If a set of agents
collaboratively read the observations of a system, then agents can also per-
form regression analysis collaboratively, i.e. by sharing all their readings and
then forming an estimate. But in some systems, the agents are unwilling
to share their observations due to privacy reasons, but do not mind sharing
estimates with others for getting similar estimates in return. To this end, a
RLS based distributed regression algorithm is developed for wireless sensor
networks (WSNs) whereby sensors exchange local estimates with one-hop
neighbors to consent on the network-wide estimates adaptively and thereby
converging to the estimates formed by a regression analysis based on all ob-
servations. The proposed algorithm has been obtained through recasting of
the weighted linear least-squares cost into a separable form. Numerical sim-
ulations demonstrate that for higher number of 1 hop neighbors, this RLS
based distributed regression algorithm implementation by an agent performs
similar to the RLS algorithm based on full information.
