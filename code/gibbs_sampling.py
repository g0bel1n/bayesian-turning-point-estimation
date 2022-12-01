import numpy as np
from scipy.stats import invwishart
from scipy.stats import invgamma



#####  hyperpriors

m0= np.array([-0.25, 0.25])
M0 = 2.22*np.identity(2)
k= 0.25
#p=?
#N=?
#g=?
c0=12
C0= 0.1* np.identity(p+2)
g0=1
G0=1
nu= 8
tau=20

##### prior distributions


# beta

# phi is a vector of size k
def I(mu1, mu2, phi):
    if mu1 <= 1 and  mu2 > 0 and np.sum(phi) <1 :
        return 1
    else:
        return 0

Ip = np.identity(p)
phi = np.random.multivariate_normal(mean= np.zeros(p), cov= k*Ip)
mu1, mu2 = np.random.multivariate_normal(mean= m0, cov= M0)


mean_beta = np.array([m0, 0])
cov_beta = np.array([[ M0, 0], [0, np.dot(k*Ip)]])
beta = np.random.multivariate_normal(mean= mean_beta, cov=cov_beta ) *I(mu1, mu2, phi)

# Q
Q = invwishart.pdf(x, df=c0, scale=C0)

#sigma
sigma = invgamma.rvs( scale = g0,  size = G0)
lambda = np.random.gamma( shape = nu/2,  scale = nu/2,N )

# gamma
gamma, gammaZ= np.random.multivariate_normal(mean= np.zeros(g), cov=tau*np.identity(g) )

#epsilon
epsilon = np.random.dirichlet( (5,2), size=None)*np.random.dirichlet( (3,7), size=None)*np.random.dirichlet( (7,3), size=None)*np.random.dirichlet( (3,7), size=None)

# transition distribution of the independent groups ???

##### posterior distributions

#let's define a bit what we work with
# mu1 = np.array([mu11, mu12])
#phi is of dimension K
# theta= np.array([mu1, ..., muK, phi,  Q1, ... QK, sigma, epsilon_etoile, epsilon, gamma, gammaZ]
# where gamma and gammaZ are of dimension K-1


#we have 4 steps we want to do
# so we need to define four functions that are proportional to the probabilities we want to compute




# The tranistion model defines how to move from sigma_current to sigma_new
transition_model = lambda x: [x[0], np.random.normal(x[1], 0.5, (1,))]


def prior(x):
    # x[0] = mu, x[1]=sigma (new or current)
    # returns 1 for all valid values of sigma. Log(1) =0, so it does not affect the summation.
    # returns 0 for all invalid values of sigma (<=0). Log(0)=-infinity, and Log(negative number) is undefined.
    # It makes the new sigma infinitely unlikely.
    if (x[1] <= 0):
        return 0
    return 1


# Computes the likelihood of the data given a sigma (new or current) according to equation (2)
def manual_log_like_normal(x, data):
    # x[0]=mu, x[1]=sigma (new or current)
    # data = the observation
    return np.sum(-np.log(x[1] * np.sqrt(2 * np.pi)) - ((data - x[0]) ** 2) / (2 * x[1] ** 2))


# Same as manual_log_like_normal(x,data), but using scipy implementation. It's pretty slow.
def log_lik_normal(x, data):
    # x[0]=mu, x[1]=sigma (new or current)
    # data = the observation
    return np.sum(np.log(scipy.stats.norm(x[0], x[1]).pdf(data)))


# Defines whether to accept or reject the new sample
def acceptance(x, x_new):
    if x_new > x:
        return True
    else:
        accept = np.random.uniform(0, 1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
        # less likely x_new are less likely to be accepted
        return (accept < (np.exp(x_new - x)))


def metropolis_hastings(likelihood_computer, prior, transition_model, param_init, iterations, data, acceptance_rule):
    # likelihood_computer(x,data): returns the likelihood that these parameters generated the data
    # transition_model(x): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # data: the data that we wish to model
    # acceptance_rule(x,x_new): decides whether to accept or reject the new sample
    x = param_init
    accepted = []
    rejected = []
    for i in range(iterations):
        x_new = transition_model(x)
        x_lik = likelihood_computer(x, data)
        x_new_lik = likelihood_computer(x_new, data)
        if (acceptance_rule(x_lik + np.log(prior(x)), x_new_lik + np.log(prior(x_new)))):
            x = x_new
            accepted.append(x_new)
        else:
            rejected.append(x_new)

    return np.array(accepted), np.array(rejected)









