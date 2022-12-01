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
gamma= np.random.multivariate_normal(mean= np.zeros(g), cov=tau*np.identity(g) )

#epsilon
epsilon = np.random.dirichlet( (5,2), size=None)*np.random.dirichlet( (3,7), size=None)*np.random.dirichlet( (7,3), size=None)*np.random.dirichlet( (3,7), size=None)

# transition distribution of the independent groups ???

##### posterior distributions











