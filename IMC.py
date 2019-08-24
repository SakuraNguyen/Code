#!/urs/bin/env python

"""
This is the integral from 0 to 1 of exp(-x)
"""
import numpy as np
import matplotlib.pyplot as plt
# UNBIASED
#M = 1000
#rnd = np.random.random(M)
#plt.plot(rnd,np.exp(-rnd))
#plt.savefig("exp.pdf")
#print np.exp(-rnd).sum()/M

# BIASED
#N = 1000
x = np.linspace(0, 1, 1000)
plt.plot(x, np.exp(-x));
pts = np.random.uniform(0,1,(1000, 2))
pts[:, 1] *= np.e
plt.scatter(pts[:, 0], pts[:, 1])
plt.xlim([0,1])
plt.ylim([0, np.e]);
plt.savefig("exp1.pdf")
