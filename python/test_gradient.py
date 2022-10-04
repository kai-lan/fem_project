import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp

binary = '../build/unit_tests'
# 0  polar decomposition
# 1  gradient of energy (force)
outLines = sp.check_output([binary, '1']).decode().split()
data = [float(x) for x in outLines]
data = np.array(data).reshape((-1, 2))

plt.loglog(data[:, 0], data[:, 1])
plt.title('Finite difference validation (force)')
plt.grid()
plt.xlabel('finite difference step size $\epsilon$')
plt.ylabel('error')
plt.savefig('plot_force.png')
plt.close()