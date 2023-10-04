import numpy as np
x_cos = [np.cos(i) for i in [0, np.pi*0.25, np.pi*0.5, np.pi*0.75, np.pi, np.pi*1.25, np.pi*1.5, np.pi*1.75]]
for i, n in enumerate(x_cos):
    if abs(n) < 0.0001:
        x_cos[i] = 0
x = np.array([3,6,9,12,15])
x_cos = np.array(x_cos)
print(0.5*(x[1:] + x[:-1]))
print(0.5*(x_cos[1:] + x_cos[:-1]))