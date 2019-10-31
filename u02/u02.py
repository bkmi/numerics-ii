import numpy as np
import matplotlib.pyplot as plt

def test_equation(x, l):
    return -l * x

def two_step_method(xk, xkminus1, l, tau):
    return xkminus1 - 2 * tau * l * xk

def solution(t, l, x0):
    return x0*np.exp(-l*t)

lamb = 1.0
x0 = 1.0
tmin = 0.0
tmax = 1.0

grid = np.linspace(tmin, tmax, 1000)
y_soln = solution(grid, lamb, x0)
plt.plot(grid, y_soln, label='Solution', linewidth=3.5)

ns = [10, 100, 1000]
for n in ns:
    grid = np.linspace(tmin, tmax, n)
    print(len(grid))
    tau = 1/n
    x1 = x0 - tau * lamb * x0
    y_twostep = [x0, x1]
    for tk in grid[2:]:
        y_twostep.append(
            two_step_method(y_twostep[-1], y_twostep[-2], lamb, tau)
        )
    plt.plot(grid, y_twostep, label=f'Two step n={n}')
plt.legend()
plt.savefig(f'compare.png')