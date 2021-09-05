import matplotlib.pyplot as plt
import numpy as np


def plot_results():
    results = np.loadtxt('results.txt', delimiter=' ')
    plt.plot(results.T[0], results.T[1], color='red',
             label='actual descent rate')
    plt.axhline(y=0.5, label='target descent rate')
    plt.xlabel('altitude(m)')
    plt.ylabel('descent rate (m/s)')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    plot_results()
