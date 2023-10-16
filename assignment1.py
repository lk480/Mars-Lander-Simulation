import numpy as np
import matplotlib.pyplot as plt

# EULER METHOD


def euler_integration_spring(m=1, k=1, x=0, v=1, t_max=100, dt=0.1, v_0=1):

    # simulation time, timestep and time
    t_array = np.arange(0, t_max, dt)

    # angular frequency of oscillation for mass-spring system
    omega = np.sqrt(k/m)

    # initialise empty lists to record trajectories
    x1_list = []
    v1_list = []

    # Euler Method Integration
    for t in t_array:
        # append current state to trajectories
        x1_list.append(x)
        v1_list.append(v)

        # calculate new position and velocity
        a = -k * x / m
        x = x + dt * v
        v = v + dt * a

    # convert trajectory lists into arrays, so they can be sliced
    x1_array = np.array(x1_list)
    v1_array = np.array(v1_list)

    return t_array, x1_array, v1_array, omega


# VERLET METHOD


def verlet_integration_spring(m=1, k=1, x=0, v=1, t_max=100, dt=0.1):
    # angular frequency of oscillation for mass-spring system
    omega = np.sqrt(k/m)
    # simulation time, timestep and time
    t_array = np.arange(0, t_max, dt)

    # initialise empty lists to record trajectories
    t, euler_position, euler_velocity, o = euler_integration_spring(
        m, k, x, v, t_max, dt)
    x2_list = [x, euler_position[1]]
    v2_list = [x, euler_velocity[1]]

    # Adjust length of t_array to match length of x2_array and v2_array
    t_array2 = t_array[2:]

    # Verlet Method Integration
    for t in t_array2:

        # calculate new position and velocity
        a = -(k * x) / m
        x = 2*(x2_list[-1]) - x2_list[-2] + a*((dt)**2)
        v = (1/dt)*(x-x2_list[-1])

        x2_list.append(x)
        v2_list.append(v)

    # convert trajectory lists into arrays, so they can be sliced
    x2_array = np.array(x2_list)
    v2_array = np.array(v2_list)

    return t_array, x2_array, v2_array, omega


def plot_spring(t, x, v, omega, v0=1):
    plt.figure(2)
    plt.clf()
    plt.xlabel('time(s)')
    plt.grid()
    plt.plot(t, x, label='x (m)')
    plt.plot(t, v, label='v (m/s)')
    plt.plot(t, (v0/omega)*np.sin(omega*t), label='x_analytical')
    plt.plot(t, np.cos(omega*t), label='v_analytical')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    t1, x1, v1, omega1 = euler_integration_spring()
    t2, x2, v2, omega2 = verlet_integration_spring()

    plot_spring(t2, x2, v2, omega2)
    plot_spring(t1, x1, v1, omega1)
