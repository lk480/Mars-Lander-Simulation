import numpy as np
import matplotlib.pyplot as plt


# EULER METHOD


def euler_integration_planetary(position, velocity, M=6.42e23,
                                G=6.67e-11, t_max=100, dt=0.1):

    # simulation time, timestep and time
    t_array = np.arange(0, t_max, dt)

    # initialise empty lists to record trajectories
    euler_position_list = []
    euler_velocity_list = []

    # Euler Method Integration
    for t in t_array:
        # append current state to trajectories
        euler_position_list.append(position)
        euler_velocity_list.append(velocity)

        # calculate new position and velocity
        mod_position = np.sqrt(np.dot(position, position))
        a = (-G * M / (mod_position)**3) * position
        position = position + np.multiply(velocity, dt)
        velocity = velocity + np.multiply(a, dt)

    # convert trajectory lists into arrays, so they can be sliced
    euler_position_array = np.array(euler_position_list)
    euler_velocity_array = np.array(euler_velocity_list)

    return t_array, euler_position_array, euler_velocity_array


# VERLET METHOD


def verlet_integration_planetary(position, velocity, M=6.42e23,
                                 G=6.67e-11, t_max=100, dt=0.1):

    # simulation time, timestep and time
    t_array = np.arange(0, t_max, dt)

    # initialise empty lists to record trajectories
    t, euler_position, euler_velocity = euler_integration_planetary(
        position, velocity, t_max, dt)
    verlet_position_list = [position, euler_position[1]]
    verlet_velocity_list = [velocity, euler_velocity[1]]

    # Adjust length of t_array to match length of x2_array and v2_array
    t_array_adjusted = t_array[2:]

    # Verlet Method Integration
    for t in t_array_adjusted:

        # calculate new position and velocity
        mod_position = np.sqrt(np.dot(position, position))
        a = (-G * M / (mod_position)**3) * position
        position = 2*(verlet_position_list[-1]) - \
            verlet_position_list[-2] + a*(dt)**2
        velocity = (1/dt)*(position-verlet_position_list[-1])

        verlet_position_list.append(position)
        verlet_velocity_list.append(velocity)

    # convert trajectory lists into arrays, so they can be sliced
    verlet_position_array = np.array(verlet_position_list)
    verlet_velocity_array = np.array(verlet_velocity_list)

    return t_array, verlet_position_array, verlet_velocity_array


if __name__ == "__main__":
    intial_p = np.array([0, 4e6, 0])
    intial_v = np.array([0, 0, 0])

    # Scenario 1 - Straight Down Descent
    t, euler_p, euler_v = euler_integration_planetary(intial_p, intial_v)
    t2, verlet_p, verlet_v = verlet_integration_planetary(intial_p, intial_v)
    height_above_surface_euler = [p[1] for p in euler_p]
    height_above_surface_verlet = [p[1] for p in verlet_p]
    # plt.plot(t, height_above_surface_euler)
    # plt.plot(t2, height_above_surface_verlet)
    # plt.show()

    # Scenario 2 - Circular Orbit
    G = 6.67e-11
    M = 6.42e23
    mod_inital_p = np.sqrt(np.dot(intial_p, intial_p))
    # correction factor for elliptical orbit
    v_orbit = np.sqrt(G*M/mod_inital_p) + 5e4
    intial_v_orbit = np.array([v_orbit, 0, 0])
    t, euler_p2, euler_v2 = euler_integration_planetary(
        intial_p, intial_v_orbit, t_max=10000)
    t2, verlet_p2, verlet_v2 = verlet_integration_planetary(
        intial_p, intial_v_orbit)
    x_coord_euler2 = [x[0] for x in euler_p2]
    y_coord_euler2 = [y[1] for y in euler_p2]
    # plt.plot(x_coord_euler2, y_coord_euler2)
    # plt.scatter(0, 0, s=100)
    # plt.show()

    # Scenario 3 - Elliptical Orbit

    G = 6.67e-11
    M = 6.42e23
    semi_major_axis = 1e6
    foci_distance = 1e6
    mod_inital_p = np.sqrt(np.dot(intial_p, intial_p))
    v_elliptical = np.sqrt(G*M/mod_inital_p)
    intial_v_elliptical = np.array([v_elliptical, 0, 0])
    t, euler_p3, euler_v3 = euler_integration_planetary(
        intial_p, intial_v_elliptical, t_max=100000)
    t2, verlet_p3, verlet_v3 = verlet_integration_planetary(
        intial_p, intial_v_elliptical)
    x_coord_euler3 = [x[0] for x in euler_p3]
    y_coord_euler3 = [y[1] for y in euler_p3]
    plt.plot(x_coord_euler3, y_coord_euler3)
    plt.scatter(0, 0, s=100)
    plt.show()

    # Scenario 4 - Hyperbolic Escape Orbit (from Circular Orbit)
    G = 6.67e-11
    M = 6.42e23
    mod_inital_p = np.sqrt(np.dot(intial_p, intial_p))
    v_escape = np.sqrt(2*G*M/mod_inital_p)
    intial_v_escape = np.array([v_escape, 0, 0])
    t, euler_p4, euler_v4 = euler_integration_planetary(
        intial_p, intial_v_escape, t_max=10000)
    t2, verlet_p4, verlet_v4 = verlet_integration_planetary(
        intial_p, intial_v_escape)
    x_coord_euler4 = [x[0] for x in euler_p4]
    y_coord_euler4 = [y[1] for y in euler_p4]
    # plt.plot(x_coord_euler2, y_coord_euler2)
    # plt.plot(x_coord_euler4, y_coord_euler4)
    # plt.scatter(0, 0, s=100)
    # plt.show()
