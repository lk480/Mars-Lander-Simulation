#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main()
{

  // declare variables
  double m, k, x, v, t_max, dt, t, a;
  vector<double> euler_t_list, euler_x_list, euler_v_list, verlet_t_list, verlet_x_list, verlet_v_list;

  // mass, spring constant, initial position and velocity
  m = 1;
  k = 1;
  x = 0;
  v = 1;

  // simulation time and timestep
  t_max = 100;
  dt = 0.1;

  // Euler integration
  for (t = 0; t <= t_max; t = t + dt)
  {

    // append current state to trajectories
    euler_t_list.push_back(t);
    euler_x_list.push_back(x);
    euler_v_list.push_back(v);

    // calculate new position and velocity
    a = -k * x / m;
    x = x + dt * v;
    v = v + dt * a;
  }

  // Verlet Integration

  verlet_x_list.push_back(x);
  verlet_x_list.push_back(euler_x_list[1]);
  verlet_v_list.push_back(x);
  verlet_v_list.push_back(euler_v_list[1]);

  for (t = 0; t <= t_max; t = t + dt)
  {
    // calculate new position and velocity
    a = -(k * x) / m;
    x = 2 * (verlet_x_list[-1]) - verlet_x_list[-2] + a * (dt * dt);
    v = (1 / dt) * (x - verlet_x_list[-1]);

    verlet_x_list.push_back(x);
    verlet_v_list.push_back(v);
  }

  // Write the trajectories to file
  ofstream fout;
  fout.open("trajectories.txt");
  if (fout)
  { // file opened successfully
    for (int i = 0; i < euler_t_list.size(); i++)
    {
      fout << euler_t_list[i] << ' ' << euler_x_list[i] << ' ' << euler_v_list[i] << endl;
    }
  }
  else
  { // file did not open successfully
    cout << "Could not open trajectory file for writing" << endl;
  }

  /* The file can be loaded and visualised in Python as follows:

  import numpy as np
  import matplotlib.pyplot as plt
  results = np.loadtxt('trajectories_euler.txt')
  plt.figure(1)
  plt.clf()
  plt.xlabel('time (s)')
  plt.grid()
  plt.plot(results[:, 0], results[:, 1], label='x (m)')
  plt.plot(results[:, 0], results[:, 2], label='v (m/s)')
  plt.legend()
  plt.show()

  Alternatively, using pandas dataframes
  import pandas as pd
  df = pd.read_csv("/Users/lohithkonathala1/git/IB-Mars-Lander-cpp/trajectories.txt", sep= ' ', index_col=0, names=['Time', 'Euler', 'Verlet'])
  %matplotlib
  df.plot()
  */
}
