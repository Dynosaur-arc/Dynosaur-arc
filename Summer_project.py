#Constants and importing packages
import math
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
h_bar = 6.6 * 10**(-34)/(2*math.pi)
m = float(input("Mass of particle:"))
k = -1*h_bar*h_bar/(2*m)
grid_length = 10**(-2)
grid_length_time = 10**(-2)




def function_maker(index):  # Must return list that contains the values of the stationary states of SHO evaluated at 500 uniformly spaced points in [-2,3)
    funct = np.empty([1,500])
    for i in range(-200,300):
        funct[0][i+200] = 1e5*math.sin(index*i/100)*math.exp(-(i/100)**2)
    return funct

def orthonormal_basis(dummy_var):
    orthonormal_basis = np.empty([500, 500])
    for i in range(0, 500):
        funct = function_maker(i)
        for j in range(0, 500):
            orthonormal_basis[i][j] = funct[0][j]
    return orthonormal_basis


orthonormal_basis = orthonormal_basis(1)

#print(orthonormal_basis)

def hamiltonian(f1,f2,f3,grid_length):
  return (k/grid_length**2) * (f1 + f3 - 2*f2)


def hamiltonian_integral(f1,f2,grid_length):
  integral = 0
  for i in range(1,499):
    integral = integral + grid_length * f1[i] * hamiltonian(f2[i-1],f2[i],f2[i+1],grid_length)
  return integral

def hamiltonian_matrix(orthonormal_basis,grid_length):
  hamiltonian_matrix = np.empty([500,500])
  for i in range(0,500):
    for j in range(0,500):
      orthonormal_basis[i][j] = hamiltonian_integral(orthonormal_basis[i],orthonormal_basis[j],grid_length)
  return hamiltonian_matrix

hamiltonian_matrix = hamiltonian_matrix(orthonormal_basis,grid_length)

e_value,coeff_matrix = np.linalg.eigh(hamiltonian_matrix)

#print(coeff_matrix)

def stationary_state(coeff_matrix,orthonormal_basis,n):
  stationary_state = np.empty([1,500])
  for i in range(0,500):

    sum = 0
    for j in range(0,500):
      sum = sum + coeff_matrix[j][n] * orthonormal_basis[j][i]
    stationary_state[0][i] = sum
  return stationary_state

n = int(input("Enter index of stationary state to be evaluated:"))
stationary_state = stationary_state(coeff_matrix,orthonormal_basis,n)

#print(stationary_state)

def time_stepper(t, stationary_state, well_potential, k, h_bar, grid_length, grid_length_time):
    wave_function = stationary_state
    if (t > 0):
        for j in range(0,int(t/grid_length_time)):
            next_state = wave_function
            for i in range(0, 498):

                if i < 200 and i > 300:
                    v = well_potential
                else:
                    v = 0
                #k = 1 #to prevent overflow
                #a = -1 * grid_length_time / (grid_length) ** 2
                #c = a
                #b = -1 * grid_length_time * (-2 * k / (grid_length ** 2) + v)

                wave_function[0][i + 1] = (1/2)*(-1 * next_state[0][i] + (2-v*grid_length_time) * next_state[0][i + 1] + -1 * next_state[0][i + 2])

    return (wave_function)


t = float(input("Enter time value:"))
well_potential = float(input("Enter value of well potential:"))
wave_function = time_stepper(t, stationary_state, well_potential, k, h_bar, grid_length, grid_length_time)
#print(wave_function)
def plotter(wave_function):
  for i in range(-199,298):
    plt.plot(i/100,wave_function[0][i+200],'r.')
  plt.show()

#print(wave_function[0][1])
plotter(wave_function)

#k = int(input("press 1 to generate plot at a different time:"))

#while(k==1):
    #t = float(input("enter time:"))
    #well_potential = float(input("enter well potential"))
    #wave_function = time_stepper(t, stationary_state, well_potential, k, h_bar, grid_length, grid_length_time)
    #plotter(wave_function)
    #k = int(input("enter 1 to continue plotting at different time values and well potentials:"))




