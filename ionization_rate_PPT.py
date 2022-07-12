import numpy as np
factorial = np.math.factorial

import scipy
gamma = scipy.math.gamma

import matplotlib.pyplot as plt

############ Ultrashort filaments of light in weakly ionizzed optically transparent media Luc Berge
########### expression of PPT
def f(l,m):
    return (2*l + 1) * factorial(l + np.abs(m)) / 2**(np.abs(m)) / factorial(np.abs(m)) / factorial(l - np.abs(m))

def phi_m(x, m):
    y = np.linspace(0, x, 10000)
    dy = y[1] - y[0]
    return np.exp(-x**2) * np.trapz((x**2 - y**2)**np.abs(m) * np.exp(y**2), y, dy)

## h bar
h_bar = 1.054e-34 # J * s
h_bar_ev = 6.582e-16 # eV * s

## electron charge and mass
e = 1.6e-19 # C
m_e = 9.1e-31 # kg

## bohr radius
a_b = 5.2e-11 # m

## Laser freq
c = 3e8 # m/s
wl = 800e-9 # m
k0 = 2 * np.pi / wl # 1/m
laser_freq = k0 * c # 1/s
energy_of_single_photon_in_ev = h_bar * laser_freq # eV

## Laser intensities
I_start = 1e12 # W/cm^2
I_end = 5e13 # W/cm^2
I_res = 100
I = np.linspace(I_start, I_end, I_res)

## from intensity to electric field
## This is from "Gaussian Beam" in wikipedia
etta = 377 # ohm = kg * m^2 * s^-1 * C^-2
etta *= 10000 # from m^2 to cm^2
E_p = np.sqrt(2 * etta * I) # N / C (SI)

## Ponderamotive energy
U_p = e**2 * E_p **2 / 4 / m_e / laser_freq**2

## Ionization potential for oxygen
U_i_eV = 12.07 # eV
U_i = U_i_eV * 1.6e-19 # J

## Coulmb field
## Coulmb force exerted on electron in hydrogen atom is 8.2e-8 N
## https://www.youtube.com/watch?v=exPlfD1dRxc&ab_channel=CyrusVandrevala
## in wolfram for hydrogen: (2 * 13.6eV)^1.5 * (9.1e-31kg)^0.5 * (h bar)^(-1)
E_0 = (2 * U_i) ** 1.5 * m_e**0.5 / h_bar / e

## required photons for MPI
nu = U_i_eV / (h_bar_ev * laser_freq)
nu_zero = np.ceil(nu)

## Keldysh Paremeter - no units
keldysh = laser_freq * np.sqrt(2 * m_e * U_i) / E_p / e

## alpha - no units
alpha = 2 * (np.arcsinh(keldysh) - keldysh / np.sqrt(1 + keldysh**2))

##  beta - no units
beta = 2 * keldysh / np.sqrt(1 + keldysh**2)


## effective quantum numbers (n star and l star)
Z = 8 # oxygen
n_star = Z / np.sqrt(2*U_i) * (a_b * m_e**0.5 / h_bar)**(-1) ## term on right is conversion from atomic units
l_star = n_star - 1
m = 0

C1 = 0
C2 = 0.683
C3 = 0
C4 = 0.033
C_nl_sqr_fl = C1/(2*U_i * h_bar**-2 * m_e * a_b**2)**(n_star/2 + 0.25)*f(1,0) + \
              C2/(2*U_i * h_bar**-2 * m_e * a_b**2)**(n_star/2 + 0.25)*f(2,0) + \
              C3/(2*U_i * h_bar**-2 * m_e * a_b**2)**(n_star/2 + 0.25)*f(3,0) + \
              C4/(2*U_i * h_bar**-2 * m_e * a_b**2)**(n_star/2 + 0.25)*f(4,0)

#### calculate the sum in the W expression
inf_sum = 100
k = np.arange(nu_zero, inf_sum)
sum_arr = np.zeros(I_res)
for i in range(I_res):
    sum_of_elements = 0
    for j in k:
        sum_of_elements += np.exp(-alpha[i] * (j - nu)) * phi_m(np.sqrt(beta[i] * (j - nu)), m)
    sum_arr[i] = sum_of_elements

####### Expression for Ionization rate PPT
W = 4 * np.sqrt(2) / np.pi  * (2 * E_0 / E_p / np.sqrt(1 + keldysh**2))**(2*n_star - 3/2 - np.abs(m)) * \
    C_nl_sqr_fl / factorial(np.abs(m)) * \
    np.exp(-2 * nu * (np.arcsinh(keldysh) - keldysh * np.sqrt(1 + keldysh**2) / (1 + 2*keldysh**2))) * \
    U_i * keldysh**2 / (1 + keldysh**2) * sum_arr / h_bar

plt.figure()
plt.plot(I, W, label = 'PPT')
## plt.plot(np.log(I), np.log(W), label = 'PPT')
plt.title('Ionization Rate vs I')
plt.xlabel('I (W/cm^2)')
plt.ylabel('W')

sigma_8 = 3e-98
plt.plot(I, sigma_8 * I**8, label = 'MPI')
plt.legend()
