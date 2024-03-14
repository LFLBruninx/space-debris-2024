import pickle
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from BreakupUtilities import compute_number
from matplotlib.ticker import PercentFormatter


# Constants

# Impact Parameters
mass_target = 1412#1630.
mass_impactor = 8.52942
impact_velocity = 4.45300*1000.
expl_flag = False

# Timing parameters
t_total = 10 # years
t_step = 1000000 # seconds
t_step_years = t_step / (365.25 * 24 * 3600)
t_array = np.arange(0, t_total, t_step_years)





def plot_objects_per_bin(m_tar, m_imp, v_imp):
    '''
    Test case to compute and plot number of objects generated in breakup

    Parameters
    ----------
    expl_flag : boolean
        flag to determine if using explosion (True) or collision (False) model

    Returns
    -------
    None.

    '''

    # Setup OneWeb case per Radtke 2017
    lc_array = np.logspace(-2, 1, 1000)
    cs = 1.

    total = 0
    N_bin = np.zeros(lc_array.shape)
    N_cum = np.zeros(lc_array.shape)
    for ii in range(len(lc_array) - 1, -1, -1):
        lc = lc_array[ii]
        N_lc = compute_number(lc, v_imp, m_tar, m_imp, expl_flag, cs)

        N_bin[ii] = math.floor(N_lc - total)
        N_cum[ii] = N_lc
        total += N_bin[ii]

    N_bin = [int(N) for N in N_bin]

    print('Total Number of Debris Particles: ', sum(N_bin))
    print('Check N cumulative', N_cum[0])

    plt.figure()
    plt.semilogx(lc_array, N_cum, 'k.')
    plt.xlabel('Characteristic Length [m]')
    plt.ylabel('Cumulative Number of Objects')

    plt.show()

    plt.figure(figsize=(5, 4))
    plt.semilogx(lc_array, N_bin, 'k.')
    plt.xlabel('Characteristic Length [m]')
    plt.ylabel('Number of Objects Per Bin')
    plt.savefig('plots\objects_per_bin.pdf')
    plt.show()

    return



def read_pkl(filename):
    with open(filename, 'rb') as pklFile:
        kep_array_initial, kep_array_final, reentry_times, latlon_array_initial, latlon_array_final = pickle.load(
            pklFile)
    return kep_array_initial, kep_array_final, reentry_times, latlon_array_initial, latlon_array_final

# Example usage
file_path_90d = 'output/breakup_data_90d.pkl'  # Replace with your actual file path
file_path_10y = 'output/breakup_data_10y.pkl'
########################################################################
# Read pickle file #####################################################
########################################################################

kep_array_initial, kep_array_90d, reentry_times_90d, latlon_array_initial, latlon_array_90d = read_pkl(file_path_90d)
kep_array_initial, kep_array_final, reentry_times_final, latlon_array_initial, latlon_array_final = read_pkl(file_path_10y)

########################################################################
# Plot number of objects per bin #######################################
########################################################################

plot_objects_per_bin(mass_target, mass_impactor, impact_velocity)

########################################################################
# Make plots for question 1 ############################################
########################################################################

SMA_initial = kep_array_initial[:, 0]
ECC_initial = kep_array_initial[:, 1]
INC_initial = kep_array_initial[:, 2]
RAAN_initial = kep_array_initial[:, 3]
AOP_initial = kep_array_initial[:, 4]
hp_initial = SMA_initial * (1 - ECC_initial) - 6378137
ha_initial = SMA_initial * (1 + ECC_initial) - 6378137

SMA_90d = kep_array_90d[:, 0]
ECC_90d = kep_array_90d[:, 1]
INC_90d = kep_array_90d[:, 2]
RAAN_90d = kep_array_90d[:, 3]
AOP_90d = kep_array_90d[:, 4]
hp_90d = SMA_90d * (1 - ECC_90d) - 6378137
ha_90d = SMA_90d * (1 + ECC_90d) - 6378137

SMA_final = kep_array_final[:, 0]
ECC_final = kep_array_final[:, 1]
INC_final = kep_array_final[:, 2]
RAAN_final = kep_array_final[:, 3]
AOP_final = kep_array_final[:, 4]
hp_final = SMA_final * (1 - ECC_final) - 6378137
ha_final = SMA_final * (1 + ECC_final) - 6378137

# how many orbits ended with a hp_final > 1000 km
num_non_deorbited = np.sum(hp_final > 1000000)

# create an array containing always the same value num_non_deorbited the size of the t_array
num_non_deorbited_array = np.ones(t_array.shape)

# plot the number of non-deorbited orbits as a function of time
# only show the y axis between 0 and 100
# add percentage to the y axis
plt.figure()
plt.plot(t_array, num_non_deorbited_array * 100, "k")
plt.ylim(0, 105)
plt.xlabel('Time [years]')
plt.ylabel('Percentage of objects remaining in orbit')
plt.gca().yaxis.set_major_formatter(PercentFormatter())
plt.savefig('plots\percentage_objects_remaining.pdf')
plt.show()


# make a Gabbard plot for the initial state plotting the h_p and h_a for each object over the SMA
plt.figure()
plt.plot(SMA_initial/1000, hp_initial/1000, "k.", label='Perigee Altitude')
plt.plot(SMA_initial/1000, ha_initial/1000, "r.", label='Apogee Altitude')
plt.xlabel('Semi-major axis [km]')
plt.ylabel('Altitude [km]')
plt.legend()
plt.title('Gabbard plot for the initial state')
plt.savefig('plots\gabbard_initial.pdf')
plt.show()

# make an incliniation versus raan plot for the initial state
plt.figure()
plt.plot(np.degrees(RAAN_initial), np.degrees(INC_initial), "k.")
plt.xlabel('RAAN [degrees]')
plt.ylabel('Inclination [degrees]')
plt.title('RAAN vs inclination for the initial state')
plt.xlim(0, 300)
plt.ylim(35, 60)
plt.savefig('plots\\raan_vs_inc_initial.pdf')
plt.show()

# make a Gabbard plot for the 90d state plotting the h_p and h_a for each object over the SMA
plt.figure()    
plt.plot(SMA_90d/1000, hp_90d/1000, "k.", label='Perigee Altitude')
plt.plot(SMA_90d/1000, ha_90d/1000, "r.", label='Apogee Altitude')
plt.xlabel('Semi-major axis [km]')
plt.ylabel('Altitude [km]')
plt.legend()
plt.title('Gabbard plot after 90 days')
plt.savefig('plots\gabbard_90d.pdf')
plt.show()

# make an incliniation versus raan plot for the 90d state
plt.figure()
plt.plot(np.degrees(RAAN_90d), np.degrees(INC_90d), "k.")
plt.xlabel('RAAN [degrees]')
plt.ylabel('Inclination [degrees]')
plt.title('RAAN vs inclination after 90 days')
plt.xlim(0, 300)
plt.ylim(35, 60)
plt.savefig('plots\\raan_vs_inc_90d.pdf')
plt.show()

# make a Gabbard plot for the final state plotting the h_p and h_a for each object over the SMA
plt.figure()
plt.plot(SMA_final/1000, hp_final/1000, "k.", label='Perigee Altitude')
plt.plot(SMA_final/1000, ha_final/1000, "r.", label='Apogee Altitude')
plt.xlabel('Semi-major axis [km]')
plt.ylabel('Altitude [km]')
plt.legend()
plt.title('Gabbard plot after 10 years')
plt.savefig('plots\gabbard_final.pdf')
plt.show()

# make an incliniation versus raan plot for the final state
plt.figure()
plt.plot(np.degrees(RAAN_final), np.degrees(INC_final), "k.")
plt.xlabel('RAAN [degrees]')
plt.ylabel('Inclination [degrees]')
plt.title('RAAN vs inclination after 10 years')
plt.xlim(0, 300)
plt.ylim(35, 60)
plt.savefig('plots\\raan_vs_inc_final.pdf')
plt.show()






print('break')

