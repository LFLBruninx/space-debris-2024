from ConjunctionUtilities import read_json_file, Pc2D_Foster, remediate_covariance
from scipy.integrate import dblquad
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BreakupUtilities import long_term_propagator, cart2kep, kep2cart
from datetime import datetime


def compute_euclidean_distance(r_sat, r_debris):
    d_eucl = np.linalg.norm(r_sat - r_debris)
    return d_eucl


def compute_mahalanobis_distance(r_sat, r_debris, P_sat, P_debris):
    Psum = P_sat + P_debris
    invP = np.linalg.inv(Psum)
    diff = r_sat - r_debris
    d_maha = float(np.sqrt(np.dot(diff.T, np.dot(invP, diff)))[0, 0])

    return d_maha



# constants

# radius of tumbling blox-wing model
HBR = 2.090746 + 0.05


file_path = 'data\\conjunction_case14.json'

dict = read_json_file(file_path)

sat = dict[36585]
sat_date = sat['UTC']
s_sat = sat['state']
P_sat = sat['covar']

debris = dict[90014]
debris_date = debris['UTC']
s_deb = debris['state']
P_deb = debris['covar']

d_eucl = np.linalg.norm(s_sat[:3] - s_deb[:3])

Pc = Pc2D_Foster(s_sat, P_sat, s_deb, P_deb, HBR)

d_maha = np.sqrt((s_sat[:3] - s_deb[:3]).T @ np.linalg.inv(P_sat[:3, :3] + P_deb[:3, :3]) @ (s_sat[:3] - s_deb[:3]))


print('break')