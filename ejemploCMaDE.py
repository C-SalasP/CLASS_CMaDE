import numpy as np
import matplotlib.pyplot as plt
from classy import Class
import time


t0 = time.perf_counter()  # desde antes de crear params

CMaDEParams = {
    "has_CMaDE": "yes",
    "save_perturbations": "no",
    "gauge": "synchronous",
    "k_c": 1.0,
    "Q_coupling": 1.0,
    "Omega0_CMaDE_dm": 0.25,
    "Omega0_CMaDE_de": 0.75,
    "Omega_g": 5.37815e-5,
    "Omega_ur": 3.71799e-5,
    "Omega_b": 0.0486773,
    "Omega_k": 0.001,
    "h": 0.6781,
    "output": "",
    "background_verbose": 0,
    "thermodynamics_verbose": 0,
}

cosmo = Class()
cosmo.set(CMaDEParams)
cosmo.compute()


z_grid = np.linspace(0.0, 3.0, 301)           
H = np.array([cosmo.Hubble(z) for z in z_grid])  



H_z1100 = cosmo.Hubble(1100.0)
print("H(z=1100) =", H_z1100) #1/mPc
rs_drag = cosmo.rs_drag()   # r_s en la época de drag (baryon drag), en Mpc
print("rs_drag =", rs_drag, "Mpc")




cosmo.struct_cleanup()
cosmo.empty()

t_end = time.perf_counter()
print(f"Tiempo de ejecución: {t_end - t0:.2f} segundos")

