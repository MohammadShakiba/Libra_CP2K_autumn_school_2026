"""NBRA nonadiabatic dynamics, adapted from NAMD_MB.py and 9_step4/tutorial.ipynb.

Run from 4_step4: python NAMD_simulation.py
Edit the assignments below. Optionally set the run folder with --prefix PATH.

1. Leave basis_size = None and sys.exit(0) active for the initial basis test.
2. Set basis_size to your chosen number of states, comment out sys.exit(0),
   select one recipe, and rerun to continue into dynamics.
3. Edit and run NAMD_plotting.py to analyze the saved data and populations.

Learning goals: set up file-based dynamics with different hopping/decoherence
methods, assess the adiabatic basis, and interpret populations and timescales.
Libra's generic_recipe is the same driver used for model-Hamiltonian dynamics;
here compute_model supplies precomputed energies and time overlaps instead.
These inputs can come from different electronic-structure codes if they use
consistent units, state ordering, and overlap conventions.

NBRA follows a precomputed nuclear trajectory without electronic back-reaction.
All matrix entries are in atomic units. DT is the spacing of input frames in fs.
The last file index fstep is excluded; states and initial conditions start at 0.
"""
import os
import sys
import argparse
import multiprocessing as mp
import time

parser = argparse.ArgumentParser(description="NBRA simulation and basis diagnostics")
parser.add_argument("--prefix", default="namd_run",
                    help="Folder for all results of this run (default: namd_run)")
args = parser.parse_args()

import numpy as np
import scipy.sparse as sp
from liblibra_core import CMATRIX, MATRIX, Random, Py2Cpp_int, pop_submatrix
from libra_py import units, data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.workflows.nbra.decoherence_times as decoherence_times
from ld_active_space import diagnose_fixed_active_spaces
from recipes import fssh2_nbra, gfsh_nbra, ida_nbra
from recipes import msdm_nbra, mash_nbra, dish_nbra, dish_rev2023_nbra
from recipes import fssh_nbra        # LD as described in our paper
from recipes import fssh_nbra_v2     # LD based on SVD procedure - more in line with Granucci and Persico
from recipes import fssh_nbra_nac_v0 # no state tracking, no scaling, no phase correction
from recipes import fssh_nbra_nac_v1 # Hungarian algorithm, no scaling, no phase correction
from recipes import fssh_nbra_nac_v2 # Hungarian algorithm, no scaling, with phase correction
from recipes import fssh_nbra_nac_v3 # Hungarian algorithm, scaling 1, with phase correction
from recipes import fssh_nbra_nac_v4 # LD for state tracking, no scaling, with phase correction



# ==================== 1. Read the input matrices ====================
# Every required file from istep through fstep-1 must exist. Both readers return
# NSTEPS = fstep-istep, NSTATES = retained basis size, and five lists of CMATRIX:
#   St_adi          overlaps between adiabatic states at successive geometries
#   Ham_adi         adiabatic electronic Hamiltonians (energies on the diagonal)
#   Hvib_adi        vibronic Hamiltonians, including nonadiabatic couplings
#   NAC_adi         time-derivative nonadiabatic coupling matrices
#   Basis_transform representation transforms (identity in this file interface)
# Store these once and reuse them at each dynamics step; avoid repeated disk
# reads and NumPy-to-CMATRIX conversions inside the model callback.
def read_files_new(data_path, istep, fstep, dt, basis_size=None):
    """Read text matrices and retain the first basis_size states.

    None retains all states. Libra reads both real entries and (real,imag)
    entries in these text files. HST gives NAC = (St - St.H())/(2*dt),
    and the vibronic Hamiltonian is Hvib = Ham - i*NAC.
    """
    nfull = np.loadtxt(f"{data_path}/st_adi_{istep}.txt", dtype=str, ndmin=2).shape[0]
    nstates = nfull if basis_size is None else basis_size
    assert 1 <= nstates <= nfull
    states = Py2Cpp_int(list(range(nstates)))
    St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = [], [], [], [], []

    for step in range(istep, fstep):
        st_full = CMATRIX(nfull, nfull)
        ham_full = CMATRIX(nfull, nfull)
        st_full.Load_Matrix_From_File(f"{data_path}/st_adi_{step}.txt")
        ham_full.Load_Matrix_From_File(f"{data_path}/ham_adi_{step}.txt")
        st, ham = CMATRIX(nstates, nstates), CMATRIX(nstates, nstates)
        pop_submatrix(st_full, st, states, states)
        pop_submatrix(ham_full, ham, states, states)
        nac = (st - st.H()) / (2.0 * dt * units.fs2au)
        hvib = ham - nac * 1j
        basis = CMATRIX(nstates, nstates)
        basis.identity()
        St_adi.append(st)
        Ham_adi.append(ham)
        Hvib_adi.append(hvib)
        NAC_adi.append(nac)
        Basis_transform.append(basis)

    return fstep - istep, nstates, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform


def read_files_old(data_path, istep, fstep, dt, basis_size=None):
    """Read the notebook's sparse NPZ format instead of text files.

    Read energies from Hvib_ci_INDEX_re.npz and overlaps from St_ci_INDEX_re.npz.
    Recompute HST NACs from these overlaps, as required by the selected recipes,
    rather than mixing them with the old Hvib_ci_INDEX_im.npz sign convention.
    """
    nfull = sp.load_npz(f"{data_path}/St_ci_{istep}_re.npz").shape[0]
    nstates = nfull if basis_size is None else basis_size
    assert 1 <= nstates <= nfull
    St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = [], [], [], [], []
    for step in range(istep, fstep):
        st = sp.load_npz(f"{data_path}/St_ci_{step}_re.npz").toarray()
        ham = sp.load_npz(f"{data_path}/Hvib_ci_{step}_re.npz").toarray()
        st = data_conv.nparray2CMATRIX(st[:nstates, :nstates].astype(complex))
        ham = data_conv.nparray2CMATRIX(np.diag(np.diag(ham)[:nstates]).astype(complex))
        nac = (st - st.H()) / (2.0 * dt * units.fs2au)
        basis = CMATRIX(nstates, nstates)
        basis.identity()
        St_adi.append(st)
        Ham_adi.append(ham)
        Hvib_adi.append(ham - nac * 1j)
        NAC_adi.append(nac)
        Basis_transform.append(basis)
    return fstep - istep, nstates, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform


# ==================== 2. Choose the data and basis ====================
# Dataset descriptions and example settings from 9_step4/tutorial.ipynb.
# Run from 4_step4: the data folders and archives are now beside this script.
# Keep exactly ONE dataset block active: comment the current block's code and
# uncomment the code for the desired example, including its reader call.
# Start with basis_size = None and leave sys.exit(0) active for the basis test.
# Then choose a size and comment out sys.exit(0) to continue into dynamics.
# The commented numerical sizes below are the notebook's example choices.

# -------------------- Example 1: old-style NPZ files --------------------
# Prepared in this tutorial:
# https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows/8_step3/1_DFT
# 11 states; dt = 1.0 fs; istep = 1200, fstep = 1399 (exclusive).
# Archive: data1.tar.bz2 (called "data.tar.bz2" in the notebook description).
# If needed, extract in a terminal: tar -xf data1.tar.bz2
#data_path = "data1"
#istep = 1200                     # First input file index
#fstep = 1399                     # Read through file 1398
#ISTATE = 10                      # Initially populated state, zero-based
#DT = 1.0                         # Input frame spacing in fs
#basis_size = None                # First run: all available states
#basis_size = 11                  # Notebook example; select after the test
#NSTEPS, NSTATES, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = read_files_old(
#    data_path, istep, fstep, DT, basis_size)

# -------------------- Example 2: (TiO2)2, Miguel Recio --------------------
# New-style text files, archive data2.tar.bz2.
# Notebook description: 10 states; dt = 1.0 fs; istep = 0, fstep = 399.
# The supplied matrices actually contain 11 states; None reads all 11, while
# basis_size = 10 reproduces the notebook's selected active-space size.
# Prepared in this tutorial:
# https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/11_program_specific_methods/3_cp2k_methods/5_namd_workflow/tutorial2-tio2.ipynb
# If needed, extract in a terminal: tar -xf data2.tar.bz2
#data_path = "data2"
#istep = 0                        # First input file index
#fstep = 399                      # Read through file 398
#ISTATE = 2                       # Initially populated state, zero-based
#DT = 1.0                         # Input frame spacing in fs
#basis_size = None                # First run: all available states
#basis_size = 10                  # Notebook example; select after the test
#NSTEPS, NSTATES, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = read_files_new(
#    data_path, istep, fstep, DT, basis_size)

# -------------------- Example 3: Mg@C60, Kosar Yasin --------------------
# New-style text files, archive data3.tar.bz2; input files are in data3/res.
# 4 states; dt = 0.5 fs; istep = 0, fstep = 49 (exclusive).
# If needed, extract in a terminal: tar -xf data3.tar.bz2
# This example is selected initially.
data_path = "data3/res"
istep = 0                        # First input file index
fstep = 49                       # Read through file 48
ISTATE = 1                       # Initially populated state, zero-based
DT = 0.5                         # Input frame spacing in fs
#basis_size = None                # First run: all available states
basis_size = 4                   # Notebook example; select after the test
NSTEPS, NSTATES, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = read_files_new(
    data_path, istep, fstep, DT, basis_size)

# -------------------- Common output for the selected dataset --------------------
output_path = args.prefix        # Use the same --prefix with NAMD_plotting.py
os.makedirs(output_path, exist_ok=True)

# NumPy copies make determinants, matrix statistics, and plotting convenient.
# Keep complex entries; take real energies or absolute couplings only when
# forming the particular observable, rather than discarding phases at input.
numpy_st_adi = np.array([data_conv.MATRIX2nparray(st) for st in St_adi])
numpy_ham_adi = np.array([data_conv.MATRIX2nparray(ham) for ham in Ham_adi])
numpy_nac_adi = np.array([data_conv.MATRIX2nparray(nac) for nac in NAC_adi])
# The plotting file reads this archive; it does not import or rerun this script.
np.savez(f"{output_path}/analysis_data.npz", overlaps=numpy_st_adi,
         hamiltonians=numpy_ham_adi, nacs=numpy_nac_adi, dt_fs=DT)


# ==================== 3. Check the adiabatic basis ====================
# A small determinant can indicate states projecting outside the chosen space,
# but determinants also decrease with dimension. Inspect singular values too.
# The rank test requires every singular value to exceed
# max(absolute_tolerance, relative_tolerance * largest singular value)
# at EVERY interval. Passing is a numerical-support test, not proof that
# populations are converged with respect to the number of states.
# The notebook calls a trial size nact_st: it is the same leading-state count
# as basis_size. Determinants need not vary monotonically with size. Local
# maxima of the minimum determinant can suggest candidate spaces, but inspect
# the singular-value support and population convergence before choosing one.
# Near-singular overlaps make LD basis transport unreliable; start with None
# so that the comparison includes every size available in the stored data.
active_sizes = list(range(1, NSTATES + 1))
res = diagnose_fixed_active_spaces(numpy_st_adi, active_sizes,
                                   relative_tolerance=1.0e-4,
                                   absolute_tolerance=1.0e-12)
print("Size   Accepted   Minimum rank   Minimum singular value   Minimum |det(St)|")
basis_summary = []
for size in active_sizes:
    determinants = np.abs(np.linalg.det(numpy_st_adi[:, :size, :size]))
    print(size, res[size]["accepted"], np.min(res[size]["effective_rank"]),
          np.min(res[size]["sigma_min"]), np.min(determinants))
    basis_summary.append([size, int(res[size]["accepted"]),
                          np.min(res[size]["effective_rank"]),
                          np.min(res[size]["sigma_min"]), np.min(determinants)])
np.savetxt(f"{output_path}/basis_diagnostics.txt", basis_summary,
           header="size accepted min_rank min_singular_value min_abs_determinant")

# Keep this line active for the initial basis test.
# After choosing basis_size, comment it out to continue into sections 4-7.
#sys.exit(0)


# ==================== 4. Energy gaps and decoherence times ====================
# Average gap fluctuations give the dephasing times/rates used by MSDM/DISH.
# Rates are in inverse atomic time, gaps in Hartree, displayed times in fs.
# energy_gaps_ave supplies |Ei(t)-Ej(t)|; the loop below averages over time,
# so gaps[i,j] = <|Ei-Ej|>, not |<Ei>-<Ej>|. ave_gaps is needed when a recipe
# enables dephasing_informed. decoherence_times_ave estimates pairwise times
# and inverse times from the trajectory's gap fluctuations.
# Computing these matrices does not itself turn on decoherence: the recipe's
# decoherence algorithm and times option determine how the rates are used.
tau, rates = decoherence_times.decoherence_times_ave([Hvib_adi], [0], NSTEPS, 0)
dE = decoherence_times.energy_gaps_ave([Ham_adi], [0], NSTEPS)
gaps = MATRIX(NSTATES, NSTATES)
for step in range(NSTEPS):
    gaps += dE[step]
gaps /= NSTEPS
deco_times = data_conv.MATRIX2nparray(tau).real * units.au2fs
np.fill_diagonal(deco_times, 0)  # Display only: no self-decoherence
np.savetxt(f"{output_path}/decoherence_times.txt", deco_times)
rates.show_matrix(f"{output_path}/decoherence_rates.txt")
gaps.show_matrix(f"{output_path}/average_gaps.txt")


# ==================== 5. Interface to the generic dynamics driver ====================
class abstr_class:
    pass


def compute_model(q, params, full_id):
    """Return the specifically named attributes expected by the Libra driver.

    q and full_id belong to the generic interface; this file-based model uses
    params['timestep'] to retrieve its data. The returned object connects the
    stored trajectory to the same propagation machinery used for model systems.
    Copies keep trajectory-specific corrections from changing the source lists.
    """
    # Libra updates timestep using icond and nfiles; do not add icond again.
    step = params["timestep"]
    obj = abstr_class()
    obj.ham_adi = CMATRIX(Ham_adi[step])
    obj.hvib_adi = CMATRIX(Hvib_adi[step])
    obj.nac_adi = CMATRIX(NAC_adi[step])
    obj.basis_transform = CMATRIX(Basis_transform[step])
    obj.time_overlap_adi = CMATRIX(St_adi[step])
    return obj


# ==================== 6. Dynamics settings and recipe selection ====================
# nsteps may exceed nfiles: Libra then repeats the stored trajectory. This is
# not new sampling; inspect the join between the last and first input frames.
# ntraj is the number of stochastic surface-hopping trajectories per icond.
# nfiles counts loaded frames; dt is DT converted from fs to atomic time.
# decoherence_rates and ave_gaps are the pairwise matrices calculated above.
# progress_frequency=0.1 reports roughly every 10% of the run (every 10 steps
# for a 100-step run). which_adi_states/which_dia_states choose state indices
# for supported state-resolved output; these do not replace basis selection.
# mem_output_level stores requested observables in memory and writes HDF5 at
# the end, avoiding repeated disk writes at the cost of memory use. BOTH the
# output level and properties_to_save must permit an observable to save it.
# Here we request time and populations. Energy, coordinate, amplitude, and
# Hamiltonian output need their appropriate saver levels and property names.
# prefix/prefix2 are assigned separately for each initial condition below.
dyn_general = {"nsteps": NSTEPS, "ntraj": 250, "nstates": NSTATES,
               "dt": DT * units.fs2au, "nfiles": NSTEPS,
               "decoherence_rates": rates, "ave_gaps": gaps,
               "Temperature": 300.0, "progress_frequency": 0.1,
               "which_adi_states": list(range(NSTATES)),
               "which_dia_states": list(range(NSTATES)),
               "mem_output_level": 2,
               "properties_to_save": ["timestep", "time", "se_pop_adi", "sh_pop_adi"],
               "isNBRA": 0, "is_nbra": 0}
# The two flag spellings cover old/new Libra drivers. Zero keeps separate
# Hamiltonians per electronic trajectory, as in the source tutorial.

# Select ONE recipe. Copy a file in recipes/ to make your own reusable recipe,
# edit its commented choices, then import it above and call its load() here.
# The notebook contrasts LD with NAC/Hvib propagation. Standard FSSH + LD uses
# tsh_method=0, decoherence_algo=-1, rep_tdse=1, electronic_integrator=2,
# state_tracking_algo=-1, and both separate phase corrections disabled.
# Its NAC-based example uses integrator=5, tracking=21, and ONE phase correction
# (do_phase_correction OR do_nac_phase_correction). Change these as a group in a
# copied recipe. FSSH2 + MSDM combines tsh_method=7, decoherence_algo=0 and
# decoherence_times_type=0 with the supplied rates; it is a separate recipe
# choice, not the default FSSH2 file's behavior.
#fssh_nbra.load(dyn_general); prf = "FSSH"
fssh_nbra_v2.load(dyn_general); prf = "FSSH"
#fssh_nbra_nac_v0.load(dyn_general); prf = "FSSH"
#fssh_nbra_nac_v1.load(dyn_general); prf = "FSSH"
#fssh_nbra_nac_v2.load(dyn_general); prf = "FSSH"
#fssh_nbra_nac_v3.load(dyn_general); prf = "FSSH"
#fssh_nbra_nac_v4.load(dyn_general); prf = "FSSH"
#fssh2_nbra.load(dyn_general); prf = "FSSH2"
#gfsh_nbra.load(dyn_general); prf = "GFSH"
#ida_nbra.load(dyn_general); prf = "IDA"
#msdm_nbra.load(dyn_general); prf = "MSDM"
#mash_nbra.load(dyn_general); prf = "MASH"
#dish_nbra.load(dyn_general); prf = "DISH"
#dish_rev2023_nbra.load(dyn_general); prf = "DISH_REV2023"

# Keep all four model entries: timestep is updated during dynamics, icond
# selects the starting geometry, model0 identifies the model for generic_recipe,
# and nstates gives the electronic dimension.
model_params = {"timestep": 0, "icond": 0, "model0": 0, "nstates": NSTATES}
# Dummy nuclear variables: NBRA reads the trajectory and computes no forces.
nucl_params = {"ndof": 1, "init_type": 3, "q": [-10.0], "p": [0.0],
               "mass": [2000.0], "force_constant": [0.01], "verbosity": -1}
# istate selects the initially occupied state. rep=1 initializes in the
# adiabatic representation (0 would be diabatic); rep_tdse in the recipe
# controls propagation. The remaining initialization choices set amplitudes
# and density matrices. ndia/nadi match the selected basis here.
elec_params = {"ndia": NSTATES, "nadi": NSTATES, "verbosity": -1,
               "init_dm_type": 0, "init_type": 1, "rep": 1, "istate": ISTATE}


# ==================== 7. Run initial conditions with a process pool ====================
# As in the notebook, one worker runs one initial-condition ensemble.
# nthreads is the number of worker PROCESSES, not electronic trajectories;
# ntraj above remains the number of trajectories per initial condition.
# Use no more workers than the CPUs allocated to this calculation.
# Each worker gets separate parameter dictionaries and a separate output folder.
# icond is an offset into the loaded files, not an absolute filename index.
# Set it in BOTH dictionaries: current Libra reads it from dyn_general.
# Compare pooled/sequential results statistically, not trajectory by trajectory.
ICONDS = [0, 5, 10]
#ICONDS = list(range(0, 25, 5))
nthreads = 4


def function1(icond):
    # Stagger random-generator creation as in the original tutorial.
    time.sleep(icond * 0.1)
    rnd = Random()
    mdl = dict(model_params)
    mdl.update({"icond": icond, "timestep": icond})
    dyn_gen = dict(dyn_general)
    prefix = f"{output_path}/{prf}_icond_{icond}"
    dyn_gen.update({"icond": icond, "prefix": prefix, "prefix2": prefix})
    tsh_dynamics.generic_recipe(dyn_gen, compute_model, mdl,
                                dict(elec_params), dict(nucl_params), rnd)
    # Do not return Libra objects: pool results would need to be pickled.


if __name__ == "__main__":
    # Explicit fork preserves the notebook's Linux workflow: workers inherit
    # the loaded CMATRIX objects without pickling or rerunning sections 1-6.
    # This version requires a platform supporting fork (intended for Linux).
    with mp.get_context("fork").Pool(nthreads) as pool:
        pool.map(function1, ICONDS)
