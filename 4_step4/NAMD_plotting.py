"""Plots for NAMD_simulation.py, adapted from 9_step4/tutorial.ipynb.

Run from 4_step4: python NAMD_plotting.py --prefix namd_run
Use the same --prefix as the simulation; figures are saved directly there.
Edit the settings and uncomment the plotting calls at the bottom as needed.
The initial selection plots input properties and basis diagnostics. After
running dynamics, also uncomment the population, decoherence, and fit calls.
This file reads saved data; it never imports or runs the simulation script.
The notebook separates post-dynamics analysis (populations, decoherence maps,
and fitted timescales) from pre-dynamics descriptions (energies, gaps, couplings,
and overlaps). The latter require only the saved input matrices, so they are
useful even before a surface-hopping calculation has been performed.
"""
import os
import argparse

parser = argparse.ArgumentParser(description="Plot saved NBRA results")
parser.add_argument("--prefix", default="namd_run",
                    help="Simulation results folder; also saves figures here (default: namd_run)")
args = parser.parse_args()

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from libra_py import units
from ld_active_space import diagnose_fixed_active_spaces


# ==================== 1. Settings to edit ====================
results_path = args.prefix      # Same folder as the simulation --prefix
plot_path = results_path        # Save figures and fit tables in the run folder
methods = ["FSSH"]
#methods = ["FSSH", "IDA", "MSDM"]
ICONDS = [0, 5, 10]             # Use initial conditions that were simulated
ISTATE = 1                     # State whose decay/recovery will be fitted
population_type = "sh_pop_adi"
#population_type = "se_pop_adi"
R2_MIN = 0.1                    # Illustrative fit-quality cutoff
os.makedirs(plot_path, exist_ok=True)


def save_plot(name):
    """Save each figure; uncomment plt.show() to also open it on screen."""
    plt.tight_layout()
    plt.savefig(f"{plot_path}/{name}.png", dpi=200)
    #plt.show()
    plt.close()


def read_populations(method, icond):
    """HDF5 populations are already averaged over SH trajectories within an IC."""
    filename = f"{results_path}/{method}_icond_{icond}/mem_data.hdf"
    with h5py.File(filename, "r") as data:
        time = np.array(data["time/data"]).reshape(-1) * units.au2fs
        se_pop = np.array(data["se_pop_adi/data"])
        sh_pop = np.array(data["sh_pop_adi/data"])
    return time - time[0], se_pop, sh_pop


# ==================== 2. Populations ====================
def plot_populations():
    """Plot SE and SH populations separately, averaged over selected ICs.

    SE populations come from electronic amplitudes; SH populations count active
    surfaces. These estimators need not agree. Compare the same time grids and
    states in each run. Each initial condition receives equal weight here.

    The notebook uses Libra's tsh_dynamics_plot.plot_dynamics for a single IC,
    selecting states through which_adi_states. Here plain Matplotlib draws the
    same saved SE/SH observables. Set ICONDS=[0] for one IC; edit the state loop
    below to display a subset. The default averages the selected ICs.
    Its sequential/pool folders have different prefixes (FSSH_NBRA versus FSSH);
    set methods/results_path to match the folders being compared.
    """
    for method in methods:
        all_se, all_sh = [], []
        for icond in ICONDS:
            time, se_pop, sh_pop = read_populations(method, icond)
            if len(all_se) == 0:
                reference_time = time
            assert time.shape == reference_time.shape and np.allclose(time, reference_time)
            all_se.append(se_pop)
            all_sh.append(sh_pop)
        mean_se = np.mean(all_se, axis=0)
        mean_sh = np.mean(all_sh, axis=0)
        for label, populations in [("SE", mean_se), ("SH", mean_sh)]:
            plt.figure(figsize=(8, 5))
            for state in range(populations.shape[1]):
                plt.plot(time, populations[:, state], label=f"S{state}")
            plt.xlabel("Time (fs)")
            plt.ylabel("Population")
            plt.title(f"{method}: {label} populations, averaged over initial conditions")
            plt.legend()
            save_plot(f"{method}_{label}_populations")


# ==================== 3. Decoherence time map ====================
def plot_decoherence_times():
    """Show a state-pair map of average decoherence times in fs.

    Each pixel corresponds to a pair of electronic states. The displayed zero
    diagonal is a plotting convention for absent self-decoherence, not a zero
    physical lifetime. A chosen color limit (e.g. 30 fs in the notebook) can
    reveal short-time structure but saturates longer times.
    """
    times = np.loadtxt(f"{results_path}/decoherence_times.txt", ndmin=2)
    plt.figure(figsize=(7, 6))
    plt.imshow(times, origin="lower", cmap="hot")
    #plt.clim(0, 30)
    plt.colorbar(label="fs")
    plt.xlabel("State index")
    plt.ylabel("State index")
    plt.title("Average decoherence times")
    save_plot("decoherence_times")


# ==================== 4. Fit population decay or recovery ====================
def exp_funct(time, tau, beta):
    # beta=1 is exponential decay; beta=2 is the notebook's Gaussian example.
    # Use the recovery expression for a state starting at zero and tending to 1.
    return np.exp(-(time / tau)**beta)
    #return 1.0 - np.exp(-(time / tau)**beta)


def plot_lifetime_fits():
    """Fit each IC, retain R^2 > R2_MIN, and summarize tau and beta.

    tau is the 1/e decay time, not generally the integrated lifetime. The model
    assumes unit initial population and zero final population (reversed for
    recovery). Choose a different model if the population has a plateau.

    The notebook's R^2 > 0.1 is a teaching example, not a convergence criterion.
    The approximate 95% error on the mean is 1.96*s/sqrt(N), using sample standard
    deviation. It assumes independent ICs; neighboring time origins may be
    correlated. A single accepted fit cannot give this error estimate.
    R^2 = 1 - sum((P-Pfit)^2)/sum((P-mean(P))^2); constant traces have no
    meaningful R^2 and are skipped. N counts accepted IC fits, not time points
    or individual hopping trajectories. Larger ntraj and longer dynamics are
    needed to assess whether these illustrative timescales are reliable.
    The notebook prose uses beta=2; its code fits beta as we do here.
    """
    for method in methods:
        taus, betas, fit_rows = [], [], []
        plt.figure(figsize=(9, 6))
        for icond in ICONDS:
            time, se_pop, sh_pop = read_populations(method, icond)
            populations = sh_pop
            if population_type == "se_pop_adi":
                populations = se_pop
            population = populations[:, ISTATE]
            plt.plot(time, population, alpha=0.5, label=f"IC {icond}")
            if np.ptp(population) < 1e-12 or len(time) < 4:
                print(method, icond, "constant or too short: skipping fit")
                continue
            try:
                popt, pcov = curve_fit(exp_funct, time, population,
                                      p0=[100.0, 1.0],
                                      bounds=([1e-12, 0.05], [np.inf, 5.0]),
                                      maxfev=10000)
            except (RuntimeError, ValueError) as error:
                print(method, icond, "fit failed:", error)
                continue
            tau, beta = popt
            fitted = exp_funct(time, tau, beta)
            r2 = 1.0 - np.sum((population - fitted)**2) / np.sum((population - population.mean())**2)
            print(f"{method}, IC {icond}: tau={tau:.3f} fs, beta={beta:.3f}, R2={r2:.3f}")
            fit_rows.append([icond, tau, beta, r2])
            plt.plot(time, fitted, "--", linewidth=1)
            if r2 > R2_MIN:
                taus.append(tau)
                betas.append(beta)
        np.savetxt(f"{plot_path}/{method}_fits.txt", np.array(fit_rows).reshape(-1, 4),
                   header="icond tau_fs beta R_squared")
        if len(taus) > 0:
            print(f"{method}: {len(taus)} accepted fits")
            print(f"Mean tau = {np.mean(taus)/1000:.6f} ps; mean beta = {np.mean(betas):.3f}")
            if len(taus) > 1:
                error_tau = 1.96 * np.std(taus, ddof=1) / np.sqrt(len(taus))
                error_beta = 1.96 * np.std(betas, ddof=1) / np.sqrt(len(betas))
                print(f"95% half-widths: tau = {error_tau/1000:.6f} ps; beta = {error_beta:.3f}")
        else:
            print(method, "no fits passed the chosen R2 cutoff")
        plt.xlabel("Time (fs)")
        plt.ylabel(f"S{ISTATE} population ({population_type})")
        plt.title(f"{method}: individual population fits")
        plt.legend()
        save_plot(f"{method}_lifetime_fits")


# ==================== 5. Input properties before dynamics ====================
def plot_excitation_energies():
    """Follow each adiabatic energy along the sampled nuclear trajectory.

    Subtract S0 to display excitation energies even if the input has an energy
    offset. Time uses the actual DT; energies are converted from Hartree to eV.
    """
    plt.figure(figsize=(8, 5))
    for state in range(NSTATES):
        plt.plot(time, energies[:, state] - energies[:, 0], label=f"S{state}")
    plt.xlabel("Time (fs)")
    plt.ylabel("Excitation energy relative to S0 (eV)")
    plt.legend()
    save_plot("excitation_energies")


def plot_gap_distributions():
    """Adjacent-state gaps keep the state ordering of the input matrices.

    Histograms describe fluctuations of E[i+1]-E[i] along the trajectory.
    Density normalization gives unit area, allowing distributions with different
    widths to be compared; this is not a histogram of electronic populations.
    """
    plt.figure(figsize=(8, 5))
    for state in range(NSTATES - 1):
        gap = energies[:, state + 1] - energies[:, state]
        plt.hist(gap, bins=50, density=True, histtype="step", label=f"S{state+1}-S{state}")
    plt.xlabel("Adjacent energy gap (eV)")
    plt.ylabel("Probability density (1/eV)")
    plt.legend()
    save_plot("gap_distributions")


def plot_nac_map():
    """Compare mean absolute couplings for different pairs of states.

    Average |NAC|, not the signed NAC, to avoid cancellation of fluctuating
    couplings. The optional diagonal band changes the displayed map only.
    """
    # In atomic units hbar=1; these time-derivative couplings are shown in meV.
    # They are not spatial derivative coupling vectors.
    average = np.mean(np.abs(nacs), axis=0) * units.au2ev * 1000.0
    np.fill_diagonal(average, 0)
    # To show only a band around the diagonal, uncomment:
    #band = 8
    #for i in range(NSTATES):
    #    for j in range(NSTATES):
    #        if abs(i - j) > band:
    #            average[i, j] = 0
    plt.figure(figsize=(7, 6))
    plt.imshow(average, origin="lower", cmap="hot")
    plt.colorbar(label="meV")
    plt.xlabel("State index")
    plt.ylabel("State index")
    plt.title("Mean absolute nonadiabatic couplings")
    save_plot("nac_map")


def plot_nac_distribution(pairs, name):
    """Unique-pair |NAC| distribution; omit exact zeros for logarithmic axes.

    All pairs use i<j to avoid counting a pair twice; adjacent pairs use j=i+1.
    Comparing the two reveals whether neighboring states represent the broader
    coupling distribution. Log-spaced bins show its range across orders of
    magnitude. Density is still per meV, not per logarithmic unit.
    """
    values = [np.abs(nacs[:, i, j]) * units.au2ev * 1000.0 for i, j in pairs]
    if len(values) == 0:
        return
    values = np.concatenate(values)
    values = values[values > 0]
    if len(values) == 0:
        print(name, "has no nonzero couplings")
        return
    low, high = values.min(), values.max()
    if low == high:
        low, high = low / 1.1, high * 1.1
    bins = np.geomspace(low, high, 80)
    hist, edges = np.histogram(values, bins=bins, density=True)
    centers = np.sqrt(edges[:-1] * edges[1:])
    keep = hist > 0
    plt.figure(figsize=(8, 5))
    plt.loglog(centers[keep], hist[keep])
    plt.xlabel("|NAC| (meV)")
    plt.ylabel("Probability density (1/meV)")
    plt.title(name)
    save_plot(name)


# ==================== 6. Time-overlap diagnostics ====================
def plot_basis_quality():
    """Plot a selected size; edit nact_st for a different leading active block.

    The determinant time trace locates intervals where the basis loses overlap
    support. Compare sizes without assuming larger determinants alone establish
    completeness. The companion singular-value/cutoff plot makes the numerical
    rank criterion explicit; use the simulation's table to compare all sizes.
    """
    nact_st = NSTATES
    #nact_st = 3
    overlaps = st[:, :nact_st, :nact_st]
    res = diagnose_fixed_active_spaces(st, [nact_st],
                                       relative_tolerance=1e-4,
                                       absolute_tolerance=1e-12)
    plt.figure(figsize=(8, 4))
    plt.plot(time, np.abs(np.linalg.det(overlaps)))
    plt.xlabel("Time (fs)")
    plt.ylabel("|det(St)|")
    plt.title(f"Active basis: {nact_st} states")
    save_plot("overlap_determinants")
    plt.figure(figsize=(8, 4))
    plt.plot(time, res[nact_st]["sigma_min"], label="Minimum singular value")
    plt.plot(time, res[nact_st]["cutoff"], "--", label="Rank cutoff")
    plt.xlabel("Time (fs)")
    plt.legend()
    save_plot("overlap_singular_values")


def plot_overlap_element():
    i, j = 1, 2                 # Edit the pair for your chosen basis
    plt.figure(figsize=(8, 4))
    plt.plot(time, np.abs(st[:, i, j]))
    plt.xlabel("Time (fs)")
    plt.ylabel(f"|St[{i},{j}]|")
    save_plot("overlap_element")


# ==================== 7. Select the plots to make ====================
# The simulation saves this archive during the initial basis-check run too.
with np.load(f"{results_path}/analysis_data.npz") as data:
    st = data["overlaps"]
    ham = data["hamiltonians"]
    nacs = data["nacs"]
    DT = float(data["dt_fs"])
NSTEPS, NSTATES = ham.shape[:2]
time = np.arange(NSTEPS) * DT
energies = np.diagonal(ham, axis1=1, axis2=2).real * units.au2ev

plot_excitation_energies()
plot_gap_distributions()
plot_nac_map()
all_pairs = [(i, j) for i in range(NSTATES) for j in range(i + 1, NSTATES)]
adjacent_pairs = [(i, i + 1) for i in range(NSTATES - 1)]
plot_nac_distribution(all_pairs, "nac_distribution_all_pairs")
plot_nac_distribution(adjacent_pairs, "nac_distribution_adjacent_pairs")
plot_basis_quality()
#plot_overlap_element()          # Choose valid i, j in the function first

# Uncomment these after running sections 4-7 of NAMD_simulation.py:
plot_populations()
plot_decoherence_times()
plot_lifetime_fits()
