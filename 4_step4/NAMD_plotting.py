"""Plot saved NBRA dynamics and analyze electronic structure.

Run from 4_step4: python NAMD_plotting.py --prefix namd_run
Use the same --prefix as the simulation; figures are saved directly there.
Edit the input folder and settings below, and the plotting calls at the bottom.
The only command-line setting is --prefix, overriding the common run/output folder.
The selected calls plot basis diagnostics and saved dynamics, including TRPES
and Shannon entropy. For a basis-only run, disable the five post-dynamics calls.
This file reads saved data; it never imports or runs the simulation script.
Pre-NAMD electronic-structure plots are provided by ../3_step3/plotting.py.
Populations, fitted timescales, entropy, and TRPES require completed dynamics.
"""
import os
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from libra_py import units
from ld_active_space import diagnose_fixed_active_spaces


# ==================== 1. Settings to edit ====================
results_path = "namd_run"      # Same folder as the simulation --prefix
plot_path = results_path        # Save figures and fit tables in the run folder
methods = ["FSSH"]
#methods = ["FSSH", "IDA", "MSDM"]
ICONDS = [0, 5, 10]             # Use initial conditions that were simulated
ISTATE = 1                     # State whose decay/recovery will be fitted
population_type = "sh_pop_adi"
#population_type = "se_pop_adi"
BASIS_SIZE = None              # None: all loaded states; or a leading-state count
R2_MIN = 0.1                    # Illustrative fit-quality cutoff

# Optional time-energy map after TSH. This is a population-weighted PES, not
# a photoelectron spectrum (ionization energies/cross sections are not supplied).
TRPES_POPULATION = "sh_pop_adi"  # or "se_pop_adi"
TRPES_STATES = None             # None: all states; e.g. [1, 2, 3] excludes S0
TRPES_EMIN = -0.5               # eV relative to S0
TRPES_EMAX = 8.0
TRPES_DE = 0.02                 # energy-grid spacing, eV
TRPES_SIGMA = 0.10              # Gaussian standard deviation, eV

# Shannon entropy is computed separately for each IC before averaging.
ENTROPY_POPULATIONS = ["sh_pop_adi", "se_pop_adi"]
ENTROPY_NORMALIZE = True        # divide by ln(number of saved states)
TRPES_NORMALIZE_PEAK = True     # tutorial-style 0..1 display; archive keeps density
TRPES_CMAP = "magma"
TRPES_MAX_POINTS = None         # None: all saved times; or a positive sample count


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


# ==================== 6. Time-overlap diagnostics ====================
def plot_basis_quality():
    """Plot the saved SVD basis check, with support for older run archives."""
    size = NSTATES if BASIS_SIZE is None else BASIS_SIZE
    if not isinstance(size, int) or not 1 <= size <= NSTATES:
        raise ValueError("BASIS_SIZE must be None or a valid leading-state count")
    filename = f"{results_path}/basis_svd.npz"
    if os.path.exists(filename):
        with np.load(filename) as data:
            if (int(data["nframes"]) != NSTEPS or int(data["nstates"]) != NSTATES
                    or not np.isclose(float(data["dt_fs"]), DT)):
                raise ValueError("Saved SVD diagnostics do not match analysis_data.npz")
            sigma_min = data[f"sigma_min_{size}"]
            cutoff = data[f"cutoff_{size}"]
    else:
        # Legacy runs have no saved SVD arrays; only they need this calculation.
        res = diagnose_fixed_active_spaces(st, [size],
                                           relative_tolerance=1e-4,
                                           absolute_tolerance=1e-12)
        sigma_min, cutoff = res[size]["sigma_min"], res[size]["cutoff"]
    plt.figure(figsize=(8, 4))
    plt.plot(time, sigma_min, label="Minimum singular value")
    plt.plot(time, cutoff, "--", label="Rank cutoff")
    plt.xlabel("Time (fs)")
    plt.title(f"Active basis: {size} states")
    plt.legend()
    save_plot("overlap_singular_values")


def plot_overlap_element():
    i, j = 1, 2                 # Edit the pair for your chosen basis
    plt.figure(figsize=(8, 4))
    plt.plot(time, np.abs(st[:, i, j]))
    plt.xlabel("Time (fs)")
    plt.ylabel(f"|St[{i},{j}]|")
    save_plot("overlap_element")


# ==================== 7. Population-weighted time-energy map ====================
def population_energy_map(energy_ev, populations, time_fs, icond, dt_fs,
                          energy_grid, sigma, states=None):
    """Broaden populations at their IC-aligned energies; return density in 1/eV.

    Saved populations precede propagation: sample k uses input (icond+k)%nfiles.
    Saved times may be strided but must lie on the input-frame time grid.
    Accumulating one state at a time avoids a time x energy x state temporary.
    """
    energy_ev, populations = np.asarray(energy_ev), np.asarray(populations)
    time_fs, energy_grid = np.asarray(time_fs), np.asarray(energy_grid)
    if energy_ev.ndim != 2 or populations.ndim != 2 or energy_ev.shape[1] != populations.shape[1]:
        raise ValueError("Energy and population arrays must have matching state columns")
    if len(time_fs) != len(populations) or len(time_fs) == 0 or len(energy_ev) == 0:
        raise ValueError("Empty input or inconsistent time/population lengths")
    if not all(np.all(np.isfinite(x)) for x in (energy_ev, populations, time_fs, energy_grid)):
        raise ValueError("TRPES inputs must be finite")
    if np.any(populations < -1e-10):
        raise ValueError("Negative populations cannot define an energy density")
    if dt_fs <= 0 or sigma <= 0 or not np.isfinite(dt_fs + sigma):
        raise ValueError("DT and Gaussian sigma must be positive and finite")
    if energy_grid.ndim != 1 or len(energy_grid) < 2 or np.any(np.diff(energy_grid) <= 0):
        raise ValueError("Energy grid must be a strictly increasing vector")
    if time_fs.ndim != 1 or np.any(np.diff(time_fs) <= 0) or not np.isclose(time_fs[0], 0):
        raise ValueError("Saved time must start at zero and increase strictly")
    if not isinstance(icond, (int, np.integer)) or not 0 <= icond < len(energy_ev):
        raise ValueError("ICOND must be an offset into the loaded input frames")
    offsets = np.rint(time_fs / dt_fs).astype(int)
    # Libra's rounded au2fs and fs2au constants are not exact reciprocals.
    if not np.allclose(time_fs, offsets * dt_fs, rtol=1e-5, atol=1e-7):
        raise ValueError("Saved times do not match DT; check propagation and input frame spacing")
    selected = list(range(energy_ev.shape[1])) if states is None else list(states)
    if not selected or len(set(selected)) != len(selected) or any(
            not isinstance(s, (int, np.integer)) or s < 0 or s >= energy_ev.shape[1] for s in selected):
        raise ValueError("TRPES_STATES must contain unique, valid state indices")
    aligned = energy_ev[(icond + offsets) % len(energy_ev)]
    aligned = aligned - aligned[:, :1]
    density = np.zeros((len(time_fs), len(energy_grid)))
    for state in selected:
        gaussian = np.exp(-0.5 * ((energy_grid[None, :] - aligned[:, state, None]) / sigma)**2)
        density += np.maximum(populations[:, state, None], 0) * gaussian / (sigma * np.sqrt(2*np.pi))
    return density


def plot_trpes():
    """Average IC-resolved population-weighted energy maps after TSH."""
    if TRPES_POPULATION not in ("sh_pop_adi", "se_pop_adi"):
        raise ValueError("Choose sh_pop_adi or se_pop_adi for TRPES_POPULATION")
    if not ICONDS or len(set(ICONDS)) != len(ICONDS):
        raise ValueError("ICONDS must be nonempty and contain no duplicates")
    if not np.isfinite(TRPES_DE) or TRPES_DE <= 0 or TRPES_EMAX <= TRPES_EMIN:
        raise ValueError("Require TRPES_DE > 0 and TRPES_EMAX > TRPES_EMIN")
    grid = np.arange(TRPES_EMIN, TRPES_EMAX + TRPES_DE * 0.5, TRPES_DE)
    for method in methods:
        total = None
        for icond in ICONDS:
            filename = f"{results_path}/{method}_icond_{icond}/mem_data.hdf"
            with h5py.File(filename, "r") as data:
                t = np.asarray(data["time/data"]).reshape(-1) * units.au2fs
                pop = np.asarray(data[f"{TRPES_POPULATION}/data"])
            if TRPES_MAX_POINTS is not None:
                if not isinstance(TRPES_MAX_POINTS, int) or TRPES_MAX_POINTS < 1:
                    raise ValueError("TRPES_MAX_POINTS must be None or a positive integer")
                t, pop = t[:TRPES_MAX_POINTS], pop[:TRPES_MAX_POINTS]
            # Preserve the time origin to detect inconsistent/truncated runs.
            if total is None:
                reference_time = t
            elif t.shape != reference_time.shape or not np.allclose(t, reference_time):
                raise ValueError("TRPES requires matching saved time grids across ICs")
            density = population_energy_map(energies, pop, t, icond, DT,
                                             grid, TRPES_SIGMA, TRPES_STATES)
            total = density if total is None else total + density
        total /= len(ICONDS)
        area = np.trapz(total, grid, axis=1)
        mean_energy = np.divide(np.trapz(total * grid[None, :], grid, axis=1), area,
                                out=np.full_like(area, np.nan), where=area > 0)
        peak = total.max()
        display = total / peak if TRPES_NORMALIZE_PEAK and peak > 0 else total
        name = f"{method}_{TRPES_POPULATION}_trpes"
        np.savez(f"{plot_path}/{name}.npz", time_fs=reference_time, energy_ev=grid,
                 density_per_ev=total, iconds=ICONDS, sigma_ev=TRPES_SIGMA,
                 display_intensity=display, mean_energy_ev=mean_energy,
                 states=np.arange(NSTATES) if TRPES_STATES is None else TRPES_STATES)
        plt.figure(figsize=(9, 5))
        plt.pcolormesh(reference_time, grid, display.T, shading="auto", cmap=TRPES_CMAP)
        plt.colorbar(label="Relative intensity" if TRPES_NORMALIZE_PEAK else "Population density (1/eV)")
        plt.plot(reference_time, mean_energy, "w-", linewidth=1.5, label="Mean energy in window")
        plt.legend()
        plt.xlabel("Time after initialization (fs)")
        plt.ylabel("Energy relative to S0 (eV)")
        plt.title(f"{method}: population-weighted PES ({TRPES_POPULATION})")
        save_plot(name)


def shannon_entropy(populations, normalize=True):
    """Population Shannon entropy in nats, optionally divided by ln(N).

    Reject materially negative or unnormalized distributions. Correct only
    roundoff drift. This is basis-dependent population entropy, not the von
    Neumann entropy of the electronic density matrix.
    """
    p = np.asarray(populations, dtype=float)
    if p.ndim != 2 or min(p.shape) < 1 or not np.all(np.isfinite(p)):
        raise ValueError("Entropy requires a finite time-by-state population array")
    if np.any(p < -1e-10) or not np.allclose(p.sum(axis=1), 1, atol=1e-6, rtol=0):
        raise ValueError("Entropy requires nonnegative populations normalized to one")
    p = np.maximum(p, 0)
    p = p / p.sum(axis=1, keepdims=True)
    logp = np.zeros_like(p)
    np.log(p, out=logp, where=p > 0)
    entropy = -np.sum(p * logp, axis=1)
    if normalize:
        entropy = entropy / np.log(p.shape[1]) if p.shape[1] > 1 else np.zeros(len(p))
    return entropy


def plot_shannon_entropy():
    """Plot mean IC entropy and entropy of the mean population separately."""
    if not ICONDS or len(set(ICONDS)) != len(ICONDS):
        raise ValueError("ICONDS must contain unique initial conditions")
    for population in ENTROPY_POPULATIONS:
        if population not in ("sh_pop_adi", "se_pop_adi"):
            raise ValueError("Entropy supports sh_pop_adi and se_pop_adi")
        plt.figure(figsize=(9, 5))
        for method in methods:
            samples, pops = [], []
            for icond in ICONDS:
                with h5py.File(f"{results_path}/{method}_icond_{icond}/mem_data.hdf", "r") as data:
                    t = np.asarray(data["time/data"]).reshape(-1) * units.au2fs
                    pop = np.asarray(data[f"{population}/data"])
                if len(t) != len(pop) or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
                    raise ValueError("Entropy time and population rows must agree and times must increase")
                if pops and (t.shape != reference_time.shape or not np.allclose(t, reference_time)
                             or pop.shape != pops[0].shape):
                    raise ValueError("Entropy requires matching time grids and state counts across ICs")
                reference_time = t
                samples.append(shannon_entropy(pop, ENTROPY_NORMALIZE))
                pops.append(pop)
            samples = np.asarray(samples)
            mean = samples.mean(axis=0)
            pooled = shannon_entropy(np.mean(pops, axis=0), ENTROPY_NORMALIZE)
            name = f"{method}_{population}_shannon_entropy"
            np.savez(f"{plot_path}/{name}.npz", time_fs=reference_time, iconds=ICONDS,
                     entropy_per_ic=samples, mean_entropy=mean, entropy_of_mean_population=pooled,
                     normalized=ENTROPY_NORMALIZE, nstates=pops[0].shape[1])
            np.savetxt(f"{plot_path}/{name}.txt", np.column_stack((reference_time, mean, pooled)),
                       header="time_fs mean_IC_entropy entropy_of_mean_population")
            line, = plt.plot(reference_time, mean, label=f"{method}: mean IC entropy")
            plt.plot(reference_time, pooled, "--", color=line.get_color(),
                     label=f"{method}: entropy of mean population")
        plt.xlabel("Time (fs)")
        plt.ylabel("Shannon entropy / ln(N)" if ENTROPY_NORMALIZE else "Shannon entropy (nats)")
        plt.legend()
        save_plot(f"{population}_shannon_entropy")


# ==================== Select the plots to make ====================
def main():
    global results_path, plot_path, st, ham, nacs, DT, NSTEPS, NSTATES, time, energies
    parser = argparse.ArgumentParser(description="Plot saved NBRA results")
    parser.add_argument("--prefix", default=results_path,
                        help="Simulation results folder; figures are saved here too")
    args = parser.parse_args()
    results_path = plot_path = args.prefix
    os.makedirs(plot_path, exist_ok=True)

    # The simulation saves this archive during the initial basis-check run too.
    with np.load(f"{results_path}/analysis_data.npz") as data:
        st = data["overlaps"]
        ham = data["hamiltonians"]
        nacs = data["nacs"]
        DT = float(data["dt_fs"])
    NSTEPS, NSTATES = ham.shape[:2]
    time = np.arange(NSTEPS) * DT
    energies = np.diagonal(ham, axis1=1, axis2=2).real * units.au2ev

    plot_basis_quality()
    #plot_overlap_element()          # Choose valid i, j in the function first

    # Enable these after running sections 4-7 of NAMD_simulation.py.
    # Comment out all five calls for a basis-only analysis.
    plot_populations()
    plot_decoherence_times()
    plot_lifetime_fits()
    plot_trpes()                    # Set TRPES_*; requires saved TSH populations
    plot_shannon_entropy()          # Set ENTROPY_*


if __name__ == "__main__":
    main()
