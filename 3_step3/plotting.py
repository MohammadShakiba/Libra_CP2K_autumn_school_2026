"""Plot pre-NAMD electronic structure after run.py; never runs CP2K or dynamics.

Run from 3_step3: python plotting.py
Edit input paths, frame selection, and PLOTS below before running.
The only command-line setting is --prefix for the output directory.
"""
import argparse
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from libra_py import units

# ==================== Settings to edit before running ====================
DATA_DIRECTORY = "ci_tio2"
OUTPUT_DIRECTORY = "plots"
ISTEP, FSTEP = 0, 20             # end excluded
DT = 1.0                       # frame spacing, fs (matrices and logs)
# Comment out entries to disable analyses; no CP2K logs are needed unless
# "excitation" or "uvvis" is selected.
PLOTS = [
    "energies",
    "gaps",
    "nacs",
    "overlaps",
    "spectra",
    "excitation",
    "uvvis",
]
# True for the aligned run.py recipe: first output is cache initialization.
SKIP_INITIAL_NAC = True
NAC_PAIRS = "all"              # "all" unique pairs or "adjacent"; one distribution
UVVIS_FWHM_EV = 0.10            # full width at half maximum, eV
UVVIS_GRID_STEP_EV = 0.005
UVVIS_STATES = None             # None: all EXC_NSTATES; or e.g. [1, 2, 3]
                               # excited-state numbers, starting at 1
DOS_BINS = 50

# Gap spectra use the saved Step 3 Hamiltonians.
SPECTRUM_PAIRS = [(0, 1)]
SPECTRUM_TEMPERATURE = 300.0     # K, for classical-to-harmonic bath conversion
SPECTRUM_MAX_WAVENUMBER = 4000.0 # cm^-1, display limit only

# CP2K log settings: select the same trajectory and state set as Step 3.
EXC_LOG_DIRECTORY = "all_logfiles"        # relative to the working directory
EXC_LOG_PATTERN = "step_{step}.log"
EXC_NSTATES = 10                # excited states only, excludes S0
EXC_CI_THRESHOLD = 0.001        # squared-amplitude cutoff
EXC_IS_UKS = False
EXC_NSDS = 4                    # leading entries/ranks to display
EXC_TOP_TRANSITIONS = 12        # most important orbital identities per state


def save_plot(name):
    plt.tight_layout()
    plt.savefig(f"{plot_path}/{name}.png", dpi=200)
    plt.close()


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


def gap_spectrum(gap_ev, dt_fs, temperature):
    """Biased centered ACF and its two-sided Fourier spectrum at omega >= 0.

    Compute in atomic units. I(omega)=integral C(t) exp(i*omega*t) dt;
    J(omega)=beta*omega*I(omega)/2 is the chosen classical harmonic-bath
    convention. Zero padding samples the transform more finely, not more accurately.
    """
    gap = np.asarray(gap_ev, dtype=float) * units.ev2Ha
    if gap.ndim != 1 or len(gap) < 3 or not np.all(np.isfinite(gap)):
        raise ValueError("A spectrum needs at least three finite gap samples")
    if dt_fs <= 0 or temperature <= 0 or not np.isfinite(dt_fs + temperature):
        raise ValueError("DT and spectrum temperature must be finite and positive")
    fluctuation = gap - gap.mean()
    n, dt_au = len(gap), dt_fs * units.fs2au
    transform = np.fft.rfft(fluctuation, n=2*n)
    power = np.abs(transform)**2 / n
    acf = np.fft.irfft(power, n=2*n)[:n]  # sum(deltaE[t]*deltaE[t+k])/N
    omega = 2*np.pi*np.fft.rfftfreq(2*n, dt_au)
    influence = dt_au * power  # two-sided PSD evaluated on nonnegative frequencies
    spectral_density = omega * influence / (2 * units.kB * temperature)
    return acf, omega, influence, spectral_density


def plot_gap_spectra():
    """Plot centered gaps, ACFs, influence spectra, and harmonic spectral density."""
    for i, j in SPECTRUM_PAIRS:
        if i == j or min(i, j) < 0 or max(i, j) >= NSTATES:
            raise ValueError(f"Invalid spectrum state pair {(i, j)}")
        gap = energies[:, j] - energies[:, i]
        acf, omega, influence, spectral = gap_spectrum(gap, DT, SPECTRUM_TEMPERATURE)
        frequency = omega * units.au2wavn
        raw_acf = acf / units.ev2Ha**2
        noise = influence / units.ev2Ha**2 / units.fs2au
        norm_acf = acf / acf[0] if acf[0] > 1e-30 else np.zeros_like(acf)
        name = f"gap_S{i}_S{j}"
        np.savez(f"{plot_path}/{name}_spectra.npz", time_fs=time, gap_ev=gap,
                 fluctuation_ev=gap-gap.mean(), acf_ev2=raw_acf, normalized_acf=norm_acf,
                 wavenumber_cm=frequency, influence_ev2_fs=noise,
                 spectral_density_ha=spectral, temperature_k=SPECTRUM_TEMPERATURE)
        fig, axes = plt.subplots(2, 2, figsize=(11, 7))
        axes[0, 0].plot(time, gap-gap.mean())
        axes[0, 0].set(xlabel="Time (fs)", ylabel="Gap fluctuation (eV)")
        axes[0, 1].plot(time, norm_acf)
        axes[0, 1].set(xlabel="Lag (fs)", ylabel="Normalized gap ACF")
        axes[1, 0].plot(frequency, noise)
        axes[1, 0].set(xlabel="Wavenumber (cm$^{-1}$)", ylabel="Influence spectrum (eV$^2$ fs)")
        axes[1, 1].plot(frequency, spectral / units.ev2Ha)
        axes[1, 1].set(xlabel="Wavenumber (cm$^{-1}$)", ylabel="Spectral density (eV)")
        for ax in axes[1]:
            ax.set_xlim(0, min(SPECTRUM_MAX_WAVENUMBER, frequency[-1]))
        fig.suptitle(f"S{i}–S{j}: gap fluctuations, T={SPECTRUM_TEMPERATURE:g} K")
        save_plot(f"{name}_spectra")


def excitation_statistics(records, nstates, nsds):
    """Track printed entries, ranked weights, and fixed (occ,vir,spin) identities.

    Normalize squared amplitudes within each retained expansion, as in the
    excitation tutorial. Also retain raw squared norms to expose truncation.
    """
    if not records or nstates < 1 or nsds < 1:
        raise ValueError("Need records and positive EXC_NSTATES/EXC_NSDS")
    shape = (len(records), nstates)
    raw_weight, count, participation = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    printed = np.zeros(shape + (nsds,))
    ranked = np.zeros_like(printed)
    transitions = [[{} for _ in range(nstates)] for _ in records]
    for frame, data in enumerate(records):
        if any(len(data[k]) < nstates for k in (1, 2, 3)):
            raise ValueError(f"Log sample {frame} has fewer than {nstates} excited states")
        for state in range(nstates):
            configs, coeffs, spins = data[1][state], data[2][state], data[3][state]
            if not (len(configs) == len(coeffs) == len(spins)):
                raise ValueError("Excitation configurations, amplitudes, and spin labels must align")
            weights = np.asarray(coeffs, dtype=float)**2
            norm = weights.sum()
            if not np.isfinite(norm) or norm <= 0:
                raise ValueError(f"Empty/invalid excitation in sample {frame}, state S{state+1}; lower cutoff")
            raw_weight[frame, state], count[frame, state] = norm, len(weights)
            weights /= norm
            participation[frame, state] = 1 / np.sum(weights**2)
            k = min(nsds, len(weights))
            printed[frame, state, :k] = weights[:k]
            ranked[frame, state, :k] = np.sort(weights)[::-1][:k]
            for (occ, vir), spin, weight in zip(configs, spins, weights):
                key = (int(occ), int(vir), str(spin))
                transitions[frame][state][key] = transitions[frame][state].get(key, 0) + weight
    return raw_weight, count, participation, printed, ranked, transitions


def plot_excitation_analysis():
    """Read only the configured CP2K logs; no dynamics or CP2K jobs are run."""
    from pathlib import Path
    import csv
    import libra_py.packages.cp2k.methods as cp2k

    if EXC_LOG_DIRECTORY is None:
        raise ValueError("Set EXC_LOG_DIRECTORY to the CP2K logs for the selected system")
    if EXC_FSTEP <= EXC_ISTEP or EXC_DT_FS <= 0 or EXC_CI_THRESHOLD < 0 or EXC_TOP_TRANSITIONS < 1:
        raise ValueError("Check excitation frame range, timestep, threshold, and transition count")
    records = []
    for step in range(EXC_ISTEP, EXC_FSTEP):
        filename = Path(EXC_LOG_DIRECTORY) / EXC_LOG_PATTERN.format(step=step)
        _, data = cp2k.read_cp2k_tddfpt_log_file({
            "logfile_name": str(filename), "number_of_states": EXC_NSTATES,
            "ci_threshold": EXC_CI_THRESHOLD, "isUKS": EXC_IS_UKS})
        records.append(data)
    raw, count, participation, printed, ranked, transitions = excitation_statistics(
        records, EXC_NSTATES, EXC_NSDS)
    t = np.arange(len(records)) * EXC_DT_FS
    np.savez(f"{plot_path}/excitation_analysis.npz", time_fs=t,
             file_steps=np.arange(EXC_ISTEP, EXC_FSTEP), raw_squared_norm=raw,
             retained_count=count, participation_ratio=participation,
             printed_weights=printed, ranked_weights=ranked, ci_threshold=EXC_CI_THRESHOLD)
    for name, values in [("printed", printed), ("ranked", ranked)]:
        mean = values.mean(axis=0)
        # Undefined uncertainty for one frame: omit bars rather than claim zero error.
        error = (1.96*values.std(axis=0, ddof=1)/np.sqrt(len(records))
                 if len(records) > 1 else np.full_like(mean, np.nan))
        rows = [[s+1, k+1, mean[s, k], error[s, k]]
                for s in range(EXC_NSTATES) for k in range(EXC_NSDS)]
        np.savetxt(f"{plot_path}/excitation_{name}_weights.txt", rows,
                   header="state entry_or_rank mean_normalized_weight naive_95pct_halfwidth")
        plt.figure(figsize=(9, 5))
        for k in range(EXC_NSDS):
            plt.errorbar(np.arange(1, EXC_NSTATES+1) + 0.05*k, mean[:, k],
                         yerr=error[:, k] if len(records) > 1 else None,
                         fmt="o", capsize=2, label=f"{name} {k+1}")
        plt.xlabel("Excited state (S1, S2, …)")
        plt.ylabel("Mean normalized squared coefficient")
        plt.legend()
        save_plot(f"excitation_{name}_weights")
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for state in range(EXC_NSTATES):
        for ax, values in zip(axes, (raw, count, participation)):
            ax.plot(t, values[:, state], label=f"S{state+1}")
    for ax, label in zip(axes, ("Retained squared norm", "Retained configurations", "Participation ratio")):
        ax.set_ylabel(label)
    axes[0].legend(ncol=min(EXC_NSTATES, 5))
    axes[-1].set_xlabel("Time from first selected log (fs)")
    save_plot("excitation_mixing")
    with open(f"{plot_path}/excitation_transitions.csv", "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["state", "occupied", "virtual", "spin", "mean_normalized_weight"])
        for state in range(EXC_NSTATES):
            identities = sorted(set().union(*(frame[state] for frame in transitions)))
            weights = np.array([[frame[state].get(key, 0) for key in identities] for frame in transitions])
            mean = weights.mean(axis=0)
            for key, weight in zip(identities, mean):
                writer.writerow([state+1, *key, weight])
            top = np.argsort(mean)[::-1][:EXC_TOP_TRANSITIONS]
            plt.figure(figsize=(10, 5))
            plt.pcolormesh(t, np.arange(len(top)), weights[:, top].T,
                           shading="auto", vmin=0, vmax=1, cmap="viridis")
            plt.yticks(np.arange(len(top)), [f"{identities[k][0]}→{identities[k][1]} ({identities[k][2]})" for k in top])
            plt.colorbar(label="Normalized squared coefficient")
            plt.xlabel("Time from first selected log (fs)")
            plt.title(f"S{state+1}: leading orbital transitions")
            save_plot(f"excitation_S{state+1}_transitions")



def read_step3_matrix(filename):
    """Read real or Libra (real,imag) text entries without running Step 3."""
    tokens = np.loadtxt(filename, dtype=str, ndmin=2)
    def value(token):
        if token.startswith("(") and token.endswith(")"):
            real, imag = token[1:-1].split(",")
            return complex(float(real), float(imag))
        return complex(token)
    matrix = np.array([[value(x) for x in row] for row in tokens])
    if matrix.shape[0] != matrix.shape[1] or not np.all(np.isfinite(matrix)):
        raise ValueError(f"Expected a finite square matrix: {filename}")
    return matrix


def plot_overlap_diagnostics():
    """Check same-frame CI orthonormality; time-overlap SVD belongs to Step 4."""
    deviation = np.max(np.abs(same_frame - np.eye(NSTATES)), axis=(1, 2))
    plt.figure(figsize=(9, 4))
    plt.plot(time, deviation)
    plt.ylabel("max |S - identity|")
    plt.xlabel("Time from first selected output (fs)")
    save_plot("overlap_orthonormality")
    np.savez(f"{plot_path}/overlap_orthonormality.npz", time_fs=time,
             orthonormality_error=deviation)


def read_optical_transitions(filename, nstates):
    """Read CP2K's eV/dipole/oscillator-strength TDDFPT table.

    Require the documented column headers; reject repeated state tables rather
    than silently combine different multiplicities or multiple calculations.
    CI coefficient squares are not oscillator strengths.
    """
    text = Path(filename).read_text()
    if "Oscillator" not in text or "energy (eV)" not in text:
        raise ValueError(f"No supported eV/oscillator-strength table in {filename}")
    rows = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) != 7 or parts[0] != "TDDFPT|" or not parts[1].isdigit():
            continue
        state = int(parts[1])
        if not 1 <= state <= nstates:
            continue
        if state in rows:
            raise ValueError(f"Multiple TDDFPT tables in {filename}; select one calculation/multiplicity")
        energy, oscillator = (float(parts[k].replace("D", "E").replace("d", "e")) for k in (2, 6))
        if not np.isfinite(energy + oscillator) or energy <= 0 or oscillator < 0:
            raise ValueError(f"Invalid optical transition S{state} in {filename}")
        rows[state] = (energy, oscillator)
    if set(rows) != set(range(1, nstates+1)):
        raise ValueError(f"Expected optical transitions S1..S{nstates} in {filename}; found {sorted(rows)}")
    return np.array([rows[i] for i in range(1, nstates+1)])


def broaden_optical_spectrum(excitation_ev, oscillator_strengths, grid_ev, fwhm_ev):
    """Average oscillator-strength density over frames on a common energy grid."""
    energy, strength = np.asarray(excitation_ev), np.asarray(oscillator_strengths)
    grid = np.asarray(grid_ev)
    if energy.ndim != 2 or min(energy.shape) < 1 or energy.shape != strength.shape:
        raise ValueError("Optical arrays must share a nonempty frame-by-state shape")
    if (not all(np.all(np.isfinite(x)) for x in (energy, strength, grid))
            or np.any(strength < 0) or not np.isfinite(fwhm_ev) or fwhm_ev <= 0):
        raise ValueError("Require finite data, nonnegative oscillator strengths, and positive FWHM")
    if grid.ndim != 1 or len(grid) < 2 or np.any(np.diff(grid) <= 0):
        raise ValueError("Optical energy grid must be strictly increasing")
    sigma = fwhm_ev / np.sqrt(8*np.log(2))
    density = np.zeros(len(grid))
    # Sum one frame at a time to keep memory independent of trajectory length.
    for e, f in zip(energy, strength):
        gaussian = np.exp(-0.5*((grid[:, None]-e[None, :])/sigma)**2)
        density += gaussian @ f / (sigma*np.sqrt(2*np.pi))
    return density / len(energy)


def plot_energy_uvvis_dos():
    """Energy traces and side-rotated excitation DOS/UV-vis with shared energy axis."""
    tables = np.stack([read_optical_transitions(
        Path(EXC_LOG_DIRECTORY) / EXC_LOG_PATTERN.format(step=step), EXC_NSTATES)
        for step in range(EXC_ISTEP, EXC_FSTEP)])
    states = list(range(1, EXC_NSTATES+1)) if UVVIS_STATES is None else list(UVVIS_STATES)
    if (not states or len(set(states)) != len(states)
            or any(not isinstance(s, int) or not 1 <= s <= EXC_NSTATES for s in states)):
        raise ValueError("UVVIS_STATES must be unique excited-state numbers from 1 to EXC_NSTATES")
    e = tables[:, np.array(states)-1, 0]
    f = tables[:, np.array(states)-1, 1]
    if (UVVIS_FWHM_EV <= 0 or UVVIS_GRID_STEP_EV <= 0
            or not np.isfinite(UVVIS_FWHM_EV + UVVIS_GRID_STEP_EV)):
        raise ValueError("UV-vis FWHM and energy-grid step must be finite and positive")
    sigma = UVVIS_FWHM_EV / np.sqrt(8*np.log(2))
    grid = np.arange(max(0, e.min()-5*sigma), e.max()+5*sigma+UVVIS_GRID_STEP_EV/2,
                     UVVIS_GRID_STEP_EV)
    spectrum = broaden_optical_spectrum(e, f, grid, UVVIS_FWHM_EV)
    # DOS counts every included excited state equally; normalized to unit area.
    dos, edges = np.histogram(e.ravel(), bins=DOS_BINS, range=(grid[0], grid[-1]), density=True)
    t = np.arange(len(tables)) * EXC_DT_FS
    np.savez(f"{plot_path}/energy_uvvis_dos.npz", time_fs=t, states=states,
             excitation_ev=e, oscillator_strengths=f, energy_grid_ev=grid,
             oscillator_strength_density_per_ev=spectrum, fwhm_ev=UVVIS_FWHM_EV,
             dos_per_ev=dos, dos_edges_ev=edges,
             file_steps=np.arange(EXC_ISTEP, EXC_FSTEP))
    np.savetxt(f"{plot_path}/uvvis_spectrum.txt", np.column_stack((grid, spectrum)),
               header="energy_ev mean_oscillator_strength_density_per_ev")
    np.savetxt(f"{plot_path}/excitation_dos.txt", np.column_stack((edges[:-1], edges[1:], dos)),
               header="energy_lower_ev energy_upper_ev probability_density_per_ev")
    fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharey=True,
                             gridspec_kw={"width_ratios": [4, 1.2, 1.2]})
    for col, state in enumerate(states):
        axes[0].plot(t, e[:, col], linewidth=1, label=f"S{state}")
    axes[0].set(xlabel="Time from first selected log (fs)", ylabel="Excitation energy (eV)")
    if len(states) <= 12:
        axes[0].legend(ncol=2, fontsize="small")
    axes[1].barh((edges[:-1]+edges[1:])/2, dos, height=np.diff(edges), color="gray", alpha=.7)
    axes[1].set(xlabel="DOS (1/eV)", title="Excited-state DOS")
    axes[2].plot(spectrum, grid, color="tab:blue")
    axes[2].set(xlabel="Strength / eV", title="UV–vis")
    axes[0].set_ylim(grid[0], grid[-1])
    fig.suptitle("Excitation energies, density of states, and optical absorption")
    save_plot("energy_uvvis_dos")
    plt.figure(figsize=(8, 4))
    plt.plot(grid, spectrum)
    plt.xlabel("Excitation energy (eV)")
    plt.ylabel("Mean oscillator-strength density (1/eV)")
    save_plot("uvvis_spectrum")


def main():
    global plot_path, EXC_LOG_DIRECTORY, EXC_ISTEP, EXC_FSTEP, EXC_DT_FS
    global st, same_frame, ham, nacs, DT, NSTEPS, NSTATES, time, energies
    parser = argparse.ArgumentParser(description="Plot Step 3 electronic structure after run.py")
    parser.add_argument("--prefix", default=OUTPUT_DIRECTORY, help="Output plot/data directory")
    args = parser.parse_args()
    if FSTEP <= ISTEP or not np.isfinite(DT) or DT <= 0:
        raise ValueError("Set FSTEP > ISTEP and finite DT > 0 in the script")
    plot_path = args.prefix
    # Use the same editable frame range and spacing for matrices and logs.
    EXC_ISTEP, EXC_FSTEP, EXC_DT_FS = ISTEP, FSTEP, DT
    plots = set(PLOTS)
    available = {"energies", "gaps", "nacs", "overlaps", "spectra", "excitation", "uvvis"}
    if not plots or plots - available:
        raise ValueError(f"PLOTS must select from {sorted(available)}")
    os.makedirs(plot_path, exist_ok=True)
    if plots & {"energies", "gaps", "nacs", "overlaps", "spectra"}:
        ham = np.stack([read_step3_matrix(Path(DATA_DIRECTORY)/f"ham_adi_{i}.txt")
                        for i in range(ISTEP, FSTEP)])
        NSTEPS, NSTATES = ham.shape[:2]
        time = np.arange(NSTEPS)*DT
        energies = np.diagonal(ham, axis1=1, axis2=2).real * units.au2ev
    if "nacs" in plots:
        st = np.stack([read_step3_matrix(Path(DATA_DIRECTORY)/f"st_adi_{i}.txt")
                       for i in range(ISTEP, FSTEP)])
        if st.shape != ham.shape:
            raise ValueError("Hamiltonian and time-overlap matrix shapes must match")
        nacs = (st - st.conj().transpose(0, 2, 1)) / (2*DT*units.fs2au)
        if SKIP_INITIAL_NAC and ISTEP == 0:
            nacs = nacs[1:]
        if "nacs" in plots and len(nacs) == 0:
            raise ValueError("No physical intervals remain for NAC analysis")
    if "energies" in plots:
        plot_excitation_energies()
    if "gaps" in plots:
        plot_gap_distributions()
    if "nacs" in plots:
        plot_nac_map()
        if NAC_PAIRS == "all":
            pairs = [(i, j) for i in range(NSTATES) for j in range(i+1, NSTATES)]
        elif NAC_PAIRS == "adjacent":
            pairs = [(i, i+1) for i in range(NSTATES-1)]
        else:
            raise ValueError('NAC_PAIRS must be "all" or "adjacent"')
        plot_nac_distribution(pairs, "nac_distribution")
    if "overlaps" in plots:
        same_frame = np.stack([read_step3_matrix(Path(DATA_DIRECTORY)/f"s_adi_{i}.txt")
                               for i in range(ISTEP, FSTEP)])
        if same_frame.shape != ham.shape:
            raise ValueError("Same-frame overlap and Hamiltonian shapes must match")
        plot_overlap_diagnostics()
    if "spectra" in plots:
        plot_gap_spectra()
    if "excitation" in plots:
        plot_excitation_analysis()
    if "uvvis" in plots:
        plot_energy_uvvis_dos()


if __name__ == "__main__":
    main()
