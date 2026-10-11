# Step 4: Nonadiabatic Molecular Dynamics within NBRA

This step simulates electronic population dynamics within the neglect of
back-reaction approximation (NBRA), using electronic energies and time overlaps
along a prescribed nuclear trajectory.

Use `NAMD_simulation.py` to check the adiabatic basis and run dynamics, then
`NAMD_plotting.py` to analyze the saved results. Edit the settings directly in
these files. The only command-line setting is `--prefix`, the common output
folder for a simulation run (default: `namd_run`).

Run from `4_step4` in a Python environment with Libra (`liblibra_core` and
`libra_py`), NumPy, SciPy, h5py, and Matplotlib available:

```bash
cd 4_step4
```

Input paths and relative prefixes are resolved from your current directory.
An absolute prefix also works. The datasets and their archives are beside the scripts; the initially selected
input is `data3/res`.

## 1. Check the basis

Section 2 of `NAMD_simulation.py` contains three dataset descriptions and
selection blocks. Comment out the current block's
code and uncomment the desired block, including its reader call. Keep exactly
one block active and leave the alternative numerical `basis_size` commented
for the first test.

| Dataset | Local input | Format | Frames (`fstep` excluded) | DT (fs) | Initial state |
| --- | --- | --- | --- | --- | --- |
| Example 1: step3 DFT tutorial, 11 states | `data1` | Sparse NPZ | 1200–1398 | 1.0 | 10 |
| Example 2: (TiO2)2, Miguel Recio | `data2` | Text | 0–398 | 1.0 | 2 |
| Example 3: Mg@C60, Kosar Yasin, 4 states | `data3/res` | Text | 0–48 | 0.5 | 1 |

Example 2 contains 11 states, all read with `basis_size = None`; setting
`basis_size = 10` selects its first 10 states. Commented `basis_size` choices
for Examples 1–3 are 11, 10, and 3, respectively. Preparation links and archive
extraction commands are included in the script comments.

For your own trajectory, edit `data_path`, `istep`, `fstep`, and `DT`. `fstep` is exclusive and `DT` is the input frame
spacing in femtoseconds. Leave `basis_size = None` to read all available states.
Leave this line active immediately after section 3:

```python
sys.exit(0)
```

Run:

```bash
python NAMD_simulation.py --prefix runs/example
```

The script prints determinant and singular-value diagnostics for every leading
active-space size, saves `analysis_data.npz` and `basis_diagnostics.txt` in the
prefix folder, and stops before decoherence calculations or dynamics.
`accepted=1` means the overlap block has full numerical rank at every interval;
it does not prove that the electronic populations are converged with basis size.

## 2. Choose the basis and run dynamics

Set `basis_size` to the chosen number of states and comment out the stop line:

```python
#sys.exit(0)
```

Sections 4–7 are already runnable. Edit `ISTATE`, the dynamics settings
(`nsteps`, `ntraj`, temperature), and `ICONDS`. States are zero-based.
Initial conditions are offsets into the loaded frames, so each must be between
`0` and `fstep - istep - 1`. Dynamics longer than the loaded trajectory repeats
the input frames; inspect the trajectory boundary before using this option.

Select exactly one recipe in section 6, for example:

```python
fssh_nbra.load(dyn_general); prf = "FSSH"
#ida_nbra.load(dyn_general); prf = "IDA"
```

Rerun with the same prefix:

```bash
python NAMD_simulation.py --prefix runs/example
```

The basis diagnostics still run first. With the stop line commented, the script
then computes decoherence statistics and runs initial conditions in a
multiprocessing pool. Edit `nthreads` in section 7 to set
the number of worker processes (default: 4); `ntraj` still controls the number
of hopping trajectories per initial condition. Match `nthreads` to your CPU
allocation. Each worker copies the parameter dictionaries, creates its own
random generator, and writes to its own initial-condition folder.

The pool explicitly uses Linux's `fork` start method so workers inherit the
already loaded Libra matrices without pickling them or rerunning the setup.
Run the script normally from the terminal; this implementation requires a
platform supporting `fork`. The pool is protected by an `if __name__ ==
"__main__"` guard and cleaned up automatically when the block finishes.

Each file in `recipes/` is a standalone, reusable recipe with active settings
and commented alternatives. To create a new recipe,
copy a file, rename it, edit its choices, import it in the simulation script,
and select its `load()` call. Keep one choice active per category. The supplied
recipes specify NAC evaluation, electronic propagation, thermal acceptance,
state tracking, phase correction, hopping, and decoherence settings. Read the
recipe comments when changing related settings together.

## 3. Plot the results

Edit `methods`, `ICONDS`, `ISTATE`, and `population_type` near the top of
`NAMD_plotting.py` to match the simulation. Set `results_path` to the saved run
folder; it is also the default output folder. Select plotting calls at the bottom.
Run `python NAMD_plotting.py`, or use `--prefix` to override the common run/output
folder. All other choices are edited in the script.
Basis plots can be made after the basis-only run. Electronic-structure plots
(energies, NACs, gap spectra, UV–vis, and excitation composition) are produced by
[Step 3 plotting.py](../3_step3/plotting.py); see the [Step 3 guide](../3_step3/README.md).
For a basis-only run, comment out the population, decoherence, and lifetime-fit
calls, plus TRPES and entropy. After dynamics, enable:

```python
plot_populations()
plot_decoherence_times()
plot_lifetime_fits()
plot_trpes()
plot_shannon_entropy()
```

Use the same prefix:

```bash
python NAMD_plotting.py --prefix runs/example
```

Figures and fit tables are saved **directly in `runs/example`**, alongside the
simulation results. Plotting never reruns dynamics. Uncomment `plt.show()` in
`save_plot()` if you also want figures displayed. On a machine without a display,
you can run `MPLBACKEND=Agg python NAMD_plotting.py --prefix runs/example`.

Fits use a stretched exponential for the selected state's decay. The commented
alternative in `exp_funct()` fits recovery. Choose a model consistent with the
initial and final populations; the default `R2_MIN = 0.1` is only an illustrative
teaching threshold.

## Additional dynamics analysis

The calls at the bottom of `NAMD_plotting.py` include TRPES and Shannon entropy.
Edit their settings near the top, then run with the same `--prefix` as the
simulation. They read saved populations and never launch TSH. For a basis-only
run, comment out all five post-dynamics calls. Electronic-structure analysis
requiring CP2K logs is in `3_step3/plotting.py`.

### Population-weighted time–energy map (TRPES)

After TSH finishes, use `plot_trpes()`. The
[TRPES tutorial](https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/6_dynamics/2_nbra_workflows/18_plotting_trpes/tutorial.ipynb)
describes the Gaussian population–energy convolution and alignment across
initial conditions. This computes a time-resolved
potential-energy-surface map in Libra's TRPES sense. It is **not a simulated
photoelectron spectrum**: no ionization energies, probe pulse, or photoionization
cross sections are supplied. It shows where the electronic population lies in
excitation energy as dynamics proceeds.

For each initial condition, the map sums normalized Gaussians centered at
`E_a(t)-E_0(t)`, weighted by the saved state populations:

$$
I(t,E)=\frac{1}{N_{\mathrm{IC}}}\sum_{r=1}^{N_{\mathrm{IC}}}
\sum_{a\in\mathcal A}p_{a,r}(t)
\frac{\exp[-(E-\epsilon_{a,r}(t))^2/(2\sigma^2)]}{\sqrt{2\pi}\sigma},
\qquad \epsilon_{a,r}(t)=E_{a,r}(t)-E_{0,r}(t).
$$

Here $r$ labels initial conditions, $\mathcal A$ is the selected state set, and
$p$ is the chosen SE or SH population. If $q_r$ is an initial-condition frame
offset and $N_f$ is the number of input frames, the energy alignment is

$$
f_r(t)=\left(q_r+\operatorname{round}\frac{t}{\Delta t}\right)\bmod N_f.
$$

The optional relative color scale and the overlaid mean energy are

$$
I_{\mathrm{rel}}(t,E)=\frac{I(t,E)}{\max_{t,E}I(t,E)},\qquad
\bar E_{\mathrm{win}}(t)=
\frac{\int_{E_{\min}}^{E_{\max}}E I(t,E)\,dE}
{\int_{E_{\min}}^{E_{\max}}I(t,E)\,dE}.
$$

The code uses trapezoidal integration on the energy grid. Zero signal remains
zero under peak normalization; a zero denominator gives an undefined (`NaN`)
mean energy. Over the full energy axis,
$\int I(t,E)dE=N_{\mathrm{IC}}^{-1}\sum_r\sum_{a\in\mathcal A}p_{a,r}(t)$.

| Setting | Meaning |
| --- | --- |
| `TRPES_POPULATION` | `"sh_pop_adi"` for active-surface populations or `"se_pop_adi"` for electronic populations. Independent of the lifetime-fit setting. |
| `TRPES_STATES` | `None` includes all states; a list such as `[1, 2, 3]` selects only excited states. State indices are zero-based. |
| `TRPES_EMIN`, `TRPES_EMAX` | Display/integration window in eV relative to S0. Include several Gaussian widths around the energies of interest. |
| `TRPES_DE` | Energy-grid spacing in eV; resolve the chosen Gaussian width. |
| `TRPES_SIGMA` | Gaussian standard deviation in eV, not FWHM; FWHM is approximately `2.355*sigma`. |
| `methods`, `ICONDS` | Completed methods and initial conditions to average with equal weight. |
| `TRPES_NORMALIZE_PEAK` | `True`: show intensity relative to the maximum, as in the tutorial; `False`: show population density in 1/eV. |
| `TRPES_CMAP` | Matplotlib color map, initially `"magma"`; alternatives include `"viridis"` and `"turbo"`. |
| `TRPES_MAX_POINTS` | `None`: all saved samples; a positive integer limits the plotted time range without changing input-energy periodicity. |

Energies come from `analysis_data.npz`; populations and actual saved times come
from each `METHOD_icond_N/mem_data.hdf`. Input frame selection is
`(icond + round(time_fs/DT)) % NSTEPS`, matching the simulation's initial offsets
and cyclic reuse. Population state columns must match the saved Hamiltonian
basis, and the initial-condition time grids must agree. Use the unmodified
archive from that simulation; a later basis-check run under the same prefix
could overwrite it. This alignment assumes the simulation's prescribed-frame
propagation and saved populations beginning at time zero.

Outputs are `METHOD_POPULATION_trpes.png` and a matching `.npz` containing
`time_fs`, `energy_ev`, and `density_per_ev`, plus the selected states, ICs, and
Gaussian width. The archive always preserves the unscaled density; its energy
integral is
the selected population when the energy window includes the Gaussian tails.
Selecting excited states therefore preserves population loss to S0. The
displayed map is peak-normalized by default. A white curve shows its
energy-weighted mean within the selected window; the archive also stores
`mean_energy_ev` and `display_intensity`. Window truncation biases this mean.
For fixed grid spacing and width, constant Gaussian prefactors cancel under
peak normalization. No temporal
broadening is applied. A finite energy window can omit signal.

### Shannon entropy of electronic populations

`plot_shannon_entropy()` calculates `S(t) = -sum_a p_a(t) ln p_a(t)`, with
`0 ln 0 = 0`. `ENTROPY_POPULATIONS` selects `sh_pop_adi`, `se_pop_adi`, or both.
The default `ENTROPY_NORMALIZE=True` divides by `ln(N)`, where `N` is the actual
number of saved state columns. A single occupied state has entropy 0; a uniform
population over all N states has normalized entropy 1. For a one-state basis,
the normalized entropy is defined as zero. With normalization disabled, entropy
is reported in nats.

Following the [Fullerene Shannon-entropy analysis](https://github.com/AkimovLab/Project_Fullerenes/blob/main/5_analysis/Shannon_entropy.ipynb),
the solid curve is the equally weighted mean of entropies computed separately
for each IC. The dashed curve is the entropy of the IC-averaged population.
These are different quantities: averaging populations first can add apparent
mixing due to differences among initial conditions. Both are saved, along with
individual IC entropies, for inspection. Populations must be finite,
nonnegative, and normalized; only roundoff-sized normalization drift is corrected.

Outputs are `POPULATION_shannon_entropy.png` comparing methods and
`METHOD_POPULATION_shannon_entropy.npz/.txt` for each method. This is a
basis-dependent population entropy, not von Neumann entropy, a direct measure
of quantum coherence, or thermodynamic entropy. Compare the same state basis,
population estimator, and initial conditions across methods.

For $N_s$ saved states, the entropy of initial condition $r$ is

$$
S_r(t)=-\sum_{a=0}^{N_s-1}p_{a,r}(t)\ln p_{a,r}(t),\qquad
\widetilde S_r(t)=\frac{S_r(t)}{\ln N_s}\quad(N_s>1).
$$

The normalized result is defined as zero for $N_s=1$. Using $0\ln0=0$,
the solid and dashed curves, respectively, are

$$
\langle S\rangle_{\mathrm{IC}}(t)=\frac{1}{N_{\mathrm{IC}}}\sum_r S_r(t),\qquad
S[\bar p](t)=-\sum_a\bar p_a(t)\ln\bar p_a(t),\qquad
\bar p_a(t)=\frac{1}{N_{\mathrm{IC}}}\sum_r p_{a,r}(t).
$$

Both curves are divided by $\ln N_s$ when normalization is enabled.
Concavity implies $S[\bar p]\geq\langle S\rangle_{\mathrm{IC}}$ for a common
state basis; their difference measures additional mixing introduced by pooling
initial conditions.

## Scientific background

The learning goals are to set up file-based NBRA dynamics with different
methods, assess the electronic basis, and interpret populations, decoherence
maps, and timescales. Libra's `generic_recipe` also supports model-Hamiltonian
calculations. Here, `compute_model` supplies stored electronic-structure data
through the same interface; the nuclear trajectory is prescribed rather than
updated in response to the evolving electronic populations.

### Input matrices and basis quality

Both readers return `NSTEPS`, `NSTATES`, and lists of Libra `CMATRIX` objects:
`St_adi` (successive-geometry overlaps), `Ham_adi` (adiabatic Hamiltonians),
`NAC_adi` (time-derivative couplings), `Hvib_adi` (vibronic Hamiltonians), and
`Basis_transform` (identity matrices in this interface). Data are read once and
reused by the model callback. NumPy copies support basis diagnostics and population–energy alignment.

The readers use HST couplings, `NAC = (St - St†)/(2*dt)`, and
`Hvib = Ham - i*NAC`, with `dt` in atomic units. The NPZ reader computes
HST couplings from the time-overlap matrices.

A nearly singular overlap can indicate projection outside the retained basis,
which is problematic for local diabatization (LD). Determinants generally
shrink with dimension but need not change monotonically: local maxima in the
minimum determinant suggest candidate sizes, not proof of completeness.
Compare the singular-value rank test and ultimately population convergence.
`nact_st` in the diagnostic plot means the same leading-state count as
`basis_size` in the simulation.

For a leading $m$-state block $T^{(m)}$ of the time-overlap matrix, the basis
plots and rank test use

$$
T^{(m)}=U\operatorname{diag}(\sigma_1,\ldots,\sigma_m)V^\dagger,\qquad
D_m(t)=|\det T^{(m)}(t)|=\prod_{a=1}^{m}\sigma_a(t),
$$

$$
s_{\min,m}(t)=\min_a\sigma_a(t),\qquad
c_m(t)=\max(\epsilon_{\mathrm{abs}},\epsilon_{\mathrm{rel}}\sigma_{\max}(t)),\qquad
r_m(t)=\sum_{a=1}^{m}\mathbf{1}[\sigma_a(t)>c_m(t)].
$$

The basis is accepted only if $r_m(t)=m$ at every interval. The plotting
function uses $\epsilon_{\mathrm{rel}}=10^{-4}$ and
$\epsilon_{\mathrm{abs}}=10^{-12}$. The optional individual-overlap plot shows
$|T_{ij}(t)|$ for the selected pair. The coupling construction is

$$
d(t)=\frac{T(t)-T^\dagger(t)}{2\Delta t},\qquad
H_{\mathrm{vib}}(t)=H_{\mathrm{adi}}(t)-i\hbar d(t),
$$

with $\hbar=1$ in the atomic-unit implementation.

### Gaps, decoherence, and dynamics controls

The average gap is `gaps[i,j] = mean(abs(Ei(t) - Ej(t)))`, rather than the
absolute difference of mean energies. `decoherence_times_ave` estimates
pairwise dephasing times and rates from gap fluctuations. The simulation saves
times in fs, rates in inverse atomic time, and gaps in Hartree. Recipes decide
whether to apply these rates; calculating them alone does not add decoherence.
`ave_gaps` is used by dephasing-informed corrections when enabled.

For the single input trajectory used here, define the absolute gap and its
population standard deviation over $N_t$ frames:

$$
g_{ab,n}=|E_a(t_n)-E_b(t_n)|,\qquad
\bar g_{ab}=\frac{1}{N_t}\sum_n g_{ab,n},\qquad
\sigma_{ab}=\sqrt{\frac{1}{N_t}\sum_n(g_{ab,n}-\bar g_{ab})^2}.
$$

The saved average-gap matrix is $\bar g_{ab}$. For $a\ne b$ and
$\sigma_{ab}>0$, Libra's estimator is

$$
\tau_{ab}=\sqrt{\frac{12}{5}}\frac{\hbar}{\sigma_{ab}},\qquad
\gamma_{ab}=\tau_{ab}^{-1}.
$$

The simulation evaluates these expressions in atomic units and converts
$\tau$ to fs for the map. The displayed diagonal is set to zero. For exactly
zero gap variance, the mathematical limit is $\tau\to\infty$ and
$\gamma\to0$; the inspected routine leaves the off-diagonal time and rate at
zero in this case. Such a zero time must not be interpreted as instantaneous
decoherence. This variance-based estimate is distinct from fitting the decay
of an electronic population.

| Setting | Meaning |
| --- | --- |
| `nfiles` / `nsteps` | Loaded frame count / propagation length; longer propagation cycles the stored frames. |
| `ntraj` | Surface-hopping trajectories per initial-condition ensemble. |
| `dt` | Propagation timestep in atomic units (`DT * units.fs2au`). |
| `progress_frequency` | Reporting fraction; 0.1 gives about ten progress reports. |
| `which_adi_states`, `which_dia_states` | State indices for supported state-resolved output, not active-basis selection. |
| `mem_output_level`, `properties_to_save` | Both must enable an observable; memory output is written to HDF5 at the end. |
| `prefix`, `prefix2` | Per-IC output folders under the command-line run prefix. |
| `istate`, `rep` in `elec_params` | Initially occupied state and initialization representation: 1 adiabatic, 0 diabatic. |

Keep `timestep`, `icond`, `model0`, and `nstates` in `model_params`: they describe
the current frame, starting frame, model identifier, and electronic dimension.
The nuclear initialization is a placeholder for this force-free NBRA example;
other electronic initialization parameters control amplitudes and density matrices.

Recipes control the hopping method, electronic propagator, and decoherence
approximation. LD handles basis rotations and phases internally. NAC/Hvib
propagation uses integrator 5; state tracking and phase correction are separate
choices, explored in the exercise below. Change related settings together when
customizing a recipe. The supplied FSSH2 recipe has no extra decoherence;
combining FSSH2 with MSDM and precomputed rates requires a customized recipe.

### Interpreting the plots and fitted times

`NAMD_plotting.py` uses Matplotlib to plot the saved SE/SH populations. `ICONDS = [0]`
selects a single IC; multiple entries give an equally weighted average. SE
populations describe electronic amplitudes, while SH populations count active
surfaces. They are distinct estimators and need not coincide.

For $N_r$ hopping trajectories within initial condition $r$, the two
population estimators are

$$
p^{\mathrm{SE}}_{a,r}(t)=\frac{1}{N_r}\sum_{\ell=1}^{N_r}|c_{a,r\ell}(t)|^2,
\qquad
p^{\mathrm{SH}}_{a,r}(t)=\frac{1}{N_r}\sum_{\ell=1}^{N_r}
\mathbf{1}[A_{r\ell}(t)=a].
$$

$c$ is the adiabatic electronic amplitude and $A$ the active surface. The HDF5
files already contain this within-IC average. The plotted curves then use

$$
\bar p_a^{X}(t)=\frac{1}{N_{\mathrm{IC}}}\sum_r p_{a,r}^{X}(t),
\qquad X\in\{\mathrm{SE},\mathrm{SH}\}.
$$

Every selected IC has equal weight, even if its number of hopping trajectories
differs. This is not a trajectory-count-weighted pooling of all runs.

Excitation energies follow the trajectory relative to S0. Adjacent-gap
histograms describe energy fluctuations. NAC maps average absolute couplings
to avoid sign cancellation; distributions compare all unique pairs with only
adjacent pairs. Log plots omit zero couplings and report density per meV.
The decoherence map shows pairwise times; its zero diagonal is a display
convention, not a zero self-decoherence lifetime.

The fitted decay model is `P(t) = exp[-(t/tau)^beta]`; recovery uses `1-P(t)`.
Both `tau` and `beta` are fitted; `beta=2` gives Gaussian decay. Fit quality is
`R² = 1 - sum((P-Pfit)^2)/sum((P-mean(P))^2)`. Constant traces have no meaningful
R² and are skipped. The mean includes only fits passing `R2_MIN`; `tau` is a
1/e decay time, not generally the integrated lifetime.

The approximate 95% half-width of the mean is `1.96*s/sqrt(N)`, where `N` is
the number of accepted IC fits and `s` is their sample standard deviation.
It is not a prediction band, and a single fit cannot estimate it. Nearby
starting frames may be correlated. Longer dynamics and more hopping trajectories
are needed to judge convergence beyond this teaching example's low R² cutoff.

For each IC, the least-squares lifetime model and its recovery alternative are

$$
P_{\mathrm{fit}}(t;\tau,\beta)=e^{-(t/\tau)^\beta},\qquad
P_{\mathrm{recovery}}(t)=1-e^{-(t/\tau)^\beta},\qquad
(\widehat\tau,\widehat\beta)=\operatorname*{arg\,min}_{\tau,\beta}
\sum_n[P(t_n)-P_{\mathrm{fit}}(t_n;\tau,\beta)]^2.
$$

The implemented bounds are $\tau\geq10^{-12}$ fs and
$0.05\leq\beta\leq5$. Fit quality is

$$
R^2=1-\frac{\sum_n[P(t_n)-P_{\mathrm{fit}}(t_n)]^2}
{\sum_n[P(t_n)-\bar P]^2},\qquad
\bar P=\frac{1}{N_t}\sum_nP(t_n).
$$

Only fits with $R^2>\texttt{R2\_MIN}$ enter the reported mean. For either
parameter $x\in\{\tau,\beta\}$ across $M$ accepted IC fits,

$$
\bar x=\frac{1}{M}\sum_{r=1}^{M}x_r,\qquad
s_x=\sqrt{\frac{1}{M-1}\sum_r(x_r-\bar x)^2},\qquad
h_x=1.96\frac{s_x}{\sqrt M}\quad(M>1).
$$

$h_x$ is the approximate 95% half-width reported for the mean, not an error
bar on an individual fitted curve. The integrated lifetime of this ideal decay
model is

$$
\int_0^\infty P_{\mathrm{fit}}(t)\,dt=\tau\,\Gamma(1+1/\beta).
$$

The script reports $\tau$, not this integral.

Each initial condition runs with separate worker parameters and an output
folder. Compare populations and fitted times across repeated runs within
sampling uncertainty; stochastic trajectories need not be identical.

## Exercise: compare methods for the four-state Mg@C60 model

**Goal:** use the calculations in the `MgC60_*` folders to investigate how
propagation, state tracking, phase correction, and hopping proposals affect
populations. Develop a physically motivated interpretation, then identify the
checks needed to distinguish numerical artifacts from genuine dynamics.

### A. Inspect the model before running dynamics

Select Example 3 (`data3/res`), initially set `basis_size = None`, and enable
`sys.exit(0)` after the basis test. Save this inspection under a separate prefix.
Examine the diagnostics for sizes 1–4, then use `basis_size = 4` for the method
comparisons. Comment out the stop line to proceed into dynamics.

Inspect excitation energies, adjacent-gap distributions, and NAC statistics
with the Step 3 analysis for the matching electronic-structure dataset. Use
the Step 4 basis plots for overlap determinants and minimum singular values. After section 4 has run, also inspect the decoherence
time map. Address these questions:

- Which states are close in energy, and which pairs have the strongest couplings?
- Which transitions would you expect to compete with population remaining on S1?
- Are the estimated decoherence times shorter or longer than the plotted dynamics?
- Does the three-state truncation have better or worse overlap support than the
  four-state space? Why is passing a rank threshold insufficient to establish
  population convergence?

Use the time traces as well as averages: a mean coupling does not specify its
sign changes, correlations, or the complete population-transfer dynamics.

### B. Reproduce the method comparisons

Use the same input frames (`istep=0`, `fstep=49`), `DT=0.5 fs`, initial state
`ISTATE=1`, four-state basis, temperature, `ICONDS=[0,5,10]`, and propagation
length. Select one recipe at a time and give every calculation a distinct
`--prefix`. Keep the number of trajectories equal for a controlled comparison.

The existing calculations provide the following comparison set. The table
reflects their saved `_dyn_params.txt` files; inspect those files when identifying
a run, since an editable recipe may have changed after the calculation.

| Existing folder | Recipe file in `recipes/` | Electronic integrator | State tracking | MK scaling | Phase correction |
| --- | --- | --- | --- | --- | --- |
| `MgC60_fssh` | `fssh_nbra.py` | LD, 2 | LD, -1 | 1 (inactive for LD tracking) | Off |
| `MgC60_fssh_v2` | `fssh_nbra_v2.py` | LD, 2 | SVD-based LD, 5 | 3 (inactive for LD tracking) | Off |
| `MgC60_fssh2` | `fssh2_nbra.py` | LD, 2 | LD, -1 | 1 (inactive for LD tracking) | Off |
| `MgC60_fssh_nac_v0` | `fssh_nbra_nac_v0.py` | NAC/Hvib, 5 | None, 0 | 3 (inactive without tracking) | Off |
| `MgC60_fssh_nac_v1` | `fssh_nbra_nac_v1.py` | NAC/Hvib, 5 | Hungarian, 21 | 3: no scaling | Off |
| `MgC60_fssh_nac_v2` | `fssh_nbra_nac_v2.py` | NAC/Hvib, 5 | Hungarian, 21 | 3: no scaling | On |
| `MgC60_fssh_nac_v3` | `fssh_nbra_nac_v3.py` | NAC/Hvib, 5 | Hungarian, 21 | 1: energy weighting | On |
| `MgC60_fssh_nac_v4` | `fssh_nbra_nac_v4.py` | NAC/Hvib, 5 | LD, -1 | 3 (inactive for LD tracking) | On |

“Phase correction” means `do_phase_correction`; `do_nac_phase_correction=0` in
all eight saved runs. All use FSSH (`tsh_method=0`) except `MgC60_fssh2`, which
uses FSSH2 (`tsh_method=7`). All disable additional decoherence
(`decoherence_algo=-1`, `decoherence_times_type=-1`): saving estimated times
does **not** mean those times were applied to propagation.

The saved `MgC60_fssh` calculation currently uses **2,500 trajectories per IC**;
the other seven use **250**. Treat this as a sampling difference, or rerun a
matched pair under new prefixes before attributing a difference solely to the
method. Each saved run has 49 steps; its plotted time points extend to 24 fs.
The nonzero IC offsets wrap through the finite input sequence during these runs.

For example, after selecting `fssh_nbra_v2.load(dyn_general); prf = "FSSH"`:

```bash
python NAMD_simulation.py --prefix runs/exercise_fssh_svd
python NAMD_plotting.py --prefix runs/exercise_fssh_svd
```

In the plotting file, match `methods`, `ICONDS`, and `ISTATE` to the simulation
and uncomment the population, decoherence-map, and lifetime-fit calls. For the
FSSH2 run, use `prf = "FSSH2"` and `methods = ["FSSH2"]`. Pool size controls
execution resources; it is not a replacement for increasing `ntraj`.

### C. Compare populations and explain the differences

Make the following pairwise comparisons, keeping the other settings fixed:

1. **LD formulation:** `fssh` versus `fssh_v2`. Compare the two treatments of
   basis transport, especially where the overlap deviates from a unitary matrix.
2. **Hop proposal:** `fssh` versus `fssh2`. Check SE and SH populations separately.
   With the same electronic propagation and no decoherence, do changes in the
   hopping proposal affect the two estimators in the same way?
3. **State tracking:** `nac_v0` versus `nac_v1`. What changes when Hungarian
   assignment is introduced without phase correction?
4. **Phase correction:** `nac_v1` versus `nac_v2`. Does correcting the phases
   change transfer, oscillations, or apparent relaxation times?
5. **Energy weighting:** `nac_v2` versus `nac_v3`. How does the assignment's
   energy weighting change dynamics within nearby versus well-separated states?
6. **Mixed propagation choices:** examine `nac_v4` against the LD and NAC runs.
   It combines integrator 5 with LD tracking; treat this as a consistency test,
   not automatically as an improved method. More than one setting differs from
   some references, so identify those differences before drawing conclusions.

For every comparison, inspect both `*_SE_populations.png` and
`*_SH_populations.png`. Record final populations, maximum transient S2/S3
populations, ground-state accumulation, and population returns or oscillations.
Check normalization and compare individual ICs (`ICONDS=[0]`, then `[5]`, then
`[10]`) with their average. Are large features common to all ICs or driven by one?

Explain the results using energies, couplings, and decoherence estimates, but
also account for the common thermal hop-acceptance rule. In particular, why
might wavefunction population on a high-energy state coexist with few accepted
hops to that state? Agreement between SE and SH alone is not a proof of accuracy,
and a visually smooth or rapidly relaxing curve need not be the best result.

### D. Assess fits and convergence

Fit the initially occupied state's population and inspect the fitted curves,
`tau`, `beta`, and R² in `*_fits.txt`. Repeat with the other population estimator.
If studying ground-state recovery, choose `ISTATE=0` in the plotting file and
the recovery expression in `exp_funct()`; do not change the simulation's initial
state for this analysis. Decide whether a monotonic decay/recovery model is
appropriate for traces with substantial population returns or incomplete decay.
Report how many fits pass the cutoff rather than quoting a lifetime alone.

Choose one disputed comparison for a limited convergence study: increase
`ntraj` at fixed ICs, then examine sensitivity to the selected starting frames.
More trajectories reduce stochastic sampling noise but cannot repair an
inconsistent propagator or basis. Longer propagation through repeated input is
not a longer independent nuclear trajectory. Changing `DT` alone while using
the same overlaps also changes the model's time scale; a timestep study needs
consistently sampled inputs or an appropriate electronic-substep study.

As an extension, compare FSSH with IDA, MSDM, or DISH using the provided recipes
and new prefixes. Identify which choices actually consume the estimated
rates, then examine whether decoherence changes oscillations, populations, and
fitted times. Keep the distinction between changing a physical approximation
and testing numerical convergence explicit.

**Deliverable:** a small comparison table of settings and population metrics,
selected overlaid or side-by-side plots, and a short argument for the most
credible dynamics. State the evidence, unresolved discrepancies, and the next
validation check; do not declare a method accurate from visual agreement alone.

## Output layout

```text
runs/example/
    analysis_data.npz
    basis_diagnostics.txt
    decoherence_times.txt
    decoherence_rates.txt
    average_gaps.txt
    FSSH_icond_0/mem_data.hdf
    FSSH_icond_5/mem_data.hdf
    FSSH_icond_10/mem_data.hdf
    overlap_determinants.png
    FSSH_sh_pop_adi_trpes.png
    sh_pop_adi_shannon_entropy.png
    ... other selected figures ...
    FSSH_fits.txt
```

Only the first two files are written before the stop line. Later files depend
on which simulation and plotting sections you run. Reusing a prefix overwrites
matching output files; use a new prefix to preserve a different calculation.
For example, keep a full-basis inspection in `runs/basis_check` and a production
run in `runs/fssh_basis3`. Both scripts also work without arguments, using
`namd_run`.

