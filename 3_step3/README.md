# Step 3: TD-DFT state overlaps and vibronic Hamiltonians

This step converts the Kohn–Sham (KS) orbital overlaps from Step 2 and the
TD-DFPT excitation amplitudes printed by CP2K into overlaps between many-electron
TD-DFT pseudowavefunctions. It then estimates time-derivative couplings and builds
vibronic Hamiltonians for subsequent nonadiabatic dynamics.

`run.py` calls `libra_py.packages.cp2k.methods.cp2k_compute_adi` for a single
trajectory. This is postprocessing: this function reads existing CP2K results;
it does not launch CP2K or propagate nuclei or electronic populations.

The explanations below follow Libra's CP2K interface implementation, including
its `citools` overlap machinery.
The exact parameter names are **`nelec_act_space`**, **`lowest_orbital`**, and
**`ci_threshold`**. Python keys are case-sensitive: `nelec_act_Space` and
`lowst_orbital` are not aliases. A misspelled active-space key is ignored by the
lookup; a missing `lowest_orbital` raises a `KeyError`.

The companion [NBRA workflow tutorial](https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/11_program_specific_methods/3_cp2k_methods/5_namd_workflow/tutorial.ipynb)
demonstrates parsing TD-DFT data, constructing CI overlaps explicitly, calling
the Libra implementation, and comparing full and truncated determinants.
Its use cases are many-body TD-DFT NACs, CI wavefunction time overlaps, and an
interface Hamiltonian for Libra dynamics. Here, “NAC” means the time-derivative
coupling along the nuclear trajectory, not a Cartesian derivative-coupling
vector. The tutorial motivates Löwdin orthogonalization by both truncation and
the finite precision of CP2K's printed coefficients.

Use the tutorial's Sections 3–5 for the workflow and convergence ideas, but
follow the current source for API details: its early parser example uses
`tolerance`, whereas the current parser reads `ci_threshold`; its prose labels
parsed energies as eV, whereas the current parser returns hartree. Its C20
orbital indices and active-electron counts must be adapted to this TiO2 dataset.

## Inputs and running

Run from `3_step3` with a Python environment containing Libra (including
`liblibra_core`), NumPy, and SciPy:

```bash
cd 3_step3
python run.py
```

Prepare these inputs from the same trajectory and Step 2 calculation:

| Input | Purpose |
| --- | --- |
| `Rutile_TiO2_MD-pos-1.xyz` | Step 1 trajectory; frames are numbered from zero. Coordinates are read in angstrom and converted to bohr. |
| `all_logfiles/step_i.log` | CP2K electron count, excitation energies, and excitation-analysis amplitudes for frame `i`. |
| `res/S_ks_i.npz` | Same-frame KS overlap matrix at frame `i`. |
| `res/St_ks_i.npz` | KS time-overlap matrix; the Step 2 writer uses `i` for the interval `i → i+1`. |

The KS matrices contain alpha and beta spin blocks and have shape
`(2*Norb, 2*Norb)`, where `Norb` is the number of stored spatial orbitals.
`E_ks_i.npz` is not read by this workflow: state energies come from the TD-DFPT
logs. Do not mix overlap files from different orbital windows.

The supplied loop processes frames `0` through `19` (`range(20)`). It creates
`ci_tio2` if needed and overwrites matching output files on reruns. Adjust the
loop and paths for your dataset. See the time-indexing issue below before using
the resulting couplings for dynamics.

## Workflow inside `cp2k_compute_adi`

1. Read the current CP2K log, convert excitation energies from eV to hartree,
   and retain excitation amplitudes satisfying `coefficient**2 > ci_threshold`.
2. Construct the reference determinant and select the active spatial orbitals.
   For this example, the log reports 48 electrons, 24 occupied orbitals, and
   10 excited states; the default reference is a closed-shell singlet.
3. Combine the retained single excitations from the previous and current calls
   into a common basis. Build spin-adapted configuration-state functions (CSFs),
   evaluate their overlaps from the KS matrices, and contract with the CI
   coefficient matrices: `T_CI = C_prev.T @ T_CSF @ C_curr`.
4. Compute the current same-frame CI overlap `S_CI`. Apply symmetric Löwdin
   orthogonalization, using `X = S_CI**(-1/2)`: the saved same-frame overlap is
   `X_curr @ S_CI @ X_curr`, and the time overlap is
   `X_prev @ T_CI @ X_curr`. The previous inverse square root is cached.
5. Set the ground-state energy to zero. For excited state `a`, use the average
   excitation energy from the previous and current calls:
   `E_a = (E_a_prev + E_a_curr)/2`. These are excitation energies relative to
   the ground state, not absolute CP2K total energies.
6. For real time overlaps `T`, form `d_ab = (T_ab - T_ba)/(2*dt)` and
   `H_vib,ab = -1j*d_ab` off the diagonal; the diagonal is `E_a`.
7. Cache the current electronic data, coordinates, and inverse overlap square
   root in `params`, then set `is_first_time[itraj] = False`.

Keep the same `params` dictionary across successive calls and process frames in
chronological order. This function does not apply CI phase correction or state
reordering. Forces and coordinate-resolved derivative couplings are placeholders
in this implementation; the time-overlap-based couplings enter `hvib_adi`.

## Parameters used by `run.py`

| Key | Supplied value | Meaning |
| --- | --- | --- |
| `atom_labels` | Labels read from the XYZ file | Required by the wrapper; currently not used in its electronic overlap calculation. |
| `timestep` | `i` | Bookkeeping in this script; the function reads the explicitly supplied filenames, not this value, to select data. |
| `dt` | `1.0*units.fs2au` | Time between adjacent input frames, in atomic units (1 fs is about 41.34 a.u.). Match the actual saved-frame spacing. Couplings scale as `1/dt`. |
| `nelec_act_space` | `None` | Use the entire stored orbital window, or choose an occupied-electron cutoff as described below. |
| `lowest_orbital` | `1` | Global, one-based CP2K index of the first spatial orbital in the stored KS matrices. |
| `ci_threshold` | `0.001` | Strict cutoff on squared excitation coefficients. |
| `nstates` | `11` | Total number of states, including ground state `0`; here, ground state plus 10 excited states. Each log must supply at least `nstates-1` excited states. |
| `is_first_time` | `{0: True}` | Initialize the cache for trajectory `0`. On this call, current CI data are also used as previous data. The function changes the flag to `False`. |
| `act_state` | `{0: 1}` | Active-state index for trajectory `0`; currently unused in the enabled calculation. It does not select which overlaps to compute. |
| `logfile_name` | `all_logfiles/step_i.log` | Current-frame TD-DFPT log. |
| `overlap_filename` | `res/S_ks_i.npz` | Current-frame KS overlap. |
| `time_overlap_filename` | `res/St_ks_i.npz` | KS overlap between two frames; the supplied script's indexing needs the adjustment discussed below. |

`full_id = Py2Cpp_int([0, 0])` selects trajectory `0` through its last element;
it is a function argument, not a parameter-dictionary entry. `q` contains its
coordinates. Optional source-level settings include `isUKS=False`,
`multiplicity=1`, `spin_projection` (defaults to `S`), and `reference_det`
(a determinant specified with positive alpha and negative beta orbital indices).
The closed-shell examples below assume the supplied restricted singlet setup.

## `lowest_orbital`: identify the stored orbital window

This key labels the first orbital already present in the KS matrices. It does
not crop those matrices. The wrapper infers the upper endpoint from their size:

```python
Norb = st_mo.shape[0] // 2
orbital_space = list(range(lowest_orbital, lowest_orbital + Norb))
highest_orbital = lowest_orbital + Norb - 1
```

For this repository:

| Step 2 dataset | Stored spatial orbitals (inclusive) | KS matrix shape | Step 3 `lowest_orbital` |
| --- | --- | --- | --- |
| Main/reference setup | `1..40` | `80 × 80` | `1` |
| `smaller_orbital_space` setup | `14..35` | `44 × 44` | `14` |

The comment in `3_step3/run.py` mentions `highest_orbital=30`, but the current
`2_step2/run_template.py` specifies `40`, and the supplied
`3_step3/res/St_ks_0.npz` is `80 × 80`. Use the actual Step 2 settings and files.

For the smaller window, use `lowest_orbital=14`, even though its first array row
has Python index zero. Setting it to `1` would mislabel the stored orbitals.
All active orbitals and retained transitions must fit in the stored window;
changing this parameter cannot recover orbitals omitted during Step 2.

## `nelec_act_space`: occupied electrons retained in the CI overlap

`None` means **all spatial orbitals in the stored overlap window**, not all
orbitals calculated by CP2K. The full log here reports 74 molecular orbitals,
whereas the supplied overlaps contain only 40 spatial orbitals. In the singlet
implementation, the effective electron count is twice the number of active
orbitals at or below the HOMO. Thus `None` with `1..40` keeps 48 electrons;
with `14..35`, it keeps 22 electrons in occupied orbitals `14..24`.

For an integer value, the wrapper sorts the absolute orbital indices in the
reference determinant (including both spins), selects
`min_indx = occupied[-nelec_act_space]`, and builds the inclusive spatial range
`min_indx..max_vir`. Here `max_vir` is the largest virtual index in the
**current log's threshold-retained excitations**. For a closed shell with HOMO
`H` and an even requested electron count `N`, this gives:

```text
first active occupied orbital = H - N/2 + 1
active spatial orbitals       = first active occupied orbital .. max_vir
```

Examples for `H=24`:

| Setting | Active occupied orbitals | Effective active electrons |
| --- | --- | --- |
| `nelec_act_space=2` | `24` | `2` |
| `nelec_act_space=4` | `23..24` | `4` |
| `nelec_act_space=22` | `14..24` | `22` |
| `nelec_act_space=48` | `1..24` | `48` |

The virtual range extends through `max_vir`; this parameter is not the number
of virtual orbitals, excited states, or excitation configurations. Use even
values for this closed-shell example. Odd values select a spatial cutoff that
includes complete occupied pairs, so the effective count can exceed the request.
The wrapper validates integer counts against `1..len(reference_det)`.

For example, with overlaps covering `14..35`:

```python
params["lowest_orbital"] = 14
params["nelec_act_space"] = 22  # occupied orbitals 14..24
```

This is valid only if all retained excitations at both endpoints fit the active
space, including the virtual bound. The active-space setting does **not**
automatically filter out incompatible excitation amplitudes: excluding an
occupied orbital still used in a retained excitation can cause an error such as
`Occupied orbital ... not present in reference determinant`.

The tutorial gives a useful lower bound. If `min_occ` is the lowest occupied
orbital involved in any retained excitation and the HOMO is `H`, then a
closed-shell active space must satisfy:

```text
nelec_act_space >= 2 * (H - min_occ + 1)
```

Equality is sufficient for this occupied-orbital coverage condition (the
tutorial's inline implementation also accepts equality). For example, a
retained excitation from orbital 18 with `H=24` requires at least 14 active
electrons, covering occupied orbitals `18..24`. Take the minimum `min_occ`
over all frames and states used in the comparison, rather than only the first
frame. Virtual-orbital coverage must also be checked. This bound establishes
whether the expansion can be represented; it does not establish NAC convergence.

Because `max_vir` is recomputed from the current log, an integer setting can
also change the virtual extent between calls. Check the union of previous and
current retained configurations. `None` provides a fixed active window when
the input matrix dimensions and labels stay fixed. Reducing the occupied space
reduces determinant sizes and changes the overlap approximation; assess
convergence of the resulting overlaps and couplings before relying on it.

## `ci_threshold`: truncate excitation weights

The parser keeps an excitation only when:

```python
coefficient**2 > ci_threshold
```

The comparison is strict and sign-independent. With `ci_threshold=0.001`,
`abs(coefficient)` must exceed approximately `0.03162`. A coefficient of `0.02`
has weight `0.0004` and is discarded; `0.04` has weight `0.0016` and is retained.
The cutoff is dimensionless and applies to each printed excitation individually,
not to a cumulative retained weight. It does not guarantee 99.9% of a state's
norm is retained.

A smaller threshold retains more configurations and usually increases the cost
of building determinant/CSF overlaps. A larger threshold can discard important
contributions or leave a state empty. The parser preserves retained coefficient
values; the wrapper subsequently orthogonalizes the CI states with Löwdin's
method. This cannot recover discarded physical contributions. Eigenvalues at
or below `1e-10` are suppressed by the inverse-square-root routine, so an empty
or linearly dependent state can produce a deficient overlap instead of an
identity matrix.

Compare, for example, `1e-3`, `1e-4`, and `1e-5` while checking retained weights,
overlap matrices, and couplings. Lowering the cutoff can introduce transitions
outside your orbital window; revisit `lowest_orbital`/Step 2 bounds and
`nelec_act_space` together. Libra can only retain amplitudes actually printed
in the CP2K logs; lowering its cutoff cannot restore unprinted contributions.

**Smaller `ci_threshold` enables more singly excited Slater-determinant (SD)
configurations to be included. These can involve deeper occupied and higher
virtual orbitals, so use larger orbital active spaces as needed.**
`nelec_act_space=None` is the ultimate untruncated choice within the available
KS overlap window: it uses all stored orbitals. If newly retained excitations
reach outside that window, even `None` is insufficient; regenerate Step 2
overlaps with a wider range and set `lowest_orbital` to its actual lower bound.

## Time indexing: align both endpoints

In the inspected Step 2 source, `St_ks_i` represents `<MO(i)|MO(i+1)>`.
However, after initialization, `cp2k_compute_adi` uses the cached CI data from
its previous call and the current log. The supplied loop loads logs `i-1` and
`i` but passes `St_ks_i`, which belongs to `i` and `i+1`. The first call also
uses frame-0 CI coefficients on both sides of a frame-0-to-1 KS overlap.
These endpoints are inconsistent.

For unrenamed files produced by this Step 2 implementation, replace the
`time_overlap_filename` assignment in the loop with:

```python
# Frame 0 initializes the cache using a same-frame overlap.
# Later calls compute the physical interval (i-1) -> i.
params["time_overlap_filename"] = (
    f"res/S_ks_{i}.npz" if i == 0 else f"res/St_ks_{i-1}.npz"
)
```

Keep `logfile_name=step_i.log` and `overlap_filename=S_ks_i.npz`. Treat output
`0` as initialization, and outputs `1..19` as the 19 physical intervals
`0→1` through `18→19`. To obtain 20 physical intervals, process frames `0..20`
and provide their logs and same-frame matrices plus `St_ks_0..19`.
This is a documented adjustment; the supplied `run.py` has not been changed.

## Outputs and checks

| Output | Contents |
| --- | --- |
| `ci_tio2/ham_adi_i.txt` | Diagonal adiabatic Hamiltonian in hartree; ground energy zero, excited energies averaged across the two calls. |
| `ci_tio2/hvib_adi_i.txt` | Complex vibronic Hamiltonian in atomic units, including imaginary off-diagonal time-derivative coupling terms. |
| `ci_tio2/s_adi_i.txt` | Real part of the Löwdin-orthogonalized same-frame CI overlap; dimensionless. |
| `ci_tio2/st_adi_i.txt` | Real part of the orthogonalized CI time overlap; dimensionless. |
| `coord-i.xyz` | Individual frames written by the XYZ reader; these files can be deleted after processing. |

With `nstates=11`, each electronic matrix is `11 × 11`, with the ground state
first. For a full-rank retained basis, `s_adi_i` should be close to identity.
The time overlap need not be identity. Check for finite entries, a Hermitian
`hvib_adi`, and stable results when the active space and CI cutoff are enlarged
and reduced, respectively. Abrupt changes can also reflect uncorrected phases
or state ordering; this routine does not resolve those issues automatically.

## Plot electronic structure after Step 3

Run `plotting.py` after `run.py` to analyze the saved electronic structure before
NA-MD. The script does not run CP2K, build new overlaps, or propagate dynamics.
It provides energy traces, gap distributions, NAC maps/distributions, overlap
diagnostics, gap spectra, excitation composition, and a combined energy/DOS/
UV–vis figure. NumPy, Matplotlib, and Libra are required; excitation parsing uses
Libra's CP2K interface. All paths are relative to your working directory.

```bash
cd 3_step3
python run.py
python plotting.py --prefix plots
```

Use the time-indexing correction described above when preparing NACs;
`plotting.py` cannot fix incorrectly paired CI and KS overlaps. By default,
initialization output 0 is omitted from NAC statistics (`SKIP_INITIAL_NAC=True`).
Set this flag to `False` only if output 0 is a genuine physical interval.
The Hamiltonian energy traces retain every selected output, including 0.

Before running, edit the settings near the top of `plotting.py`:

- `DATA_DIRECTORY`: Step 3 matrix directory (default `"ci_tio2"`).
- `EXC_LOG_DIRECTORY`: CP2K log directory (default `"all_logfiles"`).
- `ISTEP`, `FSTEP`: matrix and log file indices, end excluded (default `0, 20`).
- `DT`: frame spacing in fs (default 1), shared by matrix and log analyses.
- `OUTPUT_DIRECTORY`: default plot/data destination (default `"plots"`).
- `PLOTS`: analyses to run; comment out entries to disable them.

For example, edit `PLOTS` to one of these choices:

```python
PLOTS = ["energies", "gaps", "nacs", "overlaps"]
#PLOTS = ["spectra"]
#PLOTS = ["excitation", "uvvis"]
```

Keep one assignment active. Then run `python plotting.py`, or use
`python plotting.py --prefix plots/example` to override only the output folder.
All input and analysis choices are edited in the script; `--prefix` is the only
command-line setting.

Matrix plots read `ham_adi_i.txt`; NAC plots additionally read `st_adi_i.txt`,
and orthonormality diagnostics read `s_adi_i.txt`. Both real and Libra
`(real,imag)` text entries are supported. The NAC estimate is
`(St-St†)/(2*dt)` and is displayed as coupling energy in meV. Matrix dimensions
must agree across frames. The `overlaps` selection reads only the same-frame
overlaps and plots their deviation from identity, saving
`overlap_orthonormality.png/.npz`. Time-overlap basis diagnostics are performed
in Step 4 using SVD singular values.

### Equations for energies, gaps, NACs, and overlaps

Let $n=0,\ldots,N_t-1$ label the selected frames, $a,b$ label electronic
states, and $\Delta t$ be the saved frame spacing. The energy traces use

$$
t_n=n\Delta t,\qquad
\epsilon_a(t_n)=E_a(t_n)-E_0(t_n),\qquad
g_a(t_n)=E_{a+1}(t_n)-E_a(t_n).
$$

$\epsilon_a$ is plotted relative to S0; $g_a$ is the adjacent-state gap.
Energies are converted from hartree to eV. Each gap histogram is a probability
density: for bin $B_b=[x_b,x_{b+1})$ of width $\Delta x_b$,

$$
\rho_b=\frac{\sum_{n=0}^{N_t-1}\mathbf{1}[g_a(t_n)\in B_b]}
{N_t\,\Delta x_b},\qquad \sum_b\rho_b\Delta x_b=1.
$$

Here $\mathbf{1}$ is an indicator function; the final histogram bin includes
its right edge. Signed gaps are retained, so crossings can produce negative
values when state labels are not energy-ordered.

For time-overlap matrix $T$ and time step in atomic units, the estimated
coupling and vibronic Hamiltonian are

$$
d_{ab}(t_n)=\frac{T_{ab}(t_n)-T_{ba}^{*}(t_n)}{2\Delta t},\qquad
H^{\mathrm{vib}}_{ab}=E_a\delta_{ab}-i\hbar d_{ab}.
$$

The plotting script uses $\hbar=1$ in atomic units. If $K$ physical intervals
remain after excluding initialization, the NAC heatmap shows

$$
A_{ab}=\frac{1}{K}\sum_{n\in\mathcal I}|\hbar d_{ab}(t_n)|,
\qquad A_{aa}=0,
$$

converted to meV. The NAC distribution histograms the samples
$x_{n,ab}=|\hbar d_{ab}(t_n)|$ for either $\mathcal P=\{(a,b):a<b\}$ or
$\mathcal P=\{(a,a+1)\}$. Only strictly positive samples are retained:


$$
\rho_j^{\mathrm{NAC}} =
\frac{
\left|\left\lbrace (n,a,b) : (a,b)\in\mathcal{P},
x_{n,ab}>0,
x_{n,ab}\in B_j \right\rbrace\right|
}{
N_+ \Delta x_j
}.
$$


$N_+$ is the number of positive samples, so this is a distribution conditional
on nonzero coupling. Logarithmically spaced bins still give density per meV,
not per logarithmic unit; their plotted centers are $\sqrt{x_bx_{b+1}}$.
Set `NAC_PAIRS="all"` (default) or `NAC_PAIRS="adjacent"` in the script to
choose one pair set. One `nac_distribution.png` is written per run. Gap
distributions, the NAC map, and the NAC distribution are generated only here;
Step 4 consumes couplings for propagation but does not repeat these plots.

The same-frame overlap check plots

$$
\eta(t_n)=\max_{a,b}|S_{ab}(t_n)-\delta_{ab}|.
$$

$S$ is the saved same-frame CI overlap after Löwdin orthogonalization.
This checks orthonormality, not convergence of the NACs. The SVD time-overlap
rank test is described in the [Step 4 guide](../4_step4/README.md).

### Energy levels with sideways UV–vis and density of states

Including `"uvvis"` in `PLOTS` produces `energy_uvvis_dos.png`: energy levels versus time on
the left, with horizontally oriented excited-state DOS and UV–vis panels
sharing the vertical energy axis. It also produces `uvvis_spectrum.png`,
`uvvis_spectrum.txt`, `excitation_dos.txt`, and `energy_uvvis_dos.npz` with all
numerical arrays. This selection reads CP2K logs only.

The energy traces and optical spectrum use the **same CP2K excitation energies
and frames**, rather than combining single-frame oscillator strengths with the
midpoint energies in Step 3 Hamiltonians. The parser reads the `TDDFPT|` table
headed by excitation energy in eV and oscillator strength. It requires states
S1 through `EXC_NSTATES`, and reports missing states or duplicate state tables
explicitly. For logs containing several multiplicities or calculations, first
select a single matching optical table. CI coefficient squares are not optical
oscillator strengths and are never used as substitutes.

| Setting | Meaning |
| --- | --- |
| `UVVIS_STATES` | `None`: all `EXC_NSTATES` excited states; or a list of state numbers starting at 1. S0 is excluded from this excitation DOS. |
| `UVVIS_FWHM_EV` | Gaussian full width at half maximum in eV (default 0.10); `sigma = FWHM/sqrt(8 ln 2)`. |
| `UVVIS_GRID_STEP_EV` | Shared energy-grid spacing in eV (default 0.005). |
| `DOS_BINS` | Histogram bins for the selected excitation energies (default 50). |

Each frame's spectrum is a sum of unit-area Gaussians weighted by its oscillator
strengths. Spectra are averaged on one common energy grid, so the integrated
strength equals the mean sum of oscillator strengths, apart from truncated
Gaussian tails. No peak rescaling is applied to the saved optical spectrum.
This is an oscillator-strength spectrum, not an absolute absorption cross
section. Dark transitions contribute to the DOS but negligibly to UV–vis.

The DOS is a unit-area distribution of the selected **many-electron excitation
energies**, with all included states and frames weighted equally. It is not a
Kohn–Sham orbital DOS and has no occupation or oscillator-strength weighting.
Both side panels use excitation energy in eV, not wavelength. `EXC_CI_THRESHOLD`
affects composition analysis only; it does not modify the printed oscillator
strengths or the UV–vis spectrum.

For the combined optical plot, let $\epsilon_{na}$ and $f_{na}$ be CP2K's
excitation energy and oscillator strength for selected excited state $a$ at
frame $n$. On the same energy grid for every frame,

$$
G_\sigma(x)=\frac{e^{-x^2/(2\sigma^2)}}{\sqrt{2\pi}\sigma},\qquad
\sigma=\frac{\mathrm{FWHM}}{\sqrt{8\ln 2}},\qquad
A(E)=\frac{1}{N_t}\sum_{n=0}^{N_t-1}\sum_{a\in\mathcal A}
f_{na}G_\sigma(E-\epsilon_{na}).
$$

Consequently $\int A(E)\,dE=N_t^{-1}\sum_{n,a}f_{na}$ over the full energy
axis. The DOS histogram instead counts all selected excitations equally:

$$
\rho_b^{\mathrm{DOS}}=
\frac{\sum_{n,a\in\mathcal A}\mathbf{1}[\epsilon_{na}\in B_b]}
{N_t\,|\mathcal A|\,\Delta E_b}.
$$

The side panels plot $(\rho_b^{\mathrm{DOS}},E_b)$ and $(A(E),E)$, so energy
is vertical; rotating the display introduces no change in normalization.

### Gap fluctuations, autocorrelation, influence spectrum, and spectral density

Include `"spectra"` in `PLOTS` and edit `SPECTRUM_PAIRS`, for example
`[(0, 1), (1, 2)]`. This needs only the Step 3 Hamiltonians and frame spacing;
it can run before TSH. Each pair produces a four-panel PNG and a numerical
`gap_Si_Sj_spectra.npz` archive with the gap, centered fluctuations, raw and
normalized ACF, frequency grid, influence spectrum, and spectral density.

The calculation uses `deltaE = (E_j-E_i) - mean(E_j-E_i)` and the biased ACF
`C(k) = sum_t deltaE(t)*deltaE(t+k)/N`, corresponding to the statistician's
normalization. The influence spectrum is the two-sided Fourier transform of
this ACF evaluated at nonnegative angular frequencies, implemented as the
nonnegative periodogram `dt*abs(FFT(deltaE))**2/N`. It is reported in eV² fs;
it is not the squared Fourier amplitude of an already transformed ACF.

The harmonic-bath spectral-density convention is
`J(omega) = beta*omega*I(omega)/2`, with
`beta = 1/(kB*SPECTRUM_TEMPERATURE)`. All factors are evaluated in atomic units
before conversion: frequency is shown in cm⁻¹, and spectral density in eV
(the archive stores `spectral_density_ha`). Temperature is in kelvin and must
match the intended bath interpretation. This classical harmonic conversion is
an approximation; it does not identify a unique bath model from a short trace.

`SPECTRUM_MAX_WAVENUMBER` limits the displayed range, also bounded by the
Nyquist frequency. Zero padding only interpolates frequency samples; physical
resolution remains limited by the trajectory duration. The short teaching
trajectory cannot resolve a low-frequency vibrational spectrum.
Longer independent trajectories are needed; cyclic repetition of input frames
does not add information. A constant gap gives zero ACF and spectral signal.

For the spectral calculation, define the signed gap and its fluctuation as

$$
g_n=E_j(t_n)-E_i(t_n),\qquad
\bar g=\frac{1}{N_t}\sum_n g_n,\qquad \delta g_n=g_n-\bar g.
$$

The raw and normalized autocorrelations at lag $k\Delta t$ are

$$
C_k=\frac{1}{N_t}\sum_{n=0}^{N_t-k-1}\delta g_n\delta g_{n+k},\qquad
\widehat C_k=\frac{C_k}{C_0}.
$$

The denominator is $N_t$, not $N_t-k$. The script returns a zero normalized
ACF when $C_0\leq10^{-30}$ in squared hartree. With a zero-padded transform of
length $L=2N_t$, the sampled influence spectrum and spectral density are

$$
F_m=\sum_{n=0}^{N_t-1}\delta g_n e^{-2\pi i mn/L},\qquad
\omega_m=\frac{2\pi m}{L\Delta t},\qquad
I_m=\frac{\Delta t}{N_t}|F_m|^2,\qquad
J_m=\frac{\beta_T\omega_m}{2}I_m,\quad
\beta_T=\frac{1}{k_B T}.
$$

Only $0\leq m\leq L/2$ is displayed; $I_m$ is the two-sided spectral density
evaluated at nonnegative frequencies, with no extra factor of two.
Computations use hartree and atomic time. The horizontal axis converts angular
frequency to wavenumber, $\widetilde\nu=\omega/(2\pi c)$, and the outputs use
$[I]=\mathrm{eV}^2\,\mathrm{fs}$ and $[J]=\mathrm{eV}$ in the plots.

### Excitation composition and degree of mixing

Include `"excitation"` in `PLOTS` and set `EXC_LOG_DIRECTORY` to the logs
for the same system and trajectory used in Step 3. TSH populations and
Step 4 archives are not needed. Log-only selections do not read Step 3 matrices.

| Setting | Meaning |
| --- | --- |
| `EXC_LOG_PATTERN` | Log filename template, initially `"step_{step}.log"`. |
| `ISTEP`, `FSTEP` | Absolute log indices, end excluded; files are processed in numerical order. |
| `DT` | Time between consecutive frames in fs, also used for the log time axis. |
| `EXC_NSTATES` | Number of excited states to parse; excludes S0. |
| `EXC_CI_THRESHOLD` | Strict squared-coefficient cutoff, `c² > threshold`. |
| `EXC_IS_UKS` | Whether the CP2K output uses unrestricted spin labels. |
| `EXC_NSDS` | Number of printed entries and largest-weight ranks to display per state. |
| `EXC_TOP_TRANSITIONS` | Number of orbital identities displayed in each state's time map. |

For every retained expansion, the analysis defines normalized weights
`w_k = c_k²/sum(c²)`. It reports several complementary views:

- **Printed-entry averages:** mean weights of the first `EXC_NSDS` entries in
  each log's excitation expansion, with missing trailing entries padded by zero.
  Entry position need not identify the same orbital transition at every frame.
- **Ranked averages:** the largest, second-largest, etc. weights at each frame,
  independent of printed order. Ranks can also change orbital identity.
- **Orbital-transition maps:** weights of fixed `(occupied, virtual, spin)`
  labels through time, with absent transitions assigned zero. Each state's most
  important mean-weight transitions are displayed; all identities and their
  mean weights are exported to CSV. Orbital labels follow CP2K indices and do
  not by themselves track orbital character through crossings.
- **Mixing diagnostics:** retained raw squared norm, configuration count, and
  participation ratio `1/sum(w_k²)` versus time. A single configuration has
  ratio 1; `m` equally weighted configurations have ratio `m`.

Mean-weight plots include the illustrative half-width `1.96*s/sqrt(N)` with
sample standard deviation. Neighboring MD frames are correlated, so these are
naive independent-sample estimates, not reliable 95% confidence intervals for
correlated trajectories. One frame gives no uncertainty estimate. Empty
expansions cause an explicit error rather than being normalized to zero.
Raw squared norms expose the effect of truncation; normalized weights alone
can conceal it. Lowering the cutoff only recovers amplitudes present in the
printed logs. State and orbital crossings can change the apparent composition.

Outputs include `excitation_printed_weights.png/.txt`,
`excitation_ranked_weights.png/.txt`, `excitation_mixing.png`,
`excitation_Sn_transitions.png`, `excitation_transitions.csv`, and
`excitation_analysis.npz` with the time-dependent weights and diagnostics.
These coefficient-based quantities describe TD-DFT pseudowavefunction
composition, not oscillator strengths, transition-density charge transfer,
or spatial electron–hole separation.

For coefficients $c_{nak}$ of excited state $a$, frame $n$, and configuration
$k$, the retained set, raw weight, normalized weight, count, and participation
ratio are


$$
\mathcal{R}_{na} =
\left\lbrace k \mid |c_{nak}|^2 > \mathrm{ci\_threshold} \right\rbrace,
\qquad
Q_{na} = \sum_{k \in \mathcal{R}_{na}} |c_{nak}|^2,
\qquad
w_{nak} = \frac{|c_{nak}|^2}{Q_{na}}.
$$


$$
M_{na}=|\mathcal R_{na}|,\qquad
\mathrm{PR}_{na}=\frac{1}{\sum_{k\in\mathcal R_{na}}w_{nak}^2}.
$$

Printed-entry plots average $w$ by retained list position; ranked plots first
sort each frame's weights in descending order. Missing trailing positions
have weight zero. For either sequence $v_{nar}$, the plotted mean and error
bar are

$$
\bar v_{ar}=\frac{1}{N_t}\sum_n v_{nar},\qquad
s_{ar}^2=\frac{1}{N_t-1}\sum_n(v_{nar}-\bar v_{ar})^2,\qquad
h_{ar}=1.96\frac{s_{ar}}{\sqrt{N_t}}.
$$

These half-widths require $N_t>1$ and assume independent samples. For a fixed
orbital identity $q=(i\!\to\!j,\mathrm{spin})$, the transition maps use

$$
W_{naq}=\sum_{k\in\mathcal R_{na}}w_{nak}\mathbf{1}[q_{nak}=q],\qquad
\bar W_{aq}=\frac{1}{N_t}\sum_n W_{naq}.
$$

This sums duplicate entries with the same identity and assigns zero when an
identity is absent. The largest $\bar W_{aq}$ select the rows shown in each map.

The analysis conventions and layouts are illustrated in the
[time-resolved energies tutorial](https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/11_program_specific_methods/3_cp2k_methods/3_time_resolved_energies/tutorial.ipynb),
[excitation analysis tutorial](https://github.com/compchem-cybertraining/Tutorials_Libra/blob/master/11_program_specific_methods/3_cp2k_methods/4_excitation_analysis/tutorial.ipynb),
and [Fullerene spectra analysis](https://github.com/AkimovLab/Project_Fullerenes/blob/main/5_analysis/spectra_dynamics.py).
Population-based TRPES and Shannon entropy are calculated after NA-MD by
`4_step4/NAMD_plotting.py`, as described in the [Step 4 guide](../4_step4/README.md).

## Student exercise: convergence of TD-DFT NACs

**Objective:** find the least expensive active space and CI cutoff that give
stable many-electron NACs for this trajectory. Extend the tutorial's full versus
truncated time-overlap comparison to the couplings themselves. The two controls
are related but distinct: `ci_threshold` removes excitation configurations,
whereas `nelec_act_space` reduces the occupied-electron content of the underlying
determinants. Neither changes `nstates`.

1. **Fix the comparison conditions.** Use frames `0..19`, `nstates=11`, and
   `dt=1.0*units.fs2au`. Apply the time-indexing adjustment above and exclude
   initialization output `0`. Keep the same logs, state ordering, and KS files
   across runs. Start with the larger `1..40` overlap window and
   `lowest_orbital=1`; first verify that it covers the retained excitations.
2. **Audit the excitation expansions.** For each cutoff in
   `[1e-2, 1e-3, 1e-4, 1e-5]`, parse every log with
   `read_cp2k_tddfpt_log_file`, passing `ci_threshold` and
   `number_of_states=10` (the parser counts excited states only). Record the
   number of retained configurations and `sum(c*c for c in coefficients)` for
   each state and frame, plus the smallest occupied and largest virtual indices.
   The sum is a retained printed-coefficient weight, not a guaranteed normalized
   probability. Identify empty states and orbital-coverage failures before
   running the overlaps. From the occupied bound above, determine which active
   electron counts are admissible at each cutoff.
3. **Converge the CI cutoff using all available orbitals.** Set
   `nelec_act_space=None` and lower the cutoff through the proposed values.
   If a retained transition exceeds `1..40`, enlarge the Step 2 window and
   repeat the comparisons using that common window. Choose the smallest
   usable cutoff as a provisional reference; do not assume `1e-5` is converged
   without comparison to a smaller cutoff and adequate CP2K amplitude printing.
4. **Converge the active space at each cutoff.** Try active electron counts
   `[2, 4, 8, 12, 22, 32, 48, None]`, retaining only combinations that cover
   all necessary occupied and virtual orbitals at both endpoints. Mark rejected
   combinations as invalid, not as zero NACs. Increase the space as the cutoff
   decreases. Compare each valid integer choice to `None` at the same cutoff,
   then compare to the provisional reference from step 3. With HOMO 24,
   `48` and `None` include the same occupied orbitals, but the integer setting
   ends its virtual range at the current `max_vir`; `None` uses the full stored
   window. Account for this distinction when diagnosing differences.
5. **Keep runs independent.** Create a fresh `params` dictionary for every
   parameter pair, including fresh `{0: True}` initialization flags; do not
   carry `data_prev` or `s_ci_inv_prev` across runs. Save each run in a separate
   directory, for example `ci_N22_tol1e-3` or `ci_all_tol1e-5`, and record runtime.
6. **Compare NACs quantitatively.** Recover the real time-derivative coupling
   matrix from each saved `st_adi_i.txt` using the code below. Plot several
   representative couplings versus time, including a strongly coupled pair
   and a weakly coupled pair. Report maximum absolute and RMS errors over all
   physical intervals and unique off-diagonal state pairs. Compare coupling
   magnitudes if phases differ; signed errors require consistent phases and
   state labels. Magnitudes alone do not resolve state permutations or rotations
   within nearly degenerate subspaces.

```python
from pathlib import Path
import numpy as np
from libra_py import units

def load_nacs(directory, nframes=20, dt=1.0 * units.fs2au):
    overlaps = np.stack([
        np.loadtxt(Path(directory) / f"st_adi_{i}.txt")
        for i in range(1, nframes)  # omit initialization
    ])
    return (overlaps - overlaps.transpose(0, 2, 1)) / (2.0 * dt)

# Replace these directory names with two completed, compatible runs.
d_ref = load_nacs("ci_all_tol1e-5")
d_test = load_nacs("ci_N22_tol1e-3")
row, col = np.triu_indices(d_ref.shape[1], k=1)
error = (np.abs(d_test) - np.abs(d_ref))[:, row, col]
print("Maximum magnitude error [a.u. of inverse time]:", np.max(np.abs(error)))
print("RMS magnitude error [a.u. of inverse time]:", np.sqrt(np.mean(error**2)))
```

For the exercise, let $\mathcal P$ contain the unique off-diagonal pairs and
$K$ the number of physical intervals. The magnitude-error metrics in the code
are

$$
e_{n,ab}=|d^{\mathrm{test}}_{ab}(t_n)|-|d^{\mathrm{ref}}_{ab}(t_n)|,\qquad
\epsilon_{\max}=\max_{n,(a,b)\in\mathcal P}|e_{n,ab}|,\qquad
\epsilon_{\mathrm{RMS}}=
\sqrt{\frac{1}{K|\mathcal P|}\sum_{n,(a,b)\in\mathcal P}e_{n,ab}^2}.
$$

For consistently aligned phases and labels, replace $e$ by the signed coupling
difference to assess signed errors.

These `d_ab` values have units of inverse atomic time. The corresponding
off-diagonal coupling energy in eV is `hbar*d_ab = d_ab/units.ev2Ha`
numerically, since `hbar=1` in atomic units; `H_vib,ab=-1j*d_ab` in a.u.
Do not use `obj.nac_adi` for this exercise: it is a zero-filled placeholder in
the inspected implementation.

**Deliverables:** a table of cutoff, active-electron count, actual orbital
range, retained configuration counts/weights, runtime, maximum NAC error, and
RMS NAC error; coupling-versus-time plots; and a justified accuracy/cost choice.
Choose and state an absolute error tolerance before drawing conclusions. As an
illustrative target, require maximum coupling-energy magnitude error below
1 meV; tighter tolerances may be needed for weak couplings. Avoid relative
errors for reference couplings near zero. Check `s_adi` against identity, but
do not use orthonormality alone as evidence that NACs are converged.

As an optional extension, use the Step 4 SVD basis check for different
leading-state subsets. Small singular values can indicate leakage from the
chosen state subspace or state mixing across its boundary. This probes the
choice of electronic states, which is distinct from the orbital active space;
passing the rank test does not prove NAC convergence with `ci_threshold`.

## Libra source references

Paths below are relative to the `src` directory of a Libra source checkout,
wherever it is installed on your system:

- `libra_py/packages/cp2k/methods.py`: `read_cp2k_tddfpt_log_file` implements
  coefficient filtering and energy conversion; `cp2k_compute_adi` implements
  active-space selection, caching, orthogonalization, and Hamiltonians;
  `read_trajectory_xyz_file` extracts coordinates.
- `libra_py/citools/ci.py`: `overlap` selects `nstates-1` excited states,
  determines the effective electron count, and contracts CSF overlaps.
- `libra_py/citools/interfaces.py`: `ci_amplitudes_mtx` builds coefficient
  matrices; `sd_and_csf_overlaps` validates orbital ranges and reference
  excitations and constructs the spin-adapted overlaps.
- `libra_py/orthogonalizations.py`: `lowdin_inverse_sqrt` implements symmetric
  orthogonalization and its eigenvalue cutoff.
- `libra_py/workflows/nbra/step2.py`: the KS-overlap writer saves current-frame
  `S_ks_{step}` and previous-to-current `St_ks_{step-1}` files.

Consult these implementations when using another Libra version: some source
docstrings are less precise than the code (for example, the parser already
converts energies to hartree, and `None` uses the stored orbital window).
