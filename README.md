# Libra–CP2K Autumn School 2026

Hands-on examples for a four-step nonadiabatic dynamics workflow with CP2K and
Libra: generate nuclear motion, compute electronic structure along the
trajectory, construct state overlaps and couplings, then propagate electronic
populations. Steps 1–3 use a rutile TiO₂ unit cell; Step 4 also includes
separate example datasets.

## Four-step workflow

| Step | Task | Main scripts | Output passed onward |
| --- | --- | --- | --- |
| [1. Molecular dynamics](1_step1/README.md) | Run ground-state CP2K MD to generate nuclear geometries. | `convert.py`, `md.inp`, `submit.slm` | `Rutile_TiO2_MD-pos-1.xyz`, with the cell and saved-frame spacing. |
| [2. Electronic structure and KS overlaps](2_step2/README.md) | Run CP2K at each geometry; compute KS same-frame/time overlaps and TD-DFPT excitations. | `distribute_jobs.py`, `run_template.py`, `es_diag_temp.inp` | `res/` KS matrices and `all_logfiles/` excitation energies/amplitudes, paired with the trajectory. |
| [3. Many-electron overlaps and couplings](3_step3/README.md) | Combine KS overlaps with excitation amplitudes to construct TD-DFT state overlaps and vibronic Hamiltonians. | `run.py`, `plotting.py` | `ci_tio2/` state energies, overlaps, and Hamiltonians; electronic-structure analysis plots. |
| [4. Nonadiabatic dynamics](4_step4/README.md) | Check the electronic basis, choose a hopping/decoherence recipe, and analyze populations. | `NAMD_simulation.py`, `NAMD_plotting.py` | Basis diagnostics, saved population dynamics, figures, and fit tables. |

Step 4 uses the neglect of back-reaction approximation (NBRA): the electronic
populations evolve along prescribed nuclear motion. The electronic dynamics
does not update the nuclear trajectory from Step 1.

## Getting started

Run each script from its own step directory and follow that directory's README.
Steps 1 and 2 require CP2K and use example SLURM scripts that must be adapted
to your cluster. The Python workflow requires Libra (`libra_py` and
`liblibra_core`), NumPy, and SciPy, including the overlap support used in Step 2.
Plotting and dynamics analysis also use Matplotlib and h5py. Load the workshop
environment or your local installation before starting.

To explore without submitting CP2K calculations, use `1_step1/reference/`,
`2_step2/input/`, and `2_step2/reference/`. The Step 1 and Step 2 guides include
student exercises based on these saved outputs. Step 3 contains sample inputs;
Step 4 provides three dataset archives and instructions for selecting them.

## Connecting the steps

Keep the geometry, atom order, cell, orbital window, and frame indices consistent
as you transfer data. The short Step 1 example saves 21 frames at 1 fs spacing,
yielding 20 adjacent-frame intervals in Step 2. It demonstrates the workflow;
it is too short to establish equilibration or reliable relaxation times.

Before constructing couplings, follow the time-indexing correction documented
in the [Step 3 guide](3_step3/README.md). Its supplied loop requires that
adjustment to pair the proper TD-DFT logs with each KS time overlap.

Step 4 initially selects the separate Mg@C60 dataset (`data3/res`, 0.5 fs), so
it does not automatically consume the TiO₂ results. To use your Step 3 output,
select the text-matrix reader, point `data_path` to `../3_step3/ci_tio2`, and
set the frame range, state selection, and `DT` consistently. Exclude the
initialization output described in Step 3. Follow the Step 4 guide to enable
the basis-check stop before continuing to dynamics and plotting.
