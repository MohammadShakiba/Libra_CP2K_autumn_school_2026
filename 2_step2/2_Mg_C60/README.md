# Step 2: KS overlaps and TD-DFT excitation data

For each saved geometry from [Step 1](../1_step1/README.md), run CP2K
single-point calculations and use Libra to compute Kohn–Sham (KS) orbital
overlaps. CP2K also computes TD-DFPT excitation energies and amplitudes.
[Step 3](../3_step3/README.md) combines these outputs into many-electron
state overlaps and vibronic Hamiltonians.

Two overlap types are needed: `S_ks_i` compares orbitals at the same frame,
while `St_ks_i` compares frame `i` with frame `i+1`. The latter describes how
the orbital basis changes as the nuclei move; it is not itself a coupling.

## Files and settings

| File or folder | Purpose |
| --- | --- |
| `input/` | Precomputed Step 1 trajectory, energies, velocities, forces, and log. |
| `distribute_jobs.py` | Splits the trajectory, prepares job directories, and submits SLURM jobs. |
| `submit_template.slm` | Resource requests and environment setup for each job. |
| `run_template.py` | Libra settings and the call to `step2.run_cp2k_libint_step2`. Becomes `run.py` inside each job. |
| `es_diag_temp.inp` | CP2K energy/TD-DFPT template with orbital, Molden, and projected density-of-states output. |
| `align.py` | Optional trajectory translation/rotation alignment. |
| `reference/` | Completed calculation for spatial orbitals 1–40, including generated jobs. |
| `smaller_orbital_space/` | Alternative calculation retaining spatial orbitals 14–35. |
| `clean.sh` | Deletes `job*`, `all_*`, and `res` in the current directory; preserve needed results before using it. |

The control sequence is:

```text
distribute_jobs.py → submit_template.slm → job*/run.py → CP2K + Libra overlaps
```

The supplied settings use the same periodic cell, PBE-D3 treatment, basis sets,
and pseudopotentials as Step 1. `es_diag_temp.inp` requests 10 excited states
with `TDDFPT/NSTATES`. The example has 48 explicitly treated electrons and
24 occupied spatial orbitals, so the HOMO is orbital 24 and the LUMO is 25.

| Where to edit | Main choices |
| --- | --- |
| `distribute_jobs.py` | `istep=0`, `fstep=20`, `njobs=4`, and `run_slurm=True`. |
| `run_template.py` | Inclusive, one-based orbital window `1..40`; `nprocs=12`; CP2K executable and MPI launcher. |
| `run_template.py` and `es_diag_temp.inp` | Matching cell vectors and XYZ periodicity. |
| `es_diag_temp.inp` | SCF convergence, `ADDED_MOS=50`, TD-DFPT state count and convergence. |
| `submit_template.slm` | Cluster account, modules, wall time, memory, and CPU allocation. |

`ADDED_MOS` controls extra orbitals in the CP2K calculation; the Libra orbital
window controls which orbitals enter the saved overlap matrices. Neither is
the number of TD-DFPT excited states. Keep the window wide enough to represent
the excitations you will retain in Step 3.

## How to run

1. Enter this folder from the repository root and activate a Python environment
   with Libra, its compiled core/Libint overlap support, NumPy, and SciPy.
   CP2K and the cluster launcher must also be available in the batch jobs.
   Use your Step 1 trajectory, or copy the supplied example:

   ```bash
   cd 2_step2
   cp input/Rutile_TiO2_MD-pos-1.xyz .
   ```

   Skip the copy if your own trajectory is already in place. Preserve atom
   order and the cell definition throughout the workflow.

2. Review alignment only if needed. `python align.py` writes
   `Rutile_TiO2_MD-pos-1-aligned.xyz`, but the run template still reads the
   **original** filename. To use an aligned trajectory, change
   `trajectory_xyz_filename` explicitly. For this fixed periodic cell, keep
   the original trajectory unless you have checked that any coordinate
   transformation remains consistent with the lattice.

3. Edit the four configuration files in the table above. Match `nprocs` to the
   batch allocation and replace cluster-specific module paths. Leave the blank
   `istep` and `fstep` assignments in `run_template.py` for the distributor to
   fill; the template cannot be run directly as Python. Cube visualization is
   disabled, so VMD is not needed for the supplied calculation.

4. Submit from `2_step2`:

   ```bash
   python distribute_jobs.py
   ```

   With `run_slurm=True`, this submits jobs immediately. Set it to `False`
   to prepare and inspect the generated jobs without submitting. The reference
   jobs cover frames 0–5, 5–10, 10–15, and 15–20, sharing boundary frames so
   adjacent-frame overlaps are available. Their generated `fstep` values are
   exclusive (6, 11, 16, 21); inspect `job*/run.py` after changing the split.

5. Check each job's SLURM output and the collected CP2K logs for SCF and
   TD-DFPT convergence. Confirm that the expected files exist across the whole
   frame range, including job boundaries. Preserve successful outputs before
   rerunning, because matching names may be overwritten.

## Results and handoff

| Output | Meaning |
| --- | --- |
| `all_logfiles/step_i.log` | CP2K results for frame `i`, including excitation energies and amplitudes needed by Step 3. |
| `all_pdosfiles/` | Projected density-of-states files for electronic-structure analysis. |
| `res/E_ks_i.npz` | KS orbital energy matrices. Step 3's supplied `run.py` instead obtains excitation energies from the logs. |
| `res/S_ks_i.npz` | Same-frame KS overlaps. |
| `res/St_ks_i.npz` | KS time overlaps for interval `i → i+1`. |
| `job1/`, `job2/`, … | Generated inputs, scripts, and working files for inspection or debugging. |

The reference has 21 logs and 21 each of `E_ks` and `S_ks` (indices 0–20),
plus 20 `St_ks` files (indices 0–19). The matrices include alpha and beta spin
blocks: 40 spatial orbitals give an 80 × 80 matrix, while the 14–35 window
contains 22 spatial orbitals and gives a 44 × 44 matrix.

For a quick inspection from `2_step2`:

```python
from scipy.sparse import load_npz

S = load_npz("reference/res/S_ks_0.npz").toarray()
St = load_npz("reference/res/St_ks_0.npz").toarray()
print(S.shape, St.shape)
```

Pass the matching trajectory, `res/`, and `all_logfiles/` to Step 3. Keep each
dataset together: do not mix logs, orbital windows, or trajectories. Set
Step 3's `lowest_orbital` to 1 for the main example or 14 for the smaller
window, and use the actual saved-frame spacing. Read the **Time indexing:
align both endpoints** section of the [Step 3 guide](../3_step3/README.md)
before running its supplied loop; it documents a correction needed to match
the KS time-overlap endpoints with the TD-DFT logs.

## Student exercises

Use `reference/` and `smaller_orbital_space/` for the first three exercises;
no new CP2K runs are required. Submit a few plots or tables with short explanations.

1. **Trace one job.** Inspect `reference/job2/run.py` and its CP2K inputs.
   Identify its frame range, orbital window, and output locations. Explain why
   neighboring jobs share a frame, and count the distinct time intervals
   covered by all four jobs.
2. **Inspect overlap matrices.** Load `S_ks_0` and `St_ks_0`, report their
   dimensions, and plot their absolute values. Measure the largest absolute
   entry of `S - I` and of `S - S.conj().T`. Explain why a same-frame overlap
   should be close to identity, whereas a time overlap need not be identity
   or symmetric. Orbital signs can also change between frames.
3. **Connect orbitals to excitations.** In `reference/all_logfiles/step_0.log`,
   find the electron count, the first excitation energy, and its largest
   printed occupied-to-virtual contributions. Check whether those transitions
   fit in both orbital windows. Compare matrix dimensions and disk usage for
   the two datasets. Does fitting these few transitions establish convergence
   for every excited state and frame?
4. **Optional short calculation.** In a separate working copy, process a few
   adjacent frames with one job, then change only the orbital window. Compare
   cost and overlap blocks for the orbitals common to both runs. Explain what
   must be checked in Step 3 before using the smaller window for dynamics.
