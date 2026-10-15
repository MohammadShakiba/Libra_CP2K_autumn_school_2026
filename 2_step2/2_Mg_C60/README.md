# Step 2: Electronic-structure calculations and orbital overlaps for Mg@C60

In this tutorial, we use the molecular-dynamics (MD) trajectory of **Mg@C60**
from [Step 1](../../1_step1/2_Mg_C60/) to calculate electronic energies,
TD-DFPT excitation data, and Kohn–Sham (KS) orbital overlaps with CP2K and
Libra. These quantities are the inputs for the many-electron overlaps and
vibronic Hamiltonians constructed in [Step 3](../../3_step3/README.md).

During the tutorial session, we will align the trajectory and process only
**frames 0, 1, and 2**, giving two adjacent-frame intervals. We will then
unpack the supplied reference results to continue to Step 3 without waiting
for the full calculation.

## Files in this folder

| File | Purpose | What you need to do |
| --- | --- | --- |
| [README.md](README.md) | Instructions for the Mg@C60 Step 2 workflow. | Follow the steps below in order. |
| [align.py](align.py) | Removes whole-molecule translation and rotation from the MD coordinates and writes an aligned XYZ trajectory. | Set the input trajectory path and run this script first. |
| [run_template.py](run_template.py) | Defines the CP2K/Libra settings and calls `step2.run_cp2k_libint_step2`. The distributor turns it into `run.py` inside each job directory. | Set the aligned trajectory, executable, processor count, and orbital window. Leave `istep` and `fstep` blank. |
| [distribute_jobs.py](distribute_jobs.py) | Divides the requested frame range among jobs, fills in each job's frame limits, and launches the submission script. | Set `istep`, `fstep`, `njobs`, and the submission options. |
| [submit_template.slm](submit_template.slm) | SLURM resource requests, NESH modules, and the command that runs each generated `run.py`. | Check the reservation, resources, and environment setup. Match the CPU allocation to `nprocs`. |
| [es_diag_temp.inp](es_diag_temp.inp) | CP2K single-point template for the DFT/TD-DFPT calculation. Requests the Molden and PDOS output used by this workflow. | Review the Mg@C60 model and keep the required output sections enabled. |
| [es_ot_temp.inp](es_ot_temp.inp) | Alternative CP2K input using orbital-transformation SCF. The Libra workflow also provides an OT-template setting for its xTB route. | Keep this file available. With the supplied `isxTB=False` setting, this workflow uses `es_diag_temp.inp`. |
| [reference.tar.bz2](reference.tar.bz2) | Precomputed Mg@C60 results for frames 0–300. Contains `res/`, `all_logfiles/`, and `all_pdosfiles/`. | Unpack after the short calculation to obtain the data for Step 3. |
| [clean.sh](clean.sh) | Removes `job*`, `all_*`, and `res` from the current directory. | Use only when you intend to delete those calculation files. Save any results you need first. |

The MD trajectory comes from Step 1; it is not included in this folder's
reference archive. The single-geometry file `C60_Mg.xyz` in Step 1 is the
initial structure, whereas `C60_Mg-pos-1.xyz` is the multi-frame MD trajectory.

## 1. Load Libra on NESH

From the repository root, enter the Mg@C60 Step 2 folder:

```bash
cd 2_step2/2_Mg_C60
```

Load the tutorial Libra installation:

```bash
export PYTHONPATH=/gxfs_home/cau/suphc424/Software/libra-code/_build/src:$PYTHONPATH
source /gxfs_home/cau/suphc424/.conda/envs/libra/bin/activate
```

Check that Python can import Libra:

```bash
python -c "import liblibra_core; from libra_py.workflows.nbra import step2; print('Libra is ready')"
```

These paths refer to the NESH tutorial installation. Adapt them if you use
a different Libra build or Python environment, and use the same environment
in the batch script below. Run the commands in this guide from
`2_step2/2_Mg_C60` unless another directory is specified.

## 2. Align the Mg@C60 trajectory

Before computing the orbital overlaps, remove translation and rotation of
the molecule as a whole. Alignment retains the internal molecular motion
and preserves the atom and frame order.

Open `align.py`. Replace its hard-coded input path with the location of your
Mg@C60 MD trajectory. If the trajectory is in the Step 1 folder, use:

```python
from libra_py.md_align import align_trajectory

align_trajectory(
    "../../1_step1/2_Mg_C60/C60_Mg-pos-1.xyz",
    "C60_Mg_MD-pos-1-aligned.xyz",
)
```

If you use the precomputed Step 1 trajectory, unpack the Step 1 results
first so that the input XYZ file exists. An absolute path is also suitable
if your MD results are elsewhere.

Run the alignment:

```bash
python align.py
```

The output is `C60_Mg_MD-pos-1-aligned.xyz` in the current folder. Check
that it contains the expected 61 atoms per frame: 60 carbon atoms and one
magnesium atom. We will use this aligned file for the short calculation.

## 3. Configure the CP2K/Libra calculation

Open `run_template.py` and change the trajectory setting to the aligned file:

```python
params['trajectory_xyz_filename'] = path + '/../C60_Mg_MD-pos-1-aligned.xyz'
```

Here, `path = os.getcwd()` is evaluated inside a generated directory such
as `job1/`, so `../` points back to the Mg@C60 Step 2 folder. An absolute
path to the aligned trajectory is equally valid. Update any other
installation-specific paths in the template as needed.

Review these settings:

| Setting in `run_template.py` | Value for this tutorial | Meaning |
| --- | --- | --- |
| `params['nprocs']` | `16` | Number of MPI processes used to launch CP2K. |
| `params['mpi_executable']` | `'mpirun'` | MPI launcher; must match the CP2K module loaded in the batch job. |
| `params['cp2k_exe']` | `'cp2k.psmp'` | CP2K executable supplied by the loaded module. |
| `params['lowest_orbital']` | `90` | First spatial orbital saved in the overlap matrices, using one-based CP2K numbering. |
| `params['highest_orbital']` | `160` | Last spatial orbital saved; this upper bound is inclusive. |
| `params['isxTB']` | `False` | Use the DFT calculation with `es_diag_temp.inp`. |
| `params['isUKS']` | `False` | Restricted, closed-shell calculation. |
| `params['is_periodic']` | `False` | Treat Mg@C60 as an isolated molecule. The inactive periodic-cell example in the Python template is not used. |
| `params['cube_visualization']` | `False` | Cube rendering and VMD are not needed. |
| `params['restart_file']` | `False` for a fresh run without an external wavefunction | Disables the optional external restart-file copy. |

The supplied template currently has `params['restart_file'] = True`. For
this fresh tutorial run, change it to `False`: Libra versions that interpret
this parameter as a file path otherwise attempt to copy a file named `True`.
This setting is separate from CP2K's reuse of a preceding frame's wavefunction
when one is available.

**Leave these two assignments blank in `run_template.py`:**

```python
params['istep'] =
params['fstep'] =
```

`run_template.py` is a template and cannot be executed directly while those
values are blank. Set the overall frame range in `distribute_jobs.py`; Libra
will insert the appropriate values into each generated `job*/run.py`.

### Mg@C60 electronic-structure settings

The supplied `es_diag_temp.inp` uses PBE-D3, a nonperiodic cell of about
21.11 Å per side, and the Martyna–Tuckerman Poisson solver. It uses
`DZVP-MOLOPT-GTH`/`GTH-PBE-q4` for C and
`DZVP-MOLOPT-SR-GTH`/`GTH-PBE-q10` for Mg.

For this neutral, restricted system, the reference CP2K log reports
**250 explicitly treated electrons and 125 occupied spatial orbitals**.
Thus, the HOMO is orbital **125** and the LUMO is **126**; the old
`Here HOMO is 90` comment in the Python template is incorrect for this input.
The saved window `90..160` contains 71 spatial orbitals. Libra stores both
spin blocks, giving **142 × 142** KS matrices, as in the reference dataset.

The CP2K template requests `ADDED_MOS 60` and `TDDFPT/NSTATES 10`.
`ADDED_MOS` controls extra KS orbitals available to CP2K, the Python orbital
window selects the orbitals saved by Libra, and `NSTATES` controls the number
of excited states. Keep the orbital window wide enough to contain the
transitions retained in Step 3.

Libra generates a coordinate file and CP2K input for each frame, replacing
the template's `COORD_FILE_NAME` automatically. Select the multi-frame
trajectory through `trajectory_xyz_filename` in the Python template. Keep
the Molden and PDOS print sections enabled; the workflow reads the Molden
orbitals to calculate overlaps.

## 4. Check the NESH submission template

Open `submit_template.slm` and review the SLURM requests and module commands.
For the supplied 16-process, single-node setup, the relevant contents are:

```bash
#!/bin/bash -l
#SBATCH --reservation=autumnschool
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000

module load gcc12-env
module load cp2k-autumnschool/default
module load cp2k

export OMP_NUM_THREADS=1
export PYTHONPATH=/gxfs_home/cau/suphc424/Software/libra-code/_build/src:$PYTHONPATH
source /gxfs_home/cau/suphc424/.conda/envs/libra/bin/activate

python run.py
```

Add the Libra setup commands before `python run.py` in the supplied submission
template so that each batch job has the required environment. Keep the
`autumnschool` reservation only when it is available for your session; adapt
the account, partition, reservation, paths, and resource limits to your NESH
allocation as directed by the instructors. The commented UB-cluster examples
in the supplied file are unrelated to this NESH run.

**Keep the processor settings consistent:** this example allocates
`1 × 16 = 16` MPI tasks, matching `params['nprocs'] = 16`. With
`--cpus-per-task=1`, use `OMP_NUM_THREADS=1`. If you change the parallel
layout, update the SLURM allocation, MPI process count, and thread count
together.

## 5. Select a short frame range and submit

Open `distribute_jobs.py`. The repository settings request frames 0–300
across 50 jobs. For the tutorial session, replace them with:

```python
run_slurm = True
submit_template = 'submit_template.slm'
run_python_file = 'run_template.py'
istep = 0
fstep = 2
njobs = 1
submission_exe = 'sbatch'
```

Frame numbering starts at zero. The distributor treats `fstep` as the final
included frame, so this example processes frames **0, 1, and 2** and produces
time overlaps for **0 → 1** and **1 → 2**. In `job1/run.py`, Libra writes
`istep=0` and `fstep=3`, because the generated calculation loop excludes
its upper bound. When multiple jobs are requested, neighboring jobs share a
boundary frame so that the adjacent-frame overlaps can be computed.

Submit from the Mg@C60 Step 2 folder:

```bash
python distribute_jobs.py
```

This creates `job1/`, writes its `run.py`, copies `submit_template.slm`, and
calls `sbatch` automatically. Each batch job runs the generated `run.py`,
which launches CP2K for its assigned frames and uses Libra to compute the
orbital overlaps.

**Keep `run_slurm=True` on NESH.** In Libra's distributor, setting it to
`False` executes the submission script locally with `sh`; it is not a mode
that only prepares job files. Rerunning the distributor also recreates
matching `job*` directories, so preserve a completed run before reusing them.

Monitor your job with:

```bash
squeue -u "$USER"
```

After it finishes, inspect `job1/slurm-*.out` and the CP2K logs. Check that
both the SCF and TD-DFPT calculations converge before using the outputs.

## 6. Inspect the short-run results

The calculation creates the following directories in `2_step2/2_Mg_C60`:

| Output | Contents |
| --- | --- |
| `job1/`, `job2/`, … | Generated `run.py`, submission script, per-frame coordinates and CP2K inputs, SLURM output, and working files. |
| `all_logfiles/step_i.log` | CP2K output for frame `i`, including TD-DFPT excitation energies and amplitudes needed in Step 3. |
| `all_pdosfiles/` | Projected density-of-states output from CP2K. |
| `res/E_ks_i.npz` | KS orbital-energy matrices. |
| `res/S_ks_i.npz` | Same-frame KS overlaps, comparing orbitals at frame `i` with each other. |
| `res/St_ks_i.npz` | Time overlaps, comparing orbitals at frame `i` with those at frame `i+1`. These are inputs to the later coupling calculation. |

For the three-frame example, expect three logs and three each of `E_ks`
and `S_ks`, indexed 0–2, plus two `St_ks` files, indexed 0–1.
An empty `all_images/` directory may also be created even though cube
visualization is disabled. Molden files are removed after use when
`remove_molden=True`.

Read the sparse matrices with SciPy:

```bash
python - <<'PY'
from scipy.sparse import load_npz

S = load_npz("res/S_ks_0.npz").toarray()
St = load_npz("res/St_ks_0.npz").toarray()
print("Same-frame overlap:", S.shape)
print("Time overlap:", St.shape)
PY
```

Both shapes should be `(142, 142)` for the `90..160` window. A same-frame
overlap should be close to the identity matrix. A time-overlap matrix
need not be identity or symmetric, and orbital signs can differ between
frames.

## 7. Unpack the reference results

Once you have followed the short calculation, use the supplied archive to
continue the tutorial with the complete electronic dataset. Extract into a
separate `reference/` folder to keep the short-run results available:

```bash
mkdir -p reference
tar xvf reference.tar.bz2 -C reference
```

The archive has `res/`, `all_logfiles/`, and `all_pdosfiles/` at its root;
the command above places them under `reference/`.

| Reference data | Available indices |
| --- | --- |
| `reference/all_logfiles/step_i.log` | Frames 0–300: 301 CP2K logs. |
| `reference/res/E_ks_i.npz` | Frames 0–300: 301 KS energy matrices. |
| `reference/res/S_ks_i.npz` | Frames 0–300: 301 same-frame overlap matrices. |
| `reference/res/St_ks_i.npz` | Intervals 0 → 1 through 299 → 300: 300 time-overlap matrices. |
| `reference/all_pdosfiles/` | PDOS files for the Mg and C kinds at each frame. |

The archive contains electronic results only; keep the corresponding
Mg@C60 trajectory available for Step 3. Use logs, overlaps, and coordinates
from a matching dataset and alignment convention.

## 8. Continue to Step 3

Use the [Step 3 guide](../../3_step3/README.md) with the Mg@C60 data.
The supplied Step 3 `run.py` still contains TiO2 example paths and orbital
settings, so update them for this system before running it:

| Step 3 setting | Mg@C60 choice |
| --- | --- |
| Trajectory | The Mg@C60 XYZ trajectory corresponding to the selected electronic dataset; use the aligned file for your aligned short run. |
| Log directory | `reference/all_logfiles/` for the extracted results, using a path valid from the Step 3 working directory. |
| KS matrix directory | `reference/res/` for the extracted results, using a path valid from the Step 3 working directory. |
| `lowest_orbital` | `90` for the supplied reference window `90..160`. |
| `nstates` | `11` to include the ground state and all 10 calculated excited states. |
| `dt` | Actual time between saved MD frames. The supplied Step 1 input uses a 1 fs MD step and saves every step, giving `1.0 * units.fs2au`. |
| Frame range | A subset covered by the selected logs and matrices; the reference contains frames 0–300. |
| Output directory | Use a name appropriate for Mg@C60, such as `ci_mg_c60`. |

For an interval `i → i+1`, `St_ks_i.npz` must be paired with electronic
states from `step_i.log` and `step_{i+1}.log`. Follow the
[time-indexing instructions](../../3_step3/README.md#time-indexing-align-both-endpoints)
in the Step 3 guide when adapting its supplied loop. Keep the same orbital
window, trajectory order, and saved-frame spacing throughout the calculation.
