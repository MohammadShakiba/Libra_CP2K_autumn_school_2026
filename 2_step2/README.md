# Step 2: CP2K electronic-structure calculations and Libra orbital overlaps

This tutorial covers both **periodic rutile TiO₂** and **molecular Mg@C60**.
Starting from the MD trajectories produced in [Step 1](../1_step1/), we run
CP2K single-point calculations and use Libra to calculate Kohn–Sham (KS)
orbital overlaps. CP2K also produces TD-DFPT excitation energies and
amplitudes. These data are used in [Step 3](../3_step3/README.md) to construct
many-electron overlaps and vibronic Hamiltonians.

During the tutorial session, process only **frames 0, 1, and 2**, giving
two adjacent-frame intervals. Then use the supplied reference results to
continue to Step 3 without waiting for the full calculation.

## Choose your example

| Setting | Rutile TiO₂ | Mg@C60 |
| --- | --- | --- |
| Working folder within `2_step2/` | [1_Rutile_TiO2/](1_Rutile_TiO2/) | [2_Mg_C60/](2_Mg_C60/) |
| System | Periodic crystal; six atoms in the cell | Isolated molecule; 61 atoms |
| Trajectory alignment | **Skip alignment; use the original trajectory** | **Run `align.py` first** |
| `params['is_periodic']` | `True` | `False` |
| `params['nprocs']` | `8` | `16` |
| Saved spatial-orbital window | `1..40`, inclusive | `90..160`, inclusive |
| HOMO / LUMO indices | `24 / 25` | `125 / 126` |
| Saved KS matrix dimensions | `80 × 80` | `142 × 142` |
| Reference results | Already unpacked in `reference/` | Unpack `reference.tar.bz2` |

Orbital indices are one-based CP2K indices. Frame indices start at zero.
The matrix dimensions include both spin blocks.

## Files in the example folders

Each example has its own scripts and input templates. Edit the files inside
the folder for the case you are running.

| File or folder | Purpose |
| --- | --- |
| `run_template.py` | CP2K/Libra settings and the call to `step2.run_cp2k_libint_step2`. The distributor turns it into `run.py` inside each job directory. Edit paths, processor count, orbital window, and periodic settings here. |
| `distribute_jobs.py` | Selects the frame range and number of jobs, generates their scripts, and launches the submission template. Set `istep`, `fstep`, and `njobs` here. |
| `submit_template.slm` | SLURM resource requests, NESH modules, environment setup, and `python run.py`. Match its allocation to `nprocs`. |
| `es_diag_temp.inp` | CP2K DFT/TD-DFPT single-point template, including the cell, basis sets, pseudopotentials, and Molden/PDOS output. |
| `align.py` | Removes whole-molecule translation and rotation from an MD trajectory. Use for Mg@C60; skip for periodic TiO₂. |
| `clean.sh` | Removes `job*`, `all_*`, and `res` in the current directory. Save needed results before using it. |
| [1_Rutile_TiO2/reference/](1_Rutile_TiO2/reference/) | Completed TiO₂ calculation, including logs, KS matrices, PDOS, and generated jobs. |
| [2_Mg_C60/reference.tar.bz2](2_Mg_C60/reference.tar.bz2) | Mg@C60 reference logs, KS matrices, and PDOS. It contains electronic results only, not the MD trajectory. |
| [2_Mg_C60/es_ot_temp.inp](2_Mg_C60/es_ot_temp.inp) | Alternative CP2K input using orbital-transformation SCF. Both examples use `es_diag_temp.inp` with the supplied `isxTB=False` setting. |

The trajectories come from the corresponding Step 1 examples. Single-geometry
files such as `TiO2_unit_cell.xyz` and `C60_Mg.xyz` are initial structures;
use the multi-frame MD position trajectories for Step 2.

## 1. Enter an example folder and load Libra on NESH

From the repository root, choose **one** working folder.

For TiO₂:

```bash
cd 2_step2/1_Rutile_TiO2
```

For Mg@C60:

```bash
cd 2_step2/2_Mg_C60
```

**Run all subsequent shell commands from your chosen example folder** unless
another directory is specified. The links in this README are relative to
`2_step2/`; shell commands below assume you have entered the example folder.

Load the tutorial Libra installation:

```bash
export PYTHONPATH=/gxfs_home/cau/suphc424/Software/libra-code/_build/src:$PYTHONPATH
source /gxfs_home/cau/suphc424/.conda/envs/libra/bin/activate
```

Check the Python environment:

```bash
python -c "import liblibra_core; from libra_py.workflows.nbra import step2; print('Libra is ready')"
```

Adapt these installation paths as needed, and use the same environment in
the batch script below.

## 2. Prepare the trajectory for your case

### TiO₂: use the original periodic trajectory

Copy the supplied trajectory from [Step 1](../1_step1/1_Rutile_TiO2/README.md):

```bash
cp ../../1_step1/1_Rutile_TiO2/reference/Rutile_TiO2_MD-pos-1.xyz .
```

If you ran Step 1 yourself, use your own trajectory instead. Preserve its
atom order, frame order, and cell definition.

**No alignment is needed for this fixed-cell periodic structure. Skip
`align.py`.** The lattice supplies the reference frame. Molecular rotation
alignment is unnecessary, and rotating the coordinates alone would change
their relation to the fixed lattice vectors. Use the original, unaligned
trajectory in Step 2 and Step 3.

In this case's `run_template.py`, set:

```python
params['trajectory_xyz_filename'] = path + '/../Rutile_TiO2_MD-pos-1.xyz'
```

### Mg@C60: remove molecular translation and rotation

Use the MD position trajectory from [Step 1](../1_step1/2_Mg_C60/).
If using the precomputed Step 1 results, unpack them first so that the
trajectory exists. Open `align.py` and replace its hard-coded input path
with your trajectory location. For a trajectory in the Step 1 example folder:

```python
from libra_py.md_align import align_trajectory

align_trajectory(
    "../../1_step1/2_Mg_C60/C60_Mg-pos-1.xyz",
    "C60_Mg_MD-pos-1-aligned.xyz",
)
```

Run:

```bash
python align.py
```

The output is `C60_Mg_MD-pos-1-aligned.xyz` in the current folder. Alignment
removes translation and rotation of the molecule as a whole while retaining
its internal motion and preserving atom and frame order. Check that each
frame contains 60 C atoms and one Mg atom.

In this case's `run_template.py`, select the aligned trajectory:

```python
params['trajectory_xyz_filename'] = path + '/../C60_Mg_MD-pos-1-aligned.xyz'
```

### How the trajectory paths are resolved

In both templates, `path = os.getcwd()` is evaluated inside a generated
directory such as `job1/`. Thus, `path + '/../filename.xyz'` locates a
trajectory in the example folder. An absolute path is also valid. Update
any other installation-specific paths in the scripts accordingly.

## 3. Configure `run_template.py`

Use the settings for your chosen case:

| Parameter | TiO₂ | Mg@C60 |
| --- | --- | --- |
| `params['nprocs']` | `8` | `16` |
| `params['lowest_orbital']` | `1` | `90` |
| `params['highest_orbital']` | `40` | `160` |
| `params['is_periodic']` | `True` | `False` |

For both examples, keep `mpi_executable='mpirun'`, `cp2k_exe='cp2k.psmp'`,
`isxTB=False`, `isUKS=False`, and `cube_visualization=False`, provided these
executables match your NESH module. VMD is not required.

For a fresh run without an external wavefunction, set:

```python
params['restart_file'] = False
```

The supplied templates set this to `True`. Libra versions that interpret
the parameter as a file path otherwise attempt to copy a file named `True`.

**Leave these assignments blank in `run_template.py`:**

```python
params['istep'] =
params['fstep'] =
```

Set the overall frame range in `distribute_jobs.py`. Libra fills in the
appropriate limits in each generated `job*/run.py`. The template cannot
be executed directly while these assignments are blank.

### TiO₂: what the periodic-cell settings mean

Retain the supplied rutile lattice in the Python template:

```python
params['is_periodic'] = True

if params['is_periodic']:
    params['A_cell_vector'] = [4.6532721519, 0.0, 0.0]
    params['B_cell_vector'] = [0.0, 4.6532721519, 0.0]
    params['C_cell_vector'] = [0.0, 0.0, 2.9692029953]
    params['periodicity_type'] = 'XYZ'

    origin = [0, 0, 0]
    params['translational_vectors'] = CP2K_methods.generate_translational_vectors(
        origin, [2, 2, 2], params['periodicity_type']
    )
```

`A_cell_vector`, `B_cell_vector`, and `C_cell_vector` are Cartesian lattice
vectors in **ångströms**. They define a tetragonal cell with
`a=b=4.6532721519 Å` and `c=2.9692029953 Å`. `is_periodic=True` tells Libra
to include neighboring-cell contributions to the overlap integrals;
`periodicity_type='XYZ'` includes images in all three directions.

A neighboring cell has displacement

$$
\mathbf{T}=n_a\mathbf{A}+n_b\mathbf{B}+n_c\mathbf{C}.
$$

`origin=[0,0,0]` selects the reference cell. `[2,2,2]` generates integer
indices from −2 to +2 along each direction: 124 neighboring translations,
plus the central cell computed separately. This is the real-space image
sum used for the Γ-point overlap calculation. It controls the image range
for the integrals, not a k-point mesh or the MD simulation-cell size.
Check convergence with image range when changing the cell or basis.

**Keep the lattice identical in Step 1, `run_template.py`, and
`es_diag_temp.inp`.** The CP2K `&CELL` vectors must match those above, with
`PERIODIC XYZ`; retain `PERIODIC XYZ` and `PSOLVER PERIODIC` in `&POISSON`.
The ordinary XYZ trajectory does not contain the cell vectors, so Libra
needs them explicitly. These hard-coded vectors assume the fixed cell
used in the supplied NVT trajectory.

### Mg@C60: molecular boundary conditions

Keep `params['is_periodic']=False`; the periodic-cell block in the Python
template is inactive. The supplied CP2K input uses a cell of about
21.11 Å per side, `PERIODIC NONE`, and the Martyna–Tuckerman (`MT`) Poisson
solver for the isolated molecule. Keep the model consistent with Step 1.

### Electronic-structure settings for both cases

| Quantity | TiO₂ | Mg@C60 |
| --- | --- | --- |
| Explicitly treated electrons | `48` | `250` |
| Occupied spatial orbitals | `24` | `125` |
| HOMO / LUMO | `24 / 25` | `125 / 126` |
| `ADDED_MOS` in CP2K | `50` | `60` |
| `TDDFPT/NSTATES` | `10` | `10` |

The old `Here HOMO is 90` comment in the Mg@C60 Python template is incorrect
for the supplied input. Its saved window is `90..160`, but the HOMO is 125.

`ADDED_MOS` controls extra KS orbitals calculated by CP2K, the Python
orbital window selects those saved by Libra, and `NSTATES` is the number
of excited states. Keep the window wide enough to contain the transitions
retained in Step 3. Preserve the Molden and PDOS print sections in the
CP2K template; Libra reads the Molden orbitals to calculate overlaps.
Per-frame coordinates and CP2K input filenames are generated automatically.

## 4. Set up the NESH batch script

Edit `submit_template.slm` inside your chosen example folder. Use matching
processor counts in the Python and SLURM files:

| Resource | TiO₂ | Mg@C60 |
| --- | --- | --- |
| `params['nprocs']` | `8` | `16` |
| `--nodes` | `1` | `1` |
| `--ntasks-per-node` | `8` | `16` |
| `--cpus-per-task` | `1` | `1` |
| `OMP_NUM_THREADS` | `1` | `1` |
| Supplied wall-time request | `01:00:00` | `06:00:00` |
| Supplied memory request | `64000` MB | `64000` MB |

The following is the **TiO₂** batch example. **For Mg@C60, change
`--ntasks-per-node` to `16` and `--time` to `06:00:00`.** Add the Libra setup
commands before `python run.py` in either case:

```bash
#!/bin/bash -l
#SBATCH --reservation=autumnschool
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
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

Adapt the reservation, account, partition, paths, and limits to your NESH
allocation as directed by the instructors. The commented UB-cluster examples
in the supplied files are unrelated to this NESH run. If you change the
parallel layout, update the MPI process count, SLURM allocation, and thread
count together.

## 5. Select a few frames and submit

In either case's `distribute_jobs.py`, use:

```python
run_slurm = True
submit_template = 'submit_template.slm'
run_python_file = 'run_template.py'
istep = 0
fstep = 2
njobs = 1
submission_exe = 'sbatch'
```

The distributor includes the final frame, so this processes frames
**0, 1, and 2**, producing time overlaps for **0 → 1** and **1 → 2**.
It writes `istep=0` and `fstep=3` in `job1/run.py`, whose Python loop excludes
its upper bound. With multiple jobs, neighboring jobs share a boundary
frame so that adjacent-frame overlaps can be computed.

Submit from the chosen example folder:

```bash
python distribute_jobs.py
squeue -u "$USER"
```

The distributor creates `job1/`, writes `run.py`, copies the batch template,
and calls `sbatch`. Each job launches CP2K for its assigned geometries and
computes the overlaps with Libra.

**Keep `run_slurm=True` on NESH.** In Libra's distributor, `False` executes
the submission script locally with `sh`; it is not a preparation-only mode.
Rerunning the distributor recreates matching `job*` directories, so preserve
completed results before reusing them.

For longer calculations, the reference ranges and job splits are:

| Case | `istep` | `fstep` | `njobs` |
| --- | --- | --- | --- |
| TiO₂ reference | `0` | `20` | `4` |
| Mg@C60 reference | `0` | `300` | `50` |

Use the three-frame settings above during the session. The current TiO₂
distributor defaults to `fstep=10`, `njobs=2`; the Mg@C60 defaults to
`fstep=300`, `njobs=50`, so change these before the short run.

## 6. Inspect the short-run outputs

After the job finishes, inspect `job1/slurm-*.out` and the CP2K logs for
successful SCF and TD-DFPT convergence.

| Output within the example folder | Meaning |
| --- | --- |
| `job1/`, `job2/`, … | Generated scripts, coordinates, CP2K inputs, SLURM output, and working files. |
| `all_logfiles/step_i.log` | CP2K output for frame `i`, including excitation energies and amplitudes needed by Step 3. |
| `all_pdosfiles/` | Projected density-of-states files. |
| `res/E_ks_i.npz` | KS orbital-energy matrices. |
| `res/S_ks_i.npz` | Same-frame KS overlaps at frame `i`. |
| `res/St_ks_i.npz` | KS time overlaps between frames `i` and `i+1`; inputs to the later coupling calculation. |

For the short run, expect three logs and three each of `E_ks` and `S_ks`,
indexed 0–2, plus two `St_ks` files, indexed 0–1. Molden files are removed
after use when `remove_molden=True`. An empty `all_images/` directory may
also be created even with visualization disabled.

## 7. Use the reference data and continue to Step 3

### TiO₂ reference

The results are already in `1_Rutile_TiO2/reference/`. No extraction is
required. Keep your short-run results in the separate directories created
outside `reference/`.

### Mg@C60 reference

From `2_Mg_C60/`, extract the archive into a separate folder:

```bash
mkdir -p reference
tar xvf reference.tar.bz2 -C reference
```

The archive contains `res/`, `all_logfiles/`, and `all_pdosfiles/` at its
root; this command places them under `reference/` and keeps them separate
from the short run.

### Reference contents and checks

| Data within each case's `reference/` | TiO₂ | Mg@C60 |
| --- | --- | --- |
| `all_logfiles/step_i.log` | 21 logs, frames `0..20` | 301 logs, frames `0..300` |
| `res/E_ks_i.npz` and `res/S_ks_i.npz` | 21 of each, indices `0..20` | 301 of each, indices `0..300` |
| `res/St_ks_i.npz` | 20 matrices, indices `0..19` | 300 matrices, indices `0..299` |
| Matrix dimensions | `80 × 80` | `142 × 142` |

From either example folder, inspect the reference matrices with:

```bash
python - <<'PY'
from scipy.sparse import load_npz

S = load_npz("reference/res/S_ks_0.npz").toarray()
St = load_npz("reference/res/St_ks_0.npz").toarray()
print("Same-frame overlap:", S.shape)
print("Time overlap:", St.shape)
PY
```

Same-frame overlaps should be close to identity. Time overlaps need not
be identity or symmetric, and orbital signs can differ between frames.

### Step 3 settings

Use the [Step 3 guide](../3_step3/README.md) with matching trajectory,
logs, and overlap matrices. The following reference-data paths are relative
to a working directory of `3_step3/`:

| Setting | TiO₂ | Mg@C60 |
| --- | --- | --- |
| Reference-data root | `../2_step2/1_Rutile_TiO2/reference` | `../2_step2/2_Mg_C60/reference` |
| `lowest_orbital` | `1` | `90` |
| `nstates` | `11` | `11` |
| Example output directory | `ci_tio2` | `ci_mg_c60` |

Set the log paths to `all_logfiles/step_i.log` and matrix paths to `res/`
under the selected data root. `nstates=11` includes the ground state plus
the 10 calculated excited states. Both supplied Step 1 inputs use a 1 fs
MD step and save every step, giving `dt=1.0 * units.fs2au`; use the actual
saved-frame spacing if you change the MD settings.

Use the unaligned trajectory for TiO₂ and the aligned trajectory for your
aligned Mg@C60 run. Reference electronic results must be used with the
corresponding trajectory and coordinate convention. Keep the same orbital
window and frame ordering throughout.

For interval `i → i+1`, pair `St_ks_i.npz` with electronic states from
`step_i.log` and `step_{i+1}.log`. Follow the
[time-indexing instructions](../3_step3/README.md#time-indexing-align-both-endpoints)
before running the supplied Step 3 loop.
