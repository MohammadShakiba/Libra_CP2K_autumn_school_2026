# Step 2: Electronic-structure calculations and orbital overlaps for rutile TiO2

In this tutorial, we use the periodic rutile TiO₂ MD trajectory from
[Step 1](../../1_step1/1_Rutile_TiO2/README.md) to run CP2K single-point
calculations and calculate Kohn–Sham (KS) orbital overlaps with Libra.
CP2K also produces TD-DFPT excitation energies and amplitudes. These data
are the inputs for the many-electron overlaps and vibronic Hamiltonians
constructed in [Step 3](../../3_step3/README.md).

During the session, we will process **frames 0, 1, and 2**, giving two
adjacent-frame intervals. The supplied `reference/` folder contains results
for frames 0–20, allowing us to continue to Step 3 after the short run.

**No molecular alignment is needed for this fixed-cell periodic structure.**
Use the original MD trajectory and keep its coordinates consistent with
the lattice vectors.

## Files in this folder

| File or folder | Purpose | What you need to do |
| --- | --- | --- |
| [README.md](README.md) | Instructions for the rutile TiO₂ Step 2 workflow. | Follow the steps below. |
| [run_template.py](run_template.py) | Defines CP2K/Libra settings and calls `step2.run_cp2k_libint_step2`. Becomes `run.py` inside each job directory. | Set the trajectory path, processor count, orbital window, and periodic-cell settings. Leave `istep` and `fstep` blank. |
| [distribute_jobs.py](distribute_jobs.py) | Splits the requested frame range among jobs, fills in their frame limits, and launches the submission scripts. | Set `istep`, `fstep`, `njobs`, and the submission options. |
| [submit_template.slm](submit_template.slm) | SLURM resource requests, NESH modules, and the command that runs each generated `run.py`. | Check the allocation and environment setup; match the CPU resources to `nprocs`. |
| [es_diag_temp.inp](es_diag_temp.inp) | CP2K DFT/TD-DFPT template, including periodic boundaries, basis sets, pseudopotentials, and Molden/PDOS output. | Keep the cell consistent with Step 1 and the Python template. Preserve the output needed by Libra. |
| [align.py](align.py) | Molecular translation/rotation alignment utility retained in the folder. | **Skip this script for the periodic TiO₂ example.** |
| [reference/](reference/) | Completed calculations for spatial orbitals 1–40, including logs, overlap matrices, PDOS, and generated jobs. | Use these results for inspection and the Step 3 handoff. No extraction is required. |
| [clean.sh](clean.sh) | Removes `job*`, `all_*`, and `res` from the current directory. | Save any results you need before using it. |

The trajectory is provided in the Step 1 `reference/` folder. This Step 2
folder does not contain the `input/` or `smaller_orbital_space/` directories
mentioned in older versions of the instructions.

## 1. Load Libra on NESH

From the repository root, enter the rutile Step 2 folder:

```bash
cd 2_step2/1_Rutile_TiO2
```

Load the tutorial Libra installation:

```bash
export PYTHONPATH=/gxfs_home/cau/suphc424/Software/libra-code/_build/src:$PYTHONPATH
source /gxfs_home/cau/suphc424/.conda/envs/libra/bin/activate
```

Check the Python environment:

```bash
python -c "import liblibra_core; from libra_py.workflows.nbra import step2; print('Libra is ready')"
```

Adapt these paths if you use another installation, and use the same
environment in the batch script. Run the remaining commands from
`2_step2/1_Rutile_TiO2` unless another directory is specified.

## 2. Use the original periodic trajectory

Copy the supplied Step 1 trajectory into the current folder:

```bash
cp ../../1_step1/1_Rutile_TiO2/reference/Rutile_TiO2_MD-pos-1.xyz .
```

If you ran Step 1 yourself, use the trajectory from your own calculation
and retain its atom order, frame order, and cell definition. The supplied
rutile cell contains six atoms: two Ti and four O.

**Do not run `align.py` for this tutorial.** The periodic lattice supplies
the reference frame, so molecular translation/rotation alignment is
unnecessary. Rotating the coordinates alone would change their relation
to the fixed lattice vectors. Keep the original, unaligned trajectory
for both Step 2 and Step 3.

In `run_template.py`, replace the trajectory setting with:

```python
params['trajectory_xyz_filename'] = path + '/../Rutile_TiO2_MD-pos-1.xyz'
```

`path = os.getcwd()` is evaluated inside a generated directory such as
`job1/`, so `../` points to the Step 2 folder containing the copied XYZ.
An absolute trajectory path is also valid. This replacement corrects the
supplied template's relative path to Step 1, which is too shallow when
resolved from inside `job*/`.

## 3. Understand the periodic-cell settings

Keep the following settings in `run_template.py` for the supplied rutile cell:

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

| Setting | Meaning |
| --- | --- |
| `is_periodic=True` | Tells Libra to include neighboring-cell contributions when constructing the orbital overlaps. |
| `A_cell_vector`, `B_cell_vector`, `C_cell_vector` | The three lattice vectors in **ångströms**, expressed in Cartesian coordinates. Here they define a tetragonal cell with `a=b=4.6532721519 Å` and `c=2.9692029953 Å`. |
| `periodicity_type='XYZ'` | Includes periodic images in all three directions, as appropriate for bulk rutile. |
| `origin=[0,0,0]` | Selects the reference cell in the integer lattice-index representation. |
| `[2,2,2]` | Includes image indices from −2 to +2 along each periodic direction in the overlap calculation. |

A neighboring cell is displaced by

$$
\mathbf{T}=n_a\mathbf{A}+n_b\mathbf{B}+n_c\mathbf{C}.
$$

`generate_translational_vectors` returns the integer triples
`(n_a, n_b, n_c)`; Libra combines them with the lattice vectors to obtain
the physical displacements. With `XYZ` and `[2,2,2]`, there are
`5 × 5 × 5 − 1 = 124` neighboring translations. Libra computes the central
cell separately, giving 125 cell contributions in total.

This real-space image sum accounts for overlap between basis functions
in different cells in the Γ-point calculation used here. `[2,2,2]` controls
the image range for those integrals; it does not specify a k-point mesh
or enlarge the MD simulation cell. Its convergence can be checked by
increasing the image range, particularly when changing the cell or basis.

**Use the same lattice in Step 1, `run_template.py`, and `es_diag_temp.inp`.**
In the CP2K input, the `&CELL` vectors must match the values above, with
`PERIODIC XYZ`; the `&POISSON` section must also retain `PERIODIC XYZ` and
`PSOLVER PERIODIC`. The ordinary XYZ trajectory does not carry these cell
vectors, which is why Libra needs them explicitly. The supplied settings
assume the fixed cell used in the NVT trajectory.

## 4. Review the remaining calculation settings

| Setting in `run_template.py` | Value for this tutorial |
| --- | --- |
| `params['nprocs']` | `8`, matching the NESH submission template. |
| `params['mpi_executable']` | `'mpirun'`, compatible with the loaded CP2K module. |
| `params['cp2k_exe']` | `'cp2k.psmp'`, or the appropriate executable path for your installation. |
| `params['lowest_orbital']` | `1`. |
| `params['highest_orbital']` | `40`, inclusive. |
| `params['isxTB']` | `False`, selecting the DFT route. |
| `params['isUKS']` | `False`, for the restricted calculation. |
| `params['cube_visualization']` | `False`; VMD is not required. |
| `params['restart_file']` | `False` for a fresh run without an external wavefunction. |

The template currently sets `restart_file=True`. For the fresh run, change
it to `False`: Libra versions that interpret this parameter as a file path
otherwise attempt to copy a file named `True`.

Leave the frame assignments blank in `run_template.py`:

```python
params['istep'] =
params['fstep'] =
```

The distributor fills them in when creating `job*/run.py`. Do not execute
`run_template.py` directly while these assignments are blank.

The supplied CP2K template uses PBE-D3 with Ti and O GTH pseudopotentials.
The reference log reports **48 explicitly treated electrons and 24 occupied
spatial orbitals**, so the HOMO is orbital **24** and the LUMO is **25**.
Saving spatial orbitals `1..40` gives **80 × 80** matrices containing both
spin blocks.

`ADDED_MOS 50` sets the extra KS orbitals calculated by CP2K, the Python
orbital window selects those saved by Libra, and `TDDFPT/NSTATES 10`
requests 10 excited states. Keep these choices consistent with the states
and transitions used in Step 3. Preserve the Molden and PDOS print sections;
Libra reads the Molden orbitals to calculate overlaps. Per-frame coordinate
filenames in the CP2K input are generated automatically from the trajectory.

The Python template also contains a `cp2k_ot_input_template` setting,
although `es_ot_temp.inp` is not supplied in this folder. That file is not
read by the DFT route used here (`isxTB=False`).

## 5. Check the NESH submission template

Review `submit_template.slm`. For the supplied eight-process setup, use:

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

Add the Libra setup commands before `python run.py` in the supplied file.
Adapt the reservation, account, partition, installation paths, and resource
limits to your NESH allocation as directed by the instructors. The commented
UB-cluster examples are unrelated to this NESH run.

**Match the resources:** one node with eight tasks gives eight MPI processes,
consistent with `params['nprocs']=8`. With `--cpus-per-task=1`, retain
`OMP_NUM_THREADS=1`. Update the process count, SLURM allocation, and thread
count together if you change the parallel layout.

## 6. Run a few frames

In `distribute_jobs.py`, select the short tutorial calculation:

```python
run_slurm = True
submit_template = 'submit_template.slm'
run_python_file = 'run_template.py'
istep = 0
fstep = 2
njobs = 1
submission_exe = 'sbatch'
```

The first XYZ frame has index zero. The distributor includes the final
frame, so this processes frames **0, 1, and 2**, producing time overlaps for
**0 → 1** and **1 → 2**. It writes `istep=0` and `fstep=3` in `job1/run.py`,
whose Python loop excludes its upper bound.

Submit from the rutile Step 2 folder:

```bash
python distribute_jobs.py
squeue -u "$USER"
```

The distributor creates `job1/`, generates `run.py`, copies the submission
template, and calls `sbatch`. Each job runs CP2K for its assigned geometries
and computes the overlaps with Libra.

Keep `run_slurm=True` on NESH. In Libra's distributor, `False` executes the
submission script locally with `sh`; it is not a preparation-only mode.
Rerunning the distributor recreates matching `job*` directories, so preserve
completed results before reusing them.

For the full supplied trajectory, use `istep=0`, `fstep=20`, and `njobs=4`.
The jobs cover frames 0–5, 5–10, 10–15, and 15–20. Shared boundary frames
allow adjacent-frame overlaps across the full range. Their generated
exclusive `fstep` values are 6, 11, 16, and 21. The current distributor
defaults (`fstep=10`, `njobs=2`) cover only part of that trajectory.

## 7. Inspect the outputs and continue to Step 3

After the job finishes, inspect `job1/slurm-*.out` and the collected CP2K
logs for successful SCF and TD-DFPT convergence.

| Output | Meaning |
| --- | --- |
| `job1/`, `job2/`, … | Generated scripts, coordinates, CP2K inputs, SLURM output, and working files. |
| `all_logfiles/step_i.log` | CP2K output for frame `i`, including excitation energies and amplitudes. |
| `all_pdosfiles/` | Projected density-of-states files. |
| `res/E_ks_i.npz` | KS orbital-energy matrices. |
| `res/S_ks_i.npz` | Same-frame overlaps at frame `i`. |
| `res/St_ks_i.npz` | Time overlaps between frames `i` and `i+1`; inputs for the later coupling calculation. |

The short run should produce three logs, three `E_ks` matrices, and three
`S_ks` matrices, indexed 0–2, plus two `St_ks` matrices, indexed 0–1.
The precomputed `reference/` contains 21 logs and 21 each of `E_ks` and
`S_ks` for frames 0–20, plus 20 `St_ks` matrices for intervals 0–19.
These results are already unpacked and can be used directly.

For a quick inspection of the reference matrices:

```bash
python - <<'PY'
from scipy.sparse import load_npz

S = load_npz("reference/res/S_ks_0.npz").toarray()
St = load_npz("reference/res/St_ks_0.npz").toarray()
print("Same-frame overlap:", S.shape)
print("Time overlap:", St.shape)
PY
```

Both should be `(80, 80)` for the `1..40` orbital window. Same-frame
overlaps should be close to identity; time overlaps need not be identity
or symmetric.

For [Step 3](../../3_step3/README.md), use the matching **unaligned rutile
trajectory**, `reference/res/`, and `reference/all_logfiles/`, with paths
resolved from the Step 3 working directory. Set `lowest_orbital=1` and
`nstates=11` for the ground state plus 10 excited states. The supplied MD
input uses a 1 fs step and writes every step, giving
`dt=1.0 * units.fs2au`.

Pair `St_ks_i.npz` with electronic states from `step_i.log` and
`step_{i+1}.log`. Follow the
[time-indexing instructions](../../3_step3/README.md#time-indexing-align-both-endpoints)
before running the supplied Step 3 loop, and keep logs, orbital windows,
and trajectory frames consistent throughout the workflow.
