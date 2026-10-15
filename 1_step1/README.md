# Step 1: Molecular dynamics with CP2K

This tutorial generates ground-state molecular-dynamics (MD) trajectories
for **periodic rutile TiO₂** and **molecular Mg@C60**. At each MD step, CP2K
converges the electronic ground state and calculates the forces used to
move the nuclei. The saved geometries are the inputs for the electronic
calculations and orbital overlaps in [Step 2](../2_step2/).

During the session, run a few MD steps to learn the workflow, then use the
supplied reference trajectories for the later steps. These short examples
demonstrate the workflow; they do not establish thermal equilibration or
production-quality sampling.

## Choose your example

| Setting | Rutile TiO₂ | Mg@C60 |
| --- | --- | --- |
| Folder within `1_step1/` | [1_Rutile_TiO2/](1_Rutile_TiO2/) | [2_Mg_C60/](2_Mg_C60/) |
| Initial XYZ structure | `TiO2_unit_cell.xyz` | `C60_Mg.xyz` |
| Atoms | Two Ti and four O | 60 C and one Mg |
| Boundary conditions | Periodic in X, Y, and Z | Isolated molecule |
| CP2K project name | `Rutile_TiO2_MD` | `C60_Mg` |
| Supplied MD length | 20 steps of 1 fs | 500 steps of 1 fs |
| Target temperature / ensemble | 300 K / NVT | 300 K / NVT |
| MPI processes in `submit.slm` | `8` | `16` |
| CP2K output log | `out_md_tio2_dft.log` | `out.log` |
| Reference results | Already in `reference/` | Unpack `reference.tar.bz2` |

The requested MD length and the number of frames actually present in a
reference file are different quantities; see Section 7 for the verified
reference contents.

## Files in the example folders

| File or folder | Purpose |
| --- | --- |
| `md.inp` in each folder | CP2K electronic-structure, cell, MD, thermostat, and output settings. Review this before running. |
| `submit.slm` in each folder | SLURM resources, NESH environment, and the MPI command that runs CP2K. |
| [1_Rutile_TiO2/TiO2_unit_cell.vasp](1_Rutile_TiO2/TiO2_unit_cell.vasp) | Starting rutile structure with Cartesian coordinates and lattice vectors. |
| [1_Rutile_TiO2/TiO2_unit_cell.xyz](1_Rutile_TiO2/TiO2_unit_cell.xyz) | Six-atom structure read by the TiO₂ CP2K input. Already provided. |
| [1_Rutile_TiO2/convert.py](1_Rutile_TiO2/convert.py) | Uses Libra to convert the VASP structure to XYZ. Needed only if you regenerate the XYZ from the VASP file. |
| [1_Rutile_TiO2/reference/](1_Rutile_TiO2/reference/) | Completed 20-step TiO₂ example, including positions, energies, velocities, forces, log, and restart files. |
| [2_Mg_C60/C60_Mg.xyz](2_Mg_C60/C60_Mg.xyz) | Starting Mg@C60 geometry with 61 atoms. |
| [2_Mg_C60/reference.tar.bz2](2_Mg_C60/reference.tar.bz2) | Mg@C60 positions, energies, velocities, forces, and restart files. |
| [2_Mg_C60/out.log](2_Mg_C60/out.log) | Supplied example execution log. It includes SCF nonconvergence aborts; inspect it together with the available trajectory. |

The two `md.inp` files are available here: [TiO₂](1_Rutile_TiO2/md.inp)
and [Mg@C60](2_Mg_C60/md.inp). Their submission scripts are
[TiO₂](1_Rutile_TiO2/submit.slm) and [Mg@C60](2_Mg_C60/submit.slm).

## 1. Enter the folder for your case

From the repository root, choose **one** working folder.

For TiO₂:

```bash
cd 1_step1/1_Rutile_TiO2
```

For Mg@C60:

```bash
cd 1_step1/2_Mg_C60
```

**Run the remaining shell commands from your chosen example folder** unless
another directory is specified. The links in this README are relative to
`1_step1/`; the commands assume you have entered the appropriate subfolder.

## 2. Check the starting structure

### TiO₂: periodic rutile cell

The supplied `TiO2_unit_cell.xyz` is ready to use. If you edit or replace
`TiO2_unit_cell.vasp`, regenerate the XYZ with `convert.py`. Load Libra
on NESH for this Python conversion:

```bash
export PYTHONPATH=/gxfs_home/cau/suphc424/Software/libra-code/_build/src:$PYTHONPATH
source /gxfs_home/cau/suphc424/.conda/envs/libra/bin/activate
python convert.py
```

The script reads `TiO2_unit_cell.vasp` and writes `TiO2_unit_cell.xyz`.
This simple converter assumes Cartesian coordinates and a scale factor
of 1.0, as in the supplied VASP file. Check that the XYZ contains two Ti
and four O atoms. In `md.inp`, retain:

```text
&CELL
  A  4.6532721519  0.0           0.0
  B  0.0           4.6532721519  0.0
  C  0.0           0.0           2.9692029953
  PERIODIC XYZ
&END CELL
```

These Cartesian lattice vectors are in ångströms. They define a tetragonal
cell with `a=b=4.6532721519 Å` and `c=2.9692029953 Å`. The XYZ file contains
coordinates but does not specify this lattice, so keep `&CELL` consistent
with the source structure. The `&POISSON` section also uses `PERIODIC XYZ`
and `PSOLVER PERIODIC`.

### Mg@C60: isolated molecule

Use `C60_Mg.xyz`, containing 60 carbon atoms and one magnesium atom.
The CP2K input places the molecule in a cell of about 21.11 Å per side
and activates `CENTER_COORDINATES`.

Keep `PERIODIC NONE` in `&CELL` and `&POISSON`, with `PSOLVER MT` in
`&POISSON`. The Martyna–Tuckerman solver handles the isolated-system
electrostatics. The finite computational box does not make this a periodic
crystal.

For either case, check `COORD_FILE_NAME` in `&TOPOLOGY` and preserve the atom
order throughout the workflow. CP2K and the Python converter use paths
relative to their execution directory; use absolute paths if needed.

## 3. Review the CP2K input

Both inputs use restricted PBE with DFT-D3 dispersion and the GPW method.
The active settings are:

| Setting in `md.inp` | TiO₂ | Mg@C60 |
| --- | --- | --- |
| `RUN_TYPE` | `MD` | `MD` |
| Basis / pseudopotential | Ti: `DZVP-MOLOPT-SR-GTH` / `GTH-PBE-q12`; O: `DZVP-MOLOPT-GTH` / `GTH-PBE-q6` | C: `DZVP-MOLOPT-GTH` / `GTH-PBE-q4`; Mg: `DZVP-MOLOPT-SR-GTH` / `GTH-PBE-q10` |
| `CUTOFF` / `REL_CUTOFF` | `500 / 50` | `300 / 50` |
| SCF method | Diagonalization with Broyden mixing | Diagonalization with Broyden mixing |
| `EPS_SCF` / `MAX_SCF` | `1.0E-6 / 50` | `1.0E-6 / 500` |
| `ADDED_MOS` | `50` | `60` |
| `ENSEMBLE` / `TEMPERATURE` | `NVT / 300 K` | `NVT / 300 K` |
| `TIMESTEP` | `1.0 fs` | `1.0 fs` |
| Thermostat | `CSVR`, `TIMECON 20` | `CSVR`, `TIMECON 20` |
| Supplied `STEPS` | `20` | `500` |
| Trajectory print interval | `MD 1` | `MD 1` |

`NVT` keeps the particle number and volume fixed while the thermostat
controls the temperature. The target of 300 K does not mean that every
instantaneous temperature is exactly 300 K. The CSVR time constant is
20 fs in these inputs.

`MAX_SCF` limits the electronic iterations at each geometry; `STEPS` inside
`&MOTION / &MD` controls the number of nuclear MD steps. `ADDED_MOS` adds
unoccupied KS orbitals; it does not request excited-state calculations.
The TD-DFPT calculations are performed later in Step 2.

For the short tutorial run, change the existing `STEPS` line inside
`&MOTION / &MD` to:

```text
STEPS 2
```

With the existing 1 fs timestep and every-step trajectory output, a fresh
two-step run should save three geometries at 0, 1, and 2 fs.

For a fresh demonstration without an existing wavefunction, use
`SCF_GUESS ATOMIC`. The supplied TiO₂ input already uses it; the Mg@C60 input
currently specifies `SCF_GUESS RESTART`. This setting controls the initial
electronic guess; an MD restart also requires the nuclear and thermostat
state from the MD restart file.

CP2K must be able to locate `BASIS_MOLOPT`, `POTENTIAL`, and `dftd3.dat`.
Use the NESH module's data-file configuration or adapt the corresponding
input paths for your installation. Keep the structure, element names,
basis sets, pseudopotentials, and boundary conditions consistent.

## 4. Configure the NESH submission script

Step 1 runs CP2K directly through `submit.slm`. The Libra Python environment
is used for the optional structure conversion above and the later Step 2
workflow.

Review the resources and launch command for your case:

| Setting in `submit.slm` | TiO₂ | Mg@C60 |
| --- | --- | --- |
| `--nodes` | `1` | `1` |
| `--ntasks-per-node` | `8` | `16` |
| `--cpus-per-task` | `1` | `1` |
| `OMP_NUM_THREADS` | `1` | `1` |
| MPI launch count | `mpirun -np 8` | `mpirun -np 16` |
| Supplied wall time | `02:00:00` | `06:00:00` |
| Supplied memory | `24000` MB | `64000` MB |
| CP2K log filename | `out_md_tio2_dft.log` | `out.log` |

A TiO₂ submission script using the autumn-school module setup is:

```bash
#!/bin/bash -l
#SBATCH --reservation=autumnschool
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=24000
#SBATCH --job-name=tio2_md

module load gcc12-env
module load cp2k-autumnschool/default
module load cp2k

export OMP_NUM_THREADS=1
mpirun -np 8 cp2k.psmp -i md.inp -o out_md_tio2_dft.log
```

For **Mg@C60**, use `--ntasks-per-node=16`, `--time=06:00:00`,
`--mem=64000`, a suitable job name, and replace the final command with:

```bash
mpirun -np 16 cp2k.psmp -i md.inp -o out.log
```

Keep the `autumnschool` reservation only when it is available for your
session. Adapt account, partition, modules, and executable paths to your
NESH allocation as directed by the instructors. The Mg@C60 script already
contains the autumn-school module commands; its commented UB-cluster
examples are unrelated to this NESH run.

**Match the MPI launch count to the SLURM allocation.** These scripts use
one CPU per MPI task and `OMP_NUM_THREADS=1`. Update the allocation,
process count, and thread count together if you change the parallel layout.

## 5. Submit and check the calculation

From the chosen example folder:

```bash
sbatch submit.slm
squeue -u "$USER"
```

After the job finishes, inspect `slurm-*.out` and your case's CP2K log.
Check that SCF converges for the saved geometries and that CP2K reaches the
requested final MD step. A job disappearing from the queue or producing
an XYZ file does not establish successful completion.

Inspect the `.ener` file as well. Its header identifies time in fs,
temperature in K, and energies in atomic units. Plot temperature, potential
energy, and the printed conserved quantity against time. For NVT, the
printed conserved quantity includes the thermostat bookkeeping and should
not be identified with kinetic plus potential energy alone.

If an SCF failure occurs, retain the log and identify the failed geometry
before using the trajectory downstream. Check convergence and the actual
available frame range rather than assuming the requested MD length was
completed.

## 6. Understand the output files

The CP2K `PROJECT_NAME` controls most output filenames. In the patterns
below, `<PROJECT>` means `Rutile_TiO2_MD` for TiO₂ or `C60_Mg` for Mg@C60.

| Output | Contents and use |
| --- | --- |
| `<PROJECT>-pos-1.xyz` | Multi-frame position trajectory. This is the main input to Step 2. |
| `<PROJECT>-1.ener` | MD step, time, kinetic energy, temperature, potential energy, conserved quantity, and timing. |
| `<PROJECT>-vel-1.xyz` | Saved velocities. |
| `<PROJECT>-frc-1.xyz` | Saved forces. |
| `<PROJECT>-1.restart` | CP2K restart input with the MD state. |
| `<PROJECT>-RESTART.wfn` | Electronic wavefunction restart file. |
| `*.bak-1` | Backup restart files. |
| CP2K log and `slurm-*.out` | Convergence, MD progress, and execution information. |

Preserve outputs you need before reusing the same project name and working
directory. Keep the positions, energies, log, and restart files from each
run together.

The time between saved trajectory frames is

$$
\Delta t_{\mathrm{saved}}
=\mathrm{TIMESTEP}\times\text{trajectory output interval}.
$$

Here, `TIMESTEP 1.0 fs` and `&TRAJECTORY / &EACH / MD 1` give a saved-frame
spacing of **1 fs**. Saving every fifth MD step would give 5 fs spacing.
Record this value for the later overlap and coupling calculations.

## 7. Use the supplied reference trajectories

### TiO₂

The files are already unpacked in `1_Rutile_TiO2/reference/`. The trajectory
contains **21 geometries**, indices **0–20**, with times **0–20 fs** and
1 fs spacing. Its log records completion of the 20-step run.

This short reference is not an equilibrated 300 K trajectory: its
instantaneous temperature rises substantially above the target. Use the
energy and temperature traces to distinguish a completed calculation from
an equilibrated trajectory.

### Mg@C60

From `1_step1/2_Mg_C60/`, unpack the reference into a separate folder to
retain your short-run results:

```bash
mkdir -p reference
tar xvf reference.tar.bz2 -C reference
```

The archive stores its files at the top level, so this creates, for example,
`reference/C60_Mg-pos-1.xyz` and `reference/C60_Mg-1.ener`. It also contains
velocities, forces, and restart files; the supplied `out.log` is outside
the archive.

The archived position trajectory and energy file contain **316 frames/rows**,
indices **0–315**, with times **0–315 fs** and 1 fs spacing. Although the
input requests 500 MD steps, the archive does not contain a completed
500-step trajectory. The supplied `out.log` includes SCF nonconvergence
aborts. The available frame range covers the `0..300` range used by the
Mg@C60 Step 2 reference, but it should not be described as a completed
500-step run.

## 8. Pass the trajectory to Step 2

Follow the [combined Step 2 tutorial](../2_step2/README.md). Use the
trajectory corresponding to the run or reference dataset you selected,
and keep the atom order, geometry convention, and saved-frame spacing
consistent throughout the workflow.

### TiO₂ handoff

From `1_step1/1_Rutile_TiO2/`, copy your new trajectory:

```bash
cp Rutile_TiO2_MD-pos-1.xyz ../../2_step2/1_Rutile_TiO2/
```

To use the supplied reference instead, copy:

```bash
cp reference/Rutile_TiO2_MD-pos-1.xyz ../../2_step2/1_Rutile_TiO2/
```

**No molecular alignment is needed for the fixed periodic TiO₂ cell.** Use
the original trajectory and retain the same lattice vectors and `XYZ`
periodicity in the Step 2 CP2K and Libra settings.

### Mg@C60 handoff

From `1_step1/2_Mg_C60/`, copy your new trajectory:

```bash
cp C60_Mg-pos-1.xyz ../../2_step2/2_Mg_C60/
```

To use the extracted reference instead, copy:

```bash
cp reference/C60_Mg-pos-1.xyz ../../2_step2/2_Mg_C60/
```

In `2_step2/2_Mg_C60/`, edit `align.py` to read the copied trajectory:

```python
from libra_py.md_align import align_trajectory

align_trajectory("C60_Mg-pos-1.xyz", "C60_Mg_MD-pos-1-aligned.xyz")
```

Load Libra, run `python align.py`, and select the aligned file in Step 2's
`run_template.py`. This removes whole-molecule translation and rotation
while retaining the internal motion. Set `is_periodic=False` for this
molecular example.

For either case, Step 2 frame numbers index the geometries actually saved
in the XYZ file. A two-step practice run provides frames `0..2`; the longer
reference datasets provide the ranges listed above.
