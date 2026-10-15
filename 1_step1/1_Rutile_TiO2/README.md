# Step 1: Run molecular dynamics using CP2K

Generate a ground-state molecular dynamics (MD) trajectory of rutile TiO₂.
The nuclear geometries saved here become the input for the electronic-structure
calculations in [Step 2](../2_step2/README.md).

The supplied example is a short workflow demonstration, not an equilibrated
production trajectory. It contains two Ti and four O atoms in a periodic unit
cell of approximately 4.6533 × 4.6533 × 2.9692 Å.

## Files and settings

| File or folder | Purpose |
| --- | --- |
| `TiO2_unit_cell.vasp` | Starting coordinates and lattice vectors; an example structure is already provided. |
| `convert.py` | Uses Libra to convert the VASP structure to XYZ and print the cell information. |
| `md.inp` | CP2K electronic-structure and MD settings. |
| `submit.slm` | Example SLURM submission script; adapt it to your cluster. |
| `reference/` | Completed example, including trajectory, energies, log, and restart files. |

The main choices in `md.inp` are:

| Setting | Supplied value |
| --- | --- |
| Electronic structure | Restricted PBE with DFT-D3 dispersion; GPW method |
| Basis and pseudopotentials | Ti: DZVP-MOLOPT-SR-GTH / GTH-PBE-q12; O: DZVP-MOLOPT-GTH / GTH-PBE-q6 |
| Grid | `CUTOFF 500`, `REL_CUTOFF 50` |
| SCF | `EPS_SCF 1.0E-6`, up to 50 iterations |
| Boundary conditions | Periodic in X, Y, and Z |
| MD | NVT, target temperature 300 K, CSVR thermostat with `TIMECON 20` |
| Duration and output | 20 steps of 1 fs; positions saved every step |

## How to run

1. From the repository root, enter this folder and activate your Libra Python
   environment (for example, `conda activate libra`, if available):

   ```bash
   cd 1_step1
   python convert.py
   ```

   This reads `TiO2_unit_cell.vasp` and writes `TiO2_unit_cell.xyz`. Check that
   the XYZ file contains six atoms and that the printed lattice vectors match
   `&CELL` in `md.inp`. XYZ coordinates alone do not specify the periodic cell.
   If replacing the structure with a VESTA export, keep the filename consistent.

2. Review `md.inp`. Keep the geometry, cell, element names, basis sets, and
   pseudopotentials consistent. CP2K must be able to locate `BASIS_MOLOPT`,
   `POTENTIAL`, and `dftd3.dat` through your installation or configured paths.

3. Edit `submit.slm` for your account, partition, modules, and MPI launcher.
   The supplied paths are specific to the example cluster. Match the CPU
   allocation to the launch settings: the script currently requests one CPU
   per task but sets `OMP_NUM_THREADS=4`; either allocate four CPUs per task
   or reduce the thread count. It launches six MPI ranks.

4. Submit and inspect the CP2K log and SLURM output:

   ```bash
   sbatch submit.slm
   ```

   Check for successful SCF convergence and completion of the requested MD
   steps before using the trajectory. If CP2K is unavailable, use `reference/`
   for the analysis exercises and the supplied trajectory in Step 2.

## Results and handoff

| Output | What to inspect |
| --- | --- |
| `Rutile_TiO2_MD-pos-1.xyz` | Position trajectory; pass this to Step 2. |
| `Rutile_TiO2_MD-1.ener` | Time, kinetic energy, temperature, potential energy, conserved quantity, and timing; units are in the header. |
| `Rutile_TiO2_MD-vel-1.xyz`, `Rutile_TiO2_MD-frc-1.xyz` | Saved velocities and forces. |
| `out_md_tio2_dft.log`, `slurm-*.out` | Calculation details, convergence, and execution messages. |
| `Rutile_TiO2_MD-1.restart`, `Rutile_TiO2_MD-RESTART.wfn` | MD restart input and electronic wavefunction; `.bak-1` files are backups. |

The reference trajectory contains frames 0–20: 21 geometries spanning 20 fs,
with a saved-frame spacing of 1 fs. Record this spacing for Steps 3 and 4;
if you change the MD timestep or trajectory output frequency, update the
later analysis accordingly.

From `1_step1`, copy a successfully completed trajectory with:

```bash
cp Rutile_TiO2_MD-pos-1.xyz ../2_step2/
```

## Student exercises

Use the reference outputs if you are not running CP2K. Submit a short notebook
or report with the requested plots and a few sentences per exercise.

1. **Check the structure.** Read the VASP and XYZ files, identify the atom
   counts and cell lengths, and calculate the unit-cell volume in Å³. Explain
   which information must be supplied separately when using an XYZ file.
2. **Inspect the trajectory.** Count frames and read their times from the XYZ
   comments. Verify the 1 fs spacing and 20 fs duration. How would saving every
   fifth MD step change the number of frames and the spacing passed downstream?
3. **Assess the MD.** Plot temperature, potential energy, and the conserved
   quantity against time using `reference/Rutile_TiO2_MD-1.ener`. Compare the
   temperature with the 300 K target. Does this short trajectory support a
   claim of equilibration? For NVT, distinguish the printed conserved quantity
   from kinetic plus potential energy.
4. **Optional short rerun.** In a separate working folder, change only the
   thermostat time constant or the MD timestep. Compare the resulting traces,
   keeping the physical duration comparable. Record the changed setting and
   explain what a longer run would be needed to establish.
