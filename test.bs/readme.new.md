# MiSTI
PSMC-based Migration and Split Time Inference (MiSTI) from two genomes

## Overview
MiSTI infers split time and migration parameters using two PSMC trajectories and the joint site frequency spectrum (joint SFS) of two genomes.

High-level workflow

```text
Genome 1  ->  PSMC  ->  genome1.psmc
Genome 2  ->  PSMC  ->  genome2.psmc
Genome 1 + Genome 2  ->  joint SFS  ->  mi.sfs
genome1.psmc + genome2.psmc + mi.sfs + split time  ->  MiSTI  ->  results.mi
results.mi + psmc files  ->  MiSTIPlot  ->  plot.pdf
```

## Data preparation

### 1. Run PSMC on two genomes
Run PSMC (Li and Durbin 2011) separately on two genomes to obtain two PSMC output files.

### 2. Generate the joint site frequency spectrum (joint SFS)
Compute the joint SFS for the same two genomes.

If you use RealSFS (ANGSD), convert the RealSFS output into the MiSTI format using `utils/ANGSDSFS.py`:

```bash
./utils/ANGSDSFS.py real.sfs [genome1_label genome2_label] > mi.sfs
```

Here, `real.sfs` is the file produced by RealSFS. It is recommended to provide labels (genome or population names) to record the genome order in `real.sfs`.

Example

```bash
./utils/ANGSDSFS.py Yoruba_French.real.sfs YRI FRE > Yoruba_French.mi.sfs
```

For simulated data, compute the joint SFS from simulated chromosomes using:

Table: scripts for simulated data

| Simulator output format | Script |
|---|---|
| msHOT-lite with `-l` format | `utils/MS2JSFS.py` |
| standard ms or scrm output format | `utils/SCRM2JAF.py` |

Important constraint  
If one of the genomes is ancient, it must be the second genome in the SFS.

## About the time scale
The time axis in the MiSTI command line is defined by the PSMC time discretisation. Because discretisations typically differ between genomes, MiSTI merges the time points from the two PSMC files. Consequently, within each merged time interval, the coalescence rates for both genomes are constant.

All times provided in the MiSTI command line (except the sample date for an ancient genome) are indices of these merged time intervals.

To convert indices into generations or years, use `utils/calc_time.py`:

Example

```bash
./utils/calc_time.py genome1.psmc genome2.psmc
```

Example output

```text
0      0
1      235
2      335
3      492
4      699
5      771
...
```

Interpretation  
Setting the split time to 3 in the command line corresponds to 492 generations/years ago.

## MiSTI command line

### Minimal command line
```bash
./MiSTI.py genome1.psmc genome2.psmc mi.sfs split_time
```

Inputs and meaning

Table: required arguments

| Argument | Meaning |
|---|---|
| `genome1.psmc` | PSMC output for genome 1 |
| `genome2.psmc` | PSMC output for genome 2 |
| `mi.sfs` | joint SFS in MiSTI format |
| `split_time` | split time index (merged PSMC interval index) |

Critical requirement  
PSMC files must be supplied in the same order as the genomes appear in the SFS.

### Migration bands
Migration bands are specified using the `-mi` option. Each `-mi` is followed by five parameters (interpreted backward in time).

Table: `-mi` parameters

| Position | Parameter | Description |
|---|---|---|
| 1 | source population | 1 or 2; lineages from this population may migrate to the other population. Interpreted forward in time, this corresponds to the receiving population |
| 2 | start time | time interval index at which migration starts |
| 3 | end time | time interval index at which migration ends; must be larger than the start time |
| 4 | migration rate | initial migration rate value |
| 5 | optimised flag | 0 for fixed, 1 for optimised (included in numerical optimisation using Nelder–Mead or basin-hopping) |

Example

```bash
-mi 2 5 12 2.7 0 -mi 1 2 10 0.8 1
```

Interpretation of the example  
This command defines two migration bands. The first band is migration from population 2 between time intervals 5 and 12 with a fixed migration rate of 2.7. The second band is migration from population 1 between time intervals 2 and 10 with an initial rate of 0.8; this rate is optimised.

Schematic timeline (indices)

```text
time index:  0   2         5        10    12
            |---|---------|---------|-----|
band 2:         [=========)                 source pop 1, optimised
band 1:                   [==============)  source pop 2, fixed
```

### Pulse migration
Pulse migration is specified similarly to migration bands but represents an instantaneous event. Use `-pu` followed by four parameters.

Table: `-pu` parameters

| Position | Parameter | Description |
|---|---|---|
| 1 | source population | source population |
| 2 | time | time interval index of the event |
| 3 | fraction | fraction of lineages moving to the other population |
| 4 | optimised flag | 0 for fixed, 1 for optimised |

Constraint  
Two pulse migrations in opposite directions cannot occur at the same time. If this is required, please request support for it.

### Other command line parameters

Table: additional options

| Option | Meaning |
|---|---|
| `-o results.mi` | write an output file that can be used by the plotting script |
| `-wd PATH` | set the working directory; PSMC and SFS are read from this directory and output is written there |
| `-uf` | treat the SFS as unfolded; genomes must be polarised by ancestral state prior to SFS computation |
| `--sdate` | sampling date of the ancient genome in years (the ancient genome is always the second genome) |
| `--hetloss`, `-hl` | loss of heterozygosity for the two genomes (default 0, i.e., no loss) |
| `-rd NUM` | read PSMC round NUM; by default, the last round is used |

### Advanced and experimental command line parameters

Table: advanced options

| Option | Meaning |
|---|---|
| `-mth NUM` | mixture threshold, NUM between 0 and 1 |
| `-tol NUM` | precision parameter for numerical optimisation |
| `--cpfit` | approximate effective population sizes per interval by fitting coalescence probabilities or by fitting expected coalescence times (default) |
| `--trueEPS` | treat input coalescence rates as true effective population sizes (typically for simulated data when true EPS are known) |
| `--nosmooth` | disable smoothing |
| `--bsSize`, `-bs` | bootstrap the joint SFS to estimate the variance of the composite log-likelihood for optimised parameter values |

## Plotting results
Use `MiSTIPlot.py` to visualise results.

Minimal command line

```bash
./MiSTIPlot.py genome1.psmc genome2.psmc results.mi
```

Here, `results.mi` is produced by `MiSTI.py` using the `-o` option. For an ancient genome, provide `--sdate` as well. Additional options similar to `MiSTI.py` are supported. Run the script without arguments to list all available options.

Output  
The script generates a PDF file (default `plot.pdf`; use `-o` to change the name) containing five panels.

Panel description

1. A PSMC-like plot including the original PSMC trajectories and the effective population sizes inferred by MiSTI. A vertical line indicates the split time.
2. Probabilities that the two lineages (from the first and the second genomes) are in population 1, conditional on not having coalesced.
3. The same as panel 2, but for population 2.
4. The probability that the two lineages are in different populations.
5. The probability that the two lineages of each genome have not coalesced by the given time.

## Setting time units
In the coalescent models used by PSMC and MiSTI, all parameters, including times, are expressed in units scaled by a default effective population size N0. To rescale time to generations or years, specify:

Table: parameters required for rescaling

| Parameter | Description |
|---|---|
| mutation rate | per base pair per generation per haplotype |
| generation time | time per generation |
| bin size | the bin size used in PSMC |

You may edit `setunits.txt`, or create a new file with the same format and provide it using `--funits custom_file_name.txt`.

## Running simulations
A shell script `run_sim.sh` is provided for simulations. Install the msHOT_lite simulator (Heng Li: https://github.com/lh3/foreign/tree/master/msHOT-lite). In lines 10 and 11 of the script, set the paths to msHOT-lite and PSMC on your machine. GNU Parallel is also used (line 40) to run PSMC on two genomes in parallel.

The following command runs the ms simulation, generates PSMC and SFS files, and places them in the folder `simulated_scenario`.

Example

```bash
./run_sim.sh simulated_scenario "4 100 -t 15000 -r 1920 30000000 -l -I 2 2 2 -n 1 10 -n 2 4.5 -eN 0.025 0.2 -ej 0.045 2 1 -eN 0.175 3 -eN 0.625 1.8 -eN 3 3.2 -eN 8 5.5"
```

Then you can run MiSTI

```bash
./MiSTI.py ms2g1.psmc ms2g2.psmc sim.jafs 22 -o output.mi -uf
```

## MiSTI and GNU Parallel
If many models must be evaluated (for example, when estimating the split time), GNU Parallel is recommended: https://www.gnu.org/software/parallel/.

Example with a model containing four migration bands (migration rates change at time interval `{mc}`):

```bash
parallel --header : -j 20 ./MiSTI.py ms2g1.psmc ms2g2.psmc sim.jafs {st} -uf -o res.mi -wd data_sim/migr0_5 -mi 1 0 {mc} {mi1} 0 -mi 2 0 {mc} {mi2} 0 -mi 1 {mc} {st} {mi3} 0 -mi 2 {mc} {st} {mi4} 0 >> res.out ::: st 20 21 22 23 24 25 ::: mc 8 9 10 11 12 ::: mi1 00.0 00.5 02.0 05.0 ::: mi2 05.0 10.0 15.0 20.0 ::: mi3 00.0 00.5 02.0 05.0 ::: mi4 01.0 05.0 10.0 15.0
```

The second command sorts scenarios by the best likelihood:

```bash
grep "migration rates = " res.out | sort -t "=" -rn -k5 | less
```
