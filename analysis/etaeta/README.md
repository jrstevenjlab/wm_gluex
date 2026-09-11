# Eta Eta Flattening Helper

This directory contains helper scripts and ROOT macros for the reaction

```text
gamma p -> p eta eta
                         |
                         +-- eta -> pi0 pi+ pi-
                         |       pi0 -> gamma gamma
                         |
                         +-- eta -> gamma gamma
```

The reconstructed final state is `p pi+ pi- gamma gamma gamma gamma`. The FSRoot mode is `101_111`, with category `pi0pippimeta`, and the analysis tree is `tree_pi0pippimeta__B4_M17`.

## Data:

The analysis tree used by `writeData.sh` is the RunPeriod-2018-08 analysis version 20 tree:

```text
/cache/halld/RunPeriod-2018-08/analysis/ver20/tree_pi0pippimeta__B4_M17/merged
```

## Script: `writeData.sh`

`writeData.sh` loops over input ROOT files and runs the `FlattenForFSRoot/flatten` executable on each file.

Default settings in the script:
- `TREE=tree_pi0pippimeta__B4_M17`
- Input directory:
    - `/cache/halld/RunPeriod-2018-08/analysis/ver20/tree_pi0pippimeta__B4_M17/merged`
- Output directory:
    - `/volatile/halld/home/jrsteven/flattened/tree_pi0pippimeta__B4_M17/data`
- Flatten options:
    - `-chi2 20`
    - `-usePolarization 1`

The output file name format is:
- `tree_pi0pippimeta__B4_M17_FSROOT_<RUN>.root`

If an output file already exists for a run, the script skips that run.

## Prerequisites

1. The input files must exist under `/cache/...` (or update `INDIR` in the script).
2. The flatten executable must exist at:
     - `~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten`
3. The `FSROOT` environment variable must be set and FSRoot libraries must be available before running the ROOT macros.
     - Check with: `echo $FSROOT`
     - If this prints nothing, source the GlueX/FSRoot environment setup first.

## Run

From this directory:

```bash
bash writeData.sh
```

## Notes

If needed, edit `TREE`, `INDIR`, and `OUTDIR` at the top of `writeData.sh` before running. To use another run period or analysis launch version, update the input path and tree name consistently in `writeData.sh` and `skim.C`.

If the files do not exist on the cache disk, locate or stage the corresponding files before running the script. For example, the input files can be requested from the matching `/mss/halld/` path with `jcache get`, using the appropriate run period, version, and tree name.

## Example: Skim Then Plot (`skim.C` and `plots.C`)

The follow-up workflow reconstructs the two eta candidates in the `pi0 pi+ pi-` and `gamma gamma` decay modes:

1. Flatten the merged `tree_pi0pippimeta__B4_M17` analysis files.
2. Skim candidates with beam-energy, fit-quality, RF, and broad three-pion eta-mass cuts.
3. Create a `Chi2Rank` friend tree for the skim.
4. Plot the eta masses and diagnostic kinematic distributions.

### 1) Run the skim macro

After the flattened files exist, run:

```bash
root -l -b -q 'skim.C()'
```

The macro reads:

```text
/volatile/halld/home/jrsteven/flattened/tree_pi0pippimeta__B4_M17/data/tree_pi0pippimeta__B4_M17_FSROOT_*.root
```

It applies these general skim cuts:
- `EnPB > 8.0`
- `Chi2DOF < 20`
- `abs(RFDeltaT) < 2.004`
- `0.35 < MASS(pi0, pi+, pi-) < 0.85`

The skim output is:
- `tree_pi0pippimeta__B4_M17_GENERAL_SKIM.root`

The macro also creates the `Chi2Rank` friend tree using the RF cut. Keep this ranking output with the skim file because `plots.C` expects it.

### 2) Run the plotting macro

From the directory containing the skim and ranking files, run:

```bash
root -l -b -q 'plots.C()'
```

`plots.C` reads `tree_pi0pippimeta__B4_M17_GENERAL_SKIM.root`, attaches the `Chi2Rank` friend tree, and recreates the `plots/` directory. Its default plotting cuts include:
- `EnUnusedSh < 0.1`
- `NumUnusedTracks < 1`
- `51.2 <= ProdVz <= 78.8` cm
- `EnPB > 8.0` GeV
- `Chi2DOF < 5`
- `Chi2Rank == 1`
- `abs(-t) < 0.5` GeV$^2$
- missing-mass squared within `0.05` GeV$^2$ of zero
- `abs(MASS(eta -> gamma gamma) - 0.547) < 0.025` GeV
- `abs(MASS(eta -> pi0 pi+ pi-) - 0.547) < 0.025` GeV
- the wrong-combination `pi0` veto

The plots include unused energy, production vertex, $|t|$, beam energy, missing-mass squared, fit $\chi^2$/dof, both eta mass spectra, a two-dimensional comparison of the eta masses, and wrong-combination checks for the `pi0` veto.

The macro has an optional boolean argument for compatibility with older calls, but the current implementation does not add a separate BGGEN overlay:

```bash
root -l -b -q 'plots.C(false)'
```

## Alternative Hypotheses

This workflow is specific to the `pi0pippimeta` category and the `101_111` FSRoot mode. When using another analysis launch, run period, or final-state mode, update the tree name, input paths, category, and mode definition together in `writeData.sh`, `skim.C`, and `plots.C`.
