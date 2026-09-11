# kpkm Flattening Helper

This directory contains small helper scripts to convert GlueX analysis trees into flattened FSRoot trees.

## Data:

The analysis trees for a given run period can be found on the Private Wiki at links like the one below which refers to the RunPeriod-2018-08

https://halldweb.jlab.org/wiki-private/index.php/Fall_2018_Analysis_Launch

Typically, you want to use the latest (largest version number) of the analysis launch available, in the case of `tree_kpkm__B4` that would be version 18.

## Script: `writeData.sh`

`writeData.sh` loops over input ROOT files and runs the `FlattenForFSRoot/flatten` executable on each file.

Default settings in the script:
- `TREE=tree_kpkm__B4`
- Input directory:
  - `/cache/halld/RunPeriod-2018-08/analysis/ver18/${TREE}/merged`
- Output directory (that you can change to your own preferred location):
  - `/volatile/halld/home/jrsteven/flattened/${TREE}/data`
- Flatten options:
  - `-chi2 20`
  - `-usePolarization 1`

The output file name format is:
- `${TREE}_FSROOT_<RUN>.root`

If an output file already exists for a run, the script skips that run.

## Prerequisites

1. The input files must exist under `/cache/...` (or update `INDIR` in the script).
2. The flatten executable must exist at:
   - `~/work2/analysisGluexI/builds/hd_utilities/FlattenForFSRoot/flatten`
3. The `FSROOT` environment variable must be set (and FSRoot libraries available) before running the ROOT macros.
  - Check with: `echo $FSROOT`
  - If this prints nothing, source your GlueX/FSRoot environment setup first.

## Run

From this directory:

```bash
bash writeData.sh
```

## Notes

If needed, edit `TREE`, `INDIR`, and `OUTDIR` at the top of `writeData.sh` before running.

For example, if the files don't exist on the cache disk already, they may need to be pulled from tape.  You can search for the file stubs in the `/mss/halld/` path and then get them from tape with

`jcache get /mss/halld/RunPeriod-2018-08/analysis/ver18/tree_kpkm__B4/merged/*.root`

which will eventually place them at the same path with `mss` replaced with `cache` where they can be read from the disk.

## Example: Skim Then Plot (`skim_kpkm.C` and `plot_kpkm.C`)

The current follow-up workflow compares the `kpkm` and `pippim` hypotheses before making a `kpkm` skim:

1. Build a standard and a hybrid Chi2 ranking from both hypotheses.
2. Keep `kpkm` candidates passing the fixed Chi2 and RF cuts.
3. Compare the skimmed data and signal MC distributions.

### 1) Run the skim macro

From this directory:

```bash
root -l -b -q 'skim_kpkm.C()'
```

By default, `skim_kpkm()` calls `skim_period(5)`, corresponding to the 2018-08 sample. The input globs are currently hard-coded in `skim_period()` and use a test subset of data runs `0506*`:

- Data: `tree_kpkm__B4` and `tree_pippim__B4` under `/volatile/halld/home/jrsteven/flattened/`
- Signal MC: the corresponding flattened files under `/volatile/halld/home/jrsteven/flattened/`

The macro creates ranking friend trees for both samples. `Chi2Rank` ranks candidates by `Chi2DOF` within each `Run/Event` group, while `HybridChi2Rank` also includes the beam-energy grouping. The ranking trees provide the `Chi2Rank` and `Chi2RankGlobal` branches used by the skim and plots.

The default skim output is:

- `tree_kpkm__B4_BestChi2_SKIM_05.root`
- `tree_kpkm__B4_SIGMC_BestChi2_SKIM_05.root`

To process another period, enable the corresponding `skim_period()` call at the bottom of `skim_kpkm.C` and update the input file globs as needed.

### 2) Run the plotting macro

After the skim and ranking friend trees exist, run:

```bash
root -l -b -q 'plot_kpkm.C()'
```

`plot_kpkm.C` reads `tree_kpkm__B4_BestChi2_SKIM_*.root` and `tree_kpkm__B4_SIGMC_BestChi2_SKIM_*.root`, attaches the `Chi2Rank` friend tree, and compares data with signal MC for unused energy, production vertex, $|t|$, beam energy, missing-mass squared, fit $\chi^2$/dof, and invariant-mass distributions. Its default cuts include the RF, vertex, missing-mass, unused-energy, unused-track, beam-energy, fit-quality, and `Chi2Rank==1` selections.

The `kpkm` mass plot also shows the subset passing `Chi2RankGlobal==1`, which removes events where the competing `pippim` hypothesis has the better global fit ranking. The macro writes:

- `out_kpkm.root`, containing the selected `kpkm` and `pK-` data/MC histograms for downstream fitting.
- `plots/`, recreated at startup for plot output or interactive ROOT canvases.

## Alternative Hypotheses

The cross-hypothesis ranking requires flattened `kpkm` and `pippim` trees to be available before running the skim. Update the paths in `skim_kpkm.C` when using a full run period, a different MC sample, or a different set of runs.
