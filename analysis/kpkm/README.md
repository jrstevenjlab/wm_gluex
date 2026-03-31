# kpkm Flattening Helper

This directory contains small helper scripts to convert GlueX analysis trees into flattened FSRoot trees.

## Data:

The analysis trees for a given run period can be found on the Private Wiki at links like the one below which referes to the runPeriod-2018-08

https://halldweb.jlab.org/wiki-private/index.php/Fall_2018_Analysis_Launch

Typically, you want to use the latest (largest version number) of the analysis launch available, in the case of `tree_kpkm__B4` that would be version 13.

## Script: `writeData.sh`

`writeData.sh` loops over input ROOT files and runs the `FlattenForFSRoot/flatten` executable on each file.

Default settings in the script:
- `TREE=tree_kpkm__B4`
- Input directory:
  - `/cache/halld/RunPeriod-2018-08/analysis/ver13/${TREE}/merged`
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

`jcache get /mss/halld/RunPeriod-2018-08/analysis/ver13/tree_kpkm__B4/merged/*.root`

which will eventually place them at the same path with `mss` replaced with `cache` where they can be read from the disk.

## Example: Skim Then Plot (`skim_kpkm.C` and `plot_kpkm.C`)

An example follow-up workflow in this directory is:

1. Skim a subset of events from the flattened trees.
2. Plot and compare key distributions from the skimmed output.

### 1) Run the skim macro

From this directory:

```bash
root -l -b -q 'skim_kpkm.C()'
```

This macro is set up to produce reduced files such as:
- `tree_kpkm__B4_BestChi2_SKIM.root`
- `tree_kpkm__B4_BGGEN_BestChi2_SKIM.root`
- `tree_kpkm__B4_PSMC_BestChi2_SKIM.root`

Note: several skim steps inside `skim_kpkm.C` are currently commented out. Uncomment the `FSModeTree::skimTree(...)` and (if needed) `createChi2RankingTree(...)` lines you want to run for your sample before executing.

### 2) Run the plotting macro on the skim

After the skim files exist, run:

```bash
root -l -b -q 'plot_kpkm.C()'
```

To include BGGEN component overlays:

```bash
root -l -b -q 'plot_kpkm.C(true)'
```

This creates:
- `out_kpkm.root` with saved histograms for downstream plotting/fitting.
- A `plots/` directory (recreated each run) for generated plot outputs.
