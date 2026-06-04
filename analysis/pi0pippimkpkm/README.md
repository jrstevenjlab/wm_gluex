# pi0pippimkpkm FSRoot Helpers

This directory contains ROOT macros for the `pi0pippimkpkm` final state:

- `skim_pi0pippimkpkm.C` skims flattened FSRoot trees and makes best-chi2 skim files.
- `plot_pi0pippimkpkm.C` makes quick data/MC comparison plots from those skim files.

## Final State Numbering

FSRoot particle indices for this final state:

| Index | Particle |
| --- | --- |
| 1 | proton |
| 2 | `K+` |
| 3 | `K-` |
| 4 | `pi0` |
| 5 | `pi+` |
| 6 | `pi-` |

## Inputs

The flattened trees are already on `/sciclone`:

```bash
/sciclone/gluex10/flattened/phi2pi_gx1/tree_pi0pippimkpkm__B4/
```

The data files are read from:

```bash
/sciclone/gluex10/flattened/phi2pi_gx1/tree_pi0pippimkpkm__B4/data/tree_pi0pippimkpkm__B4_FSROOT_0*.root
```

The skim macro also has commented-out globs for BGGEN, signal MC, eta MC, and omega MC under the same `/sciclone` base. Uncomment the matching `FSModeTree::skimTree(...)`, `createChi2RankingTree(...)`, and final best-chi2 skim lines when those samples are available.

## Prerequisites

The `FSROOT` environment variable must be set, and FSRoot libraries must be available before running the ROOT macros.

Check with:

```bash
echo $FSROOT
```

If this prints nothing, source your GlueX/FSRoot environment setup first.

## Skim

Run from this directory:

```bash
root -l -b -q 'skim_pi0pippimkpkm.C()'
```

The default active skim reads the data sample and produces:

- `tree_pi0pippimkpkm__B4_GENERAL_SKIM.root`
- `tree_pi0pippimkpkm__B4_BestChi2_SKIM.root`

Optional sample outputs, once their input globs exist and their lines are uncommented, are:

- `tree_pi0pippimkpkm__B4_BGGEN_BestChi2_SKIM.root`
- `tree_pi0pippimkpkm__B4_SIGMC_BestChi2_SKIM.root`
- `tree_pi0pippimkpkm__B4_ETAMC_BestChi2_SKIM.root`
- `tree_pi0pippimkpkm__B4_OMEGAMC_BestChi2_SKIM.root`

## Plot

After the needed skim files exist, run:

```bash
root -l -b -q 'plot_pi0pippimkpkm.C()'
```

To include BGGEN component overlays:

```bash
root -l -b -q 'plot_pi0pippimkpkm.C(true)'
```

This writes `hist_pi0pippimkpkm.root` with saved histograms for downstream plotting/fitting.

Note that the plotting macro reads the data, SIGMC, ETAMC, and OMEGAMC skim filenames defined near the top of `plot_pi0pippimkpkm.C`. Create those skims first, or adjust the `FND_*` inputs before plotting a reduced set.
