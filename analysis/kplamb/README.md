# kplamb FSRoot Helpers

This directory contains ROOT macros for the `kplamb` topology, corresponding to
`gamma p -> K+ Lambda` with `Lambda -> p pi-`:

- `skim.C` skims flattened FSRoot trees by run period and makes best-chi2 skim files.
- `plot_kplamb.C` makes quick data/phase-space MC comparison plots from those skim files.
- `skimForAmpTools.C` prepares t-binned signal, sideband, and phase-space MC trees for AmpTools.

## Final State Numbering

FSRoot particle indices used by these macros:

| Index | Particle |
| --- | --- |
| 1 | `Lambda -> p pi-` |
| 1a | proton from `Lambda` |
| 1b | `pi-` from `Lambda` |
| 2 | `K+` |

The macros register the FSRoot mode code `100000000_100000` under the `kplamb`
category. Expressions such as `MASS(1)` refer to the reconstructed Lambda, while
`RMASS2(GLUEXTARGET,B,-1,-2)` uses the Lambda and `K+` as the reconstructed
final state.

## Inputs

The flattened trees are read from:

```bash
/sciclone/gluex10/flattened/kpi_gx1_pwa/ver01/tree_kplamb__B4_M18/
```

For each period, `skim.C` reads data from:

```bash
/sciclone/gluex10/flattened/kpi_gx1_pwa/ver01/tree_kplamb__B4_M18/data/tree_kplamb__B4_M18_FSROOT_%02d*.root
```

It also reads phase-space MC from:

```bash
/sciclone/gluex10/flattened/kpi_gx1_pwa/ver01/tree_kplamb__B4_M18/*phasespace*/tree_kplamb__B4_M18_FSROOT_%02d*.root
```

BGGEN input globs and output lines are present but commented out in `skim.C`.
Uncomment the matching `FSModeTree::skimTree(...)`, `createChi2RankingTree(...)`,
and best-chi2 skim lines if BGGEN skims are needed.

## Prerequisites

The `FSROOT` environment variable must be set, and FSRoot libraries must be
available before running the ROOT macros.

Check with:

```bash
echo $FSROOT
```

If this prints nothing, source your GlueX/FSRoot environment setup first.

## Skim

Run from this directory:

```bash
root -l -b -q 'skim.C()'
```

The default `skim()` function runs periods 03, 04, 05, and 07. For each period,
the active data and phase-space MC skims produce:

- `tree_kplamb__B4_M18_GENERAL_SKIM_%02d.root`
- `tree_kplamb__B4_M18_BestChi2_SKIM_%02d.root`
- `tree_kplamb__B4_M18_PSMC_GENERAL_SKIM_%02d.root`
- `tree_kplamb__B4_M18_PSMC_BestChi2_SKIM_%02d.root`

The first skim stage applies the loose skim cuts `eBeam,chi2,rf`. The second
stage creates a `Chi2Rank` friend tree and keeps only `Chi2Rank==1`.

## Plot

After the needed best-chi2 skim files exist, run:

```bash
root -l -b -q 'plot_kplamb.C()'
```

To include BGGEN component overlays:

```bash
root -l -b -q 'plot_kplamb.C(true)'
```

This writes PDF snapshots in `plots/`, including `plots/kplamb_cuts.pdf`.
When BGGEN plotting is enabled, it also writes `plots/kplamb_bggen.pdf`.

The plotting macro reads:

- data: `tree_kplamb__B4_M18_BestChi2_SKIM*.root`
- phase-space MC: `tree_kplamb__B4_M18_PSMC_BestChi2_SKIM*.root`

Update `FND_BGGEN` near the top of `plot_kplamb.C` before using the BGGEN option
with dedicated BGGEN skim files.

## AmpTools Skims

To make the default MC-only AmpTools inputs:

```bash
root -l -b -q 'skimForAmpTools.C()'
```

To also build data signal and sideband trees split by polarization:

```bash
root -l -b -q 'skimForAmpTools.C(false)'
```

Outputs are written under `out/t_<min>_<max>/` for the t bins defined in
`skimForAmpTools.C`. The signal and sideband selection uses the `Lambda2DSB`
sideband-weight definition around the reconstructed Lambda mass.
