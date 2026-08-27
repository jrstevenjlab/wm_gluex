# KHyperon Standalone Fit

This directory holds the standalone AmpTools fit binary and its companion plotter:

- `fit.cc` builds the `fit` executable
- `khyperon_plotter.cc` builds the `khyperon_plotter` executable
- `FSRootDataReader.cc` and `FSRootDataReader.h` provide the flat FSRoot reader used by the GlueX skims
- `KHyperon.cc` and `KHyperon.h` provide the custom amplitude class
- `Makefile` builds both binaries against the local AmpTools copy in `third_party/AmpTools-0.11.0`

## Environment

Before building or running, source the ROOT environment from the repository root:

```bash
source /sciclone/home/jrstevens01/wm_gluex/setup_root.csh
```

That script sets:

- `WM_GLUEX`
- ROOT 6.32.08 via `thisroot.sh`
- `FSROOT`

The standalone Makefile is already configured to use:

- ROOT 6.32.08
- the local AmpTools source tree in `third_party/AmpTools-0.11.0`
- the old libstdc++ ABI required by the AmpTools build

## Rebuild The Local AmpTools Copy

If you modify anything under `third_party/AmpTools-0.11.0`, rebuild the local AmpTools archives first:

```bash
make -C /sciclone/home/jrstevens01/wm_gluex/third_party/AmpTools-0.11.0 clean
make -C /sciclone/home/jrstevens01/wm_gluex/third_party/AmpTools-0.11.0 -j4
```

This includes the local `AmpToolsInterface::clear()` cleanup fix that avoids the post-fit bus error.

## Build The Standalone Binaries

From the repository root:

```bash
make -C analysis/kplamb/fitting/standalone clean
make -C analysis/kplamb/fitting/standalone -j4
```

This produces:

- `fit`
- `khyperon_plotter`

## How The Fit Binary Works

The `fit` executable accepts:

```bash
fit -c <config file> [-r <trials>] [-m <max iterations>] [-s <seed file>]
```

The options mean:

- `-c` config file to read
- `-r` number of randomized start-value trials
- `-m` maximum Minuit iterations
- `-s` optional seed file to write after the fit

The standalone `fit` binary loads the data and MC once, then reuses the same
`AmpToolsInterface` for every trial.  Each trial starts from a fresh randomized
parameter set: the AmpPars named by `parRange` are randomized within their
configured ranges, and the floating production coefficients are jittered as
well.  The event files are loaded only once.

The config template used by this analysis is:

- `analysis/kplamb/fitting/fit_khyperon.cfg`

It uses placeholders such as `FITNAME`, `CHANNEL`, `TREE`, `INDIR`, and `TBIN`. The driver script fills those in before each fit.

## Run A Single Fit Manually

1. Create a config file from `analysis/kplamb/fitting/fit_khyperon.cfg`.
2. Run the standalone fit:

```bash
analysis/kplamb/fitting/standalone/fit -c <your-config>.cfg -r 25 -m 1000000
```

3. Optionally write a seed file:

```bash
analysis/kplamb/fitting/standalone/fit -c <your-config>.cfg -r 25 -m 1000000 -s <seed-file>.txt
```

The fit output file name is taken from the `fit` line in the config. For the current analysis that is usually `khyperon.fit`.

When you use `-r <trials>`, the best trial is rerun once at the end so the
final `.fit` file and optional seed file are written from the winning start
point, not from a freshly reconstructed fit object.

## Inspect Fit Results

After the fit completes, use the plotter to print the summary:

```bash
analysis/kplamb/fitting/standalone/khyperon_plotter <fit-file>.fit [output.root]
```

The plotter prints:

- the likelihood
- MINUIT status
- parameter values and errors
- amplitude production coefficients and scales

It also writes a compact parameter CSV next to the fit file.  For the default
fit name this is:

```bash
khyperon_parameters.csv
```

The CSV contains `Sigma`, `Ox`, `P`, `T`, and `Oz` with their fit errors.

It also writes a ROOT file containing:

- `covariance_matrix`
- `correlation_matrix`
- `MHyp`
- `cosThetax`
- `cosThetay`
- `cosThetaz`
- `cosThetaHyp`
- `phiHyp`
- `PhiHyp`

If you do not pass an explicit output name, the plotter writes the fit-file stem with `.plots.root` appended. For example, `khyperon.fit` becomes `khyperon.plots.root`.

The plotter is built on the AmpTools `PlotGenerator` class, so the accepted-MC and generated-MC histograms are filled with the intensity weights from the fitted model instead of raw event counts.  The combined comparison histograms in the ROOT file are still named to match `analysis/kplamb/fitting/plot_plotter.C`.

The ROOT file also includes global comparison histograms for the full fit sample, using the same naming pattern as `analysis/kplamb/fitting/plot_plotter.C` but without any reaction prefix:

- `MHyperondat`
- `MHyperonbkgnd`
- `MHyperonacc`
- `MHyperongen`
- `cosThetaHypdat`
- `phiHypdat`
- `cosThetaXdat`
- `cosThetaYdat`
- `cosThetaZdat`
- `Phidat`

The accepted-MC and generated-MC histograms are scaled to the combined data yield before being written, so they overlay more cleanly in the macro.

You can draw the parameter matrices with ROOT, for example:

```bash
root -l <fit-file>.plots.root
root [0] covariance_matrix->Draw("COLZ TEXT")
root [1] correlation_matrix->Draw("COLZ TEXT")
```

To use the comparison macro, pass the reaction prefix you want to inspect, for example:

```bash
root -l -b -q 'analysis/kplamb/fitting/plot_plotter.C("khyperon.plots.root")'
```

## Run The Full t-Bin Loop

The batch driver is:

```bash
python3 analysis/kplamb/fitting/runFits.py
```

It does three things for each t bin:

- fills in the config template
- runs the standalone `fit`
- runs `khyperon_plotter` on the resulting `.fit` file
- runs `analysis/kplamb/fitting/plot_plotter.C` on the resulting `.plots.root` file to make `fit.pdf` and `fit_eff.pdf`

The plotter will also write a matching `*.plots.root` file in each fit directory.

After all t bins finish, `runFits.py` also calls:

```bash
python3 analysis/kplamb/fitting/plot_khyperon_params_vs_t.py <fit-directory>
```

This collects the per-bin `khyperon_parameters.csv` files and writes:

- `khyperon_parameters_vs_t.csv`
- `khyperon_parameters_vs_t.pdf`
- `khyperon_parameters_vs_t.root`

The combined CSV is the source table used for the PyROOT parameter-vs-t plots.
If the fit directories are named `kplambda/khyperon` instead of
`kplamb/khyperon`, pass that directory explicitly to the script.

Before running it, you may want to edit the settings near the top of `runFits.py`:

- `inDir`
- `channel`
- `tree`
- `lowT` and `highT`
- `numRand`
- `maxCalls`

The script writes its output under the current working directory, in a tree like:

```text
kplamb/khyperon/t_<low>_<high>/
```

## Notes

- The standalone reader expects the flat FSRoot branch layout used by the GlueX skims.
- If you see the old cleanup bus error again, rebuild the local AmpTools tree and then rebuild the standalone binaries.
- The local build currently uses the ROOT 6.32.08 install at `/sciclone/home/jrstevens01/builds/root/root-6.32.08`.
