# KHyperon Standalone Fit and Generator

Standalone AmpTools tools for measuring polarization observables
(Sigma, Ox, P, T, Oz) in `gamma p -> K+ Lambda`, `Lambda -> p pi-`:

- `fit`: fits polarized data with the `KHyperon` amplitude, using several randomized starting points.
- `khyperon_plotter`: summarizes a fit and writes histograms and a parameter CSV.
- `gen_khyperon`: Monte Carlo generator for the same reaction.

The only external dependencies are ROOT and AmpTools, and the Makefile can
download and build AmpTools for you.  No GlueX software is needed.

| File | Purpose |
| --- | --- |
| `fit.cc` | `fit` executable |
| `khyperon_plotter.cc` | `khyperon_plotter` executable |
| `gen_khyperon.cc` | `gen_khyperon` executable |
| `KHyperon.cc`, `KHyperon.h` | polarized K+ Lambda amplitude |
| `FSRootDataReader.cc`, `FSRootDataReader.h` | reader for flat FSRoot trees |
| `Makefile` | builds AmpTools (`make amptools`) and the executables (`make`) |

Config files are one directory up, in `analysis/kplamb/fitting/`:
`fit_khyperon.cfg` (template for the real-data fit), `gen_khyperon.cfg`
(generator) and `fit_khyperon_toy.cfg` (fit to generated samples).

## Quick start

From this directory, with ROOT installed (see [Prerequisites](#prerequisites)):

```bash
make amptools     # one time: downloads and builds AmpTools into ./AmpTools
make -j4          # builds fit, khyperon_plotter and gen_khyperon
```

Then check the whole chain with a closure test on generated events.  It needs
no data:

```bash
mkdir toy && cd toy

# "data" with known observables (values set in gen_khyperon.cfg)
../gen_khyperon -c ../../gen_khyperon.cfg -o gen_khyperon.root \
    -n 100000 -s 1 -tmin 0.1 -tmax 0.2 -t 3.0

# flat phase space used as both genmc and accmc (no detector acceptance)
../gen_khyperon -f -o gen_khyperon_flat.root \
    -n 400000 -s 2 -tmin 0.1 -tmax 0.2 -t 3.0

sed "s|TOYDIR|$PWD|g" ../../fit_khyperon_toy.cfg > toy.cfg
../fit -c toy.cfg -r 3
../khyperon_plotter khyperon.fit
cat khyperon_parameters.csv

# optional: PDF summaries (fit.pdf, fit_eff.pdf)
root -l -b -q '../../plot_plotter.C("khyperon.plots.root")'
```

This takes well under a minute.  The fitted values should agree with the
generated ones (Sigma 0.4, Ox −0.3, P −0.2, T 0.5, Oz 0.2) within the
statistical errors.  A reference run gave Sigma 0.404 ± 0.013,
Ox −0.299 ± 0.030, P −0.211 ± 0.007, T 0.492 ± 0.030 and Oz 0.186 ± 0.030.

## Prerequisites

- **A C++17 compiler, make and git.**
  - macOS: `xcode-select --install`
  - Debian/Ubuntu: `sudo apt install build-essential git`
- **ROOT 6.24 or newer** with `root-config` on your `PATH`.  Any of these works:
  - conda (Linux or macOS): `conda create -n root -c conda-forge root && conda activate root`
  - Homebrew (macOS): `brew install root`
  - a binary release or source build from <https://root.cern/install/>, then
    `source <root install>/bin/thisroot.sh`

Check with `root-config --version`.  The ROOT packages above include PyROOT,
which the Python plotting scripts in `analysis/kplamb/fitting/` use.

AmpTools needs no separate installation: `make amptools` clones
<https://github.com/mashephe/AmpTools> at a pinned release (`v0.15.3` by
default; override with `AMPTOOLS_VERSION=...`) and builds it in `./AmpTools`.

## Building

```bash
make amptools   # download and build AmpTools (one time)
make -j4        # build the executables
make clean      # remove the executables and object files
make clean-amptools   # remove ./AmpTools
```

To use an AmpTools build you already have, point the Makefile at its source
tree, the directory that contains `IUAmpTools/` and `lib/`:

```bash
make AMPTOOLS=/path/to/AmpTools/AmpTools
```

You can also set `AMPTOOLS` and `ROOT_CONFIG` as environment variables, for
example in a site setup script.  Other build variables:

| Variable | Default | Meaning |
| --- | --- | --- |
| `ROOT_CONFIG` | `root-config` | ROOT configuration script |
| `AMPTOOLS` | `./AmpTools/AmpTools` | AmpTools source tree with a built `lib/` |
| `CXXSTD` | `c++17` | C++ standard |
| `OLD_ABI` | unset | set `OLD_ABI=1` if your AmpTools was built with `-D_GLIBCXX_USE_CXX11_ABI=0` |

The Makefile links the separate AmpTools component libraries when they exist
(AmpTools 0.11) and otherwise the combined `libAmpTools.a` (newer releases).
It compiles with `-std=c++17` even when ROOT reports C++20, because AmpTools
does not build as C++20 on macOS.  ROOT may then print a warning that the C++
standard does not match its configuration; this is harmless here.

## Running the fit

The commands from here on assume the executables are on your `PATH`, or are
called with their full path.  `FITDIR` stands for `analysis/kplamb/fitting` in
your checkout:

```bash
export FITDIR=/path/to/wm_gluex/analysis/kplamb/fitting
export PATH=$FITDIR/standalone:$PATH
```

```bash
fit -c <config file> [-r <trials>] [-m <max iterations>] [-s <seed file>]
```

- `-c` config file to read
- `-r` number of randomized start-value trials
- `-m` maximum MINUIT iterations
- `-s` optional seed file to write after the fit

`fit` loads the data and MC once and reuses the same `AmpToolsInterface` for
every trial.  Each trial starts from a fresh randomized parameter set.  The
AmpTools parameters named by `parRange` lines are drawn within their ranges,
and the floating production coefficients are randomized too.  The best trial
is rerun once at the end, so the final `.fit` file and optional seed file come
from the winning start point.

The output file name comes from the `fit` line in the config, usually
`khyperon.fit`.

### Input data format

The data, background and MC samples are flat FSRoot trees, read by
`FSRootDataReader`.  A config line looks like:

```text
data <reaction> FSRootDataReader <file.root> <tree name> 3 [MC] [<weight file> <weight tree> <weight branch>]
```

- `3` is the number of final-state particles.
- `MC` reads the `MC`-prefixed (thrown) branches first.
- The optional weight file supplies per-event weights, for example sideband
  subtraction weights.

The particle order is beam, K+, p, pi-.  These correspond to FSRoot branch
labels `B`, `2`, `1a` and `1b`, so for example `EnPB` and `PxP2` for
reconstructed trees, or `MCEnPB` and `MCPxP2` for thrown trees.

### Real-data configuration

`analysis/kplamb/fitting/fit_khyperon.cfg` is the template for the GlueX
analysis.  It has four reactions, one per beam polarization orientation, with
data, sideband background, accepted MC and thrown MC in each t bin.  It uses
placeholders such as `FITNAME`, `CHANNEL`, `TREE`, `INDIR` and `TBIN` that
`runFits.py` fills in (see [Running all t bins](#running-all-t-bins)).  To run
one fit by hand, replace the placeholders yourself and run
`fit -c <your-config>.cfg -r 25 -m 1000000`.

The input trees come from `analysis/kplamb/skimForAmpTools.C`, which needs
FSRoot and access to the GlueX flattened trees.  See
`analysis/kplamb/README.md`.

## Inspecting fit results

```bash
khyperon_plotter <fit-file>.fit [output.root] [parameter.csv]
```

The plotter prints:

- the likelihood
- MINUIT status
- parameter values and errors
- amplitude production coefficients and scales

It writes `<stem>_parameters.csv` (`khyperon_parameters.csv` for the default
fit name), with `Sigma`, `Ox`, `P`, `T` and `Oz` and their fit errors.  It also
writes `<stem>.plots.root`, which contains:

- `covariance_matrix` and `correlation_matrix`, with axes labeled by parameter
- `MHyp`, `cosThetax`, `cosThetay`, `cosThetaz`, `cosThetaHyp`, `phiHyp` and `PhiHyp`
- comparison histograms named `<variable><category>`, where the category is
  `dat`, `bkgnd`, `acc` or `gen`, for example `MHyperondat` and `cosThetaYacc`

The plotter is built on the AmpTools `PlotGenerator` class.  The
accepted-MC and generated-MC histograms are weighted by the fitted intensity
and scaled to the data yield, so they overlay the data directly.

```bash
root -l khyperon.plots.root
root [0] correlation_matrix->Draw("COLZ TEXT")
```

`$FITDIR/plot_plotter.C` makes PDF summaries (`fit.pdf`, `fit_eff.pdf`) from
the ROOT file:

```bash
root -l -b -q "$FITDIR/plot_plotter.C(\"khyperon.plots.root\")"
```

## Running all t bins

`$FITDIR/runFits.py` runs the real-data fit in each t bin:

```bash
python3 $FITDIR/runFits.py [--indir <dir>] [--nrand 25] [--maxcalls 1000000]
```

`--indir` is the directory that contains `kplamb/out/t_<low>_<high>/`, as
written by `skimForAmpTools.C`.  It defaults to the `analysis/` directory of
this repository.

For each t bin the script:

- fills in the config template
- runs `fit`
- runs `khyperon_plotter` on the resulting `.fit` file
- runs `plot_plotter.C` to make `fit.pdf` and `fit_eff.pdf`

Output goes under the current working directory, in
`kplamb/khyperon/t_<low>_<high>/`.  After all bins finish, it runs
`plot_khyperon_params_vs_t.py`, which combines the per-bin CSV files into
`khyperon_parameters_vs_t.csv`, `.pdf` and `.root`.  The channel, tree names
and t binning are set near the top of `runFits.py`.

## Generating Monte Carlo

`gen_khyperon` generates `gamma p -> K+ Lambda`, `Lambda -> p pi-` events.
It follows the approach of `gen_amp` in
[halld_sim](https://github.com/JeffersonLab/halld_sim/tree/master/src/programs/Simulation),
but is self-contained and specific to this channel:

1. Throw the photon energy: flat (default), `1/E` bremsstrahlung
   (`-beam brem`), or from a flux histogram (`-flux <file> <hist>`), always
   restricted to `[-a, -b]`.
2. Throw `-t` in `[-tmin, -tmax]` and reject values that are kinematically
   forbidden at that photon energy.  The shape is a forward t-channel
   exponential plus an optional backward u-channel exponential (see
   [Production model](#production-model)).  Here
   `-t = -(Lambda - target)^2`, the same definition `skimForAmpTools.C` uses
   for its t bins.
3. Decay `Lambda -> p pi-` isotropically in the Lambda rest frame.
4. Accept or reject each event using the `KHyperon` intensity, evaluated by
   AmpTools in `kMCGeneration` mode.  With `-f` this step is skipped and the
   output is flat phase space.

```bash
gen_khyperon -o <output.root> [-c <config>] [-n <events>] [-s <seed>]
             [-a <Emin>] [-b <Emax>] [-beam flat|brem] [-flux <file> <hist>]
             [-tmin <val>] [-tmax <val>] [-t <slope>] [-u <slope>] [-uFrac <val>]
             [-tree <name>] [-batch <N>] [-f]
```

The defaults are `-a 8.2 -b 8.8 -tmin 0 -tmax 3 -t 0 -uFrac 0 -n 10000`, and
the tree name is `ntFSGlueX_100000000_100000`.  Run `gen_khyperon -h` for the
full list.

The generation config uses the same `KHyperon` amplitude syntax as the fit
and must define exactly one reaction.  Set the observables to the values to
generate, as in `analysis/kplamb/fitting/gen_khyperon.cfg`.  That file uses
polarization angle 1.77° (the PARA 0 orientation) and polarization fraction
0.35.  To generate another orientation, change the angle in the config.

The output tree uses the flat FSRoot layout with `MC`-prefixed branches
(`MCEnPB`, `MCEnP2` = K+, `MCEnP1a` = p, `MCEnP1b` = pi-, `MCEnP1` = Lambda).
`FSRootDataReader` can read it as `data`, `genmc` or `accmc`.  The file also
holds diagnostic histograms:

- production: `gen_EBeam`, `gen_t` (`-t`), `gen_u` (`-u = -(gamma - Lambda)^2`)
  and `gen_cosThetaCM` (K+ polar angle relative to the beam in the gamma p CM
  frame).  `gen_t` and `gen_u` use 0.01 GeV^2 bins over the full kinematic
  range, for comparing with published dsigma/dt, dsigma/du and dsigma/dcos(theta).
- decay: `gen_MLambda`, `gen_cosThetaX/Y/Z`, `gen_Phi` and `gen_cosThetaY_vs_Phi`.

The events are four-vectors only: no detector simulation or reconstruction is
applied.  For realistic acceptance, pass the output through a detector
simulation, or use reconstructed phase-space MC as `accmc`.

### Production model

The unpolarized production is modeled as a forward (t-channel) peak plus an
optional backward (u-channel) peak:

```text
dsigma/dt ~ exp(-b_t |t|) + uFrac * exp(-b_u |u|)
```

- `-t <b_t>`: t-channel slope in GeV^-2 (0 gives flat `-t`).
- `-u <b_u>`: u-channel slope in GeV^-2, with `u = (gamma - Lambda)^2`.
- `-uFrac <R>`: height of the backward peak at `u = 0` relative to the forward
  peak at `t = 0`.  The default 0 turns the u-channel term off.

Over a narrow beam-energy window this is equivalent to a linear-trajectory
Regge model `s^(2 alpha(t) - 2)`, whose slope is `b = 2 alpha' ln s`.  Across
8.2–8.8 GeV that slope changes by only about ±1%.  For reference, the K*(892)
and Lambda trajectories give `b_t ≈ 4.7`, `b_u ≈ 5.3` and `R ≈ 0.006`.

Only `b_t` matters for the polarization fits, because the u-channel peak sits
at `-t ≈ 14 GeV^2` at 8.5 GeV, far above the fit bins (`-t < 2.1 GeV^2`).
`b_t` sets how acceptance is averaged inside each t bin.  You can estimate it
from your own data: divide the data yield in each t bin by the phase-space MC
acceptance (accepted / thrown), then fit an exponential.

With `-uFrac 0` the generator throws `exp(-b_t|t|)` directly.  With
`-uFrac > 0` it throws `-t` flat in `[-tmin, -tmax]` and accepts or rejects on
the weight, using a ceiling at the peak of each term.  The window always
applies, so use for example `-tmax 20` to include the backward region.

```bash
# full angular range, flat decay angles, for comparing with previous measurements
gen_khyperon -f -t 4.7 -u 5.3 -uFrac 0.006 -tmax 20 -n 1000000 -o tu_full.root

# polarized sample in one fit bin with a forward slope only
gen_khyperon -c $FITDIR/gen_khyperon.cfg -t 4.7 -tmin 0.1 -tmax 0.2 -n 100000 -o bin.root
```

Check with `-t 4.7 -u 5.3 -uFrac 0.006 -tmax 20` (2M events, 8.2–8.8 GeV):
the fitted forward slope is 4.703 ± 0.004, the backward `-u` slope is
5.43 ± 0.06, and the fraction of events with `cos(theta_CM) < 0` is
0.00699 ± 0.00006 (0.00695 from numerical integration of the model).

## Troubleshooting

- **`root-config not found`**: ROOT is not set up in this shell.  Activate the
  conda environment or source `thisroot.sh`.
- **`no AmpTools build found`**: run `make amptools`, or pass
  `AMPTOOLS=/path/to/AmpTools/AmpTools`.
- **Undefined symbols mentioning `std::__cxx11`** at link time (Linux): your
  AmpTools build and these executables use different libstdc++ ABIs.  Rebuild
  both with the same setting, or add `OLD_ABI=1`.
- **Messages like `TSystem::ExpandFileName ... $SOMEVAR`** come from a
  `.rootrc` or environment that refers to an unset variable.  They are
  unrelated to this code and harmless.
- **A bus error at the end of a fit** was seen with an older AmpTools 0.11
  build and fixed by a cleanup patch in `AmpToolsInterface::clear()`.  It has
  not been seen with the pinned `v0.15.3`.
