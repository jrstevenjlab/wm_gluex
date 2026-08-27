#!/usr/bin/env python3

import argparse
import csv
import math
import re
from array import array
from pathlib import Path


PARAMETERS = ("Sigma", "Ox", "P", "T", "Oz")
T_BIN_RE = re.compile(r"^t_([+-]?\d+(?:\.\d+)?)_([+-]?\d+(?:\.\d+)?)$")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collect per-bin KHyperon parameter CSV files and plot the "
            "polarization observables versus t with PyROOT."
        )
    )
    parser.add_argument(
        "base_dir",
        nargs="?",
        help=(
            "Directory containing t_low_high subdirectories. If omitted, "
            "the script looks for kplamb/khyperon or kplambda/khyperon "
            "relative to the current directory and this script."
        ),
    )
    parser.add_argument(
        "--csv-name",
        default="khyperon_parameters.csv",
        help="Per-bin CSV filename written by khyperon_plotter.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Directory for the combined CSV, ROOT file, and PDF.",
    )
    parser.add_argument(
        "--output-csv",
        help="Combined CSV output path. Defaults to <output-dir>/khyperon_parameters_vs_t.csv.",
    )
    parser.add_argument(
        "--output-pdf",
        help="Summary PDF output path. Defaults to <output-dir>/khyperon_parameters_vs_t.pdf.",
    )
    parser.add_argument(
        "--output-root",
        help="ROOT graph output path. Defaults to <output-dir>/khyperon_parameters_vs_t.root.",
    )
    parser.add_argument("--root", help=argparse.SUPPRESS)
    parser.add_argument(
        "--x-title",
        default="t",
        help="X-axis title for the plots.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Only write the combined CSV; do not make PyROOT plots.",
    )
    return parser.parse_args()


def has_t_bin_dirs(path):
    if not path.exists() or not path.is_dir():
        return False
    return any(T_BIN_RE.match(child.name) for child in path.iterdir() if child.is_dir())


def default_base_dir():
    cwd = Path.cwd()
    script_dir = Path(__file__).resolve().parent
    candidates = (
        cwd / "kplamb" / "khyperon",
        cwd / "kplambda" / "khyperon",
        script_dir / "kplamb" / "khyperon",
        script_dir / "kplambda" / "khyperon",
        cwd,
    )

    for candidate in candidates:
        if has_t_bin_dirs(candidate):
            return candidate

    return cwd / "kplamb" / "khyperon"


def parse_float(value):
    if value is None:
        return math.nan
    value = value.strip()
    if not value:
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def format_float(value):
    if not math.isfinite(value):
        return "nan"
    return f"{value:.15g}"


def read_parameter_csv(csv_path):
    values = {}
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = row.get("parameter") or row.get("name")
            if not name:
                continue
            value = parse_float(row.get("value"))
            error = parse_float(row.get("error") or row.get("uncertainty"))
            values[name.strip()] = (value, error)
    return values


def find_parameter_csv(t_dir, csv_name):
    requested = t_dir / csv_name
    if requested.exists():
        return requested

    matches = sorted(t_dir.glob("*_parameters.csv"))
    return matches[0] if matches else requested


def collect_rows(base_dir, csv_name):
    rows = []

    for t_dir in sorted(base_dir.iterdir()):
        if not t_dir.is_dir():
            continue

        match = T_BIN_RE.match(t_dir.name)
        if not match:
            continue

        t_low = float(match.group(1))
        t_high = float(match.group(2))
        csv_path = find_parameter_csv(t_dir, csv_name)
        if not csv_path.exists():
            print(f"Skipping {t_dir}: missing {csv_path.name}")
            continue

        par_values = read_parameter_csv(csv_path)
        row = {
            "t_low": t_low,
            "t_high": t_high,
            "t_center": 0.5 * (t_low + t_high),
            "t_width": 0.5 * abs(t_high - t_low),
            "source_csv": str(csv_path),
        }

        for parameter in PARAMETERS:
            value, error = par_values.get(parameter, (math.nan, math.nan))
            row[parameter] = value
            row[f"{parameter}_error"] = error

        rows.append(row)

    return sorted(rows, key=lambda row: (row["t_low"], row["t_high"]))


def write_combined_csv(rows, output_csv):
    columns = ["t_low", "t_high", "t_center", "t_width"]
    for parameter in PARAMETERS:
        columns.extend((parameter, f"{parameter}_error"))
    columns.append("source_csv")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    column: (
                        format_float(row[column])
                        if isinstance(row[column], float)
                        else row[column]
                    )
                    for column in columns
                }
            )


def graph_points(rows, parameter):
    points = []
    for row in rows:
        value = row[parameter]
        if not math.isfinite(value):
            continue

        error = row[f"{parameter}_error"]
        points.append(
            (
                row["t_center"],
                row["t_width"],
                value,
                error if math.isfinite(error) else 0.0,
            )
        )
    return points


def y_limits(points):
    ys = [0.0]
    for _, _, value, error in points:
        ys.extend((value - error, value + error))

    low = min(ys)
    high = max(ys)
    if low == high:
        low -= 0.1
        high += 0.1
    else:
        pad = 0.12 * (high - low)
        low -= pad
        high += pad

    return low, high


def load_root():
    try:
        import ROOT
    except ImportError as exc:
        raise SystemExit(
            "Could not import PyROOT. Source the analysis environment first, "
            "for example: source /sciclone/home/jrstevens01/wm_gluex/setup_root.sh"
        ) from exc

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    return ROOT


def make_parameter_graph(ROOT, rows, parameter, x_low, x_high, x_title, color):
    points = graph_points(rows, parameter)
    if not points:
        return None, []

    x_values = array("d", [point[0] for point in points])
    x_errors = array("d", [point[1] for point in points])
    y_values = array("d", [point[2] for point in points])
    y_errors = array("d", [point[3] for point in points])

    graph = ROOT.TGraphErrors(
        len(points),
        x_values,
        y_values,
        x_errors,
        y_errors,
    )
    graph.SetName(f"graph_{parameter}")
    graph.SetTitle(f";{x_title};{parameter}")
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.9)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)

    y_low, y_high = y_limits(points)
    graph.SetMinimum(y_low)
    graph.SetMaximum(y_high)
    graph.Draw("AP")
    graph.GetXaxis().SetLimits(x_low, x_high)
    graph.GetXaxis().CenterTitle()
    graph.GetYaxis().CenterTitle()

    owned_objects = [x_values, x_errors, y_values, y_errors, graph]
    if y_low < 0.0 < y_high:
        zero_line = ROOT.TLine(x_low, 0.0, x_high, 0.0)
        zero_line.SetLineStyle(2)
        zero_line.SetLineColor(ROOT.kGray + 1)
        zero_line.Draw()
        graph.Draw("P SAME")
        owned_objects.append(zero_line)

    return graph, owned_objects


def make_parameter_plots(rows, output_pdf, output_root, x_title):
    ROOT = load_root()

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_root.parent.mkdir(parents=True, exist_ok=True)

    x_low = min(row["t_low"] for row in rows)
    x_high = max(row["t_high"] for row in rows)
    if x_low == x_high:
        x_low -= 0.1
        x_high += 0.1

    colors = {
        "Sigma": ROOT.kBlue + 1,
        "Ox": ROOT.kRed + 1,
        "P": ROOT.kGreen + 2,
        "T": ROOT.kMagenta + 1,
        "Oz": ROOT.kOrange + 7,
    }

    output_file = ROOT.TFile(str(output_root), "RECREATE")
    if output_file.IsZombie():
        raise SystemExit(f"Could not create ROOT output file: {output_root}")

    canvas = ROOT.TCanvas("canvas", "KHyperon parameters vs t", 1200, 800)
    canvas.Divide(3, 2)

    owned_objects = [canvas]
    for pad, parameter in enumerate(PARAMETERS, start=1):
        canvas.cd(pad)
        graph, graph_objects = make_parameter_graph(
            ROOT,
            rows,
            parameter,
            x_low,
            x_high,
            x_title,
            colors[parameter],
        )
        owned_objects.extend(graph_objects)

        if graph is None:
            text = ROOT.TLatex()
            text.SetNDC()
            text.SetTextSize(0.08)
            text.DrawLatex(0.18, 0.52, f"No finite {parameter} values")
            owned_objects.append(text)
            continue

        output_file.cd()
        graph.Write()

    output_file.cd()
    canvas.Write("canvas_khyperon_parameters_vs_t")
    canvas.SaveAs(str(output_pdf))
    output_file.Close()
    return owned_objects


def main():
    args = parse_args()

    base_dir = Path(args.base_dir).expanduser() if args.base_dir else default_base_dir()
    base_dir = base_dir.resolve()
    if not base_dir.exists():
        raise SystemExit(f"Base directory does not exist: {base_dir}")

    output_dir = Path(args.output_dir).expanduser() if args.output_dir else base_dir
    output_dir = output_dir.resolve()
    output_csv = (
        Path(args.output_csv).expanduser()
        if args.output_csv
        else output_dir / "khyperon_parameters_vs_t.csv"
    )
    output_pdf = (
        Path(args.output_pdf).expanduser()
        if args.output_pdf
        else output_dir / "khyperon_parameters_vs_t.pdf"
    )
    output_root = (
        Path(args.output_root).expanduser()
        if args.output_root
        else output_dir / "khyperon_parameters_vs_t.root"
    )

    rows = collect_rows(base_dir, args.csv_name)
    if not rows:
        raise SystemExit(f"No per-bin parameter CSV files found under {base_dir}")

    write_combined_csv(rows, output_csv)
    print(f"Wrote combined parameter CSV: {output_csv}")

    if args.no_plots:
        return

    make_parameter_plots(rows, output_pdf, output_root, args.x_title)
    print(f"Wrote parameter-vs-t PDF: {output_pdf}")
    print(f"Wrote parameter-vs-t ROOT file: {output_root}")


if __name__ == "__main__":
    main()
