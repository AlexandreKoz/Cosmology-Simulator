#!/usr/bin/env python3
import csv
import json
import math
import pathlib
import subprocess
import sys
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parents[2]
ART = ROOT / "validation" / "artifacts" / "research_grade" / "phase3"


def load_power(path: pathlib.Path):
    data = json.loads(path.read_text())
    bins = data.get("power_spectrum", [])
    return [
        {
            "k": float(row["k_code"]),
            "p": float(row["p_code_volume"]),
            "mode_count": int(row["mode_count"]),
        }
        for row in bins
        if float(row.get("p_code_volume", 0.0)) > 0.0
    ]


def rel_l2(lhs, rhs):
    n = min(len(lhs), len(rhs))
    if n == 0:
        return float("nan")
    num = 0.0
    den = 0.0
    for i in range(n):
        dp = lhs[i]["p"] - rhs[i]["p"]
        num += dp * dp
        den += rhs[i]["p"] * rhs[i]["p"]
    return math.sqrt(num / max(den, 1.0e-30))


def summarize_scaling(csv_path: pathlib.Path):
    rows = list(csv.DictReader(csv_path.read_text().splitlines()))
    if len(rows) != 1:
        raise RuntimeError(f"expected one row in {csv_path}")
    return rows[0]


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    (ART / "correctness").mkdir(exist_ok=True)
    (ART / "scaling").mkdir(exist_ok=True)
    (ART / "metadata").mkdir(exist_ok=True)

    diagnostics_dir = ROOT / "outputs"
    low = diagnostics_dir / "phase3_cosmo_box_ps_low" / "diagnostics" / "diag_0008_science_heavy.json"
    high = diagnostics_dir / "phase3_cosmo_box_ps_high" / "diagnostics" / "diag_0008_science_heavy.json"
    fallback_low = ROOT / "validation" / "reference" / "phase3" / "power_spectrum_low.json"
    fallback_high = ROOT / "validation" / "reference" / "phase3" / "power_spectrum_high.json"

    ps_summary = {
        "observable": "power_spectrum_consistency",
        "reference_target": "same cosmological run with higher diagnostic mesh (32 vs 24)",
        "tolerance_envelope": {"relative_l2_max": 0.15},
    }
    if low.exists() and high.exists():
        low_bins = load_power(low)
        high_bins = load_power(high)
        ps_rel_l2 = rel_l2(low_bins, high_bins)
        ps_summary.update({
            "measured": {"relative_l2": ps_rel_l2, "matched_bin_count": min(len(low_bins), len(high_bins))},
            "pass": bool(ps_rel_l2 <= 0.15),
            "source_diagnostics": [str(low.relative_to(ROOT)), str(high.relative_to(ROOT))],
            "source_kind": "runtime_heavy_diagnostics",
        })
    elif fallback_low.exists() and fallback_high.exists():
        low_bins = load_power(fallback_low)
        high_bins = load_power(fallback_high)
        ps_rel_l2 = rel_l2(low_bins, high_bins)
        ps_summary.update({
            "reference_target": "deterministic in-repo power-spectrum mesh-consistency floor (8 vs 12)",
            "tolerance_envelope": {"relative_l2_max": 0.45},
            "measured": {"relative_l2": ps_rel_l2, "matched_bin_count": min(len(low_bins), len(high_bins))},
            "pass": bool(ps_rel_l2 <= 0.45),
            "source_diagnostics": [str(fallback_low.relative_to(ROOT)), str(fallback_high.relative_to(ROOT))],
            "source_kind": "checked_in_reference_floor",
        })
    else:
        ps_summary.update({
            "pass": False,
            "blocked": "missing_diagnostics_inputs",
            "required_files": [
                str(low.relative_to(ROOT)),
                str(high.relative_to(ROOT)),
                str(fallback_low.relative_to(ROOT)),
                str(fallback_high.relative_to(ROOT)),
            ],
        })

    ps_summary.update({
    })
    (ART / "correctness" / "power_spectrum_consistency.json").write_text(json.dumps(ps_summary, indent=2) + "\n")

    scaling_entries = []
    for name in [
        "pm_only_scaling_np1.csv",
        "pm_only_scaling_np2.csv",
        "tree_only_scaling_np1.csv",
        "tree_only_scaling_np2.csv",
    ]:
        p = ROOT / "validation" / "artifacts" / name
        if p.exists():
            row = summarize_scaling(p)
            row["artifact"] = str(p.relative_to(ROOT))
            scaling_entries.append(row)
        else:
            scaling_entries.append({
                "artifact": str(p.relative_to(ROOT)),
                "status": "missing_artifact",
            })

    have_pm_np2 = any(entry.get("artifact", "").endswith("pm_only_scaling_np2.csv") and entry.get("status") != "missing_artifact"
                      for entry in scaling_entries)
    have_tree_np2 = any(entry.get("artifact", "").endswith("tree_only_scaling_np2.csv") and entry.get("status") != "missing_artifact"
                        for entry in scaling_entries)

    (ART / "scaling" / "phase2_baseline_scaling_summary.json").write_text(
        json.dumps(
            {
                "classification": "performance_evidence_only",
                "note": "np1/np2 baseline only; not full strong/weak scaling certification",
                "baseline_complete": bool(have_pm_np2 and have_tree_np2),
                "entries": scaling_entries,
            },
            indent=2,
        )
        + "\n"
    )

    try:
        git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        git_sha = "archive-no-git-metadata"
    manifest = {
        "campaign": "phase3_gravity_maturity_evidence",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_sha,
        "classification": {
            "correctness": "validation/artifacts/research_grade/phase3/correctness",
            "force_accuracy": "validation/artifacts/research_grade/phase3/force_accuracy",
            "time_integration": "validation/artifacts/research_grade/phase3/time_integration",
            "scaling": "validation/artifacts/research_grade/phase3/scaling",
            "toy_regression": "validation/artifacts/*.csv (legacy toy/perf baselines)",
        },
        "observables": [
            "cosmological_box_run_health",
            "zoom_run_health",
            "power_spectrum_consistency",
            "halo_force_potential_profiles",
            "hierarchical_stepping_time_integration_accuracy",
            "strong_scaling_baseline",
            "weak_scaling_baseline",
        ],
    }
    (ART / "metadata" / "campaign_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {(ART / 'metadata' / 'campaign_manifest.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
