#!/usr/bin/env python3
"""Assemble the compare-app comparison report and human-readable summary.

Takes the raw ASEReadCounter TSV plus the externally-counted
total_candidate_sites, computes metrics and a verdict via metrics.py, and
writes both a machine-readable JSON report and a one-page plain-text summary.
Contains no subprocess/GATK invocation logic of its own.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from metrics import ASEParseError, ThresholdPolicy, classify, parse_ase_tsv


def _decimal_or_none(value):
    return None if value is None else str(value)


def build_report(
    *,
    dna_sample_id: str | None,
    rna_sample_id: str | None,
    total_candidate_sites: int,
    ase_tsv: str,
    policy: ThresholdPolicy = ThresholdPolicy(),
) -> dict:
    metrics = parse_ase_tsv(ase_tsv, total_candidate_sites=total_candidate_sites)
    verdict = classify(metrics, policy)

    return {
        "schema_version": "swap-compare-report.v1",
        "dna_sample_id": dna_sample_id,
        "rna_sample_id": rna_sample_id,
        "verdict": verdict.value,
        "metrics": {
            "total_candidate_sites": metrics.total_candidate_sites,
            "ase_output_rows": metrics.ase_output_rows,
            "zero_read_sites": metrics.zero_read_sites,
            "covered_sites": metrics.covered_sites,
            "coverage_breadth_pct": _decimal_or_none(metrics.serialized_coverage_breadth_pct),
            "biallelic_sites": metrics.biallelic_sites,
            "biallelic_pct": _decimal_or_none(metrics.serialized_biallelic_pct),
            "ref_reads": metrics.ref_reads,
            "alt_reads": metrics.alt_reads,
            "ref_alt_ratio": _decimal_or_none(metrics.serialized_ref_alt_ratio),
        },
        "threshold_policy": {
            "min_covered_sites": policy.min_covered_sites,
            "min_coverage_breadth": str(policy.min_coverage_breadth),
            "confirmed_biallelic_fraction": str(policy.confirmed_biallelic_fraction),
            "suspected_biallelic_fraction": str(policy.suspected_biallelic_fraction),
            "min_ref_alt_ratio": str(policy.min_ref_alt_ratio),
            "max_ref_alt_ratio": str(policy.max_ref_alt_ratio),
        },
        # Every report is advisory-only. No verdict from this app ever
        # triggers, or is permitted to trigger, automated relabelling,
        # LIMS/Epic/SampleSheet/report mutation, or any clinical action.
        "validation_only_no_clinical_action": True,
        "reassignment_recommendation": False,
    }


def render_summary(report: dict) -> str:
    metrics = report["metrics"]
    lines = [
        "RNA-swap identity QC -- comparison summary",
        "===========================================",
        f"DNA sample:  {report['dna_sample_id'] or '(unlabelled)'}",
        f"RNA sample:  {report['rna_sample_id'] or '(unlabelled)'}",
        "",
        f"Verdict: {report['verdict']}",
        "",
        f"Total candidate sites (from site VCF): {metrics['total_candidate_sites']}",
        f"ASEReadCounter output rows:             {metrics['ase_output_rows']}",
        f"Sites with zero surviving reads:        {metrics['zero_read_sites']}",
        f"Covered sites (quality-pass):            {metrics['covered_sites']}",
        f"Coverage breadth:                        {metrics['coverage_breadth_pct']}%"
        if metrics["coverage_breadth_pct"] is not None
        else "Coverage breadth:                        N/A",
        f"Biallelic sites (both alleles seen):      {metrics['biallelic_sites']}"
        f" ({metrics['biallelic_pct']}% of covered)"
        if metrics["biallelic_pct"] is not None
        else f"Biallelic sites (both alleles seen):      {metrics['biallelic_sites']}",
        f"Ref/alt read ratio:                       {metrics['ref_alt_ratio']}"
        if metrics["ref_alt_ratio"] is not None
        else "Ref/alt read ratio:                       N/A (zero alt reads)",
        "",
        "This is an advisory identity-QC signal only. No automated relabelling,",
        "reassignment, or clinical action follows from this verdict -- any",
        "SUSPECTED_SWAP or AMBIGUOUS result requires human chain-of-custody review.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ase-tsv", required=True, type=Path)
    parser.add_argument("--total-candidate-sites", required=True, type=int)
    parser.add_argument("--dna-sample-id", default=None)
    parser.add_argument("--rna-sample-id", default=None)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    args = parser.parse_args(argv)

    try:
        report = build_report(
            dna_sample_id=args.dna_sample_id or None,
            rna_sample_id=args.rna_sample_id or None,
            total_candidate_sites=args.total_candidate_sites,
            ase_tsv=args.ase_tsv.read_text(),
        )
    except (ASEParseError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 3

    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    args.output_summary.write_text(render_summary(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
