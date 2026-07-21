"""Unit tests for report.py -- pure, no DNAnexus/GATK dependency."""

import json

import pytest
from report import build_report, main, render_summary

HEADER = "contig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\tlowMAPQDepth\tlowBaseQDepth\trawDepth\totherBases\timproperPairs"


def _tsv(*rows: tuple[int, int, int]) -> str:
    lines = [HEADER]
    for ref, alt, total in rows:
        lines.append(f"chr1\t1000\t.\tA\tG\t{ref}\t{alt}\t{total}\t0\t0\t{total}\t0\t0")
    return "\n".join(lines) + "\n"


def test_build_report_shape_and_advisory_markers():
    tsv = _tsv(*[(50, 50, 100)] * 80)
    report = build_report(
        dna_sample_id="26167S0040",
        rna_sample_id="26167S0040",
        total_candidate_sites=100,
        ase_tsv=tsv,
    )
    assert report["schema_version"] == "swap-compare-report.v1"
    assert report["verdict"] == "CONFIRMED"
    assert report["validation_only_no_clinical_action"] is True
    assert report["reassignment_recommendation"] is False
    assert report["metrics"]["total_candidate_sites"] == 100
    assert report["metrics"]["ase_output_rows"] == 80
    assert report["metrics"]["zero_read_sites"] == 20


def test_render_summary_contains_key_fields():
    tsv = _tsv(*[(50, 50, 100)] * 80)
    report = build_report(
        dna_sample_id="DNA1", rna_sample_id="RNA1", total_candidate_sites=100, ase_tsv=tsv
    )
    summary = render_summary(report)
    assert "DNA1" in summary
    assert "RNA1" in summary
    assert "CONFIRMED" in summary
    assert "advisory identity-QC signal only" in summary


def test_main_cli_writes_valid_json_and_summary(tmp_path):
    tsv_path = tmp_path / "ase.tsv"
    tsv_path.write_text(_tsv(*[(50, 50, 100)] * 80))
    json_path = tmp_path / "report.json"
    summary_path = tmp_path / "summary.txt"

    exit_code = main(
        [
            "--ase-tsv",
            str(tsv_path),
            "--total-candidate-sites",
            "100",
            "--dna-sample-id",
            "DNA1",
            "--rna-sample-id",
            "RNA1",
            "--output-json",
            str(json_path),
            "--output-summary",
            str(summary_path),
        ]
    )
    assert exit_code == 0
    report = json.loads(json_path.read_text())
    assert report["verdict"] == "CONFIRMED"
    assert "CONFIRMED" in summary_path.read_text()


def test_main_cli_exits_cleanly_on_malformed_tsv(tmp_path):
    tsv_path = tmp_path / "bad.tsv"
    tsv_path.write_text("not\ta\tvalid\tase\ttsv\n")
    json_path = tmp_path / "report.json"
    summary_path = tmp_path / "summary.txt"

    exit_code = main(
        [
            "--ase-tsv",
            str(tsv_path),
            "--total-candidate-sites",
            "10",
            "--output-json",
            str(json_path),
            "--output-summary",
            str(summary_path),
        ]
    )
    assert exit_code == 3
    assert not json_path.exists()


def test_main_cli_exits_cleanly_on_impossible_total_candidate_sites(tmp_path):
    tsv_path = tmp_path / "ase.tsv"
    tsv_path.write_text(_tsv(*[(50, 50, 100)] * 80))
    json_path = tmp_path / "report.json"
    summary_path = tmp_path / "summary.txt"

    exit_code = main(
        [
            "--ase-tsv",
            str(tsv_path),
            "--total-candidate-sites",
            "5",  # fewer than the 80 rows actually present -- impossible
            "--output-json",
            str(json_path),
            "--output-summary",
            str(summary_path),
        ]
    )
    assert exit_code == 3
