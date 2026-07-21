"""Unit tests for metrics.py -- pure, no DNAnexus/GATK dependency."""

from decimal import Decimal

import pytest
from metrics import (
    ASEParseError,
    ASERow,
    ThresholdPolicy,
    Verdict,
    calculate_metrics,
    classify,
    parse_ase_tsv,
)

HEADER = "contig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\tlowMAPQDepth\tlowBaseQDepth\trawDepth\totherBases\timproperPairs"


def _tsv(*rows: tuple[int, int, int]) -> str:
    lines = [HEADER]
    for ref, alt, total in rows:
        lines.append(f"chr1\t1000\t.\tA\tG\t{ref}\t{alt}\t{total}\t0\t0\t{total}\t0\t0")
    return "\n".join(lines) + "\n"


def test_parse_ase_tsv_basic_valid():
    tsv = _tsv((10, 10, 20), (5, 0, 5))
    metrics = parse_ase_tsv(tsv, total_candidate_sites=2)
    assert metrics.ase_output_rows == 2
    assert metrics.covered_sites == 2
    assert metrics.biallelic_sites == 1
    assert metrics.ref_reads == 15
    assert metrics.alt_reads == 10


def test_parse_ase_tsv_empty_raises():
    with pytest.raises(ASEParseError):
        parse_ase_tsv("", total_candidate_sites=0)


def test_parse_ase_tsv_missing_columns_raises():
    with pytest.raises(ASEParseError):
        parse_ase_tsv("contig\tposition\n1\t2\n", total_candidate_sites=1)


def test_parse_ase_tsv_alleles_exceed_total_raises():
    tsv = _tsv((10, 10, 15))  # ref+alt=20 > total=15
    with pytest.raises(ASEParseError):
        parse_ase_tsv(tsv, total_candidate_sites=1)


def test_parse_ase_tsv_wrong_field_count_raises():
    bad = HEADER + "\nchr1\t1000\t.\tA\tG\t10\t10\n"  # missing trailing columns
    with pytest.raises(ASEParseError):
        parse_ase_tsv(bad, total_candidate_sites=1)


def test_calculate_metrics_zero_read_sites_reflects_true_denominator():
    """The central regression this app exists to get right: ASEReadCounter
    silently omits zero-read loci, so total_candidate_sites (supplied
    externally) must exceed ase_output_rows whenever sites had no reads."""

    rows = (ASERow(ref_count=5, alt_count=5, total_count=10),)
    metrics = calculate_metrics(rows, total_candidate_sites=231)
    assert metrics.ase_output_rows == 1
    assert metrics.zero_read_sites == 230
    assert metrics.total_candidate_sites == 231


def test_calculate_metrics_rejects_impossible_total_candidate_sites():
    rows = (ASERow(ref_count=1, alt_count=1, total_count=2),) * 5
    with pytest.raises(ValueError):
        calculate_metrics(rows, total_candidate_sites=2)


def test_coverage_breadth_uses_total_candidate_sites_not_ase_rows():
    """Reproduces the real Phase 1 sample 26167S0040 numbers: 231 candidate
    sites, 182 ASE rows, 181 covered -- coverage breadth must be 181/231
    (78.35%), not 181/182 (99.5%, the pre-fix, misleadingly high figure)."""

    rows = (ASERow(ref_count=1, alt_count=1, total_count=2),) * 181 + (
        ASERow(ref_count=0, alt_count=0, total_count=0),
    )
    metrics = calculate_metrics(rows, total_candidate_sites=231)
    assert metrics.covered_sites == 181
    assert metrics.coverage_breadth == Decimal(181) / Decimal(231)
    assert metrics.serialized_coverage_breadth_pct == Decimal("78.35")


def test_ref_alt_ratio_and_biallelic_fraction_none_when_denominator_zero():
    metrics = calculate_metrics((), total_candidate_sites=10)
    assert metrics.coverage_breadth == Decimal(0)
    assert metrics.biallelic_fraction is None
    assert metrics.ref_alt_ratio is None


def test_classify_confirmed():
    rows = (ASERow(ref_count=50, alt_count=50, total_count=100),) * 80
    metrics = calculate_metrics(rows, total_candidate_sites=100)
    assert classify(metrics) is Verdict.CONFIRMED


def test_classify_suspected_swap():
    het_rows = (ASERow(ref_count=50, alt_count=50, total_count=100),) * 30
    ref_only_rows = (ASERow(ref_count=50, alt_count=0, total_count=50),) * 70
    metrics = calculate_metrics(het_rows + ref_only_rows, total_candidate_sites=100)
    assert metrics.biallelic_fraction == Decimal("0.30")
    assert classify(metrics) is Verdict.SUSPECTED_SWAP


def test_classify_inconclusive_low_breadth():
    rows = (ASERow(ref_count=50, alt_count=50, total_count=100),) * 20
    metrics = calculate_metrics(rows, total_candidate_sites=100)  # breadth = 20%
    assert classify(metrics) is Verdict.INCONCLUSIVE


def test_classify_inconclusive_low_covered_sites():
    rows = (ASERow(ref_count=50, alt_count=50, total_count=100),) * 40
    metrics = calculate_metrics(rows, total_candidate_sites=45)  # breadth 89%, covered 40 < 50
    assert classify(metrics) is Verdict.INCONCLUSIVE


def test_classify_ambiguous_between_thresholds():
    het_rows = (ASERow(ref_count=50, alt_count=50, total_count=100),) * 60
    ref_only_rows = (ASERow(ref_count=50, alt_count=0, total_count=50),) * 40
    metrics = calculate_metrics(het_rows + ref_only_rows, total_candidate_sites=100)
    assert metrics.biallelic_fraction == Decimal("0.60")  # between 0.55 and 0.70
    assert classify(metrics) is Verdict.AMBIGUOUS


def test_classify_custom_policy():
    rows = (ASERow(ref_count=50, alt_count=50, total_count=100),) * 60
    metrics = calculate_metrics(rows, total_candidate_sites=100)
    strict_policy = ThresholdPolicy(min_covered_sites=100)
    assert classify(metrics, strict_policy) is Verdict.INCONCLUSIVE
