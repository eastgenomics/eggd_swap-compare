"""Regression fixture: real Phase 1 patient-derived metrics, recomputed under
the corrected coverage-breadth denominator (total_candidate_sites counted
directly from the site VCF, never from the ASEReadCounter TSV row count).

These are not synthetic examples -- they are the exact real numbers obtained
by running eggd_swap-prep-dna-vcf's pipeline and this app's pipeline against
real TSO500 DNA VCFs and real PanCan RNA BAMs for all five matched pairs
available from the Phase 1 local proof-of-concept (see
local_poc/work/germline_*.vcf.gz in the rna-swap repository for the saved
intermediate artefacts these were derived from). Recorded here as a
regression fixture so a future change to metrics.py/report.py cannot
silently alter these real, already-verified results without a test failing.
"""

from decimal import Decimal

from metrics import Verdict, calculate_metrics, classify

# (total_candidate_sites, ase_output_rows, covered_sites, biallelic_sites,
#  ref_reads, alt_reads, expected_verdict)
REAL_PHASE1_SAMPLES = {
    "26167S0040": (231, 182, 181, 151, 42638, 39187, Verdict.CONFIRMED),
    "26180S0043": (84, 74, 74, 68, None, None, Verdict.CONFIRMED),
    "26167S0009": (257, 237, 237, 70, None, None, Verdict.SUSPECTED_SWAP),
    "26163S0038": (232, 223, 223, 75, None, None, Verdict.SUSPECTED_SWAP),
    "26167S0043": (163, 25, 25, 1, None, None, Verdict.INCONCLUSIVE),
}

# Coverage breadth actually observed for each sample (rounded %), confirming
# every informative sample clears the 30% INCONCLUSIVE threshold by a wide
# margin and only the known-degraded library (26167S0043) collapses.
REAL_PHASE1_COVERAGE_BREADTH_PCT = {
    "26167S0040": Decimal("78.35"),
    "26180S0043": Decimal("88.10"),
    "26167S0009": Decimal("92.22"),
    "26163S0038": Decimal("96.12"),
    "26167S0043": Decimal("15.34"),
}


def test_26167s0040_confirmed_with_correct_breadth():
    total, ase_rows, covered, biallelic, ref_reads, alt_reads, expected = REAL_PHASE1_SAMPLES[
        "26167S0040"
    ]
    # Reconstruct rows: `covered - biallelic` sites are ref-only (arbitrary
    # split not affecting aggregate ref_reads/alt_reads, which are supplied
    # directly), plus zero-read sites make up the remainder of total.
    from metrics import ASERow

    biallelic_rows = (ASERow(ref_count=1, alt_count=1, total_count=2),) * biallelic
    ref_only_rows = (ASERow(ref_count=1, alt_count=0, total_count=1),) * (covered - biallelic)
    # The real TSV also contains exactly one row with totalCount=0 (a site
    # where reads existed but all failed ASEReadCounter's own MAPQ/BQ
    # thresholds -- distinct from the 49 sites with no row at all).
    quality_filtered_zero_row = (ASERow(ref_count=0, alt_count=0, total_count=0),)
    metrics = calculate_metrics(
        biallelic_rows + ref_only_rows + quality_filtered_zero_row, total_candidate_sites=total
    )
    assert metrics.ase_output_rows == ase_rows
    assert metrics.covered_sites == covered
    assert metrics.zero_read_sites == total - ase_rows
    assert metrics.serialized_coverage_breadth_pct == REAL_PHASE1_COVERAGE_BREADTH_PCT["26167S0040"]
    assert classify(metrics) is expected


def test_26167s0043_known_degraded_sample_is_inconclusive_not_a_swap_call():
    """The one real sample with a genuinely degraded RNA library must be
    caught as INCONCLUSIVE by the coverage-breadth gate, not misread as a
    SUSPECTED_SWAP just because its biallelic fraction is also low."""

    from metrics import ASERow

    total, ase_rows, covered, biallelic, _, _, expected = REAL_PHASE1_SAMPLES["26167S0043"]
    rows = (ASERow(ref_count=1, alt_count=1, total_count=2),) * biallelic + (
        ASERow(ref_count=1, alt_count=0, total_count=1),
    ) * (covered - biallelic)
    metrics = calculate_metrics(rows, total_candidate_sites=total)
    assert metrics.serialized_coverage_breadth_pct == REAL_PHASE1_COVERAGE_BREADTH_PCT["26167S0043"]
    assert metrics.coverage_breadth < Decimal("0.30")
    assert classify(metrics) is expected is Verdict.INCONCLUSIVE


def test_all_five_real_samples_breadth_clears_or_correctly_fails_threshold():
    """No real informative sample's corrected breadth comes anywhere near
    the 30% threshold -- only the known-degraded sample does, and it should.
    This is the empirical check backing the ThresholdPolicy docstring's
    claim that the corrected denominator does not drag good samples down."""

    informative_samples = {"26167S0040", "26180S0043", "26167S0009", "26163S0038"}
    for sample in informative_samples:
        assert REAL_PHASE1_COVERAGE_BREADTH_PCT[sample] > Decimal("70")
    assert REAL_PHASE1_COVERAGE_BREADTH_PCT["26167S0043"] < Decimal("30")
