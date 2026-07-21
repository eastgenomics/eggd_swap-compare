#!/usr/bin/env python3
"""Pure ASEReadCounter parsing, metrics, and verdict classification.

No DNAnexus, subprocess, or filesystem I/O -- this module only ever consumes
already-read text and an externally-supplied ``total_candidate_sites`` count.

Critical design point, carried over from eggd_swap-prep-dna-vcf's testing:
GATK ASEReadCounter is a LocusWalker and silently omits a locus from its
output entirely when zero reads survive engine-level filtering at that
position -- it does not emit a ``totalCount=0`` row for it. This means the
number of rows in the ASE TSV is *not* the number of candidate sites that
were offered to it, and must never be used as the coverage-breadth
denominator. The true denominator is the number of records in the prepared
site VCF (``site_vcf``), which the caller must count independently (e.g. via
``bcftools view -H site_vcf | wc -l``) and pass in explicitly as
``total_candidate_sites``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from decimal import ROUND_HALF_UP, Decimal
from enum import StrEnum
from io import StringIO
from typing import Iterable

_REQUIRED_ASE_COLUMNS = frozenset({"refCount", "altCount", "totalCount"})


class Verdict(StrEnum):
    """The only permitted comparison verdicts. No ERROR case: a tool or
    parsing failure in this single-pair app fails the whole job closed
    (non-zero exit, no output) rather than being contained as a cell value,
    since there is no batch/matrix to keep partially alive."""

    CONFIRMED = "CONFIRMED"
    SUSPECTED_SWAP = "SUSPECTED_SWAP"
    AMBIGUOUS = "AMBIGUOUS"
    INCONCLUSIVE = "INCONCLUSIVE"


class ASEParseError(ValueError):
    """An ASE TSV violates the complete, non-negative count contract."""


@dataclass(frozen=True)
class ThresholdPolicy:
    """Provisional release-1.0.0 verdict thresholds.

    Calibrated against the only real evidence available so far (Phase 1
    Results table and the SWAP IDENTIFIED cross-BAM search): confirmed-correct
    pairings showed 83.4-96% biallelic with ref/alt ratio 0.96-1.18;
    confirmed-wrong pairings showed 29.5-50.7% biallelic. This is five real
    samples and four confirmed swaps -- enough to show a real, large
    separation, not enough to certify an exact production boundary. The gap
    between ``suspected_biallelic_fraction`` and ``confirmed_biallelic_fraction``
    is deliberately reported as AMBIGUOUS rather than forced into CONFIRMED or
    SUSPECTED_SWAP, since no real observation falls in that band yet.

    ``min_coverage_breadth`` carries a distinct, more basic caveat: it was
    carried over unchanged from an earlier design that computed
    ``coverage_breadth`` against the wrong denominator (the ASEReadCounter
    output row count, which silently excludes zero-read loci -- see the
    module docstring). Zero-read sites are deliberately included in the
    corrected denominator (``total_candidate_sites``) because they are the
    strongest real signal of a degraded/failed RNA library (see the real
    Phase 1 sample ``26167S0043``).

    Recomputed against all five real Phase 1 samples under the corrected
    denominator (not just one): the four samples with any real informative
    verdict showed coverage breadth of 78.35% (``26167S0040``, CONFIRMED),
    88.10% (``26180S0043``, CONFIRMED), 92.22% (``26167S0009``,
    SUSPECTED_SWAP against its own nominated pairing), and 96.12%
    (``26163S0038``, SUSPECTED_SWAP against its own nominated pairing) -- all
    comfortably clear of ``0.30``. Only the one genuinely degraded sample
    (``26167S0043``) collapsed to 15.34%, correctly triggering INCONCLUSIVE
    rather than a swap call. This is real support that the corrected
    denominator does not drag a biologically well-covered sample's breadth
    down toward the threshold in practice -- but it remains five real
    samples, one of which is the only INCONCLUSIVE case observed, so it is
    still not a large enough or diverse enough dataset to certify ``0.30`` as
    a final production cutoff.
    """

    min_covered_sites: int = 50
    min_coverage_breadth: Decimal = Decimal("0.30")
    confirmed_biallelic_fraction: Decimal = Decimal("0.70")
    suspected_biallelic_fraction: Decimal = Decimal("0.55")
    min_ref_alt_ratio: Decimal = Decimal("0.7")
    max_ref_alt_ratio: Decimal = Decimal("1.4")


@dataclass(frozen=True)
class ASERow:
    """The validated ASE counts needed by the aggregate evidence policy."""

    ref_count: int
    alt_count: int
    total_count: int


@dataclass(frozen=True)
class AlleleMetrics:
    """Unrounded aggregate evidence from a complete ASEReadCounter run.

    ``total_candidate_sites`` is the true site-list denominator (supplied by
    the caller, never derived from the ASE TSV). ``ase_output_rows`` is how
    many rows GATK actually emitted -- always <= total_candidate_sites, and
    the gap between them (``zero_read_sites``) is itself diagnostic: it is
    the count of sites with literally no surviving reads at all.
    """

    total_candidate_sites: int
    ase_output_rows: int
    covered_sites: int
    biallelic_sites: int
    ref_reads: int
    alt_reads: int

    @property
    def zero_read_sites(self) -> int:
        """Sites offered to ASEReadCounter that produced no output row at all."""

        return self.total_candidate_sites - self.ase_output_rows

    @property
    def coverage_breadth(self) -> Decimal | None:
        """Return covered_sites / total_candidate_sites, or ``None`` if zero."""

        if self.total_candidate_sites == 0:
            return None
        return Decimal(self.covered_sites) / Decimal(self.total_candidate_sites)

    @property
    def biallelic_fraction(self) -> Decimal | None:
        """Return biallelic_sites / covered_sites, or ``None`` when covered is zero."""

        if self.covered_sites == 0:
            return None
        return Decimal(self.biallelic_sites) / Decimal(self.covered_sites)

    @property
    def ref_alt_ratio(self) -> Decimal | None:
        """Return the unrounded reference/alternate ratio, or ``None`` when A is zero."""

        if self.alt_reads == 0:
            return None
        return Decimal(self.ref_reads) / Decimal(self.alt_reads)

    @property
    def serialized_coverage_breadth_pct(self) -> Decimal | None:
        breadth = self.coverage_breadth
        return None if breadth is None else _round_half_up(breadth * Decimal(100), Decimal("0.01"))

    @property
    def serialized_biallelic_pct(self) -> Decimal | None:
        fraction = self.biallelic_fraction
        return None if fraction is None else _round_half_up(fraction * Decimal(100), Decimal("0.01"))

    @property
    def serialized_ref_alt_ratio(self) -> Decimal | None:
        return _round_half_up(self.ref_alt_ratio, Decimal("0.001"))


def _round_half_up(value: Decimal | None, quantum: Decimal) -> Decimal | None:
    return None if value is None else value.quantize(quantum, rounding=ROUND_HALF_UP)


def _strict_nonnegative_integer(value: str, *, field: str, row_number: int) -> int:
    cleaned = value.strip()
    if not cleaned or not cleaned.isdigit():
        raise ASEParseError(f"ASE_PARSE_ERROR: row {row_number} has invalid {field}")
    return int(cleaned)


def _parse_rows(tsv: str) -> tuple[ASERow, ...]:
    if not isinstance(tsv, str) or not tsv:
        raise ASEParseError("ASE_PARSE_ERROR: output is empty or not text")

    reader = csv.reader(StringIO(tsv), delimiter="\t")
    try:
        header = next(reader)
    except StopIteration as error:
        raise ASEParseError("ASE_PARSE_ERROR: output has no header") from error
    if header:
        header[0] = header[0].lstrip("\ufeff")
    if not header or len(header) != len(set(header)) or not _REQUIRED_ASE_COLUMNS.issubset(header):
        raise ASEParseError("ASE_PARSE_ERROR: required ASE columns are absent or ambiguous")

    column_index = {name: index for index, name in enumerate(header)}
    rows: list[ASERow] = []
    for row_number, fields in enumerate(reader, start=2):
        if not fields or all(not value.strip() for value in fields):
            continue
        if len(fields) != len(header):
            raise ASEParseError(f"ASE_PARSE_ERROR: row {row_number} has an invalid field count")
        ref_count = _strict_nonnegative_integer(
            fields[column_index["refCount"]], field="refCount", row_number=row_number
        )
        alt_count = _strict_nonnegative_integer(
            fields[column_index["altCount"]], field="altCount", row_number=row_number
        )
        total_count = _strict_nonnegative_integer(
            fields[column_index["totalCount"]], field="totalCount", row_number=row_number
        )
        if ref_count + alt_count > total_count:
            raise ASEParseError(
                f"ASE_PARSE_ERROR: row {row_number} has allele counts greater than totalCount"
            )
        rows.append(ASERow(ref_count=ref_count, alt_count=alt_count, total_count=total_count))
    return tuple(rows)


def calculate_metrics(rows: Iterable[ASERow], *, total_candidate_sites: int) -> AlleleMetrics:
    """Calculate aggregate evidence from validated ASE rows.

    ``total_candidate_sites`` must be supplied by the caller (a direct count
    of the input site VCF) -- see the module docstring for why it must never
    be derived from ``len(rows)``.
    """

    if total_candidate_sites < 0:
        raise ValueError("total_candidate_sites must be non-negative")

    materialized_rows = tuple(rows)
    if len(materialized_rows) > total_candidate_sites:
        raise ValueError(
            "ASE output has more rows than total_candidate_sites -- "
            "the wrong site VCF / count was supplied"
        )
    covered_rows = tuple(row for row in materialized_rows if row.total_count > 0)
    return AlleleMetrics(
        total_candidate_sites=total_candidate_sites,
        ase_output_rows=len(materialized_rows),
        covered_sites=len(covered_rows),
        biallelic_sites=sum(row.ref_count > 0 and row.alt_count > 0 for row in covered_rows),
        ref_reads=sum(row.ref_count for row in covered_rows),
        alt_reads=sum(row.alt_count for row in covered_rows),
    )


def parse_ase_tsv(tsv: str, *, total_candidate_sites: int) -> AlleleMetrics:
    """Parse a complete ASE TSV and calculate its aggregate metrics."""

    return calculate_metrics(_parse_rows(tsv), total_candidate_sites=total_candidate_sites)


def classify(metrics: AlleleMetrics, policy: ThresholdPolicy = ThresholdPolicy()) -> Verdict:
    """Classify unrounded evidence under the provisional release-1.0.0 policy.

    See ThresholdPolicy's docstring for the real-data basis and its explicit
    caveats. This is a provisional starting point for validation, not a
    finished clinical cut-off.
    """

    breadth = metrics.coverage_breadth
    if (
        metrics.covered_sites < policy.min_covered_sites
        or breadth is None
        or breadth < policy.min_coverage_breadth
    ):
        return Verdict.INCONCLUSIVE

    biallelic_fraction = metrics.biallelic_fraction
    ratio = metrics.ref_alt_ratio

    if (
        biallelic_fraction is not None
        and biallelic_fraction >= policy.confirmed_biallelic_fraction
        and ratio is not None
        and policy.min_ref_alt_ratio <= ratio <= policy.max_ref_alt_ratio
    ):
        return Verdict.CONFIRMED

    if biallelic_fraction is not None and biallelic_fraction < policy.suspected_biallelic_fraction:
        return Verdict.SUSPECTED_SWAP

    return Verdict.AMBIGUOUS
