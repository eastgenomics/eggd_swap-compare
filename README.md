# eggd_swap-compare

DNAnexus app: the second of two apps in the RNA-swap identity-QC reboot. Runs
GATK `ASEReadCounter` for one prepared DNA germline site list (the `site_vcf`
output of [`eggd_swap-prep-dna-vcf`](https://github.com/eastgenomics/eggd_swap-prep-dna-vcf))
against one RNA BAM, and classifies the result as an advisory identity-QC
verdict.

Design rationale: <https://cuhbioinformatics.atlassian.net/wiki/spaces/URA/pages/4753260639>

This app is deliberately narrow in scope: **one (DNA, RNA) pair per job**, no
manifest, no candidate ID, no validation-mode flag, no coordinator
dependency. Fan-out across many pairs (expected pairs, or a cross-BAM sweep on
a suspected swap) happens outside DNAnexus, by an external driver submitting
one job per pair — never inside a coordinator job tree.

## Why coverage breadth needs a correction most implementations get wrong

GATK `ASEReadCounter` is a `LocusWalker`. If a candidate site has **zero
surviving reads** at that locus (either genuinely zero coverage, or every
overlapping read removed by an engine-level filter), it does not emit a
`totalCount=0` row for that site — it silently omits it from the output TSV
entirely. Only sites where reads survive engine-level filtering but fail
`ASEReadCounter`'s own `--min-mapping-quality`/`--min-base-quality`
thresholds get an explicit zero row.

This means **the ASE output row count is not the number of candidate sites
that were offered to it**, and using it as the coverage-breadth denominator
silently inflates the reported breadth. Confirmed both empirically (against a
real Phase 1 RNA BAM: 231 real candidate sites, only 182 rows emitted, 49
sites with literally zero reads) and independently via GATK source code.

`eggd_swap-compare` therefore requires the caller to supply
`total_candidate_sites` — counted directly from `site_vcf` (`bcftools view -H
site_vcf | wc -l`), never derived from the ASE TSV — and uses that as the true
denominator for coverage breadth. `metrics.py`'s `AlleleMetrics` makes both
numbers explicit (`ase_output_rows` and `zero_read_sites`) so the gap between
them is itself a visible diagnostic, not silently absorbed.

## Pipeline

1. Count `total_candidate_sites` directly from `site_vcf`.
2. Run GATK `ASEReadCounter` with `-V site_vcf -L site_vcf` (mandatory `-L`
   restriction — an unrestricted traversal can crash on unplaced/unlocalized
   scaffold contigs and is far slower), `--min-mapping-quality 10
   --min-base-quality 10`, and the **corrected** overlap-handling flag:
   `--count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE`
   (`--count-fragments-require-same-base` does not exist in GATK and fails
   every invocation with a `USER ERROR` before comparing a single site).
3. Compute metrics and classify a verdict (see `metrics.py`).
4. Write `comparison_report_json` (machine-readable), `comparison_summary_txt`
   (human-readable), and `ase_table` (the raw `ASEReadCounter` output, for
   audit).

## Verdict thresholds (provisional)

Calibrated against the only real evidence available so far (the Phase 1
Results table and the SWAP IDENTIFIED cross-BAM search): confirmed-correct
pairings showed 83.4–96% biallelic with ref/alt ratio 0.96–1.18;
confirmed-wrong pairings showed 29.5–50.7% biallelic. This is five real
samples and four confirmed swaps — enough to show a real, large separation,
not enough to certify an exact production boundary.

| Verdict | Rule |
|---|---|
| `INCONCLUSIVE` | Covered sites < 50, or coverage breadth < 30% |
| `CONFIRMED` | Coverage breadth ≥ 30%, biallelic % ≥ 70%, ref/alt ratio in [0.7, 1.4] |
| `SUSPECTED_SWAP` | Coverage breadth ≥ 30%, biallelic % < 55% |
| `AMBIGUOUS` | Everything else — no real observation falls in this band yet, reported explicitly rather than forced into `CONFIRMED` or `SUSPECTED_SWAP` |

**Open caveat on the `30%` coverage-breadth cutoff specifically** (distinct from the biallelic-threshold caveat above): zero-read sites are deliberately included in the `coverage_breadth` denominator, because they're the strongest real signal of a degraded/failed RNA library (see the real Phase 1 sample `26167S0043`, whose library was ~10x smaller than the others).

Recomputed against all five real Phase 1 samples under the corrected denominator:

| Sample | Total candidate sites | Covered | Coverage breadth | Verdict |
|---|---|---|---|---|
| `26167S0040` | 231 | 181 | 78.35% | CONFIRMED |
| `26180S0043` | 84 | 74 | 88.10% | CONFIRMED |
| `26167S0009` | 257 | 237 | 92.22% | SUSPECTED_SWAP (own nominated pairing) |
| `26163S0038` | 232 | 223 | 96.12% | SUSPECTED_SWAP (own nominated pairing) |
| `26167S0043` | 163 | 25 | **15.34%** | INCONCLUSIVE (known genuinely degraded library) |

Every sample with a real informative verdict clears `30%` by a wide margin (lowest is 78.35%); only the one known-degraded library collapses to 15.34% and correctly triggers `INCONCLUSIVE` rather than a swap call. This is real support that the corrected (larger) denominator does not drag a biologically well-covered sample's breadth down toward the threshold in practice. It remains five real samples, though -- not a large or diverse enough dataset to certify `0.30` as a final production cutoff, and the underlying concern (tissue-specific non-expression at some capture-region loci lowering breadth for reasons unrelated to sample identity) has not been ruled out mechanistically, only shown not to bite on this particular dataset so far.

Every report carries `validation_only_no_clinical_action: true` and
`reassignment_recommendation: false`. No verdict from this app ever triggers,
or is permitted to trigger, automated relabelling, LIMS/Epic/SampleSheet/report
mutation, or any clinical action — a `SUSPECTED_SWAP` or `AMBIGUOUS` result is
an escalation to a human for chain-of-custody investigation.

## Inputs

| Input | Class | Notes |
|---|---|---|
| `site_vcf` / `site_vcf_tbi` | file | Output of `eggd_swap-prep-dna-vcf` for one DNA sample |
| `rna_bam` / `rna_bam_bai` | file | PanCan RNA BAM (GRCh38) for the candidate/suspected specimen |
| `reference_fasta` / `_fai` / `reference_dict` | file | Must be the exact same GRCh38 reference used to prepare `site_vcf`. Must be uncompressed -- a bgzipped FASTA would also need a `.gzi` companion index, which this app has no input for |
| `gatk_jar` | file | GATK jar, supplied directly — no asset-selection indirection |
| `dna_sample_id` / `rna_sample_id` | string, optional | Labels only — never used for pairing logic |

`bcftools` comes from the team's formal, Approved `htslib_suite_asset` v1.22
(`record-J1YBvy049yKpP7kk1j4ggqxZ`, [DI-2083](https://cuhbioinformatics.atlassian.net/browse/DI-2083)),
declared under `assetDepends`.

## Outputs

| Output | Class |
|---|---|
| `comparison_report_json` | file |
| `comparison_summary_txt` | file |
| `ase_table` | file (raw `ASEReadCounter` TSV) |

## Tests

```
python3 -m pytest tests/ -v
```

Pure Python, no DNAnexus/GATK dependency — covers TSV parsing/validation,
the zero-read-sites regression this app exists to get right, and all four
verdict classification paths.
