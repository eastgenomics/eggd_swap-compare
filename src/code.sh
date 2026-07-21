#!/bin/bash
# eggd_swap-compare — run GATK ASEReadCounter for one prepared DNA site list
# against one RNA BAM, and classify the result as an advisory identity-QC
# verdict. Deliberately simple: no asset manifest, no candidate_id, no
# validation_mode, no coordinator dependency. A single (DNA, RNA) pair per
# job -- fan-out across many pairs happens outside this app.
set -eo pipefail

main() {
    local work_dir
    work_dir="$(mktemp -d)"
    trap 'rm -rf "$work_dir"' EXIT

    dx-download-all-inputs --parallel

    local site_vcf_path rna_bam_path gatk_jar_path
    site_vcf_path="$HOME/in/site_vcf/$(ls "$HOME/in/site_vcf")"
    rna_bam_path="$HOME/in/rna_bam/$(ls "$HOME/in/rna_bam")"
    gatk_jar_path="$HOME/in/gatk_jar/$(ls "$HOME/in/gatk_jar")"

    # dx-download-all-inputs keeps each declared input in its own subfolder,
    # but GATK/samtools require the .tbi/.bai to sit next to the file they
    # index and the .fai/.dict to sit next to the reference fasta with a
    # matching basename -- build both bundles by symlink.
    mkdir -p "$work_dir/site" "$work_dir/bam" "$work_dir/ref"
    ln -s "$site_vcf_path" "$work_dir/site/site.vcf.gz"
    ln -s "$HOME/in/site_vcf_tbi/$(ls "$HOME/in/site_vcf_tbi")" "$work_dir/site/site.vcf.gz.tbi"
    ln -s "$rna_bam_path" "$work_dir/bam/rna.bam"
    ln -s "$HOME/in/rna_bam_bai/$(ls "$HOME/in/rna_bam_bai")" "$work_dir/bam/rna.bam.bai"
    ln -s "$HOME/in/reference_fasta/$(ls "$HOME/in/reference_fasta")" "$work_dir/ref/reference.fasta"
    ln -s "$HOME/in/reference_fasta_fai/$(ls "$HOME/in/reference_fasta_fai")" "$work_dir/ref/reference.fasta.fai"
    ln -s "$HOME/in/reference_dict/$(ls "$HOME/in/reference_dict")" "$work_dir/ref/reference.dict"

    local site_vcf="$work_dir/site/site.vcf.gz"
    local rna_bam="$work_dir/bam/rna.bam"
    local reference_fasta="$work_dir/ref/reference.fasta"

    # --- Count the true candidate-site denominator directly from the site
    # VCF. Never derived from the ASEReadCounter TSV row count: GATK's
    # LocusWalker traversal silently omits a locus with zero surviving reads
    # rather than emitting a totalCount=0 row, so the TSV row count already
    # excludes exactly the sites with the strongest "no expression" signal.
    local total_candidate_sites
    total_candidate_sites="$(bcftools view -H "$site_vcf" | wc -l)"
    echo "Total candidate sites (from site_vcf): ${total_candidate_sites}"

    # --- Run GATK ASEReadCounter -------------------------------------------
    # -L is mandatory, not optional: without it, an unrestricted traversal
    # can crash on unplaced/unlocalized scaffold contigs and takes far longer.
    # --count-overlap-reads-handling is the corrected flag name;
    # --count-fragments-require-same-base does not exist in GATK and would
    # fail every invocation with a USER ERROR before comparing a single site.
    local ase_output="$work_dir/ase_output.tsv"
    java -XX:MaxRAMPercentage=70.0 -jar "$gatk_jar_path" ASEReadCounter \
        -R "$reference_fasta" \
        -I "$rna_bam" \
        -V "$site_vcf" \
        -L "$site_vcf" \
        -O "$ase_output" \
        --min-mapping-quality 10 \
        --min-base-quality 10 \
        --count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE

    local ase_output_rows
    ase_output_rows="$(($(wc -l <"$ase_output") - 1))"
    echo "ASEReadCounter output rows: ${ase_output_rows} (of ${total_candidate_sites} candidate sites)"

    # --- Build the comparison report and summary ----------------------------
    local report_json="$work_dir/comparison_report.json"
    local summary_txt="$work_dir/comparison_summary.txt"
    PYTHONPATH="/home/dnanexus" python3 /home/dnanexus/report.py \
        --ase-tsv "$ase_output" \
        --total-candidate-sites "$total_candidate_sites" \
        --dna-sample-id "${dna_sample_id:-}" \
        --rna-sample-id "${rna_sample_id:-}" \
        --output-json "$report_json" \
        --output-summary "$summary_txt"

    echo "--- comparison_summary.txt ---"
    cat "$summary_txt"

    local report_json_id summary_txt_id ase_table_id
    report_json_id="$(dx upload "$report_json" --brief)"
    summary_txt_id="$(dx upload "$summary_txt" --brief)"
    ase_table_id="$(dx upload "$ase_output" --brief)"

    dx-jobutil-add-output comparison_report_json "$report_json_id" --class=file
    dx-jobutil-add-output comparison_summary_txt "$summary_txt_id" --class=file
    dx-jobutil-add-output ase_table "$ase_table_id" --class=file
}

main
