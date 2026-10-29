###############################################################################
# GENERIC BIOINFORMATICS UTILITY RULES
###############################################################################

# =============================================================================
# 1. FILE FORMAT CONVERSIONS (EXCEL / TEXT)
# =============================================================================

rule tab2xlsx:
    """Convert a TSV/Tab file to Excel (.xlsx) format."""
    input: 
        "{file}"
    output: 
        "{file}.xlsx"
    conda:
        "transcript_env.yaml"
    shell: 
        "tab2xlsx < {input:q} > {output:q}"

rule gz2xlsx:
    """Decompress a gzipped tab file and export it as .gz.xlsx."""
    input: 
        "{file}.gz"
    output: 
        "{file}.gz.xlsx"
    conda:
        "transcript_env.yaml"
    shell: 
        "gzip -cd {input:q} | tab2xlsx > {output:q}"

ruleorder: gz2xlsx > tab2xlsx

rule tab2xls:
    """Convert a TSV/Tab file to the older Excel (.xls) format."""
    input: 
        "{file}"
    output: 
        "{file}.xls"
    conda:
        "transcript_env.yaml"
    shell: 
        "tab2xls < {input:q} > {output:q}"


# =============================================================================
# 2. HEADER MANIPULATION
# =============================================================================

rule add_header:
    """Add a header row to a compressed file by transposing the 2nd column."""
    input: 
        "{path}.gz"
    output: 
        "{path}.header_added.gz"
    conda:
        "transcript_env.yaml"
    shell: 
        "(bawk -M {input:q} | cut -f 2 | transpose; zcat {input:q} ) | gzip > {output:q}"

rule header_add:
    """Add a header row to a standard text file by transposing the 2nd column."""
    input: 
        "{file}"
    output: 
        "{file}.header_added"
    conda:
        "transcript_env.yaml"
    shell: 
        "(bawk -M {input:q} | cut -f 2 | transpose; cat {input:q} ) > {output:q}"

# =============================================================================
# 3. SEQUENCE MANIPULATION (FASTQ/FASTA)
# =============================================================================

rule get_fa:
    """Convert a gzipped FASTQ file to a gzipped FASTA file."""
    input: 
        "{file}.fastq.gz"
    output: 
        "{file}.fa.gz"
    conda:
        "transcript_env.yaml"
    shell: 
        "zcat {input:q} | fastq2tab | enumerate_rows | cut -f 1,3 | tab2fasta -s | gzip > {output:q}"

rule check_utils_environment:
    """Verify that the utility environment and lab PATH supply every required command."""
    conda:
        "transcript_env.yaml"
    shell:
        r"""
        for tool in tab2xlsx tab2xls bawk transpose fastq2tab enumerate_rows tab2fasta gawk cut cat gzip zcat; do
            command -v "$tool" || {{ echo "Missing utility command: $tool" >&2; exit 1; }}
        done
        """
