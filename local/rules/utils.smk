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
    shell: 
        "cat < {input} | tab2xlsx > {output}"

rule gz2xlsx:
    """Decompress a gzipped tab file and convert to Excel (.xlsx)."""
    input: 
        "{file}.gz"
    output: 
        "{file}.xlsx"
    shell: 
        "zcat < {input} | tab2xlsx > {output}"

rule tab2xls:
    """Convert a TSV/Tab file to the older Excel (.xls) format."""
    input: 
        "{file}"
    output: 
        "{file}.xls"
    shell: 
        "tab2xls < {input} > {output}"


# =============================================================================
# 2. HEADER MANIPULATION
# =============================================================================

rule add_header:
    """Add a header row to a compressed file by transposing the 2nd column."""
    input: 
        "{path}.gz"
    output: 
        "{path}.header_added.gz"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; zcat {input} ) | gzip > {output}"

rule header_add:
    """Add a header row to a standard text file by transposing the 2nd column."""
    input: 
        "{file}"
    output: 
        "{file}.header_added"
    shell: 
        "(bawk -M {input} | cut -f 2 | transpose; cat {input} ) > {output}"

# =============================================================================
# 3. SEQUENCE MANIPULATION (FASTQ/FASTA)
# =============================================================================

rule get_fa:
    """Convert a gzipped FASTQ file to a gzipped FASTA file."""
    input: 
        "{file}.fastq.gz"
    output: 
        "{file}.fa.gz"
    shell: 
        "zcat {input} | fastq2tab | enumerate_rows | cut -f 1,3 | tab2fasta -s | gzip > {output}"
