# local/rules/preflight.smk
import pandas as pd

# ── metadata.txt ─────────────────────────────────────────────────────────
_meta_ok = True

if not os.path.exists("metadata.txt"):
    logger.warning("[preflight] metadata.txt not found!")
    _meta_ok = False

elif os.path.getsize("metadata.txt") == 0:
    logger.warning("[preflight] metadata.txt is empty!")
    _meta_ok = False

if _meta_ok:
    _meta = pd.read_table("metadata.txt")

    if "sample" not in _meta.columns:
        logger.warning(
            "[preflight] metadata.txt has no 'sample' column!"
            "Columns found: %s", list(_meta.columns)
        )
    else:
        _meta_samples = set(_meta["sample"].astype(str))
        _fastq_only   = set(SAMPLES) - _meta_samples
        _meta_only    = _meta_samples - set(SAMPLES)

        if _fastq_only:
            logger.warning(
                "[preflight] Samples in fastq/ but missing from metadata.txt: %s",
                sorted(_fastq_only)
            )
        if _meta_only:
            logger.warning(
                "[preflight] Samples in metadata.txt with no fastq found: %s",
                sorted(_meta_only)
            )

        if not _fastq_only and not _meta_only:
            logger.info("[preflight] metadata OK — %d sample(s): %s", len(SAMPLES), SAMPLES)



_valid_layouts  = {"SINGLE", "PAIRED"}
_valid_trimmers = {"fastp", "trimgalore"}

if config.get("LAYOUT") not in _valid_layouts:
    raise WorkflowError(
        f"\n[preflight] Invalid LAYOUT: '{config.get('LAYOUT')}'. "
        f"Must be one of: {_valid_layouts}"
    )
if config.get("TRIMMER") not in _valid_trimmers:
    raise WorkflowError(
        f"\n[preflight] Invalid TRIMMER: '{config.get('TRIMMER')}'. "
        f"Must be one of: {_valid_trimmers}"
    )

_valid_aligners = {"star", "bwa"}
if config.get("aligner") not in _valid_aligners:
    raise WorkflowError(f"[preflight] Invalid aligner: '{config.get('aligner')}'. Must be one of: {_valid_aligners}")


_valid_star_modes = {"single_pass", "two_pass_manual", "two_pass_basic"}
if config.get("STAR_MODE") not in _valid_star_modes:
    raise WorkflowError(
        f"\n[preflight] Invalid STAR_MODE: '{config.get('STAR_MODE')}'. "
        f"Must be one of: {_valid_star_modes}"
    )

if not isinstance(config.get("USE_K2_DAEMON", False), bool):
    raise WorkflowError(
        "[preflight] USE_K2_DAEMON must be true or false (not a string). "
        f"Got: {config['USE_K2_DAEMON']!r}"
    )

if not isinstance(config.get("SAVE_UNCLASSIFIED", False), bool):
    raise WorkflowError("[preflight] SAVE_UNCLASSIFIED must be true or false")