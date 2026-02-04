#!/usr/bin/env python3
"""
ABBA_BABA_gp.py

Compute ABBA/BABA-weighted site contributions (based on population allele
frequencies) stratified by whether the test population (P2) is called
heterozygous or homozygous at the site, and (optionally) aggregate those
contributions by the genotype posterior (GP) value of the corresponding
P2 genotype.

This script is designed for biallelic SNPs in a (possibly gzipped) VCF.

Inputs
------
1) VCF/VCF.GZ with standard FORMAT fields. GT is required; GP is optional.
2) Metadata file with two whitespace-separated columns:
      <individual_id> <population_label>

You must provide four population labels (P1, P2, P3, P4). Only individuals
present in both metadata and VCF header are used.

Method (site-level)
-------------------
- For each site, compute ALT allele frequency for each population (P1..P4)
  from called genotypes (GT). Missing genotypes are ignored.
- Only biallelic SNPs are used.
- Sites are partitioned based on whether P2 genotypes are heterozygous or
  homozygous (0/0, 1/1; also supports phased calls).
- ABBA/BABA “values” follow the same expressions as in the original script,
  and are evaluated separately depending on whether P4 is fixed REF (freq=0)
  or fixed ALT (freq=1). Sites where P4 is not fixed (0<freq<1) are skipped
  to match the original behavior.
- Contributions are kept if they exceed a threshold (het vs hom).

Output
------
TSV with three columns:
    GP    Occurrences    Pattern

Where "Occurrences" is the sum of ABBA/BABA values. For sites without GP
(or if GP missing for a sample), contributions go to GP="N/A".

Notes
-----
- If multiple P2 individuals exist, each contributes independently: the same
  site may contribute multiple times (once per P2 individual) if thresholds
  are met.
- If GP is present, we use the GP corresponding to the observed GT:
      0/0 -> GP[0], 0/1 -> GP[1], 1/1 -> GP[2]
- Robustly handles missing/partial FORMAT fields.

Example
-------
python abba_baba_gp.py \
  -i data.vcf.gz -m meta.txt -o out.tsv \
  --p1 AFR --p2 TEST --p3 EUR --p4 CHIMP

"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


# ----------------------------- parsing utilities -----------------------------


def open_text_maybe_gzip(path: str):
    """Open plain text or gzipped text file for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def is_biallelic_snp(ref: str, alt: str) -> bool:
    """True iff variant is a biallelic SNP with single-base REF and ALT."""
    if ref == "." or alt == ".":
        return False
    if "," in alt:  # multiallelic
        return False
    return len(ref) == 1 and len(alt) == 1


def parse_format_indices(fmt: str) -> Dict[str, int]:
    """Map FORMAT field keys to their indices."""
    keys = fmt.split(":")
    return {k: i for i, k in enumerate(keys)}


def canonical_gt(gt: str) -> Optional[str]:
    """
    Normalize GT strings. Returns one of:
        "0/0", "0/1", "1/1"
    or None if missing/unrecognized.
    """
    if gt in (".", "./.", ".|."):
        return None
    # accept 0|1 etc.
    gt = gt.replace("|", "/")
    if gt in ("0/0", "0/1", "1/0", "1/1"):
        if gt == "1/0":
            return "0/1"
        return gt
    return None


def gt_alt_allele_count(gt: str) -> Optional[int]:
    """
    Return ALT allele count in {0,1,2} for a diploid GT, or None if missing.
    """
    cgt = canonical_gt(gt)
    if cgt is None:
        return None
    if cgt == "0/0":
        return 0
    if cgt == "0/1":
        return 1
    if cgt == "1/1":
        return 2
    return None


def selected_gp_for_gt(gt: str, gp_triplet: List[float]) -> Optional[float]:
    """
    Select GP element corresponding to normalized GT:
        0/0 -> GP[0], 0/1 -> GP[1], 1/1 -> GP[2]
    """
    cgt = canonical_gt(gt)
    if cgt is None:
        return None
    if len(gp_triplet) < 3:
        return None
    if cgt == "0/0":
        return gp_triplet[0]
    if cgt == "0/1":
        return gp_triplet[1]
    if cgt == "1/1":
        return gp_triplet[2]
    return None


# ----------------------------- core data model ------------------------------


@dataclass(frozen=True)
class PopSpec:
    """Holds the list of sample indices for a population."""
    name: str
    indices: Tuple[int, ...]


Pattern = str


def load_populations_from_metadata(
    metadata_path: str,
    pop_names: Iterable[str],
) -> Dict[str, List[str]]:
    """
    Read metadata file with columns: IND POP
    Return dict: pop -> list of individuals requested for that pop.
    """
    wanted = set(pop_names)
    pop_to_inds: Dict[str, List[str]] = {p: [] for p in pop_names}

    with open(metadata_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            ind, pop = parts[0], parts[1]
            if pop in wanted:
                pop_to_inds[pop].append(ind)

    return pop_to_inds


def map_samples_to_indices(
    header_fields: List[str],
    pop_to_inds: Dict[str, List[str]],
    logger: logging.Logger,
) -> Dict[str, PopSpec]:
    """
    Convert population -> list of individual IDs into population -> indices in VCF columns.
    Logs QC warnings for missing individuals and raises ValueError if a population has 0 samples.
    """
    pop_specs: Dict[str, PopSpec] = {}

    for pop, inds in pop_to_inds.items():
        idxs: List[int] = []
        for ind in inds:
            if ind in header_fields:
                idxs.append(header_fields.index(ind))
            else:
                logger.warning("QC: %s from %s does not appear in the VCF header.", ind, pop)

        if not idxs:
            raise ValueError(f"FATAL: no {pop} samples appear in the VCF header.")

        pop_specs[pop] = PopSpec(name=pop, indices=tuple(idxs))

    return pop_specs


def population_alt_freq(
    spline: List[str],
    pop: PopSpec,
    pos_gt: int,
    logger: logging.Logger,
) -> Optional[float]:
    """
    Compute ALT allele frequency for a population using GT calls.
    Missing GT calls are ignored; returns None if no callable genotypes.
    """
    alt_sum = 0
    allele_denom = 0

    for col_idx in pop.indices:
        fields = spline[col_idx].split(":")
        if pos_gt >= len(fields):
            continue
        gt = fields[pos_gt]
        ac = gt_alt_allele_count(gt)
        if ac is None:
            continue
        alt_sum += ac
        allele_denom += 2

    if allele_denom == 0:
        logger.debug("QC: no callable GTs for population %s at %s:%s", pop.name, spline[0], spline[1])
        return None

    return alt_sum / allele_denom


def compute_abba_baba(
    f1: float,
    f2: float,
    f3: float,
    f4: float,
) -> Optional[Tuple[float, float, str]]:
    """
    Return (abba_value, baba_value, anc_state) where anc_state is "Anc0" or "Anc1".
    Matches original branching:
      - if f4 == 0.0 -> Anc0 formulas
      - if f4 == 1.0 -> Anc1 formulas
      - else -> None (skip)
    """
    if f4 == 0.0:
        abba = (1 - f1) * f2 * f3
        baba = f1 * (1 - f2) * f3
        return abba, baba, "Anc0"
    if f4 == 1.0:
        abba = f1 * (1 - f2) * (1 - f3)
        baba = (1 - f1) * f2 * (1 - f3)
        return abba, baba, "Anc1"
    return None


def process_vcf(
    input_vcf: str,
    metadata_path: str,
    p1: str,
    p2: str,
    p3: str,
    p4: str,
    het_threshold: float,
    hom_threshold: float,
    logger: logging.Logger,
) -> Dict[Pattern, List[Tuple[float, Optional[float]]]]:
    """
    Process VCF and return a dict pattern -> list of (value, selected_gp_or_None).

    Patterns:
      ABBA_het_Anc0, BABA_het_Anc0, ABBA_het_Anc1, BABA_het_Anc1,
      ABBA_hom_Anc0, BABA_hom_Anc0, ABBA_hom_Anc1, BABA_hom_Anc1
    """
    pop_to_inds = load_populations_from_metadata(metadata_path, [p1, p2, p3, p4])

    patterns: Dict[Pattern, List[Tuple[float, Optional[float]]]] = {
        "ABBA_het_Anc0": [],
        "BABA_het_Anc0": [],
        "ABBA_het_Anc1": [],
        "BABA_het_Anc1": [],
        "ABBA_hom_Anc0": [],
        "BABA_hom_Anc0": [],
        "ABBA_hom_Anc1": [],
        "BABA_hom_Anc1": [],
    }

    with open_text_maybe_gzip(input_vcf) as fh:
        pop_specs: Optional[Dict[str, PopSpec]] = None

        for line in fh:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.strip().split()
                # VCF columns: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT samples...
                pop_specs = map_samples_to_indices(header, pop_to_inds, logger)
                logger.info("Loaded samples: %s", {k: len(v.indices) for k, v in pop_specs.items()})
                continue

            if line.startswith("#"):
                continue

            if pop_specs is None:
                raise ValueError("VCF header not found before variant lines.")

            spline = line.strip().split()
            chrom, pos, ref, alt = spline[0], spline[1], spline[3], spline[4]

            if not is_biallelic_snp(ref, alt):
                logger.debug("QC: %s:%s is not a biallelic SNP, skipping.", chrom, pos)
                continue

            fmt = spline[8]
            fmt_idx = parse_format_indices(fmt)
            if "GT" not in fmt_idx:
                logger.debug("QC: %s:%s missing GT in FORMAT, skipping.", chrom, pos)
                continue
            pos_gt = fmt_idx["GT"]
            pos_gp = fmt_idx.get("GP", -1)

            # Frequencies for P1..P4
            f1 = population_alt_freq(spline, pop_specs[p1], pos_gt, logger)
            f2 = population_alt_freq(spline, pop_specs[p2], pos_gt, logger)
            f3 = population_alt_freq(spline, pop_specs[p3], pos_gt, logger)
            f4 = population_alt_freq(spline, pop_specs[p4], pos_gt, logger)

            if any(x is None for x in (f1, f2, f3, f4)):
                # match original spirit: skip if any population cannot be computed
                continue

            abba_baba = compute_abba_baba(f1, f2, f3, f4)
            if abba_baba is None:
                continue
            abba_value, baba_value, anc_state = abba_baba

            # Evaluate per P2 individual genotype category and (optional) GP
            for col_idx in pop_specs[p2].indices:
                fields = spline[col_idx].split(":")
                if pos_gt >= len(fields):
                    continue

                gt_raw = fields[pos_gt]
                cgt = canonical_gt(gt_raw)
                if cgt is None:
                    continue

                is_het = cgt == "0/1"
                is_hom = cgt in ("0/0", "1/1")

                if not (is_het or is_hom):
                    continue

                threshold = het_threshold if is_het else hom_threshold
                category = "het" if is_het else "hom"

                # GP selection (optional)
                selected_gp: Optional[float] = None
                if pos_gp != -1 and pos_gp < len(fields):
                    gp_str = fields[pos_gp]
                    try:
                        gp_triplet = [float(x) for x in gp_str.split(",")]
                        selected_gp = selected_gp_for_gt(cgt, gp_triplet)
                    except ValueError:
                        selected_gp = None

                if abba_value > threshold:
                    patterns[f"ABBA_{category}_{anc_state}"].append((abba_value, selected_gp))
                if baba_value > threshold:
                    patterns[f"BABA_{category}_{anc_state}"].append((baba_value, selected_gp))

    return patterns


def aggregate_by_gp(
    patterns: Dict[Pattern, List[Tuple[float, Optional[float]]]]
) -> Dict[Pattern, Dict[str, float]]:
    """
    Aggregate sums per pattern keyed by GP (as string). Missing GP -> "N/A".
    Returns: pattern -> {gp_string: sum_of_values}
    """
    out: Dict[Pattern, Dict[str, float]] = {}
    for pattern, entries in patterns.items():
        sums: Dict[str, float] = defaultdict(float)
        for value, gp in entries:
            key = "N/A" if gp is None else str(gp)
            sums[key] += value
        out[pattern] = dict(sums)
    return out


def write_output_tsv(output_path: str, agg: Dict[Pattern, Dict[str, float]]) -> None:
    with open(output_path, "wt") as f:
        f.write("GP\tOccurrences\tPattern\n")
        for pattern in sorted(agg.keys()):
            gp_map = agg[pattern]
            # Write N/A last for readability
            for gp_key in sorted((k for k in gp_map.keys() if k != "N/A"), key=lambda x: float(x)):
                f.write(f"{gp_key}\t{gp_map[gp_key]}\t{pattern}\n")
            if "N/A" in gp_map:
                f.write(f"N/A\t{gp_map['N/A']}\t{pattern}\n")


# ---------------------------------- CLI -------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Aggregate ABBA/BABA site contributions by GP (optional) for P2 genotypes."
    )
    p.add_argument("-i", "--input", required=True, help="Input VCF (.vcf or .vcf.gz).")
    p.add_argument("-m", "--meta", required=True, help="Metadata file: IND POP (whitespace-separated).")
    p.add_argument("-o", "--output", required=True, help="Output TSV file.")
    p.add_argument("--p1", required=True, help="Population label for P1.")
    p.add_argument("--p2", required=True, help="Population label for P2 (test population).")
    p.add_argument("--p3", required=True, help="Population label for P3.")
    p.add_argument("--p4", required=True, help="Population label for P4 (outgroup/ancestral).")
    p.add_argument("--het-threshold", type=float, default=0.49, help="Threshold for heterozygous P2.")
    p.add_argument("--hom-threshold", type=float, default=0.99, help="Threshold for homozygous P2.")
    p.add_argument(
        "--log",
        default=None,
        help="Log file path (default: <output>.log). Use '-' for stderr only.",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging verbosity.",
    )
    return p


def setup_logger(log_path: Optional[str], level: str) -> logging.Logger:
    logger = logging.getLogger("abba_baba_gp")
    logger.setLevel(getattr(logging, level))
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")

    # Always log to stderr
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    sh.setLevel(getattr(logging, level))
    logger.addHandler(sh)

    if log_path is not None and log_path != "-":
        fh = logging.FileHandler(log_path, mode="a")
        fh.setFormatter(fmt)
        fh.setLevel(getattr(logging, level))
        logger.addHandler(fh)

    return logger


def main() -> int:
    args = build_arg_parser().parse_args()
    log_path = args.log if args.log is not None else f"{args.output}.log"
    logger = setup_logger(log_path, args.log_level)

    logger.info("Input VCF: %s", args.input)
    logger.info("Metadata:  %s", args.meta)
    logger.info("Populations: P1=%s P2=%s P3=%s P4=%s", args.p1, args.p2, args.p3, args.p4)
    logger.info("Thresholds: het=%.4f hom=%.4f", args.het_threshold, args.hom_threshold)

    patterns = process_vcf(
        input_vcf=args.input,
        metadata_path=args.meta,
        p1=args.p1,
        p2=args.p2,
        p3=args.p3,
        p4=args.p4,
        het_threshold=args.het_threshold,
        hom_threshold=args.hom_threshold,
        logger=logger,
    )

    agg = aggregate_by_gp(patterns)
    write_output_tsv(args.output, agg)

    logger.info("Done. Wrote %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
