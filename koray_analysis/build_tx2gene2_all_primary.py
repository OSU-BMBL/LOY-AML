#!/usr/bin/env python3
"""
build_tx2gene2_all_primary.py

Purpose
-------
Build the transcript-to-gene table used by tximport from the exact transcriptome
FASTA used to create the kallisto index, and automatically harmonize every
eligible human Ensembl alternative-locus gene record to its Ensembl 116 primary
assembly gene identifier. This is a genome-wide correction, not a PSMB8/PSMB9
whitelist.

Inputs
------
1. One or more Ensembl or GENCODE transcriptome FASTA files.
2. Either:
   a. the Ensembl 116 human core database, queried through a local mysql or
      mariadb command-line client, or
   b. a frozen mapping TSV produced by an earlier run.
3. Optional SYMBOL=ENSG validation checks. These checks do not restrict which
   genes are mapped.

Outputs
-------
1. A seven-column tx2gene TSV:
   transcript_id, gene_id, symbol, Chr, start, stop, strand.
2. A gene-level QC TSV documenting every source gene ID, final analysis gene ID,
   mapping status, and transcript count.
3. Optionally, a reusable Ensembl 116 all-gene alternative-to-primary map.

Author
------
Koray Dogan Kaya

Genome-wide mapping rule
------------------------
* Every Ensembl alt_allele group is inspected.
* PAR groups are retained separately to preserve biologically meaningful X/Y
  dosage information.
* A group is collapsed only when it contains exactly one gene on a canonical
  primary chromosome (1-22, X, Y, MT/M).
* If Ensembl supplies one IS_REPRESENTATIVE member, it must agree with that
  unique primary-chromosome member.
* Groups with no unique primary member, conflicting representatives, or other
  ambiguity remain unchanged and are reported in the mapping table.
* Repeated symbols outside an explicit Ensembl alt_allele group are never merged.
* Transcript IDs retain the exact FASTA versions used by kallisto.
* Gene IDs are unversioned by default, yielding stable primary human Ensembl gene
  identifiers. Pass --keep-gene-version to retain gene versions.
* Chr/start/stop/strand remain the source transcript coordinates for provenance.
  The final analysis gene identity is recorded in gene_id and the QC table.

Examples
--------
# First run. Query release-locked Ensembl 116 and save a reusable all-gene map.
python build_tx2gene2_all_primary.py \\
  --fasta Homo_sapiens.GRCh38.cdna.all.fa.gz \\
  --primary-map-out ensembl116_all_alt_to_primary.tsv \\
  --out tx2gene.tsv \\
  --qc-out tx2gene_primary_mapping_qc.tsv

# Reuse the frozen mapping without connecting to Ensembl.
python build_tx2gene2_all_primary.py \\
  --fasta Homo_sapiens.GRCh38.cdna.all.fa.gz \\
  --primary-map ensembl116_all_alt_to_primary.tsv \\
  --out tx2gene.tsv \\
  --qc-out tx2gene_primary_mapping_qc.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

VERSION_RE = re.compile(r"\.\d+$")
DEFAULT_CORE_DB = "homo_sapiens_core_116_38"
DEFAULT_PRIMARY_SEQNAMES = tuple([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])
VALID_MAPPED_STATUSES = {
    "mapped_to_unique_primary_representative",
    "mapped_to_unique_primary_member",
}


@dataclass(frozen=True)
class TranscriptRecord:
    transcript_id: str
    source_gene_id: str
    source_symbol: str
    chrom: str
    start: str
    stop: str
    strand: str
    location_type: str


@dataclass(frozen=True)
class PrimaryMapEntry:
    source_gene_id: str
    primary_gene_id: str
    primary_gene_id_versioned: str
    primary_symbol: str
    primary_seqname: str
    alt_allele_group_id: str
    status: str
    mapping_method: str


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def unversioned(identifier: str) -> str:
    return VERSION_RE.sub("", identifier or "")


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_key_value(value: str, option_name: str) -> Tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"{option_name} must have the form SYMBOL=ENSG..., received {value!r}"
        )
    key, target = (part.strip() for part in value.split("=", 1))
    if not key or not target:
        raise argparse.ArgumentTypeError(
            f"{option_name} requires non-empty values on both sides of '='"
        )
    return key, target


def parse_location(parts: Sequence[str]) -> Tuple[str, str, str, str, str]:
    """Return location type, sequence name, start, end, and strand."""
    for token in parts:
        fields = token.split(":")
        if (
            len(fields) == 6
            and fields[3].isdigit()
            and fields[4].isdigit()
            and fields[5] in {"1", "-1"}
        ):
            return (
                fields[0],
                fields[2],
                fields[3],
                fields[4],
                "+" if fields[5] == "1" else "-",
            )
    return "", "", "", "", ""


def parse_ensembl(header: str) -> TranscriptRecord:
    parts = header.split()
    transcript_id = parts[0]
    gene_id: Optional[str] = None
    symbol: Optional[str] = None
    for token in parts[1:]:
        if token.startswith("gene:"):
            gene_id = token[len("gene:") :]
        elif token.startswith("gene_symbol:"):
            symbol = token[len("gene_symbol:") :]
    location_type, chrom, start, stop, strand = parse_location(parts)
    gene_id = gene_id or transcript_id
    symbol = symbol or gene_id
    return TranscriptRecord(
        transcript_id=transcript_id,
        source_gene_id=gene_id,
        source_symbol=symbol,
        chrom=chrom,
        start=start,
        stop=stop,
        strand=strand,
        location_type=location_type,
    )


def parse_gencode(header: str) -> TranscriptRecord:
    fields = header.split("|")
    transcript_id = fields[0]
    gene_id = fields[1] if len(fields) > 1 and fields[1] else transcript_id
    symbol = fields[5] if len(fields) > 5 and fields[5] else gene_id
    return TranscriptRecord(
        transcript_id=transcript_id,
        source_gene_id=gene_id,
        source_symbol=symbol,
        chrom="",
        start="",
        stop="",
        strand="",
        location_type="",
    )


def read_fasta_records(fastas: Sequence[str]) -> Tuple[List[TranscriptRecord], Dict[str, int]]:
    records_by_transcript: MutableMapping[str, TranscriptRecord] = {}
    stats = {
        "headers": 0,
        "with_symbol": 0,
        "with_location": 0,
        "duplicate_identical": 0,
    }

    for fasta in fastas:
        if not os.path.isfile(fasta):
            raise FileNotFoundError(f"Transcriptome FASTA not found: {fasta}")
        with open_text(fasta) as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                stats["headers"] += 1
                header = line[1:].rstrip("\n\r")
                record = (
                    parse_gencode(header)
                    if "|" in header and "gene:" not in header
                    else parse_ensembl(header)
                )
                if record.source_symbol and record.source_symbol != record.source_gene_id:
                    stats["with_symbol"] += 1
                if record.chrom:
                    stats["with_location"] += 1

                previous = records_by_transcript.get(record.transcript_id)
                if previous is None:
                    records_by_transcript[record.transcript_id] = record
                elif previous == record:
                    stats["duplicate_identical"] += 1
                else:
                    raise ValueError(
                        "Conflicting FASTA headers use the same transcript ID "
                        f"{record.transcript_id!r}. Existing={previous}; new={record}"
                    )

    return list(records_by_transcript.values()), stats


def resolve_fastas(args: argparse.Namespace) -> List[str]:
    if args.fasta:
        return [str(path) for path in args.fasta]
    if not os.path.isfile(args.config):
        raise FileNotFoundError(
            f"No --fasta was supplied and configuration file was not found: {args.config}"
        )
    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required when FASTA paths are read from config.yaml. "
            "Install PyYAML or pass --fasta explicitly."
        ) from exc

    with open(args.config, "rt") as handle:
        config = yaml.safe_load(handle) or {}
    value = (config.get("reference", {}) or {}).get("transcriptome_fasta")
    if isinstance(value, str):
        return [value]
    if isinstance(value, list) and value:
        return [str(path) for path in value]
    raise KeyError(
        "config.yaml must contain reference.transcriptome_fasta as a path or list of paths"
    )


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def choose_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lookup = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        match = lookup.get(candidate.lower())
        if match:
            return match
    return None


ENSEMBL_ALT_ALLELE_SQL = r"""
SELECT
    aa.alt_allele_group_id,
    g.stable_id AS source_gene_id,
    g.version AS source_gene_version,
    COALESCE(x.display_label, '') AS source_symbol,
    sr.name AS source_seqname,
    cs.name AS coord_system,
    MAX(CASE WHEN aaa.attrib = 'IS_REPRESENTATIVE' THEN 1 ELSE 0 END) AS is_representative,
    MAX(CASE WHEN aaa.attrib = 'IS_PAR' THEN 1 ELSE 0 END) AS is_par
FROM alt_allele aa
JOIN gene g
  ON g.gene_id = aa.gene_id
JOIN seq_region sr
  ON sr.seq_region_id = g.seq_region_id
JOIN coord_system cs
  ON cs.coord_system_id = sr.coord_system_id
LEFT JOIN xref x
  ON x.xref_id = g.display_xref_id
LEFT JOIN alt_allele_attrib aaa
  ON aaa.alt_allele_id = aa.alt_allele_id
GROUP BY
    aa.alt_allele_group_id,
    g.gene_id,
    g.stable_id,
    g.version,
    x.display_label,
    sr.name,
    cs.name
ORDER BY
    aa.alt_allele_group_id,
    g.stable_id
""".strip()


def versioned_gene_id(stable_id: str, version: str) -> str:
    if stable_id and version not in {"", "0", "NULL", "None"}:
        return f"{stable_id}.{version}"
    return stable_id


def build_mapping_from_members(
    members: Iterable[Mapping[str, str]],
    primary_seqnames: Set[str],
) -> Tuple[Dict[str, PrimaryMapEntry], List[Dict[str, str]], Dict[str, int]]:
    groups: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for member in members:
        row = {str(key): str(value or "") for key, value in member.items()}
        groups[row["alt_allele_group_id"]].append(row)

    mapping: Dict[str, PrimaryMapEntry] = {}
    output_rows: List[Dict[str, str]] = []
    summary = {
        "groups_total": len(groups),
        "groups_mapped": 0,
        "groups_PAR_retained": 0,
        "groups_no_primary_member": 0,
        "groups_multiple_primary_members": 0,
        "groups_representative_conflict": 0,
    }

    for group_id, group_members in groups.items():
        par_group = any(truthy(member.get("is_par", "")) for member in group_members)
        primary_members = [
            member
            for member in group_members
            if member.get("coord_system") == "chromosome"
            and member.get("source_seqname") in primary_seqnames
        ]
        representatives = [
            member
            for member in group_members
            if truthy(member.get("is_representative", ""))
        ]

        target: Optional[Dict[str, str]] = None
        if par_group:
            status = "excluded_PAR"
            summary["groups_PAR_retained"] += 1
        elif len(primary_members) == 0:
            status = "no_primary_member"
            summary["groups_no_primary_member"] += 1
        elif len(primary_members) > 1:
            status = "multiple_primary_members"
            summary["groups_multiple_primary_members"] += 1
        else:
            candidate = primary_members[0]
            if len(representatives) > 1:
                status = "multiple_representatives"
                summary["groups_representative_conflict"] += 1
            elif len(representatives) == 1 and (
                representatives[0]["source_gene_id"] != candidate["source_gene_id"]
            ):
                status = "representative_primary_conflict"
                summary["groups_representative_conflict"] += 1
            else:
                target = candidate
                status = (
                    "mapped_to_unique_primary_representative"
                    if len(representatives) == 1
                    else "mapped_to_unique_primary_member"
                )
                summary["groups_mapped"] += 1

        primary_gene_id = target["source_gene_id"] if target else ""
        primary_gene_versioned = (
            versioned_gene_id(target["source_gene_id"], target["source_gene_version"])
            if target
            else ""
        )
        primary_symbol = target.get("source_symbol", "") if target else ""
        primary_seqname = target.get("source_seqname", "") if target else ""

        for member in group_members:
            source_gene_id = member["source_gene_id"]
            source_gene_versioned = versioned_gene_id(
                source_gene_id, member.get("source_gene_version", "")
            )
            row = {
                "alt_allele_group_id": group_id,
                "source_gene_id": source_gene_id,
                "source_gene_id_versioned": source_gene_versioned,
                "source_symbol": member.get("source_symbol", ""),
                "source_seqname": member.get("source_seqname", ""),
                "coord_system": member.get("coord_system", ""),
                "is_representative": member.get("is_representative", ""),
                "is_par": member.get("is_par", ""),
                "primary_gene_id": primary_gene_id,
                "primary_gene_id_versioned": primary_gene_versioned,
                "primary_symbol": primary_symbol,
                "primary_seqname": primary_seqname,
                "status": status,
                "mapping_method": "Ensembl_alt_allele_group",
            }
            output_rows.append(row)
            if target:
                mapping[source_gene_id] = PrimaryMapEntry(
                    source_gene_id=source_gene_id,
                    primary_gene_id=primary_gene_id,
                    primary_gene_id_versioned=primary_gene_versioned,
                    primary_symbol=primary_symbol,
                    primary_seqname=primary_seqname,
                    alt_allele_group_id=group_id,
                    status=status,
                    mapping_method="Ensembl_alt_allele_group",
                )

    return mapping, output_rows, summary


def query_ensembl_mapping(
    database: str,
    host: str,
    port: int,
    user: str,
    mysql_client: Optional[str],
    primary_seqnames: Set[str],
) -> Tuple[Dict[str, PrimaryMapEntry], List[Dict[str, str]], Dict[str, int]]:
    client = mysql_client or shutil.which("mariadb") or shutil.which("mysql")
    if not client:
        raise RuntimeError(
            "No mysql or mariadb command-line client was found. Install one, "
            "for example 'sudo apt install default-mysql-client', or use --primary-map."
        )

    command = [
        client,
        "--batch",
        "--raw",
        "--host",
        host,
        "--port",
        str(port),
        "--user",
        user,
        "--execute",
        ENSEMBL_ALT_ALLELE_SQL,
        database,
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "Ensembl MySQL query failed. Server response:\n" + completed.stderr.strip()
        )

    reader = csv.DictReader(completed.stdout.splitlines(), delimiter="\t")
    required = {
        "alt_allele_group_id",
        "source_gene_id",
        "source_gene_version",
        "source_symbol",
        "source_seqname",
        "coord_system",
        "is_representative",
        "is_par",
    }
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise RuntimeError(
            f"Unexpected Ensembl query columns: {reader.fieldnames!r}"
        )
    rows = [{str(k): str(v or "") for k, v in row.items()} for row in reader]
    if not rows:
        raise RuntimeError("The Ensembl alt_allele query returned no rows.")
    return build_mapping_from_members(rows, primary_seqnames)


def load_cached_mapping(path: str) -> Tuple[Dict[str, PrimaryMapEntry], List[Dict[str, str]]]:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Primary mapping TSV not found: {path}")

    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Mapping TSV has no header: {path}")

        source_col = choose_column(reader.fieldnames, ["source_gene_id", "alt_gene_id"])
        target_col = choose_column(
            reader.fieldnames, ["primary_gene_id", "representative_gene_id"]
        )
        if not source_col or not target_col:
            raise ValueError(
                "Mapping TSV must contain source_gene_id and primary_gene_id columns."
            )
        status_col = choose_column(reader.fieldnames, ["status", "mapping_status"])
        target_versioned_col = choose_column(
            reader.fieldnames,
            ["primary_gene_id_versioned", "representative_gene_id_versioned"],
        )
        target_symbol_col = choose_column(
            reader.fieldnames, ["primary_symbol", "representative_symbol"]
        )
        target_seqname_col = choose_column(
            reader.fieldnames, ["primary_seqname", "representative_seqname"]
        )
        group_col = choose_column(reader.fieldnames, ["alt_allele_group_id", "group_id"])
        method_col = choose_column(reader.fieldnames, ["mapping_method", "method"])

        mapping: Dict[str, PrimaryMapEntry] = {}
        all_rows: List[Dict[str, str]] = []
        for line_number, raw_row in enumerate(reader, start=2):
            row = {str(key): str(value or "") for key, value in raw_row.items()}
            all_rows.append(row)
            source = unversioned(row.get(source_col, "").strip())
            target = unversioned(row.get(target_col, "").strip())
            status = row.get(status_col, "").strip() if status_col else ""
            if not source or not target:
                continue
            if status and status not in VALID_MAPPED_STATUSES:
                continue

            entry = PrimaryMapEntry(
                source_gene_id=source,
                primary_gene_id=target,
                primary_gene_id_versioned=(
                    row.get(target_versioned_col, "").strip()
                    if target_versioned_col
                    else target
                ),
                primary_symbol=(
                    row.get(target_symbol_col, "").strip() if target_symbol_col else ""
                ),
                primary_seqname=(
                    row.get(target_seqname_col, "").strip()
                    if target_seqname_col
                    else ""
                ),
                alt_allele_group_id=(
                    row.get(group_col, "").strip() if group_col else ""
                ),
                status=status or "mapped_to_unique_primary_member",
                mapping_method=(
                    row.get(method_col, "").strip()
                    if method_col
                    else "cached_primary_map"
                ),
            )
            previous = mapping.get(source)
            if previous and previous.primary_gene_id != entry.primary_gene_id:
                raise ValueError(
                    f"Conflicting primary targets for {source} at line {line_number}: "
                    f"{previous.primary_gene_id} and {entry.primary_gene_id}"
                )
            mapping[source] = entry

    if not mapping:
        raise RuntimeError(f"No usable mapped rows were found in {path}")
    return mapping, all_rows


def write_tsv(path: str, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def validate_expectations(
    expectations: Sequence[str], final_rows: Sequence[Mapping[str, str]]
) -> None:
    for raw in expectations:
        symbol, expected_raw = parse_key_value(raw, "--expect-primary")
        expected = unversioned(expected_raw)
        observed = {
            unversioned(row["gene_id"])
            for row in final_rows
            if row["symbol"].upper() == symbol.upper()
            or row["source_symbol"].upper() == symbol.upper()
        }
        if observed != {expected}:
            raise RuntimeError(
                f"Primary-ID validation failed for {symbol}: expected only {expected}; "
                f"observed {sorted(observed) if observed else 'no matching records'}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Build tx2gene with genome-wide Ensembl 116 alternative-locus "
            "harmonization to primary human gene identifiers."
        ),
    )
    parser.add_argument(
        "--fasta",
        nargs="*",
        help="Transcriptome FASTA(s); default: config reference.transcriptome_fasta",
    )
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--out", default="tx2gene.tsv")
    parser.add_argument(
        "--qc-out",
        help="Mapping QC TSV; default: <out stem>.primary_mapping_qc.tsv",
    )
    parser.add_argument(
        "--primary-map",
        help="Frozen Ensembl all-gene alternative-to-primary mapping TSV.",
    )
    parser.add_argument(
        "--primary-map-out",
        help="Write the complete fetched Ensembl mapping table to this TSV.",
    )
    parser.add_argument(
        "--ensembl-core-db",
        default=DEFAULT_CORE_DB,
        help=f"Release-locked Ensembl human core database. Default: {DEFAULT_CORE_DB}",
    )
    parser.add_argument("--ensembl-host", default="ensembldb.ensembl.org")
    parser.add_argument("--ensembl-port", type=int, default=5306)
    parser.add_argument("--ensembl-user", default="anonymous")
    parser.add_argument(
        "--mysql-client",
        help="Path to mysql/mariadb; automatically detected when omitted.",
    )
    parser.add_argument(
        "--primary-seqnames",
        default=",".join(DEFAULT_PRIMARY_SEQNAMES),
        help="Canonical primary chromosome names used to select the unique target.",
    )
    parser.add_argument(
        "--keep-gene-version",
        action="store_true",
        help="Retain .N suffixes on gene IDs. Gene IDs are unversioned by default.",
    )
    parser.add_argument(
        "--strip-transcript-version",
        action="store_true",
        help=(
            "Remove transcript versions. Do this only when the kallisto index target "
            "IDs are also unversioned."
        ),
    )
    parser.add_argument(
        "--expect-primary",
        action="append",
        default=[],
        metavar="SYMBOL=ENSG",
        help="Optional validation only; repeatable and does not restrict mapping.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fastas = resolve_fastas(args)
    records, fasta_stats = read_fasta_records(fastas)
    if not records:
        raise RuntimeError("No transcript headers were found in the supplied FASTA file(s).")

    primary_seqnames = {
        value.strip() for value in args.primary_seqnames.split(",") if value.strip()
    }
    if args.primary_map:
        mapping, mapping_rows = load_cached_mapping(args.primary_map)
        mapping_summary: Dict[str, int] = {}
        mapping_source = f"cached TSV: {args.primary_map}"
    else:
        mapping, mapping_rows, mapping_summary = query_ensembl_mapping(
            database=args.ensembl_core_db,
            host=args.ensembl_host,
            port=args.ensembl_port,
            user=args.ensembl_user,
            mysql_client=args.mysql_client,
            primary_seqnames=primary_seqnames,
        )
        mapping_source = f"Ensembl: {args.ensembl_core_db}"

    source_full_ids_by_unversioned: Dict[str, Set[str]] = defaultdict(set)
    symbols_by_gene: Dict[str, Set[str]] = defaultdict(set)
    for record in records:
        key = unversioned(record.source_gene_id)
        source_full_ids_by_unversioned[key].add(record.source_gene_id)
        if record.source_symbol:
            symbols_by_gene[key].add(record.source_symbol)

    final_rows: List[Dict[str, str]] = []
    output_transcripts: Dict[str, Dict[str, str]] = {}
    mapped_source_genes: Set[str] = set()
    mapped_primary_genes: Set[str] = set()
    mapped_transcripts = 0

    for record in records:
        source_key = unversioned(record.source_gene_id)
        entry = mapping.get(source_key)
        if entry:
            if args.keep_gene_version:
                target_id = entry.primary_gene_id_versioned or entry.primary_gene_id
            else:
                target_id = entry.primary_gene_id
            target_symbol = entry.primary_symbol
            if not target_symbol:
                candidate_symbols = symbols_by_gene.get(entry.primary_gene_id, set())
                target_symbol = (
                    next(iter(candidate_symbols))
                    if len(candidate_symbols) == 1
                    else record.source_symbol
                )
            mapping_status = (
                f"primary:{entry.status}"
                if source_key == entry.primary_gene_id
                else f"mapped:{entry.status}"
            )
            if source_key != entry.primary_gene_id:
                mapped_transcripts += 1
                mapped_source_genes.add(source_key)
                mapped_primary_genes.add(entry.primary_gene_id)
            analysis_gene_id = target_id
            analysis_symbol = target_symbol
            group_id = entry.alt_allele_group_id
            mapping_method = entry.mapping_method
        else:
            analysis_gene_id = (
                record.source_gene_id
                if args.keep_gene_version
                else unversioned(record.source_gene_id)
            )
            analysis_symbol = record.source_symbol
            mapping_status = "unchanged"
            group_id = ""
            mapping_method = ""

        transcript_id = (
            unversioned(record.transcript_id)
            if args.strip_transcript_version
            else record.transcript_id
        )
        row = {
            "transcript_id": transcript_id,
            "gene_id": analysis_gene_id,
            "symbol": analysis_symbol,
            "Chr": record.chrom,
            "start": record.start,
            "stop": record.stop,
            "strand": record.strand,
            "source_gene_id": record.source_gene_id,
            "source_symbol": record.source_symbol,
            "mapping_status": mapping_status,
            "mapping_method": mapping_method,
            "alt_allele_group_id": group_id,
        }

        previous = output_transcripts.get(transcript_id)
        if previous is None:
            output_transcripts[transcript_id] = row
            final_rows.append(row)
        elif previous != row:
            raise RuntimeError(
                "Transcript version stripping created a conflicting transcript ID: "
                f"{transcript_id}. Preserve transcript versions or correct the FASTA inputs."
            )

    validate_expectations(args.expect_primary, final_rows)

    output_path = args.out
    qc_path = args.qc_out or str(Path(output_path).with_suffix("")) + ".primary_mapping_qc.tsv"
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "wt", newline="") as handle:
        fields = ["transcript_id", "gene_id", "symbol", "Chr", "start", "stop", "strand"]
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(final_rows)

    qc_groups: Dict[Tuple[str, str, str, str, str], Dict[str, object]] = {}
    for row in final_rows:
        key = (
            row["source_gene_id"],
            row["gene_id"],
            row["source_symbol"],
            row["symbol"],
            row["mapping_status"],
        )
        item = qc_groups.setdefault(
            key,
            {
                "source_gene_id": row["source_gene_id"],
                "source_gene_id_unversioned": unversioned(row["source_gene_id"]),
                "source_symbol": row["source_symbol"],
                "source_Chr": set(),
                "analysis_gene_id": row["gene_id"],
                "analysis_gene_id_unversioned": unversioned(row["gene_id"]),
                "analysis_symbol": row["symbol"],
                "mapping_status": row["mapping_status"],
                "mapping_method": row["mapping_method"],
                "alt_allele_group_id": row["alt_allele_group_id"],
                "n_transcripts": 0,
            },
        )
        if row["Chr"]:
            item["source_Chr"].add(row["Chr"])  # type: ignore[union-attr]
        item["n_transcripts"] = int(item["n_transcripts"]) + 1

    qc_rows: List[Dict[str, object]] = []
    for item in qc_groups.values():
        qc_rows.append(
            {
                **{key: value for key, value in item.items() if key != "source_Chr"},
                "source_Chr": ";".join(sorted(item["source_Chr"])),  # type: ignore[arg-type]
            }
        )
    qc_rows.sort(
        key=lambda row: (str(row["analysis_symbol"]), str(row["source_gene_id"]))
    )
    write_tsv(qc_path, qc_rows)

    if args.primary_map_out:
        write_tsv(args.primary_map_out, mapping_rows)

    eprint(
        f"[build_tx2gene] {fasta_stats['headers']} FASTA headers -> "
        f"{len(final_rows)} exact transcript targets"
    )
    eprint(f"[build_tx2gene] mapping source: {mapping_source}")
    eprint(
        f"[build_tx2gene] reassigned {mapped_transcripts} transcripts from "
        f"{len(mapped_source_genes)} alternative gene IDs to "
        f"{len(mapped_primary_genes)} primary gene IDs"
    )
    eprint(f"[build_tx2gene] wrote {output_path}")
    eprint(f"[build_tx2gene] wrote {qc_path}")
    if args.primary_map_out:
        eprint(f"[build_tx2gene] wrote {args.primary_map_out}")
    for key, value in mapping_summary.items():
        eprint(f"[build_tx2gene] {key}: {value}")
    if fasta_stats["with_location"] == 0:
        eprint(
            "[build_tx2gene] NOTE: FASTA headers had no source coordinates. "
            "Mapping remains valid because the target relationships come from Ensembl."
        )
    if fasta_stats["duplicate_identical"]:
        eprint(
            f"[build_tx2gene] NOTE: ignored {fasta_stats['duplicate_identical']} "
            "identical duplicate transcript headers"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        eprint(f"[build_tx2gene] ERROR: {exc}")
        sys.exit(1)
