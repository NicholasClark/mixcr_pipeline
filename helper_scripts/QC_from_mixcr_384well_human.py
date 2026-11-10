#!/usr/bin/env python3
"""
Build a 384-well QC workbook directly from MiXCR output files.

Input:
  - <mixcr_output_dir> containing:
      *.align.report.txt
      *.TCRa_UD##.assemble.report.txt
      *.TCRb_UD##.assemble.report.txt

Output:
  - <output_prefix>_TCR_QC_384_format.xlsx
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.formatting.rule import ColorScale, FormatObject, Rule
from openpyxl.utils import get_column_letter

ROW_NAMES = list("ABCDEFGHIJKLMNOP")
COL_NAMES = list(range(1, 25))
WELL_ORDER = [f"{r}{c}" for r in ROW_NAMES for c in COL_NAMES]
ID_LIST = [f"UD{i:02d}" for i in range(1, 33)]

AB_TCR_METRICS = [
    "TCRa_rawReads",
    "TCRa_rawReads_frac",
    "TCRb_rawReads",
    "TCRb_rawReads_frac",
    "TCRa_cloneReads",
    "TCRa_cloneReads_frac",
    "TCRb_cloneReads",
    "TCRb_cloneReads_frac",
    "TCRa_clones",
    "TCRb_clones",
]

PARAS1 = [
    "Sample",
    "Total_sequencing_reads",
    "Successfully_aligned_reads",
    "Successfully_aligned_reads_frac",
    "TRA_chains",
    "TRA_chains_frac",
    "TRB_chains",
    "TRB_chains_frac",
]

PARAS2 = [f"{ud}_{metric}" for ud in ID_LIST for metric in AB_TCR_METRICS]
ALL_METRICS = PARAS1 + PARAS2

ALIGN_PATTERNS = {
    "Total_sequencing_reads": re.compile(r"^Total sequencing reads:\s+(\d+)"),
    "Successfully_aligned_reads": re.compile(r"^Successfully aligned reads:\s+(\d+)\s+\(([0-9.]+)%\)"),
    "TRA_chains": re.compile(r"^TRA chains:\s+(\d+)\s+\(([0-9.]+)%\)"),
    "TRB_chains": re.compile(r"^TRB chains:\s+(\d+)\s+\(([0-9.]+)%\)"),
    "ud_tcra": re.compile(r"^\s*TCRa_UD(\d{2}):\s+(\d+)\s+\(([0-9.]+)%\)"),
    "ud_tcrb": re.compile(r"^\s*TCRb_UD(\d{2}):\s+(\d+)\s+\(([0-9.]+)%\)"),
}

ASSEMBLE_READS_RE = re.compile(r"^Reads used in clonotypes, percent of total:\s+(\d+)\s+\(([0-9.]+)%\)")
ASSEMBLE_CLONES_RE = re.compile(r"^Final clonotype count:\s+(\d+)")
WELL_RE = re.compile(r"_([A-P](?:[1-9]|1[0-9]|2[0-4]))_")


def parse_align_report(path: Path) -> dict:
    parsed = {
        "Total_sequencing_reads": 0.0,
        "Successfully_aligned_reads": 0.0,
        "Successfully_aligned_reads_frac": 0.0,
        "TRA_chains": 0.0,
        "TRA_chains_frac": 0.0,
        "TRB_chains": 0.0,
        "TRB_chains_frac": 0.0,
        "ud_tcra_raw": {},
        "ud_tcrb_raw": {},
    }

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")

            m = ALIGN_PATTERNS["Total_sequencing_reads"].match(line)
            if m:
                parsed["Total_sequencing_reads"] = float(m.group(1))
                continue

            m = ALIGN_PATTERNS["Successfully_aligned_reads"].match(line)
            if m:
                parsed["Successfully_aligned_reads"] = float(m.group(1))
                parsed["Successfully_aligned_reads_frac"] = float(m.group(2)) / 100.0
                continue

            m = ALIGN_PATTERNS["TRA_chains"].match(line)
            if m:
                parsed["TRA_chains"] = float(m.group(1))
                parsed["TRA_chains_frac"] = float(m.group(2)) / 100.0
                continue

            m = ALIGN_PATTERNS["TRB_chains"].match(line)
            if m:
                parsed["TRB_chains"] = float(m.group(1))
                parsed["TRB_chains_frac"] = float(m.group(2)) / 100.0
                continue

            m = ALIGN_PATTERNS["ud_tcra"].match(line)
            if m:
                ud = f"UD{m.group(1)}"
                parsed["ud_tcra_raw"][ud] = (float(m.group(2)), float(m.group(3)) / 100.0)
                continue

            m = ALIGN_PATTERNS["ud_tcrb"].match(line)
            if m:
                ud = f"UD{m.group(1)}"
                parsed["ud_tcrb_raw"][ud] = (float(m.group(2)), float(m.group(3)) / 100.0)
                continue

    return parsed


def parse_assemble_report(path: Path) -> tuple[float, float, float]:
    clone_reads = 0.0
    clone_reads_frac = 0.0
    clones = 0.0

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = ASSEMBLE_READS_RE.match(line)
            if m:
                clone_reads = float(m.group(1))
                clone_reads_frac = float(m.group(2)) / 100.0
                continue

            m = ASSEMBLE_CLONES_RE.match(line)
            if m:
                clones = float(m.group(1))

    return clone_reads, clone_reads_frac, clones


def apply_three_color_formatting(
    writer: pd.ExcelWriter,
    sheet_name: str,
    min_val: float,
    mid_val: float,
    max_val: float,
    number_format: str,
    startrow: int = 0,
    startcol: int = 0,
) -> None:
    worksheet = writer.sheets[sheet_name]

    # Skip the index column, color only the 24 plate-value columns.
    first_col = startcol + 2
    last_col = startcol + 25
    first_row = startrow + 2
    last_row = startrow + 17

    range_str = (
        f"{get_column_letter(first_col)}{first_row}:"
        f"{get_column_letter(last_col)}{last_row}"
    )

    color_scale = ColorScale(
        cfvo=[
            FormatObject(type="num", val=float(min_val)),
            FormatObject(type="num", val=float(mid_val)),
            FormatObject(type="num", val=float(max_val)),
        ],
        color=["FF9494", "FFFFFF", "6AB4FF"],
    )
    worksheet.conditional_formatting.add(
        range_str, Rule(type="colorScale", colorScale=color_scale)
    )

    for row in worksheet.iter_rows(
        min_row=first_row,
        max_row=last_row,
        min_col=first_col,
        max_col=last_col,
    ):
        for cell in row:
            cell.number_format = number_format


def plate_matrix(series: pd.Series) -> pd.DataFrame:
    return pd.DataFrame(
        [[series[f"{r}{c}"] for c in COL_NAMES] for r in ROW_NAMES],
        index=ROW_NAMES,
        columns=COL_NAMES,
    )


def build_qc_dataframe(mixcr_dir: Path) -> pd.DataFrame:
    rows = []
    align_reports = sorted(mixcr_dir.glob("*.align.report.txt"))

    for align_report in align_reports:
        sample = align_report.name.replace(".align.report.txt", "")
        m_well = WELL_RE.search(sample)
        if not m_well:
            continue

        well = m_well.group(1)
        row = {metric: 0.0 for metric in ALL_METRICS if metric != "Sample"}
        row["Sample"] = sample
        row["Well"] = well

        align = parse_align_report(align_report)
        row["Total_sequencing_reads"] = align["Total_sequencing_reads"]
        row["Successfully_aligned_reads"] = align["Successfully_aligned_reads"]
        row["Successfully_aligned_reads_frac"] = align["Successfully_aligned_reads_frac"]
        row["TRA_chains"] = align["TRA_chains"]
        row["TRA_chains_frac"] = align["TRA_chains_frac"]
        row["TRB_chains"] = align["TRB_chains"]
        row["TRB_chains_frac"] = align["TRB_chains_frac"]

        for ud in ID_LIST:
            tcra_report = mixcr_dir / f"{sample}.TCRa_{ud}.assemble.report.txt"
            tcrb_report = mixcr_dir / f"{sample}.TCRb_{ud}.assemble.report.txt"
            if not tcra_report.exists() or not tcrb_report.exists():
                continue

            tcra_raw, tcra_raw_frac = align["ud_tcra_raw"].get(ud, (0.0, 0.0))
            tcrb_raw, tcrb_raw_frac = align["ud_tcrb_raw"].get(ud, (0.0, 0.0))
            row[f"{ud}_TCRa_rawReads"] = tcra_raw
            row[f"{ud}_TCRa_rawReads_frac"] = tcra_raw_frac
            row[f"{ud}_TCRb_rawReads"] = tcrb_raw
            row[f"{ud}_TCRb_rawReads_frac"] = tcrb_raw_frac

            tcra_clone_reads, tcra_clone_frac, tcra_clones = parse_assemble_report(tcra_report)
            tcrb_clone_reads, tcrb_clone_frac, tcrb_clones = parse_assemble_report(tcrb_report)

            row[f"{ud}_TCRa_cloneReads"] = tcra_clone_reads
            row[f"{ud}_TCRa_cloneReads_frac"] = tcra_clone_frac
            row[f"{ud}_TCRb_cloneReads"] = tcrb_clone_reads
            row[f"{ud}_TCRb_cloneReads_frac"] = tcrb_clone_frac
            row[f"{ud}_TCRa_clones"] = tcra_clones
            row[f"{ud}_TCRb_clones"] = tcrb_clones

        rows.append(row)

    if not rows:
        raise RuntimeError(
            f"No 384-well sample align reports found in {mixcr_dir} "
            "(expected names like *_A1_* through *_P24_*)"
        )

    qc_df = pd.DataFrame(rows)

    existing_wells = set(qc_df["Well"].tolist())
    missing_wells = [w for w in WELL_ORDER if w not in existing_wells]
    for well in missing_wells:
        empty_row = {metric: 0.0 for metric in ALL_METRICS if metric != "Sample"}
        empty_row["Sample"] = well
        empty_row["Well"] = well
        qc_df = pd.concat([qc_df, pd.DataFrame([empty_row])], ignore_index=True)

    qc_df["Well"] = pd.Categorical(qc_df["Well"], categories=WELL_ORDER, ordered=True)
    qc_df = qc_df.sort_values("Well")
    qc_df = qc_df.set_index("Well")
    if not qc_df.index.is_unique:
        raise RuntimeError("Duplicate wells detected after QC table assembly.")

    return qc_df


def write_qc_workbook(qc_df: pd.DataFrame, output_prefix: Path) -> Path:
    summary_rows = []
    for ud in ID_LIST:
        for metric in AB_TCR_METRICS:
            col = f"{ud}_{metric}"
            mean_val = qc_df[col].mean()
            failed = (qc_df[col] < 0.5 * mean_val).sum()
            summary_rows.append(
                {
                    "UD": ud,
                    "Parameter": metric,
                    "Average": mean_val,
                    "QC_failed_wells": failed,
                }
            )
    summary_df = pd.DataFrame(summary_rows)

    id_rank = (
        qc_df.filter(regex=r"_TCRb_clones$", axis=1)
        .sum()
        .sort_values(ascending=False)
        .index.str.replace("_TCRb_clones", "", regex=False)
        .tolist()
    )

    out_file = Path(f"{output_prefix}_TCR_QC_384_format.xlsx")
    with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
        for idx, stat in enumerate(["Average", "QC_failed_wells"]):
            pivot = summary_df.pivot(index="UD", columns="Parameter", values=stat)
            pivot = pivot.reindex(columns=AB_TCR_METRICS)
            pivot.to_excel(
                writer,
                sheet_name="Summary",
                index_label=stat,
                startrow=0,
                startcol=0 if idx == 0 else 14,
            )

        for metric in [
            "Sample",
            "Total_sequencing_reads",
            "Successfully_aligned_reads",
            "Successfully_aligned_reads_frac",
        ]:
            df_plate = plate_matrix(qc_df[metric])
            df_plate.to_excel(writer, sheet_name=metric)
            if metric != "Sample":
                number_format = "0.0%" if metric.endswith("_frac") else "#,##0"
                col_values = qc_df[metric].values.astype(float)
                apply_three_color_formatting(
                    writer,
                    sheet_name=metric,
                    min_val=np.min(col_values),
                    mid_val=np.percentile(col_values, 50),
                    max_val=np.max(col_values),
                    number_format=number_format,
                )

        chain_layout = {
            "TRA_chains": (0, 0),
            "TRA_chains_frac": (0, 30),
            "TRB_chains": (22, 0),
            "TRB_chains_frac": (22, 30),
        }
        for metric, (startrow, startcol) in chain_layout.items():
            df_plate = plate_matrix(qc_df[metric])
            df_plate.to_excel(
                writer,
                sheet_name="TRA_TRB_chains",
                index_label=metric,
                startrow=startrow,
                startcol=startcol,
            )
            number_format = "0.0%" if metric.endswith("_frac") else "#,##0"
            col_values = qc_df[metric].values.astype(float)
            apply_three_color_formatting(
                writer,
                sheet_name="TRA_TRB_chains",
                min_val=np.min(col_values),
                mid_val=np.percentile(col_values, 50),
                max_val=np.max(col_values),
                number_format=number_format,
                startrow=startrow,
                startcol=startcol,
            )

        metric_layout = {
            "TCRa_rawReads": (0, 0),
            "TCRa_rawReads_frac": (0, 30),
            "TCRb_rawReads": (22, 0),
            "TCRb_rawReads_frac": (22, 30),
            "TCRa_cloneReads": (44, 0),
            "TCRa_cloneReads_frac": (44, 30),
            "TCRb_cloneReads": (66, 0),
            "TCRb_cloneReads_frac": (66, 30),
            "TCRa_clones": (88, 0),
            "TCRb_clones": (110, 0),
        }

        for ud in id_rank:
            sheet_name = f"abTCR_{ud}"
            for metric, (startrow, startcol) in metric_layout.items():
                col = f"{ud}_{metric}"
                df_plate = plate_matrix(qc_df[col])
                df_plate.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index_label=col,
                    startrow=startrow,
                    startcol=startcol,
                )
                number_format = "0.0%" if metric.endswith("_frac") else "#,##0"
                col_values = qc_df[col].values.astype(float)
                apply_three_color_formatting(
                    writer,
                    sheet_name=sheet_name,
                    min_val=np.min(col_values),
                    mid_val=np.percentile(col_values, 50),
                    max_val=np.max(col_values),
                    number_format=number_format,
                    startrow=startrow,
                    startcol=startcol,
                )

    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate 384-well QC workbook directly from MiXCR output files."
    )
    parser.add_argument("mixcr_output_dir", help="Directory with MiXCR output reports")
    parser.add_argument(
        "output_prefix",
        help="Output prefix (produces <output_prefix>_TCR_QC_384_format.xlsx)",
    )
    args = parser.parse_args()

    mixcr_dir = Path(args.mixcr_output_dir).resolve()
    if not mixcr_dir.exists():
        raise FileNotFoundError(f"MiXCR output directory not found: {mixcr_dir}")

    qc_df = build_qc_dataframe(mixcr_dir)
    out_file = write_qc_workbook(qc_df, Path(args.output_prefix))

    print(f"Input samples in QC table: {qc_df.shape[0]}")
    print(f"Output workbook: {out_file}")


if __name__ == "__main__":
    main()
