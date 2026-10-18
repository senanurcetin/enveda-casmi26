"""Safely render NEXUS read-only CASMI26 BigQuery SQL templates.

No credentials or datasets are accessed here. Execute via bq CLI after
inspecting the rendered SQL and a BigQuery dry run.

  python scripts/render_bigquery_sql.py \
    --project my-cloud-project --dataset casmi26 --table train_spectra \
    --output runs/library_quality.sql

  bq query --use_legacy_sql=false --dry_run < runs/library_quality.sql
  bq query --use_legacy_sql=false --format=prettyjson < runs/library_quality.sql
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "sql" / "bigquery" / "library_quality.sql"
PROJECT_PATTERN = re.compile(r"[a-z][a-z0-9-]{4,61}[a-z0-9]\Z")
IDENTIFIER_PATTERN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")


def render(template: str, project: str, dataset: str, table: str) -> str:
    if not PROJECT_PATTERN.fullmatch(project):
        raise ValueError("Unsafe or invalid GCP project ID")
    if not IDENTIFIER_PATTERN.fullmatch(dataset):
        raise ValueError("Unsafe or invalid BigQuery dataset ID")
    if not IDENTIFIER_PATTERN.fullmatch(table):
        raise ValueError("Unsafe or invalid BigQuery table ID")
    rendered = (template.replace("{{PROJECT_ID}}", project)
                .replace("{{DATASET_ID}}", dataset)
                .replace("{{TABLE_NAME}}", table))
    if "{{" in rendered or "}}" in rendered:
        raise ValueError("Unresolved SQL template tokens")
    return rendered


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--table", default="train_spectra")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    sql = render(TEMPLATE.read_text(encoding="utf-8"), args.project,
                 args.dataset, args.table)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(sql, encoding="utf-8")
        print(f"Rendered {args.output}")
    else:
        print(sql)


if __name__ == "__main__":
    main()
