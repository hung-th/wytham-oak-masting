"""
Upload large GWAS files to Zenodo.

Files uploaded (all gitignored due to >100 MB):
  gwas/lfmm2_res_WWO_37WGS_pruned_imputed.rds   (~152 MB)
  gwas/WWO_37WGS_pruned_snmf_sub.snmf/           (~160 MB, zipped before upload)
  gwas/output/n37/lfmm.z.tsv                    (~179 MB)
  gwas/output/n37/lfmm.p.tsv                    (~178 MB)
  gwas/output/n37/lfmm.q.tsv                    (~176 MB)

Usage
-----
1. Get a Zenodo access token:
     https://zenodo.org/account/settings/applications/tokens/new/
   Tick scope: deposit:actions, deposit:write
   For sandbox testing first:
     https://sandbox.zenodo.org/account/settings/applications/tokens/new/

2. Install dependency:
     pip install requests

3. Run (sandbox first):
     python scripts/zenodo_upload.py --token YOUR_TOKEN --sandbox

4. Once happy, publish for real:
     python scripts/zenodo_upload.py --token YOUR_TOKEN

   The script stops before publishing — review the record URL printed,
   then re-run with --publish to mint the DOI.

5. Copy the DOI into README.md and update .gitignore notes.
"""

import argparse
import json
import os
import sys
import zipfile
from pathlib import Path

import requests

# ── Metadata ─────────────────────────────────────────────────────────────────
METADATA = {
    "title": (
        "Wytham Woods oak masting GWAS — large genomic data files"
    ),
    "description": (
        "<p>Large genomic data files from a genome-wide association study (GWAS) "
        "of acorn production variation in <em>Quercus robur</em> at Wytham Woods, "
        "Oxford. These files are too large to host on GitHub and accompany the "
        "repository at https://github.com/hung-th/wytham-oak-masting.</p>"
        "<p>Files included:</p><ul>"
        "<li>lfmm2_res_WWO_37WGS_pruned_imputed.rds — fitted LFMM2 model object</li>"
        "<li>WWO_37WGS_pruned_snmf_sub.snmf.zip — sNMF project directory</li>"
        "<li>lfmm.z.tsv — LFMM2 Z-scores (all SNPs × phenotypes)</li>"
        "<li>lfmm.p.tsv — LFMM2 p-values (all SNPs × phenotypes)</li>"
        "<li>lfmm.q.tsv — LFMM2 q-values (all SNPs × phenotypes)</li>"
        "</ul>"
    ),
    "upload_type": "dataset",
    "access_right": "open",
    "license": "cc-by-4.0",
    "creators": [
        {"name": "Hung, Tin Hang", "affiliation": "University of Oxford"}
    ],
    "related_identifiers": [
        {
            "identifier": "https://github.com/hung-th/wytham-oak-masting",
            "relation": "isSupplementTo",
            "scheme": "url",
        }
    ],
    "keywords": [
        "Quercus robur", "masting", "GWAS", "LFMM2", "sNMF",
        "Wytham Woods", "oak", "acorn production"
    ],
}

# ── Files to upload ───────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).parent.parent

FILES = [
    REPO_ROOT / "gwas" / "lfmm2_res_WWO_37WGS_pruned_imputed.rds",
    REPO_ROOT / "gwas" / "output" / "n37" / "lfmm.z.tsv",
    REPO_ROOT / "gwas" / "output" / "n37" / "lfmm.p.tsv",
    REPO_ROOT / "gwas" / "output" / "n37" / "lfmm.q.tsv",
]

SNMF_DIR  = REPO_ROOT / "gwas" / "WWO_37WGS_pruned_snmf_sub.snmf"
SNMF_ZIP  = REPO_ROOT / "gwas" / "WWO_37WGS_pruned_snmf_sub.snmf.zip"


def zip_snmf():
    """Zip the sNMF directory if not already done."""
    if SNMF_ZIP.exists():
        print(f"  Using existing {SNMF_ZIP.name}")
        return
    print(f"  Zipping {SNMF_DIR.name} ...", end=" ", flush=True)
    with zipfile.ZipFile(SNMF_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in SNMF_DIR.rglob("*"):
            zf.write(f, f.relative_to(SNMF_DIR.parent))
    print(f"done ({SNMF_ZIP.stat().st_size / 1e6:.0f} MB)")


def create_deposition(base_url, token):
    r = requests.post(
        f"{base_url}/deposit/depositions",
        params={"access_token": token},
        json={},
        headers={"Content-Type": "application/json"},
    )
    r.raise_for_status()
    dep = r.json()
    print(f"  Created deposition: {dep['id']}")
    return dep


def upload_file(bucket_url, token, path):
    size_mb = path.stat().st_size / 1e6
    print(f"  Uploading {path.name} ({size_mb:.0f} MB) ...", end=" ", flush=True)
    with open(path, "rb") as f:
        r = requests.put(
            f"{bucket_url}/{path.name}",
            params={"access_token": token},
            data=f,
        )
    r.raise_for_status()
    print("done")


def set_metadata(base_url, token, dep_id):
    r = requests.put(
        f"{base_url}/deposit/depositions/{dep_id}",
        params={"access_token": token},
        json={"metadata": METADATA},
        headers={"Content-Type": "application/json"},
    )
    r.raise_for_status()
    print("  Metadata set.")


def publish(base_url, token, dep_id):
    r = requests.post(
        f"{base_url}/deposit/depositions/{dep_id}/actions/publish",
        params={"access_token": token},
    )
    r.raise_for_status()
    doi = r.json()["doi"]
    print(f"\n  Published! DOI: {doi}")
    return doi


def main():
    parser = argparse.ArgumentParser(description="Upload large GWAS files to Zenodo.")
    parser.add_argument("--token",   required=True, help="Zenodo access token")
    parser.add_argument("--sandbox", action="store_true",
                        help="Use sandbox.zenodo.org for testing")
    parser.add_argument("--publish", action="store_true",
                        help="Publish the deposition and mint the DOI")
    args = parser.parse_args()

    base_url = (
        "https://sandbox.zenodo.org/api"
        if args.sandbox
        else "https://zenodo.org/api"
    )
    env = "SANDBOX" if args.sandbox else "PRODUCTION"
    print(f"\n── Zenodo upload [{env}] ──────────────────────────────")

    # Check all source files exist
    all_files = FILES + [SNMF_ZIP if SNMF_ZIP.exists() else SNMF_DIR]
    missing = [f for f in FILES if not f.exists()]
    if not SNMF_DIR.exists() and not SNMF_ZIP.exists():
        missing.append(SNMF_DIR)
    if missing:
        print("ERROR: missing files:")
        for f in missing:
            print(f"  {f}")
        sys.exit(1)

    # Zip snmf directory
    print("\n[1/4] Preparing files")
    zip_snmf()

    # Create deposition
    print("\n[2/4] Creating deposition")
    dep = create_deposition(base_url, args.token)
    dep_id     = dep["id"]
    bucket_url = dep["links"]["bucket"]
    html_url   = dep["links"]["html"]

    # Upload files
    print("\n[3/4] Uploading files")
    for path in FILES + [SNMF_ZIP]:
        upload_file(bucket_url, args.token, path)

    # Set metadata
    print("\n[4/4] Setting metadata")
    set_metadata(base_url, args.token, dep_id)

    print(f"\n  Review your deposition at:\n  {html_url}\n")

    if args.publish:
        print("── Publishing ────────────────────────────────────────")
        doi = publish(base_url, args.token, dep_id)
        print(f"\n  Add this DOI to README.md:\n  https://doi.org/{doi}")
    else:
        print("  Run with --publish to mint the DOI when ready.")


if __name__ == "__main__":
    main()
