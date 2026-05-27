import json
import os
import urllib.request
from Bio import PDB

OUTPUT_DIR = "pdb_files"
os.makedirs(OUTPUT_DIR, exist_ok=True)

existing = {f[:-4] for f in os.listdir(OUTPUT_DIR) if f.endswith(".pdb")}
print(f"Already on disk: {len(existing)} structures")

SEARCH_URL   = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{}.pdb"

query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": "adenylate kinase"},
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.experimental_method",
                    "operator": "exact_match",
                    "value": "X-ray",
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": 2.5,
                },
            },
        ],
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {"start": 0, "rows": 500},
        "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
    },
}

print("Querying RCSB for adenylate kinase X-ray structures (resolution ≤ 2.5 Å)…")
data = json.dumps(query).encode()
req  = urllib.request.Request(SEARCH_URL, data=data, headers={"Content-Type": "application/json"})
with urllib.request.urlopen(req) as resp:
    result = json.loads(resp.read())

all_ids     = [e["identifier"] for e in result["result_set"]]
to_download = [pid for pid in all_ids if pid not in existing]

print(f"  RCSB total matching: {result['total_count']}")
print(f"  In query window:     {len(all_ids)}")
print(f"  Already downloaded:  {len(existing)}")
print(f"  To download:         {len(to_download)}\n")

parser     = PDB.PDBParser(QUIET=True)
downloaded = []
failed     = []

for pid in to_download:
    path = os.path.join(OUTPUT_DIR, f"{pid}.pdb")
    try:
        urllib.request.urlretrieve(DOWNLOAD_URL.format(pid), path)
        struct     = parser.get_structure(pid, path)
        n_atoms    = sum(1 for _ in struct.get_atoms())
        resolution = struct.header.get("resolution") or "N/A"
        res_str    = f"{resolution:.2f}" if isinstance(resolution, float) else str(resolution)
        print(f"  {pid:<8} {res_str:<8} Å   {n_atoms:>6} atoms")
        downloaded.append(pid)
    except Exception as exc:
        print(f"  {pid:<8} FAILED: {exc}")
        failed.append((pid, str(exc)))
        if os.path.exists(path):
            os.remove(path)

total = len(existing) + len(downloaded)
print(f"\n{'─'*50}")
print(f"  Downloaded:  {len(downloaded):>4} new")
print(f"  Failed:      {len(failed):>4}")
print(f"  Total:       {total:>4} structures in {OUTPUT_DIR}/")
