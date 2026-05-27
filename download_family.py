import json
import os
import time
import urllib.request
from Bio import PDB

OUTPUT_DIR  = "pdb_files"
TARGET      = 1000
BATCH_SIZE  = 500
os.makedirs(OUTPUT_DIR, exist_ok=True)

existing = {f[:-4] for f in os.listdir(OUTPUT_DIR) if f.endswith(".pdb")}
print(f"Already on disk: {len(existing)} structures")

SEARCH_URL   = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{}.pdb"

def rcsb_query(start, rows):
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type":"terminal","service":"full_text",
                 "parameters":{"value":"adenylate kinase"}},
                {"type":"terminal","service":"text",
                 "parameters":{"attribute":"rcsb_entry_info.experimental_method",
                               "operator":"exact_match","value":"X-ray"}},
                {"type":"terminal","service":"text",
                 "parameters":{"attribute":"rcsb_entry_info.resolution_combined",
                               "operator":"less_or_equal","value":2.5}},
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                      "direction": "asc"}],
        },
    }
    data = json.dumps(query).encode()
    req  = urllib.request.Request(SEARCH_URL, data=data,
                                   headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())

# Paginate to collect up to TARGET unique IDs
print(f"Querying RCSB for up to {TARGET} adenylate kinase X-ray ≤2.5Å structures…")
all_ids   = []
seen_set  = set()
start_row = 0
total_count = None

while len(all_ids) < TARGET:
    result      = rcsb_query(start_row, BATCH_SIZE)
    if total_count is None:
        total_count = result["total_count"]
        print(f"  RCSB total matching: {total_count}")
    batch = result.get("result_set", [])
    if not batch:
        break
    for e in batch:
        pid = e["identifier"]
        if pid not in seen_set:
            seen_set.add(pid)
            all_ids.append(pid)
    if len(all_ids) >= TARGET or start_row + BATCH_SIZE >= total_count:
        break
    start_row += BATCH_SIZE
    time.sleep(0.3)   # be polite to RCSB

all_ids     = all_ids[:TARGET]
to_download = [pid for pid in all_ids if pid not in existing]

print(f"  In query window:     {len(all_ids)}")
print(f"  Already downloaded:  {len(existing) - len([p for p in existing if p not in set(all_ids)])}")
print(f"  To download:         {len(to_download)}\n")

parser     = PDB.PDBParser(QUIET=True)
downloaded = []
failed     = []

for i, pid in enumerate(to_download):
    path = os.path.join(OUTPUT_DIR, f"{pid}.pdb")
    try:
        urllib.request.urlretrieve(DOWNLOAD_URL.format(pid), path)
        struct  = parser.get_structure(pid, path)
        n_atoms = sum(1 for _ in struct.get_atoms())
        res     = struct.header.get("resolution") or "N/A"
        res_str = f"{res:.2f}" if isinstance(res, float) else str(res)
        print(f"  {pid:<8} {res_str:<8} Å   {n_atoms:>6} atoms")
        downloaded.append(pid)
    except Exception as exc:
        print(f"  {pid:<8} FAILED: {exc}")
        failed.append((pid, str(exc)))
        if os.path.exists(path):
            os.remove(path)
    if (i + 1) % 50 == 0:
        print(f"  — {i + 1}/{len(to_download)} processed —")

total = len({f[:-4] for f in os.listdir(OUTPUT_DIR) if f.endswith(".pdb")})
print(f"\n{'─'*50}")
print(f"  Downloaded:  {len(downloaded):>4} new")
print(f"  Failed:      {len(failed):>4}")
print(f"  Total:       {total:>4} structures in {OUTPUT_DIR}/")
if failed:
    print("\nFailed IDs:")
    for pid, err in failed:
        print(f"  {pid:<8} {err[:80]}")
