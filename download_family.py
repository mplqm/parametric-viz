import json
import os
import urllib.request
from Bio import PDB

OUTPUT_DIR = "pdb_files"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DOWNLOAD_URL = "https://files.rcsb.org/download/{}.pdb"

query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": "adenylate kinase"
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.experimental_method",
                    "operator": "exact_match",
                    "value": "X-ray"
                }
            }
        ]
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {"start": 0, "rows": 20},
        "sort": [
            {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
        ]
    }
}

print("Querying RCSB for adenylate kinase X-ray structures...")
data = json.dumps(query).encode("utf-8")
req = urllib.request.Request(
    SEARCH_URL,
    data=data,
    headers={"Content-Type": "application/json"}
)
with urllib.request.urlopen(req) as resp:
    result = json.loads(resp.read())

pdb_ids = [entry["identifier"] for entry in result["result_set"]]
print(f"Found {result['total_count']} total structures — downloading top {len(pdb_ids)} by resolution\n")
print(f"{'PDB ID':<8} {'Resolution (Å)':<16} {'Atoms':>6}")
print("-" * 32)

parser = PDB.PDBParser(QUIET=True)

for pdb_id in pdb_ids:
    pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}.pdb")
    try:
        urllib.request.urlretrieve(DOWNLOAD_URL.format(pdb_id), pdb_file)
        structure = parser.get_structure(pdb_id, pdb_file)
        atom_count = sum(1 for _ in structure.get_atoms())
        resolution = structure.header.get("resolution") or "N/A"
        res_str = f"{resolution:.2f}" if isinstance(resolution, float) else str(resolution)
        print(f"{pdb_id:<8} {res_str:<16} {atom_count:>6}")
    except Exception as e:
        print(f"{pdb_id:<8} ERROR: {e}")
