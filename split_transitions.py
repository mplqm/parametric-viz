import json
import os

SRC     = "transitions.json"
OUT_DIR = "transitions"

os.makedirs(OUT_DIR, exist_ok=True)

print(f"Reading {SRC}…")
with open(SRC) as f:
    data = json.load(f)

index_transitions = []

for t in data["transitions"]:
    tid      = f"C{t['from_cluster']}-C{t['to_cluster']}"
    out_path = os.path.join(OUT_DIR, f"{tid}.json")
    with open(out_path, "w") as f:
        json.dump(t, f, separators=(",", ":"))

    index_transitions.append({
        "id":           tid,
        "from_cluster": t["from_cluster"],
        "to_cluster":   t["to_cluster"],
        "from_members": t["from_members"],
        "to_members":   t["to_members"],
        "n_residues":   t["n_residues"],
        "res_seq":      t["res_seq"],
    })

index = {
    "reference":    data["reference"],
    "n_clusters":   data["n_clusters"],
    "n_frames":     data["n_frames"],
    "interpolation": data["interpolation"],
    "transitions":  index_transitions,
}

idx_path = os.path.join(OUT_DIR, "index.json")
with open(idx_path, "w") as f:
    json.dump(index, f, separators=(",", ":"))

idx_kb  = os.path.getsize(idx_path) / 1024
avg_kb  = sum(
    os.path.getsize(os.path.join(OUT_DIR, f"C{t['from_cluster']}-C{t['to_cluster']}.json"))
    for t in data["transitions"]
) / max(len(data["transitions"]), 1) / 1024

print(f"Index : {idx_path}  ({idx_kb:.0f} KB)")
print(f"Files : {len(data['transitions'])} transitions  (avg {avg_kb:.0f} KB each)")
print(f"Done  → {OUT_DIR}/")
