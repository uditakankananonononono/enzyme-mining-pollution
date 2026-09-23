"""GATE3 vetting (cluster + canon-collision audit) + figures, industrial slice."""
import csv, json, os, re, gzip, glob
from collections import Counter, defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BASE = os.path.dirname(os.path.abspath(__file__)) + "/.."
DATA, CORPUS, RES = f"{BASE}/data", f"{BASE}/data/corpus", f"{BASE}/results"
os.makedirs(f"{RES}/figures", exist_ok=True)

cands = list(csv.DictReader(open(f"{RES}/candidates.csv")))
refs = list(csv.DictReader(open(f"{DATA}/reference_set.csv")))
def read_fasta_text(txt):
    seqs, name, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs
corpus = {}
for fn in os.listdir(CORPUS):
    p = f"{CORPUS}/{fn}"
    if fn.endswith(".gz"):
        for n, s in read_fasta_text(gzip.open(p, "rt").read()).items(): corpus["DECOY|" + n] = s
    elif fn.startswith("bg_"):
        for n, s in read_fasta_text(open(p).read()).items(): corpus["BG|" + n] = s
for p in glob.glob(f"{CORPUS}/mgnify/*.faa"):
    mag = os.path.basename(p)[:-4]
    for n, s in read_fasta_text(open(p).read()).items(): corpus[f"{mag}|{n}"] = s

def kset(s, k=3): return set(s[i:i+k] for i in range(len(s)-k+1))
# greedy clustering of candidates at jacc >= 0.5
clusters = []
for r in cands:
    s = corpus.get(r["id"])
    if not s: continue
    ks = kset(s); placed = False
    for c in clusters:
        j = len(ks & c["ks"]) / len(ks | c["ks"])
        if j >= 0.5:
            c["members"].append(r["id"]); placed = True; break
    if not placed: clusters.append({"rep": r["id"], "ks": ks, "members": [r["id"]], "cls": r["cls"], "margin": float(r["margin"])})
clusters.sort(key=lambda c: -c["margin"])
vet = []
for i, c in enumerate(clusters, 1):
    rep = c["members"][0]
    row = next(r for r in cands if r["id"] == rep)
    status = "NOVEL" if float(row["max_jaccard_vs_canon"]) < 0.5 else "HOMOLOG-OF-CANON"
    vet.append({"cluster": f"IND-C{i:03d}", "rep_id": rep, "cls": c["cls"], "n_members": len(c["members"]),
                "margin": c["margin"], "score": row["score"], "source": row["source"],
                "rieske_motif": row["rieske_motif"], "max_jaccard_vs_canon": row["max_jaccard_vs_canon"],
                "gate3_status": status})
w = csv.DictWriter(open(f"{RES}/vetted_candidates.csv", "w", newline=""), fieldnames=list(vet[0].keys()))
w.writeheader(); w.writerows(vet)
print("clusters:", len(vet), Counter(v["gate3_status"] for v in vet), Counter(v["cls"] for v in vet))

# --- figures ---
val1 = json.load(open(f"{RES}/validation_v1.json")); val2 = json.load(open(f"{RES}/validation.json"))
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
# fig1: discriminator separation - margins of heldout VERIFIED vs eval negatives (v2)
split = json.load(open(f"{DATA}/split.json"))
# recompute margins from candidates? use validation summary instead: plot AUROC/precision bars
ax = axes[0]
ax.bar(["GATE1 recall\n(bar 0.90)", "GATE2 precision\n(bar 0.90)", "GATE2 AUROC\n(bar 0.90)"],
       [val2["gate1_recall_verified_heldout"]["recall"], val2["gate2"]["precision"], val2["gate2"]["auroc"]],
       color=["#2b6cb0", "#2f855a", "#2f855a"])
ax.axhline(0.9, color="red", ls="--", lw=1); ax.set_ylim(0, 1.05); ax.set_title("Frozen-model gate metrics (v2)")
ax = axes[1]
ax.bar(["v1 (contaminated\nreference set)", "v2 (curated\nalpha-only)"],
       [val1["gate1_recall_verified_heldout"]["recall"], val2["gate1_recall_verified_heldout"]["recall"]],
       color=["#c53030", "#2b6cb0"])
ax.axhline(0.9, color="red", ls="--", lw=1); ax.set_ylim(0, 1.05); ax.set_title("GATE1 recall: v1 vs v2")
fig.tight_layout(); fig.savefig(f"{RES}/figures/fig1_gates.png", dpi=150); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
ax = axes[0]
fam_counts = Counter((r["cls"], r["tier"]) for r in refs)
labels = ["HLD VERIFIED", "HLD THIN", "RHD VERIFIED", "RHD THIN"]
vals = [fam_counts[("HLD","VERIFIED")], fam_counts[("HLD","THIN")], fam_counts[("RHD","VERIFIED")], fam_counts[("RHD","THIN")]]
ax.bar(labels, vals, color=["#2b6cb0", "#90cdf4", "#2f855a", "#9ae6b4"]); ax.set_title("Reference set composition (v2)")
ax.tick_params(axis='x', rotation=20)
ax = axes[1]
src = Counter(r["source"].split(":")[0] for r in cands)
ax.bar(range(len(src)), list(src.values()), color="#6b46c1")
ax.set_xticks(range(len(src))); ax.set_xticklabels([k.replace(".fasta","").replace("bg_","bg:\n") for k in src], fontsize=8)
ax.set_title("Candidate sources"); ax.set_ylabel("candidates")
fig.tight_layout(); fig.savefig(f"{RES}/figures/fig2_sets_sources.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 4.2))
cls_counts = Counter(r["cls"] for r in cands)
novel = Counter(r["cls"] for r in cands if float(r["max_jaccard_vs_canon"]) < 0.5)
x = range(len(cls_counts))
ax.bar(x, [cls_counts[c] for c in cls_counts], color="#b7791f", label="called")
ax.bar(x, [novel[c] for c in cls_counts], color="#2f855a", label="novel (jacc<0.5)")
ax.set_xticks(list(x)); ax.set_xticklabels(list(cls_counts)); ax.legend()
ax.set_title(f"Candidates by class ({len(cands)} total, {len(vet)} clusters)")
fig.tight_layout(); fig.savefig(f"{RES}/figures/fig3_candidates.png", dpi=150); plt.close(fig)
print("figures written")
