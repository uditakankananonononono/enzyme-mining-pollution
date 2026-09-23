"""Figures for the pesticides slice paper. v2-level scoring of labeled seqs (fast, representative)."""
import csv, json, os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyfamsa, pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
ALPHA = Alphabet.amino(); BUILDER = Builder(ALPHA); BG = Background(ALPHA)
DATA, RESULTS = "../data", "../results"
os.makedirs(f"{RESULTS}/figures", exist_ok=True)

def read_fasta_text(txt):
    seqs, name, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

refs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/reference_set.csv"))}
negs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/negative_set.csv"))}
seqs = {}
for fn in os.listdir(f"{DATA}/fasta"):
    seqs[fn[:-6]] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
split = json.load(open(f"{DATA}/split.json"))

def build_hmm(name, accs):
    fs = [pyfamsa.Sequence(a.encode(), seqs[a].encode()) for a in accs]
    if len(fs) >= 2:
        msa = pyfamsa.Aligner(threads=2).align(fs)
        tmsa = TextMSA(name=name.encode(), sequences=[pyhmmer.easel.TextSequence(name=s.id, sequence=s.sequence.decode()) for s in msa])
        hmm, _, _ = BUILDER.build_msa(tmsa.digitize(ALPHA), BG)
    else:
        hmm, _, _ = BUILDER.build(pyhmmer.easel.TextSequence(name=name.encode(), sequence=seqs[accs[0]]).digitize(ALPHA), BG)
    hmm.name = name.encode(); return hmm

fam_accs = defaultdict(list)
for acc in split["train"]: fam_accs[f"{refs[acc]['class']}:{refs[acc]['family']}"].append(acc)
pos = {k: build_hmm(k, a) for k, a in fam_accs.items()}
pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "PTE_NEIGHBOR"]
gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "GENERIC_ESTERASE"]
negh = {"NEG:PLL": build_hmm("NEG:PLL", pll), "NEG:EST": build_hmm("NEG:EST", gest)}
names = list(seqs)
block = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=a.encode(), sequence=seqs[a]).digitize(ALPHA) for a in names])
scores = {a: defaultdict(float) for a in names}
for hmm_obj, hits in zip(list(pos.values()) + list(negh.values()), pyhmmer.hmmsearch(list(pos.values()) + list(negh.values()), block, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        t = h.name.decode() if isinstance(h.name, bytes) else h.name
        scores[t][hn] = max(scores[t][hn], h.score)

def posmax(a): return max(scores[a].get(k, 0.0) for k in pos)
def negmax(a): return max(scores[a].get("NEG:PLL", 0.0), scores[a].get("NEG:EST", 0.0))
groups = {
    "train pos": ([a for a in split["train"]], "#1f77b4", "o"),
    "heldout pos": ([a for a in split["heldout"]], "#2ca02c", "s"),
    "tune neg (PLL/EST)": ([a for a in split["negatives_tune"]], "#d62728", "x"),
    "eval neg (PLL/EST)": ([a for a in split["negatives_eval"]], "#ff7f0e", "^"),
}
fig, ax = plt.subplots(figsize=(7.2, 5.2))
for label, (accs, c, m) in groups.items():
    ax.scatter([negmax(a) for a in accs], [posmax(a) for a in accs], c=c, marker=m, s=42, alpha=0.85, label=f"{label} (n={len(accs)})")
ax.plot([0, 1200], [0, 1200], "k--", lw=0.8, alpha=0.5, label="pos = neg")
ax.set_xlabel("best negative-profile score (bits)"); ax.set_ylabel("best positive-profile score (bits)")
ax.set_title("Discriminator separation: labeled positives vs documented false-positive modes (v2 scoring)")
ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{RESULTS}/figures/fig1_discriminator.png", dpi=160); plt.close(fig)

# Fig 2: recall boundary
bd = json.load(open(f"{RESULTS}/boundary_diagnosis.json"))["jaccard_vs_train"]
fig, ax = plt.subplots(figsize=(7.2, 4.6))
import random
random.seed(7)
for status, c, y in [("recovered", "#2ca02c", 1), ("MISSED", "#d62728", 0)]:
    xs = [r["jaccard"] for r in bd if r["status"] == status]
    ys = [y + random.uniform(-0.08, 0.08) for _ in xs]
    ax.scatter(xs, ys, c=c, s=70, label=f"{status} (n={len(xs)})", zorder=3)
    for r, yy in zip([r for r in bd if r["status"] == status], ys):
        ax.annotate(r["acc"], (r["jaccard"], yy), textcoords="offset points", xytext=(5, 4), fontsize=7)
ax.axvline(0.15, color="k", ls="--", lw=1, alpha=0.6)
ax.text(0.155, 0.5, "coverage boundary\n(jacc ~0.15-0.28)", fontsize=8)
ax.set_yticks([0, 1]); ax.set_yticklabels(["MISSED", "recovered"])
ax.set_xlabel("3-mer Jaccard similarity to nearest TRAINING enzyme")
ax.set_title("Recall is bounded by fold coverage of the characterized canon (held-out, v4)")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f"{RESULTS}/figures/fig2_recall_boundary.png", dpi=160); plt.close(fig)

# Fig 3: version metrics
vers = ["v1", "v2", "v3", "v4"]
recall = [0.733, 0.692, 0.462, 0.692]; auroc = [0.912, 0.916, 0.911, 0.906]
pc = [29/32, 1.0, 1.0, 1.0]
fig, ax = plt.subplots(figsize=(7.2, 4.2))
x = range(len(vers)); w = 0.27
ax.bar([i - w for i in x], recall, w, label="GATE 1 recall (held-out)", color="#1f77b4")
ax.bar(list(x), auroc, w, label="GATE 2 AUROC", color="#2ca02c")
ax.bar([i + w for i in x], pc, w, label="positive control (train recovery)", color="#9467bd")
ax.axhline(0.90, color="r", ls="--", lw=1); ax.text(2.4, 0.905, "gate threshold 0.90", color="r", fontsize=8)
ax.set_xticks(list(x)); ax.set_xticklabels(vers); ax.set_ylim(0, 1.05); ax.set_ylabel("metric")
ax.set_title("Validation across pipeline designs (thresholds frozen pre-measurement)")
ax.legend(fontsize=8, loc="lower right"); ax.grid(axis="y", alpha=0.25)
fig.tight_layout(); fig.savefig(f"{RESULTS}/figures/fig3_versions.png", dpi=160); plt.close(fig)

# Fig 4: reference set composition
from collections import Counter
cnt = Counter((r["class"], r["family"], r["tier"]) for r in refs.values())
fams = sorted({f for _, f, _ in cnt})
ver = [cnt[("OP", f, "VERIFIED")] for f in fams]; thin_op = [cnt[("OP", f, "THIN")] for f in fams]
vpyr = [cnt[("PYR", f, "VERIFIED")] for f in fams]; tpyr = [cnt[("PYR", f, "THIN")] for f in fams]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
x = range(len(fams))
b1 = ax.bar(x, ver, label="OP verified", color="#1f77b4")
b2 = ax.bar(x, thin_op, bottom=ver, label="OP thin", color="#9ecae1")
b3 = ax.bar(x, vpyr, bottom=[a + b for a, b in zip(ver, thin_op)], label="PYR verified", color="#2ca02c")
b4 = ax.bar(x, tpyr, bottom=[a + b + c for a, b, c in zip(ver, thin_op, vpyr)], label="PYR thin", color="#a1d99b")
ax.set_xticks(list(x)); ax.set_xticklabels(fams, rotation=30, ha="right", fontsize=8)
ax.set_ylabel("enzymes"); ax.set_title("Reference canon: 76 enzymes across OP and pyrethroid classes (tiered)")
ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.25)
fig.tight_layout(); fig.savefig(f"{RESULTS}/figures/fig4_reference_set.png", dpi=160); plt.close(fig)

# Fig 5: the intrinsic collision - Pye3 vs YeiG
fig, ax = plt.subplots(figsize=(6.4, 3.6))
items = ["ACJ07038.1\n(Pye3, true\npyrethroid est.)", "B6VG94\n(YeiG, S-formylglutathione\nhydrolase, FP mode)"]
vals = [posmax("ACJ07038.1"), posmax("B6VG94")]
bars = ax.bar(items, vals, color=["#2ca02c", "#d62728"], width=0.5)
for b, v in zip(bars, vals): ax.text(b.get_x() + b.get_width()/2, v + 5, f"{v:.1f} bits", ha="center", fontsize=10)
ax.set_ylabel("best positive-profile score (bits)")
ax.set_title("Documented intrinsic collision: sequence-only scores cannot\nseparate Pye3 from YeiG-family esterases")
ax.grid(axis="y", alpha=0.25)
fig.tight_layout(); fig.savefig(f"{RESULTS}/figures/fig5_collision.png", dpi=160); plt.close(fig)
print("figures written:", os.listdir(f"{RESULTS}/figures"))
