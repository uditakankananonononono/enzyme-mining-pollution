"""Industrial-pollutant enzyme mining pipeline (builder 21 replacement slice).
Adapted from validated pesticides/ platform (mine_v5.py). Single frozen run:
Stage 1: family profile HMMs (training VERIFIED+THIN references, split frozen pre-mining).
Stage 2: bg-corpus-expanded HMMs (FAMX) + negative-mode HMMs (tune negatives only).
Stage 3: positive control, GATE1 held-out recall, GATE2 discriminator precision/AUROC,
         then candidate calls + GATE3 novelty proxy on the same frozen model.
Classes: HLD (HALOALKANE/HALOACID), RHD (RHD_PAH/RHD_BPH/RHD_BTX/RHD_OTHER)."""
import csv, json, os, re, gzip, glob, sys
from collections import defaultdict, Counter
import pyfamsa, pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
ALPHA = Alphabet.amino(); BUILDER = Builder(ALPHA); BG = Background(ALPHA)
BASE = os.path.dirname(os.path.abspath(__file__)) + "/.."
DATA, CORPUS, RESULTS = f"{BASE}/data", f"{BASE}/data/corpus", f"{BASE}/results"
os.makedirs(RESULTS, exist_ok=True)

def read_fasta_text(txt):
    seqs, name, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def build_hmm(name, keys, allseq):
    fs = [pyfamsa.Sequence(k.encode(), allseq[k].encode()) for k in keys]
    if len(fs) >= 2:
        msa = pyfamsa.Aligner(threads=2).align(fs)
        tmsa = TextMSA(name=name.encode(), sequences=[pyhmmer.easel.TextSequence(name=s.id, sequence=s.sequence.decode()) for s in msa])
        hmm, _, _ = BUILDER.build_msa(tmsa.digitize(ALPHA), BG)
    else:
        hmm, _, _ = BUILDER.build(pyhmmer.easel.TextSequence(name=name.encode(), sequence=allseq[keys[0]]).digitize(ALPHA), BG)
    hmm.name = name.encode(); return hmm

refs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/reference_set.csv"))}
negs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/negative_set.csv"))}
LABELED = set(refs) | set(negs)
seqs = {}
for fn in os.listdir(f"{DATA}/fasta"):
    seqs[fn[:-6]] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
split = json.load(open(f"{DATA}/split.json")); train, heldout = split["train"], split["heldout"]

def bg_acc(n):
    return n.split("|")[1] if n.startswith(("sp|", "tr|")) and n.count("|") >= 2 else n
bg = {}
for fn in os.listdir(CORPUS):
    if fn.startswith("bg_"):
        for n, s in read_fasta_text(open(f"{CORPUS}/{fn}").read()).items():
            if bg_acc(n) not in LABELED: bg[n] = s
allseq = dict(seqs); allseq.update(bg)

fam_accs = defaultdict(list)
for acc in train: fam_accs[f"{refs[acc]['cls']}:{refs[acc]['family']}"].append(acc)
fam_hmms = {k: build_hmm(k, a, allseq) for k, a in fam_accs.items()}

bgnames = list(bg)
bgblock = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=bg[n]).digitize(ALPHA) for n in bgnames])
labblock = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=a.encode(), sequence=seqs[a]).digitize(ALPHA) for a in seqs])
famscore_bg = {n: {} for n in bgnames}; famscore_lab = {a: {} for a in seqs}
fh = list(fam_hmms.values())
for hmm_obj, hits in zip(fh, pyhmmer.hmmsearch(fh, bgblock, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        t = h.name.decode() if isinstance(h.name, bytes) else h.name
        famscore_bg[t][hn] = max(famscore_bg[t].get(hn, 0.0), h.score)
for hmm_obj, hits in zip(fh, pyhmmer.hmmsearch(fh, labblock, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        t = h.name.decode() if isinstance(h.name, bytes) else h.name
        famscore_lab[t][hn] = max(famscore_lab[t].get(hn, 0.0), h.score)
fam_T0 = {f: 0.8 * min(famscore_lab[a].get(f, 0.0) for a in m) for f, m in fam_accs.items()}
expanded = {}
for f, members in fam_accs.items():
    hits = sorted(((famscore_bg[n].get(f, 0.0), n) for n in bgnames), reverse=True)
    add = [n for s, n in hits if s >= fam_T0[f]][:60]
    expanded["FAMX:" + f] = build_hmm("FAMX:" + f, members + add, allseq)

# negative-mode HMMs from TUNE negatives only
hldn = [a for a in split["negatives_tune"] if negs[a]["neg_class"] in ("HAD_PHOSPHATASE", "EPOXIDE_HYDROLASE")]
rhdn = [a for a in split["negatives_tune"] if negs[a]["neg_class"] in ("RIESKE_FERREDOXIN", "RING_CLEAVAGE_DO")]
neg_hmms = {"NEG:HLD": build_hmm("NEG:HLD", hldn, allseq), "NEG:RHD": build_hmm("NEG:RHD", rhdn, allseq)}

corpus, source = {}, {}
for fn in os.listdir(CORPUS):
    p = f"{CORPUS}/{fn}"
    if fn.endswith(".gz"):
        for n, s in read_fasta_text(gzip.open(p, "rt").read()).items(): corpus["DECOY|" + n] = s; source["DECOY|" + n] = fn
    elif fn.startswith("bg_"):
        for n, s in read_fasta_text(open(p).read()).items():
            if bg_acc(n) not in LABELED: corpus["BG|" + n] = s; source["BG|" + n] = fn
for p in glob.glob(f"{CORPUS}/mgnify/*.faa"):
    mag = os.path.basename(p)[:-4]
    for n, s in read_fasta_text(open(p).read()).items():
        key = f"{mag}|{n}"; corpus[key] = s; source[key] = f"mgnify:{mag}"
for acc, s in allseq.items(): corpus[f"LAB|{acc}"] = s; source[f"LAB|{acc}"] = "labeled"

names = list(corpus)
block = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=corpus[n]).digitize(ALPHA) for n in names])
prof_members = {("FAMX:" + f): m for f, m in fam_accs.items()}
all_hmms = list(expanded.values()) + list(neg_hmms.values())
scores = {n: defaultdict(float) for n in names}
for hmm_obj, hits in zip(all_hmms, pyhmmer.hmmsearch(all_hmms, block, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        t = h.name.decode() if isinstance(h.name, bytes) else h.name
        scores[t][hn] = max(scores[t][hn], h.score)
print("corpus:", len(names), "HMMs:", len(all_hmms), flush=True)

def negmax(n):
    return max(scores[n].get("NEG:HLD", 0.0), scores[n].get("NEG:RHD", 0.0))
T, D = {}, {}
for pid, members in prof_members.items():
    own = [scores[f"LAB|{a}"].get(pid, 0.0) for a in members]
    T[pid] = 0.8 * min(own)
    deltas = sorted((scores[f"LAB|{a}"].get(pid, 0.0) - negmax(f"LAB|{a}") for a in split["negatives_tune"]), reverse=True)
    D[pid] = max((deltas[1] if len(deltas) > 1 else 0.0) + 10.0, 0.0)

def call(n):
    best, bp = -1, None
    for pid in prof_members:
        s = scores[n].get(pid, 0.0)
        if s >= T[pid] and (s - negmax(n)) >= D[pid] and s > best: best, bp = s, pid
    return bp

tr = [a for a in train if call(f"LAB|{a}")]
ho_all = [a for a in heldout if call(f"LAB|{a}")]
ho_ver = [a for a in heldout if refs[a]["tier"] == "VERIFIED"]
ho_ver_hit = [a for a in ho_ver if call(f"LAB|{a}")]
fpe = [a for a in split["negatives_eval"] if call(f"LAB|{a}")]
prec = len(ho_ver_hit) / (len(ho_ver_hit) + len(fpe)) if (len(ho_ver_hit) or len(fpe)) else 1.0
# AUROC on margin: heldout VERIFIED (pos) vs eval negatives (neg)
def margin(a):
    pid_scores = {pid: scores[f"LAB|{a}"].get(pid, 0.0) for pid in prof_members}
    best = max(pid_scores.values()) if pid_scores else 0.0
    return best - negmax(f"LAB|{a}")
pairs = [(margin(a), 1) for a in ho_ver] + [(margin(a), 0) for a in split["negatives_eval"]]
pos = [m for m, y in pairs if y == 1]; neg = [m for m, y in pairs if y == 0]
wins = sum(1 if p > n else 0.5 if p == n else 0 for p in pos for n in neg)
auroc = wins / (len(pos) * len(neg)) if pos and neg else float("nan")
print(f"PC {len(tr)}/{len(train)} | GATE1 recall(VERIFIED heldout) {len(ho_ver_hit)}/{len(ho_ver)} | GATE2 prec {prec:.3f} AUROC {auroc:.3f}", flush=True)

VER = [a for a, r in refs.items() if r["tier"] == "VERIFIED"]
verk = {a: set(seqs[a][i:i+3] for i in range(len(seqs[a])-2)) for a in VER}
def max_jac(s):
    A = set(s[i:i+3] for i in range(len(s)-2))
    best = 0.0
    for a, K in verk.items():
        j = len(A & K) / len(A | K)
        if j > best: best = j
    return best

RIESKE = re.compile(r"C.H.{15,20}C..H")
rows = []
for n in names:
    if n.startswith("LAB|") or n.startswith("DECOY|"): continue
    pid = call(n)
    if not pid: continue
    s = scores[n][pid]; nm = negmax(n)
    cls = pid.split(":")[1]
    rows.append({"id": n, "source": source[n], "length": len(corpus[n]),
                 "profile": pid, "cls": cls,
                 "score": round(s, 1), "neg_score": round(nm, 1), "margin": round(s - nm, 1),
                 "rieske_motif": bool(RIESKE.search(corpus[n])) if cls == "RHD" else "",
                 "max_jaccard_vs_canon": round(max_jac(corpus[n]), 3)})
rows.sort(key=lambda r: -r["margin"])
with open(f"{RESULTS}/candidates.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id","source","length","profile","cls","score","neg_score","margin","rieske_motif","max_jaccard_vs_canon"])
    w.writeheader(); w.writerows(rows)
validation = {"seed": 20260923,
    "positive_control": {"recovered": len(tr), "total": len(train)},
    "gate1_recall_verified_heldout": {"recovered": len(ho_ver_hit), "total": len(ho_ver),
        "recall": round(len(ho_ver_hit)/len(ho_ver), 3) if ho_ver else None,
        "misses": [a for a in ho_ver if a not in ho_ver_hit]},
    "gate1_all_heldout": {"recovered": len(ho_all), "total": len(heldout)},
    "gate2": {"precision": round(prec, 3), "auroc": round(auroc, 3),
        "eval_negatives_called": fpe, "n_eval_negatives": len(split["negatives_eval"])},
    "thresholds_T": {k: round(v,1) for k,v in T.items()}, "thresholds_D": {k: round(v,1) for k,v in D.items()}}
json.dump(validation, open(f"{RESULTS}/validation.json", "w"), indent=1)
summary = {"corpus_total": len(names),
    "mgnify_proteins": sum(1 for n in names if n.startswith("MGYG")),
    "bg_proteins": sum(1 for n in names if n.startswith("BG|")),
    "decoy_proteins": sum(1 for n in names if n.startswith("DECOY|")),
    "candidates": len(rows),
    "candidates_by_class": dict(Counter(r["cls"] for r in rows)),
    "candidates_by_source": dict(Counter(r["source"].split(":")[0] for r in rows)),
    "novel_jacc_lt_0.5": sum(1 for r in rows if r["max_jaccard_vs_canon"] < 0.5),
    "validation": validation}
json.dump(summary, open(f"{RESULTS}/mining_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "validation"}, indent=1))
