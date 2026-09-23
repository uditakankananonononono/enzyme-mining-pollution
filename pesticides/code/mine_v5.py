"""Mining v2: v4 platform (validated), bg corpus excludes LABELED accessions (canon-duplicate leak fix). Original v4 docstring preserved below.
Final mining run: v4 platform (validated) applied to MGnify MAG proteomes + UniProt bg corpus.
Candidates: CONFIDENT calls, excluding labeled/decoys, GATE 3 novelty (jacc<0.5 proxy for <80% id vs VERIFIED canon, flagged).
Outputs: results/candidates_v2.csv + results/mining_summary.json."""
import csv, json, os, re, gzip, glob
from collections import defaultdict
import pyfamsa, pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
ALPHA = Alphabet.amino(); BUILDER = Builder(ALPHA); BG = Background(ALPHA)
DATA, CORPUS, RESULTS = "../data", "../data/corpus", "../results"

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
for acc in train: fam_accs[f"{refs[acc]['class']}:{refs[acc]['family']}"].append(acc)
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
expanded, exp_count = {}, {}
for f, members in fam_accs.items():
    hits = sorted(((famscore_bg[n].get(f, 0.0), n) for n in bgnames), reverse=True)
    add = [n for s, n in hits if s >= fam_T0[f]][:60]
    exp_count[f] = len(add)
    expanded["FAMX:" + f] = build_hmm("FAMX:" + f, members + add, allseq)

v3 = json.load(open(f"{RESULTS}/validation_v3.json"))
sub_members = {cid: v["train_members"] for cid, v in v3["subfamilies"].items()}
sub_hmms = {cid: build_hmm(cid, m, allseq) for cid, m in sub_members.items()}
pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "PTE_NEIGHBOR"]
gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "GENERIC_ESTERASE"]
neg_hmms = {"NEG:PLL": build_hmm("NEG:PLL", pll, allseq), "NEG:EST": build_hmm("NEG:EST", gest, allseq)}

# corpus with source tags
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
prof_members = {("FAMX:" + f): m for f, m in fam_accs.items()}; prof_members.update(sub_members)
all_hmms = list(expanded.values()) + list(sub_hmms.values()) + list(neg_hmms.values())
scores = {n: defaultdict(float) for n in names}
for hmm_obj, hits in zip(all_hmms, pyhmmer.hmmsearch(all_hmms, block, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        t = h.name.decode() if isinstance(h.name, bytes) else h.name
        scores[t][hn] = max(scores[t][hn], h.score)
print("corpus:", len(names), "HMMs:", len(all_hmms))

def negmax(n): return max(scores[n].get("NEG:PLL", 0.0), scores[n].get("NEG:EST", 0.0))
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

# re-validate gates on this exact frozen model
tr = [a for a in train if call(f"LAB|{a}")]
ho = [a for a in heldout if call(f"LAB|{a}")]
fpe = [a for a in split["negatives_eval"] if call(f"LAB|{a}")]
prec = len(ho) / (len(ho) + len(fpe)) if (len(ho) or len(fpe)) else 1.0
print(f"re-check: PC {len(tr)}/{len(train)}, GATE1 {len(ho)}/{len(heldout)}, GATE2 prec {prec:.3f}")

# novelty: max 3mer-jaccard vs VERIFIED references
VER = [a for a, r in refs.items() if r["tier"] == "VERIFIED"]
verk = {a: set(seqs[a][i:i+3] for i in range(len(seqs[a])-2)) for a in VER}
def max_jac(s):
    A = set(s[i:i+3] for i in range(len(s)-2))
    best = 0.0
    for a, K in verk.items():
        j = len(A & K) / len(A | K)
        if j > best: best = j
    return best

GXSXG = re.compile(r"G.S.G")
rows = []
for n in names:
    if n.startswith("LAB|") or n.startswith("DECOY|"): continue
    pid = call(n)
    if not pid: continue
    s = scores[n][pid]; nm = negmax(n)
    rows.append({"id": n, "source": source[n], "length": len(corpus[n]),
                 "profile": pid, "class": pid.split(":")[1] if pid.startswith("FAMX:") else refs.get(sub_members.get(pid, ["?"])[0], {}).get("class", "?"),
                 "score": round(s, 1), "neg_score": round(nm, 1), "margin": round(s - nm, 1),
                 "gxsxg": bool(GXSXG.search(corpus[n])),
                 "max_jaccard_vs_canon": round(max_jac(corpus[n]), 3)})
rows.sort(key=lambda r: -r["margin"])
with open(f"{RESULTS}/candidates_v2.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["id"]); w.writeheader(); w.writerows(rows)
summary = {"corpus_total": len(names),
           "mgnify_proteins": sum(1 for n in names if n.startswith("MGYG")),
           "candidates": len(rows),
           "candidates_by_class": dict(__import__("collections").Counter(r["class"] for r in rows)),
           "candidates_by_source": dict(__import__("collections").Counter(r["source"].split(":")[0] for r in rows)),
           "novel_jacc_lt_0.5": sum(1 for r in rows if r["max_jaccard_vs_canon"] < 0.5),
           "revalidation": {"pc": f"{len(tr)}/{len(train)}", "gate1": f"{len(ho)}/{len(heldout)}", "gate2_precision": round(prec, 3)}}
json.dump(summary, open(f"{RESULTS}/mining_summary_v2.json", "w"), indent=2)
print(json.dumps(summary, indent=1))
