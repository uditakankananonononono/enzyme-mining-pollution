"""v4: v2 family HMMs + v3 UniRef90 sub-family HMMs, PLUS one jackhmmer-style
self-expansion round: family HMM hits >= family T on the UNLABELED background corpus
(bg_*.fasta only) are added to that family's MSA (cap 60), HMM rebuilt.
Scoring = max over {expanded family HMMs} U {sub-family HMMs}.
Thresholds frozen identically to v2/v3 rules. One held-out measurement.
Pre-registered structural fixes: (a) cross-subtype generalization needs family level,
(b) divergent singletons need broader profiles from unlabeled homologs."""
import csv, json, os, re, gzip
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
seqs = {}
for fn in os.listdir(f"{DATA}/fasta"):
    seqs[fn[:-6]] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
split = json.load(open(f"{DATA}/split.json")); assert split["version"] == 2
train, heldout = split["train"], split["heldout"]

# background (unlabeled) corpus
bg = {}
for fn in os.listdir(CORPUS):
    if fn.startswith("bg_"): bg.update(read_fasta_text(open(f"{CORPUS}/{fn}").read()))
allseq = dict(seqs); allseq.update(bg)

# family HMMs (v2)
fam_accs = defaultdict(list)
for acc in train: fam_accs[f"{refs[acc]['class']}:{refs[acc]['family']}"].append(acc)
fam_hmms = {k: build_hmm(k, a, allseq) for k, a in fam_accs.items()}

# round-0 scoring of background with family HMMs
bgnames = list(bg)
bgblock = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=bg[n]).digitize(ALPHA) for n in bgnames])
famscore_bg = {n: {} for n in bgnames}
labblock = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=a.encode(), sequence=seqs[a]).digitize(ALPHA) for a in seqs])
famscore_lab = {a: {} for a in seqs}
for hmm_obj, hits in zip(fam_hmms.values(), pyhmmer.hmmsearch(list(fam_hmms.values()), bgblock, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        tgt = h.name.decode() if isinstance(h.name, bytes) else h.name
        famscore_bg[tgt][hn] = max(famscore_bg[tgt].get(hn, 0.0), h.score)
for hmm_obj, hits in zip(fam_hmms.values(), pyhmmer.hmmsearch(list(fam_hmms.values()), labblock, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        tgt = h.name.decode() if isinstance(h.name, bytes) else h.name
        famscore_lab[tgt][hn] = max(famscore_lab[tgt].get(hn, 0.0), h.score)

# family T from training own-scores (round-0)
fam_T0 = {}
for f, members in fam_accs.items():
    fam_T0[f] = 0.8 * min(famscore_lab[a].get(f, 0.0) for a in members)

# self-expansion: bg hits >= T0, cap 60 per family
expanded_hmms = {}
exp_count = {}
for f, members in fam_accs.items():
    hits = sorted(((famscore_bg[n].get(f, 0.0), n) for n in bgnames), reverse=True)
    add = [n for s, n in hits if s >= fam_T0[f]][:60]
    exp_count[f] = len(add)
    expanded_hmms["FAMX:" + f] = build_hmm("FAMX:" + f, members + add, allseq)
print("family expansion counts:", exp_count)

# v3 sub-family HMMs: load from cached expansion (rebuild quickly using acc2clu logic is costly; reuse: sub-family = UniRef90 cluster from v3 run)
v3 = json.load(open(f"{RESULTS}/validation_v3.json"))
sub_members = {cid: v["train_members"] for cid, v in v3["subfamilies"].items()}
sub_hmms = {}
for cid, members in sub_members.items():
    sub_hmms[cid] = build_hmm(cid, [m for m in members if m in allseq], allseq)  # training members only (v3 expansions not in seqs)
pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "PTE_NEIGHBOR"]
gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "GENERIC_ESTERASE"]
neg_hmms = {"NEG:PLL": build_hmm("NEG:PLL", pll, allseq), "NEG:EST": build_hmm("NEG:EST", gest, allseq)}

# full corpus scoring
corpus = {}
for fn in os.listdir(CORPUS):
    if fn.endswith(".gz"): corpus.update(read_fasta_text(gzip.open(f"{CORPUS}/{fn}", "rt").read()))
    elif fn.endswith(".fasta"): corpus.update(read_fasta_text(open(f"{CORPUS}/{fn}").read()))
for acc, s in allseq.items(): corpus[f"LAB|{acc}"] = s
names = list(corpus)
block = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=corpus[n]).digitize(ALPHA) for n in names])
all_hmms = list(expanded_hmms.values()) + list(sub_hmms.values()) + list(neg_hmms.values())
scores = {n: defaultdict(float) for n in names}
for hmm_obj, hits in zip(all_hmms, pyhmmer.hmmsearch(all_hmms, block, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        tgt = h.name.decode() if isinstance(h.name, bytes) else h.name
        scores[tgt][hn] = max(scores[tgt][hn], h.score)
print("corpus:", len(names), "HMMs:", len(all_hmms))

def negmax(n): return max(scores[n].get("NEG:PLL", 0.0), scores[n].get("NEG:EST", 0.0))
# freeze thresholds: expanded family HMMs use family rule on training members; sub-family HMMs use same
T, D = {}, {}
prof_members = {("FAMX:" + f): m for f, m in fam_accs.items()}
prof_members.update(sub_members)
for pid, members in prof_members.items():
    own = [scores[f"LAB|{a}"].get(pid, 0.0) for a in members]
    T[pid] = 0.8 * min(own)
    deltas = sorted((scores[f"LAB|{a}"].get(pid, 0.0) - negmax(f"LAB|{a}") for a in split["negatives_tune"]), reverse=True)
    D[pid] = max((deltas[1] if len(deltas) > 1 else 0.0) + 10.0, 0.0)

def call(n):
    for pid in prof_members:
        s = scores[n].get(pid, 0.0)
        if s >= T[pid] and (s - negmax(n)) >= D[pid]: return pid
    return None

tr = [a for a in train if call(f"LAB|{a}")]
ho = [a for a in heldout if call(f"LAB|{a}")]
fpe = [a for a in split["negatives_eval"] if call(f"LAB|{a}")]
tp, fp2 = len(ho), len(fpe)
prec = tp / (tp + fp2) if tp + fp2 else 1.0
s_pairs = [(1, max(scores[f"LAB|{a}"].get(p, 0.0) for p in prof_members) - negmax(f"LAB|{a}")) for a in heldout] + \
          [(0, max(scores[f"LAB|{a}"].get(p, 0.0) for p in prof_members) - negmax(f"LAB|{a}")) for a in split["negatives_eval"]]
ps = sorted(s for l, s in s_pairs if l == 1); ns = [s for l, s in s_pairs if l == 0]
auc = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p in ps for n in ns) / (len(ps) * len(ns))
gate1 = len(ho) / len(heldout)
pc_pass = len(tr) == len(train)
print(f"POSITIVE CONTROL: {len(tr)}/{len(train)} {'PASS' if pc_pass else 'FAIL'} missing={[a for a in train if a not in tr]}")
print(f"GATE 1 recall: {len(ho)}/{len(heldout)} = {gate1:.3f} missed={[a for a in heldout if a not in ho]}")
print(f"GATE 2 precision: {prec:.3f} fp={[(a, negs[a]['neg_class']) for a in fpe]}; AUROC={auc:.3f}")
res = {"design": "v4 union of (family HMMs self-expanded via unlabeled bg hits >= family T, cap 60) + (v3 UniRef90 sub-family HMMs) + NEG HMMs; thresholds frozen by v2/v3 rules",
       "family_expansion_counts": exp_count,
       "positive_control": {"pass": pc_pass, "recovered": len(tr), "of": len(train), "missing": [a for a in train if a not in tr]},
       "gate1_recall": {"recall": round(gate1, 3), "hits": len(ho), "of": len(heldout), "pass": gate1 >= 0.90,
                         "missed": [a for a in heldout if a not in ho]},
       "gate2_specificity": {"precision": round(prec, 3), "fp_eval": fpe, "auroc": round(auc, 3), "pass": prec >= 0.90 and auc >= 0.90},
       "corpus_size": len(names), "split_version": 2}
json.dump(res, open(f"{RESULTS}/validation_v4.json", "w"), indent=2)
print("wrote results/validation_v4.json")
