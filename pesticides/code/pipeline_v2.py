"""v2: per-family thresholds. T_f = 0.8*min(training own-family score);
D_f = 2nd-highest tune-negative delta + 10 (<=1 tune FP allowance, disclosed:
documented irreducible homology collision Pye3 vs YeiG-family esterases)."""
import csv, json, os, re, gzip
from collections import defaultdict
import pyfamsa, pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
ALPHA = Alphabet.amino(); BUILDER = Builder(ALPHA); BG = Background(ALPHA)
DATA, CORPUS, RESULTS = "../data", "../data/corpus", "../results"
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

refs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/reference_set.csv"))}
negs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/negative_set.csv"))}
seqs = {}
for fn in os.listdir(f"{DATA}/fasta"):
    seqs[fn[:-6]] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
split = json.load(open(f"{DATA}/split.json"))
assert split["version"] == 2

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
for acc in split["train"]:
    fam_accs[f"{refs[acc]['class']}:{refs[acc]['family']}"].append(acc)
pos_hmms = {k: build_hmm(k, a) for k, a in fam_accs.items()}
pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "PTE_NEIGHBOR"]
gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "GENERIC_ESTERASE"]
neg_hmms = {"NEG:PLL": build_hmm("NEG:PLL", pll), "NEG:EST": build_hmm("NEG:EST", gest)}

corpus = {}
for fn in os.listdir(CORPUS):
    if fn.endswith(".gz"): corpus.update(read_fasta_text(gzip.open(f"{CORPUS}/{fn}", "rt").read()))
    elif fn.endswith(".fasta"): corpus.update(read_fasta_text(open(f"{CORPUS}/{fn}").read()))
for acc, s in seqs.items(): corpus[f"LAB|{acc}"] = s
names = list(corpus)
block = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=corpus[n]).digitize(ALPHA) for n in names])
all_hmms = list(pos_hmms.values()) + list(neg_hmms.values())
scores = {n: defaultdict(float) for n in names}
for hmm_obj, hits in zip(all_hmms, pyhmmer.hmmsearch(all_hmms, block, cpus=2, E=100.0)):
    hn = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        tgt = h.name.decode() if isinstance(h.name, bytes) else h.name
        scores[tgt][hn] = max(scores[tgt][hn], h.score)
print("corpus:", len(names), "| HMMs:", len(all_hmms))

def famscore(n, f): return scores[n].get(f, 0.0)
def negmax(n): return max(scores[n].get("NEG:PLL", 0.0), scores[n].get("NEG:EST", 0.0))

# freeze per-family T and D on training + tune negatives
fams = list(pos_hmms)
T, D = {}, {}
for f in fams:
    members = fam_accs[f]
    own = [famscore(f"LAB|{a}", f) for a in members]
    T[f] = 0.8 * min(own)
    deltas = sorted((famscore(f"LAB|{a}", f) - negmax(f"LAB|{a}") for a in split["negatives_tune"]), reverse=True)
    D[f] = (deltas[1] if len(deltas) > 1 else 0.0) + 10.0
    D[f] = max(D[f], 0.0)
print("frozen T:", {k: round(v,1) for k,v in T.items()})
print("frozen D:", {k: round(v,1) for k,v in D.items()})

def call(n):
    best, bf = -1, None
    for f in fams:
        s = famscore(n, f)
        if s >= T[f] and (s - negmax(n)) >= D[f] and s > best: best, bf = s, f
    return bf

def evaluate():
    train_hit = [a for a in split["train"] if call(f"LAB|{a}")]
    hold_hit = [a for a in split["heldout"] if call(f"LAB|{a}")]
    fp_tune = [a for a in split["negatives_tune"] if call(f"LAB|{a}")]
    fp_eval = [a for a in split["negatives_eval"] if call(f"LAB|{a}")]
    tp, fp2 = len(hold_hit), len(fp_eval)
    prec = tp / (tp + fp2) if tp + fp2 else 1.0
    scored = [(1, max(famscore(f"LAB|{a}", f) for f in fams) - negmax(f"LAB|{a}")) for a in split["heldout"]] + \
             [(0, max(famscore(f"LAB|{a}", f) for f in fams) - negmax(f"LAB|{a}")) for a in split["negatives_eval"]]
    ps = sorted(s for l, s in scored if l == 1); ns = [s for l, s in scored if l == 0]
    auc = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p in ps for n in ns) / (len(ps) * len(ns))
    return train_hit, hold_hit, fp_tune, fp_eval, prec, auc

tr, ho, fpt, fpe, prec, auc = evaluate()
pc_pass = len(tr) == len(split["train"])
gate1 = len(ho) / len(split["heldout"])
print(f"POSITIVE CONTROL: {len(tr)}/{len(split['train'])} -> {'PASS' if pc_pass else 'FAIL'} (missing: {[a for a in split['train'] if a not in tr]})")
print(f"GATE 1 recall: {len(ho)}/{len(split['heldout'])} = {gate1:.3f} (need >=0.90) missed: {[a for a in split['heldout'] if a not in ho]}")
print(f"GATE 2 precision: {prec:.3f} (fp={fpe} {[ (a, negs[a]['neg_class']) for a in fpe]}); AUROC={auc:.3f}")
res = {"design": "v2 per-family T=0.8*min(train own), D=2nd-max tune-neg delta + 10 (<=1 tune FP allowance, disclosed)",
       "frozen_T": {k: round(v, 2) for k, v in T.items()}, "frozen_D": {k: round(v, 2) for k, v in D.items()},
       "positive_control": {"pass": pc_pass, "recovered": len(tr), "of": len(split["train"]), "missing": [a for a in split["train"] if a not in tr]},
       "gate1_recall": {"recall": round(gate1, 3), "hits": len(ho), "of": len(split["heldout"]), "pass": gate1 >= 0.90,
                         "missed": [a for a in split["heldout"] if a not in ho]},
       "gate2_specificity": {"precision": round(prec, 3), "fp_eval": fpe, "auroc": round(auc, 3), "pass": prec >= 0.90 and auc >= 0.90,
                              "fp_detail": [(a, negs[a]["neg_class"]) for a in fpe],
                              "tune_fp_disclosed": fpt},
       "corpus_size": len(names), "split_version": 2}
json.dump(res, open(f"{RESULTS}/validation_v2.json", "w"), indent=2)
print("wrote results/validation_v2.json")
