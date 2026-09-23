"""Pesticide-enzyme mining pipeline (builder 20 slice).
Stage 1: family profile HMMs (training VERIFIED references + negative families).
Stage 2: candidate scoring over corpus; thresholds frozen on training+tune-negatives only.
Stage 3: positive control (training recovery), held-out recall (GATE 1), discriminator (GATE 2).
"""
import csv, json, os, re, gzip, sys
from collections import defaultdict
import pyfamsa
import pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, SequenceFile, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, CORPUS, RESULTS = f"{ROOT}/data", f"{ROOT}/data/corpus", f"{ROOT}/results"
os.makedirs(RESULTS, exist_ok=True)
ALPHA = Alphabet.amino()
BUILDER = Builder(ALPHA)
BG = Background(ALPHA)

def read_fasta_text(txt):
    seqs, name, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def load_labeled():
    refs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/reference_set.csv"))}
    negs = {r["accession"]: r for r in csv.DictReader(open(f"{DATA}/negative_set.csv"))}
    seqs = {}
    for fn in os.listdir(f"{DATA}/fasta"):
        acc = fn[:-6]
        seqs[acc] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
    split = json.load(open(f"{DATA}/split.json"))
    return refs, negs, seqs, split

def build_hmm(name, accs, seqs):
    fs = [pyfamsa.Sequence(acc.encode(), seqs[acc].encode()) for acc in accs if acc in seqs]
    if len(fs) >= 2:
        msa = pyfamsa.Aligner(threads=2).align(fs)
        tmsa = TextMSA(name=name.encode(),
                       sequences=[pyhmmer.easel.TextSequence(name=s.id, sequence=s.sequence.decode()) for s in msa])
        hmm, _, _ = BUILDER.build_msa(tmsa.digitize(ALPHA), BG)
    else:
        s = pyhmmer.easel.TextSequence(name=name.encode(), sequence=seqs[accs[0]])
        hmm, _, _ = BUILDER.build(s.digitize(ALPHA), BG)
    hmm.name = name.encode()
    return hmm

def main():
    refs, negs, seqs, split = load_labeled()
    fam_accs = defaultdict(list)
    for acc in split["train"]:
        fam_accs[(refs[acc]["class"], refs[acc]["family"])].append(acc)
    pos_hmms = {f"{c}:{f}": build_hmm(f"{c}:{f}", a, seqs) for (c, f), a in fam_accs.items()}
    pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "PTE_NEIGHBOR"]
    gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "GENERIC_ESTERASE"]
    neg_hmms = {"NEG:PLL": build_hmm("NEG:PLL", pll, seqs), "NEG:EST": build_hmm("NEG:EST", gest, seqs)}
    print("positive HMMs:", list(pos_hmms), "| negative HMMs:", list(neg_hmms))

    # corpus: decoys + background + ALL labeled seqs (references + negatives embedded)
    corpus = {}
    for fn in os.listdir(CORPUS):
        if fn.endswith(".gz"): corpus.update(read_fasta_text(gzip.open(f"{CORPUS}/{fn}", "rt").read()))
        elif fn.endswith(".fasta"): corpus.update(read_fasta_text(open(f"{CORPUS}/{fn}").read()))
    for acc, s in seqs.items():
        corpus[f"LAB|{acc}"] = s
    print("corpus sequences:", len(corpus))

    names = list(corpus)
    block = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=corpus[n]).digitize(ALPHA) for n in names])
    all_hmms = list(pos_hmms.values()) + list(neg_hmms.values())
    scores = {n: defaultdict(float) for n in names}
    for hmm_obj, hits in zip(all_hmms, pyhmmer.hmmsearch(all_hmms, block, cpus=2, E=100.0)):
        hname = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
        for h in hits:
            hn = h.name.decode() if isinstance(h.name, bytes) else h.name
            scores[hn][hname] = max(scores[hn][hname], h.score)
    print("scored.")

    GXSXG = re.compile(r"G.S.G")
    def features(n):
        sc = scores[n]
        op = max((sc[k] for k in pos_hmms if k.startswith("OP:")), default=0.0)
        pyr = max((sc[k] for k in pos_hmms if k.startswith("PYR:")), default=0.0)
        neg = max(sc["NEG:PLL"], sc["NEG:EST"])
        return op, pyr, neg, bool(GXSXG.search(corpus[n]))

    # --- threshold freeze on training VERIFIED + tune negatives ONLY ---
    train = split["train"]; tune_neg = split["negatives_tune"]
    def candidates(T, D, subset):
        out = []
        for n in subset:
            op, pyr, neg, g = features(n)
            cls = "OP" if op >= pyr else "PYR"; best = max(op, pyr)
            if best >= T and (best - neg) >= D: out.append((n, cls, best, best - neg))
        return out
    grid_T = [0, 20, 40, 60, 80, 100, 125, 150, 200, 250, 300]
    grid_D = [0, 5, 10, 20, 40, 60]
    lab_train = [f"LAB|{a}" for a in train]; lab_tune = [f"LAB|{a}" for a in tune_neg]
    best_cfg = None
    for T in grid_T:
        for D in grid_D:
            rec_hits = {n for n, _, _, _ in candidates(T, D, lab_train)}
            fp = len(candidates(T, D, lab_tune))
            rec = len(rec_hits) / len(lab_train)
            prec = len(rec_hits) / (len(rec_hits) + fp) if (len(rec_hits) or fp) else 1.0
            if prec >= 0.98 and (best_cfg is None or (rec, -T, -D) > (best_cfg[0], -best_cfg[1], -best_cfg[2])):
                best_cfg = (rec, T, D, fp)
    rec_train, T, D, fp = best_cfg
    print(f"frozen: T={T} D={D} (training recall={rec_train:.3f}, tune FP={fp})")

    # --- positive control ---
    pc = {n for n, _, _, _ in candidates(T, D, lab_train)}
    pc_ok = len(pc) == len(lab_train)
    print(f"POSITIVE CONTROL: {len(pc)}/{len(lab_train)} training recovered -> {'PASS' if pc_ok else 'FAIL'}")

    # --- GATE 1: held-out recall ---
    lab_hold = [f"LAB|{a}" for a in split["heldout"]]
    hold_hits = {n for n, _, _, _ in candidates(T, D, lab_hold)}
    gate1 = len(hold_hits) / len(lab_hold)
    print(f"GATE 1 recall: {len(hold_hits)}/{len(lab_hold)} = {gate1:.3f} (need >=0.90)")

    # --- GATE 2: discriminator vs eval negatives ---
    lab_eval = [f"LAB|{a}" for a in split["negatives_eval"]]
    fp_eval = candidates(T, D, lab_eval)
    tp = len(hold_hits); fp2 = len(fp_eval)
    prec = tp / (tp + fp2) if tp + fp2 else 1.0
    # AUROC over heldout positives + eval negatives, score = best - neg
    pairs = [(1, features(n)[2] and None) for n in []]  # placeholder removed
    scored = []
    for n in lab_hold:
        op, pyr, neg, _ = features(n); scored.append((1, max(op, pyr) - neg))
    for n in lab_eval:
        op, pyr, neg, _ = features(n); scored.append((0, max(op, pyr) - neg))
    pos_s = sorted(s for l, s in scored if l == 1)
    neg_s = [s for l, s in scored if l == 0]
    auc = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p in pos_s for n in neg_s) / (len(pos_s) * len(neg_s))
    print(f"GATE 2 precision: {tp}/{tp}+{fp2} = {prec:.3f} (need >=0.90); AUROC={auc:.3f} (need >=0.90)")
    fp_detail = [(n.split("|")[1], negs[n.split('|')[1]]['neg_class'], round(b,1), round(d,1)) for n, _, b, d in fp_eval]

    res = {"frozen": {"T": T, "D": D, "train_recall": rec_train, "tune_fp": fp},
           "positive_control": {"pass": pc_ok, "recovered": len(pc), "of": len(lab_train)},
           "gate1_recall": {"recall": gate1, "hits": len(hold_hits), "of": len(lab_hold),
                            "pass": gate1 >= 0.90, "missed": [n.split("|")[1] for n in lab_hold if n not in hold_hits]},
           "gate2_specificity": {"precision": prec, "fp": fp2, "auroc": auc,
                                  "pass": prec >= 0.90 and auc >= 0.90, "fp_detail": fp_detail},
           "corpus_size": len(names)}
    json.dump(res, open(f"{RESULTS}/validation.json", "w"), indent=2)
    print("wrote results/validation.json")

if __name__ == "__main__":
    main()
