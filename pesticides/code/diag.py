import sys, os, csv, json, gzip
from collections import defaultdict
sys.path.insert(0, ".")
import importlib.util
spec = importlib.util.spec_from_file_location("pipe", "pipeline.py")
# import pieces without running main
import pyfamsa, pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
ALPHA = Alphabet.amino(); BUILDER = Builder(ALPHA); BG = Background(ALPHA)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath("pipeline.py")))
DATA = "../data"

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
    acc = fn[:-6]
    seqs[acc] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
split = json.load(open(f"{DATA}/split.json"))

def build_hmm(name, accs):
    fs = [pyfamsa.Sequence(a.encode(), seqs[a].encode()) for a in accs]
    if len(fs) >= 2:
        msa = pyfamsa.Aligner(threads=2).align(fs)
        tmsa = TextMSA(name=name.encode(), sequences=[pyhmmer.easel.TextSequence(name=s.id, sequence=s.sequence.decode()) for s in msa])
        hmm, _, _ = BUILDER.build_msa(tmsa.digitize(ALPHA), BG)
    else:
        s = pyhmmer.easel.TextSequence(name=name.encode(), sequence=seqs[accs[0]])
        hmm, _, _ = BUILDER.build(s.digitize(ALPHA), BG)
    hmm.name = name.encode(); return hmm

fam_accs = defaultdict(list)
for acc in split["train"]:
    fam_accs[(refs[acc]["class"], refs[acc]["family"])].append(acc)
print("train strata sizes:", {f"{c}:{f}": len(a) for (c,f),a in fam_accs.items()})
print("heldout:", [(a, refs[a]['class'], refs[a]['family']) for a in split['heldout']])
hmms = {f"{c}:{f}": build_hmm(f"{c}:{f}", a) for (c,f),a in fam_accs.items()}
pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"]=="PTE_NEIGHBOR"]
gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"]=="GENERIC_ESTERASE"]
hmms["NEG:PLL"] = build_hmm("NEG:PLL", pll); hmms["NEG:EST"] = build_hmm("NEG:EST", gest)

allacc = split["train"] + split["heldout"] + split["negatives_tune"] + split["negatives_eval"]
block = DigitalSequenceBlock(ALPHA, [pyhmmer.easel.TextSequence(name=a.encode(), sequence=seqs[a]).digitize(ALPHA) for a in allacc])
sc = {a: {} for a in allacc}
for hmm_obj, hits in zip(hmms.values(), pyhmmer.hmmsearch(list(hmms.values()), block, cpus=2, E=100.0)):
    hname = hmm_obj.name.decode() if isinstance(hmm_obj.name, bytes) else hmm_obj.name
    for h in hits:
        hn = h.name.decode() if isinstance(h.name, bytes) else h.name
        sc[hn][hname] = max(sc[hn].get(hname, 0.0), h.score)

print("\n== training members: own-family score, class max, neg max ==")
for a in split["train"]:
    famkey = f"{refs[a]['class']}:{refs[a]['family']}"
    own = sc[a].get(famkey, 0.0)
    cls = refs[a]['class']
    clsmax = max((sc[a].get(k,0.0) for k in hmms if k.startswith(cls+":")), default=0.0)
    negmax = max(sc[a].get("NEG:PLL",0.0), sc[a].get("NEG:EST",0.0))
    flag = "LOW" if own < 300 else ""
    print(f"  {a:12s} {famkey:16s} own={own:7.1f} clsmax={clsmax:7.1f} negmax={negmax:6.1f} {flag}")
print("\n== tune negatives: max positive-class score ==")
for a in split["negatives_tune"]:
    opmax = max((sc[a].get(k,0.0) for k in hmms if k.startswith("OP:")), default=0.0)
    pyrmax = max((sc[a].get(k,0.0) for k in hmms if k.startswith("PYR:")), default=0.0)
    negmax = max(sc[a].get("NEG:PLL",0.0), sc[a].get("NEG:EST",0.0))
    if max(opmax,pyrmax) > 40: print(f"  {a:12s} {negs[a]['neg_class']:18s} op={opmax:6.1f} pyr={pyrmax:6.1f} neg={negmax:6.1f}")
