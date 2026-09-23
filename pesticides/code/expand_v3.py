"""v3: UniRef90 homolog expansion of TRAINING-only seeds -> sub-family HMMs.
Leak guards: expansion excludes every accession in reference_set.csv + negative_set.csv.
Thresholds frozen by same rules as v2 (T=0.8*min training own-score; D=2nd-max tune-neg delta+10)."""
import csv, json, os, re, gzip, time, urllib.request, urllib.parse
from collections import defaultdict
import pyfamsa, pyhmmer
from pyhmmer.easel import Alphabet, TextMSA, DigitalSequenceBlock
from pyhmmer.plan7 import Builder, Background
ALPHA = Alphabet.amino(); BUILDER = Builder(ALPHA); BG = Background(ALPHA)
UA = {"User-Agent": "enzyme-mining-research/1.0"}
DATA, CORPUS, RESULTS = "../data", "../data/corpus", "../results"

def uget(url):
    req = urllib.request.Request(url, headers=UA)
    for i in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r: return r.read().decode()
        except Exception: time.sleep(2 + i)
    return None

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
LABELED = set(refs) | set(negs)
seqs = {}
for fn in os.listdir(f"{DATA}/fasta"):
    seqs[fn[:-6]] = list(read_fasta_text(open(f"{DATA}/fasta/{fn}").read()).values())[0]
split = json.load(open(f"{DATA}/split.json")); assert split["version"] == 2
train = split["train"]

# 1) map training accessions to UniRef90 clusters
acc2clu = {}
for acc in train:
    j = uget("https://rest.uniprot.org/uniref/search?" + urllib.parse.urlencode({"query": acc, "size": 10}))
    cid = None
    if j:
        try:
            for r in json.loads(j).get("results", []):
                if r.get("entryType") == "UniRef90": cid = r["id"]; break
                if not cid:
                    rm = r.get("representativeMember", {}).get("uniref90Id")
                    if rm: cid = rm
        except Exception: pass
    acc2clu[acc] = cid or f"SINGLETON:{acc}"
    time.sleep(0.25)
clusters = defaultdict(list)
for acc, cid in acc2clu.items(): clusters[cid].append(acc)
print("sub-families:", len(clusters))

# 2) fetch expansion members (fasta per cluster), excluding labeled
expansion = {}  # member_name -> seq ; record provenance
expand_of = defaultdict(list)
for cid, members in clusters.items():
    if cid.startswith("SINGLETON:"): continue
    body = uget("https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": f"(uniref_cluster_90:{cid})", "format": "fasta", "size": 40}))
    if not body: continue
    got = read_fasta_text(body)
    for hdr, s in got.items():
        # uniprot fasta header id looks like sp|P0A433|OPD_SPHFU ; extract accession
        m = re.match(r"(?:sp|tr)\|([^|]+)\|", hdr)
        acc = m.group(1) if m else hdr
        if acc in LABELED: continue
        key = f"EXP|{acc}"
        if key not in expansion:
            expansion[key] = s; expand_of[cid].append(key)
    time.sleep(0.3)
print("expansion seqs:", len(expansion))

# 3) build sub-family HMMs
allseq = dict(seqs); allseq.update(expansion)
def build_hmm(name, keys):
    fs = [pyfamsa.Sequence(k.encode(), allseq[k].encode()) for k in keys]
    if len(fs) >= 2:
        msa = pyfamsa.Aligner(threads=2).align(fs)
        tmsa = TextMSA(name=name.encode(), sequences=[pyhmmer.easel.TextSequence(name=s.id, sequence=s.sequence.decode()) for s in msa])
        hmm, _, _ = BUILDER.build_msa(tmsa.digitize(ALPHA), BG)
    else:
        hmm, _, _ = BUILDER.build(pyhmmer.easel.TextSequence(name=name.encode(), sequence=allseq[keys[0]]).digitize(ALPHA), BG)
    hmm.name = name.encode(); return hmm

pos_hmms, sub_members = {}, {}
for cid, members in clusters.items():
    keys = members + expand_of.get(cid, [])[:40]
    pos_hmms[cid] = build_hmm(cid, keys)
    sub_members[cid] = members
pll = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "PTE_NEIGHBOR"]
gest = [a for a in split["negatives_tune"] if negs[a]["neg_class"] == "GENERIC_ESTERASE"]
neg_hmms = {"NEG:PLL": build_hmm("NEG:PLL", pll), "NEG:EST": build_hmm("NEG:EST", gest)}

# 4) score corpus
corpus = {}
for fn in os.listdir(CORPUS):
    if fn.endswith(".gz"): corpus.update(read_fasta_text(gzip.open(f"{CORPUS}/{fn}", "rt").read()))
    elif fn.endswith(".fasta"): corpus.update(read_fasta_text(open(f"{CORPUS}/{fn}").read()))
for acc, s in allseq.items(): corpus[f"LAB|{acc}"] = s
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

def negmax(n): return max(scores[n].get("NEG:PLL", 0.0), scores[n].get("NEG:EST", 0.0))
T, D = {}, {}
for cid, members in sub_members.items():
    own = [scores[f"LAB|{a}"].get(cid, 0.0) for a in members]
    T[cid] = 0.8 * min(own)
    deltas = sorted((scores[f"LAB|{a}"].get(cid, 0.0) - negmax(f"LAB|{a}") for a in split["negatives_tune"]), reverse=True)
    D[cid] = max((deltas[1] if len(deltas) > 1 else 0.0) + 10.0, 0.0)

def call(n):
    best, bf = -1, None
    for cid in pos_hmms:
        s = scores[n].get(cid, 0.0)
        if s >= T[cid] and (s - negmax(n)) >= D[cid] and s > best: best, bf = s, cid
    return bf

tr = [a for a in train if call(f"LAB|{a}")]
ho = [a for a in split["heldout"] if call(f"LAB|{a}")]
fpe = [a for a in split["negatives_eval"] if call(f"LAB|{a}")]
tp, fp2 = len(ho), len(fpe)
prec = tp / (tp + fp2) if tp + fp2 else 1.0
s_pairs = [(1, max(scores[f"LAB|{a}"].get(c, 0.0) for c in pos_hmms) - negmax(f"LAB|{a}")) for a in split["heldout"]] + \
          [(0, max(scores[f"LAB|{a}"].get(c, 0.0) for c in pos_hmms) - negmax(f"LAB|{a}")) for a in split["negatives_eval"]]
ps = sorted(s for l, s in s_pairs if l == 1); ns = [s for l, s in s_pairs if l == 0]
auc = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p in ps for n in ns) / (len(ps) * len(ns))
gate1 = len(ho) / len(split["heldout"])
pc_pass = len(tr) == len(train)
print(f"POSITIVE CONTROL: {len(tr)}/{len(train)} {'PASS' if pc_pass else 'FAIL'} missing={[a for a in train if a not in tr]}")
print(f"GATE 1 recall: {len(ho)}/{len(split['heldout'])} = {gate1:.3f} missed={[a for a in split['heldout'] if a not in ho]}")
print(f"GATE 2 precision: {prec:.3f} fp={[(a, negs[a]['neg_class']) for a in fpe]}; AUROC={auc:.3f}")
res = {"design": "v3 UniRef90-expanded sub-family HMMs (training-only seeds; labeled accessions excluded from expansion); T=0.8*min(train own); D=2nd-max tune-neg delta+10",
       "subfamilies": {cid: {"train_members": m, "n_expansion": len(expand_of.get(cid, []))} for cid, m in sub_members.items()},
       "positive_control": {"pass": pc_pass, "recovered": len(tr), "of": len(train), "missing": [a for a in train if a not in tr]},
       "gate1_recall": {"recall": round(gate1, 3), "hits": len(ho), "of": len(split["heldout"]), "pass": gate1 >= 0.90,
                         "missed": [a for a in split["heldout"] if a not in ho]},
       "gate2_specificity": {"precision": round(prec, 3), "fp_eval": fpe, "auroc": round(auc, 3), "pass": prec >= 0.90 and auc >= 0.90},
       "corpus_size": len(names), "split_version": 2}
json.dump(res, open(f"{RESULTS}/validation_v3.json", "w"), indent=2)
print("wrote results/validation_v3.json")
