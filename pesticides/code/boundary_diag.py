import csv, json, os
from collections import Counter
def read_fasta_text(txt):
    seqs, name, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs
seqs = {}
for fn in os.listdir("../data/fasta"):
    seqs[fn[:-6]] = list(read_fasta_text(open(f"../data/fasta/{fn}").read()).values())[0]
split = json.load(open("../data/split.json"))
v4 = json.load(open("../results/validation_v4.json"))
refs = {r["accession"]: r for r in csv.DictReader(open("../data/reference_set.csv"))}

def kmers(s, k=3): return set(s[i:i+k] for i in range(len(s)-k+1))
def jac(a, b):
    A, B = kmers(a), kmers(b)
    return len(A & B) / len(A | B) if A | B else 0.0

train_k = {a: kmers(seqs[a]) for a in split["train"]}
def max_jac(acc):
    A = kmers(seqs[acc])
    best, ba = 0.0, None
    for t, T in train_k.items():
        j = len(A & T) / len(A | T)
        if j > best: best, ba = j, t
    return best, ba

print("=== held-out 3mer-Jaccard vs nearest training member (missed vs recovered) ===")
rows = []
for a in split["heldout"]:
    j, nb = max_jac(a)
    status = "MISSED" if a in v4["gate1_recall"]["missed"] else "recovered"
    rows.append((a, refs[a]["family"], status, round(j, 3), nb, refs[a]["name"][:40]))
    print(f"  {a:12s} {refs[a]['family']:10s} {status:9s} jacc={j:.3f} nearest={nb} ({refs[a]['name'][:35]})")
json.dump({"jaccard_vs_train": [{"acc": a, "family": f, "status": s, "jaccard": j, "nearest_train": n}
                                 for a, f, s, j, n, nm in rows]}, open("../results/boundary_diagnosis.json", "w"), indent=2)
misses = [r for r in rows if r[2] == "MISSED"]
hits = [r for r in rows if r[2] == "recovered"]
print(f"\nmissed n={len(misses)} max-jacc mean={sum(r[3] for r in misses)/len(misses):.3f}")
print(f"recovered n={len(hits)} max-jacc mean={sum(r[3] for r in hits)/len(hits):.3f}")
