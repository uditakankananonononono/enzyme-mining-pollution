"""Assemble reference_set.csv + negative_set.csv for industrial/ slice from raw UniProt pulls.
Tiers per B20 convention: VERIFIED = reviewed w/ target activity evidence; THIN = annotation-only (unreviewed)."""
import csv, random, re, os
D = "/home/sandbox/industrial/data"
random.seed(20260923)

def load(fn):
    rows = list(csv.DictReader(open(f"{D}/{fn}"), delimiter="\t"))
    return [r for r in rows if r.get("Entry")]

def fam_hld(r):
    return "HALOALKANE" if "3.8.1.5" in (r.get("EC number") or "") else "HALOACID"

def fam_rhd(r):
    s = ((r.get("Protein names") or "") + " " + (r.get("Gene Names") or "")).lower()
    if re.search(r"naphthalene|phenanthrene|nahac|ndob|nida|phda|phnac|nagac|paah", s): return "RHD_PAH"
    if re.search(r"biphenyl|bpha", s): return "RHD_BPH"
    if re.search(r"toluene|benzene|todc|bnza|xylx|bedc", s): return "RHD_BTX"
    return "RHD_OTHER"

out = []
# HLD VERIFIED
for r in load("raw_haloalkane_reviewed.tsv") + load("raw_haloacid_reviewed.tsv"):
    out.append(dict(accession=r["Entry"], cls="HLD", family=fam_hld(r), source_db="uniprot",
        name=(r["Protein names"] or "").split("(")[0].strip(), organism=r["Organism"][:60],
        length=r["Length"], reviewed="True", tier="VERIFIED",
        evidence_url=f"https://www.uniprot.org/uniprotkb/{r['Entry']}/entry",
        evidence_note="Swiss-Prot reviewed, EC "+(r.get("EC number") or ""), func_snippet=(r["Protein names"] or "")[:160]))
# RHD VERIFIED
for r in load("raw_rhd_reviewed.tsv"):
    out.append(dict(accession=r["Entry"], cls="RHD", family=fam_rhd(r), source_db="uniprot",
        name=(r["Protein names"] or "").split("(")[0].strip(), organism=r["Organism"][:60],
        length=r["Length"], reviewed="True", tier="VERIFIED",
        evidence_url=f"https://www.uniprot.org/uniprotkb/{r['Entry']}/entry",
        evidence_note="Swiss-Prot reviewed, EC "+(r.get("EC number") or ""), func_snippet=(r["Protein names"] or "")[:160]))
# dedupe by accession
seen, refs = set(), []
for r in out:
    if r["accession"] not in seen: seen.add(r["accession"]); refs.append(r)
# THIN: seeded sample 12 per class
hld_thin = [r for r in load("raw_hld_thin.tsv")]
rhd_thin = [r for r in load("raw_rhd_thin.tsv")]
random.shuffle(hld_thin); random.shuffle(rhd_thin)
for r in hld_thin[:12]:
    refs.append(dict(accession=r["Entry"], cls="HLD", family=fam_hld(r), source_db="uniprot",
        name=(r["Protein names"] or "").split("(")[0].strip(), organism=r["Organism"][:60],
        length=r["Length"], reviewed="False", tier="THIN",
        evidence_url=f"https://www.uniprot.org/uniprotkb/{r['Entry']}/entry",
        evidence_note="TrEMBL annotation-only propagation", func_snippet=(r["Protein names"] or "")[:160]))
for r in rhd_thin[:12]:
    refs.append(dict(accession=r["Entry"], cls="RHD", family=fam_rhd(r), source_db="uniprot",
        name=(r["Protein names"] or "").split("(")[0].strip(), organism=r["Organism"][:60],
        length=r["Length"], reviewed="False", tier="THIN",
        evidence_url=f"https://www.uniprot.org/uniprotkb/{r['Entry']}/entry",
        evidence_note="TrEMBL annotation-only propagation", func_snippet=(r["Protein names"] or "")[:160]))

fn = f"{D}/reference_set.csv"
cols = ["accession","cls","family","source_db","name","organism","length","reviewed","tier","evidence_url","evidence_note","func_snippet"]
w = csv.DictWriter(open(fn,"w",newline=""), fieldnames=cols); w.writeheader(); w.writerows(refs)

# negatives
negs = []
def add_negs(fn_, neg_class, n, lenlo=100, lenhi=700, note=""):
    rows = [r for r in load(fn_) if r.get("Length") and lenlo <= int(r["Length"]) <= lenhi]
    random.shuffle(rows)
    for r in rows[:n]:
        negs.append(dict(accession=r["Entry"], neg_class=neg_class,
            name=(r["Protein names"] or "").split("(")[0].strip(), organism=r["Organism"][:60],
            length=r["Length"], evidence_url=f"https://www.uniprot.org/uniprotkb/{r['Entry']}/entry", note=note))
add_negs("raw_neg_had.tsv", "HAD_PHOSPHATASE", 20, note="HAD-fold hydrolase, actual phosphatase/nucleotidase activity")
add_negs("raw_neg_eh.tsv", "EPOXIDE_HYDROLASE", 15, note="alpha/beta hydrolase fold, no dehalogenation")
add_negs("raw_neg_rieske_fd.tsv", "RIESKE_FERREDOXIN", 20, note="Rieske [2Fe-2S] electron carrier, no oxygenase")
add_negs("raw_neg_ringcleavage.tsv", "RING_CLEAVAGE_DO", 15, note="intra/extradiol ring-cleavage dioxygenase, not ring-hydroxylating")
ncols = ["accession","neg_class","name","organism","length","evidence_url","note"]
w = csv.DictWriter(open(f"{D}/negative_set.csv","w",newline=""), fieldnames=ncols); w.writeheader(); w.writerows(negs)

from collections import Counter
print("refs:", len(refs), Counter((r["cls"], r["tier"]) for r in refs))
print("families:", Counter((r["cls"], r["family"]) for r in refs))
print("negs:", len(negs), Counter(n["neg_class"] for n in negs))
