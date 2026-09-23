import csv, os, subprocess, sys
D = "/home/sandbox/industrial/data"
accs = [r["accession"] for r in csv.DictReader(open(f"{D}/reference_set.csv"))]
accs += [r["accession"] for r in csv.DictReader(open(f"{D}/negative_set.csv"))]
accs = sorted(set(accs)); print("total accessions:", len(accs))
os.makedirs(f"{D}/fasta", exist_ok=True)
for i in range(0, len(accs), 40):
    chunk = accs[i:i+40]
    q = " OR ".join(chunk)
    url = "https://rest.uniprot.org/uniprotkb/stream"
    out = subprocess.run(["curl","-s","-m","120","--get",url,"--data-urlencode",f"query=accession:({q})",
        "--data-urlencode","format=fasta"], capture_output=True, text=True)
    cur, buf = None, []
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            if cur: open(f"{D}/fasta/{cur}.fasta","w").write("\n".join(buf)+"\n")
            parts = line[1:].split("|"); cur = parts[1] if len(parts)>1 else line[1:].split()[0]; buf=[line]
        else: buf.append(line)
    if cur: open(f"{D}/fasta/{cur}.fasta","w").write("\n".join(buf)+"\n")
    print(f"chunk {i//40}: got so far {len(os.listdir(D+'/fasta'))}")
missing = [a for a in accs if not os.path.exists(f"{D}/fasta/{a}.fasta")]
print("missing:", missing)
