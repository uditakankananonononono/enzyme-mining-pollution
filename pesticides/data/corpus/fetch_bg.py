import urllib.request, urllib.parse, time
UA = {"User-Agent": "enzyme-mining-research/1.0"}
queries = [
    ("bg_esterase.fasta", 'protein_name:esterase AND taxonomy_id:2 AND length:[200 TO 600]', 7000),
    ("bg_amidohydrolase.fasta", '(amidohydrolase OR lactonase OR deaminase) AND taxonomy_id:2 AND length:[200 TO 600]', 5000),
    ("bg_euk_ce.fasta", 'carboxylesterase AND taxonomy_id:2759 AND length:[300 TO 800]', 3000),
]
for fname, query, want in queries:
    cursor = None; got = 0
    with open(fname, "w") as out:
        while got < want:
            params = {"query": query, "format": "fasta", "size": 500}
            if cursor: params["cursor"] = cursor
            u = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
            req = urllib.request.Request(u, headers=UA)
            try:
                with urllib.request.urlopen(req, timeout=40) as r:
                    body = r.read().decode()
                    link = r.headers.get("Link", "")
            except Exception as ex:
                print(fname, "fetch err", str(ex)[:60]); time.sleep(3); continue
            n = body.count("\n>")
            out.write(body)
            got += n
            if 'rel="next"' not in link: break
            import re
            m = re.search(r'cursor=([^&>]+)', link)
            if not m: break
            cursor = urllib.parse.unquote(m.group(1))
            time.sleep(0.3)
    print(fname, "sequences:", got)
