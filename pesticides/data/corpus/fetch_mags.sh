cd ~/enzyme-mining-pollution/pesticides/data/corpus/mgnify
for id in MGYG000535630 MGYG000535629 MGYG000535628 MGYG000535627 MGYG000535626 MGYG000535625 MGYG000535624 MGYG000535622; do
  curl -s --max-time 20 "https://www.ebi.ac.uk/metagenomics/api/v1/genomes/$id" | python3 -c "
import json,sys
try:
    d=json.load(sys.stdin)['data']
    a=d['attributes']; b=d['relationships']
    cat=b.get('catalogue',{}).get('data',{}).get('id','?') if b.get('catalogue',{}).get('data') else '?'
    bio=b.get('biome',{}).get('data',{}).get('id','?') if b.get('biome',{}).get('data') else '?'
    print(f\"{d['id']}\t{cat}\t{bio}\tlen={a.get('length','?')}\")
except Exception as e: print('ERR', e)
" >> provenance.tsv
  curl -s --max-time 60 -o "$id.faa" "https://www.ebi.ac.uk/metagenomics/api/v1/genomes/$id/downloads/$id.faa"
  echo "$id $(grep -c '^>' $id.faa 2>/dev/null) proteins"
done
