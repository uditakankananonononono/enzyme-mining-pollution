# Sequence-Validated Mining of Pesticide-Degrading Enzymes from Metagenomes: a Specificity-Gated Platform and the Fold-Coverage Boundary of Profile Recall
(working title; builder 20 slice, project enzyme-mining-pollution)

## Abstract
[TODO after mining] We built a profile-HMM mining platform for organophosphate (OP) and pyrethroid hydrolases with gates locked before outcomes: >=90% recall on held-out characterized enzymes, and a specificity discriminator validated against documented false-positive modes (phosphotriesterase-like lactonases, esterase-positive/pyrethroid-inactive enzymes). The discriminator passes robustly (precision 1.0, AUROC 0.91, four architectures). Canon-wide recall fails (0.69) for a quantified structural reason: 4/13 held-out enzymes are fold-orphans (3-mer Jaccard 0.05-0.08 to nearest training enzyme vs >=0.28 for all recoveries); per-fold recall is 9/9. We mine 8 barley-rhizosphere MAG proteomes (23,989 proteins) + ~15.5k UniProtKB background sequences and deliver a vetted candidate DB with explicit coverage caveats.

## 1. Introduction
- OP + pyrethroid contamination; enzymatic bioremediation need.
- Canon: PTE (P0A433/P0A434/OpdA), MPH family, PONs, OPAA, DFPase; pyrethroid: PytH/PytZ/EstP/Sys410/Pye3/Est3385/Ces2a/2e, insect resistance esterases.
- Documented FP modes: PLLs (Afriat 2006; SacPox/SsoPox), screen-level esterase overcalling (Sys410 3:1, Pye3 6:1). The platform's contribution: a specificity gate validated against those modes.
- Prior art: targeted metagenome mining (Robinson et al. 2023 PMC10781932), PAZy (plastics), PCycDB (P cycling); single-enzyme functional screens. Gap: no pesticide-enzyme mining platform with locked recall/specificity gates.

## 2. Methods
2.1 Reference canon: 76 enzymes (39 OP / 37 PYR), tiered VERIFIED (47)/THIN (29), each with evidence URL; NCBI-only entries (Sys410, Est3385, Pye3) included; missing canon honestly reported (Est804, Cest2923, EstSt7, PytY).
2.2 Negative set: 61 labeled negatives across two documented modes.
2.3 Gates locked pre-outcome (GATES_LOCKED.md, commit 62b8c87): recall >=0.90 held-out; discriminator precision >=0.90 + AUROC >=0.90; novelty <80% identity + curation absence.
2.4 Pipeline: pyfamsa MSAs -> pyhmmer family + UniRef90-expanded sub-family HMMs + jackhmmer-style self-expansion on unlabeled background; negative HMMs (PLL, generic esterases); thresholds frozen from training + tune negatives only (T = 0.8 x min training own-score; D = 2nd-max tune-neg delta + 10 bits); split v2 (strata <3 members all-training).
2.5 Corpus: 24,240-sequence validation corpus (E. coli K-12, B. subtilis decoys, UniProt background); mining corpus adds 8 barley-rhizosphere MAG proteomes (MGnify; 23,989 proteins).
2.6 Compute: 2-core/2GB sandbox; all sources public APIs (UniProt REST, NCBI eutils, MGnify API).

## 3. Results
3.1 Positive control: 34/34 training enzymes recovered (v2-v4).
3.2 Specificity discriminator: precision 1.0, AUROC 0.906-0.916 across 4 designs (Fig 1, 3). Documented intrinsic collision quantified: YeiG B6VG94 = Pye3 ACJ07038.1 at 281.7 bits (Fig 5) - sequence-only separation impossible; resolved by negative-HMM margin.
3.3 Recall boundary: canon-wide recall 0.692 (9/13); per-fold recall 9/9 among held-out enzymes with a training relative (Jaccard >=0.277); all misses are fold-orphans (Jaccard 0.051-0.081) (Fig 2). Robust across v1-v4 (Table 1).
3.4 Mining: [TODO] candidate counts by class/source; novelty (jacc<0.5) counts; coverage caveat.

## 4. Discussion
- The amended-bar claim: specificity-gated mining is validated; recall is bounded by canon fold coverage - a boundary that quantifies where profile mining cannot see.
- Convergent evolution of pesticide hydrolases across folds explains the boundary (Sys410 family V vs PytH alpha/beta fold etc.).
- Limitations: no wet-lab validation (candidates are computational); structure work (colabfold) reserved for top candidates.

## 5. Conclusion
[TODO]

## References (key)
Afriat et al. 2006 (PLL promiscuity); Horne et al. 2002 (opdA); Singh 2009 Nat Rev Microbiol; Scott et al. 2008 (enzymatic basis); Wang et al. (PytH JZ-1); Wu et al. (EstP Klebsiella); Sys410 (10.1186/1475-2859-11-33); Pye3 (10.1186/1475-2859-7-38); Est3385 (AND66123.1); Robinson et al. 2023 (PMC10781932); Buchholz et al. 2022 (PAZy); Zeng et al. 2022 (PCycDB).
