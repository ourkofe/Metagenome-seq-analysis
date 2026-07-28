# Metagenome 분석 실습 기록 (16S rRNA / QIIME2)

2026 바이오 빅데이터(유전체) 분석 취업희망자 과정 – Metagenome Day 실습 복기 노트.
결과 이미지는 `figures/` 폴더에 캡쳐해서 넣고, 아래 해당 섹션에 링크.

> ⚠️ 원본 FASTQ, SILVA/Greengenes2 DB, 대용량 `.qza`는 이 레포에 올리지 않습니다. (`.gitignore` 참고)
> ⚠️ 아래 커맨드의 서버 경로는 실제 계정 경로 대신 `<WORKDIR>`로 표기했습니다.

## 실습 개요
- 샘플: Mouse fecal 16개 (Control 8 / Diabetes 8) — 제2형 당뇨병과 장내 미생물 연관성 연구
- Target region: 16S rDNA V3–V4 (Forward 341F `CCTACGGGNGGCWGCAG` / Reverse 806R `GACTACHVGGGTATCTAATCC`)
- Run: MiSeq 301 PE
- 환경: QIIME2 `quay.io/qiime2/amplicon:2024.10` (Docker)

## 진행 상황
- [x] 전처리 (import → primer trimming)
- [x] ASV 구축 (DADA2 → Rarefy)
- [x] Taxonomic assignment (Greengenes2 / SILVA)
- [x] Phylogenetic tree
- [x] Alpha diversity
- [ ] Beta diversity ← **내일 진행**
- [ ] ANCOM-BC
- [ ] Shotgun metagenome (MAGs)

---

## 1. 전처리 (Preprocessing)

```bash
# Paired-end fastq import (filepath manifest 사용)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path filepath \
  --output-path Raw.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Primer trimming (cutadapt)
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences Raw.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences Trim1_Primer.qza

qiime demux summarize \
  --i-data Trim1_Primer.qza --p-n 50000 \
  --o-visualization Trim1_Primer.qzv
```

**결과**: Q20 기준 truncate 지점 → Forward ~250bp, Reverse ~220bp (overlap ≈43bp, 필요 최소 12bp 충족)

![Trim1_Primer](figures/trim1_primer.png)

---

## 2. ASV 구축 (DADA2 → Rarefy)

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs Trim1_Primer.qza \
  --p-trunc-len-f 250 --p-trunc-len-r 220 \
  --p-max-ee-f 2 --p-max-ee-r 4 \
  --p-n-reads-learn 100000 \
  --p-chimera-method pooled \
  --o-table ASVtable.qza \
  --o-representative-sequences ASVrep-seqs.qza \
  --o-denoising-stats ASVstats.qza

# Rarefy (최소 샘플 depth로 정규화)
qiime feature-table rarefy \
  --i-table ASVtable.qza \
  --p-sampling-depth 25347 \
  --o-rarefied-table ASVtable_rarefied.qza
```

**결과**: 840 ASVs → rarefy 후 834 ASVs / 16 샘플 / min depth 25,347 reads (남은 reads 비율 50% 이상 확인)

![ASV table](figures/asv_table_summary.png)

---

## 3. Taxonomic Assignment

```bash
# A. classify-sklearn (Greengenes2, 사전학습 분류기)
qiime feature-classifier classify-sklearn \
  --i-classifier gg2_202409.qza \
  --i-reads ASVrep-seqs.qza \
  --o-classification Taxonomy_gg2.qza

# B. classify-consensus-blast (SILVA 138.99)
qiime feature-classifier classify-consensus-blast \
  --i-query ASVrep-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --i-reference-reads silva-138-99-seqs.qza \
  --o-classification Taxonomy_silva.qza

qiime taxa barplot \
  --i-table ASVtable_rarefied.qza \
  --i-taxonomy Taxonomy_gg2.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization Taxa_Plot.qzv
```

**결과**: Phylum~Species level 상대풍부도 barplot 확인 (Greengenes2 vs SILVA 비교)

![Taxa barplot](figures/taxa_barplot.png)

---

## 4. Phylogenetic Tree

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ASVrep-seqs.qza \
  --o-alignment align-rep-seqs.qza \
  --o-masked-alignment masked-align-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```
MAFFT(multiple alignment) → FastTree(ML tree 추정). Rooted tree는 이후 Faith's PD, UniFrac 계산에 사용.

---

## 5. Alpha Diversity (샘플 내 다양성)

```bash
qiime diversity alpha --i-table ASVtable_rarefied.qza \
  --p-metric shannon --o-alpha-diversity adiv_shannon.qza

qiime diversity alpha-phylogenetic \
  --i-table ASVtable_rarefied.qza \
  --i-phylogeny rooted-tree.qza \
  --p-metric faith_pd --o-alpha-diversity adiv_faith_pd.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity adiv_shannon.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization sig_shannon.qzv

qiime diversity alpha-rarefaction \
  --i-table ASVtable_rarefied.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 25347 \
  --m-metadata-file metadata.tsv \
  --o-visualization rarefaction.qzv
```

**사용 지표**: observed_features, chao1, shannon, simpson, simpson_e, faith_pd, singles

**결과 요약** (실제 수치로 채우기):

| Group | ASVs 평균 | Shannon 평균 | Faith_PD 평균 |
|---|---|---|---|
| Control | | | |
| Diabetes | | | |

Kruskal-Wallis (Shannon 기준): H = __, p-value = __

![Alpha diversity boxplot](figures/alpha_diversity_boxplot.png)
![Rarefaction curve](figures/rarefaction_curve.png)

---

## 6. Beta Diversity 🚧 (내일 진행 예정)

```bash
qiime diversity-lib weighted-unifrac \
  --i-table ASVtable_rarefied.qza \
  --i-phylogeny rooted-tree.qza \
  --o-distance-matrix wUniFrac-dm.qza

qiime diversity pcoa \
  --i-distance-matrix wUniFrac-dm.qza \
  --o-pcoa wUniFrac-pcoa.qza

qiime emperor plot \
  --i-pcoa wUniFrac-pcoa.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization wUniFrac-pcoa-plot.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix wUniFrac-dm.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Group \
  --p-method anosim \
  --o-visualization wUniFrac-anosim.qzv
```

- [ ] weighted / unweighted UniFrac, Bray-Curtis distance 비교
- [ ] PCoA plot으로 Control vs Diabetes 분리 여부 확인
- [ ] ANOSIM R / p-value 기록

---

## 7. ANCOM-BC (그룹 간 차등 풍부도) 🚧 예정

```bash
qiime taxa collapse --i-table ASVtable_rarefied.qza \
  --i-taxonomy Taxonomy_gg2.qza --p-level 6 \
  --o-collapsed-table ASVtable-l6.qza

qiime composition ancombc \
  --i-table ASVtable-l6.qza \
  --m-metadata-file metadata.tsv \
  --p-formula Group \
  --p-reference-levels Group::Control \
  --o-differentials l6-ancombc.qza
```

---

## 8. Shotgun Metagenome (MAGs) 🚧 추후 진행
Raw → QC/host DNA 제거 → Assembly(contigs) → Binning(MAGs) → Dereplication → Taxonomic/Functional profiling

---

## 참고
- QIIME2 docs: https://docs.qiime2.org/
- QIIME2 view (qza/qzv 확인): https://view.qiime2.org/
