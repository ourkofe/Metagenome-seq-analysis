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

## 8. Shotgun Metagenome (MAGs)
## 8. Shotgun Metagenome

Shotgun metagenome 분석은 QIIME2 기반 **MOSHPIT**(MOdular SHotgun metagenome Pipelines with Integrated provenance Tracking) 플러그인으로 진행함.  
파이프라인 흐름: `전처리 → Assembly → Binning → Dereplication → Taxonomic/Functional Annotation`

- 샘플: `sample1`, `sample2` (실습 제공 mock community paired-end FASTQ)
- 실행 환경: `quay.io/qiime2/moshpit:2026.4` (docker)

---

### 8-1. 실습 데이터 및 환경 준비

```bash
mkdir META
cd META/
wget https://data.macrogen.com/~meta/BioEdu/202607/META.zip
unzip META.zip
```

데이터 위치: `<WORKDIR>/META/Shotgun/Raw/`

```
Sample1_00_L001_R1_001.fastq.gz
Sample1_00_L001_R2_001.fastq.gz
Sample2_00_L001_R1_001.fastq.gz
Sample2_00_L001_R2_001.fastq.gz
```

MOSHPIT docker 이미지 실행:

```bash
docker run -v <WORKDIR>/META/Shotgun:/data -v <WORKDIR>/META_Ref:/ref -it quay.io/qiime2/moshpit:2026.4
```

버전 확인:

```bash
qiime info
```

```
QIIME 2 release: 2026.4
QIIME 2 version: 2026.4.0
q2cli version: 2026.4.0
```

---

### 8-2. 전처리 (Preprocessing) - Importing data

```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Raw --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path Raw.qza
```

- Semantic Type: `SampleData` (샘플 단위 분류) + `PairedEndSequences` (양방향 시퀀싱) + `WithQuality` (quality 포함)
- Input format: `CasavaOneEightSingleLanePerSampleDirFmt` (Casava 1.8, single lane, demultiplexing 완료된 per-sample 파일)

결과 확인:

```bash
qiime tools peek Raw.qza
```

```
Type:        SampleData[PairedEndSequencesWithQuality]
Data format: SingleLanePerSamplePairedEndFastqDirFmt
```

---

### 8-3. Assembly (MEGAHIT) 및 평가 (QUAST)

**Assembly**

```bash
mkdir Assembly
mosh assembly assemble-megahit --i-reads Raw.qza --p-presets meta-sensitive --p-num-cpu-threads 4 --p-min-contig-len 500 --o-contigs Assembly/contigs.qza --verbose
```

- `--p-presets meta-sensitive`: k-mer 간격을 촘촘히 늘려가며 low-abundance 종까지 최대한 contig화 (연산량↑, 정밀도↑)

**Assembly 품질 평가 (QUAST)**

```bash
mosh assembly evaluate-quast --i-contigs Assembly/contigs.qza --i-references /ref/reference-genomes.qza --p-threads 4 --p-min-contig 500 --o-reference-genomes Assembly/ref-genomes.qza --o-results-table Assembly/quast-results.qza --o-visualization Assembly/contigs-qc-quast.qzv --verbose
```

**결과 (Statistics without reference)**

| 지표 | sample1_contigs | sample2_contigs |
|---|---|---|
| # contigs | 9,997 | 10,200 |
| Largest contig | 10,884 | 77,570 |
| Total length | 10,050,842 | 20,130,423 |
| N50 | 1,040 | 2,941 |
| N90 | 570 | 770 |
| L50 | 2,753 | 1,363 |
| L90 | 8,115 | 6,932 |

- **N50**: contig를 긴 순서로 정렬 시, 누적 길이 50%를 넘기는 지점의 contig 길이 (클수록 좋음)
- **L50**: N50 조건을 만족하는 데 필요한 최소 contig 개수 (작을수록 좋음)

시각화: `figures/contigs-qc-quast.png` (QUAST 리포트 캡처)
샘플별 균주 genome fraction(%) 예시:

| Organism | sample1_contigs | sample2_contigs |
|---|---|---|
| B_subtilis | 87.724 | 99.598 |
| E_coli_K-12 | 19.273 | 14.058 |
| E_coli_O157_H7 | 56.959 | 93.927 |
| L_monocytogenes_EGD-e | 40.961 | 92.979 |
| M_tuberculosis_H37Rv | 19.631 | 76.603 |
| P_aeruginosa_PAO1 | 0.737 | 8.938 |
| S_aureus_NCTC8325 | 0.925 | 8.932 |
| S_enterica_LT2 | 24.230 | 80.680 |

---

### 8-4. Binning (MetaBAT2) 및 평가 (BUSCO)

**Contigs Indexing**

```bash
mkdir Binning
mosh assembly index-contigs --i-contigs Assembly/contigs.qza --p-threads 4 --p-seed 100 --o-index Binning/contigs-index.qza --verbose
```

**Mapping (Raw → Contigs)**

```bash
mosh assembly map-reads --i-index Binning/contigs-index.qza --i-reads Raw.qza --p-threads 4 --p-seed 100 --o-alignment-maps Binning/Raw-to-contigs-aln.qza --verbose
```

**Binning**

```bash
mosh annotate bin-contigs-metabat --i-contigs Assembly/contigs.qza --i-alignment-maps Binning/Raw-to-contigs-aln.qza --p-num-threads 4 --p-seed 100 --o-mags Binning/mags.qza --o-contig-map Binning/contig-map.qza --o-unbinned-contigs Binning/unbinned-contigs.qza --verbose
```

**Binning 품질 평가 (BUSCO)**

```bash
mosh annotate evaluate-busco --i-mags Binning/mags.qza --i-db /ref/busco-db-bacteria.qza --i-unbinned-contigs Binning/unbinned-contigs.qza --p-lineage-dataset bacteria_odb12 --p-cpu 2 --o-results Binning/busco-results.qza --o-visualization Binning/mags.qzv --verbose
```

**결과 (BUSCO summary statistics, Total MAGs = 7)**

| 지표 | Min | Median | Mean | Max |
|---|---|---|---|---|
| % single | 0 | 25.9 | 32.5 | 94 |
| % duplicated | 0 | 0 | 0.9 | 2.6 |
| % fragmented | 0 | 5.2 | 6.5 | 13.8 |
| % missing | 0.9 | 64.7 | 60.1 | 100 |
| % complete | 0 | 26.7 | 33.4 | 96.6 |
| % completeness | 0 | 35.3 | 39.9 | 99.1 |
| % contamination | 0 | 1.4 | 1.8 | 5 |
| % unbinned contigs | 70.9 | 70.9 | 73.5 | 88.9 |

- Completeness = (100 - Missing) / Total genes
- Contamination = 100 × Duplicated / Complete genes

시각화: `figures/mags-busco-qzv.png`

---

### 8-5. Dereplication

**Low quality MAG 필터링**

```bash
mkdir dRep
mosh annotate filter-mags --i-mags Binning/mags.qza --m-metadata-file Binning/busco-results.qza --p-where 'completeness>50 AND contamination<10' --p-on 'mag' --o-filtered-mags dRep/mags-filtered.qza --verbose
```

**MinHash signature 계산 (sourmash)**

```bash
qiime sourmash compute --i-sequence-file dRep/mags-filtered.qza --p-ksizes 105 --p-scaled 100 --o-min-hash-signature dRep/min-hash.qza --verbose
```

**유사도 비교**

```bash
qiime sourmash compare --i-min-hash-signature dRep/min-hash.qza --p-ksize 105 --o-compare-output dRep/min-hash-compare.qza --verbose
```

**중복 제거 (대표 MAG 선정)**

```bash
qiime annotate dereplicate-mags --i-mags dRep/mags-filtered.qza --i-distance-matrix dRep/min-hash-compare.qza --m-metadata-file Binning/busco-results.qza --p-metadata-column completeness --p-threshold 0.9 --o-dereplicated-mags dRep/mags-derep.qza --o-table dRep/table.qza --verbose
```

- Completeness 기준 대표 MAG 선정, 거리(유사도) threshold 0.9 이상인 MAG끼리 하나로 묶어 중복 제거
- Filtering 조건: `completeness > 50 AND contamination < 10` 통과한 MAG만 dereplication 대상에 포함

---

### 8-6. Taxonomic Classification (Kraken2)

```bash
mkdir Taxonomy
mosh annotate classify-kraken2 --i-seqs dRep/mags-derep.qza --i-db /ref/kraken2-db.qza --p-threads 4 --p-memory-mapping --o-reports Taxonomy/kraken2-reports-mags.qza --o-outputs Taxonomy/kraken2-hits-mags.qza --verbose
```

**MAG 단위 데이터 변환**

```bash
mosh annotate kraken2-to-mag-features --i-reports Taxonomy/kraken2-reports-mags.qza --i-outputs Taxonomy/kraken2-hits-mags.qza --o-taxonomy Taxonomy/mags-taxonomy.qza --verbose
```

결과 예시 (`mags-taxonomy.qza`):

| Feature ID | Taxon |
|---|---|
| 3b937762-... | d__Bacteria;k__Bacteria |
| 7e3f38da-... | d__Bacteria;...;o__Enterobacterales;f__Enterobacteriaceae |
| ebc4e5d1-... | d__Bacteria;...;f__Bacillaceae;g__Bacillus |

---

### 8-7. Abundance Estimation (TPM)

**MAG indexing → mapping → 길이 계산 → abundance 계산**

```bash
mkdir Abundance
mosh assembly index-derep-mags --i-mags dRep/mags-derep.qza --p-threads 4 --p-seed 100 --o-index Abundance/mags-derep-index.qza --verbose

mosh assembly map-reads --i-index Abundance/mags-derep-index.qza --i-reads Raw.qza --p-threads 4 --p-seed 100 --o-alignment-maps Abundance/Raw-to-mags-aln.qza --verbose

mosh annotate get-feature-lengths --i-features dRep/mags-derep.qza --o-lengths Abundance/mags-derep-lengths.qza --verbose

mosh annotate estimate-abundance --i-alignment-maps Abundance/Raw-to-mags-aln.qza --i-feature-lengths Abundance/mags-derep-lengths.qza --p-metric tpm --p-min-mapq 42 --p-threads 4 --o-abundances Abundance/mags-abundances.qza
```

- `--p-metric tpm`: RPKM과 달리 정규화 후 샘플 간 총합이 동일(=100만)하여 샘플 간 직접 비교 가능
- `--p-min-mapq 42`: 완벽 일치 매핑만 사용 (균주 수준 정밀 정량)

**Taxa barplot 시각화**

```bash
mosh taxa barplot --i-table Abundance/mags-abundances.qza --i-taxonomy Taxonomy/mags-taxonomy.qza --o-visualization Abundance/mags-taxa-barplot.qzv --verbose
```

Level 7 (genus) 기준 검출된 주요 속(genus):

- *Bacillus*
- Enterobacteriaceae (미분류)
- *Salmonella*
- *Listeria*
- *Mycobacterium*
- *Staphylococcus*
- *Pseudomonas*

전체 4개 샘플에서 genus-level 상대 비율 구성이 유사하게 나타남 (`figures/mags-taxa-barplot.png` 참고).

---

### 8-8. Functional Annotation (eggNOG / DIAMOND)

**Ortholog 검색 (DIAMOND)**

```bash
mkdir eggNOG
mosh annotate search-orthologs-diamond --i-seqs dRep/mags-derep.qza --i-db /ref/diamond_db.qza --p-num-cpus 4 --p-db-in-memory --o-eggnog-hits eggNOG/eggNOG_hits --o-table eggNOG/eggNOG_ft --o-loci eggNOG/eggNOG_loci --verbose
```

**기능 정보 매핑 (COG / KEGG / CAZy 등)**

```bash
mosh annotate map-eggnog --i-eggnog-hits eggNOG/eggNOG_hits.qza --i-db /ref/eggnog_db.qza --p-num-cpus 4 --p-db-in-memory --o-ortholog-annotations eggNOG/eggNOG_annotations --verbose
```

**특정 주석 추출 (KEGG pathway 빈도)**

```bash
mosh annotate extract-annotations --i-ortholog-annotations eggNOG/eggNOG_annotations.qza --p-annotation kegg_pathway --p-max-evalue 0.001 --o-annotation-frequency eggNOG/kegg_pathway_frequency --verbose
```

**MAG 풍부도 × 기능 주석 결합**

```bash
mosh annotate multiply-tables --i-table1 Abundance/mags-abundances.qza --i-table2 eggNOG/kegg_pathway_frequency.qza --o-result-table eggNOG/eggNOG_kegg_pathway_abundance --verbose
```

MAG별 KEGG pathway 빈도 예시:

| KEGG pathway | MAG1 | MAG2 | MAG3 |
|---|---|---|---|
| map01100 (Metabolic pathways) | 308 | 790 | 656 |
| map01110 | 144 | 333 | 305 |
| map01120 | 101 | 318 | 206 |
| map02010 (ABC transporters) | 79 | 197 | 153 |

→ 샘플별 총 pathway 존재량으로 환산 (`eggNOG_kegg_pathway_abundance`)까지 산출 완료.

---

### 8-9. 전체 파이프라인 요약

| 단계 | 사용 도구 | 주요 output |
|---|---|---|
| 전처리 | `qiime tools import` | Raw.qza |
| Assembly | MEGAHIT (`mosh assembly assemble-megahit`) | contigs.qza |
| Assembly 평가 | QUAST | N50, L50, genome fraction |
| Binning | MetaBAT2 (`mosh annotate bin-contigs-metabat`) | mags.qza |
| Binning 평가 | BUSCO | completeness / contamination |
| Dereplication | sourmash (MinHash) | mags-derep.qza |
| Taxonomic classification | Kraken2 | mags-taxonomy.qza |
| Abundance | TPM 기반 정량 | mags-abundances.qza |
| Functional annotation | eggNOG / DIAMOND | KEGG/COG/CAZy pathway 빈도 |

최종 dereplicated MAG 7개 기준, genus-level taxonomy 및 KEGG pathway 기능 주석까지 완료.

---

## 참고
- QIIME2 docs: https://docs.qiime2.org/
- QIIME2 view (qza/qzv 확인): https://view.qiime2.org/
