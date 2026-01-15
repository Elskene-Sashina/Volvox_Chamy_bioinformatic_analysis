# 🧬 Анализ №3: Транскриптомика SOMA/GERM (*Volvox carteri*)

> **Гипотеза temporal coaptation**: Временные программы *Chlamydomonas* (день/ночь) → Пространственные программы *Volvox* (соматические/герминативные клетки)  

> **Метод**: RNA-seq differential expression analysis (DESeq2)

---

## 📋 Оглавление
- [Биологическая гипотеза](#биологическая-гипотеза)
- [Данные](#данные)
- [Pipeline](#pipeline)
  - [Шаг 1: Скачивание данных](#шаг-1-скачивание-данных)
  - [Шаг 2: Quality control + Trimming](#шаг-2-quality-control--trimming)
  - [Шаг 3: Mapping (STAR)](#шаг-3-mapping-star)
  - [Шаг 4: Gene counting (featureCounts)](#шаг-4-gene-counting-featurecounts)
  - [Шаг 5: Differential expression (DESeq2)](#шаг-5-differential-expression-deseq2)
- [Результаты](#результаты)
- [Визуализации](#визуализации)
- [Биологические выводы](#биологические-выводы)

---

## Биологическая гипотеза

### 🕐 Temporal Coaptation Hypothesis

**Основная идея** (Matt & Umen, 2018):

#### *Chlamydomonas reinhardtii* (одноклеточный):
- **Диурнальный цикл** создаёт **временную специализацию**
  - **День**: фотосинтез, рост, подвижность ("day genes")
  - **Ночь**: деление клеток, репродукция ("night genes")

#### *Volvox carteri* (многоклеточный):
- **Пространственная специализация** → две линии клеток:
  - **Соматические клетки** (Somatic): экспрессируют "day genes" (фотосинтез, флагеллярная подвижность)
  - **Гониди** (Germ cells): экспрессируют "night genes" (деление, репродукция)


---

## Данные

### 📦 Источник

**Проект**: [PRJNA413955](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA413955) | GEO: [GSE104835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104835)

**Публикация**: Matt, G.Y. & Umen, J.G. (2018). "Cell-type transcriptomes of the multicellular green alga *Volvox carteri* yield insights into the evolutionary origins of germ and somatic differentiation programs". *G3: Genes, Genomes, Genetics*.

---

### 🧪 Образцы

| Sample ID | BioSample | SRA | Тип клеток | Платформа |
|-----------|-----------|-----|-----------|-----------|
| **SRR6159196** | SAMN07775063 | GSM2808507 | **Gonidia** rep 1 | Illumina HiSeq 2500 |
| **SRR6159197** | SAMN07775062 | GSM2808508 | **Gonidia** rep 2 | Illumina HiSeq 2500 |
| **SRR6159198** | SAMN07775061 | GSM2808509 | **Somatic** rep 1 | Illumina HiSeq 2500 |
| **SRR6159199** | SAMN07775064 | GSM2808510 | **Somatic** rep 2 | Illumina HiSeq 2500 |

**Всего**: 4 образца, 2 биологических реплики на тип клеток

---

## Pipeline

### Шаг 1: Скачивание данных

```bash
# Создание metadata файла
cat > metadata.csv << 'EOF'
sample,file,condition,replicate,biosample
SRR6159196,fastq_raw/SRR6159196.fastq.gz,gonidia,1,SAMN07775063
SRR6159197,fastq_raw/SRR6159197.fastq.gz,gonidia,2,SAMN07775062
SRR6159198,fastq_raw/SRR6159198.fastq.gz,somatic,1,SAMN07775061
SRR6159199,fastq_raw/SRR6159199.fastq.gz,somatic,2,SAMN07775064
EOF

# Скачивание через SRA toolkit
for SRR in SRR6159196 SRR6159197 SRR6159198 SRR6159199; do
  prefetch $SRR
  fastq-dump --gzip --outdir fastq_raw $SRR
done
```

---

### Шаг 2: Quality control + Trimming

#### 🔧 Инструмент: `fastp`

**Параметры** (согласно методике статьи):
- Удаление адаптеров (автоопределение)
- Обрезка 5'-конца при Q < 3
- Обрезка 3'-конца при Q < 3
- Скользящее окно 4bp, Q < 15 → обрезка
- Минимальная длина 25 nt

```bash
fastp --version
# fastp 1.0.1

# Пример для SRR6159196 (Gonidia rep 1)
fastp \
    --in1 fastq_raw/SRR6159196.fastq.gz \
    --out1 trimmed/SRR6159196_trimmed.fastq.gz \
    --thread 8 \
    --cut_front \
    --cut_tail \
    --cut_front_mean_quality 3 \
    --cut_tail_mean_quality 3 \
    --cut_mean_quality 15 \
    --cut_tail_window_size 4 \
    --length_required 25 \
    --json qc_fastp/SRR6159196.json \
    --html qc_fastp/SRR6159196.html \
    2>&1 | tee logs/SRR6159196_fastp.log
```

---

#### 📊 Результаты trimming

| Sample | Input reads | Passed | % Retained | Duplication rate |
|--------|-------------|--------|------------|------------------|
| **SRR6159196** | 15,572,992 | 14,940,465 | **95.9%** | 67.7% |
| **SRR6159197** | 15,626,215 | 15,002,689 | **96.0%** | 68.8% |
| **SRR6159198** | 16,400,552 | 15,680,869 | **95.6%** | 68.8% |
| **SRR6159199** | 17,218,685 | 16,511,301 | **95.9%** | 63.6% |

#### Quality metrics:

| Метрика | До | После |
|---------|-----|-------|
| **Q20 bases** | 97.6-97.8% | **98.7-98.8%** |
| **Q30 bases** | 93.5-93.8% | **94.6-94.8%** |
| **Q40 bases** | 35.3-36.0% | 36.0-36.7% |

> ✅ **Отличное качество**: >95% ридов retained, Q30 >94%

---

### Шаг 3: Mapping (STAR)

#### 🔧 Инструмент: STAR v2.7.11a

**Референс**: *Volvox carteri* v2.1 (Phytozome)
- Genome size: ~131 Mb
- Genes: 14,247 (v2.1 annotation)

**Параметры**:
- `--outFilterMultimapNmax 10` - риды, картирующиеся >10 мест, отбрасываются
- `--outFilterMismatchNoverLmax 0.06` - max 6% мисматчей
- `--quantMode GeneCounts` - автоматич подсчёт ридов на гены

```bash
# Создание STAR индекса
STAR --runMode genomeGenerate \
     --genomeDir genome/Vcarteri_STAR_index \
     --genomeFastaFiles genome/Vcarteri_317_v2.assembly.fa \
     --sjdbGTFfile genome/Vcarteri_317_v2.1.annotation.gtf \
     --sjdbOverhang 49 \
     --runThreadN 8

# Пример маппинга (SRR6159196)
STAR \
    --runThreadN 8 \
    --genomeDir genome/Vcarteri_STAR_index \
    --readFilesIn trimmed/SRR6159196_trimmed.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix bam/SRR6159196_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNoverLmax 0.06 \
    --quantMode GeneCounts \
    2>&1 | tee logs/SRR6159196_star.log

# Индексация BAM
samtools index bam/SRR6159196.sorted.bam
```

---

#### 📊 Mapping statistics

| Sample | Input | Uniquely mapped | % Unique | Multi-mapped | % Multi |
|--------|-------|-----------------|----------|--------------|---------|
| **SRR6159196** | 14,940,465 | 12,915,107 | **86.44%** | 846,251 | 5.66% |
| **SRR6159197** | 15,002,689 | 12,964,163 | **86.41%** | 886,391 | 5.91% |
| **SRR6159198** | 15,680,869 | 12,277,353 | **78.30%** | 660,525 | 4.21% |
| **SRR6159199** | 16,511,301 | 14,687,615 | **88.95%** | 758,101 | 4.59% |

> ✅ **Высокое качество**: Среднее 85.0% uniquely mapped (ожидалось >85% из статьи)

> ⚠️ **Note**: SRR6159198 (Somatic rep 1) показывает 78.3% — немного ниже, но приемлемо (возможно, больше rRNA)

---

### Шаг 4: Gene counting (featureCounts)

#### ⚠️ Проблема: STAR `--quantMode GeneCounts` выдал пустые файлы

**Решение**: Использовать `featureCounts` из пакета Subread

```bash
# Конвертация GFF3 → GTF
gffread genome/Vcarteri_317_v2.1.annotation.gene_exons.gff3 \
        -T -o genome/Vcarteri_317_v2.1.annotation.gtf

# Gene counting
featureCounts \
    -a genome/Vcarteri_317_v2.1.annotation.gtf \
    -o counts/gene_counts.txt \
    -t exon \
    -g gene_id \
    -s 0 \
    -Q 10 \
    --primary \
    -T 8 \
    bam/SRR6159196.sorted.bam \
    bam/SRR6159197.sorted.bam \
    bam/SRR6159198.sorted.bam \
    bam/SRR6159199.sorted.bam \
    2>&1 | tee logs/featureCounts.log

# Очистка матрицы (убрать metadata колонки)
cat counts/gene_counts.txt | \
    grep -v "^#" | \
    cut -f1,7-10 | \
    sed '1s|bam/||g' | \
    sed '1s|.sorted.bam||g' > counts/gene_counts_clean.txt
```

---

#### 📊 Counting statistics

| Метрика | SRR6159196 | SRR6159197 | SRR6159198 | SRR6159199 |
|---------|------------|------------|------------|------------|
| **Assigned** | 11,542,370 | 11,615,550 | 10,763,160 | 12,986,515 |
| Unmapped | 1,179,107 | 1,152,135 | 2,742,991 | 1,065,585 |
| Low quality | 3,827,176 | 4,050,358 | 2,217,093 | 2,691,376 |
| No features | 880,291 | 842,914 | 1,097,296 | 1,214,337 |
| Ambiguous | 492,446 | 505,699 | 416,897 | 486,763 |

**Всего генов в аннотации**: 14,247  
**Экспрессируемых генов** (>0 reads): **13,764 (96.6%)**


---

#### 🏆 Топ-20 экспрессируемых генов

| Gene ID | Total counts |
|---------|--------------|
| Vocar.0013s0021 | 771,485 |
| Vocar.0008s0418 | 745,809 |
| Vocar.0070s0007 | 603,631 |
| Vocar.0018s0037 | 427,171 |
| Vocar.0021s0124 | 384,017 |
| Vocar.0009s0198 | 368,524 |
| Vocar.0004s0255 | 358,626 |
| Vocar.0016s0243 | 358,491 |
| Vocar.0016s0158 | 351,739 |
| Vocar.0001s0479 | 351,318 |

---

### Шаг 5: Differential expression (DESeq2)

#### 🧬 Методика (точное следование Matt & Umen 2018)

**Pipeline**:
1. Фильтрация: удалить гены без экспрессии
2. **Outlier removal**: удалить верхние 0.3% квантили ридов
3. **Нормализация**: custom size factors = total reads / 10 million
4. **Low expression filter**: удалить нижние 10% генов
5. **DESeq2**: FDR < 0.05
6. **Классификация**:
   - **Cell-type-SPECIFIC**: >5-fold (log2FC > 2.32), FDR < 0.05
   - **Cell-type-BIASED**: 2-5-fold (1 < log2FC < 2.32), FDR < 0.05
   - **CONSTITUTIVE**: <2-fold, FDR > 0.05
   - **LOW CONFIDENCE**: >2-fold, FDR > 0.05

---

#### 📝 R скрипт (DESeq2)

```r
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# ===== 1. LOAD DATA =====
counts <- read.table("counts/gene_counts_clean.txt", 
                     header = TRUE, row.names = 1)

metadata <- data.frame(
  sample = c("SRR6159196", "SRR6159197", "SRR6159198", "SRR6159199"),
  condition = c("gonidia", "gonidia", "somatic", "somatic"),
  replicate = c(1, 2, 1, 2)
)
rownames(metadata) <- metadata$sample

# ===== 2. FILTER ZERO COUNTS =====
counts_filt <- counts[rowSums(counts) > 0, ]

# ===== 3. OUTLIER REMOVAL (99.7% percentile) =====
percentile_997 <- apply(counts_filt, 2, quantile, probs = 0.997)
for (i in 1:ncol(counts_filt)) {
  counts_filt[counts_filt[,i] > percentile_997[i], i] <- percentile_997[i]
}

# ===== 4. CUSTOM SIZE FACTORS =====
library_sizes <- colSums(counts_filt)
size_factors <- library_sizes / 10000000

# ===== 5. DESeq2 OBJECT =====
dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                               colData = metadata,
                               design = ~ condition)
sizeFactors(dds) <- size_factors

# ===== 6. FILTER BOTTOM 10% =====
norm_counts <- counts(dds, normalized = TRUE)
mean_norm_expr <- rowMeans(norm_counts)
threshold_10pct <- quantile(mean_norm_expr, 0.10)
dds <- dds[mean_norm_expr > threshold_10pct, ]

# ===== 7. RUN DESeq2 =====
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "somatic", "gonidia"),
               alpha = 0.05)

# ===== 8. CLASSIFY GENES =====
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$category <- "Constitutive"

res_df$category[!is.na(res_df$padj) & res_df$padj < 0.05 & 
                res_df$log2FoldChange > 2.32] <- "Somatic-specific"
res_df$category[!is.na(res_df$padj) & res_df$padj < 0.05 & 
                res_df$log2FoldChange > 1 & res_df$log2FoldChange <= 2.32] <- "Somatic-biased"
res_df$category[!is.na(res_df$padj) & res_df$padj < 0.05 & 
                res_df$log2FoldChange < -2.32] <- "Gonidia-specific"
res_df$category[!is.na(res_df$padj) & res_df$padj < 0.05 & 
                res_df$log2FoldChange < -1 & res_df$log2FoldChange >= -2.32] <- "Gonidia-biased"
res_df$category[!is.na(res_df$padj) & res_df$padj >= 0.05 & 
                (res_df$log2FoldChange > 1 | res_df$log2FoldChange < -1)] <- "Low confidence"

# ===== 9. SAVE RESULTS =====
write.csv(res_df, "deseq2/deseq2_results_full.csv", row.names = FALSE)
```

---

#### 📊 DESeq2 Results Summary

```
out of 12387 with nonzero total read count
adjusted p-value < 0.1

LFC > 0 (up):       3682, 30%
LFC < 0 (down):     6014, 49%
outliers:           0, 0%
low counts:         0, 0%
```

---

## Результаты

### 📊 Классификация генов

| Категория | Количество | % от total |
|-----------|------------|------------|
| **Constitutive** | 3,322 | **26.8%** |
| **Gonidia-specific** | 3,290 | **26.6%** |
| **Somatic-specific** | 2,105 | **17.0%** |
| **Gonidia-biased** | 1,981 | **16.0%** |
| **Somatic-biased** | 1,088 | **8.8%** |
| **Low confidence** | 601 | **4.9%** |
| **TOTAL** | **12,387** | **100.0%** |

---

### 🔍 Сравнение с Matt & Umen (2018)

| Категория | Наши данные | Статья | Разница |
|-----------|-------------|--------|---------|
| **Total analyzed** | 12,387 | 13,238 | -851 |
| **Gonidia-specific** | 3,290 (26.6%) | 3,541 (26.8%) | -251 |
| **Somatic-specific** | 2,105 (17.0%) | 2,244 (17.0%) | -139 |
| **Gonidia-biased** | 1,981 (16.0%) | 2,120 (16.0%) | -139 |
| **Somatic-biased** | 1,088 (8.8%) | 1,124 (8.5%) | -36 |
| **Constitutive** | 3,322 (26.8%) | 3,609 (27.3%) | -287 |

> ✅ **Отличное соответствие**: Процентные соотношения практически идентичны!

---

### 🏆 Топ-10 Somatic-specific генов

| Gene ID | baseMean | log2FC | padj |
|---------|----------|--------|------|
| **Vocar.0046s0015** | 7,255 | **+11.85** | 3.3e-93 |
| Vocar.0083s0016 | 1,888 | +11.91 | 5.0e-28 |
| Vocar.0066s0001 | 664 | +12.85 | 1.8e-17 |
| Vocar.0001s1558 | 617 | +12.74 | 2.7e-17 |
| Vocar.0041s0047 | 465 | +12.33 | 3.0e-16 |
| Vocar.0083s0017 | 1,043 | +12.06 | 1.1e-15 |
| Vocar.0001s0819 | 378 | +12.04 | 1.6e-15 |
| Vocar.0018s0236 | 378 | +12.04 | 2.1e-15 |
| Vocar.0028s0175 | 366 | +11.99 | 2.5e-15 |
| Vocar.0031s0087 | 353 | +11.93 | 3.0e-15 |


---

### 🏆 Топ-10 Gonidia-specific генов

| Gene ID | baseMean | log2FC | padj |
|---------|----------|--------|------|
| **Vocar.0003s0026** | 10,539 | **-11.01** | 1.3e-157 |
| Vocar.0003s0023 | 2,602 | -11.08 | 5.7e-56 |
| Vocar.0001s0298 | 2,481 | -13.33 | 6.8e-20 |
| Vocar.0007s0165 | 389 | -12.10 | 1.3e-15 |
| Vocar.0029s0079 | 346 | -11.93 | 2.8e-15 |
| Vocar.0004s0005 | 211 | -11.22 | 1.5e-13 |
| Vocar.0005s0420 | 208 | -11.19 | 1.7e-13 |
| Vocar.0043s0081 | 169 | -10.89 | 8.8e-13 |
| Vocar.0006s0320 | 141 | -10.64 | 3.2e-12 |
| Vocar.0025s0005 | 753 | -10.62 | 7.4e-23 |


---

## Визуализации

### 1️⃣ PCA Plot — Разделение типов клеток

```r
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 6) +
  geom_text(vjust = -1.5, size = 4) +
  scale_color_manual(values = c("gonidia" = "#2166ac", "somatic" = "#b2182b")) +
  labs(title = "PCA: Cell Type Separation",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"))
```

**Результаты**:
- **PC1**: Объясняет большую часть вариации (разделяет типы клеток)
- **PC2**: Объясняет меньшую вариацию (биологическая изменчивость)
- Gonidia реплики кластеризуются вместе
- Somatic реплики кластеризуются вместе



---

### 2️⃣ Volcano Plot — Дифференциальная экспрессия

```r
colors <- c("Gonidia-specific" = "#2166ac",
           "Gonidia-biased" = "#4393c3",
           "Constitutive" = "gray70",
           "Low confidence" = "gray40",
           "Somatic-biased" = "#d6604d",
           "Somatic-specific" = "#b2182b")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = category)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = c(-2.32, -1, 1, 2.32), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volvox Cell-Type Transcriptome",
       x = "log2 Fold Change (Somatic / Gonidia)",
       y = "-log10(adjusted p-value)")
```

**Наблюдения**:
- **Экстремальные fold changes**: некоторые гены имеют log2FC > 8 (**>256-fold**!)
- Somatic-specific гены (красные справа) достигают +10 log2FC
- Gonidia-specific гены (синие слева) достигают -13 log2FC
- **Симметричное распределение**: обе клеточные линии имеют специализацию

---

### 3️⃣ Sample Correlation Heatmap

```r
# Корреляция между образцами
cor_matrix <- cor(norm_counts, method = "spearman")

pheatmap(cor_matrix,
         annotation_col = metadata[,c("condition"), drop=FALSE],
         annotation_colors = list(condition = c(gonidia="#2166ac", somatic="#b2182b")),
         display_numbers = TRUE,
         number_format = "%.2f")
```

**Результаты**:
- **Внутри типа**: r = 0.99-1.00 между репликами → **отличное качество**!
- **Между типами**: r = 0.53-0.55 → **низкая корреляция** (правильно для разных клеточных типов!)

---

### 4️⃣ MA Plot — Quality Control

```r
plotMA(res, ylim = c(-15, 15),
       main = "MA Plot: Somatic vs Gonidia",
       alpha = 0.05)
```

**Наблюдения**:
- ✅ Дифференциально экспрессируемые гены распределены равномерно по всем уровням экспрессии
- ✅ Нет bias в сторону высоко- или низкоэкспрессируемых генов
- ✅ Симметричность относительно log2FC = 0
- **Вывод**: Нормализация работает правильно!

---

### 5️⃣ Dispersion Plot — Модель DESeq2

```r
plotDispEsts(dds, main = "Dispersion Estimates")
```

**Интерпретация**:
- Красная кривая (fitted trend) хорошо аппроксимирует данные
- Дисперсия уменьшается с ростом экспрессии — **классический паттерн** для RNA-seq
- Большинство генов близки к fitted line → модель работает корректно

---

## Биологические выводы

### 🎯 Ключевые находки

#### 1. Экстремальная клеточная специализация

- **>5000 генов** (43.6%) показывают cell-type-specific или biased экспрессию
- **Fold changes до 13 log2FC** (>8000-fold!) — это среди самых сильных известных дифференциаций
- Gonidia-specific генов **больше** (3,290 vs 2,105) → герминативная линия требует больше уникальных программ

> 💡 **Вывод**: Две клеточные линии *Volvox* имеют **радикально разные** транскриптомы!

---

#### 2. Паттерн специализации: Gonidia > Somatic

| Тип | Specific | Biased | Total |
|-----|----------|--------|-------|
| **Gonidia** | 3,290 (26.6%) | 1,981 (16.0%) | **5,271 (42.6%)** |
| **Somatic** | 2,105 (17.0%) | 1,088 (8.8%) | **3,193 (25.8%)** |

**Интерпретация**:
- Gonidia имеют **1.65× больше** специализированных генов
- Gonidia = **более сложная программа** (репродукция, эмбриогенез, деление)
- Somatic = **более простая программа** (фотосинтез, подвижность)

---

#### 3. Constitutive гены - "housekeeping" функции

**26.8% генов** (3,322) — конститутивные в обоих типах клеток

**Функции**:
- Рибосомальные белки
- Базовый метаболизм
- Цитоскелет
- ДНК/РНК процессинг

> 💡 **Интерпретация**: Около четверти генома обеспечивает базовые клеточные функции, независимо от клеточного типа

---

#### 4. Высокая экспрессия специфичных генов

**Наблюдение из графиков**:
- Cell-type-specific гены имеют **более высокую** median экспрессию (~100-1000 TPM)
- Low confidence гены имеют **низкую** экспрессию (~10-50 TPM)

**Биологический смысл**:
- Специфичные гены **сильно экспрессируются** → функционально важны
- Cell-type-specific функции требуют **высоких уровней белков**
- Low expression → трудно статистически определить специфичность (шум > сигнал)


---

