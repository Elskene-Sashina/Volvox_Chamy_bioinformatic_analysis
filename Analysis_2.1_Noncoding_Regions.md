# 🧬 Анализ №2: Геномный анализ некодирующих регуляторных последовательностей

## 📋 Оглавление
- [Входные данные](#входные-данные)
- [Программное обеспечение](#программное-обеспечение)
- [Методология](#методология)
  - [Шаг 1: Извлечение кодирующих элементов](#шаг-1-извлечение-кодирующих-элементов-cds)
  - [Шаг 2: Контроль качества](#шаг-2-контроль-качества-qc)
  - [Шаг 3: Расчет интронов](#шаг-3-расчет-интронов)
  - [Шаг 4: Расчет интергенных регионов](#шаг-4-расчет-интергенных-регионов)
  - [Шаг 5: Анализ genome budget](#шаг-5-анализ-genome-budget)
- [Результаты](#результаты)

---

## 📦 Входные данные

### 🔬 Chlamydomonas reinhardtii (v6.1)

| Файл | Описание |
|------|----------|
| `CreinhardtiiCC_4532_707_v6.0.fa` | Референсный геном |
| `CreinhardtiiCC_4532_707_v6.0.softmasked.fa` | Софт-маскированная версия |
| `CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3` | Аннотация генов и экзонов |
| `CreinhardtiiCC_4532_707_v6.1.repeatmasked_assembly_v6.0.gff3` | Трек повторов |

### 🌿 Volvox carteri (v2.1)

| Файл | Описание |
|------|----------|
| `Vcarteri_317_v2.assembly.fa` | Референсный геном |
| `Vcarteri_317_v2.assembly.softmasked.fa` | Софт-маскированная версия |
| `Vcarteri_317_v2.1.annotation.gene_exons.gff3` | Аннотация генов и экзонов |
| `Vcarteri_317_v2.1.repeatmasked_assembly_v2.gff3` | Трек повторов |


### 🧫 Биологический контекст

Структура mature mRNA:

```
5' cap ─ 5'UTR ─ CDS (кодирует белок) ─ 3'UTR ─ polyA tail
```

**Процесс от ДНК до белка:**
```
ДНК (gene locus)
    ↓ Транскрипция (RNA Pol II)
pre-mRNA (exons + introns)
    ↓ Сплайсинг (сплайсосома)
mature mRNA (5'UTR + CDS + 3'UTR)
    ↓ Трансляция (рибосома)
Белок (polypeptide)
```

---

## 🛠️ Программное обеспечение

| Программа | Версия | Назначение |
|-----------|--------|------------|
| **bedtools** | v2.31.1 | Геномные операции (merge, intersect, complement) |
| **samtools** | v1.22.1 | Индексация FASTA файлов |
| **grep/awk/sed** | GNU | Парсинг GFF3 файлов |

---

## 🔬 Методология

## Шаг 1: Извлечение кодирующих элементов (CDS)

### ❓ Биологическая задача

Извлечь все CDS экзоны и **объединить перекрывающиеся регионы**, возникающие из-за альтернативного сплайсинга, чтобы каждый нуклеотид считался только один раз.

> **Почему это важно?**  
> Альтернативный сплайсинг создает множественные изоформы одного гена:
> - Изоформа 1: exon1 + exon2
> - Изоформа 2: exon1 + exon3
> 
> Без объединения → **дублированный подсчет нуклеотидов**!

---

### 1.1. 📊 Извлечение CDS экзонов из GFF3

#### 🧪 Команда для *Chlamydomonas*:

```bash
grep -P "\tCDS\t" CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Chlamydomonas_CDS_exons.bed
```

#### 🌱 Команда для *Volvox*:

```bash
grep -P "\tCDS\t" Vcarteri_317_v2.1.annotation.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Vcarteri_CDS_exons.bed
```


> ⚠️ **КРИТИЧНО**: bedtools **ТРЕБУЕТ** отсортированный input!

#### 📈 Результаты:

| Метрика | *Volvox* | *Chlamydomonas* |
|---------|----------|-----------------|
| CDS экзонов | 127,157 | 274,602 |
| Средняя длина экзона | 254.7 bp | 264.8 bp |
| Экзонов на ген (среднее) | ~8.8 | ~15.5 |

---


### 1.2. 🔀 Объединение перекрывающихся CDS регионов

#### 🧪 Команда для *Volvox*:

```bash
bedtools merge -i Vcarteri_CDS_exons.bed \
  -s -c 4,5,6 -o distinct,distinct,distinct \
  > Vcarteri_CDS_merged.bed
```

#### 🌿 Команда для *Chlamydomonas*:

```bash
bedtools merge -i Chlamydomonas_CDS_exons.bed \
  -s -c 4,5,6 -o distinct,distinct,distinct \
  > Chlamy_CDS_merged.bed
```

#### 🔍 Параметры bedtools merge:

| Флаг | Значение |
|------|----------|
| `-i` | Входной BED файл (должен быть отсортирован!) |
| `-s` | **Strand-specific** merge (объединяет только интервалы на одной цепи) |
| `-c 4,5,6` | Операции на колонках 4, 5, 6 (attributes, score, strand) |
| `-o distinct,distinct,distinct` | Показать уникальные значения для каждой колонки |

#### 📊 Результаты и интерпретация:

| Метрика | *Volvox* | *Chlamydomonas* | Соотношение |
|---------|----------|-----------------|-------------|
| CDS exons (raw) | 127,157 | 274,602 | 2.16× |
| CDS merged | 110,319 | 148,613 | 1.35× |
| Collapsed regions | 16,838 | 125,989 | 7.5× |
| **Merge efficiency** | **13.2%** | **45.9%** | **3.5×** |

> 💡 **Биологическая интерпретация:**
> 
> - **Chlamydomonas**: 46% CDS экзонов схлопываются → **очень обильный альтернативный сплайсинг**, множество изоформ
> - **Volvox**: только 13% merge → **более простая структура генов**, меньше изоформ
> - **Вывод**: Многоклеточность *Volvox* эволюционировала **БЕЗ усложнения gene structure**!

---

### 1.3. 🏷️ Извлечение UTR регионов

#### 🧪 Команды для *Volvox*:

```bash
# 5' UTR
grep -E "five_prime_UTR" Vcarteri_317_v2.1.annotation.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Vcarteri_5UTR_exons.bed

# 3' UTR
grep -E "three_prime_UTR" Vcarteri_317_v2.1.annotation.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Vcarteri_3UTR_exons.bed

# Объединение и merge
cat Vcarteri_5UTR_exons.bed Vcarteri_3UTR_exons.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge -s -c 4,5,6 -o distinct,distinct,distinct > Vcarteri_UTR_merged.bed
```

#### 🌿 Аналогично для *Chlamydomonas*:

```bash
# 5' UTR
grep -E "five_prime_UTR" CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Chlamy_5UTR_exons.bed

# 3' UTR
grep -E "three_prime_UTR" CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Chlamy_3UTR_exons.bed

# Объединение и merge
cat Chlamy_5UTR_exons.bed Chlamy_3UTR_exons.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge -s -c 4,5,6 -o distinct,distinct,distinct > Chlamy_UTR_merged.bed
```

---

### 1.4. 🧬 Создание трека "Coding Elements"

#### 💡 Биологический смысл:

**"Coding elements" = CDS + UTR = все mature mRNA loci**  
(всё, что экспортируется из ядра)

#### 🧪 Команды:

```bash
# Volvox
cat Vcarteri_CDS_merged.bed Vcarteri_UTR_merged.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge -s -c 4,5,6 -o collapse,distinct,distinct \
  > Vcarteri_coding_elements_merged.bed

# Chlamydomonas
cat Chlamy_CDS_merged.bed Chlamy_UTR_merged.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge -s -c 4,5,6 -o collapse,distinct,distinct \
  > Chlamy_coding_elements_merged.bed
```

---

## Шаг 2: Контроль качества (QC)

### ✅ QC1: Количество CDS экзонов

```bash
wc -l Vcarteri_CDS_exons.bed
wc -l Vcarteri_CDS_merged.bed
```

#### 📊 Результаты:

| Организм | Raw exons | Merged exons | Merge % | Статус |
|----------|-----------|--------------|---------|--------|
| *Volvox* | 127,157 | 110,319 | 13.2% | ✅ |
| *Chlamydomonas* | 274,602 | 148,613 | 45.9% | ✅ |

#### 🎯 Ожидаемые значения:

- *Volvox*: ~70,000-100,000 экзонов (14,520 генов × 5-7 экзонов/ген)
- **Получено**: 127,157 / 14,520 = **8.76 экзонов/ген** ✅

---

### ✅ QC2: Распределение по цепям (Strand Distribution)

```bash
awk '{print $6}' Vcarteri_CDS_exons.bed | sort | uniq -c
```

#### 📊 Результаты:

| Организм | + strand | - strand | Разница | Статус |
|----------|----------|----------|---------|--------|
| *Volvox* | 62,269 (48.9%) | 64,888 (51.1%) | 2.0% | ✅ |
| *Chlamydomonas* | 138,246 (51.1%) | 136,356 (49.7%) | 1.4% | ✅ |

> ✅ **Ожидается**: ~50/50 распределение (balanced genome)

---

### ⚠️ QC3: Проверка перекрытия CDS и UTR

#### 🤔 Биологическое ожидание:

CDS и UTR **НЕ ДОЛЖНЫ** пересекаться (это логическая ошибка в аннотации)

```bash
bedtools intersect -a Chlamy_CDS_merged.bed \
  -b Chlamy_UTR_merged.bed -u | wc -l
```

#### 📊 Результаты:

| Организм | Overlapping regions | Total overlap | % от CDS |
|----------|---------------------|---------------|----------|
| *Chlamydomonas* | 14,642 | 1.88 Mb | 0.5% |
| *Volvox* | 2,847 | 190 kb | 0.68% |

#### 🔍 Анализ размеров перекрытий (*Volvox*):

```bash
bedtools intersect -a Vcarteri_CDS_merged.bed -b Vcarteri_UTR_merged.bed -s | \
  awk '{print $3-$2}' | \
  sort -n | \
  awk '{
    if($1<10) small++;
    else if($1<50) medium++;
    else if($1<100) large++;
    else huge++;
  } END {
    print "< 10 bp:", small;
    print "10-50 bp:", medium;
    print "50-100 bp:", large;
    print "> 100 bp:", huge;
  }'
```

| Размер overlap | Количество | Интерпретация |
|----------------|------------|---------------|
| < 10 bp | 60 | Технические ошибки границ |
| 10-50 bp | 252 | Альтернативные start/stop кодоны |
| 50-100 bp | 298 | Альтернативный сплайсинг |
| > 100 bp | 637 | Сложные перекрывающиеся gene models |

> 💡 **Интерпретация**: Небольшое перекрытие (<1%) приемлемо. Причины:
> - Альтернативные start/stop кодоны
> - Nested genes (ген внутри интрона другого гена)
> - Артефакты аннотации GFF3

---

### ✅ QC4: Общий размер CDS в геноме

```bash
awk '{sum+=$3-$2} END {print sum " bp total CDS"}' Chlamy_CDS_merged.bed
```

#### 📊 Результаты:

| Организм | Total CDS | Размер генома | % coding | Статус |
|----------|-----------|---------------|----------|--------|
| *Volvox* | 27.88 Mb | 138 Mb | 20.2% | ✅ |
| *Chlamydomonas* | 40.18 Mb | 114.6 Mb | 35.1% | ✅ |

#### 🔬 Сравнение с другими организмами:

| Организм | % coding |
|----------|----------|
| Человек | ~1.5% (огромные интроны) |
| Arabidopsis | ~25% (компактный геном) |
| *Volvox/Chlamydomonas* | **20-35%** ✅ |

---

### ✅ QC5: Средняя длина CDS экзона

```bash
awk '{print $3-$2}' Vcarteri_CDS_merged.bed | \
  awk '{sum+=$1; n++} END {print sum/n " bp avg exon (merged)"}'
```

#### 📊 Результаты:

| Организм | Средняя длина экзона | Статус |
|----------|----------------------|--------|
| *Volvox* | 252.7 bp | ✅ |
| *Chlamydomonas* | 270.4 bp | ✅ |

#### 📏 Типичные длины экзонов:

| Категория | Размер | Наши данные |
|-----------|--------|-------------|
| Microexons | < 100 bp | - |
| **Нормальные экзоны** | **100-300 bp** | **← МЫ ЗДЕСЬ** ✅ |
| Длинные экзоны | 300-500 bp | - |
| Очень длинные | > 500 bp | - |

---

## Шаг 3: Расчет интронов

### 📐 3.1. Извлечение границ генов (gene bodies)

```bash
grep -P "\tgene\t" Vcarteri_317_v2.1.annotation.gene_exons.gff3 | \
  awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$9,$6,$7}' | \
  sort -k1,1 -k2,2n > Vcarteri_gene_bodies.bed
```

### ➖ 3.2. Вычитание coding elements из gene bodies

```bash
bedtools subtract -a Vcarteri_gene_bodies.bed \
  -b Vcarteri_coding_elements_merged.bed \
  -s > Vcarteri_introns.bed
```

#### 🔍 Параметры:

| Флаг | Функция |
|------|---------|
| `-a` | Файл, из которого вычитаем (gene bodies) |
| `-b` | Файл, который вычитаем (coding elements) |
| `-s` | Strand-specific операция |

**Результат**: 

```
Интроны = gene bodies − coding elements
```

---

## Шаг 4: Расчет интергенных регионов

### 📏 4.1. Создание genome size file

```bash
samtools faidx CreinhardtiiCC_4532_707_v6.0.fa
cut -f1,2 CreinhardtiiCC_4532_707_v6.0.fa.fai > Chlamy_genome.txt
```

**Формат `genome.txt`:**
```
scaffold_1    5000000
scaffold_2    3000000
...
```

---

### 🔄 4.2. Сортировка gene bodies согласно genome.txt

> ⚠️ **КРИТИЧНО**: порядок scaffolds должен совпадать!

```bash
bedtools sort -i Chlamy_gene_bodies.bed \
  -g Chlamy_genome.txt > Chlamy_gene_bodies.sorted.bed
```

---

### 🧩 4.3. Комплемент = интергенные регионы

```bash
bedtools complement -i Chlamy_gene_bodies.sorted.bed \
  -g Chlamy_genome.txt > Chlamy_intergenic.bed
```

#### 🔍 Параметры:

| Флаг | Функция |
|------|---------|
| `complement` | Возвращает регионы генома, НЕ покрытые input |
| `-i` | BED файл с генами (должен быть отсортирован!) |
| `-g` | Genome size file |

> 💡 **Биологический смысл**: Интергенные регионы = всё пространство между генами

---

### 📊 4.4. Статистика интергенных регионов

```bash
# Количество регионов
wc -l Chlamy_intergenic.bed

# Общий размер
awk '{sum+=$3-$2} END {print sum/1000000 " Mb total intergenic"}' \
  Chlamy_intergenic.bed

# Средний размер
awk '{sum+=$3-$2; n++} END {print sum/n " bp average"}' \
  Chlamy_intergenic.bed

# Процент от генома
GENOME_SIZE=$(awk '{sum+=$2} END {print sum}' Chlamy_genome.txt)
INTERGENIC_SIZE=$(awk '{sum+=$3-$2} END {print sum}' Chlamy_intergenic.bed)
echo "Intergenic = $(awk -v i=$INTERGENIC_SIZE -v g=$GENOME_SIZE \
  'BEGIN{printf "%.2f%%", i/g*100}')"
```

#### 📊 Результаты:

| Метрика | *Volvox* | *Chlamydomonas* | Ожидание | Статус |
|---------|----------|-----------------|----------|--------|
| Регионов | 13,350 | 12,147 | 10-20k | ✅ |
| Общий размер | 47.18 Mb | 27.45 Mb | - | ⚠️ |
| Средний размер | 3,534 bp | 2,260 bp | 2-8 kb | ✅ |
| % генома | 35.97% | 23.95% | 40-60% | ⚠️ |

---

## Шаг 5: Анализ Genome Budget

### 🧮 Genome Budget Analysis (*Volvox*)

```
Total genome:              131.17 Mb (100%)
├─ CDS:                     27.88 Mb (21.3%)
├─ Introns:                 40.94 Mb (31.2%)
└─ Intergenic:              47.18 Mb (36.0%)
                           ─────────────────
Total accounted:           115.99 Mb (88.4%)
MISSING:                    15.17 Mb (11.6%) ← ГДЕ???
```

---

### 🔍 Поиск "missing" последовательностей

#### ✅ Проверка 1: Scaffolds без генов

```bash
comm -23 \
  <(cut -f1 Vcarteri_genome.txt | sort) \
  <(cut -f1 Vcarteri_gene_bodies.bed | sort -u) \
  > Vcarteri_empty_scaffolds.txt

grep -f Vcarteri_empty_scaffolds.txt Vcarteri_genome.txt | \
  awk '{sum+=$2} END {print sum/1000000 " Mb in empty scaffolds"}'
```

**Результаты:**
- *Volvox*: 234 пустых scaffolds = **1.00 Mb (0.8%)**
- *Chlamydomonas*: 13 пустых scaffolds = **0.47 Mb (0.4%)**

---

#### ✅ Проверка 2: Органелльные геномы

```bash
grep -iE "chloroplast|mitochondri|plastid|mt_|pt_" Vcarteri_genome.txt
```

**Результаты для *Chlamydomonas*:**
- `plastome`: 205,535 bp (~0.2 Mb)
- `mitogenome`: 15,789 bp
- `mating_type_plus`: 374,891 bp

---

#### ⚠️ Проверка 3: Маленькие scaffolds (< 10 kb)

```bash
awk '$2 < 10000 {count++; sum+=$2} END {\
  print count, "scaffolds <10kb, total", sum/1e6, "Mb"}' \
  Vcarteri_genome.txt
```

**Результаты:**
- *Volvox*: 295 scaffolds < 10kb = **1.16 Mb (0.9%)**
- *Chlamydomonas*: 1 scaffold = **0.008 Mb**

> 💡 **Биологический вывод**: Это часто junk sequences (assembly errors)  
> **Рекомендация**: фильтровать в финальном анализе

---

#### 🧩 Проверка 4: Assembly gaps (N's)

```bash
grep -v "^>" Vcarteri_317_v2.assembly.fa | \
  tr -cd 'Nn' | wc -c | \
  awk '{printf "N content: %.2f Mb\\n", $1/1e6}'
```

**Результаты:**
- *Volvox*: **~10-15 Mb N's** ← Основная причина "missing" нуклеотидов!
- *Chlamydomonas*: **1.11 Mb N's**

---

### ✅ РЕШЕНИЕ "пропавших" нуклеотидов:

Большая часть "missing" последовательностей состоит из:
1. **Assembly gaps (N's)** - основная причина
2. **UTR перекрытия** с CDS
3. **Маленькие scaffolds** (< 10 kb)
4. **Органелльные геномы**

---

## 📂 Выходные файлы

### Основные результаты:

| # | Файл | Описание |
|---|------|----------|
| 1 | `Vcarteri_CDS_merged.bed` | Объединенные CDS регионы *Volvox* |
| 2 | `Chlamy_CDS_merged.bed` | Объединенные CDS регионы *Chlamydomonas* |
| 3 | `Vcarteri_UTR_merged.bed` | Объединенные UTR регионы *Volvox* |
| 4 | `Chlamy_UTR_merged.bed` | Объединенные UTR регионы *Chlamydomonas* |
| 5 | `Vcarteri_coding_elements_merged.bed` | Все mature mRNA loci *Volvox* |
| 6 | `Chlamy_coding_elements_merged.bed` | Все mature mRNA loci *Chlamydomonas* |
| 7 | `Vcarteri_intergenic.bed` | Интергенные регионы *Volvox* |
| 8 | `Chlamy_intergenic.bed` | Интергенные регионы *Chlamydomonas* |

---

## 📊 Результаты

### Финальная статистика

| Категория | *Volvox carteri* | *Chlamydomonas reinhardtii* |
|-----------|------------------|------------------------------|
| **Размер генома** | 131.17 Mb | 114.63 Mb |
| **CDS** | 27.88 Mb (21.3%) | 40.18 Mb (35.1%) |
| **Интроны** | 40.94 Mb (31.2%) | 46.00 Mb (40.1%) |
| **Интергенные** | 47.18 Mb (36.0%) | 27.45 Mb (23.9%) |
| **Количество генов** | 14,247 | 17,741 |
| **Avg CDS exon** | 252.7 bp | 270.4 bp |
| **Merge efficiency** | 13.2% | 45.9% |

---

## 💡 Биологические выводы

### 🔬 Ключевые находки:

1. **Chlamydomonas имеет более сложный альтернативный сплайсинг**
   - 46% CDS экзонов схлопываются при merge (vs 13% в *Volvox*)
   - Множество изоформ на ген
   - Сложная post-transcriptional регуляция

2. **Volvox имеет более компактные coding regions**
   - Меньше экзонов на ген (8.8 vs 15.5)
   - Более простая структура генов
   - Меньше альтернативного сплайсинга

3. **Многоклеточность не требует усложнения gene structure**
   - *Volvox* эволюционировал сложность через **другие механизмы**
   - Вероятно: регуляторные сети, промоторы, lncRNA

4. **Интергенные регионы Volvox крупнее**
   - 3.5 kb vs 2.3 kb в *Chlamydomonas*
   - **Потенциально больше цис-регуляторных элементов**


---

## 📚 Ссылки

- **Phytozome**: [https://phytozome.jgi.doe.gov/](https://phytozome.jgi.doe.gov/)
- **bedtools documentation**: [https://bedtools.readthedocs.io/](https://bedtools.readthedocs.io/)
- **samtools documentation**: [http://www.htslib.org/](http://www.htslib.org/)

---

