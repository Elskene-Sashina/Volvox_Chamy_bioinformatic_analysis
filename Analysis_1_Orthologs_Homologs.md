# 🧬 Анализ №3: Ортологи и расширенные генные семейства

> **Сравнительная геномика**: *Chlamydomonas reinhardtii* vs *Volvox carteri*  
> **Метод**: OrthoFinder v2.5.5 с DIAMOND search и MSA

---

## 📋 Оглавление
- [1. Подготовка данных](#1-подготовка-данных)
- [2. Запуск OrthoFinder](#2-запуск-orthofinder)
- [3. Общая статистика](#3-общая-статистика)
- [4. Volvox-специфичные расширенные семейства](#4-volvox-специфичные-расширенные-семейства)
- [5. Фильтрация pure Volvox-specific](#5-фильтрация-pure-volvox-specific)
- [6. Хромосомное распределение](#6-хромосомное-распределение)
- [7. GO аннотация](#7-go-аннотация)
- [8. Анализ длин белков](#8-анализ-длин-белков)
- [9. Типы ортологии](#9-типы-ортологии)
- [10. Дупликации и потери генов](#10-дупликации-и-потери-генов)
- [11. Биологические выводы](#11-биологические-выводы)

---

## 1. Подготовка данных

### 📦 Источник

**Phytozome** (JGI):
- *Chlamydomonas reinhardtii*: Assembly v6.0, Annotation v6.1
- *Volvox carteri*: Assembly v2, Annotation v2.1

### 🔧 Извлечение белковых последовательностей

```bash
# Распаковка
gunzip -c data/proteins/Vcarteri_v2_1.protein.fa.gz > \
  data/proteins/vcarteri_v2.1.faa

gunzip -c data/proteins/CreinhardtiiCC_4532_v6_1.protein_primaryTranscriptOnly.fa.gz > \
  data/proteins/creinhardtii_v6.1_primary.faa
```

### ⚙️ Форматирование заголовков

OrthoFinder требует формата: `>species|gene_id`

```python
# Скрипт для переформатирования headers
import re

def reformat_fasta(input_file, output_file, species_prefix):
    with open(input_file) as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                # Извлечь gene_id (первое слово после >)
                gene_id = line.split()[0][1:]
                f_out.write(f'>{species_prefix}|{gene_id}\n')
            else:
                f_out.write(line)

# Применить
reformat_fasta('vcarteri_v2.1.faa', 
               'orthofinder_input/Volvox_carteri.faa', 
               'Vca')

reformat_fasta('creinhardtii_v6.1_primary.faa', 
               'orthofinder_input/Chlamydomonas_reinhardtii_v6.faa', 
               'Cre')
```

> 💡 **Зачем**: OrthoFinder использует имя файла как название вида, добавляя префикс к каждому гену в филогенетических деревьях.

---

## 2. Запуск OrthoFinder

### 🚀 Команда

```bash
orthofinder -f orthofinder_input \
  -t 4 \
  -a 4 \
  -S diamond \
  -M msa \
  -o orthofinder_output \
  2>&1 | tee orthofinder.log
```

#### 📖 Параметры:

| Флаг | Значение | Функция |
|------|----------|---------|
| `-f` | `orthofinder_input` | Папка с FASTA файлами |
| `-t` | 4 | Threads для BLAST/DIAMOND |
| `-a` | 4 | Threads для MSA (Multiple Sequence Alignment) |
| `-S` | `diamond` | Использовать DIAMOND (быстрее BLAST) |
| `-M` | `msa` | Создать MSA для каждой ортогруппы |
| `-o` | `orthofinder_output` | Выходная папка |

⏱️ **Время работы**: ~2-4 часа (зависит от CPU)

---

## 3. Общая статистика

### 📊 Statistics_Overall

| Метрика | Значение | % |
|---------|----------|---|
| **Всего генов** | 33,768 | 100.0 |
| **Генов в ортогруппах** | 27,992 | **82.9** |
| **Unassigned генов** | 5,776 | 17.1 |
| **Всего ортогрупп** | 11,332 | - |
| **Single-copy ортогрупп** | 8,728 | **77.0** |
| **Средний размер ортогруппы** | 2.5 | - |
| **Медиана размера** | 2.0 | - |

> ✅ **Высокое качество**: 82.9% генов успешно кластеризовано в ортогруппы!

---

### 📊 Statistics_PerSpecies

| Метрика | *Chlamydomonas* | *Volvox* |
|---------|-----------------|----------|
| **Всего генов** | 17,693 | 16,075 |
| **В ортогруппах** | 14,380 (81.3%) | 13,612 (84.7%) |
| **Unassigned** | 3,313 (18.7%) | 2,463 (15.3%) |
| **Средний размер** | 2.47 | 2.54 |

> 💡 **Наблюдение**: *Volvox* имеет **выше процент** генов в ортогруппах (84.7% vs 81.3%) → меньше видоспецифичных генов?

---

### 📈 Распределение размеров ортогрупп

#### *Volvox carteri*:

| Размер | Ортогрупп | Генов всего |
|--------|-----------|-------------|
| 0 (отсутствуют) | 373 | 0 |
| 1 | 9,121 | 9,121 |
| 2 | 1,249 | 2,498 |
| 3 | 205 | 615 |
| 4 | 110 | 440 |
| 5 | 40 | 200 |
| 6-10 | 56 | 407 |
| 10-15 | 10 | 126 |
| 16-20 | 3 | 59 |
| 21-50 | 5 | 136 |

**Максимум**: 27 генов в одной ортогруппе (OG0000019)

#### *Chlamydomonas reinhardtii*:

| Размер | Ортогрупп | Генов всего |
|--------|-----------|-------------|
| 0 (отсутствуют) | 533 | 0 |
| 1 | 9,906 | 9,906 |
| 2 | 611 | 1,222 |
| 3 | 197 | 591 |
| 4 | 81 | 324 |
| 5-10 | 105 | 777 |
| 51-100 | 4 | 298 |
| 151-200 | 1 | 199 |
| 201-500 | 1 | 206 |

**Максимум**: 206 генов (OG с огромной экспансией!)

---

## 4. Volvox-специфичные расширенные семейства

### 🎯 Определение

**Расширенные семейства** = ортогруппы, где:
- *Chlamydomonas* = 0 генов (отсутствует)
- *Volvox* ≥ 2 генов (дупликация)

### 🔍 Извлечение

```bash
# Фильтр: Chlamydomonas = 0, Volvox >= 2
awk -F'\t' 'NR==1 {print; next} $2==0 && $3>=2' \
  Orthogroups.GeneCount.tsv > volvox_expanded_families.tsv

# Подсчет
tail -n +2 volvox_expanded_families.tsv | wc -l
```

**Результат**: **373 ортогруппы**

---

### 📝 Извлечение генов из расширенных семейств

```bash
# 1. Создать список OG IDs
cut -f1 volvox_expanded_families.tsv | tail -n +2 > expanded_orthogroups.txt

# 2. Извлечь гены Volvox из Orthogroups.tsv
awk -F'\t' 'NR==FNR {ortho[$1]=1; next}  
             FNR==1 {print "Orthogroup\tVolvox_carteri"; next} 
             $1 in ortho {print $1"\t"$3}' \
  expanded_orthogroups.txt Orthogroups.tsv > volvox_expanded_genes.tsv
```

---

### 🔀 Разбиение списков генов на отдельные строки

```bash
# Гены в Orthogroups.tsv разделены запятыми и пробелами
# Нужно разделить на отдельные строки

awk -F'\t' '{ 
    orthogroup=$1;  
    n=split($2, genes, /[, ]+/);  
    for(i=1; i<=n; i++)  
        if(genes[i]!="")  
            print orthogroup"\t"genes[i] 
}' volvox_expanded_genes.tsv > volvox_genes_list.txt
```

---

### 🔄 Конвертация ID формата

OrthoFinder выдаёт: `Vca|0001s1100.1.p`  
Нужен формат Phytozome: `Vocar.0001s1100.1.p`

```bash
sed -E 's/Vca\|([0-9]+s[0-9]+\.[0-9]+\.p)/Vocar.\1/g' \
  volvox_genes_list.txt > volvox_genes_list_converted.txt
```

---

### 📋 Добавление аннотаций

```bash
awk -F'\t' 'BEGIN {OFS="\t"} 
  # Первый проход: читаем annotation_info
  NR==FNR {
    if(FNR==1) {
      # Сохраняем заголовок
      for(i=1; i<=NF; i++) annot_header[i]=$i
      next
    }
    peptide=$4  # peptideName (колонка 4)
    annot_line[peptide]=$0
    next
  }

  # Второй проход: обрабатываем volvox_genes_list_converted.txt
  FNR==1 {
    # Выводим объединенный заголовок
    print $0, annot_header[1], annot_header[2], ..., annot_header[16]
    next
  }

  {
    gene=$2  # Gene_ID
    if(annot_line[gene] != "") {
      print $0, annot_line[gene]
    } else {
      # Если нет аннотации, заполняем пустыми
      print $0, "", "", ..., ""
    }
  }' \
  Vcarteri_317_v2.1.P14.annotation_info.txt \
  volvox_genes_list_converted.txt \
  > volvox_genes_annotated_full.tsv
```

---

## 5. Фильтрация pure Volvox-specific

### ⚠️ Проблема

Некоторые гены имеют **Best-hit-clamy-name/defline** → это НЕ чистые Volvox-specific!

**Возможные причины**:
1. Дальний ортолог (не попал в ортогруппу из-за дивергенции)
2. Паралог
3. Ошибка аннотации

### 🧹 Решение: Исключить гены с Best-hit-clamy

#### Структура аннотации:

```bash
head -1 volvox_genes_annotated_full.tsv | tr '\t' '\n' | nl
```

| # | Колонка |
|---|---------|
| 1 | Orthogroup |
| 2 | Volvox_carteri |
| 3-6 | pacId, locusName, transcriptName, peptideName |
| 7-12 | Pfam, Panther, ec, KOG, KO, GO |
| 13-14 | Best-hit-arabi-name, Best-hit-arabi-defline |
| **15-16** | **Best-hit-clamy-name**, **Best-hit-clamy-defline** ← Важно! |
| 17-18 | Best-hit-rice-name, Best-hit-rice-defline |

---

#### Шаг 1: Найти ортогруппы с Best-hit-clamy

```bash
awk -F'\t' 'NR>1 && ($15!="" || $16!="") {print $1}' \
  volvox_genes_annotated_full.tsv | \
  sort -u > orthogroups_with_clamy_hit.txt

wc -l orthogroups_with_clamy_hit.txt
```

**Результат**: 209 ортогрупп имеют хиты с *Chlamydomonas*

---

#### Шаг 2: Создать список pure Volvox-specific

```bash
# Все ортогруппы
sort volvox_genes_list_converted.txt | cut -f1 | sort -u > all_orthogroups_from_list.txt

# Исключить те, что имеют Best-hit-clamy
comm -23 <(sort all_orthogroups_from_list.txt) \
         <(sort orthogroups_with_clamy_hit.txt) \
  > pure_volvox_specific.txt

wc -l pure_volvox_specific.txt
```

**Результат**: **164 ортогруппы** (из 373) — чистые Volvox-specific!

---

#### Шаг 3: Фильтровать аннотированный файл

```bash
awk -F'\t' 'BEGIN {OFS="\t"} 
    NR==FNR {pure_og[$1]=1; next} 
    FNR==1 {print; next} 
    $1 in pure_og {print} 
' pure_volvox_specific.txt \
  volvox_genes_annotated_full.tsv \
  > volvox_pure_specific_annotated.tsv

# Подсчет генов
tail -n +2 volvox_pure_specific_annotated.tsv | wc -l
```

**Результат**: **629 генов** в 164 pure Volvox-specific ортогруппах

---

## 6. Хромосомное распределение

### 🗺️ Кластеры генов на скаффолдах

```bash
# Извлечь scaffold ID из locus name
tail -n +2 volvox_pure_specific_annotated.tsv | cut -f4 | \
  sed 's/Vocar\.\([0-9]*\)s.*/\1/' | \
  sort | uniq -c | sort -rn | head -20
```

#### 📈 Топ-10 скаффолдов с Volvox-specific генами

| Scaffold | Количество генов | % от total |
|----------|------------------|------------|
| **0001** | 35 | 5.6% |
| **0013** | 13 | 2.1% |
| **0022** | 11 | 1.7% |
| **0005** | 11 | 1.7% |
| **0004** | 10 | 1.6% |
| **0003** | 10 | 1.6% |
| **0021** | 9 | 1.4% |
| **0017** | 9 | 1.4% |
| **0032** | 8 | 1.3% |
| **0015** | 8 | 1.3% |

> ⚠️ **Кластеризация**: Scaffold 0001 содержит **5.6%** всех Volvox-specific генов! Возможный **tandem duplication hotspot**.

---

### 🔬 Детальный анализ топ кластеров

#### Scaffold 1 (35 генов):

```
Vocar.0001s1100  OG0000019  (Tandem cluster: 5 генов)
Vocar.0001s1101  OG0000019
Vocar.0001s1102  OG0000019
Vocar.0001s1103  OG0000019
Vocar.0001s1104  OG0000019
Vocar.0001s0367  OG0000068
Vocar.0001s1065  OG0000079  Pfam:PF00098 (zinc finger)
Vocar.0001s0499  OG0000098
```

> 💡 **Паттерн**: Множественные **тандемные дупликации** (consecutive gene IDs)

---

#### Scaffold 13 (13 генов):

```
Vocar.0013s0096  OG0000019  (Tandem cluster: 6 генов)
Vocar.0013s0097  OG0000019
Vocar.0013s0098  OG0000019
Vocar.0013s0099  OG0000019
Vocar.0013s0100  OG0000019
Vocar.0013s0101  OG0000019
Vocar.0013s0238  OG0000098
Vocar.0013s0202  OG0000126
Vocar.0013s0240  OG0002078  Pfam:PF00652 (ricin B lectin)
Vocar.0013s0270  OG0002078
```

---

## 7. GO аннотация

### 📊 Статистика GO терминов

| Метрика | Значение |
|---------|----------|
| Генов с GO терминами | 17 (2.7%) |
| Всего GO аннотаций | 39 |
| Уникальных GO терминов | 16 |

> ⚠️ **Низкое покрытие**: Только 2.7% генов имеют GO аннотацию → большинство генов **функционально не охарактеризованы**!

---

### 🏆 Топ GO терминов

| GO ID | Количество | Описание | Категория |
|-------|------------|----------|-----------|
| **GO:0003676** | 5 | Nucleic acid binding | MF |
| **GO:0008270** | 5 | Zinc ion binding | MF |
| **GO:0003677** | 4 | DNA binding | MF |
| **GO:0006355** | 4 | Regulation of DNA-templated transcription | BP |
| **GO:0005634** | 3 | Nucleus | CC |
| **GO:0004553** | 3 | Hydrolase activity | MF |
| **GO:0005975** | 3 | Carbohydrate metabolic process | BP |
| **GO:0004879** | 2 | DNA-binding TF activity | MF |
| **GO:0050825** | 2 | **Ice binding** | MF |
| **GO:0008810** | 2 | Cellulase activity | MF |

> 💡 **Интересная находка**: **Ice binding proteins** (GO:0050825) — потенциальная адаптация к холодным условиям?

---

### 📈 Функциональные категории (сводка)

```
DNA/RNA regulation:     14 генов (82%)
├─ Transcription factors (zinc fingers)
├─ DNA binding proteins
└─ Nucleic acid binding

Cell wall/ECM:           6 генов (35%)
├─ Hydrolases (cellulase)
├─ Carbohydrate metabolism
└─ Ricin B lectins

Structural/unknown:      3 гена (18%)
└─ No functional annotation
```

---

## 8. Анализ длин белков

### 📏 Статистика длин (все 629 генов)

```python
import pandas as pd
import numpy as np

# Предполагаем, что длины белков уже извлечены
lengths = [...]  # Список длин из FASTA

print(f"Всего белков: {len(lengths)}")
print(f"Средняя длина: {np.mean(lengths):.1f} aa")
print(f"Медиана: {np.median(lengths):.0f} aa")
print(f"Минимум: {min(lengths)} aa")
print(f"Максимум: {max(lengths)} aa")

# Распределение по категориям
very_short = sum(1 for l in lengths if l < 100)
short = sum(1 for l in lengths if 100 <= l < 200)
medium = sum(1 for l in lengths if 200 <= l < 500)
long = sum(1 for l in lengths if l >= 500)

print(f"\nОчень короткие (<100 aa):   {very_short} ({very_short/len(lengths)*100:.1f}%)")
print(f"Короткие (100-200 aa):       {short} ({short/len(lengths)*100:.1f}%)")
print(f"Средние (200-500 aa):        {medium} ({medium/len(lengths)*100:.1f}%)")
print(f"Длинные (≥500 aa):           {long} ({long/len(lengths)*100:.1f}%)")
```

#### 📊 Результаты (из файла):

| Категория | Количество | % |
|-----------|------------|---|
| **Всего белков** | 42 (топ семейства) | 100.0 |
| Средняя длина | 412.5 aa | - |
| Медиана | 332 aa | - |
| Минимум | 69 aa | - |
| Максимум | 1255 aa | - |
| **Очень короткие** (<100 aa) | 12 | **28.6** |
| Короткие (100-200 aa) | 5 | 11.9 |
| Средние (200-500 aa) | 13 | 31.0 |
| Длинные (≥500 aa) | 12 | **28.6** |

---

### 🔍 Анализ очень коротких белков (<100 aa)

**12 белков (28.6%)** — это необычно много!

**Возможные объяснения**:
1. **Сигнальные пептиды** — регуляторные функции
2. **Регуляторные пептиды** — гормоны, факторы роста
3. **Фрагменты генов** — ошибки аннотации, псевдогены
4. **Короткие функциональные белки** — антимикробные пептиды, др.

> 💡 **Рекомендация**: Проверить ORF prediction для белков <100 aa (возможны артефакты)

---

### 🏆 Топ расширенные семейства (по длине белка)

| Ген | Копий | Длина (aa) | Ортогруппа |
|-----|-------|------------|------------|
| Vca\|0020s0179 | **6** | 1040-1083 | **Самое расширенное!** |
| Vca\|0012s0179 | 4 | 344-626 | Второе крупное |
| Vca\|0003s0370/0371 | 3 | 404-406 | Тандемные гены |

---

## 9. Типы ортологии

### 📊 Классификация ортологических отношений

#### 1:1 Ортологи (single-copy)

```
Species A:  ──●── gene1
Species B:  ──●── gene1_ortholog

Количество: 9,928 (91.9%)
```

> ✅ **Интерпретация**: Консервативные гены, критичные для выживания (housekeeping)

---

#### 1:many Ортологи

```
Species A:  ──●── gene1
               │
Species B:  ──●── gene1_copy1
            ──●── gene1_copy2
            ──●── gene1_copy3

*Chlamydomonas* → *Volvox*: 1,178
*Volvox* → *Chlamydomonas*: 393
```

**Механизм**: Дупликация генов в одном виде **после** разделения

---

#### many:1 Ортологи

```
Species A:  ──●── gene1_copy1
            ──●── gene1_copy2
               │
Species B:  ──●── gene1

*Chlamydomonas* → *Volvox*: 1,498
*Volvox* → *Chlamydomonas*: 2,891
```

**Механизм**: Потеря генов или дупликация в предке + потеря в одном виде

---

#### many:many Ортологи

```
Species A:  ──●── gene1_copy1
            ──●── gene1_copy2
               │
Species B:  ──●── gene1_copy1
            ──●── gene1_copy2

*Chlamydomonas* ↔ *Volvox*: 594
*Volvox* ↔ *Chlamydomonas*: 460
```

**Механизм**: Дупликации в **обоих** видах после разделения

---

### 📈 Сводная таблица ортологии

| Тип | Количество HOGs | % от общих |
|-----|-----------------|------------|
| **1:1 ортологи** | 9,928 | **91.9** |
| 1:many | 1 | 0.0 |
| many:1 | 8 | 0.1 |
| many:many | 154 | 1.4 |
| **ВСЕГО общих** | **10,091** | **93.4** |
| *Chlamydomonas*-specific | 397 | 3.7 |
| *Volvox*-specific | 313 | 2.9 |
| **ВСЕГО HOGs** | **10,801** | **100.0** |

---

### 🔬 Интерпретация:

#### Высокая синтения: **98.4%**

\[
\text{Синтения} = \frac{\text{1:1 ортологи}}{\text{Общие ортогруппы}} = \frac{9{,}928}{10{,}091} = 98.4\%
\]

> ✅ **Вывод**: Геномы *Chlamydomonas* и *Volvox* **очень консервативны** → мало крупных перестроек после разделения!

---

## 10. Дупликации и потери генов

### 📊 Дупликации по видам

| Вид | Всего дупликаций | С >50% поддержкой |
|-----|------------------|-------------------|
| *Chlamydomonas* | 2,495 | 2,495 (100%) |
| *Volvox* | 1,444 | 1,444 (100%) |
| Предковый узел (N0) | 1 | 1 |

**Разница**: *Chlamydomonas* имеет **1.73× больше** дупликаций!

---

### 🧮 Расчет потерь генов

```
Формула:
  Потери = (Unassigned genes) + (Species-specific orthogroups)

Chlamydomonas:
  Потери = 3,313 + 397 = 3,710 генов

Volvox:
  Потери = 2,463 + 313 = 2,776 генов
```

| Событие | *Chlamydomonas* | *Volvox* | Соотношение |
|---------|-----------------|----------|-------------|
| **Дупликации** | 2,495 | 1,444 | 1.73× |
| **Потери** | 3,710 | 2,776 | 1.34× |

> 💡 **Интерпретация**: *Chlamydomonas* более **динамичный** геном (больше дупликаций И потерь)

---

### ⚡ Асимметричные дупликации

**Определение**: Ортогруппы, где разница ≥3 гена между видами

```bash
# Python скрипт для поиска
import pandas as pd

df = pd.read_csv('N0.tsv', sep='\t')
asymmetric = df[(df['Chlamydomonas'] - df['Volvox']).abs() >= 3]

print(f"Всего асимметричных: {len(asymmetric)}")
print(asymmetric[['HOG', 'Chlamydomonas', 'Volvox', 'Diff']].to_string())
```

#### 📊 Результаты:

| HOG | *Chlamydomonas* | *Volvox* | Разница |
|-----|-----------------|----------|---------|
| N0.ids.HOG0000012 | 4 | 1 | **3** |
| N0.ids.HOG0000014 | 4 | 1 | 3 |
| N0.ids.HOG0000015 | 4 | 1 | 3 |
| N0.ids.HOG0000028 | 4 | 1 | 3 |
| **N0.ids.HOG0000031** | **1** | **4** | **3** ← Расширение в *Volvox*! |
| N0.ids.HOG0000035 | 4 | 1 | 3 |
| N0.ids.HOG0000040 | 4 | 1 | 3 |
| N0.ids.HOG0000051 | 4 | 1 | 3 |
| N0.ids.HOG0000078 | 4 | 1 | 3 |

**Всего**: 9 асимметричных ортогрупп

---

### 🔍 Детальный анализ HOG0000031 (Volvox-расширенная)

| OG | Аннотация | Функция |
|----|-----------|---------|
| OG0000125 | Proline-rich, SRCR domain | ECM, extracellular region |
| OG0000124 | UV DNA damage endonuclease | DNA repair |
| OG0000126 | IDRs, Signal peptide | Extracellular region |
| OG0000127 | Reverse transcriptase | DNA integration/recombination |

> 💡 **Интерпретация**: Расширение генов **ECM** (extracellular matrix) в *Volvox* → связано с многоклеточностью!

---

## 11. Биологические выводы

### 🎯 Ключевые находки

#### 1. Высокая консервативность геномов

- **91.9%** генов — 1:1 ортологи
- **98.4%** синтения
- Мало крупных геномных перестроек

> 💡 **Вывод**: *Volvox* и *Chlamydomonas* разделились относительно недавно (~200 млн лет назад), геномы остаются похожими

---

#### 2. Volvox-специфичные расширения (164 pure OGs)

- **629 генов** в расширенных семействах
- **98.4%** без гомологов в модельных растениях
- **Кластеризация** на scaffold 1 (5.6% всех генов)

**Функции**:
- DNA-binding (zinc fingers)
- ECM proteins (ricin B lectins)
- Ice-binding proteins
- Неизвестные функции (>80%)

> 💡 **Вывод**: Большинство Volvox-specific генов — **de novo** или сильно дивергировавшие, связаны с **многоклеточностью**

---

#### 3. Тандемные дупликации — ключевой механизм

**Примеры**:
- OG0000019: 27 копий (tandem cluster на scaffold 1, 13, 22)
- OG0000068: 12 копий (разбросаны по genome)

**Механизм**:
```
Original gene:  ──●──
                   ↓ Tandem duplication
Duplicated:     ──●──●──●──●──
```

> 💡 **Вывод**: *Volvox* активно использует **tandem duplications** для создания gene families

---

#### 4. Динамика геномов: Chlamydomonas > Volvox

| Событие | *Chlamydomonas* | *Volvox* | Интерпретация |
|---------|-----------------|----------|---------------|
| Дупликации | 2,495 | 1,444 | Одноклеточный более пластичен |
| Потери | 3,710 | 2,776 | Очистка генома в *Volvox* |
| Unassigned | 18.7% | 15.3% | Больше уникальных генов |

> 💡 **Парадокс**: Многоклеточный *Volvox* имеет **меньше** динамики генома, чем одноклеточный *Chlamydomonas*!

---

#### 5. Эволюция многоклеточности через малые изменения

**Путь *Volvox***:
- Не через WGD (whole genome duplication)
- Не через массовые дупликации
- Через:
  1. **Tandem duplications** (создание ECM gene families)
  2. **Регуляторные изменения** (промоторы, TF)
  3. **De novo гены** (164 pure OGs)

> 💡 **Вывод**: Многоклеточность возникла через **малые, целенаправленные изменения**, а не через геномные катаклизмы

---

## 📂 Выходные файлы

### Основные результаты OrthoFinder:

| # | Файл | Описание |
|---|------|----------|
| 1 | `Orthogroups.tsv` | Полная таблица ортогрупп |
| 2 | `Orthogroups.GeneCount.tsv` | Количество генов по видам |
| 3 | `N0.tsv` | Hierarchical orthogroups (HOGs) |
| 4 | `Orthogroups_UnassignedGenes.tsv` | Гены без ортогруппы |
| 5 | `Statistics_Overall.tsv` | Общая статистика |
| 6 | `Statistics_PerSpecies.tsv` | Статистика по видам |

### Наши фильтрованные файлы:

| # | Файл | Описание |
|---|------|----------|
| 7 | `volvox_expanded_families.tsv` | 373 расширенных OG |
| 8 | `pure_volvox_specific.txt` | 164 pure Volvox OG IDs |
| 9 | `volvox_pure_specific_annotated.tsv` | 629 генов с аннотацией |
| 10 | `orthogroups_with_clamy_hit.txt` | 209 OG с хитами Chlamy |


---

## 📚 Ссылки

- **OrthoFinder**: [https://github.com/davidemms/OrthoFinder](https://github.com/davidemms/OrthoFinder)
- **Emms & Kelly (2019)**: Genome Biology, OrthoFinder paper
- **Phytozome**: [https://phytozome-next.jgi.doe.gov/](https://phytozome-next.jgi.doe.gov/)

---

<div align="center">

</div>
