下面是一份你可以直接放到 GitHub 的简单 `README.md`，风格比较简洁、应付审稿人也够用 😄
（你自己把文件名 / 路径里 HDOM/LDOM 改成通用描述就行）

````markdown
# Reactomics PMD Annotation Script

This repository contains a Python script for simple reactomics analysis based on FT-ICR-MS molecular formula data.  
The workflow calculates paired mass differences (PMDs) between two samples and annotates potential reaction types and KEGG enzyme classes.

## 1. Overview

The script:

1. Reads two FT-ICR-MS result files (e.g. `0-LDOM-0.xlsx` and `0-LDOM-14.xlsx`).
2. Parses molecular formulas (`sumFormula`) into elemental counts of C, H, N, O, S.
3. Calculates all pairwise elemental differences (ΔC, ΔH, ΔN, ΔO, ΔS) between sample 0 and sample 14.
4. Counts the frequency of each unique PMD pattern.
5. Classifies PMDs into reaction types (e.g. reduction, oxidation, alkylation, decarboxylation, etc.).
6. Suggests possible KEGG enzyme classes (EC number, enzyme category, brief description).
7. Filters out unmatched PMDs (no EC assigned).
8. Exports the annotated PMD table as an Excel file.

## 2. Requirements

- Python 3.8+  
- Required packages:
  - `pandas`
  - `openpyxl` (for reading/writing Excel files)

You can install them via:

```bash
pip install pandas openpyxl
````

## 3. Input Data

The script expects two Excel files with at least one column:

* `sumFormula` – molecular formula strings, e.g. `C10H12O3`, `C15H18O5N1`, etc.

In the example code, the input files are:

```python
sample1.xlsx
sample2.xlsx
```

You can change these paths to your own FT-ICR-MS result files.

## 4. How to Run

1. Edit the script and update the input file paths:

```python
df_hdom0 = pd.read_excel(r"PATH/TO/sample1.xlsx")
df_hdom14 = pd.read_excel(r"PATH/TO/sample2.xlsx")
```

2. Run the script:

```bash
python reactomics_pmd_annotation.py
```

3. After successful execution, an output file will be generated:

```text
Reactomics_Final_Annotated.xlsx
```

## 5. Output Description

The output Excel file contains one row per **unique PMD** (elemental difference), with the following fields:

* `Frequency` – how many times this PMD occurs between the two samples
* `ΔC`, `ΔH`, `ΔN`, `ΔO`, `ΔS` – elemental differences between sample 0 and sample 14
* `反应类型` – inferred reaction type in Chinese (e.g. `+2H（还原）`, `-CO2（脱羧/呼吸）`)
* `EC编号` – suggested KEGG EC class (e.g. `EC 1.3.1.-`, `EC 4.1.1.-`)
* `酶类别` – enzyme category (e.g. oxidoreductase, transferase, lyase)
* `功能描述` – brief English description of the enzyme function

Rows where no EC class could be assigned (`EC编号 = "NA"`) are removed before export.

## 6. Notes & Limitations

* The mapping from PMD to reaction type and KEGG class is **rule-based and approximate**; it is intended for exploratory interpretation, not for strict mechanistic proof.
* Only C, H, N, O, and S elements are considered.
* All pairwise combinations between the two samples are used, so runtime and memory usage grow with the number of formulas.

## 7. Citation / Acknowledgement

If you use or adapt this script in a publication, you can simply acknowledge it as:

> Reactomics PMD analysis was performed using an in-house Python script for PMD-based reaction and enzyme class annotation.
