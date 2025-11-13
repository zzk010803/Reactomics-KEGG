#  1. 导入必要的库
import pandas as pd
import re
from itertools import product
from collections import Counter

#  2. 读取 HDOM-0 和 HDOM-14 的 FT-ICR-MS 数据
df_hdom0 = pd.read_excel(r"0-LDOM-0.xlsx")       # 修改为本地路径
df_hdom14 = pd.read_excel(r"0-LDOM-14.xlsx")         # 修改为本地路径

#  3. 分子式解析函数（返回元素计数）
def parse_formula(formula):
    elements = ['C', 'H', 'N', 'O', 'S']
    pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    counts = dict.fromkeys(elements, 0)

    text = str(formula).replace(" ", "")
    for m in pattern.finditer(text):
        elem = m.group(1)
        num_str = m.group(2)
        if elem in counts:
            counts[elem] += int(num_str) if num_str else 1

    return counts

#  4. 提取样本中存在的唯一分子式并解析元素
#hdom0_formulas = df_hdom0['sumFormula'].dropna().unique()
#hdom14_formulas = df_hdom14['sumFormula'].dropna().unique()

#hdom0_parsed = {f: parse_formula(f) for f in hdom0_formulas}
#hdom14_parsed = {f: parse_formula(f) for f in hdom14_formulas}

#  4. 直接对每一行分子式进行解析（不再按分子式去重）
hdom0_series = df_hdom0['sumFormula'].dropna()
hdom14_series = df_hdom14['sumFormula'].dropna()

hdom0_parsed_list = [parse_formula(f) for f in hdom0_series]
hdom14_parsed_list = [parse_formula(f) for f in hdom14_series]

#  5. PMD（元素差异）计算：HDOM0 vs. HDOM14 所有组合
pmd_results = []
for h_counts in hdom0_parsed_list:
    for l_counts in hdom14_parsed_list:
        diff = {e: h_counts[e] - l_counts[e] for e in ['C', 'H', 'N', 'O', 'S']}
        pmd_results.append(tuple(diff.values()))

#  6. 统计频率，转换为 DataFrame
pmd_counter = Counter(pmd_results)
pmd_df = pd.DataFrame(pmd_counter.items(), columns=["ΔCHNOS", "Frequency"])
pmd_df[["ΔC", "ΔH", "ΔN", "ΔO", "ΔS"]] = pd.DataFrame(pmd_df["ΔCHNOS"].tolist(), index=pmd_df.index)
pmd_df = pmd_df.drop(columns="ΔCHNOS")
top_pmd_df = pmd_df.sort_values("Frequency", ascending=False).copy()

#  7. 反应类型识别函数（增强逻辑）
def classify_reaction_final(row):
    c, h, n, o, s = row["ΔC"], row["ΔH"], row["ΔN"], row["ΔO"], row["ΔS"]

    if o == n == s == 0 and h > 0 and h % 2 == 0:
        return f"+{h}H（还原）"
    elif o == n == s == 0 and h < 0 and h % 2 == 0:
        return f"{h}H（氧化）"
    elif c > 0 and abs(h - 2 * c) <= 2 and o == 0:
        return f"+{abs(c)}CH₂（烷基化）"
    elif c < 0 and abs(h - 2 * c) <= 2 and o == 0:
        return f"{c}CH₂（脱烷基）"
    elif o > 0 and c == 0:
        return f"+{o}O（氧化/羟基化）"
    elif o < 0 and c == 0:
        return f"{o}O（脱氧）"
    elif n > 0 and abs(h - 2 * n) <= 2:
        return f"+{n}NH₂（胺化）"
    elif n < 0 and abs(h - 2 * n) <= 2:
        return f"{n}NH₂（脱胺）"
    elif c > 0 and o == 2 * c and h == 0:
        return f"+{c}COOH（羧基化）"
    elif c < 0 and o == 2 * c and h == 0:
        return f"{c}COOH（脱羧）"
    elif c == 1 and h == 2 and o == 1:
        return "+CH2O（醇/醛化）"
    elif c == 2 and h == 4 and o == 2:
        return "+C2H4O2（乙酸化）"
    elif c == 0 and h == 1 and o == 1:
        return "+OH（羟基化）"
    elif c == 1 and h == 3 and o == 1:
        return "+CH3O（甲氧基化）"
    elif c == 2 and h == 3 and o == 2:
        return "+C2H3O2（酯化/乙酰化）"
    elif c == 1 and h == 1 and o == 2:
        return "+CO2（羧基相关）"
    elif c == -1 and h == -1 and o == -2:
        return "-CO2（脱羧/呼吸）"
    elif c == 2 and h == 2 and o == 1:
        return "+C2H2O（酮基添加）"
    elif h == 0 and o == 3:
        return "+O3（臭氧化）"
    else:
        return "未知或复合反应"

#  8. 更精细的 KEGG 酶类推测函数
def suggest_enzyme_detailed(reaction_type):
    mapping = {
        "还原": ("EC 1.3.1.-", "氧化还原酶", "Oxidoreductase (on CH-CH group)"),
        "氧化": ("EC 1.1.1.-", "氧化还原酶", "Oxidoreductase (on CH-OH group)"),
        "羟基化": ("EC 1.14.13.-", "氧化还原酶", "Monooxygenase"),
        "烷基化": ("EC 2.1.1.-", "甲基转移酶", "Methyltransferase"),
        "脱烷基": ("EC 3.5.1.-", "水解酶", "Dealkylase"),
        "胺化": ("EC 2.6.1.-", "氨基转移酶", "Aminotransferase"),
        "脱胺": ("EC 3.5.4.-", "脱胺酶", "Deaminase"),
        "羧基化": ("EC 6.4.1.-", "羧化酶", "Carboxylase"),
        "脱羧": ("EC 4.1.1.-", "裂解酶", "Decarboxylase"),
        "醇/醛化": ("EC 1.1.1.-", "氧化还原酶", "Alcohol dehydrogenase / Aldehyde oxidase"),
        "乙酸化": ("EC 2.3.1.-", "酰基转移酶", "Acetyltransferase"),
        "羟基": ("EC 1.14.-.-", "单加氧酶", "Hydroxylase"),
        "甲氧基": ("EC 2.1.1.6", "甲基转移酶", "O-Methyltransferase"),
        "酯化": ("EC 2.3.1.-", "酰基转移酶", "Esterase / Acetyltransferase"),
        "呼吸": ("EC 4.1.1.-", "脱羧酶", "Decarboxylase"),
        "酮基添加": ("EC 1.2.1.-", "脱氢酶", "Keto group dehydrogenase"),
        "臭氧化": ("EC 1.13.11.-", "双加氧酶", "Dioxygenase"),
    }

    for keyword, (ec, category, desc) in mapping.items():
        if keyword in reaction_type:
            return ec, category, desc
    return "NA", "暂未匹配", "Unknown reaction"

#  9. 应用 KEGG 分类、反应类型注释
top_pmd_df["反应类型"] = top_pmd_df.apply(classify_reaction_final, axis=1)
top_pmd_df[["EC编号", "酶类别", "功能描述"]] = top_pmd_df["反应类型"].apply(
    lambda x: pd.Series(suggest_enzyme_detailed(x))
)
# 过滤掉无法匹配反应的 PMD（删除 NA）
top_pmd_df = top_pmd_df[top_pmd_df["EC编号"] != "NA"]

#  10. 保存最终结果
top_pmd_df.to_excel("Reactomics_Final_Annotated.xlsx", index=False)
print("分析完成，结果保存为：Reactomics_Final_Annotated.xlsx")
