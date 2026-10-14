# Rebuild Application result CSVs as bordered HTML tables (wide form where clearer).
# Run: python Code/Application/format_app_tables.py

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(r"F:/ising model thesis/-Fourth-Corner-Ising-Model-")
TAB = ROOT / "Application_SOCPR" / "tables"
OUT = TAB / "app_tables_formatted.html"

ENVS = ["(const)", "salinity", "water_temp", "PAR"]
GROUPS = [
    "Amphipod",
    "Appendicularian",
    "Chaetognath",
    "Copepod",
    "Foraminifera",
    "Krill",
    "Other",
    "Pteropod",
]


def read_csv(name: str) -> list[dict[str, str]]:
    with (TAB / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def fmt_num(x, digits=4, keep_int=False):
    if x is None or x == "":
        return ""
    if isinstance(x, str) and x in ("TRUE", "FALSE"):
        return "是" if x == "TRUE" else "否"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if keep_int or abs(v - round(v)) < 1e-12:
        if abs(v) >= 1 or v == 0:
            return str(int(round(v)))
    if abs(v) > 0 and abs(v) < 1e-3:
        return f"{v:.4g}"
    return f"{v:.{digits}f}"


def td(val, cls=""):
    c = f' class="{cls}"' if cls else ""
    return f"<td{c}>{val}</td>"


def th(val):
    return f"<th>{val}</th>"


def table(headers, rows, caption, note=""):
    head = "".join(th(h) for h in headers)
    body = []
    for row in rows:
        cells = "".join(td(c) for c in row)
        body.append(f"<tr>{cells}</tr>")
    note_html = f'<p class="note">{note}</p>' if note else ""
    return (
        f'<section class="block">'
        f"<h2>{caption}</h2>{note_html}"
        f'<div class="wrap"><table>'
        f"<thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody>"
        f"</table></div></section>"
    )


def pivot_lk(rows, model):
    grid = {e: {g: "" for g in GROUPS} for e in ENVS}
    sel = {e: {g: False for g in GROUPS} for e in ENVS}
    ref = {e: {g: False for g in GROUPS} for e in ENVS}
    for r in rows:
        if r["model"] != model:
            continue
        grid[r["env"]][r["group"]] = r["estimate"]
        sel[r["env"]][r["group"]] = r["selected"] == "TRUE"
        if r.get("is_reference_group") == "TRUE":
            ref[r["env"]][r["group"]] = True
    out = []
    for e in ENVS:
        cells = [e]
        for g in GROUPS:
            lab = fmt_num(grid[e][g], 4)
            if ref[e][g]:
                lab = f"{lab}†"
            elif not sel[e][g]:
                lab = "0"
            cells.append(lab)
        out.append(cells)
    return out


def pivot_beta(rows):
    species = []
    seen = set()
    cell = {}
    meta = {}
    for r in rows:
        sp = r["species"]
        if sp not in seen:
            seen.add(sp)
            species.append(sp)
            meta[sp] = (r["group"], r["prevalence"])
        cell[(sp, r["env"])] = (r["estimate"], r["selected"] == "TRUE")
    out = []
    for sp in species:
        grp, prev = meta[sp]
        row = [sp, grp, fmt_num(prev, 4)]
        for e in ENVS:
            est, selected = cell[(sp, e)]
            row.append(fmt_num(est, 4) if selected else "0")
        out.append(row)
    return out


overview = read_csv("app_model_overview.csv")
A = read_csv("app_A_fourth_corner.csv")
B = read_csv("app_B_fourth_corner.csv")
beta = read_csv("app_FCIR_I_Beta.csv")
sg = read_csv("app_Theta_same_vs_cross_group.csv")
yr = read_csv("app_species_recording_windows.csv")
edges = read_csv("app_FCIR_M_Theta_edges.csv")

sections = []

sections.append(
    table(
        ["模型", "主效应", "交互", "参数数", "非零数", "非零比例(%)", "λ_min", "glm 别名数", "用时(s)", "准完全分离"],
        [
            [
                r["model"],
                r["main_effect"],
                r["interaction"],
                fmt_num(r["n_params"], keep_int=True),
                fmt_num(r["n_nonzero"], keep_int=True),
                fmt_num(r["pct_nonzero"], 1),
                fmt_num(r["lambda_min"], 6),
                fmt_num(r["glm_aliased"], keep_int=True),
                fmt_num(r["elapsed_s"], 1),
                "是" if r["quasi_separation"] == "TRUE" else "否",
            ]
            for r in overview
        ],
        "表 1  三个模型拟合概览",
        "来源：app_model_overview.csv。N = 29,983。数值按原始 CSV 四舍五入显示。",
    )
)

lk_headers = ["环境"] + GROUPS
sections.append(
    table(
        lk_headers,
        pivot_lk(B, "FCIR"),
        "表 2a  第四角主效应 B（FCIR）",
        "来源：app_B_fourth_corner.csv。格子为估计值；0 表示 lasso 恰好压成 0。"
        "† Pteropod 是参照组，B[, K] ≡ 0，不是被选掉的。",
    )
)
sections.append(
    table(
        lk_headers,
        pivot_lk(B, "FCIR_M"),
        "表 2b  第四角主效应 B（FCIR_M）",
        "来源：app_B_fourth_corner.csv。与表 2a 同一套长表，按模型拆成宽表。",
    )
)
sections.append(
    table(
        lk_headers,
        pivot_lk(A, "FCIR"),
        "表 3a  第四角交互 A（FCIR）",
        "来源：app_A_fourth_corner.csv。A 只作用于跨功能群配对；同群 99 对由 α₀ 吸收。",
    )
)
sections.append(
    table(
        lk_headers,
        pivot_lk(A, "FCIR_I"),
        "表 3b  第四角交互 A（FCIR_I）",
        "来源：app_A_fourth_corner.csv。",
    )
)

sections.append(
    table(
        ["物种", "功能群", "出现率", "(const)", "salinity", "water_temp", "PAR"],
        pivot_beta(beta),
        "表 4  FCIR_I 物种级主效应 β_j",
        "来源：app_FCIR_I_Beta.csv。0 表示恰好为零。出现率 = 该物种在 29,983 条记录中的出现比例。",
    )
)

sections.append(
    table(
        ["配对类型", "总边数", "被选中", "选中比例(%)", "|θ| 中位数", "θ 均值"],
        [
            [
                "跨功能群" if r["same_group"] == "FALSE" else "同功能群",
                fmt_num(r["n_pairs"], keep_int=True),
                fmt_num(r["n_selected"], keep_int=True),
                fmt_num(r["pct_selected"], 1),
                fmt_num(r["median_abs_theta"], 3),
                fmt_num(r["mean_theta"], 3),
            ]
            for r in sg
        ],
        "表 5  Θ_int 同群 vs 跨群",
        "来源：app_Theta_same_vs_cross_group.csv。",
    )
)

sections.append(
    table(
        ["物种", "功能群", "出现率", "被记录年数", "首年", "末年"],
        [
            [
                r["species"],
                r["group"],
                fmt_num(r["prevalence"], 4),
                fmt_num(r["n_years_recorded"], keep_int=True),
                fmt_num(r["first_year"], keep_int=True),
                fmt_num(r["last_year"], keep_int=True),
            ]
            for r in yr
        ],
        "表 6  各物种记录窗口",
        "来源：app_species_recording_windows.csv。调查年份共 24 年（1991–2016）。",
    )
)

edge_rows = []
for i, r in enumerate(edges, 1):
    edge_rows.append(
        [
            i,
            r["sp1"],
            r["sp2"],
            r["grp1"],
            r["grp2"],
            "同群" if r["same_group"] == "TRUE" else "跨群",
            fmt_num(r["theta"], 4),
            fmt_num(r["prev1"], 4),
            fmt_num(r["prev2"], 4),
            fmt_num(r["n11"], keep_int=True),
            fmt_num(r["n10"], keep_int=True),
            fmt_num(r["n01"], keep_int=True),
            fmt_num(r["n00"], keep_int=True),
            fmt_num(r["min_cell"], keep_int=True),
            "是" if r["selected"] == "TRUE" else "否",
            fmt_num(r["years_shared"], keep_int=True),
            fmt_num(r["years_jaccard"], 3),
        ]
    )

sections.append(
    table(
        [
            "#",
            "物种 1",
            "物种 2",
            "群 1",
            "群 2",
            "配对",
            "θ",
            "出现率1",
            "出现率2",
            "n11",
            "n10",
            "n01",
            "n00",
            "min cell",
            "选中",
            "共有年数",
            "年份 Jaccard",
        ],
        edge_rows,
        "表 7  FCIR_M 全部 435 条 Θ_int 边（按 |θ| 降序）",
        "来源：app_FCIR_M_Theta_edges.csv。n11 = 共同出现次数；min cell = 2×2 表最小格，0 表示完全分离。"
        "年份 Jaccard = 两物种被记录年份集合的交/并。",
    )
)

html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>南大洋 CPR 应用结果表</title>
<style>
  body {{
    font-family: "Times New Roman", "SimSun", serif;
    margin: 24px auto;
    max-width: 1400px;
    color: #111;
    background: #fff;
  }}
  h1 {{ font-size: 22px; margin: 0 0 8px; }}
  .lead {{ color: #333; margin: 0 0 20px; font-size: 14px; }}
  nav {{
    border: 1px solid #222;
    padding: 10px 14px;
    margin-bottom: 28px;
    background: #f7f7f7;
  }}
  nav a {{ margin-right: 14px; color: #0645ad; }}
  .block {{ margin: 0 0 36px; page-break-inside: avoid; }}
  h2 {{ font-size: 16px; margin: 0 0 8px; }}
  .note {{ font-size: 12px; color: #333; margin: 0 0 8px; }}
  .wrap {{ overflow-x: auto; }}
  table {{
    border-collapse: collapse;
    width: 100%;
    font-size: 12.5px;
    border: 1.5px solid #111;
  }}
  th, td {{
    border: 1px solid #111;
    padding: 4px 7px;
    vertical-align: middle;
  }}
  th {{
    background: #e8e8e8;
    font-weight: 700;
    text-align: center;
  }}
  td {{ text-align: center; }}
  td:first-child, th:first-child {{ text-align: left; }}
  tbody tr:nth-child(even) {{ background: #f4f4f4; }}
  @media print {{
    nav {{ display: none; }}
    body {{ margin: 12px; }}
    table {{ font-size: 10px; }}
  }}
</style>
</head>
<body>
<h1>南大洋 CPR 应用结果：原始估计值表</h1>
<p class="lead">
由 Application_SOCPR/tables 中的 CSV 重排而成。
第四角矩阵改成环境 × 功能群宽表；β_j 改成物种 × 环境宽表；Θ_int 保留全部 435 条边。
0 表示 adaptive lasso 恰好压成零，不是四舍五入。
</p>
<nav>
  <a href="#t1">概览</a>
  <a href="#t2a">B-FCIR</a>
  <a href="#t2b">B-FCIR_M</a>
  <a href="#t3a">A-FCIR</a>
  <a href="#t3b">A-FCIR_I</a>
  <a href="#t4">β_j</a>
  <a href="#t5">同群/跨群</a>
  <a href="#t6">记录窗口</a>
  <a href="#t7">Θ_int 全表</a>
</nav>
{''.join(sections)}
<script>
  const caps = document.querySelectorAll('section.block');
  const ids = ['t1','t2a','t2b','t3a','t3b','t4','t5','t6','t7'];
  caps.forEach((el, i) => {{ if (ids[i]) el.id = ids[i]; }});
</script>
</body>
</html>
"""

OUT.write_text(html, encoding="utf-8")
print(f"Wrote {OUT}")
