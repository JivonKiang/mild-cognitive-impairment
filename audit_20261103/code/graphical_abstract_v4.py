#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
C-2 Graphical Abstract：v4 转译重定位版四闸审计框架
四闸：development stability -> independent test -> decision utility -> cell context
风格：Graphical_Abstract_豆包协作说明.md 规范（暖白底/低饱和青绿+珊瑚红+金/细箭头/无3D）
输出：figures/GraphicalAbstract_v4_nature_grid.png (300dpi) + .tiff (600dpi)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib.lines import Line2D

plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

OUT = r'D:/20241004_MCI/20261103_MCI_DCA_NRI/figures'
os.makedirs(OUT, exist_ok=True)

# ---- 色板（低饱和四色体系） ----
BG    = '#FBFAF7'   # 暖白
INK   = '#2B2E33'   # 深墨灰
GRAY  = '#7E8287'   # 辅助灰
TEAL  = '#4F8A85'   # 低饱和青绿：开发/稳定/预期
CORAL = '#C0776A'   # 低饱和珊瑚红：反转/边界
GOLD  = '#B3923C'   # 低饱和金：决策效用强调
HAIR  = '#D8D5CE'   # 发丝线

fig = plt.figure(figsize=(16, 9), dpi=100)
fig.patch.set_facecolor(BG)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 16); ax.set_ylim(0, 9)
ax.axis('off')
ax.set_facecolor(BG)

def txt(x, y, s, size=9, color=INK, weight='normal', ha='left', va='center', style='normal'):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, style=style, zorder=5)

def thin_arrow(x1, x2, y, color=GRAY, lw=1.0):
    ax.annotate('', xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=lw,
                                shrinkA=0, shrinkB=0, mutation_scale=9),
                zorder=4)

# ---- 标题 ----
txt(0.55, 8.42, 'A 12-GENE BLOOD SCORE UNDER A PREDEFINED TRANSPORTABILITY AUDIT',
    size=15.5, weight='semibold')

# ---- 四闸区域几何 ----
ZW, GAP = 3.30, 0.55
X0 = [0.50, 0.50 + ZW + GAP, 0.50 + 2 * (ZW + GAP), 0.50 + 3 * (ZW + GAP)]
TOP = 7.45

def gate_frame(i, label, title, source, accent):
    x0 = X0[i]
    ax.add_line(Line2D([x0, x0 + ZW], [TOP, TOP], color=accent, lw=1.2, zorder=3))
    txt(x0, TOP - 0.42, label, size=8, color=GRAY, weight='semibold')
    txt(x0, TOP - 0.80, title, size=13, weight='bold')
    txt(x0, TOP - 1.28, source, size=8.5, color=GRAY)

gate_frame(0, 'GATE 1', 'DEVELOPMENT', 'AddNeuroMed · blood expression microarray', TEAL)
gate_frame(1, 'GATE 2', 'INDEPENDENT TEST', 'GSE249477 · blood RNA-seq', CORAL)
gate_frame(2, 'GATE 3', 'DECISION UTILITY', 'decision curve analysis + NRI, two cohorts', GOLD)
gate_frame(3, 'GATE 4', 'CELL CONTEXT', 'single-cell immune localization', GRAY)

# ---- 闸间细箭头 ----
for i in range(3):
    thin_arrow(X0[i] + ZW + 0.06, X0[i + 1] - 0.06, 4.55, color=GRAY, lw=1.1)

# ================= GATE 1：registry -> locked score，稳定 =================
g1 = X0[0]
# 95-gene registry：淡灰点阵（10 x 4 示意，收窄避让标签）
rng = np.random.default_rng(12)
for r in range(4):
    for c in range(10):
        ax.add_patch(Circle((g1 + 0.30 + c * 0.26, 6.30 - r * 0.20), 0.035,
                            facecolor=HAIR, edgecolor='none', alpha=0.75, zorder=2))
txt(g1 + 0.30, 5.44, '95-gene registry', size=8, color=GRAY)
ax.annotate('', xy=(g1 + 0.62, 5.02), xytext=(g1 + 0.62, 5.30),
            arrowprops=dict(arrowstyle='-|>', color=GRAY, lw=0.9, mutation_scale=8))
# locked 12 条稳定表达线（青绿，整齐；压低高度避让上方箭头）
heights = np.array([0.56, 0.61, 0.58, 0.62, 0.59, 0.64, 0.60, 0.63, 0.58, 0.61, 0.59, 0.62])
for c, h in enumerate(heights):
    ax.add_patch(Rectangle((g1 + 0.30 + c * 0.215, 4.30), 0.115, h,
                           facecolor=TEAL, edgecolor='none', alpha=0.88, zorder=3))
ax.add_line(Line2D([g1 + 0.24, g1 + 3.05], [4.30, 4.30], color=GRAY, lw=0.8, zorder=2))
txt(g1 + 0.30, 4.05, 'locked 12-gene score', size=8, color=INK)
txt(g1, 3.35, 'stable signal', size=10, color=TEAL, weight='semibold')
txt(g1, 2.95, 'in the development system', size=9, color=INK)

# ================= GATE 2：10/12 mapped + 方向反转 =================
g2 = X0[1]
# 12 圆点：10 实（映射成功）+ 2 空
for c in range(12):
    filled = c < 10
    ax.add_patch(Circle((g2 + 0.34 + c * 0.225, 5.95), 0.075,
                        facecolor=(TEAL if filled else 'none'),
                        edgecolor=(TEAL if filled else GRAY), lw=0.9,
                        alpha=0.9 if filled else 0.8, zorder=3))
txt(g2 + 0.34, 5.62, '10 / 12 entries mapped', size=8.5, color=INK)
# 方向对比
txt(g2, 4.95, 'EXPECTED', size=8.5, color=GRAY, weight='semibold')
thin_arrow(g2 + 1.05, g2 + 1.85, 4.95, color=TEAL, lw=1.2)
txt(g2, 4.42, 'OBSERVED', size=8.5, color=GRAY, weight='semibold')
ax.annotate('', xy=(g2 + 1.05, 4.42), xytext=(g2 + 1.85, 4.42),
            arrowprops=dict(arrowstyle='-|>', color=CORAL, lw=1.2, mutation_scale=9))
txt(g2, 3.70, 'directions reversed', size=10, color=CORAL, weight='semibold')
txt(g2, 3.28, 'transfer to an independent', size=9, color=INK)
txt(g2, 2.95, 'RNA-seq cohort fails', size=9, color=INK)

# ================= GATE 3：mini-DCA 示意 + 数字事实 =================
g3 = X0[2]
bx, by, bw, bh = g3 + 0.42, 3.90, 2.30, 1.75
# 发丝坐标框
ax.add_patch(Rectangle((bx, by), bw, bh, facecolor='none', edgecolor=HAIR, lw=0.9, zorder=2))
xx = np.linspace(0, 1, 60)
def mapx(t): return bx + 0.12 + t * (bw - 0.24)
def mapy(v): return by + 0.14 + v * (bh - 0.30)
# treat-all（虚线，递减）、treat-none（虚线，水平）
ax.add_line(Line2D(mapx(xx), mapy(0.72 - 0.62 * xx), color=GRAY, lw=0.9, ls=(0, (3, 2)), zorder=3))
ax.add_line(Line2D([mapx(0), mapx(1)], [mapy(0.10), mapy(0.10)], color=GRAY, lw=0.9, ls=(0, (3, 2)), zorder=3))
# 模型曲线（金，中段高于两条默认策略）
nb = 0.14 + 0.58 * np.exp(-((xx - 0.42) ** 2) / 0.055) - 0.18 * xx
ax.add_line(Line2D(mapx(xx), mapy(np.clip(nb, 0.06, 0.92)), color=GOLD, lw=1.5, zorder=4))
txt(bx, by - 0.24, 'threshold', size=7.5, color=GRAY, ha='center')
txt(bx, by + bh + 0.16, 'net benefit', size=7.5, color=GRAY, ha='center')
txt(g3, 3.35, 'net benefit above default', size=9, color=INK)
txt(g3, 3.02, 'strategies across 10-45%', size=9, color=INK)
txt(g3, 2.69, 'risk thresholds', size=9, color=INK)
txt(g3, 2.30, 'adj AUC 0.690 (n=41) · 0.619 (n=766)', size=8.5, color=GOLD, weight='semibold')

# ================= GATE 4：节点图 + UCP2 开关 =================
g4 = X0[3]
cx, cy = g4 + 1.62, 4.90
angles = np.deg2rad([90, 150, 210, 270, 330, 30])
sats = [(cx + 0.78 * np.cos(a), cy + 0.78 * np.sin(a)) for a in angles]
for (sx, sy) in sats:
    ax.add_line(Line2D([cx, sx], [cy, sy], color=HAIR, lw=0.9, zorder=2))
for (sx, sy) in sats:
    ax.add_patch(Circle((sx, sy), 0.075, facecolor=HAIR, edgecolor='none', alpha=0.9, zorder=3))
ax.add_patch(Circle((cx, cy), 0.10, facecolor=GOLD, edgecolor='none', alpha=0.92, zorder=4))
# UCP2 开关：MCI 期 B 细胞（青绿下调）→ AD 期 NK（珊瑚上调）
ax.add_patch(Circle((sats[3][0], sats[3][1]), 0.085, facecolor='none', edgecolor=TEAL, lw=1.3, zorder=4))
ax.add_patch(Circle((sats[0][0], sats[0][1]), 0.085, facecolor='none', edgecolor=CORAL, lw=1.3, zorder=4))
txt(g4, 3.35, 'UCP2 stage-specific switch:', size=9, color=INK)
txt(g4, 3.02, 'B cells lower in MCI, NK cells', size=9, color=INK)
txt(g4, 2.69, 'higher in AD', size=9, color=INK)
txt(g4, 2.30, 'hypothesis-generating context', size=8.5, color=GRAY, style='italic')

# ================= 底部证据边界条 =================
ax.add_line(Line2D([0.50, 15.50], [1.55, 1.55], color=HAIR, lw=1.0, zorder=2))
txt(0.55, 1.12, 'EVIDENCE BOUNDARY   stable in development does not equal independent reproduction',
    size=9.5, color=INK)
txt(15.45, 1.12, 'not a validated biomarker', size=9.5, color=CORAL, weight='semibold', ha='right')

fig.savefig(os.path.join(OUT, 'GraphicalAbstract_v4_nature_grid.png'), dpi=300,
            facecolor=BG, bbox_inches=None)
fig.savefig(os.path.join(OUT, 'GraphicalAbstract_v4_nature_grid.tiff'), dpi=600,
            facecolor=BG, bbox_inches=None, pil_kwargs={'compression': 'tiff_lzw'})
print('GA saved:', OUT)
