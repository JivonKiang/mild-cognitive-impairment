#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
B-2 细胞来源定位：locked12 基因在 PBMC 各细胞类型的表达定位
数据：razer_figures_20261016/source_scRNA_locked12_expression_context.csv（12基因×9细胞×3组全网格）
产出：双面板热图（MCI vs HC delta z / AD vs HC delta z）Nature 网格
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = r'D:/20241004_MCI'
OUT_F = os.path.join(ROOT, '20261103_MCI_DCA_NRI/figures')
os.makedirs(OUT_F, exist_ok=True)

rows = list(csv.DictReader(open(os.path.join(ROOT, 'razer_figures_20261016/source_scRNA_locked12_expression_context.csv'))))

# 基因按 ledger 顺序（首现顺序）
genes = ['C22ORF39','EIF3B','GCA','GPI','LAT1-3TM','PHF15','PRDM10','RP2','UCP2','USP38','VAT1','ZSWIM6']
cts = ['Naive CD8 T cells','Central memory CD8 T cells','Effector memory CD8 T cells',
       'Terminal effector CD8 T cells','MAIT cells','Natural killer cells',
       'Non-switched memory B cells','Th2 cells','Neutrophil']
ct_short = {'Naive CD8 T cells':'Naive CD8','Central memory CD8 T cells':'CM CD8',
            'Effector memory CD8 T cells':'EM CD8','Terminal effector CD8 T cells':'Term eff CD8',
            'MAIT cells':'MAIT','Natural killer cells':'NK',
            'Non-switched memory B cells':'NSM B','Th2 cells':'Th2','Neutrophil':'Neutrophil'}

# mean_expr[pair][(gene,ct)]
me = {}
for r in rows:
    try:
        me.setdefault((r['disease_group'], r['panel_gene'], r['cell_type']), float(r['mean_expr']))
    except ValueError:
        pass  # NA 值跳过

def delta(g1, g0):
    M = np.zeros((len(genes), len(cts)))
    for i, g in enumerate(genes):
        for j, ct in enumerate(cts):
            M[i, j] = me.get((g1, g, ct), np.nan) - me.get((g0, g, ct), np.nan)
    return M

dM = delta('MCI', 'HC')   # MCI vs HC
dA = delta('AD', 'HC')    # AD vs HC

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.6})
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6))
vmax = np.nanmax(np.abs(np.concatenate([dM.ravel(), dA.ravel()])))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
for ax, M, pid, lab in zip(axes, [dM, dA], ['a', 'b'], ['MCI vs HC', 'AD vs HC']):
    im = ax.imshow(M.T, aspect='auto', cmap='RdBu_r', norm=norm)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=6)
    ax.set_yticks(range(len(cts)))
    ax.set_yticklabels([ct_short[c] for c in cts], fontsize=6.5)
    ax.text(0.0, 1.06, f'{pid}  {lab}', transform=ax.transAxes, fontsize=7, fontweight='bold')
    for i in range(len(genes)):
        for j in range(len(cts)):
            v = M[i, j]
            if not np.isnan(v):
                ax.text(i, j, f'{v:+.2f}', ha='center', va='center',
                        fontsize=4.6, color='white' if abs(v) > vmax*0.55 else '#222222')
cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
cb.set_label('Δ mean expression (log-norm)', fontsize=6.5)
cb.ax.tick_params(labelsize=6)
fig.savefig(os.path.join(OUT_F, 'Fig_celllocal_nature_grid.png'), dpi=600, bbox_inches='tight')
fig.savefig(os.path.join(OUT_F, 'Fig_celllocal_nature_grid.tiff'), dpi=600, bbox_inches='tight')
print('saved Fig_celllocal_nature_grid.png')
# 打印最强信号 cell×gene 供正文引用
order = np.dstack(np.unravel_index(np.argsort(-np.abs(dM), axis=None), dM.shape))[0]
print('Top MCI-vs-HC localizations:')
for i, j in order[:6]:
    print(f'  {genes[i]} @ {ct_short[cts[j]]}: {dM[i,j]:+.3f}')
