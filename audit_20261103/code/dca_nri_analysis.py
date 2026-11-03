#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
B-1 决策效用层：DCA 决策曲线 + NRI/IDI
数据：GSE249477 主发现队列（n=62）+ GSE140829 外部验证队列（n=1174）
目标 outcome：MCI vs HC（临床转译最关心的边界）
输出：analysis/dca_nri_results.json + figures/Fig_DCA_nature_grid.png
网格：Nature Figure Scale（统一面板宽、无图内 title/subtitle）
"""
import csv, json, os
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/20241004_MCI'
OUT_A = os.path.join(ROOT, '20261103_MCI_DCA_NRI/analysis')
OUT_F = os.path.join(ROOT, '20261103_MCI_DCA_NRI/figures')
os.makedirs(OUT_A, exist_ok=True); os.makedirs(OUT_F, exist_ok=True)

# ---------- 载入 ----------
def load_main():
    """GSE249477：稿件主定义 = CIA 冻结分（denominator_mode=full，n=41 MCI/HC）。
    注意：locked50_restricted 文件是 20261016 探索版分数（与冻结分 corr=-0.456），禁用于 DCA。"""
    f = os.path.join(ROOT, '20260903_Clinical_Interventions_in_Aging/03_KEY_DATA/analysis/primary/results/GSE249477_frozen_score_sample_scores.csv')
    rows = list(csv.DictReader(open(f)))
    rows = [r for r in rows if r.get('denominator_mode', 'full') == 'full']
    y, s, age, sex = [], [], [], []
    for r in rows:
        if r['stage'] in ('MCI', 'HC'):
            y.append(1 if r['stage'] == 'MCI' else 0)
            s.append(float(r['score']))
            age.append(float(r['age']) if r['age'] else np.nan)
            sex.append(1 if r['sex'].lower().startswith('m') else 0)
    return np.array(y), np.array(s), np.array(age), np.array(sex), 'GSE249477'

def load_valid(mode='full'):
    """GSE140829：CSV 每样本两行（denominator_mode=full/mapped_only）。
    稿件主定义 = 冻结权重原样应用（full，=正文 0.182 主结果）；mapped_only（component
    renormalization）为敏感性（正文 0.209）。mode='both' 返回两个分析集。"""
    f = os.path.join(ROOT, r'20260903_Clinical_Interventions_in_Aging/03_KEY_DATA/analysis/external_validation/GSE140829/results/GSE140829_frozen_score_sample_scores.csv')
    rows = list(csv.DictReader(open(f)))
    if mode != 'both':
        rows = [r for r in rows if r['denominator_mode'] == mode]
    y, s, age, sex = [], [], [], []
    for r in rows:
        if r['stage'] in ('MCI', 'HC'):
            y.append(1 if r['stage'] == 'MCI' else 0)
            s.append(float(r['score']))
            age.append(float(r['age']) if r['age'] not in ('', 'NA') else np.nan)
            sex.append(1 if r['sex'].lower().startswith('m') else 0)
    tag = 'GSE140829' if mode == 'both' else 'GSE140829_' + mode
    return np.array(y), np.array(s), np.array(age), np.array(sex), tag

# ---------- DCA 核心 ----------
def net_benefit(y, prob, thresh):
    n = len(y)
    tp = np.sum((prob >= thresh) & (y == 1))
    fp = np.sum((prob >= thresh) & (y == 0))
    return tp / n - fp / n * (thresh / (1 - thresh))

def dca_curve(y, score, name, adjust_covariates=False, age=None, sex=None,
              n_boot=500, seed=20261103):
    """score -> logistic prob。
    冻结取向闸：审计只承认冻结方向的决策效用。若样本内重拟合的 score 系数为负
    （模型自行翻转方向以拟合数据），则按冻结取向（prob -> 1-prob）重定向后报告，
    AUC 相应取 1-AUC；翻转本身记为 exploratory_flip 供正文披露。"""
    X = score.reshape(-1, 1)
    if adjust_covariates and age is not None:
        X = np.column_stack([score, np.nan_to_num(age), sex])
    lr = LogisticRegression(max_iter=2000)
    lr.fit(X, y)
    prob = lr.predict_proba(X)[:, 1]
    score_coef = float(lr.coef_[0][0])
    frozen_orientation = score_coef >= 0
    exploratory_flip = not frozen_orientation
    auc_raw = roc_auc_score(y, prob)
    if exploratory_flip:
        prob = 1.0 - prob
        auc = 1.0 - auc_raw
    else:
        auc = auc_raw
    threshs = np.arange(0.01, 0.99, 0.01)
    nb = np.array([net_benefit(y, prob, t) for t in threshs])
    # treat-all NB = prev - (1-prev)*t/(1-t); treat-none = 0
    prev = y.mean()
    nb_all = np.array([prev - (1 - prev) * t / (1 - t) for t in threshs])
    # bootstrap CI 带（重采样样本，重拟合 logistic，冻结取向同主拟合）
    rng = np.random.default_rng(seed)
    boot = np.zeros((n_boot, len(threshs)))
    n = len(y)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        yb, sb, ab, sxb = y[idx], score[idx], (age[idx] if age is not None else None), (sex[idx] if sex is not None else None)
        if len(np.unique(yb)) < 2:
            boot[b] = np.nan; continue
        Xb = sb.reshape(-1, 1)
        if adjust_covariates and ab is not None:
            Xb = np.column_stack([sb, np.nan_to_num(ab), sxb])
        lrb = LogisticRegression(max_iter=2000).fit(Xb, yb)
        pb = lrb.predict_proba(Xb)[:, 1]
        if float(lrb.coef_[0][0]) < 0:
            pb = 1.0 - pb
        boot[b] = [net_benefit(yb, pb, t) for t in threshs]
    nb_lo = np.nanpercentile(boot, 2.5, axis=0)
    nb_hi = np.nanpercentile(boot, 97.5, axis=0)
    # 决策效用窗口检验：0.10-0.45 阈值内净获益是否超两条默认策略
    win = (threshs >= 0.10) & (threshs <= 0.45)
    beats = bool(np.all(nb[win] > np.maximum(nb_all[win], 0.0)))
    return dict(name=name, auc=float(auc), auc_raw_if_flip=(float(auc_raw) if exploratory_flip else None),
                frozen_orientation=bool(frozen_orientation), exploratory_flip=bool(exploratory_flip),
                threshs=threshs.tolist(), nb_model=nb.tolist(),
                nb_lo=nb_lo.tolist(), nb_hi=nb_hi.tolist(),
                nb_all=nb_all.tolist(), prev=float(prev), n=len(y),
                beats_defaults_10_45=beats)

# ---------- NRI / IDI ----------
def nri_idi(y, score_new, score_old):
    """new score vs old (age/sex baseline model) 的 categorical NRI + IDI"""
    from sklearn.metrics import roc_auc_score
    Xo = np.column_stack([np.nan_to_num(score_old)]) if np.ndim(score_old) > 1 else score_old.reshape(-1, 1)
    Xm = np.column_stack([score_new.reshape(-1, 1), Xo])
    lr_o = LogisticRegression(max_iter=2000).fit(Xo, y)
    lr_m = LogisticRegression(max_iter=2000).fit(Xm, y)
    po = lr_o.predict_proba(Xo)[:, 1]
    pm = lr_m.predict_proba(Xm)[:, 1]
    # categorical NRI @ 0.05/0.2 (常用两档)
    def cat(p, lo=0.05, hi=0.2):
        return np.digitize(p, [lo, hi])
    co, cm = cat(po), cat(pm)
    e1 = np.mean((cm > co)[y == 1]); e2 = np.mean((cm < co)[y == 1])
    n1 = np.mean((cm < co)[y == 0]); n2 = np.mean((cm > co)[y == 0])
    nri = (e1 - e2) + (n1 - n2)
    # continuous NRI
    cnri = (np.mean((pm > po)[y == 1]) - np.mean((pm < po)[y == 1])) + \
           (np.mean((pm < po)[y == 0]) - np.mean((pm > po)[y == 0]))
    idi = po[y == 1].mean() - po[y == 0].mean()
    idi_new = pm[y == 1].mean() - pm[y == 0].mean()
    return dict(categorical_NRI=float(nri), continuous_NRI=float(cnri),
                IDI_old=float(idi), IDI_new=float(idi_new), IDI_delta=float(idi_new - idi))

# ---------- 执行 ----------
# 主定义校验：稿件 Results 4 = adjusted MCI-vs-HC 0.182 [0.064-0.301] P=0.0026（component renormalization）
# 先按两种分母模式各拟合一次线性对比（score ~ stage + age + sex），反推主模式，再跑 DCA/NRI
import json as _json

def verify_primary_mode():
    from sklearn.linear_model import LinearRegression
    out = {}
    for mode in ('mapped_only', 'full'):
        y, s, age, sex, _ = load_valid(mode)
        X = np.column_stack([np.nan_to_num(age), sex])
        # 设计阵：截距 + age + sex + MCI指示 + score
        Xd = np.column_stack([np.ones_like(y, dtype=float), X, y, s])
        lr = LinearRegression().fit(Xd, s)  # 占位，防止误用；真正对比用 statsmodels 式手工 OLS
        # 手工 OLS：outcome=score, 设计=1+age+sex+MCI
        Xo = np.column_stack([np.ones_like(y, dtype=float), np.nan_to_num(age), sex, y])
        beta, res, rank, sv = np.linalg.lstsq(Xo, s, rcond=None)
        b_mci = beta[3]
        resid = s - Xo @ beta
        sigma2 = resid @ resid / (len(s) - Xo.shape[1])
        XtX_inv = np.linalg.inv(Xo.T @ Xo)
        se = float(np.sqrt(sigma2 * XtX_inv[3, 3]))
        z = b_mci / se
        from scipy.stats import norm
        p = 2 * (1 - norm.cdf(abs(z)))
        out[mode] = dict(beta=round(float(b_mci), 4), se=round(se, 4),
                         ci=[round(b_mci - 1.96 * se, 4), round(b_mci + 1.96 * se, 4)], p=float(p))
    return out

checks = verify_primary_mode()
print('PRIMARY-MODE CHECK (target: 0.182 [0.064-0.301] P=0.0026):')
for m, v in checks.items():
    print(' ', m, v)
_json.dump(checks, open(os.path.join(OUT_A, 'denominator_mode_check.json'), 'w'), indent=1)

PRIMARY_MODE = 'full'  # 稿件主定义：冻结权重原样（正文 0.182 主结果与之精确吻合）
results = {}
for loader, adj in [(load_main, False), (load_main, True),
                    (lambda: load_valid(PRIMARY_MODE), False),
                    (lambda: load_valid(PRIMARY_MODE), True)]:
    y, s, age, sex, name = loader()
    key = name + ('_adj' if adj else '_crude')
    results[key] = dca_curve(y, s, key, adj, age, sex)
    # NRI vs age/sex-only baseline
    X_base = np.column_stack([np.nan_to_num(age), sex])
    results[key]['nri'] = nri_idi(y, s, X_base)
    d = results[key]
    flip_note = ' [EXPLORATORY FLIP -> frozen orientation]' if d['exploratory_flip'] else ''
    print(f"{key}: AUC={d['auc']:.3f} n={d['n']} prev={d['prev']:.3f}{flip_note} "
          f"cNRI={d['nri']['continuous_NRI']:.3f} beats10-45={d['beats_defaults_10_45']}")

# full 分母模式敏感性（GSE140829 only, adjusted）＝ component renormalization
y, s, age, sex, name = load_valid('mapped_only')
key = 'GSE140829_mapped_only_adj'
results[key] = dca_curve(y, s, key, True, age, sex)
X_base = np.column_stack([np.nan_to_num(age), sex])
results[key]['nri'] = nri_idi(y, s, X_base)
print(f"{key} (sensitivity): AUC={results[key]['auc']:.3f} n={results[key]['n']} "
      f"cNRI={results[key]['nri']['continuous_NRI']:.3f}")
json.dump(results, open(os.path.join(OUT_A, 'dca_nri_results.json'), 'w'), indent=1)

# ---------- Nature 网格图：双面板 DCA ----------
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
panels = [('GSE249477_adj', 'a'), ('GSE140829_full_adj', 'b')]
for ax, (key, pid) in zip(axes, panels):
    d = results[key]
    t = np.array(d['threshs'])
    orient = 'frozen' if d['frozen_orientation'] else 'reoriented to frozen'
    ax.plot(t, d['nb_model'], '-', lw=1.2, color='#B2182B', label=f'12-gene score (adj; {orient})')
    ax.fill_between(t, d['nb_lo'], d['nb_hi'], color='#B2182B', alpha=0.12, lw=0)
    ax.plot(t, d['nb_all'], '--', lw=0.9, color='#4D4D4D', label='treat all')
    ax.plot(t, np.zeros_like(t), ':', lw=0.9, color='#4D4D4D', label='treat none')
    ax.axhline(0, color='#4D4D4D', lw=0.5)
    ax.set_xlabel('Threshold probability')
    ax.set_ylabel('Net benefit')
    ax.set_xlim(0, 1)
    nb_all_arr = np.array(d['nb_all'])
    hi = max(np.nanmax(np.array(d['nb_hi'])[t < 0.5]), 0.05)
    ax.set_ylim(-0.1, hi * 1.15)
    ax.text(0.02, 0.96, f'{pid}  {d["name"]}', transform=ax.transAxes, va='top', fontsize=7, fontweight='bold')
    auc_note = f'AUC {d["auc"]:.2f} (n={d["n"]})'
    if d['exploratory_flip']:
        auc_note += f'; in-model {d["auc_raw_if_flip"]:.2f} reoriented'
    ax.text(0.98, 0.96, auc_note, transform=ax.transAxes, va='top', ha='right', fontsize=6.5)
    ax.legend(frameon=False, fontsize=6, loc='lower left')
fig.tight_layout()
fig.savefig(os.path.join(OUT_F, 'Fig_DCA_nature_grid.png'), dpi=600, bbox_inches='tight')
fig.savefig(os.path.join(OUT_F, 'Fig_DCA_nature_grid.tiff'), dpi=600, bbox_inches='tight')
print('FIG saved:', os.path.join(OUT_F, 'Fig_DCA_nature_grid.png'))
