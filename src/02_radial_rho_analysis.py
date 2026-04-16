from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import spearmanr

print("RUNNING BPT-CUT SCRIPT")

# =========================
# 0. 路径
# =========================
BASE = Path(r"C:\Users\30126\Documents\MaNGA_reproduction")
FN = BASE / "data" / "manga-7443-12703-MAPS-HYB10-MILESHC-MASTARHC2.fits"
OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

# =========================
# 1. 读取 MAPS 文件
# =========================
with fits.open(FN) as hdul:
    rre = hdul["SPX_ELLCOO"].data[1].astype(float)   # channel 2 => Python index 1
    gflux = hdul["EMLINE_GFLUX"].data.astype(float)
    givar = hdul["EMLINE_GFLUX_IVAR"].data.astype(float)
    gmask = hdul["EMLINE_GFLUX_MASK"].data
    gsigma = hdul["EMLINE_GSIGMA"].data.astype(float)

# =========================
# 2. 通道编号
# =========================
CH_OII_3727 = 0
CH_OII_3729 = 1
CH_HB       = 14
CH_OIII_5008 = 16
CH_HA       = 23
CH_NII_6585 = 24
CH_SII_6718 = 25
CH_SII_6732 = 26

# =========================
# 3. 提取二维图
# =========================
oii = gflux[CH_OII_3727] + gflux[CH_OII_3729]
hb = gflux[CH_HB]
oiii = gflux[CH_OIII_5008]
ha = gflux[CH_HA]
nii = gflux[CH_NII_6585]
sii = gflux[CH_SII_6718] + gflux[CH_SII_6732]
sigma_ha = gsigma[CH_HA]

# =========================
# 4. 由 IVAR 估计误差与 S/N
# =========================
def variance_from_ivar(ivar):
    out = np.full_like(ivar, np.nan, dtype=float)
    m = np.isfinite(ivar) & (ivar > 0)
    out[m] = 1.0 / ivar[m]
    return out

var_oii = variance_from_ivar(givar[CH_OII_3727]) + variance_from_ivar(givar[CH_OII_3729])
var_sii = variance_from_ivar(givar[CH_SII_6718]) + variance_from_ivar(givar[CH_SII_6732])

sn_oii = np.full_like(oii, np.nan, dtype=float)
m_oii = np.isfinite(var_oii) & (var_oii > 0)
sn_oii[m_oii] = oii[m_oii] / np.sqrt(var_oii[m_oii])

sn_sii = np.full_like(sii, np.nan, dtype=float)
m_sii = np.isfinite(var_sii) & (var_sii > 0)
sn_sii[m_sii] = sii[m_sii] / np.sqrt(var_sii[m_sii])

sn_hb = np.full_like(hb, np.nan, dtype=float)
m_hb = np.isfinite(givar[CH_HB]) & (givar[CH_HB] > 0)
sn_hb[m_hb] = hb[m_hb] * np.sqrt(givar[CH_HB][m_hb])

sn_ha = np.full_like(ha, np.nan, dtype=float)
m_ha = np.isfinite(givar[CH_HA]) & (givar[CH_HA] > 0)
sn_ha[m_ha] = ha[m_ha] * np.sqrt(givar[CH_HA][m_ha])

sn_oiii = np.full_like(oiii, np.nan, dtype=float)
m_oiii = np.isfinite(givar[CH_OIII_5008]) & (givar[CH_OIII_5008] > 0)
sn_oiii[m_oiii] = oiii[m_oiii] * np.sqrt(givar[CH_OIII_5008][m_oiii])

sn_nii = np.full_like(nii, np.nan, dtype=float)
m_nii = np.isfinite(givar[CH_NII_6585]) & (givar[CH_NII_6585] > 0)
sn_nii[m_nii] = nii[m_nii] * np.sqrt(givar[CH_NII_6585][m_nii])

# =========================
# 5. 简化版 mask 过滤
# =========================
good_mask = (
    (gmask[CH_HB] == 0) &
    (gmask[CH_HA] == 0) &
    (gmask[CH_OIII_5008] == 0) &
    (gmask[CH_NII_6585] == 0) &
    (gmask[CH_OII_3727] == 0) &
    (gmask[CH_OII_3729] == 0) &
    (gmask[CH_SII_6718] == 0) &
    (gmask[CH_SII_6732] == 0)
)

# =========================
# 6. S/N cuts
# =========================
good_lines = (
    np.isfinite(rre) &
    good_mask &
    np.isfinite(ha) & np.isfinite(hb) &
    np.isfinite(oiii) & np.isfinite(oii) &
    np.isfinite(nii) & np.isfinite(sii) &
    np.isfinite(sigma_ha) &
    (ha > 0) & (hb > 0) &
    (oiii > 0) & (oii > 0) &
    (nii > 0) & (sii > 0) &
    (sigma_ha > 0) &
    np.isfinite(sn_ha) & (sn_ha > 5) &
    np.isfinite(sn_hb) & (sn_hb > 5) &
    np.isfinite(sn_oii) & (sn_oii > 3) &
    np.isfinite(sn_oiii) & (sn_oiii > 3) &
    np.isfinite(sn_nii) & (sn_nii > 3) &
    np.isfinite(sn_sii) & (sn_sii > 3)
)

print("Number of valid spaxels after S/N cut =", np.sum(good_lines))

# =========================
# 7. 构造 BPT 所需线比
# =========================
log_o3hb = np.full_like(ha, np.nan, dtype=float)
m = good_lines & (oiii > 0) & (hb > 0)
log_o3hb[m] = np.log10(oiii[m] / hb[m])

log_n2ha = np.full_like(ha, np.nan, dtype=float)
m = good_lines & (nii > 0) & (ha > 0)
log_n2ha[m] = np.log10(nii[m] / ha[m])

log_s2ha = np.full_like(ha, np.nan, dtype=float)
m = good_lines & (sii > 0) & (ha > 0)
log_s2ha[m] = np.log10(sii[m] / ha[m])

# =========================
# 8. 简化版 BPT 分类
# =========================
# [NII]-BPT: 用经典 AGN 分界（Kewley 2001 style）
nii_bpt_agn = (
    np.isfinite(log_o3hb) &
    np.isfinite(log_n2ha) &
    (log_n2ha < 0.47) &
    (log_o3hb > (0.61 / (log_n2ha - 0.47) + 1.19))
)

# [SII]-BPT: 先分出 AGN，再用 Seyfert/LINER 分界
sii_bpt_agn = (
    np.isfinite(log_o3hb) &
    np.isfinite(log_s2ha) &
    (log_s2ha < 0.32) &
    (log_o3hb > (0.72 / (log_s2ha - 0.32) + 1.30))
)

sii_bpt_seyfert = (
    sii_bpt_agn &
    (log_o3hb > (1.89 * log_s2ha + 0.76))
)

sii_bpt_liner = (
    sii_bpt_agn &
    (~sii_bpt_seyfert)
)

# =========================
# 9. 近似中心 3 arcsec 区域
#    MaNGA: 0.5"/spaxel => 3" ~ 6 spaxels
# =========================
iy0, ix0 = np.unravel_index(np.nanargmin(rre), rre.shape)

yy, xx = np.indices(rre.shape)
rpix = np.sqrt((xx - ix0)**2 + (yy - iy0)**2)
r_arcsec = 0.5 * rpix

central_3arcsec = r_arcsec <= 3.0

# =========================
# 10. 论文思路的简化版 BPT 剔除
#     - 全部剔除 [SII]-BPT Seyfert
#     - 中心 3" 内，额外剔除 [SII]-BPT LINER
#     - 中心 3" 内，额外剔除 [NII]-BPT AGN
# =========================
exclude_bpt = (
    sii_bpt_seyfert |
    (central_3arcsec & sii_bpt_liner) |
    (central_3arcsec & nii_bpt_agn)
)

print("Number of Seyfert-like spaxels removed =", np.sum(sii_bpt_seyfert))
print("Number of central LINER-like spaxels removed =", np.sum(central_3arcsec & sii_bpt_liner))
print("Number of central NII-BPT AGN-like spaxels removed =", np.sum(central_3arcsec & nii_bpt_agn))

good_final = good_lines & (~exclude_bpt)
print("Number of valid spaxels after simplified BPT cut =", np.sum(good_final))

# =========================
# 11. 构造最终物理量
# =========================
k_hb = 3.61
k_ha = 2.53

ebv_gas = np.full_like(ha, np.nan, dtype=float)
ha_hb = np.full_like(ha, np.nan, dtype=float)
ha_hb[good_final] = ha[good_final] / hb[good_final]

valid = good_final & np.isfinite(ha_hb) & (ha_hb > 0)
ebv_gas[valid] = 2.5 / (k_hb - k_ha) * np.log10(ha_hb[valid] / 2.86)
ebv_gas[ebv_gas < 0] = np.nan

valid = np.isfinite(ebv_gas)

log_o3o2 = np.full_like(ha, np.nan, dtype=float)
m = valid & np.isfinite(oiii) & np.isfinite(oii) & (oiii > 0) & (oii > 0)
log_o3o2[m] = np.log10(oiii[m] / oii[m])

log_n2s2 = np.full_like(ha, np.nan, dtype=float)
m = valid & np.isfinite(nii) & np.isfinite(sii) & (nii > 0) & (sii > 0)
log_n2s2[m] = np.log10(nii[m] / sii[m])

log_sigma_ha = np.full_like(ha, np.nan, dtype=float)
m = valid & np.isfinite(sigma_ha) & (sigma_ha > 0)
log_sigma_ha[m] = np.log10(sigma_ha[m])

# =========================
# 12. 半径 bins
# =========================
rbins = [(0.0, 0.3), (0.3, 0.6), (0.6, 0.9), (0.9, 1.2), (1.2, 1.5), (1.5, 2.5)]

def collect_binned_data(xmap, ymap, rmap, rbins):
    results = []
    for r0, r1 in rbins:
        m = (
            np.isfinite(xmap) &
            np.isfinite(ymap) &
            np.isfinite(rmap) &
            (rmap >= r0) & (rmap < r1)
        )

        x = xmap[m]
        y = ymap[m]

        if len(x) > 5:
            rho, p = spearmanr(x, y)
        else:
            rho, p = np.nan, np.nan

        results.append({
            "rbin": (r0, r1),
            "x": x,
            "y": y,
            "rho": rho,
            "p": p,
            "n": len(x),
        })
    return results

res_o3o2 = collect_binned_data(log_o3o2, ebv_gas, rre, rbins)
res_n2s2 = collect_binned_data(log_n2s2, ebv_gas, rre, rbins)
res_sig = collect_binned_data(log_sigma_ha, ebv_gas, rre, rbins)

# =========================
# 13. 打印 rho
# =========================
print("\n=== Spearman rho: E(B-V)_gas vs log([OIII]/[OII]) ===")
for r in res_o3o2:
    print(f"R/Re {r['rbin']}: rho={r['rho']:.3f}, N={r['n']}")

print("\n=== Spearman rho: E(B-V)_gas vs log([NII]/[SII]) ===")
for r in res_n2s2:
    print(f"R/Re {r['rbin']}: rho={r['rho']:.3f}, N={r['n']}")

print("\n=== Spearman rho: E(B-V)_gas vs log(sigma_Ha) ===")
for r in res_sig:
    print(f"R/Re {r['rbin']}: rho={r['rho']:.3f}, N={r['n']}")

# =========================
# 14. 作图函数
# =========================
rng = np.random.default_rng(42)

def subsample_xy(x, y, max_points=4000):
    n = len(x)
    if n <= max_points:
        return x, y
    idx = rng.choice(n, size=max_points, replace=False)
    return x[idx], y[idx]

def plot_panel_set(results, xlabel, outname, ylabel="E(B-V)_gas"):
    fig, axes = plt.subplots(1, len(results), figsize=(20, 3.8), sharey=True)

    for ax, r in zip(axes, results):
        x, y = subsample_xy(r["x"], r["y"], max_points=4000)
        ax.scatter(x, y, s=3, alpha=0.15)

        if len(x) > 30:
            nbins = 12
            edges = np.linspace(np.nanmin(x), np.nanmax(x), nbins + 1)
            xc, ym = [], []
            for i in range(nbins):
                mm = (x >= edges[i]) & (x < edges[i + 1])
                if np.sum(mm) > 5:
                    xc.append(np.nanmedian(x[mm]))
                    ym.append(np.nanmedian(y[mm]))
            if len(xc) > 1:
                ax.plot(xc, ym, color="black", lw=2)

        r0, r1 = r["rbin"]
        ax.set_title(f"R/Re: [{r0}, {r1})\nrho={r['rho']:.3f}")
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.2)

    axes[0].set_ylabel(ylabel)
    plt.tight_layout()
    savepath = OUT / outname
    plt.savefig(savepath, dpi=200, bbox_inches="tight")
    print(f"Saved: {savepath}")
    plt.show()

# =========================
# 15. 画三组图
# =========================
plot_panel_set(res_o3o2, "log([OIII]/[OII])", "rho_o3o2_bpt.png")
plot_panel_set(res_n2s2, "log([NII]/[SII])", "rho_n2s2_bpt.png")
plot_panel_set(res_sig, "log(sigma_Ha)", "rho_sigma_ha_bpt.png")