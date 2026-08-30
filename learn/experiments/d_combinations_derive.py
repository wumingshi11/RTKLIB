#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块 D：载波相位与组合 —— 用真实数据一步步推导（D1~D3）
================================================
沿用固定数据（G03），从 RINEX 的 L1/L2/C1/P2 出发：
  D1 载波相位与整周模糊度
  D2 线性组合：GF（几何无关）、MW（Melbourne-Wübbena）
  D3 周跳：人为给 L1 加周跳，看 GF/MW 响应

运行：python3 learn/experiments/d_combinations_derive.py
"""
import re
import math
from pathlib import Path

C = 299792458.0
F1 = 1.57542e9
F2 = 1.22760e9
LAM1 = C / F1
LAM2 = C / F2
LAMW = C / (F1 - F2)   # 宽巷波长

ROOT = Path(__file__).resolve().parents[2]
OBS = ROOT / 'test/data/rinex/07590920.05o'


def read_g03_epochs(path, n=10):
    """读取 G03 前 n 个历元的 L1/C1/L2/P2。"""
    lines = path.read_text().splitlines()
    i = 0
    obs_types = None
    while 'END OF HEADER' not in lines[i]:
        if 'TYPES OF OBSERV' in lines[i]:
            parts = lines[i].split()
            obs_types = parts[1:1 + int(parts[0])]
        i += 1
    i += 1
    pat = re.compile(r'^\s*(\d{2})\s+\d+\s+\d+\s+\d+\s+\d+\s+[\d.]+\s+\d+\s+\d+[A-Z]')
    out = []
    while i < len(lines) and len(out) < n:
        l = lines[i]
        m = pat.match(l)
        if m:
            nsat = int(l[28:31])
            st = l[33:].rstrip()
            sats = []
            for k in range(0, len(st), 3):
                c = st[k:k + 3].strip()
                if c:
                    sats.append(int(c[:-1]) if c[-1].isalpha() else int(c))
            i += 1
            data = {}
            for prn in sats:
                toks = lines[i].split()
                i += 1
                data[prn] = {t: float(v.replace('D', 'E'))
                             for t, v in zip(obs_types, toks)}
            if 3 in data:
                out.append(data[3])
        else:
            i += 1
    return out


def gf_meters(l1_cyc, l2_cyc):
    """GF = λ1·L1 - λ2·L2（几何无关，米）。"""
    return l1_cyc * LAM1 - l2_cyc * LAM2


def mw_cycles(l1_cyc, l2_cyc, c1, p2):
    """MW（Melbourne-Wübbena），宽巷周数。
    公式：MW = [ (f1·L1m - f2·L2m)/(f1-f2) - (f1·P1+f2·P2)/(f1+f2) ] / λw
    其中 L1m=λ1·L1, L2m=λ2·L2。"""
    L1m = l1_cyc * LAM1
    L2m = l2_cyc * LAM2
    mw_m = (F1 * L1m - F2 * L2m) / (F1 - F2) - (F1 * c1 + F2 * p2) / (F1 + F2)
    return mw_m / LAMW


def main():
    print('=' * 78)
    print('模块 D：载波相位与组合 —— 真实数据推导（D1~D3）')
    print('=' * 78)
    eps = read_g03_epochs(OBS, n=10)
    d0 = eps[0]
    print(f'\n[输入] G03 前 {len(eps)} 个历元（30 s 采样）')
    print(f'       波长：λ1={LAM1:.6f} m, λ2={LAM2:.6f} m, λw={LAMW:.6f} m')

    # ---------- D1 ----------
    print('\n' + '-' * 78)
    print('D1 载波相位与整周模糊度（真实数字）')
    print('-' * 78)
    L1m = d0['L1'] * LAM1
    L2m = d0['L2'] * LAM2
    print(f'  第 1 个历元 G03：')
    print(f'    L1 = {d0["L1"]:.3f} cycle → λ1·L1 = {L1m:.3f} m')
    print(f'    L2 = {d0["L2"]:.3f} cycle → λ2·L2 = {L2m:.3f} m')
    print(f'    C1 = {d0["C1"]:.3f} m（伪距，无模糊）')
    print(f'    P2 = {d0["P2"]:.3f} m')
    # 用模块 B/C 的几何距离作参照
    rho_ref = 24873907.226  # 来自 B/C 固定数据（G03 含 Sagnac）
    print(f'    参考几何距离 ρ ≈ {rho_ref:.3f} m（来自 B/C 固定数据）')
    print(f'    λ1·L1 - ρ = {L1m - rho_ref:.3f} m ← 这就是未知整周模糊度 λ1·N1 的量级')
    print(f'    λ2·L2 - ρ = {L2m - rho_ref:.3f} m ← 同理')
    print('  → 载波相位观测值含有未知整数 N，不能直接当几何距离用；')
    print('    但它的噪声是毫米级，比伪距（分米~米级）精密得多。')

    # ---------- D2 ----------
    print('\n' + '-' * 78)
    print('D2 线性组合：GF 与 MW（真实序列）')
    print('-' * 78)
    print(f'  {"历元":<6}{"L1(cyc)":>14}{"L2(cyc)":>14}{"GF(m)":>12}{"MW(cyc)":>16}')
    gfs = []
    mws = []
    for k, d in enumerate(eps):
        g = gf_meters(d['L1'], d['L2'])
        m = mw_cycles(d['L1'], d['L2'], d['C1'], d['P2'])
        gfs.append(g)
        mws.append(m)
        print(f'  {k:<6}{d["L1"]:>14.3f}{d["L2"]:>14.3f}{g:>12.3f}{m:>16.3f}')
    print(f'\n  GF 变化范围：{min(gfs):.3f} ~ {max(gfs):.3f} m（约 {max(gfs)-min(gfs):.3f} m）')
    print(f'  MW 变化范围：{min(mws):.3f} ~ {max(mws):.3f} cyc（约 {max(mws)-min(mws):.3f} cyc）')
    print('  → GF 主要剩电离层+模糊度，随时间缓慢变化；')
    print('  → MW 几乎常数（宽巷模糊度+码噪声），是周跳/模糊度的好工具。')

    # ---------- D3 ----------
    print('\n' + '-' * 78)
    print('D3 周跳：人为给 L1 加周跳，看 GF/MW 响应')
    print('-' * 78)
    print('  以第 3 个历元为基准，给 L1 加 +1 周、+5 周：')
    base = eps[2]
    g0 = gf_meters(base['L1'], base['L2'])
    m0 = mw_cycles(base['L1'], base['L2'], base['C1'], base['P2'])
    print(f'  原始：GF={g0:.3f} m, MW={m0:.3f} cyc')
    for slip in [1, 5]:
        l1s = base['L1'] + slip
        gs = gf_meters(l1s, base['L2'])
        ms = mw_cycles(l1s, base['L2'], base['C1'], base['P2'])
        print(f'  L1+{slip} 周：GF={gs:.3f} m (Δ={gs-g0:+.3f}), '
              f'MW={ms:.3f} cyc (Δ={ms-m0:+.3f})')
    print('  → GF 对 L1 单频周跳敏感：+1 周 ≈ +0.190 m；')
    print('  → MW 对周跳也敏感（+1 周 L1 → MW 也跳约 1 周量级）。')
    print('  → 实际中常用“GF 历元间差分”或“MW 与上一历元比较”检测周跳。')


if __name__ == '__main__':
    main()
