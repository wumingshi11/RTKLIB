#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块 E：差分与 RTK 浮点解 —— 用真实数据一步步推导（E1~E5）
================================================
用 rover(0759) 和 base(3040) 第一个历元、共视卫星：
  E1 单差/双差构造
  E2 双差协方差 R（参考星导致互相关）
  E3 RTK 状态布局（浮点状态）
  E4/E5 用双差相位+码做一次浮点最小二乘（简化 Kalman 更新）

运行：python3 learn/experiments/e_dd_derive.py
"""
import re
import math
from pathlib import Path

C = 299792458.0
MU = 3.9860050E14
OMEGA_E = 7.2921151467E-5
F = -4.442807633E-10
LAM1 = C / 1.57542e9

ROOT = Path(__file__).resolve().parents[2]
NAV = ROOT / 'test/data/rinex/07590920.05n'
OBS_R = ROOT / 'test/data/rinex/07590920.05o'   # rover
OBS_B = ROOT / 'test/data/rinex/30400920.05o'   # base
TOW0 = 518400.0

BASE_POS = (-3978242.4348, 3382841.1715, 3649902.7667)
ROVER_APPROX = (-3976210.472, 3382362.102, 3652503.783)  # 来自 C 模块 SPP 解


# ---------------- RINEX 解析 ----------------
def parse_tail_fields(line):
    s = line[3:]
    out = []
    for j in range(0, len(s), 19):
        f = s[j:j + 19].strip()
        if f:
            out.append(float(f.replace('D', 'E')))
    return out


def read_nav(path):
    lines = path.read_text().splitlines()
    i = 0
    while 'END OF HEADER' not in lines[i]:
        i += 1
    i += 1
    nav = {}
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        sat = int(line[0:2].strip())
        block = [lines[i + k] for k in range(8) if i + k < len(lines)]
        i += 8
        if len(block) < 8:
            continue
        def fld(a, b):
            return block[0][a:b].strip()
        af0 = float(fld(22, 41).replace('D', 'E'))
        af1 = float(fld(41, 60).replace('D', 'E'))
        af2 = float(fld(60, 79).replace('D', 'E'))
        v = []
        for b in block[1:]:
            v += parse_tail_fields(b)
        eph = {'IODE': v[0], 'Crs': v[1], 'deln': v[2], 'M0': v[3],
               'Cuc': v[4], 'e': v[5], 'Cus': v[6], 'sqrtA': v[7],
               'toe': v[8], 'Cic': v[9], 'OMG0': v[10], 'Cis': v[11],
               'i0': v[12], 'Crc': v[13], 'omg': v[14], 'OMGd': v[15],
               'idot': v[16], 'codes': v[17], 'week': v[18], 'L2flag': v[19],
               'svacc': v[20], 'svh': v[21], 'tgd': v[22], 'iodc': v[23],
               'ttm': v[24], 'af0': af0, 'af1': af1, 'af2': af2}
        nav.setdefault(sat, []).append((eph['toe'], eph))
    return nav


def read_first_epoch(path):
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
    while i < len(lines):
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
            return sats, data
        i += 1
    return None, None


# ---------------- 卫星位置（B 模块一步法） ----------------
def kepler(M, e):
    E = M
    for _ in range(30):
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < 1e-14:
            break
    return E


def eph2pos(eph, t):
    A = eph['sqrtA'] ** 2
    n0 = math.sqrt(MU / A ** 3)
    tk = t - eph['toe']
    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800
    n = n0 + eph['deln']
    M = eph['M0'] + n * tk
    E = kepler(M, eph['e'])
    nu = math.atan2(math.sqrt(1 - eph['e'] ** 2) * math.sin(E),
                    math.cos(E) - eph['e'])
    phi = nu + eph['omg']
    du = eph['Cus'] * math.sin(2 * phi) + eph['Cuc'] * math.cos(2 * phi)
    dr = eph['Crs'] * math.sin(2 * phi) + eph['Crc'] * math.cos(2 * phi)
    di = eph['Cis'] * math.sin(2 * phi) + eph['Cic'] * math.cos(2 * phi)
    u = phi + du
    r = A * (1 - eph['e'] * math.cos(E)) + dr
    i = eph['i0'] + di + eph['idot'] * tk
    xp = r * math.cos(u)
    yp = r * math.sin(u)
    O = eph['OMG0'] + (eph['OMGd'] - OMEGA_E) * tk - OMEGA_E * eph['toe']
    x = xp * math.cos(O) - yp * math.cos(i) * math.sin(O)
    y = xp * math.sin(O) + yp * math.cos(i) * math.cos(O)
    z = yp * math.sin(i)
    dts = eph['af0'] + eph['af1'] * tk + eph['af2'] * tk * tk
    dts += F * eph['e'] * math.sqrt(A) * math.sin(E)
    return (x, y, z), dts


def satpos(eph, P, t_rx):
    t1 = t_rx - P / C
    _, dts1 = eph2pos(eph, t1)
    t_emit = t1 - dts1
    rs, dts = eph2pos(eph, t_emit)
    return rs, dts


def geodist(rs, rr):
    dx, dy, dz = rs[0] - rr[0], rs[1] - rr[1], rs[2] - rr[2]
    rho = math.sqrt(dx * dx + dy * dy + dz * dz)
    rho += OMEGA_E / C * (rs[0] * rr[1] - rs[1] * rr[0])
    return rho


def ecef2pos(r):
    x, y, z = r
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    e2 = (1 / 298.257223563) * (2 - 1 / 298.257223563)
    lat = math.atan2(z, p * (1 - e2))
    h = 0.0
    for _ in range(5):
        N = 6378137.0 / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    return lat, lon, h


def xyz2enu(r, lat, lon):
    sinp, cosp = math.sin(lat), math.cos(lat)
    sinl, cosl = math.sin(lon), math.cos(lon)
    return (-sinl * r[0] + cosl * r[1],
            -sinp * cosl * r[0] - sinp * sinl * r[1] + cosp * r[2],
            cosp * cosl * r[0] + cosp * sinl * r[1] + sinp * r[2])


def solve(A, b):
    """高斯消元 + 回代解线性方程组 A x = b。"""
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        pv = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[pv] = M[pv], M[col]
        piv = M[col][col]
        for c in range(col, n + 1):
            M[col][c] /= piv
        for r in range(n):
            if r == col:
                continue
            f = M[r][col]
            for c in range(col, n + 1):
                M[r][c] -= f * M[col][c]
    return [M[i][n] for i in range(n)]


def main():
    print('=' * 78)
    print('模块 E：差分与 RTK 浮点解 —— 真实数据推导（E1~E5）')
    print('=' * 78)
    nav = read_nav(NAV)
    sats_r, obs_r = read_first_epoch(OBS_R)
    sats_b, obs_b = read_first_epoch(OBS_B)
    common = sorted(set(sats_r) & set(sats_b))
    print(f'\n[输入] rover 可见星 = {sats_r}')
    print(f'       base  可见星 = {sats_b}')
    print(f'       共视星 = {common}')

    # 卫星位置（用 rover 伪距一步法，共视星同一位置）
    rs = {}
    for prn in common:
        eph = min(nav[prn], key=lambda x: abs(x[0] - TOW0))[1]
        r, dts = satpos(eph, obs_r[prn]['C1'], TOW0)
        rs[prn] = r

    # 高度角（用 rover 近似位置选参考星）
    pos_r = ecef2pos(ROVER_APPROX)
    el = {}
    for prn in common:
        e = ((rs[prn][0] - ROVER_APPROX[0]) / math.dist(rs[prn], ROVER_APPROX),
             (rs[prn][1] - ROVER_APPROX[1]) / math.dist(rs[prn], ROVER_APPROX),
             (rs[prn][2] - ROVER_APPROX[2]) / math.dist(rs[prn], ROVER_APPROX))
        enu = xyz2enu(e, pos_r[0], pos_r[1])
        el[prn] = math.asin(enu[2])
    ref = max(common, key=lambda p: el[p])
    print(f'       参考星（最高仰角）= G{ref:02d}, el={math.degrees(el[ref]):.2f}°')

    # ---------- E1 ----------
    print('\n' + '-' * 78)
    print('E1 单差 / 双差构造（真实数字，L1 相位转米）')
    print('-' * 78)
    print(f'  {"星":<6}{"SD码(m)":>14}{"SD相(m)":>14}{"DD码(m)":>14}{"DD相(m)":>14}{"DD几何(m)":>14}')
    dd_code = {}
    dd_phase = {}
    dd_geo = {}
    for prn in common:
        sd_code = obs_r[prn]['C1'] - obs_b[prn]['C1']
        sd_phase = LAM1 * (obs_r[prn]['L1'] - obs_b[prn]['L1'])
        if prn == ref:
            print(f'  G{prn:<4}{sd_code:>14.3f}{sd_phase:>14.3f}{"-":>14}{"-":>14}{"-":>14}  (参考星)')
            continue
        dc = sd_code - (obs_r[ref]['C1'] - obs_b[ref]['C1'])
        dp = sd_phase - LAM1 * (obs_r[ref]['L1'] - obs_b[ref]['L1'])
        rho_r = geodist(rs[prn], ROVER_APPROX)
        rho_b = geodist(rs[prn], BASE_POS)
        rho_rr = geodist(rs[ref], ROVER_APPROX)
        rho_br = geodist(rs[ref], BASE_POS)
        dg = (rho_r - rho_b) - (rho_rr - rho_br)
        dd_code[prn] = dc
        dd_phase[prn] = dp
        dd_geo[prn] = dg
        print(f'  G{prn:<4}{sd_code:>14.3f}{sd_phase:>14.3f}{dc:>14.3f}{dp:>14.3f}{dg:>14.3f}')

    # ---------- E2 ----------
    print('\n' + '-' * 78)
    print('E2 双差协方差 R（参考星导致互相关）')
    print('-' * 78)
    targets = [p for p in common if p != ref]
    n = len(targets)
    print(f'  非参考星 = {["G%02d" % p for p in targets]}（{n} 条双差）')
    print('  若单差噪声 σ²，则 R = σ²·D·Dᵀ，D 为 (n)×(n+1) 差分矩阵')
    # D matrix: rows = target - ref
    D = [[0.0] * (n + 1) for _ in range(n)]
    for i, p in enumerate(targets):
        D[i][0] = -1.0       # ref 列
        D[i][i + 1] = 1.0    # target 列
    # R = D D^T
    R = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            R[i][j] = sum(D[i][k] * D[j][k] for k in range(n + 1))
    print('  R/σ² 矩阵（对角=2，非对角=1）：')
    print('       ' + ''.join(f'{"G%02d" % p:>8}' for p in targets))
    for i, p in enumerate(targets):
        print(f'  G{p:02d}  ' + ''.join(f'{R[i][j]:>8.1f}' for j in range(n)))
    print('  → 同参考星导致双差之间互相关：对角=2σ²，非对角=σ²')

    # ---------- E3/E4/E5 ----------
    print('\n' + '-' * 78)
    print('E3/E4/E5 双差浮点最小二乘（简化：位置3 + 浮点模糊度n）')
    print('-' * 78)
    # 状态：[dx,dy,dz, N1..Nn]（N 为双差相位模糊度，周）
    m = n
    # 观测：每条双差 1 个码 + 1 个相位 = 2n 个
    H = []
    v = []
    W = []
    for i, p in enumerate(targets):
        # 方向导数：e_ref - e_s
        es = ((rs[p][0] - ROVER_APPROX[0]) / math.dist(rs[p], ROVER_APPROX),
              (rs[p][1] - ROVER_APPROX[1]) / math.dist(rs[p], ROVER_APPROX),
              (rs[p][2] - ROVER_APPROX[2]) / math.dist(rs[p], ROVER_APPROX))
        er = ((rs[ref][0] - ROVER_APPROX[0]) / math.dist(rs[ref], ROVER_APPROX),
              (rs[ref][1] - ROVER_APPROX[1]) / math.dist(rs[ref], ROVER_APPROX),
              (rs[ref][2] - ROVER_APPROX[2]) / math.dist(rs[ref], ROVER_APPROX))
        g = (er[0] - es[0], er[1] - es[1], er[2] - es[2])
        # 码观测
        row = [0.0] * (3 + m)
        row[0], row[1], row[2] = g
        H.append(row)
        v.append(dd_code[p] - dd_geo[p])
        w_code = 1.0 / max(math.sin(el[p]) ** 2, 1e-6)
        W.append([w_code])
        # 相位观测（含模糊度）
        row = [0.0] * (3 + m)
        row[0], row[1], row[2] = g
        row[3 + i] = LAM1
        H.append(row)
        v.append(dd_phase[p] - dd_geo[p])
        w_phase = 1e4 * w_code   # 相位比码精很多
        W.append([w_phase])
    Wm = [[0.0] * len(v) for _ in range(len(v))]
    for i in range(len(v)):
        Wm[i][i] = W[i][0]
    # 法方程
    A = [[0.0] * (3 + m) for _ in range(3 + m)]
    b = [0.0] * (3 + m)
    for i in range(len(v)):
        wi = Wm[i][i]
        for r in range(3 + m):
            b[r] += H[i][r] * wi * v[i]
            for c in range(3 + m):
                A[r][c] += H[i][r] * wi * H[i][c]
    dx = solve(A, b)
    print(f'  状态数 = 3 位置 + {m} 浮点模糊度')
    print(f'  基线修正 = ({dx[0]:.3f}, {dx[1]:.3f}, {dx[2]:.3f}) m')
    rr_new = [ROVER_APPROX[i] + dx[i] for i in range(3)]
    lat, lon, h = ecef2pos(rr_new)
    print(f'  新 rover ECEF = ({rr_new[0]:.3f}, {rr_new[1]:.3f}, {rr_new[2]:.3f}) m')
    print(f'  新 rover LLH  = {math.degrees(lat):.6f}, {math.degrees(lon):.6f}, {h:.3f} m')
    base_llh = ecef2pos(BASE_POS)
    print(f'  base     LLH  = {math.degrees(base_llh[0]):.6f}, {math.degrees(base_llh[1]):.6f}, {base_llh[2]:.3f} m')
    print(f'  浮点模糊度 N（周）：')
    for i, p in enumerate(targets):
        print(f'    G{p:02d}-G{ref:02d}: {dx[3+i]:.3f} 周')
    print(f'  （这是浮点解；F 模块用 LAMBDA 把这些 N 固定成整数）')


if __name__ == '__main__':
    main()
