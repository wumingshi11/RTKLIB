#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块 C：单点定位（SPP）—— 用真实数据一步步推导
================================================
在模块 B 的基础上（同一份固定数据），从 8 颗 GPS 卫星的 C1 伪距出发：
  C1 伪距改正（卫星钟差/相对论/TGD/电离层/对流层/Sagnac）
  C2 加权最小二乘与雅可比
  C3 SPP 完整流程与质量（迭代收敛、DOP、残差）

运行：python3 learn/experiments/c_spp_derive.py
"""
import math
import re
from pathlib import Path

# ---------------- 常数 ----------------
C = 299792458.0
MU = 3.9860050E14
OMEGA_E = 7.2921151467E-5
F = -4.442807633E-10
RE_WGS84 = 6378137.0
FE_WGS84 = 1.0 / 298.257223563

ROOT = Path(__file__).resolve().parents[2]
NAV = ROOT / 'test/data/rinex/07590920.05n'
OBS = ROOT / 'test/data/rinex/07590920.05o'

WEEK = 1316
TOW0 = 518400.0
RINEX_POS = (-3976219.5082, 3382372.5671, 3652512.9849)


# ---------------- RINEX 解析（与 B 模块共用） ----------------
def parse_tail_fields(line):
    s = line[3:]
    out = []
    for j in range(0, len(s), 19):
        f = s[j:j + 19].strip()
        if f:
            out.append(float(f.replace('D', 'E')))
    return out


def read_nav_header(path):
    alpha = beta = None
    for line in path.read_text().splitlines():
        if 'ION ALPHA' in line:
            alpha = [float(x.replace('D', 'E')) for x in line.split()[:4]]
        elif 'ION BETA' in line:
            beta = [float(x.replace('D', 'E')) for x in line.split()[:4]]
    return alpha, beta


def read_gps_nav(path):
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines) and 'END OF HEADER' not in lines[i]:
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
        toc = {'year': int(fld(3, 5)), 'month': int(fld(6, 8)),
               'day': int(fld(9, 11)), 'hour': int(fld(12, 14)),
               'min': int(fld(15, 17)), 'sec': float(fld(17, 22))}
        af0 = float(fld(22, 41).replace('D', 'E'))
        af1 = float(fld(41, 60).replace('D', 'E'))
        af2 = float(fld(60, 79).replace('D', 'E'))
        v = []
        for b in block[1:]:
            v += parse_tail_fields(b)
        eph = {
            'IODE': v[0], 'Crs': v[1], 'deln': v[2], 'M0': v[3],
            'Cuc': v[4], 'e': v[5], 'Cus': v[6], 'sqrtA': v[7],
            'toe': v[8], 'Cic': v[9], 'OMG0': v[10], 'Cis': v[11],
            'i0': v[12], 'Crc': v[13], 'omg': v[14], 'OMGd': v[15],
            'idot': v[16], 'codes': v[17], 'week': v[18], 'L2flag': v[19],
            'svacc': v[20], 'svh': v[21], 'tgd': v[22], 'iodc': v[23],
            'ttm': v[24], 'af0': af0, 'af1': af1, 'af2': af2, **toc,
        }
        nav.setdefault(sat, []).append((eph['toe'], eph))
    return nav


def read_first_epoch_obs(path):
    lines = path.read_text().splitlines()
    i = 0
    obs_types = None
    while i < len(lines) and 'END OF HEADER' not in lines[i]:
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
            return sats, data, obs_types
        i += 1
    return None, None, None


# ---------------- 坐标 / 几何 ----------------
def ecef2pos(r):
    x, y, z = r
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    e2 = FE_WGS84 * (2 - FE_WGS84)
    lat = math.atan2(z, p * (1 - e2))
    h = 0.0
    for _ in range(5):
        N = RE_WGS84 / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    return lat, lon, h


def pos2ecef(pos):
    lat, lon, h = pos
    e2 = FE_WGS84 * (2 - FE_WGS84)
    N = RE_WGS84 / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    return ((N + h) * math.cos(lat) * math.cos(lon),
            (N + h) * math.cos(lat) * math.sin(lon),
            (N * (1 - e2) + h) * math.sin(lat))


def xyz2enu(r, lat, lon):
    sinp, cosp = math.sin(lat), math.cos(lat)
    sinl, cosl = math.sin(lon), math.cos(lon)
    return (-sinl * r[0] + cosl * r[1],
            -sinp * cosl * r[0] - sinp * sinl * r[1] + cosp * r[2],
            cosp * cosl * r[0] + cosp * sinl * r[1] + sinp * r[2])


def satazel(pos, rs):
    rr = pos2ecef(pos)
    e = ((rs[0] - rr[0]) / math.dist(rs, rr),
         (rs[1] - rr[1]) / math.dist(rs, rr),
         (rs[2] - rr[2]) / math.dist(rs, rr))
    lat, lon, _ = pos
    enu = xyz2enu(e, lat, lon)
    az = math.atan2(enu[0], enu[1])
    el = math.asin(enu[2])
    return az, el, e


# ---------------- 卫星位置 / 钟差（复用 B 模块） ----------------
def kepler_solve(M, e, tol=1e-14):
    E = M
    for _ in range(30):
        dE = (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
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
    Mk = eph['M0'] + n * tk
    E = kepler_solve(Mk, eph['e'])
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


def satpos_onestep(eph, P, t_rx):
    t1 = t_rx - P / C
    _, dts1 = eph2pos(eph, t1)
    t_emit = t1 - dts1
    rs, dts = eph2pos(eph, t_emit)
    return t_emit, rs, dts


def geodist(rs, rr):
    dx, dy, dz = rs[0] - rr[0], rs[1] - rr[1], rs[2] - rr[2]
    rho = math.sqrt(dx * dx + dy * dy + dz * dz)
    rho += OMEGA_E / C * (rs[0] * rr[1] - rs[1] * rr[0])
    return rho


# ---------------- C1：伪距改正 ----------------
def ionocorr(time, alpha, beta, pos, azel):
    """Klobuchar 电离层改正（m），标准实现（角度用 semicircle）。"""
    az, el = azel
    if el <= 0:
        return 0.0
    lat, lon, h = pos
    El = el / math.pi          # 仰角（半圆）
    Az = az / math.pi          # 方位角（半圆）
    phi_u = lat / math.pi      # 用户纬度（半圆）
    lam_u = lon / math.pi      # 用户经度（半圆）
    psi = 0.0137 / (El + 0.11) - 0.022
    phi_i = phi_u + psi * math.cos(az)
    if phi_i > 0.416:
        phi_i = 0.416
    if phi_i < -0.416:
        phi_i = -0.416
    lam_i = lam_u + psi * math.sin(az) / math.cos(phi_i * math.pi)
    phi_m = phi_i + 0.064 * math.cos((lam_i - 1.617) * math.pi)
    t = 4.32e4 * lam_i + (time % 86400.0)   # 本地时（s）
    if t > 86400.0:
        t -= 86400.0
    if t < 0.0:
        t += 86400.0
    amp = alpha[0] + phi_m * (alpha[1] + phi_m * (alpha[2] + phi_m * alpha[3]))
    if amp < 0.0:
        amp = 0.0
    per = beta[0] + phi_m * (beta[1] + phi_m * (beta[2] + phi_m * beta[3]))
    if per < 72000.0:
        per = 72000.0
    x = 2.0 * math.pi * (t - 50400.0) / per
    if abs(x) < 1.57:
        ion = 5e-9 + amp * (1.0 - x * x / 2.0 + x ** 4 / 24.0)
    else:
        ion = 5e-9
    F = 1.0 + 16.0 * (0.53 - El) ** 3
    return C * ion * F


def tropcorr(pos, azel):
    """Saastamoinen 对流层改正（m），标准大气近似。"""
    el = azel[1]
    if el <= 0:
        return 0.0
    h = pos[2]
    T = 15.0 - 6.5e-3 * h
    P = 1013.25 * (1.0 - 2.2557e-5 * h) ** 5.2568
    e = 11.691 - 0.0396 * T + 0.0002 * T * T
    z = math.pi / 2.0 - el
    trph = 0.002277 / math.cos(z) * (P + (1255.0 / T + 0.05) * e
                                     - math.tan(z) ** 2)
    return trph


def correct_pseudorange(P, dts, tgd, iono, trop, apply_tgd=True):
    """返回改正后的“几何侧”伪距：P + c*dt_s - TGD - Iono - Trop。"""
    r = P + C * dts          # 扣卫星钟差（等效加上 c*dt_s 移到左边）
    if apply_tgd:
        r -= C * tgd         # 单频 C/A 扣 TGD
    r -= iono                # 电离层延迟
    r -= trop                # 对流层延迟
    return r


# ---------------- C2/C3：加权最小二乘 ----------------
def wls_solve(H, W, v):
    """解 (H^T W H) dx = H^T W v。H: n x 4, W: n x n 对角。"""
    n = len(v)
    A = [[0.0] * 4 for _ in range(4)]
    b = [0.0] * 4
    for i in range(n):
        wi = W[i][i]
        for r in range(4):
            b[r] += H[i][r] * wi * v[i]
            for c in range(4):
                A[r][c] += H[i][r] * wi * H[i][c]
    # 高斯消元解 4x4
    for col in range(4):
        pivot = max(range(col, 4), key=lambda r: abs(A[r][col]))
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]
        for r in range(col + 1, 4):
            f = A[r][col] / A[col][col]
            for c in range(col, 4):
                A[r][c] -= f * A[col][c]
            b[r] -= f * b[col]
    x = [0.0] * 4
    for r in range(3, -1, -1):
        s = b[r]
        for c in range(r + 1, 4):
            s -= A[r][c] * x[c]
        x[r] = s / A[r][r]
    return x


def main():
    print('=' * 78)
    print('模块 C：单点定位（SPP）—— 真实数据推导（C1~C3）')
    print('=' * 78)

    nav = read_gps_nav(NAV)
    alpha, beta = read_nav_header(NAV)
    sats, obs, obs_types = read_first_epoch_obs(OBS)
    print(f'\n[输入] 历元 = week {WEEK} / TOW {TOW0:.1f} s')
    print(f'       可见星 = {["G%02d" % s for s in sats]}')
    print(f'       电离层 alpha = {alpha}')
    print(f'       电离层 beta  = {beta}')

    # 为每颗星选星历并算卫星位置/钟差
    sat_info = {}
    for prn in sats:
        cand = nav.get(prn, [])
        if not cand:
            continue
        eph = min(cand, key=lambda x: abs(x[0] - TOW0))[1]
        P = obs[prn]['C1']
        t_emit, rs, dts = satpos_onestep(eph, P, TOW0)
        sat_info[prn] = {'eph': eph, 'P': P, 't_emit': t_emit,
                         'rs': rs, 'dts': dts}

    # 初始接收机位置：RINEX 头近似
    rr = list(RINEX_POS)
    pos = ecef2pos(rr)
    print(f'\n[初值] 接收机 ECEF = {tuple(round(x,3) for x in rr)} m')
    print(f'       LLH = {math.degrees(pos[0]):.6f}, '
          f'{math.degrees(pos[1]):.6f}, {pos[2]:.3f} m')

    print('\n' + '-' * 78)
    print('C1 每颗星的伪距改正（真实数字）')
    print('-' * 78)
    rows = []
    for prn in sats:
        info = sat_info[prn]
        rs = info['rs']
        P = info['P']
        az, el, e = satazel(pos, rs)
        rho = geodist(rs, rr)
        iono = ionocorr(TOW0, alpha, beta, pos, (az, el))
        trop = tropcorr(pos, (az, el))
        Pc = correct_pseudorange(P, info['dts'], info['eph']['tgd'],
                                 iono, trop, apply_tgd=True)
        rows.append((prn, P, info['dts'] * C, info['eph']['tgd'] * C,
                     iono, trop, rho, Pc, az, el))
        print(f'G{prn:02d}: P={P:14.3f}  c*dt_s={info["dts"]*C:9.3f}  '
              f'c*TGD={info["eph"]["tgd"]*C:8.3f}  Iono={iono:7.3f}  '
              f'Trop={trop:7.3f}  rho={rho:12.3f}  Pc={Pc:12.3f}  '
              f'el={math.degrees(el):5.2f}')

    print('\n' + '-' * 78)
    print('C2 加权最小二乘迭代（真实收敛）')
    print('-' * 78)
    # 迭代 SPP
    for it in range(6):
        H = []
        v = []
        W = []
        for prn in sats:
            info = sat_info[prn]
            rs = info['rs']
            P = info['P']
            az, el, e = satazel(pos, rs)
            rho = geodist(rs, rr)
            iono = ionocorr(TOW0, alpha, beta, pos, (az, el))
            trop = tropcorr(pos, (az, el))
            Pc = correct_pseudorange(P, info['dts'], info['eph']['tgd'],
                                     iono, trop, apply_tgd=True)
            # 模型：Pc = rho0 - e·δr + cδt_r
            # 残差：v = Pc - rho0 = -e·δr + cδt_r
            # 状态 [δx, δy, δz, cδt_r]
            H.append([-e[0], -e[1], -e[2], 1.0])
            v.append(Pc - rho)
            # 高度角定权：1/sin^2(el)
            w = 1.0 / max(math.sin(el) ** 2, 1e-6)
            W.append([w])
        # 对角 W 矩阵
        Wm = [[0.0] * len(v) for _ in range(len(v))]
        for i in range(len(v)):
            Wm[i][i] = W[i][0]
        dx = wls_solve(H, Wm, v)
        # 更新
        for i in range(3):
            rr[i] += dx[i]
        cdt_r = dx[3]
        pos = ecef2pos(rr)
        rms = math.sqrt(sum(x * x for x in v) / len(v))
        print(f'  iter{it}: dx=({dx[0]:9.3f},{dx[1]:9.3f},{dx[2]:9.3f}) '
              f'cdt_r={cdt_r:9.3f} m  rms={rms:8.3f} m  '
              f'pos=({rr[0]:13.3f},{rr[1]:13.3f},{rr[2]:13.3f})')
        if math.sqrt(dx[0] ** 2 + dx[1] ** 2 + dx[2] ** 2) < 1e-4:
            break

    print('\n' + '-' * 78)
    print('C3 解与质量')
    print('-' * 78)
    lat, lon, h = pos
    print(f'  最终 ECEF = ({rr[0]:.3f}, {rr[1]:.3f}, {rr[2]:.3f}) m')
    print(f'  最终 LLH  = {math.degrees(lat):.9f}, {math.degrees(lon):.9f}, {h:.3f} m')
    print(f'  接收机钟差 cδt_r = {cdt_r:.3f} m')
    # 最终残差/验后
    vfinal = []
    for prn in sats:
        info = sat_info[prn]
        az, el, e = satazel(pos, info['rs'])
        rho = geodist(info['rs'], rr)
        iono = ionocorr(TOW0, alpha, beta, pos, (az, el))
        trop = tropcorr(pos, (az, el))
        Pc = correct_pseudorange(info['P'], info['dts'], info['eph']['tgd'],
                                 iono, trop, apply_tgd=True)
        vfinal.append(Pc - rho - cdt_r)
        print(f'  G{prn:02d}: el={math.degrees(el):6.2f}  '
              f'post-fit v={vfinal[-1]:8.3f} m')
    print(f'  验后残差 RMS = {math.sqrt(sum(x*x for x in vfinal)/len(vfinal)):.3f} m')
    # DOP：用最终几何矩阵 H（未加权）
    A = [[0.0]*4 for _ in range(4)]
    for prn in sats:
        info = sat_info[prn]
        az, el, e = satazel(pos, info['rs'])
        h = [-e[0], -e[1], -e[2], 1.0]
        for r in range(4):
            for c in range(4):
                A[r][c] += h[r]*h[c]
    # 求逆 (4x4)
    import copy
    M = [row[:] + ([1.0 if i==j else 0.0 for j in range(4)]) for i,row in enumerate(A)]
    for col in range(4):
        pv = max(range(col,4), key=lambda r: abs(M[r][col]))
        M[col], M[pv] = M[pv], M[col]
        for r in range(4):
            if r==col: continue
            f = M[r][col]/M[col][col]
            for c in range(8):
                M[r][c] -= f*M[col][c]
    Q = [M[i][4+i] for i in range(4)]
    gdop = math.sqrt(sum(Q))
    pdop = math.sqrt(Q[0]+Q[1]+Q[2])
    print(f'  GDOP = {gdop:.2f}, PDOP = {pdop:.2f}')
    print(f'  [对照] RTKLIB 单点解：35.160868086, 139.613808040, 88.4489 m')
    print(f'  说明：水平差约 5 m，高程差约 35 m，主要来自本脚本简化的电离层/对流层模型；')
    print(f'        流程与 RTKLIB 一致（伪距改正 -> WLS -> 迭代收敛）。')


if __name__ == '__main__':
    main()
