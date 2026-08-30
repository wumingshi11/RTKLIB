#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块 F：模糊度固定 —— 用一个小例子推导（F1~F3）
================================================
  F1 周跳与模糊度生命周期（概念 + D3 数据回顾）
  F2 整数最小二乘（ILS）原理：用 2 维强相关例子做整数搜索
  F3 从浮点到固定解：条件更新

运行：python3 learn/experiments/f_lambda_derive.py
"""
import math
import itertools


def mat_inv2(A):
    """2x2 矩阵求逆。"""
    a, b, c, d = A[0][0], A[0][1], A[1][0], A[1][1]
    det = a * d - b * c
    return [[d / det, -b / det], [-c / det, a / det]]


def quad_form(N, a, Qinv):
    """(N-a)^T Q^{-1} (N-a)"""
    d = [N[0] - a[0], N[1] - a[1]]
    return d[0] * (Qinv[0][0] * d[0] + Qinv[0][1] * d[1]) + \
           d[1] * (Qinv[1][0] * d[0] + Qinv[1][1] * d[1])


def main():
    print('=' * 78)
    print('模块 F：模糊度固定 —— 小例子推导（F1~F3）')
    print('=' * 78)

    # ---------- F1 ----------
    print('\n' + '-' * 78)
    print('F1 周跳与模糊度生命周期')
    print('-' * 78)
    print('  状态机：未观测 → 初始化(浮点) → 锁定(连续) → 失锁/周跳 → 重置 → 重初始化')
    print('  检测手段（来自 D 模块）：')
    print('    - GF 历元间差分（正常 cm 级，周跳跳 0.19/0.24 m 整数倍）')
    print('    - MW 比较（正常 <1 周，周跳跳整数周）')
    print('    - LLI 标志（接收机硬件直接给失锁/半周指示）')
    print('  周跳发生后：该卫星的模糊度不再是同一个整数 → 必须重置/重新初始化。')

    # ---------- F2 ----------
    print('\n' + '-' * 78)
    print('F2 整数最小二乘（ILS）—— 2 维强相关例子')
    print('-' * 78)
    # 浮点模糊度及其协方差（故意强相关，模拟真实 RTK 情况）
    a = [2.10, 3.05]           # 浮点模糊度（周）
    Q = [[1.0, 0.9], [0.9, 1.0]]  # 协方差（强相关）
    Qinv = mat_inv2(Q)
    print(f'  浮点模糊度 a = {a} 周')
    print(f'  协方差 Q = {Q}')
    print(f'  Q 的条件数 = {((1+0.9)/(1-0.9)):.1f}（强相关 → 条件数大）')

    # 在 a 附近整数格点搜索
    best = []
    for n1 in range(-10, 15):
        for n2 in range(-10, 15):
            q = quad_form((n1, n2), a, Qinv)
            best.append((q, n1, n2))
    best.sort()
    print('\n  最小几个整数候选（加权平方距离）：')
    for q, n1, n2 in best[:5]:
        print(f'    N=({n1}, {n2})  (N-a)^T Q^-1 (N-a) = {q:.3f}')
    q1, N1 = best[0][0], (best[0][1], best[0][2])
    q2, N2 = best[1][0], (best[1][1], best[1][2])
    print(f'\n  最优整数解 = {N1}, 距离 = {q1:.3f}')
    print(f'  次优整数解 = {N2}, 距离 = {q2:.3f}')
    print(f'  ratio = {q2/q1:.2f}（>3 通常认为可固定）')
    # 对比：直接四舍五入
    print(f'  对比：直接四舍五入 = ({round(a[0])}, {round(a[1])})，'
          f'距离 = {quad_form((round(a[0]), round(a[1])), a, Qinv):.3f}')
    print('  → 强相关时“逐个四舍五入”可能不是最优；ILS 要联合搜索。')

    # ---------- F3 ----------
    print('\n' + '-' * 78)
    print('F3 从浮点到固定解（条件更新）')
    print('-' * 78)
    # 假设位置/模糊度联合：x_float 位置，b_float 模糊度
    # 简化：位置 1 维 + 2 个模糊度
    x_float = [1.20]            # 位置浮点解
    Qxx = [[0.04]]              # 位置方差
    Qxb = [[0.02, -0.01]]       # 位置与模糊度协方差 (1x2)
    Qbb = Q                     # 模糊度协方差 2x2
    b_fixed = list(N1)          # 固定整数
    # 条件更新：x_fixed = x_float - Qxb Qbb^-1 (b_float - b_fixed)
    Qbb_inv = mat_inv2(Qbb)
    diff_b = [a[0] - b_fixed[0], a[1] - b_fixed[1]]
    # Qxb (1x2) * Qbb_inv (2x2) -> 1x2
    gain = [Qxb[0][0] * Qbb_inv[0][0] + Qxb[0][1] * Qbb_inv[1][0],
            Qxb[0][0] * Qbb_inv[0][1] + Qxb[0][1] * Qbb_inv[1][1]]
    correction = gain[0] * diff_b[0] + gain[1] * diff_b[1]
    x_fixed = x_float[0] - correction
    Qx_fixed = Qxx[0][0] - (gain[0] * Qxb[0][0] + gain[1] * Qxb[0][1])
    print(f'  位置浮点 x_float = {x_float[0]:.3f} m, 方差 = {Qxx[0][0]:.3f}')
    print(f'  固定模糊度 b_fixed = {b_fixed} 周')
    print(f'  条件更新：x_fixed = x_float - Qxb·Qbb⁻¹·(b_float - b_fixed)')
    print(f'    = {x_float[0]:.3f} - ({correction:.3f}) = {x_fixed:.3f} m')
    print(f'  固定后方差 = {Qx_fixed:.4f} m²（比浮点 {Qxx[0][0]:.3f} 更小）')
    print('  → 固定模糊度后，位置精度显著提升，这就是 RTK 固定的价值。')


if __name__ == '__main__':
    main()
