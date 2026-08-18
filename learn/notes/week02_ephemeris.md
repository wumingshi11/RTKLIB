# 第 2 周笔记：如何计算卫星位置与钟差

> 目标：讲清 RTKLIB 从广播星历出发，得到"信号发射时刻"的卫星 ECEF 位置 `rs[]` 和钟差 `dts[]` 的完整流程。
>
> 对应源码：`src/ephemeris.c`、`src/pntpos.c`、`src/rtklib.h`

---

## 0. 一张调用图

```text
pntpos()                              src/pntpos.c
  └─ satposs(sol->time, obs, n, nav, opt_.sateph, rs, dts, var, svh)
       │  ① 由伪距反推信号发射时刻   time[i] = obs[i].time - P/c
       │  ② ephclk() 用广播星历算卫星钟差 → 修正发射时刻
       └─ satpos(time[i], teph, sat, ephopt, nav, rs, dts, var, svh)
            按 ephopt 分发：
            EPHOPT_BRDC   → ephpos()            ← 广播星历（最常用）
            EPHOPT_SBAS   → satpos_sbas()
            EPHOPT_SSRAPC → satpos_ssr()
            EPHOPT_SSRCOM → satpos_ssr()
            EPHOPT_PREC   → peph2pos()          ← 精密星历
```

`ephpos()` 内部再按卫星系统分发：

```text
ephpos()
  ├─ GPS / GAL / QZS / CMP / IRN → seleph() + eph2pos()   开普勒解析法
  ├─ GLO                          → selgeph() + geph2pos() 数值积分法
  └─ SBS                          → selseph() + seph2pos()
  最后：t 与 t+1ms 各算一次，差分求 速度 rs[3:6] 和 钟漂 dts[1]
```

---

## 1. 数据流与输出约定

`satposs()` 的输出（注释里写得很清楚）：

| 变量 | 含义 | 单位 |
|---|---|---|
| `rs[0:3]+i*6` | 第 i 颗卫星位置 {x,y,z} | m (ECEF) |
| `rs[3:6]+i*6` | 第 i 颗卫星速度 {vx,vy,vz} | m/s |
| `dts[0]+i*2` | 卫星钟差 bias | s |
| `dts[1]+i*2` | 卫星钟漂 drift | s/s |
| `var[i]` | 位置+钟差方差 | m² |
| `svh[i]` | 健康标志（-1: 无星历可用） | — |

**关键约定（易错）**：
- 位置/钟差都是**信号发射时刻**的值，不是接收时刻。
- 位置参考**天线相位中心**（广播星历给的就是质心，RTKLIB 后续在 `satantoff()` 里改正）。
- 钟差不含码偏差（TGD/BGD），由 `obsd_t` 观测的码需要单独做码偏差改正。
- 矩阵**列主序**：`rs` 用 `mat(6,n)`，所以第 i 颗卫星是 `rs[j + i*6]`。

### 为什么必须先算发射时刻？

观测方程中的伪距 $P_r^s$ 是接收时刻 $t_r$ 记录下的观测量，但卫星位置必须取**信号发射时刻** $t^s$ 的值：

$$t^s = t_r - \frac{P_r^s}{c}$$

RTKLIB 分两步：
1. `time[i] = timeadd(obs[i].time, -pr/CLIGHT)` —— 先用伪距粗估发射时刻
2. `ephclk(time[i], ...)` 算出卫星钟差 `dt`，再 `time[i] = timeadd(time[i], -dt)` —— 扣除卫星钟差

> 这也是学习计划里"常见误区：把接收时刻当发射时刻"的源码落点。
> 注意：任何伪距都没有时，`satposs()` 直接跳过该卫星（`if (j>=NFREQ) continue;`）。

---

## 2. 广播星历 `eph_t`：用到的开普勒参数

`src/rtklib.h` 中 `eph_t` 结构体（GPS/GAL/QZS/CMP/IRN 通用）：

```c
double A;      /* 轨道长半轴 √A 的平方 A (m)   */
double e;      /* 偏心率                         */
double i0;     /* 参考时刻轨道倾角 (rad)         */
double OMG0;   /* 参考时刻升交点赤经 Ω0 (rad)    */
double omg;    /* 近地点幅角 ω (rad)             */
double M0;     /* 参考时刻平近点角 M0 (rad)      */
double deln;   /* 平均角速度改正 Δn (rad/s)      */
double OMGd;   /* 升交点赤经变化率 Ω̇ (rad/s)    */
double idot;   /* 轨道倾角变化率 İ (rad/s)       */
double crc,crs,cuc,cus,cic,cis; /* 调和改正系数   */
double toes;   /* Toe 在周内秒 (s)               */
double f0,f1,f2; /* 卫星钟差系数 af0,af1,af2     */
double toe,toc;  /* 星历参考时刻 / 钟参考时刻     */
```

**RINEX 导航电文 → eph_t** 由 `rinex.c` 的 `decode_eph()` + `add_eph()` 完成（第 2 周另一块精读内容）。

---

## 3. `eph2pos()`：GPS 广播星历 → 卫星位置的 11 步

这是最核心的一段，把下面的公式和 `src/ephemeris.c` 的代码逐行对照。

### 第 1 步：时间差

$$t_k = t - t_{oe}$$

代码：`tk = timediff(time, eph->toe);`

### 第 2 步：选系统常数

| 系统 | 引力常数 μ (m³/s²) | 地球自转角速度 Ω̇e (rad/s) |
|---|---|---|
| GPS | MU_GPS = 3.986005e14 | OMGE = 7.2921151467e-5 |
| GAL | MU_GAL | OMGE_GAL |
| CMP | MU_CMP | OMGE_CMP |

### 第 3 步：平均角速度 + 平近点角

$$n = \sqrt{\frac{\mu}{A^3}} + \Delta n$$

$$M = M_0 + n\, t_k$$

代码：`M = eph->M0 + (sqrt(mu/(A*A*A)) + eph->deln) * tk;`

### 第 4 步：开普勒方程牛顿迭代求偏近点角 E

$$E - e\sin E = M,\qquad E_{k+1} = E_k - \frac{E_k - e\sin E_k - M}{1 - e\cos E_k}$$

代码：

```c
for (n=0, E=M, Ek=0.0; fabs(E-Ek)>RTOL_KEPLER && n<MAX_ITER_KEPLER; n++) {
    Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
}
if (n>=MAX_ITER_KEPLER) { trace(2,"kepler iteration overflow..."); return; }
```

> 收敛容差 `RTOL_KEPLER=1E-13`，最多 `MAX_ITER_KEPLER=30` 次。迭代不收敛说明星历数据有问题，直接 `return`（此时 `rs/dts` 保持 0）。

### 第 5~7 步：纬度幅角 u、矢径 r、倾角 i（原始值）

$$\tan u = \frac{\sqrt{1-e^2}\sin E}{\cos E - e},\qquad u = \arctan2(\cdots) + \omega$$

$$r = A\,(1 - e\cos E)$$

$$i = i_0 + \dot{i}\, t_k$$

代码：

```c
u  = atan2(sqrt(1.0-e*e)*sinE, cosE-e) + eph->omg;
r  = eph->A*(1.0-eph->e*cosE);
i  = eph->i0 + eph->idot*tk;
```

### 第 8 步：二阶调和改正（6 个系数）

$$u \leftarrow u + C_{us}\sin 2u + C_{uc}\cos 2u$$

$$r \leftarrow r + C_{rs}\sin 2u + C_{rc}\cos 2u$$

$$i \leftarrow i + C_{is}\sin 2u + C_{ic}\cos 2u$$

代码：

```c
sin2u=sin(2.0*u); cos2u=cos(2.0*u);
u += eph->cus*sin2u + eph->cuc*cos2u;
r += eph->crs*sin2u + eph->crc*cos2u;
i += eph->cis*sin2u + eph->cic*cos2u;
```

### 第 9 步：轨道平面内坐标

$$x = r\cos u,\qquad y = r\sin u$$

### 第 10 步：旋转到 ECEF（含地球自转 Sagnac 改正）

升交点赤经（注意地球自转的两种写法，GEO 卫星特殊处理）：

$$\Omega = \Omega_0 + (\dot\Omega - \dot\Omega_e)\,t_k - \dot\Omega_e\, t_{oe}$$

ECEF 坐标：

$$X = x\cos\Omega - y\cos i\sin\Omega$$

$$Y = x\sin\Omega + y\cos i\cos\Omega$$

$$Z = y\sin i$$

代码：

```c
O    = eph->OMG0 + (eph->OMGd - omge)*tk - omge*eph->toes;
rs[0]= x*cosO - y*cosi*sinO;
rs[1]= x*sinO + y*cosi*cosO;
rs[2]= y*sin(i);
```

> 为什么有 $-\dot\Omega_e\,t_k$？因为 ECEF 系随地球自转，而开普勒轨道是在**惯性系**里描述的，必须扣除地球在 $t_k$ 内的自转角度，否则视线方向会有几十米级误差（赤道处约 $R_e\,\dot\Omega_e\,t_k \approx$ 数 m~数十 m）。

> 北斗 GEO 卫星（`prn<=5 || prn>=59`）走另一个分支，还要再做一次绕 z 轴的旋转（`COS_5/SIN_5`，对应 5° 倾角的 GEO 轨道面旋转），这是 BDS 特有的处理。

### 第 11 步：卫星钟差 + 相对论改正

广播钟差多项式（此时 $t_k = t - t_{oc}$）：

$$dt^s = a_{f0} + a_{f1}\,t_k + a_{f2}\,t_k^2$$

相对论改正（由轨道偏心率引起，约为 -7 µs ~ +7 µs 的周期项）：

$$\Delta dt_{rel} = -\frac{2\sqrt{\mu A}}{c^2}\,e\sin E$$

代码：

```c
tk   = timediff(time, eph->toc);
*dts = eph->f0 + eph->f1*tk + eph->f2*tk*tk;
*dts -= 2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);  /* relativity */
```

方差：`*var = var_uraeph(sys, eph->sva);` —— 由 URA 指数换算成 m²。

---

## 4. 速度与钟漂：1 ms 差分近似（`ephpos()`）

`ephpos()` 对同一条星历算两次：

```c
eph2pos(time,        eph, rs,  dts, var);
time = timeadd(time, tt);                    /* tt = 1E-3 s */
eph2pos(time,        eph, rst, dtst, var);
/* 卫星速度 */
for (i=0;i<3;i++) rs[i+3] = (rst[i]-rs[i])/tt;
/* 钟漂 */
dts[1] = (dtst[0]-dts[0])/tt;
```

数值差分即可满足定位需求，不需要解析求导。

---

## 5. 其他系统的卫星位置

### GLONASS：`geph2pos()` —— 4 阶 Runge-Kutta 数值积分

GLONASS 广播星历不直接给开普勒参数，而是给**参考时刻的位置/速度/加速度**（`geph_t.pos[3], vel[3], acc[3]`），需要积分运动方程：

$$m\ddot{\mathbf r} = -\mu\frac{\mathbf r}{r^3} + \text{J}_2\text{ 项} + \text{科里奥利项（地球旋转）}$$

`deq()` 里实现右端：

```c
a = 1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3;   /* 3/2*J2*μ*Re²/r⁵ */
b = 5.0*x[2]*x[2]/r2;                        /* 5z²/r² */
c = -MU_GLO/r3 - a*(1.0-b);
xdot[3] = (c+omg2)*x[0] + 2.0*OMGE_GLO*x[4] + acc[0];
xdot[4] = (c+omg2)*x[1] - 2.0*OMGE_GLO*x[3] + acc[1];
xdot[5] = (c-2.0*a)*x[2] + acc[2];
```

`glorbit()` 是经典 RK4：`k1+2k2+2k3+k4`，从 `geph->toe` 积分到目标时刻。

### 精密星历：`peph2pos()` —— 插值

精密星历是**表格化**的（每 5/15/30 min 一档位置+钟差），需要用插值多项式内插到目标时刻（RTKLIB 用 Lagrange 多项式）。这里只提入口，细节在 `peph2pos` 与 `interppol()`。

---

## 6. 调用 `satpos` 之后

`pntpos()` 拿到 `rs/dts` 后：

1. `estpos()`：算几何距离 `geodist(rs, rr)`（含地球自转 Sagnac 改正）、视线方向、伪距残差，做加权最小二乘 → 接收机位置。
2. `estvel()`：用多普勒 + 卫星速度算接收机速度。

所以卫星位置是**整个 SPP 的第一个前置环节**，位置错了后面全错。

---

## 7. 卫星位置/钟差计算的偏差来源与量级（学习重点：不用抠开普勒，但要懂误差）

### 7.1 一句话定位

开普勒解析法是 ICD 里的成熟算法，按公式照做即可；**学习重点不是公式本身，而是"每一步漏掉会引入多大偏差"以及这些偏差如何进入最终定位解**。

### 7.2 三类偏差

**A. 广播星历自身质量（单机内无法消除）**

| 误差源 | 典型量级 | 说明 |
|---|---|---|
| 广播轨道误差 | 0.5~2 m | 地面监测站预测+注入，用 URA 指数 `sva` 描述 |
| 广播钟差误差 | ~1~2 ns（0.3~0.6 m） | af0 拟合残差；`var_uraeph()` 只给统计方差，不给真值 |

**B. 模型改正项（漏掉就是系统偏差，`eph2pos()` 内已处理）**

| 改正项 | 量级 | 处理位置 |
|---|---|---|
| 相对论周期项 | 最大 ~7 m（±23 ns） | `eph2pos()` 末尾 `-2√(μA)·e·sinE/c²` |
| 地球自转 Sagnac | 视线方向可达数十 m | 升交点赤经 $-Ω̇_e·t_k$ 项 + `geodist()` |

值得会手算的量级：
- 相对论周期项：$2\sqrt{\mu A}\,e/c^2 \approx \dfrac{2\sqrt{3.986\times10^{14}\times2.656\times10^{7}}\times0.01}{(3\times10^8)^2} \approx 23\,\text{ns} \approx 7\,\text{m}$
- Sagnac 项：$R_e\,\dot\Omega_e\,\tau \approx 6.37\times10^{6}\times7.29\times10^{-5}\times0.07 \approx 32\,\text{m}$（信号传播约 70 ms，赤道附近切线方向）

**C. 传播与硬件类偏差（进入伪距预测，不是卫星位置本身）**

| 项 | 量级 | 处理位置 |
|---|---|---|
| 天线相位中心 PCO/PCV | 0.5~1.5 m | `satantoff()`（质心→相位中心） |
| 电离层延迟 | 天顶 5~15 m，低仰角更大 | 模型/双频消（`estpos` 内） |
| 对流层延迟 | 天顶 ~2.5 m，低仰角可达 20 m+ | 模型/状态估计 |
| 多路径 | 米级、难建模 | 天线/接收机设计 |

### 7.3 偏差如何进入定位解：几何因子 DOP

定位误差 ≈ 观测误差 × DOP（取决于卫星几何构型）。卫星位置/钟差误差沿视线方向投影后，与接收机钟差耦合求解，最终决定 SPP 精度（水平约 2~5 m）。

### 7.4 为什么这是"差分定位的起点"（衔接第 4~6 周）

- 星历误差、卫星钟差、电离层延迟是**空间相关的**：相距不远的两个测站看同一颗卫星，这些误差几乎相同。
- 多路径、热噪声是**空间不相关**的，差分消不掉。
- → 单差/双差（第 4~6 周）的核心就是消掉这些公共误差。**理解"哪类误差空间相关、哪类不相关"，比背公式重要得多。**

### 7.5 电离层/对流层模型：算的是"积分/拟合"，不是"厚度"

**电离层**（色散介质，延迟 ∝ 1/f²）：
- 关键量是 **TEC = ∫Ne·ds**（电子密度沿路径积分），不是层厚。1 TECU = 10¹⁶ e⁻/m²，L1 上 1 TECU ≈ 0.163 m。
- Klobuchar 模型（`ionmodel()`，`IONOOPT_BRDC`）：用广播星历 8 参数 {α0~α3, β0~β3} 拟合一条**余弦曲线**近似 TEC 随地方时的变化 + 5 ns 常数底，是**经验拟合**不是物理厚度。
- `ionmapf()` 的"单层模型"：假设电离层集中在 350 km 薄壳（`HION=350`）做天顶→斜向投影，是**几何简化**，不是算厚度。

**对流层**（非色散，与频率无关）：
- 延迟 = ∫(n−1)·ds（折射率路径积分），分干/湿。
- Saastamoinen（`tropmodel()`，`TROPOPT_SAAS`）：标准大气算气压/温度/湿度 → 天顶干延迟 ZHD ≈ 2.3 m（看气压，好估）+ 湿延迟 ZWD（看水汽，难估）。也不是厚度。

**残余量（SPP 视角）**

| 项 | 模型 | 消除比例 | 残余 |
|---|---|---|---|
| 电离层 | Klobuchar（广播 8 参） | ~50–60% | 白天 2~10 m，夜间 <1~2 m（SPP 头号误差源） |
| 电离层 | 双频无电离层组合 IF | 一阶 >99.9% | 二阶 ~0.1~1 cm |
| 电离层 | IONEX 格网 | 高 | cm~几十 cm |
| 对流层 | Saastamoinen+映射 | 干延迟 ~90–95% | 天顶 ~0.1~0.3 m（湿延迟主导），低仰角放大 |

> 关键认知：Klobuchar 只能消一半多电离层 → 这就是差分定位/双频消必须存在的理由；它们不是靠"更好的模型"，而是靠"消除空间相关误差"（差分）或"频率相关性"（IF 组合）。

---

## 8. 验收 checklist（对应学习计划 2.4 数值验收）

1. 选测试数据 `test/data/rinex/07590920.05o` + `07590920.05n` 里同一颗 GPS 卫星。
2. 记录 `satposs()` 输出：发射时刻、`rs[0:6]`、`dts[0:2]`、`var`、`svh`。
3. 独立小程序复算 `eph2pos`：核对中间量 `tk, n, M, E, u, r, i, Ω, rs`。
4. 若不一致，优先检查：时间系统（GPST）、角度单位（rad）、`toe` 取整、调和改正的符号。
5. 分别关闭/保留卫星钟差、Sagnac 项，观察伪距预测变化量（应分别是 ~10⁵ m 量级和 ~数 m 量级）。

## 9. 常见误区

- **把接收时刻当发射时刻** → 卫星沿轨道移动 + 钟差，伪距预测误差可达数百米。
- **`dts` 忘了相对论改正** → 数米级常数偏差（GPS 广播星历已在参数中隐含部分，这里补的是周期项）。
- **`toe` 与 `toc` 混用** → 轨道用 `toe`，钟差用 `toc`，两个时刻不同！
- **忽略北斗 GEO 分支** → CMP PRN≤5 / ≥59 走旋转分支，用普通公式会错。
- **角度单位**：RINEX 里是弧度；`toes` 是秒，别拿周数当秒用。
