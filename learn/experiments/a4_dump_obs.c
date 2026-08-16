/* a4_dump_obs.c
 *
 * 实验目的（A4 观测数据流）：
 *   直接调用 readrnx()，把 RINEX 观测文件解码成 obsd_t，
 *   打印首个历元的观测值，验证：
 *     1) 单位：P->m, L->cycle, D->Hz, SNR->0.001dBHz
 *     2) code[] 是 CODE_??? 索引（用 code2obs() 转字符串）
 *     3) sat 是 RTKLIB 内部卫星号（用 satno2id() 转 "G03"）
 *     4) L*lambda 与 P 的差应接近"整数周 × 波长"（含整周模糊度）
 *
 * 编译（复用 rnx2rtkp 已编译对象，避免重复编译整个库）：
 *   gcc -Wall -O2 -DTRACE -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN -DNFREQ=5 \
 *       -I src -o learn/experiments/a4_dump_obs learn/experiments/a4_dump_obs.c \
 *       $(ls app/consapp/rnx2rtkp/gcc/*.o | grep -v rnx2rtkp.o) \
 *       lib/iers/gcc/iers.a -lgfortran -lm -lrt
 *
 * 运行：
 *   learn/experiments/a4_dump_obs test/data/rinex/07590920.05o \
 *       test/data/rinex/07590920.05n
 */
#include "rtklib.h"
#include <stdio.h>
#include <string.h>

/* postpos.o 需要的回调（rnx2rtkp.c 中也这样实现） */
extern int showmsg(const char *format, ...) { return 0; }
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

static void print_obs(const obs_t *obs, const nav_t *nav, int maxepoch)
{
    gtime_t t0 = obs->data[0].time;
    int i, f, cnt = 0;
    char id[8], code[8];

    for (i = 0; i < obs->n; i++) {
        if (timediff(obs->data[i].time, t0) != 0.0) {
            if (++cnt >= maxepoch) break;
            t0 = obs->data[i].time;
            printf("\n=== 下一个历元 ===\n");
        }
        if (cnt == 0) {
            printf("=== 历元 %s GPST ===\n", time_str(t0, 3));
            cnt++; /* 首个历元标记，避免重复 */
        }
        satno2id(obs->data[i].sat, id);
        printf("sat=%s (内部号 %2d)\n", id, obs->data[i].sat);
        for (f = 0; f < NFREQ + NEXOBS; f++) {
            double lam, nxlam;
            if (obs->data[i].P[f] == 0.0 && obs->data[i].L[f] == 0.0) continue;
            strcpy(code, code2obs(obs->data[i].code[f]));
            lam = sat2freq(obs->data[i].sat, obs->data[i].code[f], nav) > 0.0 ?
                  CLIGHT / sat2freq(obs->data[i].sat, obs->data[i].code[f], nav) : 0.0;
            nxlam = lam > 0.0 && obs->data[i].P[f] > 0.0 ?
                    (obs->data[i].P[f] - obs->data[i].L[f] * lam) / lam : 0.0;
            printf("  f%d code=%-3s  P=%13.3f m  L=%15.3f cyc  "
                   "D=%10.1f Hz  SNR=%6.1f dBHz  LLI=%d\n",
                   f, code, obs->data[i].P[f], obs->data[i].L[f],
                   obs->data[i].D[f], obs->data[i].SNR[f] * 0.001,
                   obs->data[i].LLI[f]);
            if (lam > 0.0) {
                printf("           lambda=%.4f m  L*lambda=%13.3f m  "
                       "  (P - L*lambda)/lambda = %.4f (应接近整数=整周N)\n",
                       lam, obs->data[i].L[f] * lam, nxlam);
            }
        }
    }
}

int main(int argc, char **argv)
{
    obs_t obs = {0};
    nav_t nav = {0};
    sta_t sta;
    int i;

    if (argc < 2) {
        fprintf(stderr, "usage: %s obsfile [navfile ...]\n", argv[0]);
        return -1;
    }
    for (i = 1; i < argc; i++) {
        /* rcv=1: 观测文件；nav 文件时 obs=NULL */
        if (readrnx(argv[i], 1, "", &obs, &nav, &sta) <= 0) {
            fprintf(stderr, "readrnx error: %s\n", argv[i]);
            return -1;
        }
    }
    sortobs(&obs);
    printf("total obs records = %d\n", obs.n);

    print_obs(&obs, &nav, 2);

    freeobs(&obs);
    freenav(&nav, 1);
    return 0;
}
