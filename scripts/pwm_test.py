import csv, math, random
from collections import Counter

D = '/home/hackathon/hackathon-data/her2-trastuzumab/'
BASE = 'WGGDGFYAMD'
AAs = 'ACDEFGHIKLMNPQRSTVWY'

# ---- build PWM log-odds (high vs low) from FACS ----
hi = [Counter() for _ in range(10)]
lo = [Counter() for _ in range(10)]
facs_class = {}
for r in csv.DictReader(open(D + 'facs_unz/tz_her2_facs_3classes_524346.csv')):
    s, c = r['cdrh3'], r['binding_class']
    facs_class[s] = c
    if c == 'high':
        for i, a in enumerate(s):
            hi[i][a] += 1
    elif c == 'low':
        for i, a in enumerate(s):
            lo[i][a] += 1

W = []
for i in range(10):
    th, tl = sum(hi[i].values()), sum(lo[i].values())
    W.append({a: math.log(((hi[i][a] + 1) / (th + 20)) / ((lo[i][a] + 1) / (tl + 20))) for a in AAs})


def pwm(s):
    return sum(W[i].get(a, 0.0) for i, a in enumerate(s)) if len(s) == 10 else None


def ed(a, b):
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else 0.0


# ---- SPR set: trastuzumab-framework 10-mers with a KD, excluding the 8 reference rows ----
spr = list(csv.DictReader(open(D + 'trastuzumab_her2_spr_all_768.csv')))
rows = [r for r in spr if r['kd_nM'] and len(r['cdrh3']) == 10 and r['cdrh3'] != BASE
        and 'GFNIKDTYIH' in r['vh_sequence']]
print('SPR 10-mer variants with KD (reference excluded):', len(rows))

seqs = [r['cdrh3'] for r in rows]
logkd = [math.log10(float(r['kd_nM'])) for r in rows]
p = [pwm(s) for s in seqs]
d = [ed(s, BASE) for s in seqs]

print('\nSpearman vs log10 KD (negative = good, predicts tighter binding):')
print('  PWM log-odds       %+.3f' % spearman(p, logkd))
print('  edit distance      %+.3f' % spearman(d, logkd))
# combined: residualise PWM against distance is overkill; just report within-distance
print('\nwithin edit-distance strata (PWM only, where n>=8):')
from collections import defaultdict
st = defaultdict(list)
for s, k, pp, dd in zip(seqs, logkd, p, d):
    st[dd].append((pp, k))
for dd in sorted(st):
    v = st[dd]
    if len(v) >= 8:
        print('  dist %2d  n=%2d  rho=%+.3f' % (dd, len(v), spearman([a for a, _ in v], [b for _, b in v])))

# ---- how well does PWM separate the top decile? ----
order = sorted(range(len(seqs)), key=lambda i: -p[i])
top = [10 ** logkd[i] for i in order[:19]]
bot = [10 ** logkd[i] for i in order[-19:]]
print('\ntop decile by PWM:    median KD %8.1f nM   best %.3f' % (sorted(top)[len(top) // 2], min(top)))
print('bottom decile by PWM: median KD %8.1f nM   best %.3f' % (sorted(bot)[len(bot) // 2], min(bot)))

# ---- enrichment check: does PWM rank the 7 known sub-0.18nM binders highly? ----
good = {'WGSQGFYAYD', 'WGSDGWYAMD', 'WGGQGFYAMD', 'WGGSGPYLFD', 'WGGEGFYAMD', 'WGGDGFYLMD', 'WGGDGFYEMD'}
ranked = [seqs[i] for i in order]
print('\nrank (of %d) of the 7 sub-reference binders by PWM:' % len(ranked))
for g in good:
    if g in ranked:
        print('   %-12s rank %3d   dist %d' % (g, ranked.index(g) + 1, ed(g, BASE)))

# ---- novelty-window subset: >=4 edits, the designs that would pass the filter ----
nov = [(s, k) for s, k in zip(seqs, logkd) if ed(s, BASE) >= 4]
print('\nnovelty-passing (>=4 edits) SPR variants with KD: %d' % len(nov))
kds = sorted(10 ** k for _, k in nov)
print('   median %.1f nM, best %.3f nM, n under 1 nM: %d' % (kds[len(kds) // 2], kds[0], sum(1 for v in kds if v < 1)))
pn = [pwm(s) for s, _ in nov]
print('   PWM Spearman within this window: %+.3f' % spearman(pn, [k for _, k in nov]))
