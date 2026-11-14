import csv, statistics as st
from collections import Counter, defaultdict

BASE = 'WGGDGFYAMD'
D = '/home/hackathon/hackathon-data/her2-trastuzumab/'

# ---------- SPR ----------
spr = list(csv.DictReader(open(D + 'trastuzumab_her2_spr_all_768.csv')))
print('SPR rows:', len(spr))
print('binding_result:', Counter(r['binding_result'] for r in spr).most_common())
print('cdrh3 len:', Counter(len(r['cdrh3']) for r in spr).most_common())
print('unique cdrh3:', len({r['cdrh3'] for r in spr}))
vh_counts = Counter(r['vh_sequence'][:0] or ('tz' if 'GFNIKDTYIH' in r['vh_sequence'] else 'other') for r in spr)
print('framework:', vh_counts.most_common())

kd = [(float(r['kd_nM']), r['cdrh3']) for r in spr if r['kd_nM']]
print('rows with KD:', len(kd))
vals = sorted(v for v, _ in kd)
qs = [vals[int(p * (len(vals) - 1))] for p in (0, .25, .5, .75, 1)]
print('KD nM min/q1/med/q3/max: %.3f %.2f %.2f %.1f %.1f' % tuple(qs))
print('KD under 1 nM:', sum(1 for v in vals if v < 1), ' under 0.2:', sum(1 for v in vals if v < 0.2))

ref = [v for v, c in kd if c == BASE]
print('trastuzumab (WGGDGFYAMD) rows:', sum(1 for r in spr if r['cdrh3'] == BASE), 'KDs:', sorted(ref))


def ed(a, b):
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))


# KD vs edit distance, trastuzumab-framework 10-mers only
buckets = defaultdict(list)
for v, c in kd:
    if len(c) == 10:
        buckets[ed(c, BASE)].append(v)
print('\nKD by CDRH3 edit distance from trastuzumab (10-mers with KD):')
print('  dist  n   median KD    best KD')
for d in sorted(buckets):
    b = sorted(buckets[d])
    print('  %4d %3d  %9.3f  %9.3f' % (d, len(b), st.median(b), b[0]))

# how often does a measured variant beat trastuzumab's worst reference KD?
if ref:
    thr = max(ref)
    better = [(v, c) for v, c in kd if v < thr and c != BASE]
    print('\nvariants with KD below trastuzumab reference max (%.2f nM): %d' % (thr, len(better)))
    for v, c in sorted(better)[:8]:
        print('   %-18s %.3f nM  (edit dist %d)' % (c, v, ed(c, BASE)))

# ---------- FACS ----------
print('\n' + '=' * 60)
cls = Counter()
dist = Counter()
by = defaultdict(Counter)
pos = {c: [Counter() for _ in range(10)] for c in ('high', 'medium', 'low')}
facs_seqs = set()
vls = set()
n = 0
with open(D + 'facs_unz/tz_her2_facs_3classes_524346.csv') as f:
    for r in csv.DictReader(f):
        n += 1
        c, s = r['binding_class'], r['cdrh3']
        cls[c] += 1
        d = int(r['edit_distance'])
        dist[d] += 1
        by[d][c] += 1
        for i, aa in enumerate(s[:10]):
            pos[c][i][aa] += 1
        facs_seqs.add(s)
        if n <= 50000:
            vls.add(r['vl_sequence'])
print('FACS rows:', n, 'unique cdrh3:', len(facs_seqs), 'unique VL in first 50k:', len(vls))
print('classes:', cls.most_common())
print('\nedit distance -> class share:')
print('  dist   n      high%  med%   low%')
for d in sorted(dist):
    t = dist[d]
    print('  %4d %7d   %5.1f %5.1f %5.1f' % (d, t, 100 * by[d]['high'] / t, 100 * by[d]['medium'] / t, 100 * by[d]['low'] / t))

print('\nper-position enrichment, high vs low (top residue by log-ratio):')
print('  pos  tz  high-favoured      low-favoured')
for i in range(10):
    h, l = pos['high'][i], pos['low'][i]
    th, tl = sum(h.values()), sum(l.values())
    sc = {}
    for aa in set(h) | set(l):
        ph = (h.get(aa, 0) + 1) / (th + 20)
        pl = (l.get(aa, 0) + 1) / (tl + 20)
        sc[aa] = ph / pl
    up = sorted(sc.items(), key=lambda x: -x[1])[:3]
    dn = sorted(sc.items(), key=lambda x: x[1])[:3]
    print('  %3d  %s   %-18s %s' % (i + 1, BASE[i],
          ' '.join('%s%.1fx' % (a, v) for a, v in up),
          ' '.join('%s%.2fx' % (a, v) for a, v in dn)))

# overlap
spr_seqs = {r['cdrh3'] for r in spr}
print('\nSPR CDRH3s also in FACS:', len(spr_seqs & facs_seqs))
kd_in = [(v, c) for v, c in kd if c in facs_seqs]
print('of those, with a KD:', len(kd_in))
