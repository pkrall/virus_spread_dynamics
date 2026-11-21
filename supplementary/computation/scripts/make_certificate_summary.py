"""Certificate summary table (box radius r, iterate power P, bound on the derivative of the P-step map over the box,
step at which the certified trajectory enters the box) for every certified trajectory, from summary.txt."""
import os, re, csv, math
up = lambda x: f'{math.ceil(float(x) * 1e4) / 1e4:.4f}'   # outward (upward) rounding of an upper bound
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')); C = os.path.join(ROOT, 'outputs', 'certificates')
txt = open(os.path.join(C, 'summary.txt')).read()
box = {int(s): (r, P, n) for s, r, P, n in re.findall(r'scenario (\d): fixed point .*?box radius ([\d.]+), P = (\d+), \|DF\^P\| <= ([\dE.+-]+)', txt)}
rows = []
for s, step in re.findall(r'Table 1 scenario (\d): box entered at step (\d+)', txt):
    r, P, n = box[int(s)]; rows.append([f'Table 2, scenario {s}', s, r, P, up(n), step])
for tab, ramp, step in re.findall(r'Table (S[12]), (\d+)-step ramp: box entered at step (\d+)', txt):
    s = 4 if tab == 'S1' else 3; r, P, n = box[s]
    rows.append([f'Table {tab}, {ramp}-step ramp', s, r, P, up(n), step])
with open(os.path.join(C, 'certificate_summary.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['certified trajectory', 'box of scenario', 'box radius r', 'iterate power P', 'bound on |DF^P| over box', 'step box entered']); w.writerows(rows)
print('wrote certificate_summary.csv')
