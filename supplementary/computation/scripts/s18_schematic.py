"""Fig. 1 of the article: model schematic (an aggregate of turnover-renewed transmission units seeding a clustered,
subcritical population).  X and Y are seeded separately, at rates kappa_X(t) and kappa_Y(t), and do not compete in the
seeded population; their clusters are drawn with the same pattern.  Drawing only; no model quantities are computed."""
import math
from _common import *
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
out = os.path.join(FIG, 'fig1.png')
plt.rcParams.update({'font.family': 'DejaVu Serif', 'mathtext.fontset': 'dejavuserif', 'font.size': 8})
CX, CY = '#1f77b4', '#ff7f0e'
Y0, Y1 = 0.55, 7.45                                     # vertical extent of the drawing (0.5 inch per unit)
fig = plt.figure(figsize=(6.6, (Y1 - Y0) / 2)); ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 13.2); ax.set_ylim(Y0, Y1); ax.axis('off')


def box(x, y, w, h, fc='white', ec='0.35', lw=0.9, r=0.12, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw, zorder=z))


def arrow(x0, y0, x1, y1, lw=0.9, c='0.25', ls='-', ms=7, z=3, cs='arc3,rad=0'):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=ms, lw=lw, color=c, ls=ls, zorder=z,
                                 connectionstyle=cs, shrinkA=0, shrinkB=0))


# ---------------- aggregate of transmission units
B = 2.75                                                # lower edge of the front unit
ax.text(2.95, 7.15, 'Aggregate of units (TRTA)', ha='center', va='center', fontsize=9, weight='bold')
for dx, dy in [(0.36, 0.30), (0.18, 0.15)]:
    box(0.35 + dx, B + dy, 4.9, 3.6, fc='0.97', ec='0.6', lw=0.7, z=1)
box(0.35, B, 4.9, 3.6, fc='white', ec='0.3', lw=1.0, z=2)
ax.text(0.55, B + 3.35, 'transmission unit (TRTU)', ha='left', va='center', fontsize=8, style='italic', zorder=3)
cw, ch, cyy = 0.78, 0.62, B + 1.5
xs = {'S': 0.75, 'E': 1.95, 'I': 3.15, 'R': 4.35}
for lab, x in xs.items():
    if lab in ('E', 'I'):
        box(x, cyy + ch / 2, cw, ch / 2, fc=CX if lab == 'I' else '#c6dbef', ec='0.3', lw=0.8, r=0.05, z=4)
        box(x, cyy, cw, ch / 2, fc=CY if lab == 'I' else '#fdd0a2', ec='0.3', lw=0.8, r=0.05, z=4)
        tc = 'white' if lab == 'I' else 'black'
        ax.text(x + cw / 2, cyy + 3 * ch / 4, lab + r'$_X$', ha='center', va='center', fontsize=7.5, color=tc, zorder=5)
        ax.text(x + cw / 2, cyy + ch / 4, lab + r'$_Y$', ha='center', va='center', fontsize=7.5, color=tc, zorder=5)
    else:
        box(x, cyy, cw, ch, fc='#e5f5e0' if lab == 'S' else '0.88', ec='0.3', lw=0.8, r=0.05, z=4)
        ax.text(x + cw / 2, cyy + ch / 2, lab, ha='center', va='center', fontsize=8.5, zorder=5)
for a, b in (('S', 'E'), ('E', 'I'), ('I', 'R')):
    arrow(xs[a] + cw + 0.04, cyy + ch / 2, xs[b] - 0.04, cyy + ch / 2, z=5)
ax.text(1.95, B + 2.7, 'X and Y compete for\nrenewed susceptibles', ha='left', va='center', fontsize=7, zorder=5)
arrow(1.9, B + 2.7, 1.3, cyy + ch + 0.05, lw=0.6, c='0.4', ms=5, z=5, cs='arc3,rad=0.25')
# departures from every compartment
ybar = B + 0.7
for lab, x in xs.items():
    ax.plot([x + cw / 2 + (0.18 if lab == 'S' else 0)] * 2, [cyy - 0.03, ybar], color='0.55', lw=0.7, zorder=4)
ax.plot([xs['S'] + cw / 2 + 0.18, xs['R'] + cw / 2], [ybar, ybar], color='0.55', lw=0.7, zorder=4)
arrow(xs['R'] + cw / 2, ybar, xs['R'] + cw / 2, B - 0.4, lw=0.9, c='0.55', ms=6, z=4)
ax.text(xs['R'] + cw / 2 - 0.1, B + 0.35, 'departures\n(turnover $m$)', ha='right', va='center', fontsize=7, color='0.4', zorder=5)
# susceptible arrivals
arrow(xs['S'] + cw / 2 - 0.15, B - 0.4, xs['S'] + cw / 2 - 0.15, cyy - 0.04, z=5, c='#2ca02c', lw=1.1)
ax.text(xs['S'] + cw / 2 - 0.02, B + 0.35, 'susceptible\narrivals', ha='left', va='center', fontsize=7, color='#2ca02c', zorder=5)
ax.text(2.8, B - 0.7, 'several units, seeded at different times', ha='center', va='center', fontsize=7, color='0.35')
# ---------------- seeded population: one block of clusters per strain, drawn with the same pattern
R = 0.47
cols = [8.3 + 1.28 * i for i in range(4)]
rows = {CX: (6.2, 5.05), CY: (3.7, 2.55)}              # upper and lower row of the block of each strain
ax.text((cols[0] + cols[-1]) / 2, 7.15, 'Seeded population (CSCP)', ha='center', va='center', fontsize=9, weight='bold')
pattern = {(0, 0): 4, (1, 0): 3, (1, 1): 2}              # (column, row in block): exposed members
links = [((0, 0), (1, 0)), ((0, 0), (1, 1))]             # between-cluster transmission


def centre(col, i, j):
    return cols[i], rows[col][j]


for col in (CX, CY):
    for i in range(4):
        for j in range(2):
            cx, cy = centre(col, i, j); ninf = pattern.get((i, j), 0)
            ax.add_patch(Circle((cx, cy), R, fc='white', ec='0.3' if ninf else '0.6', lw=0.9 if ninf else 0.7, zorder=2))
            pts = [(cx + 0.25 * math.cos(2 * math.pi * k / 5 + 0.3), cy + 0.25 * math.sin(2 * math.pi * k / 5 + 0.3)) for k in range(5)]
            for k, (px, py) in enumerate(pts):
                ax.add_patch(Circle((px, py), 0.07, fc=(col if k < ninf else 'white'), ec='0.3', lw=0.5, zorder=4))
            for k in range(ninf - 1):
                (x0, y0), (x1, y1) = pts[k], pts[k + 1]
                ax.plot([x0, x1], [y0, y1], color=col, lw=1.6, zorder=3, solid_capstyle='round')
    for a, b in links:
        (x0, y0), (x1, y1) = centre(col, *a), centre(col, *b)
        d = math.hypot(x1 - x0, y1 - y0); ux, uy = (x1 - x0) / d, (y1 - y0) / d
        arrow(x0 + ux * R, y0 + uy * R, x1 - ux * R, y1 - uy * R, lw=0.9, c=col, ls=(0, (3, 2)), ms=6, z=3)
# ---------------- seeding: one arrow per strain, ending at the first seeded cluster of its block
xa0, xa1 = 5.72, cols[0] - R - 0.16
xm = (xa0 + xa1) / 2
for col, lab in ((CX, r'$\kappa_X(t)=\alpha\, p_{I,X}(t)$'), (CY, r'$\kappa_Y(t)=\alpha\, p_{I,Y}(t)$')):
    y = rows[col][0]
    arrow(xa0, y, xa1, y, lw=1.8, c=col, ms=12, z=3)
    ax.text(xm, y - 0.42, lab, ha='center', va='center', fontsize=8)
ax.text(xm, (rows[CX][0] + rows[CY][0]) / 2 - 0.2, 'seeding', ha='center', va='center', fontsize=8)
# legend
ly = 1.45
ax.plot([7.75, 8.25], [ly, ly], color='0.3', lw=1.6, solid_capstyle='round')
ax.text(8.4, ly, 'within clusters: strong transmission', ha='left', va='center', fontsize=7)
arrow(7.75, ly - 0.45, 8.25, ly - 0.45, lw=0.9, c='0.3', ls=(0, (3, 2)), ms=6)
ax.text(8.4, ly - 0.45, 'between clusters: weak (subcritical)', ha='left', va='center', fontsize=7)
fig.savefig(out, dpi=300); plt.close(fig); print('wrote fig1.png')
