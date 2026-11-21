import os, sys, csv
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
TAB = os.path.join(ROOT, 'outputs', 'tables'); FIG = os.path.join(ROOT, 'outputs', 'figures')
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 8, 'axes.titlesize': 8.5, 'legend.fontsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.spines.top': False, 'axes.spines.right': False})
SCEN = {1: (0.05, 8, 0.025, 3, 4), 2: (0.05, 20, 0.035, 3, 4), 3: (0.12, 10, 0.035, 3, 4), 4: (0.12, 20, 0.035, 3, 4), 5: (0.12, 20, 0.020, 4, 5)}

def write_csv(name, header, rows):
    with open(os.path.join(TAB, name), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print('wrote', name)

def savefig(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(FIG, name), dpi=300, bbox_inches='tight'); plt.close(fig)
    print('wrote', name)
