"""Ordering of Alpha disappearance and Delta establishment in EU countries.

Input: ECDC "Data on SARS-CoV-2 variants in the EU/EEA" (country-week records), comma- or
semicolon-separated (e.g. after saving from a spreadsheet program). Only GISAID-derived records of the 27 EU member states are used.

Definitions (main specification; see PARAMS):
  shares          variant detections / sequences with known variant (number_sequenced_known_variant)
  Alpha           B.1.1.7 and B.1.1.7+E484K;  Delta: B.1.617.2 and AY.4.2
  qualifying week at least MIN_SEQ sequences with known variant
  Alpha established  first of K consecutive qualifying weeks with Alpha share >= ALPHA_EST
  Delta established  first of K consecutive qualifying weeks (from Alpha establishment on) with Delta share >= DELTA_EST
  endpoint        first of two consecutive calendar weeks, both qualifying, with Alpha share < ALPHA_END
  violation       a qualifying week strictly between Alpha and Delta establishment with Alpha share < ALPHA_END

Output: table_s7.csv (Table S7 of Online Resource 1) and sensitivity.csv.
Usage: python alpha_delta_ordering.py strain_data.csv outdir
"""
import sys, itertools
import numpy as np, pandas as pd

EU = ['Austria', 'Belgium', 'Bulgaria', 'Croatia', 'Cyprus', 'Czechia', 'Denmark', 'Estonia', 'Finland', 'France',
      'Germany', 'Greece', 'Hungary', 'Ireland', 'Italy', 'Latvia', 'Lithuania', 'Luxembourg', 'Malta', 'Netherlands',
      'Poland', 'Portugal', 'Romania', 'Slovakia', 'Slovenia', 'Spain', 'Sweden']
ALPHA = ['B.1.1.7', 'B.1.1.7+E484K']
DELTA = ['B.1.617.2', 'AY.4.2']
PARAMS = dict(MIN_SEQ=50, ALPHA_EST=0.50, DELTA_EST=0.10, ALPHA_END=0.01, K=2, FIRST='2020-40', LAST='2021-52')


def weekly(path):
    with open(path, encoding='utf-8-sig') as f:
        sep = ';' if f.readline().count(';') > 3 else ','
    d = pd.read_csv(path, sep=sep, encoding='utf-8-sig', low_memory=False)
    d = d[(d.source == 'GISAID') & d.country.isin(EU)]
    key = ['country', 'year_week']
    base = d.groupby(key)['number_sequenced_known_variant'].first().rename('n')
    a = d[d.variant.isin(ALPHA)].groupby(key).number_detections_variant.sum().rename('alpha')
    de = d[d.variant.isin(DELTA)].groupby(key).number_detections_variant.sum().rename('delta')
    w = pd.concat([base, a, de], axis=1).fillna({'alpha': 0, 'delta': 0}).reset_index()
    w = w.sort_values(key).reset_index(drop=True)
    w['pa'] = np.where(w.n > 0, w.alpha / w.n.clip(lower=1), np.nan)
    w['pd'] = np.where(w.n > 0, w.delta / w.n.clip(lower=1), np.nan)
    return w


def first_run(mask, k, start=0):
    run = 0
    for i in range(start, len(mask)):
        run = run + 1 if mask[i] else 0
        if run >= k:
            return i - k + 1
    return None


def analyse(w, MIN_SEQ, ALPHA_EST, DELTA_EST, ALPHA_END, K, FIRST, LAST):
    rows = []
    for c in EU:
        g = w[(w.country == c) & w.year_week.between(FIRST, LAST)].reset_index(drop=True)
        q = (g.n >= MIN_SEQ).values
        pa = g.pa.fillna(0).values; pdl = g.pd.fillna(0).values
        tA = first_run(q & (pa >= ALPHA_EST), K)
        tD = first_run(q & (pdl >= DELTA_EST), K, tA) if tA is not None else None
        tE = None
        if tA is not None:
            low = q & (pa < ALPHA_END)
            tE = next((i for i in range(tA + 1, len(g) - 1) if low[i] and low[i + 1]), None)
        between = [i for i in range(tA + 1, tD) if q[i]] if (tA is not None and tD is not None) else []
        yw = lambda i: g.year_week[i] if i is not None else ''
        rows.append(dict(country=c, qualifying_weeks=int(q.sum()), alpha_established=yw(tA), delta_established=yw(tD),
                         qualifying_weeks_between=len(between),
                         min_alpha_share_between=(round(100 * min(pa[i] for i in between), 1) if between else np.nan),
                         violations=sum(pa[i] < ALPHA_END for i in between), endpoint=yw(tE),
                         weeks_delta_to_endpoint=(tE - tD if (tE is not None and tD is not None) else np.nan)))
    return pd.DataFrame(rows)


def summary(r):
    est = r.alpha_established != ''
    return dict(alpha_established=int(est.sum()), endpoint=int((r.endpoint != '').sum()),
                endpoint_after_delta=int((r.weeks_delta_to_endpoint >= 0).sum()),
                countries_with_violation=int((r.violations > 0).sum()))


if __name__ == '__main__':
    path, out = sys.argv[1], sys.argv[2]
    w = weekly(path)
    r = analyse(w, **PARAMS)
    r.to_csv(f'{out}/table_s7.csv', index=False)
    print('main specification', PARAMS, summary(r))
    sens = []
    for n, a, b in itertools.product([20, 30, 50, 100], [0.25, 0.5], [0.05, 0.1, 0.5]):
        p = dict(PARAMS, MIN_SEQ=n, ALPHA_EST=a, DELTA_EST=b)
        rr = analyse(w, **p)
        sens.append(dict(min_seq=n, alpha_est=a, delta_est=b, **summary(rr),
                         excluded=', '.join(rr.country[rr.alpha_established == '']),
                         violation_countries=', '.join(rr.country[rr.violations > 0])))
    pd.DataFrame(sens).to_csv(f'{out}/sensitivity.csv', index=False)
    print(pd.DataFrame(sens).to_string())
