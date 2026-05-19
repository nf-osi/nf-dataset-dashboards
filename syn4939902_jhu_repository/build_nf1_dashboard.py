#!/usr/bin/env python3
"""
Build script for the JHU NF1 Biospecimen Repository Dashboard.

Queries the three curated Synapse datasets, processes specimen/patient
metadata, and produces a fully self-contained nf1_dashboard.html with
Chart.js and Plotly inlined (no external dependencies).

Usage:
    python build_nf1_dashboard.py

Requirements:
    pip install synapseclient pandas
"""

import ast
import json
import urllib.request
from pathlib import Path

import pandas as pd
import synapseclient

# ── Configuration ──────────────────────────────────────────────────────────────

DATASETS = {
    'syn74376449': 'WES',
    'syn74376463': 'RNA-seq',
    'syn74376440': 'WGS',
}

# Pinned to current stable versions (update if datasets are re-versioned)
DS_VERSIONS = {
    'syn74376449': 'syn74376449.2',
    'syn74376463': 'syn74376463.2',
    'syn74376440': 'syn74376440.2',
}

DS_COLORS = {
    'WES':     '#1a5fd4',   # cobalt blue
    'RNA-seq': '#b45309',   # amber
    'WGS':     '#047857',   # emerald
}

TUMOR_COLORS = {
    'Plexiform Neurofibroma':                 '#4e79a7',
    'Malignant Peripheral Nerve Sheath Tumor': '#e15759',
    'Not Applicable':                          '#6b7280',
    'Unknown':                                 '#4b5563',
    'Cutaneous Neurofibroma':                  '#76b7b2',
    'Diffuse Infiltrating Neurofibroma':       '#59a14f',
    'Neurofibroma':                            '#f28e2b',
    'Atypical Neurofibroma':                   '#c084fc',
    'Neurofibroma with Degenerative Atypia':   '#9c755f',
    'Nodular Neurofibroma':                    '#edc948',
    'Localized Neurofibroma':                  '#b07aa1',
    'Massive Soft Tissue Neurofibroma':        '#fabfd2',
}

# Prefix symbols encoding tumor type in the sunburst
TUMOR_SYMBOLS = {
    'Plexiform Neurofibroma':                 '◎',
    'Malignant Peripheral Nerve Sheath Tumor': '▲',
    'Cutaneous Neurofibroma':                  '○',
    'Diffuse Infiltrating Neurofibroma':       '≋',
    'Atypical Neurofibroma':                   '△',
    'Nodular Neurofibroma':                    '●',
    'Neurofibroma':                            '◦',
    'Localized Neurofibroma':                  '□',
    'Not Applicable':                          '·',
    'Unknown':                                 '·',
}

CDN_PLOTLY  = 'https://cdn.plot.ly/plotly-2.32.0.min.js'
CDN_CHARTJS = 'https://cdn.jsdelivr.net/npm/chart.js@4.4.3/dist/chart.umd.min.js'

HERE      = Path(__file__).parent
TEMPLATE  = HERE / 'nf1_dashboard_template.html'
OUT_HTML  = HERE / 'nf1_dashboard.html'
LIB_CACHE = HERE / '.lib_cache'


# ── Synapse ────────────────────────────────────────────────────────────────────

def query_datasets(syn):
    frames = []
    for ds_id, assay in DATASETS.items():
        print(f'  Querying {ds_id} ({assay})…')
        df = syn.tableQuery(f'SELECT * FROM {ds_id}').asDataFrame()
        df['assay_norm'] = assay
        df['dataset_id'] = ds_id
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


# ── Data helpers ───────────────────────────────────────────────────────────────

def parse_list_col(val):
    """Parse Synapse STRING_LIST values (stored as Python list literals)."""
    if pd.isna(val):
        return []
    if isinstance(val, list):
        return val
    try:
        result = ast.literal_eval(str(val))
        return result if isinstance(result, list) else [str(val)]
    except Exception:
        return [str(val)]


def explode_tumor_types(df):
    """Return a long DataFrame with one row per (specimenID, tumorType)."""
    rows = []
    for _, row in df.iterrows():
        types = parse_list_col(row.get('tumorType')) or ['Unknown']
        for t in types:
            rows.append({'specimenID': row['specimenID'],
                         'assay_norm': row['assay_norm'],
                         'tumorType':  t})
    return pd.DataFrame(rows)


# ── Data processing ────────────────────────────────────────────────────────────

def build_data(df):
    specs_df = df.drop_duplicates('specimenID').copy()
    pats_df  = df.drop_duplicates('individualID').copy()

    spec_counts    = specs_df.groupby('individualID').size()
    multi_pat_ids  = set(spec_counts[spec_counts > 1].index)

    # --- Summary ---
    summary = {
        'total_files':             int(len(df)),
        'total_patients':          int(pats_df['individualID'].nunique()),
        'total_specimens':         int(specs_df['specimenID'].nunique()),
        'multi_specimen_patients': int(len(multi_pat_ids)),
    }

    # --- Per-dataset counts ---
    files_by_ds = {k: int(v) for k, v in df.groupby('assay_norm').size().items()}
    specs_by_ds = {k: int(v) for k, v in
                   specs_df.groupby('assay_norm')['specimenID'].nunique().items()}
    pats_by_ds  = {k: int(v) for k, v in
                   (df.drop_duplicates(['individualID', 'assay_norm'])
                      .groupby('assay_norm')['individualID'].nunique()).items()}

    # --- Sex distribution (patient level) ---
    sex_patients = {k: int(v) for k, v in
                    pats_df['sex'].fillna('Unknown').value_counts().items()
                    if k != 'Unknown'}

    # --- Age at diagnosis ---
    ages   = pd.to_numeric(pats_df['age'], errors='coerce').dropna()
    bins   = [0, 10, 20, 30, 40, 50, 60, 70, 200]
    labels = ['0-10', '11-20', '21-30', '31-40', '41-50', '51-60', '61-70', '71+']
    age_series = pd.Series(pd.cut(ages, bins=bins, labels=labels, right=True)).value_counts().sort_index()
    age_dist = {str(k): int(v) for k, v in age_series.items() if v > 0}

    # --- Tumor types (specimen level, exploded) ---
    tt_df   = explode_tumor_types(specs_df)
    tt_counts = {k: int(v) for k, v in tt_df['tumorType'].value_counts().items()}
    tumor_order = list(tt_counts.keys())

    # --- Tumor × dataset ---
    all_assays = list(DATASETS.values())
    tumor_modality = {assay: {t: 0 for t in tumor_order} for assay in all_assays}
    for assay in all_assays:
        sub = specs_df[specs_df['assay_norm'] == assay]
        for t, cnt in explode_tumor_types(sub)['tumorType'].value_counts().items():
            if t in tumor_modality[assay]:
                tumor_modality[assay][t] = int(cnt)

    # --- Tissue counts ---
    tissue_counts = {k: int(v) for k, v in
                     specs_df['tissue'].fillna('Unknown').value_counts().items()}

    # --- Sample type (tissue) × dataset ---
    tissue_order = list(tissue_counts.keys())
    tissue_modality = {assay: {t: 0 for t in tissue_order} for assay in all_assays}
    for assay in all_assays:
        sub = specs_df[specs_df['assay_norm'] == assay]
        for t, cnt in sub['tissue'].fillna('Unknown').value_counts().items():
            if t in tissue_modality[assay]:
                tissue_modality[assay][t] = int(cnt)

    # --- Multi-specimen patient breakdown (replaces sunburst) ---
    multi_pat_rows = []
    multi = specs_df[specs_df['individualID'].isin(multi_pat_ids)]
    for pat, grp in multi.groupby('individualID'):
        counts = grp.groupby('assay_norm').size().to_dict()
        multi_pat_rows.append({
            'patient': pat,
            'total':   int(grp.shape[0]),
            'WES':     int(counts.get('WES', 0)),
            'RNA-seq': int(counts.get('RNA-seq', 0)),
            'WGS':     int(counts.get('WGS', 0)),
        })
    multi_pat_rows.sort(key=lambda x: -x['total'])

    return {
        'summary':        summary,
        'files_by_ds':    files_by_ds,
        'specs_by_ds':    specs_by_ds,
        'pats_by_ds':     pats_by_ds,
        'sex_patients':   sex_patients,
        'age_dist':       age_dist,
        'tumor_types':    tt_counts,
        'tumor_order':    tumor_order,
        'tumor_modality': tumor_modality,
        'tissue_counts':   tissue_counts,
        'tissue_order':    tissue_order,
        'tissue_modality': tissue_modality,
        'tumor_colors':   TUMOR_COLORS,
        'ds_colors':      DS_COLORS,
        'multi_pat_data': multi_pat_rows,
        'sunburst':       build_sunburst(specs_df, multi_pat_ids),
    }


def build_sunburst(specs_df, multi_pat_ids):
    """Plotly sunburst: ring 1 = patients, ring 2 = specimens.
    Color encodes assay; symbol prefix encodes tumor type (blood overrides)."""
    ids     = ['root']
    labels  = ['All Patients']
    parents = ['']
    values  = [0]
    colors  = ['#f0ece4']   # root node: warm light

    multi = (specs_df[specs_df['individualID'].isin(multi_pat_ids)]
                     .sort_values(['individualID', 'specimenID']))

    for pat, grp in multi.groupby('individualID'):
        pat_node = f'pat_{pat}'
        ids.append(pat_node)
        labels.append(pat)
        parents.append('root')
        values.append(0)
        colors.append('#dde4f0')  # patient node: light blue-gray

        for i, (_, row) in enumerate(grp.iterrows()):
            spec_id = row['specimenID']
            assay   = row['assay_norm']
            tumors  = parse_list_col(row.get('tumorType'))
            tissue  = str(row.get('tissue', '')).lower()

            if tissue == 'blood':
                sym = '◆'
            else:
                primary = (tumors[0] if tumors else 'Unknown')
                sym = TUMOR_SYMBOLS.get(primary, '·')

            ids.append(f'spec_{pat}_{i}')
            labels.append(f'{sym} {spec_id}')
            parents.append(pat_node)
            values.append(1)
            colors.append(DS_COLORS.get(assay, '#94a3b8'))

    return {
        'ids': ids, 'labels': labels, 'parents': parents,
        'values': values, 'colors': colors, 'hovertext': [],
    }


# ── HTML assembly ──────────────────────────────────────────────────────────────

def fetch_or_cache(url, filename):
    LIB_CACHE.mkdir(exist_ok=True)
    path = LIB_CACHE / filename
    if path.exists():
        print(f'  Using cached {filename}')
    else:
        print(f'  Downloading {url}…')
        with urllib.request.urlopen(url) as resp:
            path.write_bytes(resp.read())
    return path.read_text(encoding='utf-8')


def assemble_html(data):
    html = TEMPLATE.read_text(encoding='utf-8')

    # Inject data
    html = html.replace('{{DATA}}', json.dumps(data, ensure_ascii=False))

    # Inline JS libraries (cached after first download)
    plotly_js  = fetch_or_cache(CDN_PLOTLY,  'plotly-2.32.0.min.js')
    chartjs_js = fetch_or_cache(CDN_CHARTJS, 'chart.umd.min.js')

    html = html.replace(
        f'<script src="{CDN_PLOTLY}"></script>',
        f'<script>{plotly_js}</script>',
    )
    html = html.replace(
        f'<script src="{CDN_CHARTJS}"></script>',
        f'<script>{chartjs_js}</script>',
    )
    return html


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print('Logging into Synapse…')
    syn = synapseclient.Synapse()
    syn.login(silent=True)

    print('Querying datasets…')
    df = query_datasets(syn)
    print(f'  {len(df):,} total files')

    print('Processing data…')
    data = build_data(df)
    s = data['summary']
    print(f"  {s['total_patients']} patients · "
          f"{s['total_specimens']} specimens · "
          f"{s['multi_specimen_patients']} multi-specimen")

    print('Assembling HTML…')
    html = assemble_html(data)
    OUT_HTML.write_text(html, encoding='utf-8')
    print(f'Done → {OUT_HTML}  ({len(html) // 1024} KB)')


if __name__ == '__main__':
    main()
