import json
import os
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from ..core.utils import log_message, get_output_path, ensure_directory_exists


def _compute_iqr_threshold(distances: np.ndarray) -> float:
    if distances.size == 0:
        return float("inf")
    q1 = np.percentile(distances, 25)
    q3 = np.percentile(distances, 75)
    iqr = q3 - q1
    return float(q3 + 1.5 * iqr)


def _load_gdf_any(gdf_path: str) -> gpd.GeoDataFrame:
    try:
        gdf = gpd.read_parquet(gdf_path)
    except Exception:
        gdf = gpd.read_file(gdf_path)
    required = {"id", "x", "y", "geometry"}
    missing = required.difference(gdf.columns)
    if missing:
        raise ValueError(f"gdf missing required columns: {sorted(missing)}")
    return gdf


def _infer_sample_names(files_config: Sequence[Mapping[str, str]]) -> Sequence[str]:
    unique_samples = sorted(list(set(
        f['name'].replace('_TRA', '').replace('_TRB', '') for f in files_config
    )))
    return unique_samples


def _load_env_gdf_map() -> Dict[str, str]:
    path = os.environ.get('SPATIALTCR_GDF_MAP')
    if path and os.path.exists(path):
        with open(path, 'r') as f:
            data = json.load(f)
        if isinstance(data, dict):
            return {str(k): str(v) for k, v in data.items()}
    return {}


def _compute_st_include(df: pd.DataFrame, mask_xy: np.ndarray, radius: float = 0.1) -> pd.DataFrame:
    if mask_xy.size == 0 or df.empty:
        df = df.copy()
        df['st_include'] = True
        return df
    from scipy.spatial import cKDTree
    tree = cKDTree(mask_xy)
    x_col = 'major_XPOS' if 'major_XPOS' in df.columns else 'XPOS'
    y_col = 'major_YPOS' if 'major_YPOS' in df.columns else 'YPOS'
    pts = np.column_stack([df[x_col].astype(float).values, df[y_col].astype(float).values])
    dists, _ = tree.query(pts, k=1, distance_upper_bound=radius)
    include = np.isfinite(dists)
    out = df.copy()
    out['st_include'] = include
    return out


def _spatial_join_map(points_df: pd.DataFrame, gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    pts = gpd.GeoDataFrame(points_df.copy(), geometry=[Point(xy) for xy in zip(points_df['major_XPOS'], points_df['major_YPOS'])], crs=gdf.crs)
    try:
        joined = gpd.sjoin(pts, gdf[['id', 'x', 'y', 'geometry']], how='left', predicate='within')
    except Exception:
        joined = gpd.sjoin(pts, gdf[['id', 'x', 'y', 'geometry']], how='left')
    out = points_df.copy()
    out['id'] = joined['id'].astype(str)
    out['x'] = joined['x']
    out['y'] = joined['y']
    out['geometry'] = joined['geometry_right'] if 'geometry_right' in joined.columns else joined['geometry']
    out['id'] = out['id'].fillna('NotMatched')
    return out


def _remap_not_matched(df: pd.DataFrame, gdf: gpd.GeoDataFrame, distance_method: str = "iqr", max_distance: Optional[float] = None) -> pd.DataFrame:
    from scipy.spatial import cKDTree

    id_series = df['id'].astype(str).fillna("").str.strip()
    not_matched_mask = id_series.str.lower().isin([v.lower() for v in ("Not Matched", "NotMatched", "NA", "", "not matched", "Not matched")])
    if not not_matched_mask.any():
        return df

    gdf_coords = np.column_stack([gdf['x'].astype(float).values, gdf['y'].astype(float).values])
    if len(gdf_coords) == 0:
        return df
    tree = cKDTree(gdf_coords)

    candidates = df.loc[not_matched_mask].copy()
    pts = np.column_stack([
        candidates['major_XPOS'].astype(float).values,
        candidates['major_YPOS'].astype(float).values,
    ])
    distances, nn_idx = tree.query(pts, k=1)

    if max_distance is not None:
        threshold = float(max_distance)
    else:
        m = (distance_method or 'iqr').lower()
        if m == 'iqr':
            threshold = _compute_iqr_threshold(distances)
        elif m == 'zscore':
            mu = float(np.mean(distances)); sigma = float(np.std(distances))
            threshold = mu + 3.0 * sigma
        elif m == 'none':
            threshold = float('inf')
        else:
            threshold = _compute_iqr_threshold(distances)

    within = distances <= threshold
    if within.any():
        matched_cells = gdf.iloc[nn_idx[within]]
        target_index = candidates.index[within]
        df.loc[target_index, 'id'] = matched_cells['id'].astype(str).values
        df.loc[target_index, 'x'] = matched_cells['x'].astype(float).values
        df.loc[target_index, 'y'] = matched_cells['y'].astype(float).values
        df.loc[target_index, 'geometry'] = matched_cells['geometry'].apply(lambda g: getattr(g, 'wkt', g)).values
        out_dist = pd.Series(distances, index=candidates.index, name='nn_distance')
        df['nn_distance'] = out_dist
    else:
        df['nn_distance'] = np.nan

    return df


def _load_major_second_source(base_dir: str, chain: str) -> Optional[pd.DataFrame]:
    path = os.path.join(base_dir, '06_filtered_results', f'major_second_{chain}.csv')
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
        if 'chain_type' not in df.columns:
            df['chain_type'] = chain
        if 'sample' not in df.columns:
            df['sample'] = None
        return df
    except Exception as e:
        log_message(f"Mapper: failed reading 06 file {os.path.basename(path)}: {e}")
        return None


def _prepare_mapping_frame(df: pd.DataFrame) -> pd.DataFrame:
    mapped = df.copy()
    # Prefer major coords; fall back to generic XPOS/YPOS if present
    if 'major_XPOS' not in mapped.columns and 'XPOS' in mapped.columns:
        mapped['major_XPOS'] = mapped['XPOS']
    if 'major_YPOS' not in mapped.columns and 'YPOS' in mapped.columns:
        mapped['major_YPOS'] = mapped['YPOS']
    for col in ['major_XPOS', 'major_YPOS']:
        if col not in mapped.columns:
            mapped[col] = np.nan
        mapped[col] = pd.to_numeric(mapped[col], errors='coerce')
    mapped = mapped.dropna(subset=['major_XPOS', 'major_YPOS'])
    if 'id' not in mapped.columns:
        mapped['id'] = 'NotMatched'
    if 'x' not in mapped.columns:
        mapped['x'] = np.nan
    if 'y' not in mapped.columns:
        mapped['y'] = np.nan
    if 'geometry' not in mapped.columns:
        mapped['geometry'] = None
    return mapped


def _write_cell_observations(mapped_frames: List[pd.DataFrame], obs_dir: str) -> None:
    if not mapped_frames:
        log_message("Obs builder: no mapped frames; skipping")
        return
    os.makedirs(obs_dir, exist_ok=True)
    all_df = pd.concat(mapped_frames, ignore_index=True)
    # Keep matched cells only
    all_df = all_df[all_df['id'].astype(str).str.lower() != 'notmatched']
    if all_df.empty:
        log_message("Obs builder: no matched cells to summarize")
        return

    def _top_two_for_chain(cell_df: pd.DataFrame, chain: str) -> Dict[str, Any]:
        records: List[Dict[str, Any]] = []
        for _, row in cell_df[cell_df['chain_type'] == chain].iterrows():
            records.append({
                'tcr_id': row.get('major_tcr_id'),
                'abundance': pd.to_numeric(row.get('major_abundance'), errors='coerce'),
                'v_call': row.get('major_v_call'),
                'j_call': row.get('major_j_call'),
                'cdr3aa': row.get('major_junction_aa')
            })
            records.append({
                'tcr_id': row.get('second_tcr_id'),
                'abundance': pd.to_numeric(row.get('second_abundance'), errors='coerce'),
                'v_call': row.get('second_v_call'),
                'j_call': row.get('second_j_call'),
                'cdr3aa': row.get('second_junction_aa')
            })
        # Remove rows with no TCR id or calls
        rec_df = pd.DataFrame(records)
        rec_df = rec_df.dropna(subset=['tcr_id', 'v_call', 'j_call', 'cdr3aa'], how='all')
        if rec_df.empty:
            return {}
        rec_df['abundance'] = rec_df['abundance'].fillna(0)
        rec_df = rec_df.groupby('tcr_id').agg({
            'abundance': 'max',
            'v_call': 'first',
            'j_call': 'first',
            'cdr3aa': 'first'
        }).reset_index()
        rec_df = rec_df.sort_values(by='abundance', ascending=False)
        out: Dict[str, Any] = {}
        for idx in range(2):
            if idx < len(rec_df):
                r = rec_df.iloc[idx]
                out[f'chain_{idx+1}_{chain}_abundance'] = r['abundance']
                out[f'chain_{idx+1}_{chain}_v_call'] = r['v_call']
                out[f'chain_{idx+1}_{chain}_j_call'] = r['j_call']
                out[f'chain_{idx+1}_{chain}_cdr3aa'] = r['cdr3aa']
                out[f'chain_{idx+1}_{chain}_tcrid'] = r['tcr_id']
            else:
                out[f'chain_{idx+1}_{chain}_abundance'] = None
                out[f'chain_{idx+1}_{chain}_v_call'] = None
                out[f'chain_{idx+1}_{chain}_j_call'] = None
                out[f'chain_{idx+1}_{chain}_cdr3aa'] = None
                out[f'chain_{idx+1}_{chain}_tcrid'] = None
        return out

    for sample, sample_df in all_df.groupby('sample'):
        if sample_df.empty:
            continue
        rows: List[Dict[str, Any]] = []
        for cell_id, cell_df in sample_df.groupby('id'):
            row: Dict[str, Any] = {
                'sample': sample,
                'cell_id': cell_id,
                'x': float(cell_df['x'].dropna().median()) if cell_df['x'].notna().any() else None,
                'y': float(cell_df['y'].dropna().median()) if cell_df['y'].notna().any() else None,
            }
            if 'geometry' in cell_df.columns:
                geom = cell_df['geometry'].dropna()
                row['geometry'] = geom.iloc[0] if not geom.empty else None
            row.update(_top_two_for_chain(cell_df, 'TRA'))
            row.update(_top_two_for_chain(cell_df, 'TRB'))
            rows.append(row)
        if not rows:
            log_message(f"Obs builder: no rows for sample {sample}")
            continue
        obs_df = pd.DataFrame(rows)
        out_path = os.path.join(obs_dir, f'cell_obs_{sample}.parquet')
        obs_df.to_parquet(out_path, index=False)
        log_message(f"Obs builder: wrote {len(obs_df)} cells -> {os.path.basename(out_path)}")


def run_cid_mapper(
    base_dir: str,
    files_config: Sequence[Mapping[str, str]],
    gdf_map: Optional[Mapping[str, str]] = None,
    distance_method: str = "iqr",
    max_distance: Optional[float] = None
) -> None:
    log_message("### Step 5: Mapping CIDs to Cells (Mapper) ###")

    samples = _infer_sample_names(files_config)
    env_map = _load_env_gdf_map()
    gdf_map_final: Dict[str, str] = {}
    if gdf_map:
        gdf_map_final.update({str(k): str(v) for k, v in gdf_map.items()})
    gdf_map_final.update(env_map)

    major_sources = {
        chain: _load_major_second_source(base_dir, chain) for chain in ['TRA', 'TRB']
    }

    def _load_st_mask_for_sample(sample: str) -> np.ndarray:
        p = os.environ.get(f"SPATIALTCR_ST_MASK_{sample}")
        if not p or not os.path.exists(p):
            return np.empty((0, 2), dtype=float)
        if p.endswith(".parquet"):
            m = pd.read_parquet(p)
        else:
            m = pd.read_csv(p)
        if not {'x', 'y'}.issubset(m.columns):
            return np.empty((0, 2), dtype=float)
        return np.column_stack([m['x'].astype(float).values, m['y'].astype(float).values])

    mapped_out_dir = os.path.join(base_dir, '07_mapped_results')
    obs_out_dir = os.path.join(base_dir, '08_cell_obs')
    ensure_directory_exists(mapped_out_dir)
    ensure_directory_exists(obs_out_dir)

    mapped_frames: List[pd.DataFrame] = []

    for sample in samples:
        gdf_path = gdf_map_final.get(sample)
        if not gdf_path or not os.path.exists(gdf_path):
            log_message(f"Mapper: skipping sample '{sample}' - missing gdf path. Provide via gdf_map or SPATIALTCR_GDF_MAP.")
            continue
        try:
            gdf = _load_gdf_any(gdf_path)
        except Exception as e:
            log_message(f"Mapper: failed to load gdf for {sample}: {e}")
            continue

        st_mask_xy = _load_st_mask_for_sample(sample)

        for chain in ['TRA', 'TRB']:
            source_df = major_sources.get(chain)
            if source_df is not None:
                df = source_df.copy()
                if df['sample'].notna().any():
                    df = df[df['sample'] == sample]
            else:
                in_path = get_output_path(base_dir, 'final_tcr', f'final_tcr_{sample}_{chain}.parquet')
                if not os.path.exists(in_path):
                    log_message(f"Mapper: input missing for {sample}_{chain}, skipping")
                    continue
                try:
                    df = pd.read_parquet(in_path)
                    df['chain_type'] = chain
                    df['sample'] = sample
                except Exception as e:
                    log_message(f"Mapper: failed reading {os.path.basename(in_path)}: {e}")
                    continue

            if df.empty:
                log_message(f"Mapper: empty TCR table for {sample}_{chain}, skipping")
                continue

            map_df = _prepare_mapping_frame(df)
            if map_df.empty:
                log_message(f"Mapper: no valid coordinates for {sample}_{chain}")
                continue

            map_df['chain_type'] = chain
            map_df['sample'] = sample

            map_df = _compute_st_include(map_df, st_mask_xy, radius=0.1)
            map_df = map_df[map_df['st_include'] == True].copy()
            if map_df.empty:
                log_message(f"Mapper: {sample}_{chain} no points within ST mask; skipping")
                continue

            mapped_df = _spatial_join_map(map_df, gdf)
            mapped_df = _remap_not_matched(mapped_df, gdf, distance_method=distance_method, max_distance=max_distance)

            out_path = os.path.join(mapped_out_dir, f'cell_mapping_{sample}_{chain}.csv')
            if 'geometry' in mapped_df.columns:
                mapped_df = mapped_df.copy()
                mapped_df['geometry'] = mapped_df['geometry'].apply(lambda g: getattr(g, 'wkt', g))
            mapped_df.to_csv(out_path, index=False)
            mapped_frames.append(mapped_df.assign(chain_type=chain, sample=sample))
            log_message(f"Mapper: wrote {len(mapped_df)} rows -> {os.path.basename(out_path)}")

    _write_cell_observations(mapped_frames, obs_out_dir)
    log_message("### Step 5: Completed ###")


