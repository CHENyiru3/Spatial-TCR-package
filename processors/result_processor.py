import glob
import os
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from ..core.utils import log_message, get_output_path, validate_file_exists
from ..core.config import FILE_PATTERNS, VALIDATION_THRESHOLDS


def parse_tcr_chain(chain_str: str) -> Optional[Dict[str, str]]:
    """
    Parse a single TCR chain string from TRUST4 report.
    
    Args:
        chain_str: Chain string from TRUST4 output
        
    Returns:
        Parsed chain data or None if invalid
    """
    if not isinstance(chain_str, str) or chain_str == '*':
        return None
    
    parts = chain_str.split(',')
    # Expected format: V-gene,*,J-gene,C-gene,CDR3-nt,CDR3-aa,abundance,contig_id,identity,flags
    if len(parts) < 10 or parts[5] in ['out_of_frame', '*']:
        return None  # Non-productive or malformed chain
    
    return {
        'v_call': parts[0],
        'j_call': parts[2], 
        'c_call': parts[3],
        'junction': parts[4],
        'junction_aa': parts[5],
        'abundance': float(parts[6]),
        'contig_id': parts[7],
    }


def finalize_cid_mode_results(
    base_dir: str, 
    files_config: Sequence[Mapping[str, str]]
) -> None:

    log_message("### Step 4: Finalizing TRUST4 Results (CID Mode) ###")
    
    unique_samples = set(
        f['name'].replace('_TRA', '').replace('_TRB', '') 
        for f in files_config
    )
    
    for sample_prefix in unique_samples:
        for chain_type in ['TRA', 'TRB']:
            _process_sample_chain_cid_mode(base_dir, sample_prefix, chain_type)
    
    log_message("### Step 4: Completed ###")


def finalize_cellular_mode_results(
    base_dir: str,
    files_config: Sequence[Mapping[str, str]], 
    cell_data: pd.DataFrame
) -> None:
    """
    Process TRUST4 results for cellular mode (cell-level analysis).
    
    Args:
        base_dir: Base output directory
        files_config: List of file configurations
        cell_data: Cell coordinate/geometry data
    """
    log_message("### Step 4: Finalizing TRUST4 Results (Cellular Mode) ###")
    
    unique_samples = sorted(list(set(
        f['name'].replace('_TRA', '').replace('_TRB', '') 
        for f in files_config
    )))
    
    for sample_prefix in unique_samples:
        _process_sample_cellular_mode(base_dir, sample_prefix, cell_data)
    
    log_message("### Step 4: Completed ###")


def _process_sample_chain_cid_mode(
    base_dir: str, 
    sample_prefix: str, 
    chain_type: str
) -> None:
    """Process a single sample and chain type for CID mode."""
    log_message(f"Processing: {sample_prefix}_{chain_type}")
    
    report_path = get_output_path(base_dir, 'trust4_raw', f'{sample_prefix}_{chain_type}_barcode_report.tsv')
    coords_path = get_output_path(base_dir, 'processed', f'{sample_prefix}_{chain_type}.coordinates.tsv')
    
    # Check if required files exist
    if not os.path.exists(report_path):
        log_message(f"WARNING: Report file not found, skipping: {os.path.basename(report_path)}")
        return
    
    if not os.path.exists(coords_path):
        log_message(f"WARNING: Coordinates file not found, skipping: {os.path.basename(coords_path)}")
        return
    

    # Load TRUST4 report
    report_df = pd.read_csv(report_path, sep='\t')
    report_df.rename(columns={'#barcode': 'CID'}, inplace=True)
    
    # Parse TCR chains
    all_tcrs = []
    primary_col = 'chain2' if chain_type == 'TRA' else 'chain1'
    secondary_col = 'secondary_chain2' if chain_type == 'TRA' else 'secondary_chain1'
    
    for _, row in report_df.iterrows():
        cid = row['CID']
        
        # Parse primary chain
        primary_tcr = parse_tcr_chain(row[primary_col])
        if primary_tcr:
            primary_tcr['CID'] = cid
            primary_tcr['is_primary'] = True
            all_tcrs.append(primary_tcr)
        
        # Parse secondary chain
        secondary_tcr = parse_tcr_chain(row[secondary_col])
        if secondary_tcr:
            secondary_tcr['CID'] = cid
            secondary_tcr['is_primary'] = False
            all_tcrs.append(secondary_tcr)
    
    if not all_tcrs:
        log_message(f"No productive TCRs found in {os.path.basename(report_path)}")
        return
    
    tcrs_df = pd.DataFrame(all_tcrs)
    
    # Load coordinates and merge
    coords_df = pd.read_csv(coords_path, sep='\t')
    coords_df.drop_duplicates(subset=['CID'], inplace=True)
    
    final_df = pd.merge(tcrs_df, coords_df, on='CID', how='left')
    
    # Generate clonotype IDs
    final_df.sort_values(by=['v_call', 'j_call', 'junction_aa'], inplace=True)
    final_df['tcr_id'] = final_df.groupby(['v_call', 'j_call', 'junction_aa']).ngroup()
    
    # Save results
    output_path = get_output_path(base_dir, 'final_tcr', f'final_tcr_{sample_prefix}_{chain_type}.parquet')
    final_df.to_parquet(output_path, index=False)
    
    log_message(f"Processed {len(tcrs_df)} TCRs for {len(final_df['CID'].unique())} CIDs. "
                f"Saved to {os.path.basename(output_path)}")


def _process_sample_cellular_mode(
    base_dir: str,
    sample_prefix: str, 
    cell_data: pd.DataFrame
) -> None:
    """Process a single sample for cellular mode (both TRA and TRB)."""
    log_message(f"Processing cellular results for: {sample_prefix}")
    
    all_cellular_tcrs = []
    
    for chain_type in ['TRA', 'TRB']:
        # Check for barcode mapping files
        barcode_map_path = get_output_path(base_dir, 'processed', f'{sample_prefix}_{chain_type}.cell_barcode_map.tsv')
        report_path = get_output_path(base_dir, 'trust4_raw', f'{sample_prefix}_{chain_type}_barcode_report.tsv')
        
        if not os.path.exists(barcode_map_path) or not os.path.exists(report_path):
            log_message(f"WARNING: Missing files for {sample_prefix}_{chain_type}, skipping")
            continue
    
        # Load barcode mapping and TRUST4 report
        barcode_map_df = pd.read_csv(barcode_map_path, sep='\t')
        report_df = pd.read_csv(report_path, sep='\t')
        report_df.rename(columns={'#barcode': 'pseudo_barcode'}, inplace=True)
        
        # Debug barcode intersection
        map_barcodes = set(barcode_map_df['pseudo_barcode'].unique())
        report_barcodes = set(report_df['pseudo_barcode'].unique())
        intersection = map_barcodes & report_barcodes
        
        log_message(f"[DEBUG] {sample_prefix}_{chain_type}: "
                    f"barcode_map.tsv count: {len(map_barcodes)}, "
                    f"barcode_report.tsv count: {len(report_barcodes)}")
        log_message(f"[DEBUG] {sample_prefix}_{chain_type}: intersection count: {len(intersection)}")
        
        if len(intersection) == 0:
            log_message(f"[DEBUG] {sample_prefix}_{chain_type}: No overlapping barcodes!")
            continue
        
        # Merge data
        merged_df = pd.merge(report_df, barcode_map_df, on='pseudo_barcode', how='inner')
        
        if merged_df.empty:
            continue
        
        # Parse TCR chains
        primary_col = 'chain2' if chain_type == 'TRA' else 'chain1'
        secondary_col = 'secondary_chain2' if chain_type == 'TRA' else 'secondary_chain1'
        
        for _, row in merged_df.iterrows():
            cell_id = row['cell_id']
            
            # Parse chains
            primary_tcr = parse_tcr_chain(row[primary_col])
            if primary_tcr:
                primary_tcr.update({
                    'cell_id': cell_id,
                    'chain_type': chain_type,
                    'is_primary': True
                })
                all_cellular_tcrs.append(primary_tcr)
            
            secondary_tcr = parse_tcr_chain(row[secondary_col])
            if secondary_tcr:
                secondary_tcr.update({
                    'cell_id': cell_id,
                    'chain_type': chain_type,
                    'is_primary': False
                })
                all_cellular_tcrs.append(secondary_tcr)
        
    
    # Combine all TCRs and add cell coordinates
    final_df = pd.DataFrame(all_cellular_tcrs)
    
    # Merge with cell coordinate data
    cell_coords_df = cell_data[['id']].copy()
    cell_coords_df['XPOS'] = cell_data.geometry.centroid.x
    cell_coords_df['YPOS'] = cell_data.geometry.centroid.y
    cell_coords_df['cell_id'] = cell_coords_df['id'].astype(str)
    
    final_df = pd.merge(final_df, cell_coords_df, on='cell_id', how='left')
    
    # Generate clonotype IDs
    final_df.sort_values(by=['v_call', 'j_call', 'junction_aa'], inplace=True)
    final_df['tcr_id'] = final_df.groupby(['v_call', 'j_call', 'junction_aa']).ngroup()
    
    # Save results
    output_path = get_output_path(base_dir, 'final_tcr', f'final_cellular_tcr_{sample_prefix}.parquet')
    final_df.to_parquet(output_path, index=False)
    
    log_message(f"Processed {len(final_df)} TCRs for {len(final_df['cell_id'].unique())} cells. "
               f"Saved to {os.path.basename(output_path)}")


def validate_results(base_dir: str, mode: str) -> bool:
    """
    Validate pipeline results.
    
    Args:
        base_dir: Base output directory
        mode: Pipeline mode ('cid' or 'cellular')
        
    Returns:
        True if validation passes
    """
    log_message(f"Validating {mode} mode results...")
    

    final_tcr_dir = get_output_path(base_dir, 'final_tcr', '')
    
    if not os.path.exists(final_tcr_dir):
        log_message("ERROR: Final TCR directory not found")
        return False
    
    parquet_files = [f for f in os.listdir(final_tcr_dir) if f.endswith('.parquet')]
    
    if not parquet_files:
        log_message("ERROR: No final TCR files found")
        return False
    
    log_message(f"Found {len(parquet_files)} final TCR files")
    
    # Validate each file
    total_tcrs = 0
    for parquet_file in parquet_files:
        file_path = os.path.join(final_tcr_dir, parquet_file)
        df = pd.read_parquet(file_path)
        total_tcrs += len(df)
        log_message(f"  {parquet_file}: {len(df)} TCRs")
    
    if total_tcrs < VALIDATION_THRESHOLDS['min_reads_per_sample']:
        log_message(f"WARNING: Low total TCR count: {total_tcrs}")
    
    log_message(f"Validation passed. Total TCRs: {total_tcrs}")
    return True
    

# =============================
# Format organizer (06 outputs)
# =============================

def _merge_major_second(df: pd.DataFrame) -> pd.DataFrame:
    """Organize per-chain TCRs into major/second format with exact column order."""
    # Normalize key column name: prefer 'CID', fallback to 'cell_id'
    if 'CID' not in df.columns:
        if 'cell_id' in df.columns:
            df = df.copy()
            df['CID'] = df['cell_id']
        else:
            raise KeyError("Input DataFrame must contain either 'CID' or 'cell_id' column")

    major = df[df['is_primary'] == True].copy().drop_duplicates(subset=['CID'])
    second = df[df['is_primary'] == False].copy().drop_duplicates(subset=['CID'])

    major_cids = set(major['CID'])
    promotable_seconds = second[~second['CID'].isin(major_cids)].copy()
    if not promotable_seconds.empty:
        promotable_seconds.loc[:, 'is_primary'] = True

    true_seconds = second[second['CID'].isin(major_cids)].copy()

    new_major_df = pd.concat([major, promotable_seconds], ignore_index=True)
    new_major_df = new_major_df.add_prefix('major_').rename(columns={'major_CID': 'CID'})
    true_seconds = true_seconds.add_prefix('second_').rename(columns={'second_CID': 'CID'})

    merged = pd.merge(new_major_df, true_seconds, on='CID', how='left')

    # Map contig ids
    if 'major_cell_id' in merged.columns:
        merged = merged.rename(columns={'major_cell_id': 'major_contig_id'})
    elif 'major_contig_id' not in merged.columns:
        merged['major_contig_id'] = merged.get('major_CID', merged['CID'])

    if 'second_cell_id' in merged.columns:
        merged = merged.rename(columns={'second_cell_id': 'second_contig_id'})
    elif 'second_contig_id' not in merged.columns:
        merged['second_contig_id'] = merged.get('second_CID', merged['CID'])

    desired_columns = [
        'sample',
        'major_v_call','major_j_call','major_c_call','major_junction','major_junction_aa',
        'major_abundance','major_contig_id','CID','major_is_primary','major_XPOS','major_YPOS','major_tcr_id',
        'second_v_call','second_j_call','second_c_call','second_junction','second_junction_aa',
        'second_abundance','second_contig_id','second_is_primary','second_XPOS','second_YPOS','second_tcr_id'
    ]

    for col in desired_columns:
        if col not in merged.columns:
            merged[col] = pd.NA

    return merged[desired_columns]


def organize_final_outputs(base_dir: str) -> None:
    """Produce 06_filtered_results/major_second_{TRA,TRB}.csv from 05_validated_tcr parquets.

    - Supports both combined cellular parquets with chain_type and per-chain CID parquets.
    - If chain_type not present, infer from filename suffix (_TRA/_TRB).
    """
    final_dir = get_output_path(base_dir, 'final_tcr', '')
    out_dir = os.path.join(base_dir, '06_filtered_results')
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(final_dir):
        log_message(f"Format Organizer: final dir missing {final_dir}")
        return

    parquet_files = sorted(
        [os.path.join(final_dir, f) for f in os.listdir(final_dir) if f.endswith('.parquet')]
    )
    if not parquet_files:
        log_message("Format Organizer: no parquet files found to organize")
        return

    frames: list[pd.DataFrame] = []
    for path in parquet_files:
        try:
            df = pd.read_parquet(path)
            # Track sample from filename if possible
            name = os.path.basename(path)
            sample_name = None
            if name.startswith('final_tcr_') and ('_TRA' in name or '_TRB' in name):
                # final_tcr_{sample}_{chain}.parquet
                trimmed = name[len('final_tcr_'):]
                sample_name = trimmed.rsplit('_', 1)[0]
            elif name.startswith('final_cellular_tcr_'):
                # final_cellular_tcr_{sample}.parquet
                sample_name = name[len('final_cellular_tcr_'):-len('.parquet')]
            if sample_name:
                df = df.copy()
                df['sample'] = sample_name
            # Infer chain if missing
            if 'chain_type' not in df.columns:
                chain = 'TRA' if '_TRA' in name else ('TRB' if '_TRB' in name else None)
                if chain is None:
                    # leave unspecified; will be dropped later
                    pass
                else:
                    df = df.copy()
                    df['chain_type'] = chain
            frames.append(df)
        except Exception as e:
            log_message(f"Format Organizer: failed reading {os.path.basename(path)}: {e}")

    if not frames:
        log_message("Format Organizer: no readable parquet frames")
        return

    df_all = pd.concat(frames, ignore_index=True)

    # Basic QC filter
    for chain in ['TRA', 'TRB']:
        sub = df_all[df_all.get('chain_type') == chain].copy()
        if sub.empty:
            log_message(f"Format Organizer: {chain} empty; skipping")
            continue
        # Filters
        sub = sub.dropna(subset=['XPOS', 'YPOS'])
        sub = sub[(sub['v_call'] != '*') & (sub['j_call'] != '*') & (sub['c_call'] != '*')]

        organized = _merge_major_second(sub)
        out_path = os.path.join(out_dir, f'major_second_{chain}.csv')
        organized.to_csv(out_path, index=False)
        log_message(f"Format Organizer: wrote {len(organized)} rows -> {os.path.basename(out_path)}")
    