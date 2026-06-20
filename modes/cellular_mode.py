import gzip
import json
import os
import random
from typing import Dict, Mapping, Optional, Sequence, Tuple

import geopandas as gpd
import pandas as pd
import pyfastx
from shapely.geometry import Point
from shapely.affinity import scale

from ..core.utils import log_message, setup_output_directories, validate_file_exists, get_output_path, remove_file_if_exists
from ..core.config import CONFIG
from ..processors.umi_processor import extract_umi_and_cid
from ..processors.spatial_processor import decode_h5_and_map_coordinates
from ..processors.trust4_processor import run_trust4_annotation, validate_fastq_format
from ..processors.result_processor import finalize_cellular_mode_results, validate_results, organize_final_outputs


def run_cellular_pipeline(input_file: str, output_dir: str, gdf_path: str) -> bool:

    log_message("=== Starting SpatialTCR Cellular Mode Pipeline ===")
    
    # Load and validate inputs
    files_config = _load_and_validate_config(input_file)
    gdf = _load_and_validate_gdf(gdf_path)
    
    # Setup output directories
    setup_output_directories(output_dir)
    
    # Step 1: UMI-tools extraction
    extract_umi_and_cid(output_dir, files_config)
    
    # Step 2: H5 decoding and spatial mapping
    decode_h5_and_map_coordinates(output_dir, files_config)
    
    # Step 3: Cellular mapping and TRUST4 preparation
    fastq_preparations = _prepare_cellular_trust4_data(output_dir, files_config, gdf)
    
    # Step 4: TRUST4 annotation with cellular barcodes
    run_trust4_annotation(output_dir, files_config, fastq_preparations)
    
    # Step 5: Result finalization
    finalize_cellular_mode_results(output_dir, files_config, gdf)
    # Step 5.5: Organize final outputs into 06_filtered_results
    organize_final_outputs(output_dir)
    
    # Step 6: Validation
    success = validate_results(output_dir, 'cellular')
    
    return success
    

def _load_and_validate_config(input_file: str) -> Sequence[Mapping[str, str]]:
    """Load and validate the input configuration file."""
    validate_file_exists(input_file, "Input configuration file")
    
    with open(input_file, 'r') as f:
        files_config = json.load(f)
    
    if not isinstance(files_config, list):
        raise ValueError("Configuration file must contain a list of sample configurations")
    
    # Validate each file configuration  
    required_keys = ['name', '1', '2', 'h5']
    for i, file_config in enumerate(files_config):
        for key in required_keys:
            if key not in file_config:
                raise ValueError(f"Missing required key '{key}' in file configuration {i}")
        
        # Validate file paths
        for file_key in ['1', '2', 'h5']:
            validate_file_exists(file_config[file_key], f"Input file {file_key} for {file_config['name']}")
    
    log_message(f"Loaded configuration for {len(files_config)} samples")
    return files_config


def _load_and_validate_gdf(gdf_path: str) -> gpd.GeoDataFrame:
    """Load and validate the GeoDataFrame file."""
    validate_file_exists(gdf_path, "GeoDataFrame file")
    
    gdf = gpd.read_parquet(gdf_path)
    
    if 'id' not in gdf.columns:
        raise ValueError("GDF must contain an 'id' column for cell identification")
    
    if gdf.empty:
        raise ValueError("GDF file is empty")
    
    log_message(f"Loaded GDF with {len(gdf)} cell boundaries")
    return gdf


def _prepare_cellular_trust4_data(
    base_dir: str,
    files_config: Sequence[Mapping[str, str]], 
    gdf: gpd.GeoDataFrame
) -> Dict[str, Tuple[str, str]]:
    log_message("### Step 3: Preparing Cellular TRUST4 Data ###")
    
    fastq_preparations = {}
    
    for file_config in files_config:
        base_name = file_config['name']
        
        # Step 3a: Map CIDs to cells
        mapping_file = _map_cids_to_cells(base_dir, base_name, gdf)
        if not mapping_file:
            log_message(f"WARNING: No cellular mapping for {base_name}, skipping")
            continue
        
        # Step 3b: Generate TRUST4-ready FASTQ files
        r1_path, r2_path = _create_cellular_trust4_fastqs(base_dir, base_name, mapping_file)
        if r1_path and r2_path:
            fastq_preparations[base_name] = (r1_path, r2_path)
    log_message(f"Prepared {len(fastq_preparations)} samples for cellular TRUST4")
    return fastq_preparations


def _map_cids_to_cells(base_dir: str, base_name: str, gdf: gpd.GeoDataFrame) -> Optional[str]:
    """Map CIDs to cells using spatial boundaries."""
    log_message(f"Cellular Mapping: Starting for {base_name}")
    
    coords_path = get_output_path(base_dir, 'processed', f"{base_name}.coordinates.tsv")
    output_path = get_output_path(base_dir, 'processed', f'{base_name}.cid_to_cell_mapping.tsv')
    
    if not os.path.exists(coords_path):
        log_message(f"WARNING: Coordinates file not found for {base_name}")
        return None
    
    coords_df = pd.read_csv(coords_path, sep='\t')
    if coords_df.empty:
        log_message(f"WARNING: Empty coordinates file for {base_name}")
        return None
    
    # Scale cell boundaries by 2x for better capture
    gdf_scaled = gdf.copy()
    gdf_scaled['geometry'] = gdf_scaled['geometry'].apply(
        lambda geom: scale(geom, xfact=2, yfact=2, origin=(0, 0))
    )
    
    # Log coordinate ranges for diagnostics
    log_message(f"  CID coordinate range - X: [{coords_df['XPOS'].min():.2f}, {coords_df['XPOS'].max():.2f}], "
                f"Y: [{coords_df['YPOS'].min():.2f}, {coords_df['YPOS'].max():.2f}]")
    log_message(f"  Cell boundaries - X: [{gdf_scaled.total_bounds[0]:.2f}, {gdf_scaled.total_bounds[2]:.2f}], "
                f"Y: [{gdf_scaled.total_bounds[1]:.2f}, {gdf_scaled.total_bounds[3]:.2f}]")
    
    # Create spatial points
    coords_gdf = gpd.GeoDataFrame(
        coords_df,
        geometry=[Point(xy) for xy in zip(coords_df['XPOS'], coords_df['YPOS'])],
        crs=gdf_scaled.crs
    )
    
    # Perform spatial join
    result = gpd.sjoin(coords_gdf, gdf_scaled, how='left', predicate='within')
    
    # Report mapping statistics
    unmapped_count = result['id'].isna().sum()
    total_count = len(result)
    log_message(f"  Total CID coordinates: {total_count}")
    log_message(f"  Unmapped CIDs: {unmapped_count} ({unmapped_count/total_count*100:.1f}%)")
    
    # Keep only mapped entries
    mapped_result = result[['CID', 'id']].rename(columns={'id': 'cell_id'}).dropna(subset=['cell_id'])
    
    if mapped_result.empty:
        log_message(f"WARNING: No CIDs mapped to cells for {base_name}")
        return None
    
    mapped_result.to_csv(output_path, sep='\t', index=False)
    log_message(f"Cellular Mapping: Mapped {len(mapped_result)} reads to {mapped_result['cell_id'].nunique()} cells.")
    return output_path


def _generate_cellular_barcode(idx: int, length: int = 20) -> str:
    """Generate a pseudo-barcode for cellular analysis."""
    random.seed(CONFIG['barcode_generation']['seed_base'] + idx)
    return ''.join(random.choices(CONFIG['barcode_generation']['alphabet'], k=length))


def _create_cellular_trust4_fastqs(
    base_dir: str, 
    base_name: str, 
    mapping_file_path: str
) -> Tuple[Optional[str], Optional[str]]:
    """Create TRUST4-compatible FASTQ files with cellular barcodes."""
    barcode_length = CONFIG['barcode_generation']['default_length']
    umi_length = 10
    expected_r1_length = barcode_length + umi_length
    
    log_message(f"TRUST4 Prep: Rewriting FASTQ for {base_name} - "
                f"Direct CID replacement (AGCT barcode, {barcode_length}bp)")
    
    # Input and output file paths
    r1_in_path = get_output_path(base_dir, 'umi', f"{base_name}_1.fq.gz")
    r2_in_path = get_output_path(base_dir, 'umi', f"{base_name}_mapped_2.fq.gz")
    r1_out_path = get_output_path(base_dir, 'processed', f"{base_name}.cellular_trust4_R1.fq.gz")
    r2_out_path = get_output_path(base_dir, 'processed', f"{base_name}.cellular_trust4_R2.fq.gz")
    
    # Remove existing output files
    for out_path in [r1_out_path, r2_out_path]:
        remove_file_if_exists(out_path)
    
    # Validate input files
    if not os.path.exists(r1_in_path) or not os.path.exists(r2_in_path):
        log_message(f"ERROR: Input files missing for {base_name}")
        return None, None
    
    # Load mapping data
    mapping_df = pd.read_csv(mapping_file_path, sep='\t', dtype={'cell_id': str})
    cid_to_cell_dict = pd.Series(mapping_df.cell_id.values, index=mapping_df.CID).to_dict()
    
    # Generate cellular barcodes
    unique_cells = sorted(mapping_df['cell_id'].unique())
    cell_to_barcode = {
        cell_id: _generate_cellular_barcode(idx, barcode_length) 
        for idx, cell_id in enumerate(unique_cells)
    }
    cid_to_barcode = {
        cid: cell_to_barcode[cell_id] 
        for cid, cell_id in cid_to_cell_dict.items() 
        if cell_id in cell_to_barcode
    }
    
    log_message(f"Generated {len(cell_to_barcode)} unique cellular barcodes ({barcode_length}bp each)")
    
    reads_written = 0
    unmapped_reads = 0

    fq1 = pyfastx.Fastq(r1_in_path, build_index=False)
    fq2 = pyfastx.Fastq(r2_in_path, build_index=False)
    
    with gzip.open(r1_out_path, 'wt') as f_out1, gzip.open(r2_out_path, 'wt') as f_out2:
        for (h1, s1, q1), (h2, s2, q2) in zip(fq1, fq2):
            # Extract CID from header
            cid = h1.split('_')[1][:25] if '_' in h1 else None
            
            if cid and cid in cid_to_barcode:
                pseudo_cell_barcode = cid_to_barcode[cid]
                
                # Extract UMI from original sequence
                umi = s1[25:35] if len(s1) >= 35 else s1[-10:].ljust(10, 'N')
                
                # Update headers
                def replace_cid(header, old_cid, new_barcode):
                    return header.replace(old_cid, new_barcode, 1)
                
                h1_new = replace_cid(h1, cid, pseudo_cell_barcode)
                h2_new = replace_cid(h2, cid, pseudo_cell_barcode)
                
                # Create new R1 sequence (barcode + UMI)
                new_r1_seq = pseudo_cell_barcode + umi
                new_r1_seq = new_r1_seq[:expected_r1_length].ljust(expected_r1_length, 'N')
                new_r1_qual = 'F' * expected_r1_length
                new_r2_qual = 'F' * len(s2)
                
                # Write FASTQ records
                f_out1.write(f"@{h1_new}\n{new_r1_seq}\n+\n{new_r1_qual}\n")
                f_out2.write(f"@{h2_new}\n{s2}\n+\n{new_r2_qual}\n")
                
                reads_written += 1
            else:
                unmapped_reads += 1

    
    log_message(f"TRUST4 Prep: Wrote {reads_written} reads for {base_name} ({unmapped_reads} unmapped).")
    
    # Save barcode mapping files
    _save_barcode_mappings(base_dir, base_name, cell_to_barcode, cid_to_barcode)
    
    # Validate output
    if reads_written == 0:
        log_message(f"ERROR: No reads written for {base_name}")
        return None, None
    
    # Validate FASTQ format
    if not validate_fastq_format(r1_out_path, expected_r1_length):
        log_message(f"ERROR: FASTQ format validation failed for {base_name}")
        return None, None
    
    return r1_out_path, r2_out_path


def _save_barcode_mappings(
    base_dir: str,
    base_name: str,
    cell_to_barcode: Dict[str, str],
    cid_to_barcode: Dict[str, str]
) -> None:
    """Save barcode mapping files for reference."""
    # Cell to barcode mapping
    barcode_map_path = get_output_path(base_dir, 'processed', f"{base_name}.cell_barcode_map.tsv")
    with open(barcode_map_path, 'w') as f:
        f.write("cell_id\tpseudo_barcode\n")
        for cell_id, barcode in cell_to_barcode.items():
            f.write(f"{cell_id}\t{barcode}\n")
    
    # CID to barcode mapping  
    cid_barcode_map_path = get_output_path(base_dir, 'processed', f"{base_name}.cid_barcode_map.tsv")
    with open(cid_barcode_map_path, 'w') as f:
        f.write("CID\tpseudo_barcode\n")
        for cid, barcode in cid_to_barcode.items():
            f.write(f"{cid}\t{barcode}\n")


def run_full_pipeline(input_file: str, output_dir: str, gdf_path: str) -> bool:

    return run_cellular_pipeline(input_file, output_dir, gdf_path)