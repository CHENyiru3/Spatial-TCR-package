import csv
import gzip
import os
import re
from typing import Mapping, Sequence

import pandas as pd

from ..core.utils import log_message, get_output_path, execute_command, validate_file_exists, TOOL_PATHS, ensure_directory_exists


def decode_h5_and_map_coordinates(base_dir: str, files: Sequence[Mapping[str, str]]) -> None:
    log_message("### Step 2: H5 Decoding and Spatial Mapping ###")
    
    for file_config in files:
        base_name = file_config['name']
        log_message(f"H5 Decode: Starting spatial mapping for {base_name}")
        
        # Validate required files
        validate_file_exists(file_config['h5'], "H5 file")
        r1_path = get_output_path(base_dir, 'umi', f"{base_name}_1.fq.gz")
        r2_path = get_output_path(base_dir, 'umi', f"{base_name}_2.fq.gz")
        validate_file_exists(r1_path, "UMI-processed R1 file")
        validate_file_exists(r2_path, "UMI-processed R2 file")
        
        # Run ST_BarcodeMap for R2 (TCR sequences with coordinates)
        _run_barcodemap(
            base_dir=base_dir,
            base_name=base_name,
            h5_file=file_config['h5'],
            input_r1=r1_path,
            input_r2=r2_path,
            output_file=get_output_path(base_dir, 'umi', f"{base_name}_mapped_2.fq.gz"),
            description="R2 mapping"
        )
        
        # Run ST_BarcodeMap for R1 (coordinates reference)
        _run_barcodemap(
            base_dir=base_dir,
            base_name=base_name, 
            h5_file=file_config['h5'],
            input_r1=r1_path,
            input_r2=r1_path,  # Use R1 as both inputs
            output_file=get_output_path(base_dir, 'umi', f"{base_name}_mapped_1.fq.gz"),
            description="R1 mapping"
        )
        
        # Extract coordinates from mapped reads
        extract_coordinates_from_mapped_reads(base_dir, base_name)
        log_message(f"H5 Decode: Completed spatial mapping for {base_name}")
    
    log_message("### Step 2: Completed ###")


def _run_barcodemap(
    base_dir: str,
    base_name: str,
    h5_file: str,
    input_r1: str,
    input_r2: str,
    output_file: str,
    description: str
) -> None:
    """Run ST_BarcodeMap for spatial coordinate mapping."""
    command = [
        TOOL_PATHS['ST_BARCODEMAP'],
        "--in", h5_file,
        "--in1", input_r1,
        "--in2", input_r2,
        "--out", output_file,
        "--mismatch", "1",
        "--umiStart", "25",
        "--thread", "12"
    ]
    
    log_file = get_output_path(base_dir, 'logs', f'barcode_map_{description}_{base_name}.log')
    execute_command(command, log_file, f"ST_BarcodeMap {description} for {base_name}")


def extract_coordinates_from_mapped_reads(base_dir: str, base_name: str) -> None:
    """
    Extract spatial coordinates from ST_BarcodeMap output.
    
    Args:
        base_dir: Base output directory
        base_name: Sample name
    """
    mapped_fq_path = get_output_path(base_dir, 'umi', f"{base_name}_mapped_1.fq.gz")
    coordinate_path = get_output_path(base_dir, 'processed', f"{base_name}.coordinates.tsv")
    mapping_path = get_output_path(base_dir, 'processed', f"{base_name}_readid2cid.tsv")
    
    log_message(f"Extracting coordinates from {os.path.basename(mapped_fq_path)} to TSV")
    
    validate_file_exists(mapped_fq_path, "Mapped FASTQ file")
    
    # Ensure output directory exists
    ensure_directory_exists(os.path.dirname(coordinate_path))
    
    # Regex pattern for coordinate extraction
    pattern = re.compile(r'CB:Z:(\d+_\d+)')
    coordinate_count = 0


    with open(coordinate_path, 'w', newline='') as coord_file, \
            open(mapping_path, 'w') as map_file:
        
        coord_writer = csv.writer(coord_file, delimiter='\t')
        coord_writer.writerow(["CID", "XPOS", "YPOS"])
        
        # Process FASTQ file
        with gzip.open(mapped_fq_path, 'rt') as fq_file:
            lines = fq_file.readlines()
            
            for i in range(0, len(lines), 4):
                if i < len(lines) and lines[i].startswith('@'):
                    header = lines[i].strip()
                    
                    # Extract coordinate from CB:Z tag
                    coord_match = pattern.search(header)
                    
                    # Extract CID from read name
                    read_name = header.split()[0][1:]  # Remove '@'
                    name_parts = read_name.split('_')
                    
                    if len(name_parts) >= 2:
                        cid = name_parts[1][:25]  # First 25 chars as CID
                        read_id = name_parts[0]
                        map_file.write(f"{read_id}\t{cid}\n")
                        
                        if coord_match:
                            coord_str = coord_match.group(1)
                            x, y = coord_str.split('_')
                            coord_writer.writerow([cid, x, y])
                            coordinate_count += 1

    
    # Remove duplicates
    if coordinate_count > 0:
        df = pd.read_csv(coordinate_path, sep='\t')
        original_count = len(df)
        df.drop_duplicates(subset=['CID', 'XPOS', 'YPOS'], inplace=True)
        df.to_csv(coordinate_path, sep='\t', index=False)
        
        log_message(f"Extracted {coordinate_count} coordinate entries, "
                    f"{len(df)} unique entries for {base_name}")
    else:
        log_message(f"WARNING: No coordinates extracted for {base_name}")
