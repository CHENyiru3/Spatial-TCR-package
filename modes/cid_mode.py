import json
import os
from typing import Mapping, Sequence

from ..core.utils import log_message, setup_output_directories, validate_file_exists
from ..processors.umi_processor import extract_umi_and_cid
from ..processors.spatial_processor import decode_h5_and_map_coordinates
from ..processors.trust4_processor import run_trust4_annotation
from ..processors.mapper_processor import run_cid_mapper
from ..processors.result_processor import finalize_cid_mode_results, validate_results, organize_final_outputs


def run_cid_pipeline(input_file: str, output_dir: str) -> bool:

    log_message("=== Starting SpatialTCR CID Mode Pipeline ===")

    # Load and validate configuration
    files_config = _load_and_validate_config(input_file)
    
    # Setup output directories
    setup_output_directories(output_dir)
    
    # Step 1: UMI-tools extraction
    extract_umi_and_cid(output_dir, files_config)
    
    # Step 2: H5 decoding and spatial mapping
    decode_h5_and_map_coordinates(output_dir, files_config)
    
    # Step 3: TRUST4 annotation
    run_trust4_annotation(output_dir, files_config)
    
    # Step 4: Result finalization
    finalize_cid_mode_results(output_dir, files_config)
    # Step 4.5: Organize final outputs into 06_filtered_results
    organize_final_outputs(output_dir)
    # Step 5: CID->Cell mapping and nearest remap
    run_cid_mapper(output_dir, files_config)
    
    # Step 6: Validation
    success = validate_results(output_dir, 'cid')
    
    if success:
        log_message("=== CID Mode Pipeline Completed Successfully ===")
    else:
        log_message("=== CID Mode Pipeline Completed with Warnings ===")
        
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


def run_full_pipeline(input_file: str, output_dir: str) -> bool:

    return run_cid_pipeline(input_file, output_dir)