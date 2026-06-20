import os
import subprocess
from typing import List, Mapping, Sequence

import joblib

from ..core.utils import log_message, get_output_path, execute_command, validate_file_exists, UMI_PATTERNS, TOOL_PATHS
from ..core.config import get_config_value


def extract_umi_and_cid(base_dir: str, files: Sequence[Mapping[str, str]]) -> None:
    """
    Process raw FASTQ files with UMI-tools to extract CID and UMI into read headers.
    
    Args:
        base_dir: Base output directory
        files: List of file configurations with 'name', '1', '2' keys
    """
    log_message("### Step 1: UMI-tools Extraction ###")
    
    def process_single_file(file_config: Mapping[str, str]) -> None:
        """Process a single sample with UMI-tools."""
        base_name = file_config['name']
        log_message(f"UMI-tools: Starting extraction for {base_name}")
        
        # Validate input files
        validate_file_exists(file_config['1'], "R1 input file")
        validate_file_exists(file_config['2'], "R2 input file")
        
        # Build UMI-tools command
        bc_pattern = "".join([
            "^",
            UMI_PATTERNS['DISCARD_MGI_ADAPTOR'],
            UMI_PATTERNS['CELL_1'],
            UMI_PATTERNS['DISCARD_INTEM_SEQ'],
            UMI_PATTERNS['UMI_1'],
            UMI_PATTERNS['DISCARD_OLIGO_DT'],
            ".*$"
        ])
        
        command = [
            "umi_tools", "extract", "--extract-method=regex",
            f"--bc-pattern={bc_pattern}",
            f"--log={get_output_path(base_dir, 'logs', f'{base_name}.extract.log')}",
            f"--stdin={file_config['1']}",
            f"--read2-in={file_config['2']}",
            f"--stdout={get_output_path(base_dir, 'umi', f'{base_name}_1.fq.gz')}",
            f"--read2-out={get_output_path(base_dir, 'umi', f'{base_name}_2.fq.gz')}",
            f"--filtered-out={get_output_path(base_dir, 'umi', f'{base_name}_1.failed.fq.gz')}",
            f"--filtered-out2={get_output_path(base_dir, 'umi', f'{base_name}_2.failed.fq.gz')}",
        ]
        
        execute_command(
            command, 
            get_output_path(base_dir, 'logs', f'{base_name}.umi_extract.log'),
            f"UMI-tools extraction for {base_name}"
        )
        log_message(f"UMI-tools: Completed extraction for {base_name}")
        
        # Post-process R1 for spatial mapping
        rewrite_r1_for_mapping(base_dir, base_name)
        

    
        # Use os.cpu_count() as in original script
    n_jobs = os.cpu_count()
    joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(process_single_file)(file_config) 
        for file_config in files
    )
    
    log_message("### Step 1: Completed ###")


def rewrite_r1_for_mapping(base_dir: str, base_name: str) -> None:
    """
    Rewrite R1 FASTQ to contain only CID+UMI sequence for spatial mapping.
    
    Args:
        base_dir: Base output directory
        base_name: Sample name
    """
    log_message(f"Post-processing: Rewriting R1 for {base_name} (CID+UMI only)")
    
    r1_path = get_output_path(base_dir, 'umi', f"{base_name}_1.fq.gz")
    r1_temp_path = r1_path + ".tmp"
    
    validate_file_exists(r1_path, "UMI-tools output R1 file")
    validate_file_exists(TOOL_PATHS['REWRITE_R1_SCRIPT'], "R1 rewrite script")
    
    command = ["bash", TOOL_PATHS['REWRITE_R1_SCRIPT'], r1_path, r1_temp_path]

    result = execute_command(
        command,
        get_output_path(base_dir, 'logs', f'{base_name}.r1_rewrite.log'),
        f"R1 rewrite for {base_name}",
        capture_output=True
    )
    
    if result.stdout:
        log_message(f"R1 rewrite output: {result.stdout.strip()}")
    
    # Validate output and replace original
    validate_file_exists(r1_temp_path, "R1 rewrite output")
    os.replace(r1_temp_path, r1_path)
    
    log_message(f"Post-processing: Successfully rewrote R1 for {base_name}")
