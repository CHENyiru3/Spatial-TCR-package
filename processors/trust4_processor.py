import os
import subprocess
from typing import List, Mapping, Optional, Sequence, Tuple

import joblib

from ..core.utils import log_message, get_output_path, execute_command, validate_file_exists, get_file_size, TOOL_PATHS, TRUST4_READ_FORMATS
from ..core.config import get_config_value


def run_trust4_annotation(
    base_dir: str,
    files: Sequence[Mapping[str, str]],
    fastq_preparations: Optional[dict] = None
) -> None:

    log_message("### Step 3: TRUST4 TCR Annotation ###")
    
    def run_single_trust4_task(file_config: Mapping[str, str]) -> Tuple[bool, str, str]:
        """Run TRUST4 for a single sample."""
        base_name = file_config['name']
    
        # Determine input files
        if fastq_preparations and base_name in fastq_preparations:
            r1_path, r2_path = fastq_preparations[base_name]
            read_format = TRUST4_READ_FORMATS['CELLULAR_MODE']  # For cellular mode with 20bp barcodes
        else:
            # Standard CID mode
            r1_path = get_output_path(base_dir, 'umi', f"{base_name}_mapped_1.fq.gz")
            r2_path = get_output_path(base_dir, 'umi', f"{base_name}_mapped_2.fq.gz") 
            read_format = TRUST4_READ_FORMATS['CID_MODE']  # For CID mode with 25bp barcodes
        
        # Validate input files
        validate_file_exists(r1_path, f"R1 file for {base_name}")
        validate_file_exists(r2_path, f"R2 file for {base_name}")
        
        # Build TRUST4 command
        trust4_output_dir = get_output_path(base_dir, 'trust4_raw', '')
        log_file = get_output_path(base_dir, 'logs', f'{base_name}.trust4_run.log')
        
        command = [
            "run-trust4",
            "-1", r1_path,
            "-2", r2_path,
            "--barcode", r1_path,
            "-f", TOOL_PATHS['TRUST4_TCR_REF'],
            "--ref", TOOL_PATHS['TRUST4_IMGT_REF'],
            "--od", trust4_output_dir,
            "-o", base_name,
            "-t", "16"
        ]
        
        # Add read format for cellular mode
        if read_format:
            command.extend(["--readFormat", read_format])
        
        log_message(f"TRUST4: Running for {base_name}")
        log_message(f"TRUST4 CMD: {' '.join(command)}")
        
        execute_command(command, log_file, f"TRUST4 for {base_name}")
        
        # Verify output
        report_path = os.path.join(trust4_output_dir, f'{base_name}_barcode_report.tsv')
        report_size = get_file_size(report_path)
        
        log_message(f"TRUST4: Completed for {base_name}. Report size: {report_size} bytes")
        return (True, base_name, f"Success: {report_size} bytes")
        
    # Run TRUST4 jobs in parallel with original script settings
    if fastq_preparations:
        # Cellular mode: Use cpu_count()//4 as in original cellular_level_annotation.py
        n_jobs = max(1, min(len(files), os.cpu_count() // 4))
    else:
        # CID mode: Use 4 jobs as in original raw_run_umi_tool.py  
        n_jobs = 4
    
    log_message(f"Running TRUST4 with {n_jobs} parallel jobs")
    
    results = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(run_single_trust4_task)(file_config) for file_config in files
    )
    
    # Report results
    successful = sum(1 for success, _, _ in results if success)
    total = len(results)
    
    log_message(f"TRUST4: Completed {successful}/{total} jobs successfully")
    
    # Log any failures
    for success, base_name, message in results:
        if not success:
            log_message(f"TRUST4 FAILED - {base_name}: {message}")
    
    log_message("### Step 3: Completed ###")


def validate_fastq_format(fastq_path: str, expected_seq_len: Optional[int] = None) -> bool:
    """
    Validate FASTQ file format and sequence lengths.
    
    Args:
        fastq_path: Path to FASTQ file
        expected_seq_len: Expected sequence length (optional)
        
    Returns:
        True if validation passes
    """
    import pyfastx
    
    fq = pyfastx.Fastq(fastq_path, build_index=False)
    
    # Check first few sequences
    count = 0
    for name, seq, qual in fq:
        if expected_seq_len and len(seq) != expected_seq_len:
            log_message(f"WARNING: Sequence length mismatch in {fastq_path}. "
                        f"Expected: {expected_seq_len}, Got: {len(seq)}")
            return False
        
        count += 1
        if count >= 10:  # Check first 10 sequences
            break
    
    return True
        
