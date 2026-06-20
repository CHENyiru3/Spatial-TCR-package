import datetime
import os
import subprocess
from typing import Any, List, Optional

# Standardized output directory structure
OUTPUT_STRUCTURE = {
    'subdirs': {
        'logs': 'logs',
        'umi': "02_processed",
        'processed': '02_processed', 
        'trust4_raw': '03_trust4_results',
        'final_tcr': '05_validated_tcr'
    }
}

# UMI-tools regex patterns for read extraction
UMI_PATTERNS = {
    'DISCARD_MGI_ADAPTOR': r".+",
    'CELL_1': r"(?P<cell_1>.{25})",
    'DISCARD_INTEM_SEQ': r"(TTGTCTTCCTAAGAC){e<=2}",
    'UMI_1': r"(?P<umi_1>.{10})",
    'DISCARD_OLIGO_DT': r"T+",
}

# Tool paths - update these based on your environment
TOOL_PATHS = {
    'ST_BARCODEMAP': "/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/SpatialTCR_data/upstream_results/05_scripts/ST_BarcodeMap-0.0.1",
    'REWRITE_R1_SCRIPT': "/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/SpatialTCR_data/upstream_results/05_scripts/rewrite_r1_cid_umi.sh",
    'TRUST4_SCRIPT': "/mnt/volume4/spatialTCR/spatialTCR_tool/spatialTCR_upstream/trust4_pipeline.py",
    'TRUST4_TCR_REF': "/mnt/volume3/trn/04_benchmark/01_trust4/refdata/tcr.fa",
    'TRUST4_IMGT_REF': "/mnt/volume3/trn/04_benchmark/01_trust4/refdata/IMGT+C_tcr.fa"
}

# Constants
UMI_LENGTH = 10

# TRUST4 read formats for different modes
TRUST4_READ_FORMATS = {
    'CID_MODE': "bc:0:24,um:25:34,r2:0:-1",  # For CID mode (25 char barcode + 10 char UMI)
    'CELLULAR_MODE': "bc:0:19,um:20:29",     # For cellular mode (20 char barcode + 10 char UMI)
}


def log_message(message: str) -> None:
    """Print a message with timestamp."""
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")


def get_output_path(base_dir: str, subdir_key: str, filename: str) -> str:
    """Construct an output path using the standardized directory structure."""
    if subdir_key not in OUTPUT_STRUCTURE['subdirs']:
        raise ValueError(f"Unknown subdirectory key: {subdir_key}")
    return os.path.join(base_dir, OUTPUT_STRUCTURE['subdirs'][subdir_key], filename)


def ensure_directory_exists(directory_path: str) -> None:
    """Create directory if it doesn't exist."""
    os.makedirs(directory_path, exist_ok=True)


def setup_output_directories(base_dir: str) -> None:
    """Create all required output directories."""
    for subdir in OUTPUT_STRUCTURE['subdirs'].values():
        ensure_directory_exists(os.path.join(base_dir, subdir))


def execute_command(
    command: List[str], 
    log_file: str, 
    description: Optional[str] = None,
    capture_output: bool = False
) -> subprocess.CompletedProcess:
    if description:
        log_message(f"Executing: {description}")
    
    log_message(f"CMD: {' '.join(command)}")
    

    if capture_output:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
    else:
        with open(log_file, "w") as f:
            result = subprocess.run(command, check=True, stdout=f, stderr=f, text=True)
    
    return result
        


def validate_file_exists(file_path: str, description: Optional[str] = None) -> None:
    """Validate that a file exists, raising an error if not."""
    if not os.path.exists(file_path):
        desc = description or "File"
        raise FileNotFoundError(f"{desc} not found: {file_path}")


def get_file_size(file_path: str) -> int:
    """Get file size in bytes."""
    return os.path.getsize(file_path) if os.path.exists(file_path) else 0


def remove_file_if_exists(file_path: str) -> bool:
    """Remove file if it exists. Returns True if file was removed."""
    if os.path.exists(file_path):
        os.remove(file_path)
        log_message(f"Removed old file: {file_path}")
        return True
    return False