from typing import Dict, Any

CONFIG = {
    'version': '1.0.0',
    'parallel_jobs': {
        'umi_jobs': 'cpu_count',  # Use os.cpu_count() for UMI processing
        'trust4_jobs': 4,         # Use 4 jobs for TRUST4 (from raw_run_umi_tool.py)
        'cellular_trust4_jobs': 'cpu_count_div_4',  # Use cpu_count()//4 for cellular mode TRUST4
    },
    'validation': {
        'min_tcr_length': 50,
        'max_tcr_length': 500,
        'expected_fastq_seq_len': 30, 
    },
    'barcode_generation': {
        'default_length': 20,
        'alphabet': 'ACGT',
        'seed_base': 42, 
    }
}

FILE_PATTERNS = {
    'fastq_extensions': ['.fq', '.fastq', '.fq.gz', '.fastq.gz'],
    'output_extensions': {
        'coordinates': '.coordinates.tsv',
        'mapping': '.cid_to_cell_mapping.tsv',
        'barcode_map': '.cell_barcode_map.tsv',
        'cid_barcode_map': '.cid_barcode_map.tsv',
        'trust4_r1': '.cellular_trust4_R1.fq.gz',
        'trust4_r2': '.cellular_trust4_R2.fq.gz',
        'final_tcr': '.parquet'
    }
}

# Validation thresholds
VALIDATION_THRESHOLDS = {
    'min_reads_per_sample': 1000,
    'min_cells_per_sample': 10,
    'max_unmapped_rate': 0.95,  
    'min_tcr_abundance': 0.1
}

def get_config_value(key_path: str, default: Any = None) -> Any:
    """
    Get a configuration value using dot notation.
    
    Args:
        key_path: Dot-separated path to config value (e.g., 'parallel_jobs.trust4_jobs')
        default: Default value if key not found
        
    Returns:
        Configuration value or default
    """
    keys = key_path.split('.')
    value = CONFIG
    
    for key in keys:
        value = value[key]
    return value
