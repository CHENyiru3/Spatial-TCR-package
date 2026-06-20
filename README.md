# SpatialTCR Tool - Modular Architecture

A comprehensive tool for spatial T-cell receptor (TCR) analysis with support for both CID-based and cellular-level annotation modes. Copyright@LabW, Zhejiang University


1. **Run CID mode analysis:**
   ```bash
   python spatialTCR_cli.py cid --input config.json --output results/
   ```

2. **Run cellular mode analysis:**
   ```bash
   python spatialTCR_cli.py cellular --input config.json --output results/ --gdf cells.gdf
   ```

## Architecture Overview

```
spatialTCR_tool/
├── core/                     # Core utilities and configuration
│   ├── utils.py             # Logging and utility functions
│   └── config.py            # Default configuration settings
├── processors/              # Processing modules
│   ├── umi_processor.py     # UMI extraction and CID processing
│   ├── spatial_processor.py # Spatial coordinate mapping
│   ├── trust4_processor.py  # TCR annotation processing
│   └── result_processor.py  # Result aggregation and output
├── modes/                   # Analysis mode implementations
│   ├── cid_mode.py         # CID-based analysis pipeline
│   └── cellular_mode.py    # Cellular-level analysis pipeline
├── examples/               # Usage examples and backward compatibility
│   ├── legacy_cid_api.py   # Original CID mode API
│   ├── legacy_cellular_api.py # Original cellular mode API
│   └── usage_examples.py   # Modern API examples
├── spatialTCR_cli.py       # Main CLI interface
├── setup_check.py          # Setup validation script
└── README.md              # This file
```

## Dependencies

### Python Packages
- **pandas**: Data manipulation and analysis
- **geopandas**: Spatial data operations
- **joblib**: Parallel processing support
- **shapely**: Geometric operations
- **pyfastx**: FASTQ file processing (optional for basic testing)

### External Tools
- **UMI-tools**: UMI extraction and deduplication
- **TRUST4**: TCR assembly and annotation
## 📖 Usage Modes

### CID Mode
Processes UMI data and performs CID-based TCR analysis:
```bash
python spatialTCR_cli.py cid --input config.json --output results/ [--threads 4] [--debug]
```

### Cellular Mode
This mode will annotate the cell unique id to CIDs first, so that we can regard them as single cells.

Performs cellular-level TCR annotation with spatial mapping:
```bash
python spatialTCR_cli.py cellular --input config.json --output results/ --gdf cells.gdf [--threads 4] [--debug]
```


## 🔧 Advanced Usage

### Programmatic API
```python
from modes.cid_mode import CIDModeRunner
from modes.cellular_mode import CellularModeRunner

# CID mode
runner = CIDModeRunner()
results = runner.run(config, output_dir, threads=4)

# Cellular mode
runner = CellularModeRunner()
results = runner.run(config, output_dir, gdf_file, threads=4)
```
