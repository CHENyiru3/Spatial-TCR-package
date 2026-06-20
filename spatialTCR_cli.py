#!/usr/bin/env python3
"""
SpatialTCR Tool - Command Line Interface

Provides command-line access to both CID and cellular mode pipelines.
"""

import argparse
import sys
from pathlib import Path

from .modes.cid_mode import run_cid_pipeline
from .modes.cellular_mode import run_cellular_pipeline
from .core.utils import log_message


def main():
    """Main command-line interface."""
    parser = argparse.ArgumentParser(
        description="SpatialTCR Tool - TCR Analysis Pipeline for Spatial Transcriptomics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # CID mode - direct CID-level analysis
  python spatialTCR_cli.py cid --input config.json --output results/

  # Cellular mode - cell-level analysis with geometric boundaries  
  python spatialTCR_cli.py cellular --input config.json --output results/ --gdf cells.gdf

  # Show help for specific mode
  python spatialTCR_cli.py cid --help
  python spatialTCR_cli.py cellular --help
        """
    )
    
    # Add version argument
    parser.add_argument(
        '--version', 
        action='version', 
        version='SpatialTCR Tool v1.0.0'
    )
    
    # Create subparsers for different modes
    subparsers = parser.add_subparsers(
        dest='mode',
        help='Analysis mode',
        metavar='MODE'
    )
    
    # CID mode subparser
    cid_parser = subparsers.add_parser(
        'cid',
        help='CID-level TCR analysis (direct CID processing)',
        description='Analyze TCRs at the CID level without cellular aggregation'
    )
    
    cid_parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input JSON configuration file with sample information'
    )
    
    cid_parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for results'
    )
    
    # Cellular mode subparser
    cellular_parser = subparsers.add_parser(
        'cellular',
        help='Cellular-level TCR analysis (with cell boundary mapping)',
        description='Analyze TCRs at the cellular level using geometric cell boundaries'
    )
    
    cellular_parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input JSON configuration file with sample information'
    )
    
    cellular_parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for results'
    )
    
    cellular_parser.add_argument(
        '--gdf', '-g',
        required=True,
        help='Path to GeoDataFrame file (.gdf) containing cell boundary data'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Show help if no mode is specified
    if not args.mode:
        parser.print_help()
        sys.exit(1)
    
    # Validate input files exist
    if not Path(args.input).exists():
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)
    
    if args.mode == 'cellular' and not Path(args.gdf).exists():
        print(f"ERROR: GDF file not found: {args.gdf}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    Path(args.output).mkdir(parents=True, exist_ok=True)
    
    # Run the appropriate pipeline
    try:
        if args.mode == 'cid':
            log_message(f"Starting CID mode pipeline")
            log_message(f"Input: {args.input}")
            log_message(f"Output: {args.output}")
            
            success = run_cid_pipeline(args.input, args.output)
            
        elif args.mode == 'cellular':
            log_message(f"Starting cellular mode pipeline")
            log_message(f"Input: {args.input}")
            log_message(f"Output: {args.output}")
            log_message(f"GDF: {args.gdf}")
            
            success = run_cellular_pipeline(args.input, args.output, args.gdf)
        
        # Exit with appropriate code
        if success:
            log_message("Pipeline completed successfully!")
            sys.exit(0)
        else:
            log_message("Pipeline completed with errors or warnings.")
            sys.exit(1)
            
    except KeyboardInterrupt:
        log_message("Pipeline interrupted by user.")
        sys.exit(130)
        
    except Exception as e:
        log_message(f"Pipeline failed with error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()