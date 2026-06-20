#!/usr/bin/env python3
"""
CID Mode Example - Backward compatible with original raw_run_umi_tool.py

This script maintains the same interface as the original raw_run_umi_tool.py
while using the new modular architecture.
"""

from spatialTCR_tool.modes.cid_mode import run_full_pipeline

if __name__ == "__main__":
    # Example configurations - uncomment the one you want to run
    
    # Test with thymus subsampling data
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/no_filter_test/test_thymus_subsampling/testing.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/no_filter_test/test_thymus_subsampling/version3_output"
    # )

    # Expression data samples
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_D03854D1.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_results/D03854D1_rerun"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_tonsil.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_results/tonsil_rerun"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_D03854D2.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_results/D03854D2_rerun"
    # )
    
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_thymus.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_results/thymus_rerun"
    # )

    # Control data samples
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_D03854D2.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/CID_level_result/control_D03854D2_rerun"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_thymus.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/CID_level_result/control_thymus_rerun"
    # )

    run_full_pipeline(
        input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_tonsil.json",
        output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/CID_level_result/control_tonsil_rerun"
    )