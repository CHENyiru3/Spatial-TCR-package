#!/usr/bin/env python3
"""
Cellular Mode Example - Backward compatible with original cellular_level_annotation.py
"""

from spatialTCR_tool.modes.cellular_mode import run_full_pipeline


if __name__ == "__main__":
    # Example configurations - uncomment the one you want to run

    # Expression data with cellular analysis
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_tonsil.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/tonsil_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D04919G2.DAPI.cell_mask_nofilter.gdf"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_D03854D2.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/D03854D2_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D03854D2.DAPI.cell_mask.gdf"
    # )
    
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_thymus.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/thymus_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/C04895C1.DAPI.cell_mask.gdf"
    # )

    # Control data with cellular analysis
    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_D03854D2.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/control_D03854D2_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D03854D2.DAPI.cell_mask.gdf"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_thymus.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/control_thymus_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/C04895C1.DAPI.cell_mask.gdf"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_tonsil.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/control_tonsil_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D04919G2.DAPI.cell_mask_nofilter.gdf"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_D03854D1.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/pipeline_test",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D03854D1.DAPI_Image.cell_mask.gdf"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/expr_json/expr_D03854D1.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/D03854D1_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D03854D1.DAPI_Image.cell_mask.gdf"
    # )

    # run_full_pipeline(
    #     input_file="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/control_json/control_D03854D2.json",
    #     output_dir="/mnt/volume4/spatialTCR/spatialTCR_tool/upstream_data/Cell_level_result/control_D03854D2_rerun",
    #     gdf_path="/mnt/volume4/spatialTCR/spatialTCR_tool/SpatialT/stero-tcr/downstream/pipeline_outputs/gdfs/D03854D2.DAPI.cell_mask.gdf"
    # )