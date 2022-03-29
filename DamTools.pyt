# -*- coding: utf-8 -*-

##########################################
# Automated Dam Reservoir estimation Toolbox
# Original Author: Frithjof POLLMÜLLER (2019, Master Thesis)
# Updated to process multiple dams: Leonhard Urs BÜRGER (2021, Master Project)
# Updated to Python3 : Leonhard Urs BÜRGER (2022)
# Supervisor: Jürgen BERLEKAMP
# Version: March 2022
# for ArcGIS Pro 2.8.0 (Esri Inc.)
##########################################

import arcpy
import numpy as np
import math
import pandas as pd
import scipy.optimize as sciopt
import shutil

import time
import os
import matplotlib
import csv
import traceback

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

from collections import deque
import zipfile
import shutil
import sys

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "MScProject"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Enumerate_river_segments, Automated_reservoir_estimation, Combine_dem_rasters,
         Merge_results]

class Automated_reservoir_estimation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Automated reservoir estimation (2022.03)"
        self.description = "A tool to estimate potential hydropower reservoirs (Bürger, Pollmueller + Berlekamp)"
        self.canRunInBackground = False
        arcpy.env.overwriteOutput = True # override

        ## Change Variables here to match the Field names in the dam feature class
        # given dam height [m]
        self.height_field = "g_hgt_m" #"DAM_HGT_M"#"Height"# "Height_m" #"Height"
        # given hydropower capacity [MW]
        self.cap_field = "g_cap_mw"#"Capacity__MW_"
        # unique ID
        self.id_field = "HPPD_ID" #"GRAND_ID"#"DAM_ID" #"GRAND_ID" #"DAM_ID"
        # latitude of the dam position
        self.lat_field = "lat"#"LAT_DD" #"Lat_2016" #"LAT_DD" #"Lat_2015"      # former: "Current_LAT"

        ## can be left unchanged when working with the HydroATLAS and SRTMGL3 data
        self.basin_field = "HYBAS_ID"
        self.discharge_field = "dis_m3_pmx" # used for estimation
        self.discharge_pmx = "dis_m3_pmx"
        self.discharge_pmn = "dis_m3_pmn"
        self.discharge_pyr = "dis_m3_pyr"
        self.hydroid_field = "HYRIV_ID"
        self.nextdown_field = "NEXT_DOWN"
        self.dem_cell_size = 92  # cell dimensions at the equator

        ## output field names
        self.out_vol_field = "r_vol_mcm" #"FHReD_vol_mcm"
        self.out_area_field = "r_area_msm" #"FHReD_area_msm"
        self.out_hyid_field = "near_hy_ID"
        self.out_est_key = "d_est_key"

        ## Min and max values for possible dam heights in meters
        self.min_dam_height = 5 ## former 20
        self.max_dam_height = 300  # meters
        self.max_dam_width = 8000  # meters
        self.max_width_x_height = 300000
        self.max_reservoir_area = 10000  # KM^2
        self.max_reservoir_volume = 250  # KM^3
        self.max_dam_interruption = 0.2  # Share of profile points above dam line segments

        # Variables to be used by the process
        self.logger = Logger()
        self.project_folder = None
        self.project_name = None
        self.project_path = None
        self.graphs_path = None
        self.scratch_path = None
        self.out_gdb = None
        self.dam_positions = None
        self.dam_layer_name = None
        self.dam_pos_basin = None
        self.dam_dict = None
        self.river_seg_dict = None
        self.dem = None
        self.dem_basin = None
        self.dem_watershed_clip = None
        self.river_seg = None
        self.river_seg_basin = None
        #self.res_overlapping = False
        self.basins = None
        self.current_basin_id = None
        self.height_step = None
        self.id_list = None
        self.id_list_basin = None
        self.continent = None
        self.keep_tmp_files = False
        self.not_calculated_string = "not calculated"
        self.snap_range = 0.001

        self.upriver_segments = None
        self.rs_buffer_dist = "10000 Meters"

        self.initial_dam_snap_range = "1000 Meters"



    def getParameterInfo(self):
        """Define parameter definitions"""
        # Dam Positions - Points of dam locations
        param_dam_positions = arcpy.Parameter(
            displayName="Dam Positions",
            name="in_dam_positions",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")

        # DEM - The digital elevation model
        param_dem = arcpy.Parameter(
            displayName="SRTMGL3 DEM",
            name="in_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        # RiverATLAS river segments
        param_river_segment = arcpy.Parameter(
            displayName = "RiverATLAS river segments",
            name = "in_river_segments",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Input")

        # Height Step - Define at which intervals potential reservoirs should be validated. Impacts runtime.
        param_height_step = arcpy.Parameter(
            displayName="Height Step [Meters]",
            name="in_height_step",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param_height_step.value = 2

        param_upriver_segments = arcpy.Parameter(
            displayName = "Number of upriver river segments defining the initial search area",
            name="in_upriver_segments",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )
        param_upriver_segments.value = 700

        param_reservoir_overlapping = arcpy.Parameter(
            displayName = "Try to minimize reservoir overlapping? (max. 20 km)",
            name="in_reservoir_overlapping",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        #param_reservoir_overlapping.value = False

        param_keep_tmp_files = arcpy.Parameter(
            displayName = "Kepp all temporary files for debugging?",
            name="in_keep_tmp_files",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        )
        param_keep_tmp_files.value = False

        # Dam ID List - Comma seperated list of dam IDs to only run tool for certian dams
        param_dam_id_list = arcpy.Parameter(
            displayName="Dam ID List",
            name="in_dam_id_list",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        #param_dam_id_list.value = "3103"#"4031, 4030"

        # Project Name - Used for naming output folders and files
        param_project_name = arcpy.Parameter(
            displayName="Project Name",
            name="in_project_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param_project_name.value = "0110_"

        now = time.strftime("%c")
        now = now.replace(":", "-")
        now = now.replace(".", "")
        now = now.replace(" ", "-")
        # param_project_name.value = "ReservoirEstimator_%s" % now
        #param_project_name.value = "FWRE_Tests_Settings_Log"

        # Result Folder - Path to a folder for saving the scripts CSV and graph outputs
        param_output_folder = arcpy.Parameter(
            displayName="Output Folder",
            name="out_output_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        param_basins = arcpy.Parameter(
            displayName="BasinATLAS (lev03 recommended)",
            name="in_basins",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input"
        )

        # Just for development, set defaults for required parameters
        #param_dem.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\dem_con\\af_con_3s_rect.tif"
        #param_dem.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\SRTMGL3\\SRTMGL3_WORLD.tif"
        #param_dam_positions.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\v2_0\\v2_0.gdb"
        #param_output_folder.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\v2_0\\Tool_run"
        #param_river_segment.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\RiverATLAS_v10.gdb\\RiverATLAS_v10"
        #param_basins.value ="C:\\Users\\leobu\\Documents\\GIS\M.Sc.Projekt\\Daten\\BasinATLAS_v10.gdb\\BasinATLAS_v10_lev03"

        params = [param_dam_positions, param_dem, param_river_segment, param_height_step, param_upriver_segments,
                  param_dam_id_list, param_project_name, param_output_folder,  param_basins,
                  param_keep_tmp_files] #param_reservoir_overlapping,
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        arcpy.env.overwriteOutput = True
        arcpy.env.parallelProcessingFactor = "100%"

        # Get current extent and workspaces to reset after the program is done
        start_workspace = arcpy.env.workspace
        start_scratch = arcpy.env.scratchWorkspace
        start_extent = arcpy.env.extent

        sum_area = 0    # My Temp additions (JB)
        sum_vol  = 0    # My Temp additions (JB)

        try:
            # Read and verify user inputs
            b_user_input = self.read_user_inputs(parameters)
            if b_user_input:
                self.logger.log_msg("All user inputs seem to be okay")
                b_create_folders = self.create_folder_structure()
                if b_create_folders:
                    self.logger.log_msg("Created output folders and logfile")

                    # Write settings used to logfile
                    self.write_settings_to_log()

                    # Sum important parameters
                    sum_area = 0
                    sum_vol = 0

                    failed_to_process = 0
                    failed_dam_ids = []
                    count_dams = 0
                    basin_list = self.intersect_basins()

                    for basin_id in basin_list:
                        self.current_basin_id = basin_id
                        try:
                            self.set_basin_variables(basin_id)
                        except Exception:
                            self.logger.log_msg("Basin is not in the DEM raster and will be skipped", "Warning")
                            continue
                        # discharge to the FHReD
                        self.build_dam_dict()

                        # river segment dict with links to upriver segments
                        self.build_river_seg_dict()
                        # For every ID in the ID-List
                        for idx, dam_id in enumerate(self.id_list_basin):
                            try:

                                arcpy.env.extent = "MAXOF"
                                count_dams += 1
                                self.logger.log_msg("Progress[%i/%i] # current dam %i" % (count_dams, len(self.id_list) , dam_id))
                                #ram_info = psutil.virtual_memory()
                                #self.loger.log_msg("Of % RAM % MB \% used" % (ram_info[0] / (1024.0**2),ram_info[2]))
                                self.logger.log_msg("Dams failed to process: %i" % failed_to_process)
                                self.logger.log_msg("Failed dam IDs: " + str(failed_dam_ids))

                                # Get watershed polygon and line
                                watershed_poly, watershed_line = self.get_watershed(dam_id)
                                # Get profile of watershed line
                                profile = self.get_profile_as_dict(watershed_line)
                                profile = self.clean_profile_2(profile)
                                #self.export_as_graph(profile["dists"], profile["heights"], "profile_%i" % dam_id,
                                #                     "Full Watershed Profiles\\Original")

                                # Get centered profile
                                profile_centered = self.center_profile(profile)
                                self.export_as_graph(profile_centered["dists"], profile_centered["heights"],
                                                     "profile_%i" % dam_id, "Centered Watersehd Profiles")

                                # Get full info for all possible heights of current dam
                                self.logger.log_msg("Getting dam info for all heights")
                                dam_info = self.get_full_dam_info(dam_id, watershed_poly, profile_centered)

                                # Write dam info to CSV
                                self.logger.log_msg("Exporting results as CSV")
                                #self.save_csv(dam_info)

                                # Save best reservoir as raster and polygon
                                best_dam = self.get_best_dam_2(dam_info)
                                self.save_reservoir_at_h_as_raster_and_poly(watershed_poly, best_dam)

                                # Save best dam as graph as well
                                best_height = best_dam["abs_height"]
                                self.save_best_dam_as_graph(profile_centered, dam_info, best_height)
                                self.logger.log_msg("Exporting selected dam as CSV")
                                #self.save_selected_dam_csv(best_dam)
                                self.save_csv2(dam_info, best_dam)

                                sum_area += best_dam["area"]
                                sum_vol += best_dam["volume"]


                            except Exception:
                                self.logger.log_msg("An Error occurred on dam %i: %s. It will be skipped."
                                                 % (dam_id, str(traceback.format_exc())), "Error")
                                failed_to_process += 1
                                failed_dam_ids.append(dam_id)
                                arcpy.env.extent = "MAXOF"
                            self.delete_features([self.dem_watershed_clip],silent=True)
                        self.delete_features([self.dam_pos_basin,self.dem_basin, self.river_seg_basin],silent=True)
                else:
                    self.logger.log_msg("Output folders could not be created, the tool will terminate now")
            else:
                self.logger.log_msg("User inputs could not be verified, the tool will terminate now")
            self.merge_results()
        finally:
            self.log_summary(sum_area, sum_vol)
            self.logger.log_msg("Resetting workspace environments and extent")
            arcpy.env.scratchWorkspace = start_scratch
            arcpy.env.workspace = start_workspace
            arcpy.env.extent = start_extent

        return

    ####################################################################################################################
    def set_basin_variables(self,basin_id):

        if basin_id == "No Basins given":

            self.dem_basin = self.dem
            self.river_seg_basin = self.river_seg
            self.dam_pos_basin = self.dam_positions

        else:
            self.logger.log_msg("Perparations for basin: {}".format(int(basin_id)))
            basin_lyr = arcpy.env.workspace + "\\basin_lyr"  # was in memory
            arcpy.MakeFeatureLayer_management(self.basins, basin_lyr)

            # Select the basin by id
            arcpy.SelectLayerByAttribute_management(basin_lyr, "NEW_SELECTION", "%s = %i" % (self.basin_field, basin_id))
            dem_basin_clip = arcpy.env.scratchWorkspace + "\\dem_basin_clip_%i.tif" % basin_id
            arcpy.Clip_management(self.dem,"#",dem_basin_clip,basin_lyr,"","ClippingGeometry")
            # First Fill DEM to avoid sinks    # added by Juergen 27.02.20
            self.dem_basin = dem_basin_clip
            #self.dem_basin = arcpy.env.scratchWorkspace + "\\dem_basin_%i.tif" % basin_id
            self.river_seg_basin = arcpy.env.workspace + "\\river_seg_basin_%i" % basin_id
            arcpy.Clip_analysis(self.river_seg,basin_lyr,self.river_seg_basin)
            self.dam_pos_basin = arcpy.env.workspace + "\\dam_pos_basin_%i" % basin_id
            arcpy.Clip_analysis(self.dam_positions,basin_lyr,self.dam_pos_basin)
        cursor = arcpy.da.SearchCursor(self.dam_pos_basin,self.id_field)
        self.id_list_basin = sorted({row[0] for row in cursor}) # unique DAM_IDs in the current basin

    def intersect_basins(self):
        self.logger.log_msg("Copies and snaps the dams to the river segments and builds the basin list")
        dam_lyr = arcpy.env.workspace + "\\dam_lyr"
        arcpy.MakeFeatureLayer_management(self.dam_positions,dam_lyr)
        selection_str = ""
        for dam_id in self.id_list:
            selection_str += " ," + str(dam_id)
        #self.logger.log_msg( "%s IN (%s) " % (self.id_field,selection_str[2:]))
        arcpy.SelectLayerByAttribute_management(dam_lyr, "NEW_SELECTION", "%s IN (%s) " % (self.id_field,selection_str[2:]))
        dam_pos = self.out_gdb + "\\" + self.dam_layer_name + "_dams"
        arcpy.CopyFeatures_management(dam_lyr, dam_pos)
        self.dam_positions = dam_pos
        #arcpy.Snap_edit(self.dam_positions,[[self.river_seg,"EDGE",self.initial_dam_snap_range]])
        self.AddField_helper(self.dam_positions, self.out_hyid_field, "LONG")

        if not self.basins == "None":
            dams_isect_basins = arcpy.env.workspace + "\\dams_isect_basins"
            arcpy.Intersect_analysis([self.dam_positions,self.basins],dams_isect_basins)
            cursor = arcpy.da.SearchCursor(dams_isect_basins,self.basin_field)
            unique_basin_list = sorted({row[0] for row in cursor})
            self.delete_features([dams_isect_basins],silent=True)
            return unique_basin_list
            #self.logger.log_msg(unique_basin_list)
        else:
            return ["No Basins given"]

    def merge_results(self):
        arcpy.env.workspace = self.project_folder + "\\" + self.project_name + "\\Workspace.gdb"
        self.logger.log_msg("Merging the results to new output datasets")
        all_files = arcpy.ListFeatureClasses()
        #arcpy.AddMessage(all_files)
        if all_files:
            res_polys = []
            #res_rasters = arcpy.ListRasters()
            watershed_polys = []
            watershed_lines = []
            for file in all_files:
                file = str(file)
                if "reservoir_polygon" in file:
                    res_polys.append(arcpy.env.workspace + "\\" + file)
                #elif "watershed_line" in file:
                #    watershed_lines.append(arcpy.env.workspace + "\\" + file)
                elif "watershed_poly_projected" in file:
                    tmp_file = arcpy.env.workspace + "\\" + file
                    # Dissolves the watershed polygon
                    tmp_file_diss = arcpy.env.workspace + "\\tmp_diss_" + file
                    arcpy.Dissolve_management(tmp_file, tmp_file_diss)
                    watershed_polys.append(tmp_file_diss)
            all_res_polys = self.out_gdb + "\\"+ self.dam_layer_name + "_reservoirs"
            #all_res_rasters = self.out_gdb + "\\reservoir_rasters"
            all_watershed_polys = self.out_gdb +  "\\"+ self.dam_layer_name + "_watersheds"
            #all_watershed_lines = self.out_gdb + "\\watershed_lines_projected"
            #arcpy.AddMessage(str(res_polys) + str(all_res_polys))
            #arcpy.Merge_management(res_rasters, all_res_rasters)
            try:
                arcpy.Merge_management(res_polys, all_res_polys)
                arcpy.Merge_management(watershed_polys, all_watershed_polys)
                #arcpy.Merge_management(watershed_lines, all_watershed_lines)

                basin_level = str(self.basins).split("lev")[1]
                self.AddField_helper(self.dam_positions, self.out_vol_field, "DOUBLE" )
                self.AddField_helper(self.dam_positions, self.out_area_field, "DOUBLE")
                self.AddField_helper(self.dam_positions, self.out_est_key, "SHORT")

                #self.AddField_helper(self.dam_positions, "in_BasinATLAS_lev" + str(basin_level), "LONG")
                dam_pos_lyr = arcpy.env.workspace + "\\dam_pos_lyr"  # was in memory
                arcpy.MakeFeatureLayer_management(self.dam_positions, dam_pos_lyr)
                arcpy.AddJoin_management(dam_pos_lyr, self.id_field, all_res_polys, self.id_field)
                arcpy.CalculateField_management(dam_pos_lyr, self.out_vol_field,"!{}_reservoirs.r_vol_mcm!".format(self.dam_layer_name),"PYTHON3")
                arcpy.CalculateField_management(dam_pos_lyr, self.out_area_field,"!{}_reservoirs.r_area_msm!".format(self.dam_layer_name),"PYTHON3")
                arcpy.CalculateField_management(dam_pos_lyr, self.out_est_key,"!{}_reservoirs.d_est_key!".format(self.dam_layer_name),"PYTHON3")
                arcpy.RemoveJoin_management(dam_pos_lyr, "{}_reservoirs".format(self.dam_layer_name))
                #self.delete_features(dam_pos_lyr,silent=True)
            except Exception:
                self.logger.log_msg("Merging failed! Most likely due to missing output data", "Warning")
            if self.keep_tmp_files == False:
                self.logger.log_msg("Deletes all temporary files and Workspace GDB & folder")
                arcpy.env.workspace = self.project_folder + "\\" + self.project_name + "\\Workspace.gdb"
                all_files = arcpy.ListFeatureClasses()
                arcpy.Delete_management(all_files)
                arcpy.Delete_management(self.project_folder + "\\" + self.project_name + "\\Workspace.gdb")
                arcpy.Delete_management(self.project_folder + "\\" + self.project_name + "\\ScratchWorkspace")

    def build_dam_dict(self):

        self.logger.log_msg("Builds the dam dictionary")


        dams_isect_rivseg = arcpy.env.workspace + "\\dams_isect_rivseg_per_basin"
        #self.logger.log_msg(dams_isect_rivseg)
        #arcpy.Snap_edit(self.dam_pos_basin,[[self.river_seg_basin,"EDGE","500 Meters"]])
        dams_buffer = arcpy.env.workspace + "\\dams_buffer_per_basin"

        arcpy.Buffer_analysis(self.dam_pos_basin,dams_buffer,self.initial_dam_snap_range)
        if self.discharge_field in [f.name for f in arcpy.ListFields(dams_buffer)]:
            arcpy.DeleteField_management(dams_buffer, self.discharge_field)
        arcpy.Intersect_analysis([dams_buffer,self.river_seg],dams_isect_rivseg) #,cluster_tolerance="10 Meters"

        #self.logger.log_msg([f.name for f in arcpy.ListFields(dams_isect_rivseg)])
        cursor = arcpy.da.SearchCursor(dams_isect_rivseg,[self.id_field,self.discharge_field, self.hydroid_field,
                                        self.discharge_pmx, self.discharge_pmn, self.discharge_pyr]) # , self.elevation_field
        #self.logger.log_msg(cursor.fields)
        self.dam_dict = {}
        for row in cursor:
            #self.logger.log_msg("row " + str(row[0]))
            if not row[0] in self.dam_dict:
                self.dam_dict[row[0]] = {self.discharge_field:row[1], self.hydroid_field:row[2], self.discharge_pmx:row[3],
                                            self.discharge_pmn:row[4], self.discharge_pyr:row[5]} #, self.elevation_field:row[3]
            else:
                if row[1] > self.dam_dict[row[0]][self.discharge_field]:
                    self.dam_dict[row[0]] = {self.discharge_field:row[1], self.hydroid_field:row[2], self.discharge_pmx:row[3],
                                                self.discharge_pmn:row[4], self.discharge_pyr:row[5]}
        self.delete_features([dams_isect_rivseg,dams_buffer], silent=True)
        self.AddField_helper(self.dam_positions, self.discharge_pmx, "FLOAT")
        self.AddField_helper(self.dam_positions, self.discharge_pmn, "FLOAT")
        self.AddField_helper(self.dam_positions, self.discharge_pyr, "FLOAT")
        cursor = arcpy.da.UpdateCursor(self.dam_positions,[self.id_field, self.out_hyid_field, self.discharge_pmx, self.discharge_pmn, self.discharge_pyr])
        for row in cursor:
            if row[0] in self.dam_dict:
                row[1] = self.dam_dict[row[0]][self.hydroid_field]
                row[2] = self.dam_dict[row[0]][self.discharge_pmx]
                row[3] = self.dam_dict[row[0]][self.discharge_pmn]
                row[4] = self.dam_dict[row[0]][self.discharge_pyr]
                cursor.updateRow(row)
        #self.logger.log_msg(self.dam_dict)
        if not self.dam_dict:
            self.logger.log_msg("Intersect of dams and river segments failed, check if they span the same area","Error")

    def build_river_seg_dict(self):
        self.logger.log_msg("Builds the river segment dictionary")
        try:

            self.river_seg_dict = {}
            cursor = arcpy.da.SearchCursor(self.river_seg_basin,[self.hydroid_field, self.nextdown_field])
            for row in cursor:
                self.river_seg_dict[row[0]]={self.nextdown_field:row[1],"NEXT_UP":[],"DAM":None}
            for key in self.river_seg_dict.keys():
                next_down = self.river_seg_dict[key][self.nextdown_field]
                if next_down and next_down in self.river_seg_dict:
                    self.river_seg_dict[next_down]["NEXT_UP"].append(key)
            for key in self.dam_dict.keys():
                hydroid = self.dam_dict[key][self.hydroid_field]
                self.river_seg_dict[hydroid]["DAM"] = key
        except Exception:
            self.logger.log_msg("An Error occured building the river segment dictionary, this can be caused by an overused memory and fixed by restarting the PC","Error")
            self.logger.log_msg("The given Error was: " + str(traceback.format_exc()), "Error")

    def log_summary(self, sum_area, sum_volume):
        self.logger.log_msg("################SUMMARY################")
        self.logger.log_msg("Area Sum [km^2]: %.4f" % sum_area)
        self.logger.log_msg("Volume Sum [km^3]: %.4f" % sum_volume)
        self.logger.log_msg("#######################################")

    def save_best_dam_as_graph(self, profile, dam_info, height):
        # Get at height index of best dam
        best_dam_info = {}
        dam_id = dam_info["dam_id"]
        for info in dam_info["per_height"]:
            if info["abs_height"] == height:
                best_dam_info = info

        abs_h = best_dam_info["abs_height"]
        cut_indices = self.get_cut_indices_for_height_line(abs_h, profile)
        segments = self.get_dam_line_segments(profile, cut_indices, abs_h)

        # TODO: find better solution
        dam_info["best_dam"] = best_dam_info

        # Cut profile to only show relevant information for dam line
        x_left = cut_indices[0][0]
        x_right = cut_indices[-1][1]
        dists = profile["dists"][x_left: x_right+1]
        heights = profile["heights"][x_left: x_right+1]

        self.export_as_graph_with_dam_line(dists, heights, segments,
                                           abs_h, "dam_%i_%i" % (dam_id, abs_h), dam_info,
                                           "Results")

    def get_best_dam_2(self, dam_info):
        """
        A more structured approach for dam selection:
        1. Dams are filtered out by following critria:
            -Max height x meters
            -Max width x meters (hard to validate this threshold)
            -Max reservoir area x square kilometers
            -Max reservoir volume x cubic kilometers
            -Max dam line interruption value

        2. If drop height <= 300m:
            -Choose dam closest to drop height, if a dam with a height difference of less than 50% is available
           Else:
            -Choose dam with max vol/cost from remaining dams

        3. If no dams remain after filtering, filter again only by dam line interruption and  max width + 30%
        """

        # create list of dams that fulfill all criteria
        self.logger.log_msg("Filtering out implausible dams")
        #self.logger.log_msg("xxxdam info: " +str(dam_info["per_height"]))
        possible_dams = []
        additional_comment = ""
        for dam in dam_info["per_height"]:
            if (dam["rel_height"] <= self.max_dam_height
                    and dam["rel_height"] >= self.min_dam_height
                    and dam["total_length"] <= self.max_dam_width
                    and dam["area"] <= self.max_reservoir_area
                    and dam["volume"] <= self.max_reservoir_volume
                    and dam["share_above_dam_line"] <= self.max_dam_interruption)\
                    and dam["total_length"] * dam["rel_height"] <= self.max_width_x_height:
                possible_dams.append(dam)
        self.logger.log_msg("%i possible dams out of %i dams remain for selection"
                         % (len(possible_dams), len(dam_info["per_height"])))
        #self.logger.log_msg("xxxxxxxpos dams: " + str( possible_dams))
        #self.logger.log_msg(possible_dams)
        if len(possible_dams) < 1:
            additional_comment = " (Soft filters used)"
            self.logger.log_msg("Using softer filters now to get a result")
            for dam in dam_info["per_height"]:
                if dam["share_above_dam_line"] <= self.max_dam_interruption \
                        and dam["total_length"] <= 1.3 * self.max_dam_width:
                    possible_dams.append(dam)
            self.logger.log_msg("%i possible dams out of %i dams remain for selection after soft filtering"
                         % (len(possible_dams), len(dam_info["per_height"])))
        if len(possible_dams) < 1:
            raise Exception("Using soft filters also resulted in 0 usable dams, so no dam can be selected")
        self.logger.log_msg("Checking if given or estimated drop height can be used for dam selection")
        est_drop_height = dam_info["estimated_drop_height"]
        est_drop_height_usable = False
        if est_drop_height and self.max_dam_height >= est_drop_height > 0:
            est_drop_height_usable = True
            est_drop_height = int(est_drop_height)
        #self.logger.log_msg("drop_usable: " + str(est_drop_height_usable) + " height " + str(est_drop_height))
        # Add 3m freeboard as safety factor to drop height
        # est_drop_height += 3
        given_height = dam_info["given_dam_height"]
        given_height_usable = False

        #self.logger.log_msg("xxxxxx given height: " + str(given_height))
        if given_height and self.max_dam_height >= given_height > 0:
            given_height_usable = True
            given_height = int(given_height)
        elif given_height:
            self.logger.log_msg("Given dam height of %i is out of max. dam height bounds of %i" % (given_height, self.max_dam_height))

        if given_height_usable:
            self.logger.log_msg("Given dam height of %i meters is within acceptable bounds and will be tested for dam selection"
                             % int(given_height))
            # Select dam closest to drop height
            closest_to_given_height = None
            closest_to_given_height_idx = None
            #self.logger.log_msg("xxxx"+str(possible_dams)) #########################################
            for i in range(0, len(possible_dams)):
                if closest_to_given_height is None or abs(given_height - possible_dams[i]["rel_height"]) < closest_to_given_height:
                    closest_to_given_height = abs(given_height - possible_dams[i]["rel_height"])
                    closest_to_given_height_idx = i

            # if the difference to the closest dam is more than 50% of the calculated drop height
            #self.logger.log_msg("xxxx" + str(closest_to_given_height) + " xxxx" + str(given_height)) ############
            try:
                delta_h = float(closest_to_given_height) / float(given_height)
            except Exception as e:
                self.logger.log_msg("No suitable dam was found","Error")
            self.logger.log_msg("The closest dam deviates %.2f percent in height (%i/%i)." % (delta_h*100,int(possible_dams[closest_to_given_height_idx]["rel_height"]),given_height))
            if delta_h > 0.5 and given_height > self.min_dam_height: # second part only filters too big dams nor too small ones
                self.logger.log_msg("There is no dam close enough to the given dam height. "
                                 "A dam will be selected by a capacity derived dam height")
            else:
                message = "Selected a dam of height %i meters for a given dam height of %i meters" \
                    % (int(possible_dams[closest_to_given_height_idx]["rel_height"]), given_height)
                decision = 1
                self.logger.log_msg(message)
                possible_dams[closest_to_given_height_idx]["comment"] = message + additional_comment
                possible_dams[closest_to_given_height_idx]["decision"] = decision
                possible_dams[closest_to_given_height_idx]["deviation"] = abs(possible_dams[closest_to_given_height_idx]["rel_height"]/given_height-1)
                return possible_dams[closest_to_given_height_idx]

        if est_drop_height_usable:
            self.logger.log_msg("Estimated drop height of %i meters is within acceptable bounds and will be tested for dam selection"
                             % int(est_drop_height))
            # Select dam closest to drop height
            closest_to_drop_height = None
            closest_to_drop_height_idx = None
            #self.logger.log_msg("xxxx"+str(possible_dams)) #########################################
            for i in range(0, len(possible_dams)):
                if closest_to_drop_height is None or abs(est_drop_height - possible_dams[i]["rel_height"]) < closest_to_drop_height:
                    closest_to_drop_height = abs(est_drop_height - possible_dams[i]["rel_height"])
                    closest_to_drop_height_idx = i

            # if the difference to the closest dam is more than 50% of the calculated drop height
            #self.logger.log_msg("xxxx" + str(closest_to_drop_height) + " xxxx" + str(est_drop_height)) ############
            error_calc_h = False
            try:
                delta_h = float(closest_to_drop_height) / float(est_drop_height)
                self.logger.log_msg("The closest dam deviates %.2f percent in height (%i/%i)." % (delta_h*100,int(possible_dams[closest_to_drop_height_idx]["rel_height"]),est_drop_height))
            except Exception as e:
                self.logger.log_msg("No dam was found using the estimated drop height")
                error_calc_h = True


            if error_calc_h or (delta_h > 0.5 and est_drop_height > self.min_dam_height): # second part only filters too big dams nor too small ones
                if not error_calc_h:
                    self.logger.log_msg("There is no dam close enough to the estimated drop height. ")
                self.logger.log_msg("A dam will be selected by the maximum of reservoir volume per cost")
                # Get dam with max vol/cost
                try:
                    max_vol_cost = None
                    max_vol_cost_idx = None
                    for i in range(0, len(possible_dams)):
                        modified_cost = np.power(possible_dams[i]["cost"], 1.1)  # experimental
                        vol_cost = float(possible_dams[i]["volume"]) / float(modified_cost)
                        if max_vol_cost is None or vol_cost > max_vol_cost:
                            max_vol_cost = vol_cost
                            max_vol_cost_idx = i
                    message = "Selected dam based on max volume per cost, because there was no dam close enough" \
                              " to the estimated drop height available. Resulting dam height is %i meters" \
                              % int(possible_dams[max_vol_cost_idx]["rel_height"])
                    decision = 3
                    self.logger.log_msg(message)
                    possible_dams[max_vol_cost_idx]["comment"] = message + additional_comment
                    possible_dams[max_vol_cost_idx]["decision"] = decision
                    possible_dams[max_vol_cost_idx]["deviation"] = None
                    return possible_dams[max_vol_cost_idx]
                except Exception:
                    raise Exception("Calucaltion based on max volume per cost failed, so no reservoir could be calculated!")
            else:
                message = "Selected a dam of height %i meters for an estimated drop height of %i meters" \
                    % (int(possible_dams[closest_to_drop_height_idx]["rel_height"]), est_drop_height)
                decision = 2
                self.logger.log_msg(message)
                possible_dams[closest_to_drop_height_idx]["comment"] = message + additional_comment
                possible_dams[closest_to_drop_height_idx]["decision"] = decision
                possible_dams[closest_to_drop_height_idx]["deviation"] = abs(possible_dams[closest_to_drop_height_idx]["rel_height"]/est_drop_height-1)
                return possible_dams[closest_to_drop_height_idx]
        else:
            self.logger.log_msg("Drop height of %i meters is not usable, dam selection will be based on"
                             " the maximum of reservoir volume per estimated cost" % int(est_drop_height))
            # Get dam with max vol/cost
            try:
                max_vol_cost = None
                max_vol_cost_idx = None
                for i in range(0, len(possible_dams)):
                    modified_cost = np.power(possible_dams[i]["cost"], 1.1)  # experimental
                    vol_cost = float(possible_dams[i]["volume"]) / float(modified_cost)
                    if max_vol_cost is None or vol_cost > max_vol_cost:
                        max_vol_cost = vol_cost
                        max_vol_cost_idx = i
                message = "Selected dam based on max volume per cost. Resulting dam height is %i meters" \
                          % int(possible_dams[max_vol_cost_idx]["rel_height"])
                decision = 3
                self.logger.log_msg(message)
                possible_dams[max_vol_cost_idx]["comment"] = message + additional_comment
                possible_dams[max_vol_cost_idx]["decision"] = decision
                possible_dams[max_vol_cost_idx]["deviation"] = None
                return possible_dams[max_vol_cost_idx]
            except Exception:
                raise Exception("Calucaltion based on max volume per cost failed, so no reservoir could be calculated!")

    def get_best_dam(self, dam_info):

        """
        This is the first decision function, as discussed with Berlekamp friday, the 14.06.2019
        It uses drop height and width x height of the worlds largest dams as the main criteria.
        If no dam fulfills all criteria, the last one segment dam before multi segment dams will be used
        """
        # Check, if drop height is reasonable
        est_drop_height = dam_info["estimated_drop_height"]
        if self.min_dam_height <= est_drop_height <= self.max_dam_height:
            self.logger.log_msg("Drop Height of %i m. seems reasonable" % int(est_drop_height))
            # Find dam closest to calculated drop height
            delta_drop_height_list = []
            for dam_at_h in dam_info["per_height"]:
                delta_drop_height_list.append(abs(est_drop_height - dam_at_h["rel_height"]))

            # Find index of minimal delta to drop height
            min_delta_idx = np.argmin(delta_drop_height_list)
            self.logger.log_msg("Dam closest to drop height is at %i m."
                             % int(dam_info["per_height"][min_delta_idx]["rel_height"]))

            # Check width x height of the selected dam
            w_x_h = dam_info["per_height"][min_delta_idx]["total_length"] * \
                dam_info["per_height"][min_delta_idx]["rel_height"]

            # If width x height ist not withing bounds
            if w_x_h > self.max_width_x_height:
                self.logger.log_msg("Width x height of selected dam is above the average of the worlds largest dams."
                                 " A dam not exceeding this threshold will be selected instead")
                for idx in range(min_delta_idx, 0, -1):
                    w_x_h = dam_info["per_height"][idx]["total_length"] * \
                        dam_info["per_height"][idx]["rel_height"]
                    if w_x_h <= self.max_width_x_height:
                        self.logger.log_msg("The first dam not exceeding the maximum width x height"
                                         " criteria was found at %i meters"
                                         % int(dam_info["per_height"][idx]["rel_height"]))
                        dam_info["per_height"][idx]["comment"] = "Selected width x heigth criteria," \
                                                            " because drop height selection resulted in a too huge dam"
                        return dam_info["per_height"][idx]
            else:
                dam_info["per_height"][min_delta_idx]["comment"] = "Selected by drop height calculation"
                return dam_info["per_height"][min_delta_idx]
        else:
            self.logger.log_msg("Drop Height of %i m. is not within acceptable bounds, "
                             "or no dam within the width x height criteria was found."
                             " Will choose dam based on segment count now" % int(est_drop_height))

            # naive decision: find last one segment dam before first multi segment dam
            # find first one segment dam first to start from
            first_one_segment_dam_idx = 0
            for idx, dam in enumerate(dam_info["per_height"]):
                if dam["segment_count"] == 1:
                    first_one_segment_dam_idx = idx
                    break

            best_idx = -1
            for idx in range(first_one_segment_dam_idx, len(dam_info["per_height"])):
                dam = dam_info["per_height"][idx]
                if dam["segment_count"] > 1 or idx == len(dam_info["per_height"])-1:
                    best_idx = idx-1
                    break

            # check, if best index is a valid index, else return 0
            if best_idx < 0 or best_idx > len(dam_info["per_height"])-1:
                self.logger.log_msg("No valid best dam could be selected", "Warning")
                dam_info["per_height"][0]["comment"] = "No valid dam could be selected"
                return dam_info["per_height"][0]
            else:
                # Check, if retrieved dam is within max width x height bounds
                # Check width x height of the selected dam
                w_x_h = dam_info["per_height"][best_idx]["total_length"] * \
                        dam_info["per_height"][best_idx]["rel_height"]
                # If width x height ist not withing bounds
                if w_x_h > self.max_width_x_height:
                    self.logger.log_msg("Width x height of selected one-segment dam is above the average "
                            "of the worlds largest dams. A dam not exceeding this threshold will be selected instead")
                    for idx in range(best_idx, 0, -1):
                        w_x_h = dam_info["per_height"][idx]["total_length"] * \
                                dam_info["per_height"][idx]["rel_height"]
                        if w_x_h <= self.max_width_x_height:
                            self.logger.log_msg("The first dam not exceeding the maximum width x height"
                                             " criteria was found at %i meters"
                                             % int(dam_info["per_height"][idx]["rel_height"]))
                            dam_info["per_height"][idx]["comment"] = "Selected width x heigth criteria," \
                                                        " because drop height selection resulted in a too huge dam"
                            return dam_info["per_height"][idx]
                else:
                    self.logger.log_msg("Best Height for dam %i is %i" % (dam_info["dam_id"],
                                                                       dam_info["per_height"][best_idx]["rel_height"]))
                    dam_info["per_height"][best_idx]["comment"] = "Selected largest one-segment dam before" \
                                                                  " multi segment dams occurred"
                    return dam_info["per_height"][best_idx]

    def get_full_dam_info(self, dam_id, watershed_polygon, profile):
        """Collects and returns information for all possible dam heights for a given dam
           This includes: relative height, absolute height, area, volume, cost estimation, dam length,
                          number of dam segments, valley factor, drop height from discharge and capacity,
                          theoretical capacity at every height.
                          Results are returned in a easily expandable dictionary
        """

        # Get the raster values within the watershed polygon as a numpy array
        #self.logger.log_msg(arcpy.env.workspace + "\\" + str(watershed_polygon))
        '''self.logger.log_msg("Clipping DEM for dam %i" % dam_id)
        dem_clip = arcpy.env.scratchWorkspace + "\\dem_clip"
        arcpy.Clip_management(self.dem_fil, "#", dem_clip, arcpy.env.workspace + "\\" + str(watershed_polygon), clipping_geometry="ClippingGeometry")
        dem_array = arcpy.RasterToNumPyArray(dem_clip, nodata_to_value=20000)'''

        dem_array = arcpy.RasterToNumPyArray(self.dem_watershed_clip, nodata_to_value=20000)
        # Get min and max height of dem clip as boundaries
        heights = profile["heights"]
        dists = profile["dists"]
        min_h = int(np.min(dem_array)) # floor from actual filled dem in the watersehd
        max_h = int(np.max(heights)) + 1 # roof from the height profile
        # self.logger.log_msg("min_h %i max_h %i" % (min_h, max_h)) ###################################
        # Limit maximum dam height to max dam height meters
        if max_h - min_h > self.max_dam_height:
            max_h = min_h + self.max_dam_height

        # Estimate relative drop height from capacity and discharge
        est_drop_height, given_capacity, given_dam_height = self.estimate_drop_height(dam_id)
        #self.logger.log_msg("xxxxxx" + str(given_dam_height))

        results = {
                        "dam_id": dam_id,
                        "estimated_drop_height": est_drop_height,
                        "given_capacity": given_capacity,
                        "given_dam_height": given_dam_height,
                        "per_height": []
                   }

        # Get area calculation correction based on latitude of the dam position
        lat = self.get_dam_lat(dam_id)
        self.logger.log_msg("Latitude: " + str(lat))

        # Abweitung bezogen auf ein Längengrad
        abweitung = (40030.0 / 360.0) * np.cos(np.radians(np.absolute(lat)))

        lat_corr_fac = abweitung / 111.319 # km per 1 degree lat at the equator

        self.logger.log_msg("Latitude correction: " + str(lat_corr_fac))

        # For every height between max and min with height_step step size
        self.logger.log_msg("Calculating all possible dam heights and their reservoirs")
        height_step_counter = 0

        for h in range(min_h + self.height_step, max_h + self.height_step, self.height_step):
            #self.logger.log_msg("Dam %i height %i meters" % (dam_id, h))

            rel_h = h-min_h  # Relative height of the dam

            area = self.get_area_at_h(dem_array, h, lat_corr_fac)

            vol = self.get_volume_at_h(dem_array, h, lat_corr_fac)

            electric_output = self.estimate_electric_output(rel_h, dam_id)

            cost, segment_count, total_length, fragmentation, segments = self.calculate_cost(profile, h, dam_id)
            if cost == -1:
                continue # if no intersection is found, skip this height
            if not cost:
                break  # Stop at this dam height

            # Get percent of points above dam height in profile between outmost dam line segements
            share_above_dam_line = self.get_share_above_dam_line(profile, h)

            # Add to results
            results["per_height"].append({
                "dam_id": dam_id,  # redundant but useful
                "estimated_drop_height": est_drop_height,  # redundant but useful
                "given_capacity": given_capacity,  # redundant but useful
                "abs_height": h,
                "rel_height": rel_h,
                "area": area,
                "volume": vol,
                "electric_output": electric_output,
                "cost": cost,
                "vol_cost": vol/cost,
                "segment_count": segment_count,
                "segments": segments,
                "total_length": total_length,
                "fragmentation": fragmentation,
                "share_above_dam_line": share_above_dam_line
            })

            #self.logger.log_msg("results" + str(results["best_dam"])) #################################################

            # Export graph of dam with additional information
            '''height_step_counter += self.height_step
            if height_step_counter >= 10:
                height_step_counter = 0
                self.export_as_graph_with_dam_line(profile["dists"], profile["heights"], segments,
                                                   h, "dam_%i_%i" % (dam_id, h), results,
                                                   "All 10m Heights Dam Lines\\%i" % dam_id)'''

        # Delete no longer needed data
        #self.delete_features([dem_watershed_clip,watershed_raster], True) ####################################gggggggg

        return results

    def get_area_at_h(self, data, h, lat_corr_fac):
        """Calculates the area at a given height in million square meters"""
        # Copy array, to not manipulate original data
        data_copy = data.copy()

        # Set all  values > h to zero
        data_copy[data_copy > h] = 0

        # Count non-zero cells
        non_zero_elements_count = np.count_nonzero(data_copy)

        # Calculate area in km^2
        area = (non_zero_elements_count * self.dem_cell_size * self.dem_cell_size * lat_corr_fac) / 1000000.0
        return area

    def get_volume_at_h(self, data, h, lat_corr_fac):
        """Calculates the volume at a given height in cubic kilometers"""
        # Copy array, to not manipulate original data
        data_copy = data.copy()

        # Set all  values > h to zero
        data_copy[data_copy > h] = 0

        # Count non-zero cells
        non_zero_elements_count = np.count_nonzero(data_copy)

        at_height = non_zero_elements_count * h
        below_height = data_copy.sum()
        delta_height = at_height - below_height
        patch_size = np.int64(self.dem_cell_size * self.dem_cell_size * lat_corr_fac)
        vol = np.float64(delta_height * patch_size)
        vol_ckm = vol / 1000000000.0  # cubic kilometers

        return vol_ckm

    def estimate_electric_output(self, height, dam_id):
        """Calculates theoretical output from discharge and relative height"""
        q = self.dam_dict[dam_id][self.discharge_field]

        # calculate electric power at given height
        p_el = 8.5 * q * height
        p_el = p_el / 1000.0  # Results in megawatt

        return p_el

    def estimate_drop_height(self, dam_id):
        """
            Calculates theoretical output from discharge and relative height.
            Also returns the given capacity of the dam.
        """

        capacity, given_dam_height = self.get_dam_cap_height(dam_id)
        #self.logger.log_msg("cap: " + str(capacity) + " given_height " + str(given_dam_height))
        try:
            p_el = float(capacity) * 1000 # Change here to capacity * 0.4 to use avg dis data
        except TypeError:
            self.logger.log_msg("No capacity was provided, setting drop height to -1 meters")
            return -1, -1, -1

        # read discharge from dam info dictionary
        #self.logger.log_msg("xxxx p_el: " + str(p_el))
        q = self.dam_dict[dam_id][self.discharge_field]
        #self.logger.log_msg("xxxx q: " + str(q))
        # calculate drop height
        h = p_el / q / 8.5
        #self.logger.log_msg("xxxx h: " + str(h))
        if h < 1:
            h = 1
        return h, capacity, given_dam_height

    def get_dam_xy(self, dam_id):
        with arcpy.da.SearchCursor(self.dam_pos_basin, [self.id_field, "SHAPE@XY"]) as cursor:
            for row in cursor:
                curr_id = int(row[0])
                dam_pos_xy = row[1]
                if curr_id == dam_id:
                    return dam_pos_xy

    def get_dam_lat(self, dam_id):
        with arcpy.da.SearchCursor(self.dam_pos_basin, [self.id_field, self.lat_field]) as cursor:
            for row in cursor:
                curr_id = int(row[0])
                lat = row[1]
                if curr_id == dam_id:
                    lat = float(str(lat).replace(',', '.'))
                    return lat

    def get_dam_cap_height(self, dam_id):
        with arcpy.da.SearchCursor(self.dam_pos_basin, [self.id_field, self.cap_field, self.height_field]) as cursor:
            for row in cursor:
                if row[0] == dam_id:
                    #self.logger.log_msg(str(row))
                    return row[1] , row[2]
                    """
                curr_id = int(row[0])
                cap = row[1]
                height = row[2]
                if curr_id == dam_id:
                    return cap"""

    def save_reservoir_at_h_as_raster_and_poly(self, watershed, dam):
        height = dam["abs_height"]
        dam_id = dam["dam_id"]
        self.logger.log_msg("Saving reservoir raster at %i" % height)
        #dem_clip = arcpy.env.scratchWorkspace + "\\dem_clip"
        #arcpy.Clip_management(self.dem_fil, "#", dem_clip, arcpy.env.workspace + "\\" + str(watershed), clipping_geometry="ClippingGeometry") # dem_basin
        #dem_clip = self.dem_watershed_clip
        reservoir_raster = arcpy.env.scratchWorkspace + "\\reservoir_raster_without_height_%i.tif" % dam_id
        reservoir_raster_wh = arcpy.env.scratchWorkspace + "\\reservoir_raster_%i.tif" % dam_id
        reservoir_raster_res = arcpy.sa.Raster(self.dem_watershed_clip) < height
        #self.logger.log_msg(height)
        #reservoir_raster_res.save(arcpy.env.scratchWorkspace + "\\test_%i.tif" % dam_id)
        #reservoir_raster_res_maj_filter = arcpy.sa.MajorityFilter(reservoir_raster_res,"FOUR","MAJORITY")
        reservoir_raster_res_boundary_clean = arcpy.sa.BoundaryClean(reservoir_raster_res,"NO_SORT","ONE_WAY")
        reservoir_raster_res_with_no_data = arcpy.sa.SetNull(reservoir_raster_res_boundary_clean == 0, reservoir_raster_res_boundary_clean)
        reservoir_raster_res_with_no_data.save(reservoir_raster)

        reservoir_raster_res_with_height = arcpy.sa.Raster(reservoir_raster) * int(height)
        reservoir_raster_res_with_height.save(reservoir_raster_wh)

        # Smoothing result polygon

        # Resample to four times the original resolution
        cell_size_result = arcpy.GetRasterProperties_management(reservoir_raster_wh, "CELLSIZEX")
        cell_size = cell_size_result.getOutput(0).replace(",", ".")

        self.logger.log_msg("Resampling reservoir raster to higher resolution")
        try:
            resampled = arcpy.env.scratchWorkspace + "\\resampled.tif"
            arcpy.Resample_management(reservoir_raster_wh, resampled, float(cell_size)/4.0, "NEAREST") # former NEAREST
        except Exception as e:
            self.logger.log_msg("Error trying to resample the raster, skipping the resample", "Warning")
            resampled = reservoir_raster_wh

        # Lowpass filter
        self.logger.log_msg("Executing lowpass filter")
        low_pass_res = arcpy.env.scratchWorkspace + "\\low_pass.tif"
        #filter_out = arcpy.sa.Filter(resampled, "LOW", "DATA")
        filter_out = arcpy.sa.FocalStatistics(resampled) # 3x3 rect, mean
        #filter_out = arcpy.sa.MajorityFilter(resampled,"")
        # Save the output
        filter_out.save(low_pass_res)

        # From float to int
        self.logger.log_msg("Converting from float to integer representation")
        low_pass_int_res = arcpy.env.scratchWorkspace + "\\low_pass_int.tif"
        out_int = arcpy.sa.Int(low_pass_res)

        try:
            # Save the output
            out_int.save(low_pass_int_res)
        except Exception as e:
            self.logger.log_msg("Error saving the output, most likely a memory error. Try to rerun the script in an empty ArcMap file.\n %s"
                                    % str(traceback.format_exc()),"Error")

        # Raster to poly
        self.logger.log_msg("Saving reservoir polygon at %i" % height)
        reservoir_poly_raw = "%s\\reservoir_raw_polygon_%i" % (arcpy.env.workspace, dam_id)
        arcpy.RasterToPolygon_conversion(low_pass_int_res, reservoir_poly_raw)

        # Smooth polygon with PEAK and 500m tolerance
        self.logger.log_msg("Smoothing resulting polygon")
        reservoir_poly = "%s\\reservoir_polygon_%i" % (arcpy.env.workspace, dam_id)
        arcpy.cartography.SmoothPolygon(reservoir_poly_raw, reservoir_poly, "PAEK", "500 Meters")

        self.delete_features([reservoir_raster, resampled, low_pass_res, reservoir_poly_raw, low_pass_int_res, reservoir_raster_res_boundary_clean, reservoir_raster_wh],silent=True)



        # Add fields and information to polygon
        self.logger.log_msg("Writing dam info to polygon")
        self.AddField_helper(reservoir_poly, self.id_field, "LONG")
        #basin_level = str(self.basins).split("lev")[1]
        #self.AddField_helper(reservoir_poly, "Dam_in_BasinATLAS_lev" + str(basin_level), "DOUBLE")
        self.AddField_helper(reservoir_poly, "r_sf_asl_m", "LONG")
        self.AddField_helper(reservoir_poly, "d_hgt_m", "LONG")
        self.AddField_helper(reservoir_poly, "r_area_msm", "DOUBLE")
        #self.AddField_helper(reservoir_poly, "Reservoir_area_msm", "DOUBLE")
        self.AddField_helper(reservoir_poly, "r_vol_mcm", "DOUBLE")
        self.AddField_helper(reservoir_poly, "d_len_m", "LONG")
        self.AddField_helper(reservoir_poly, "d_seg_no", "LONG")
        self.AddField_helper(reservoir_poly, "e_cap_mw", "DOUBLE")
        self.AddField_helper(reservoir_poly, "e_drop_m", "DOUBLE")
        self.AddField_helper(reservoir_poly, "e_cost_maud", "DOUBLE")
        #self.AddField_helper(reservoir_poly, "Estimated_volume_per_cost", "DOUBLE")
        self.AddField_helper(reservoir_poly, "r_est_com", "TEXT")
        self.AddField_helper(reservoir_poly, "d_est_key", "SHORT")


        update_cursor = arcpy.UpdateCursor(reservoir_poly)
        for update_row in update_cursor:
            update_row.setValue(self.id_field, int(dam["dam_id"]))
            #update_row.setValue("Dam_in_BasinATLAS_lev" + str(basin_level), float(self.current_basin_id))
            update_row.setValue("r_sf_asl_m", int(dam["abs_height"]))
            update_row.setValue("d_hgt_m", int(dam["rel_height"]))
            update_row.setValue("r_area_msm", float(dam["area"]))
            update_row.setValue("r_vol_mcm", float(dam["volume"]*1000))
            update_row.setValue("d_len_m", int(dam["total_length"]))
            update_row.setValue("d_seg_no", int(dam["segment_count"]))
            update_row.setValue("e_cap_mw", float(dam["electric_output"]))
            update_row.setValue("e_drop_m", float(dam["estimated_drop_height"]))
            update_row.setValue("e_cost_maud", float(dam["cost"]))
            #update_row.setValue("Estimated_volume_per_cost", float(dam["vol_cost"]))
            update_row.setValue("r_est_com", str(dam["comment"]))
            update_row.setValue("d_est_key", str(dam["decision"]))

            update_cursor.updateRow(update_row)

        #self.delete_features([dem_clip],silent=True)

    def calculate_cost(self, profile, height, dam_id):
        #self.logger.log_msg("xxxxxxxheight: " + str(height))
        #self.logger.log_msg("xxxxxxprofile" + str(profile))
        # If the height given is higher than one edge of the profile data, return False
        if height >= profile["heights"][0] or height >= profile["heights"][-1]:
            return False, False, False, False, False

        cut_indices = self.get_cut_indices_for_height_line(height, profile)
        #self.logger.log_msg("xxxxxxcut_index" + str(cut_indices))
        dam_line_segments = self.get_dam_line_segments(profile, cut_indices, height)
        #self.logger.log_msg("xxxxxxdam_line_segments" + str(dam_line_segments))
        total_length = self.get_segments_length(dam_line_segments)
        try:
            fragmentation = self.get_dam_fragmentation(dam_line_segments)
        except Exception:
            return -1, False, False, False, False
        cost = 0
        for segment in dam_line_segments:
            segment_length = segment["x_right"] - segment["x_left"]
            segment_min_h = segment["height"]
            if segment_length > 0:
                cost += 0.0039 * pow(height-segment_min_h, 1.5681) * pow(segment_length, 0.6148)
                cost += 1

        return cost, len(dam_line_segments), total_length, fragmentation, dam_line_segments

    def get_cut_indices_for_height_line(self, line_height, profile):
        """
        returns index pairs where a horizontal line of a given height input_height
        (dam height in meters) cuts through dam_line_data
        index pairs are ordered from left to right
        """
        dam_heights = profile["heights"]
        dam_dists = profile["dists"]
        index_pairs = []

        for i in range(0, len(dam_dists) - 1):
            lower_point = min(dam_heights[i], dam_heights[i + 1])
            higher_point = max(dam_heights[i], dam_heights[i + 1])
            if lower_point <= line_height <= higher_point:
                index_pairs.append([i, i + 1])

        return index_pairs

    def get_dam_line_segments(self, profile, index_pairs, dam_height):
        """get inclusion line segments from dam height, index-pairs and dam data"""
        dam_heights = profile["heights"]
        dam_dists = profile["dists"]

        dam_line_segments = []

        # outer search
        for outer_i in range(0, len(index_pairs) - 1):
            # if pair goes down to right
            if dam_heights[index_pairs[outer_i][0]] > dam_heights[index_pairs[outer_i][1]]:
                # inner search
                inner_i = outer_i + 1

                # calculate crossover x-positions between height horizontal
                # line and outer and inner index pair connections
                x_left = [dam_dists[index_pairs[outer_i][0]], dam_dists[index_pairs[outer_i][1]]]
                y_left = [dam_heights[index_pairs[outer_i][0]], dam_heights[index_pairs[outer_i][1]]]
                left_coefficients = np.polyfit(x_left, y_left, 1)

                x_right = [dam_dists[index_pairs[inner_i][0]], dam_dists[index_pairs[inner_i][1]]]
                y_right = [dam_heights[index_pairs[inner_i][0]], dam_heights[index_pairs[inner_i][1]]]
                right_coefficients = np.polyfit(x_right, y_right, 1)

                # left crossover x-value = (y-b)/m
                x_left = (dam_height - left_coefficients[1]) / left_coefficients[0]
                x_right = (dam_height - right_coefficients[1]) / right_coefficients[0]

                # find lowest point between current index pair
                # todo: do not use loop
                min_segment_h = 99999
                for pos in range(index_pairs[outer_i][1], index_pairs[inner_i][1]):
                    if dam_heights[pos] < min_segment_h:
                        min_segment_h = dam_heights[pos]

                dam_line_segments.append({"x_left": x_left,
                                          "x_right": x_right,
                                          "height": min_segment_h})

        return dam_line_segments

    def get_share_above_dam_line(self, profile, dam_line_height_abs):
        cut_indices = self.get_cut_indices_for_height_line(dam_line_height_abs, profile)
        x_start = cut_indices[0][1]
        x_end = cut_indices[-1][0]

        total_profile_points = x_end - x_start
        points_above = 0

        for h in profile["heights"][x_start: x_end]:
            if h > dam_line_height_abs:
                points_above += 1

        if points_above > 0:
            share_above = float(points_above) / float(total_profile_points)
            return share_above
        else:
            return 0.0

    def get_segments_length(self, dam_line_segments):
        dam_length = 0
        #self.logger.log_msg("xxxxxxxxx" + str(dam_line_segments))
        for segment in dam_line_segments:
            dam_length += segment["x_right"] - segment["x_left"]
        #self.logger.log_msg("xxxxxxxx" + str(dam_length))
        return dam_length

    def get_dam_fragmentation(self, dam_line_segments):
        """Calculates, how much of the dams total length is covered by the largest segment"""
        total_length = self.get_segments_length(dam_line_segments)

        # Find largest segment's length
        largest_segment = 0
        for segment in dam_line_segments:
            if segment["x_right"] - segment["x_left"] > largest_segment:
                largest_segment = segment["x_right"] - segment["x_left"]

        factor = 1.0 - (float(largest_segment)/float(total_length))
        return factor

    def get_watershed(self, dam_id):
        """Returns the watershed polygon and line for a given dam id in area true projection"""

        start_extent = arcpy.env.extent
        self.logger.log_msg("Select the next %s upriver river segments" % self.upriver_segments)
        rs_lyr = arcpy.env.workspace + "\\rs_lyr"
        arcpy.MakeFeatureLayer_management(self.river_seg_basin, rs_lyr)

        dam_pos_lyr = arcpy.env.workspace + "\\dam_pos_lyr"  # was in memory
        arcpy.MakeFeatureLayer_management(self.dam_pos_basin, dam_pos_lyr)

        # Select the dam by id
        arcpy.SelectLayerByAttribute_management(dam_pos_lyr, "NEW_SELECTION", "%s = %i" % (self.id_field, dam_id))

        # Create buffer around dam position to limit the watershed line
        self.logger.log_msg("Creating buffer around dam position")
        dam_pos_buffer = arcpy.env.workspace + "\\dam_pos_buffer_%i" % dam_id  # was in memory
        arcpy.Buffer_analysis(in_features=dam_pos_lyr, out_feature_class=dam_pos_buffer,
                              buffer_distance_or_field="20 Kilometers")

        selected_rs = []
        next_rs = deque() # a FIFO queue
        if not dam_id in self.dam_dict:
            raise Exception("The position of dam %i seems to be ambigious (i.e. further than %s away from HydroATLAS segments), it need to be checked by hand!" % (dam_id,self.initial_dam_snap_range))
        tmp_rs = self.dam_dict[dam_id][self.hydroid_field]
        #selected_rs.append(tmp_rs)
        selection_str = str(tmp_rs)
        for elem in self.river_seg_dict[tmp_rs]["NEXT_UP"]:
            next_rs.append(elem)
        counter = 0
        while counter < self.upriver_segments and next_rs:

            tmp_rs = next_rs.popleft()
            #if self.res_overlapping and self.river_seg_dict[tmp_rs]["DAM"]: # tries to prevent reservoir overlapping but doesnt work
            #    continue
            selection_str += ", " + str(tmp_rs)
            #selected_rs.append(tmp_rs)
            for elem in self.river_seg_dict[tmp_rs]["NEXT_UP"]:
                next_rs.append(elem)
            counter += 1
        #self.logger.log_msg(selection_str)
        arcpy.SelectLayerByAttribute_management(rs_lyr, "NEW_SELECTION", "%s IN (%s) " % (self.hydroid_field,selection_str))

        rs_buffer = arcpy.env.scratchWorkspace + "\\rs_buffer_%i.shp" % dam_id
        arcpy.Buffer_analysis(rs_lyr, rs_buffer, self.rs_buffer_dist)

        desc = arcpy.Describe(rs_buffer)
        arcpy.env.extent = rs_buffer

        self.logger.log_msg("Clipping the raster to the potentional dam watershed")
        dem_clipped = arcpy.env.scratchWorkspace + "\\dem_clip_%i.tif" % dam_id
        extent_str = "{} {} {} {}".format(desc.extent.XMin,desc.extent.YMin,desc.extent.XMax,desc.extent.YMax)
        dem_clipped_res = arcpy.Clip_management(self.dem_basin,extent_str,dem_clipped,rs_buffer,"9999","ClippingGeometry") # ,rs_buffer,"0","ClippingGeometry"

        #self.logger.log_msg("NEW Convert to float Raster")
        #self.logger.log_msg("NEW Resample to cubic float raster")
        #dem_float = arcpy.env.scratchWorkspace + "\\dem_float_%i.tif" % dam_id
        #arcpy.Resample_management(dem_clipped_res, dem_float, dem_clipped,"CUBIC")


        # First Fill DEM to avoid sinks    # added by Juergen 27.02.20

        self.logger.log_msg("Filling sinks")
        dem_fil = arcpy.env.scratchWorkspace + "\\dem_fil_%i.tif" % dam_id
        dem_fil_res = arcpy.sa.Fill(dem_clipped_res)
        dem_fil_res.save(dem_fil)

        '''
        # set all to high vaues to NoData
        self.logger.log_msg("Removing all areas higher than %i meters above the dam" % (self.max_dam_height + 100))
        dem_low = arcpy.env.scratchWorkspace + "\\dem_low_%i.tif" % dam_id
        max_height = self.dam_dict[dam_id][self.elevation_field] + self.max_dam_height + 100
        dem_low_res = arcpy.sa.SetNull(dem_fil_res,dem_fil_res,"VALUE > %i" % (max_height))
        dem_low_res.save(dem_low)'''

        # Create flow directions within extent of buffer
        self.logger.log_msg("Calculating flow directions")
        flow_dir = arcpy.env.scratchWorkspace + "\\flow_dir_%i.tif" % dam_id
        flow_dir_res = arcpy.sa.FlowDirection(dem_fil_res)
        flow_dir_res.save(flow_dir)

        # Create accumulation raster from flow directions
        self.logger.log_msg("Calculating flow accumulation")
        flow_acc = arcpy.env.scratchWorkspace + "\\flow_acc_%i.tif" % dam_id
        flow_acc_res = arcpy.sa.FlowAccumulation(flow_dir_res)
        flow_acc_res.save(flow_acc)

        # Clip buffer around dam point from flow acc raster
        self.logger.log_msg("Clips a buffer from the flow accumulation raster")
        dam_pos_buffer_flow = arcpy.env.workspace + "\\dam_pos_buffer_flow_%i" % dam_id  # was in memory
        arcpy.Buffer_analysis(in_features=dam_pos_lyr, out_feature_class=dam_pos_buffer_flow,
                              buffer_distance_or_field="1500 Meters",line_end_type="FLAT")

        flow_acc_clip = arcpy.env.scratchWorkspace + "\\flow_acc_clip_%i.tif" % dam_id
        flow_acc_clip_res = arcpy.Clip_management(flow_acc_res,"#",flow_acc_clip,dam_pos_buffer_flow)

        # differs between 'high' and 'low' flow values and snaps the dam to the clostest 'high' value
        self.logger.log_msg("Differs between 'high' and 'low' flow values and snaps the dam to the clostest 'high' value")
        arr = arcpy.RasterToNumPyArray(flow_acc_clip,nodata_to_value=0)
        arr_sort = np.sort(arr,axis=None)
        arr_diff = arr_sort[1:] - arr_sort[:-1] # all elements n - (n-1)
        max_index = np.argmax(arr_diff[:-2])
        thresh = int(arr_sort[max_index]+1)
        try:
            arr_diff_high = arr_sort[max_index+2:] - arr_sort[max_index+1:-1]
            max_diff_high = np.amax(arr_diff_high)
            if max_diff_high > 0.2 * np.amax(arr_diff):
                #self.logger.log_msg("XXX reached " + str(np.amax(arr_diff)) + " " + str(max_diff_high))
                max_index_high = np.argmax(arr_diff_high)
                potential_thresh = int(arr_sort[max_index+1+max_index_high]+1)
                over_thresh = arr_sort[arr_sort > potential_thresh]
                #self.logger.log_msg("xxxxxx" + str(over_thresh))
                if len(over_thresh) > 1:
                    thresh = potential_thresh
                #self.logger.log_msg("XXX thresh " + str(thresh))
        except:
            self.logger.log_msg("Something went wrong ")
        flow_acc_thresh = arcpy.env.scratchWorkspace + "\\flow_acc_thresh_%i.tif" % dam_id
        ras = arcpy.sa.Raster(flow_acc_clip)
        flow_acc_thresh_res = arcpy.sa.Con(ras > thresh, 1,0) # 1 for high values, 0 for low
        flow_acc_thresh_res.save(flow_acc_thresh)

        flow_acc_line = arcpy.env.workspace +"\\flow_acc_line_%i" % dam_id

        #tmp_snap_range = self.snap_range

        arcpy.RasterToPolyline_conversion(flow_acc_thresh_res,flow_acc_line)
        arcpy.Snap_edit(dam_pos_lyr,[[flow_acc_line,"EDGE","1500 Meters"]])

        # move the real geometry
        xy = None
        for row in arcpy.da.SearchCursor(dam_pos_lyr, ["SHAPE@XY"]):
            # Print x,y coordinates of each point feature
            xy = row[0]
            print("{}, {}".format(xy[0], xy[1]))

        cursor = arcpy.da.UpdateCursor(self.dam_positions, ["SHAPE@XY", self.id_field])
        for row in cursor:
            if dam_id == row[1]:
                row[0] = xy
                cursor.updateRow(row)





        # Snap dam position to accumulation cell with highest value in given snap range
        self.logger.log_msg("Snapping pour point")
        pour_point_raster = arcpy.env.scratchWorkspace + "\\pour_point_%i.tif" % dam_id
        arcpy.PointToRaster_conversion(dam_pos_lyr, self.id_field, pour_point_raster, cellsize=flow_acc)
        snap_pour = arcpy.env.scratchWorkspace + "\\snap_pour_%i.tif" % dam_id
        snap_pour_res = arcpy.sa.SnapPourPoint(pour_point_raster, flow_acc_res, self.snap_range)
        snap_pour_res.save(snap_pour)

        # Calculate watershed from flow directions and pour point
        self.logger.log_msg("Calculating watershed")
        watershed = arcpy.env.scratchWorkspace + "\\watershed_%i.tif" % dam_id
        watershed_res = arcpy.sa.Watershed(flow_dir, snap_pour)
        watershed_clean_res = arcpy.sa.BoundaryClean(watershed_res,"NO_SORT","ONE_WAY")
        watershed_clean_res.save(watershed)

        #watershed_raster = arcpy.env.scratchWorkspace + "\\watershed_raster_%i.tif" % dam_id
        self.dem_watershed_clip = arcpy.env.scratchWorkspace + "\\dem_watershed_clip_%i.tif" % dam_id
        watershed_raster_res = arcpy.sa.Raster(watershed)
        watershed_raster_res /= watershed_raster_res

        #watershed_raster_res.save(watershed_raster)
        dem_watershed_clip_res = arcpy.sa.Raster(dem_fil) * watershed_raster_res
        dem_watershed_clip_res.save(self.dem_watershed_clip)

        # Create watershed polygon
        self.logger.log_msg("Creating watershed polygon")
        watershed_poly = arcpy.env.workspace + "\\watershed_poly_%i" % dam_id
        arcpy.RasterToPolygon_conversion(watershed, watershed_poly, "SIMPLIFY")

        # Create line from watershed polygon
        self.logger.log_msg("Creating line from watershed polygon")
        watershed_line = arcpy.env.workspace + "\\watershed_line_%i" % dam_id  # was in memory
        arcpy.PolygonToLine_management(watershed_poly, watershed_line)

        # to get only the area aroud the dam as profile
        watershed_line_clip = arcpy.env.workspace + "\\watershed_line_clip_%i" % dam_id
        arcpy.Clip_analysis(watershed_line,dam_pos_buffer,watershed_line_clip)


        # Project results
        self.logger.log_msg("Projecting the watershed polygon and line for accurate area calculations")
        utm_projection = self.get_utm_projection(watershed_poly)
        #self.logger.log_msg("UTM: " + utm_projection)
        #projection = self.get_projections_strings(self.continent)
        #self.logger.log_msg("Projecting Results to %s" % projection[1])
        watershed_polygon_prj = "watershed_poly_projected_%i" % dam_id
        #arcpy.Project_management(watershed_poly, watershed_polygon_prj, arcpy.SpatialReference(projection[1]))
        arcpy.Project_management(watershed_poly, watershed_polygon_prj, utm_projection)
        watershed_line_prj = "watershed_line_projected_%i" % dam_id
        #arcpy.Project_management(watershed_line, watershed_line_prj, arcpy.SpatialReference(projection[1]))
        arcpy.Project_management(watershed_line_clip, watershed_line_prj, utm_projection)

        self.logger.log_msg("Adding dam and basin infos to the watershed polygon")
        self.AddField_helper(watershed_polygon_prj, "DAM_ID", "LONG")
        basin_level = str(self.basins).split("lev")[1]
        self.AddField_helper(watershed_polygon_prj, "Dam_in_BasinATLAS_lev" + str(basin_level), "DOUBLE")


        update_cursor = arcpy.UpdateCursor(watershed_polygon_prj)
        for update_row in update_cursor:
            update_row.setValue("DAM_ID", int(dam_id))
            update_row.setValue("Dam_in_BasinATLAS_lev" + str(basin_level), float(self.current_basin_id))

            update_cursor.updateRow(update_row)
        # Delete no longer needed files
        self.logger.log_msg("Deleting no longer needed files")
        self.delete_features([dam_pos_lyr,  flow_dir, snap_pour, flow_acc, flow_dir, watershed,
                             watershed_poly, watershed_line, flow_acc_line, dem_fil, #dem_low,
                              rs_lyr, rs_buffer, dem_clipped, rs_lyr, dam_pos_buffer,# dem_float,
                              dam_pos_lyr, watershed_line_clip, watershed_line,
                              dam_pos_buffer_flow, flow_acc_clip, flow_acc_thresh, pour_point_raster
                              ], silent=True) # rs_buffer,

        arcpy.env.extent = start_extent
        return [watershed_polygon_prj, watershed_line_prj]

    def get_utm_projection(self, geom):
        # Add field to write projection string to
        self.AddField_helper(geom, "X_UTM_ZONE", "TEXT", field_length=1000)
        arcpy.CalculateUTMZone_cartography(geom, "X_UTM_ZONE")
        with arcpy.da.SearchCursor(geom, ["X_UTM_ZONE"]) as cursor:
            for row in cursor:
                utm_zone = row[0]
                return utm_zone

    def read_user_inputs(self, params):
        """
        Verify user inputs and set them with correct data types to the tools instance variables.
        Returns True if everything is okay, False otherwise
        """
        dam_positions = params[0].valueAsText
        if not arcpy.Exists(dam_positions):
            self.logger.log_msg("No dam positions provided or failed exist check", "Error")
            return False
        else:
            self.dam_positions = str(dam_positions)
            self.dam_layer_name = str(dam_positions).split("\\")[-1].split(".")[0]
            self.logger.log_msg("Set dam positions to %s" % self.dam_positions)

        dem = params[1].value
        if not arcpy.Exists(dem):
            self.logger.log_msg("No DEM provided or failed exist check", "Error")
            return False
        else:
            self.dem = str(dem)
            self.logger.log_msg("Set DEM to %s" % self.dem)

        river_segment = params[2].valueAsText
        if not arcpy.Exists(river_segment):
            self.logger.log_msg("No HydroATLAS river segments provided or failed exist check", "Error")
            return False
        else:
            self.river_seg = str(river_segment)
            self.logger.log_msg("Set HydroATLAS river segments to %s" % self.river_seg)
        '''
        dis = params[2].value
        if dis is None:
            self.logger.log_msg("No discharge raster or failed exist check. "
                               "Will continue without drop height and electric output calculations", "Warning")
        elif arcpy.Exists(dis):
            self.discharge = str(dis)
            self.logger.log_msg("Set discharge raster to %s" % self.discharge)
        '''
        height_step = int(params[3].valueAsText)
        if height_step < 0:
            self.logger.log_msg("Height step has to be a positive integer", "Error")
            return False
        else:
            self.height_step = int(height_step)
            self.logger.log_msg("Set heights step to %i" % self.height_step)

        #self.continent = str(params[5].valueAsText)
        #self.logger.log_msg("Set continent to %s" % self.continent)

        self.project_folder = str(params[7].valueAsText)
        self.logger.log_msg("Set output folder to %s" % self.project_folder)

        self.project_name = str(params[6].valueAsText)
        self.logger.log_msg("Set project name to %s" % self.project_name)

        upriver_segments = int(params[4].valueAsText)
        if upriver_segments < 0:
            self.logger.log_msg("Number of upriver segments has to be a positive integer", "Error")
            return False
        else:
            self.upriver_segments = int(upriver_segments)
            self.logger.log_msg("Set number of upriver segments to %s" % self.upriver_segments)

        id_list = params[5].valueAsText
        #self.logger.log_msg(id_list)
        if id_list is None:
            self.logger.log_msg("No ID list was provided, calculations will be done for all dam positions")
            self.id_list = self.get_id_list_of_dam_positions()
            self.logger.log_msg("ID list now contains: %s" % str(self.id_list))
        else:
            id_list_split = id_list.replace(" ", "").split(",")
            self.id_list = [int(item) for item in id_list_split]
            #self.id_list = id_list
            self.logger.log_msg("ID list contains: %s" % str(self.id_list))
        #res_overlapping = params[8].valueAsText
        #if res_overlapping:
        #    self.res_overlapping = True

        if params[8]:
            self.basins = str(params[8].valueAsText)
        # if no errors occured, return True
        keep_tmp_files = params[9].valueAsText
        #self.logger.log_msg(keep_tmp_files)
        if keep_tmp_files == "true":
            self.keep_tmp_files = True

        return True

    def get_id_list_of_dam_positions(self):
        dam_id_list = []
        with arcpy.da.SearchCursor(self.dam_positions, [self.id_field]) as dam_position_cursor:
            for dam_position_row in dam_position_cursor:
                dam_id = int(dam_position_row[0])
                dam_id_list.append(dam_id)
        return dam_id_list

    def center_profile(self, profile):
        """Rearange profile data, so that the lowest point (dam) is centered on the x-axis"""

        # Find index of lowest point
        low_idx = np.argmin(profile["heights"])

        # Calculate how many places the list/array needs to be rolled
        mid_idx = int(len(profile["heights"]) / 2)
        roll = mid_idx - low_idx
        # if roll < 0:
        # roll = len(profile["heights"]) - roll
        heights_rolled = np.roll(profile["heights"], roll)
        # dists_rolled = np.roll(profile["dists"], roll)

        return {"heights": heights_rolled, "dists": profile["dists"]}

    def clean_profile(self, profile):
        """Because of possible self-intersections on lines, raw profiles sometimes contain data
           of multiple lines. This tool finds and returns the profile of the longest line of those lines
        """

        # Find all indices where dists are 0
        zero_indices = []
        for idx in range(0, len(profile["dists"])):
            if profile["dists"][idx] == 0:
                zero_indices.append(idx)

        zero_indices.append(len(profile["dists"])-1)

        if len(zero_indices) > 2:

            idx_dists = []
            for zero_idx in range(0, len(zero_indices)-1):
                idx_dists.append(zero_indices[zero_idx+1]-zero_indices[zero_idx])

            max_start = np.argmax(idx_dists)
            max_end = max_start+1

            start_idx = zero_indices[max_start]+1
            end_idx = zero_indices[max_end]-1

            cleaned_dists = profile["dists"][start_idx: end_idx]
            cleaned_heights = profile["heights"][start_idx: end_idx]

            if cleaned_dists[0] > cleaned_dists[1]:
                cleaned_dists.reverse()
            return {"heights": cleaned_heights, "dists": cleaned_dists}
        else:
            if profile["dists"][0] > profile["dists"][1]:
                profile["dists"].reverse()
            return {"heights": profile["heights"], "dists": profile["dists"]}

    def clean_profile_2(self, profile):
        dists_separated = [[]]
        heights_separated = [[]]
        for idx in range(0, len(profile["dists"])):
            if profile["dists"][idx] == 0:
                dists_separated.append([])
                heights_separated.append([])
                continue
            else:
                dists_separated[-1].append(profile["dists"][idx])
                heights_separated[-1].append(profile["heights"][idx])
        heights = []
        #self.logger.log_msg("xxxxxxxheights:" + str(heights))
        for x in heights_separated:
            if x:
                heights.append(min(x))
            else:
                heights.append(999999)
        #heights = [min(x) for x in heights_separated]
        lowest_idx = np.argmin(heights)

        lowest_dists = dists_separated[lowest_idx]
        lowest_heights = heights_separated[lowest_idx]
        lowest_heights_sort = [x for (y,x) in sorted(zip(lowest_dists,lowest_heights), key=lambda pair: pair[0])]
        lowest_dists_sort = sorted(lowest_dists)
        #self.logger.log_msg("xxxxxxxxxxlowest heights:" + str(lowest_heights_sort))
        #self.logger.log_msg("xxxxxxxxxxlowest dists:" + str(lowest_dists_sort))
        return {"heights": lowest_heights_sort, "dists": lowest_dists_sort}
        """
        # Find longest sub list
        list_sizes = [len(x) for x in dists_separated]
        longest_idx = np.argmax(list_sizes)

        longest_dists = dists_separated[longest_idx]
        longest_heights = heights_separated[longest_idx]

        # If necessary, reverse dists list
        if longest_dists[0] > longest_dists[1]:
            longest_dists.reverse()
        #arcpy.AddMessage(longest_heights)
        return {"heights": longest_heights, "dists": longest_dists}
        """

    def get_profile_as_dict(self, line):
        """Return dist heights list for line and dem"""

        # Get stackprofile from damline
        out_profile_table = "damline_profile"
        arcpy.ddd.StackProfile(line, profile_targets=[self.dem_basin], out_table=out_profile_table)

        heights = []
        dists = []

        with arcpy.da.SearchCursor(out_profile_table, ["FIRST_DIST", "FIRST_Z"]) as profile_cursor:
            for profile_row in profile_cursor:
                dists.append(profile_row[0])
                heights.append(profile_row[1])

        self.delete_features([out_profile_table], True)
        #arcpy.AddMessage(heights)
        return {"heights": heights, "dists": dists}

    def delete_features(self, features, silent=False):
        if not self.keep_tmp_files:
            for feature in features:
                """
                p = arcpy.mp.ArcGISProject("CURRENT")
                for df in p.listMaps("*"):
                    for lyr in df.listLayers("*"):
                        if lyr.name == str(feature):
                            arcpy.mp.RemoveLayer(df, lyr)
                del p
                """

                # Delete
                if arcpy.Exists(feature):
                    if not silent:
                        self.logger.log_msg("Deleting no longer needed features: %s" % str(feature))
                    try:
                        arcpy.Delete_management(feature)
                    except Exception:
                        self.logger.log_msg("Error deleting features: %s"
                                         % str(traceback.format_exc()), "Error")

    def create_folder_structure(self):
        """Create the folder structure for all script outputs, the result GDB and a logfile"""
        # Folders
        try:
            # Set workspaces for the tools runtime
            output_folder = r'%s' % self.project_folder
            os.chdir(output_folder)
            workspace = r'%s' % output_folder + "\\" + self.project_name + "\\Workspace.gdb"
            arcpy.env.workspace = workspace
            arcpy.env.scratchWorkspace = self.project_folder + "\\" + self.project_name + "\\ScratchWorkspace"

            self.project_path = self.project_folder + "/" + self.project_name
            self.graphs_path = self.project_path + "/Graphs"
            self.scratch_path = self.project_path + "/ScratchWorkspace"
            if not os.path.exists(self.project_path):
                os.makedirs(self.project_path)
            if not os.path.exists(self.graphs_path):
                os.makedirs(self.graphs_path)
            if not os.path.exists(self.scratch_path):
                os.makedirs(self.scratch_path)
        except Exception:
            self.logger.log_msg("Error creating folder structure: %s"
                             % str(traceback.format_exc()), "Error")
            return False

        self.logger = Logger(path = str(self.project_folder + "/" + self.project_name))


        # the new all incuding csv file
        tmp_file = self.project_path + "\\"+ self.dam_layer_name + "_height_steps.csv"
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        with open(tmp_file, mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            csv_file.writerow([self.id_field, "sel_hgt_m", "d_hgt_m", "r_area_msm", "r_vol_mcm",
                                "e_drop_m", "g_cap_mw", "r_sl_asl_m", "e_cap_mw", "e_cost_maud"])

        """
        # the output pc readable csv file
        tmp_file = self.project_path + "\\PC_results.csv"
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        with open(tmp_file, mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            csv_file.writerow(["dam_id", "drop_height_m", "capacity_mw", "abs_height_m", "rel_height_m",
                                "area_km2", "volume_km3", "electric_output_mw", "cost", "volume_per_cost",
                                "dam_segments", "total_length_m", "width_x_height", "fragmentation",
                                "dam_line_interruption"])

        # the output selected dam pc readable csv file
        tmp_file = self.project_path + "\\PC_selected_dams.csv"
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        with open(tmp_file, mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)

            csv_file.writerow(["dam_id", "drop_height_m", "capacity_mw", "abs_height_m", "rel_height_m",
                                "area_km2", "volume_km3", "electric_output_mw", "cost", "volume_per_cost",
                                "dam_segments", "total_length_m", "width_x_height", "fragmentation",
                                "dam_line_interruption"])
        """


        # make sure the files are empty
        tmp_file = self.project_path + "\\Human_selected_dams.csv"
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

        tmp_file = self.project_path + "\\Human_results.csv"
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

        # File Geodatabase
        try:
            if not arcpy.Exists(self.project_folder + "\\" + self.project_name + "\\Workspace.gdb"):
                arcpy.CreateFileGDB_management(self.project_folder + "\\" + self.project_name, "Workspace.gdb")
        except Exception:
            self.logger.log_msg("Error creating results GDB: %s"
                             % str(traceback.format_exc()), "Error")
            return False

        try:
            out_gdb = self.project_folder + "\\" + self.project_name + "\\Result.gdb"
            if not arcpy.Exists(out_gdb):
                arcpy.CreateFileGDB_management(self.project_folder + "\\" + self.project_name, "Result.gdb")
            self.out_gdb = str(out_gdb)
            #arcpy.AddMessage(self.out_gdb)
        except Exception:
            self.logger.log_msg("Error creating output GDB: %s"
                             % str(traceback.format_exc()), "Error")
            return False
        return True

    def export_as_graph(self, x_data, y_data, plot_name, subfolder=""):
        # Create subfolder if it not exists already
        try:
            if subfolder != "":
                subfolder_full_path = self.graphs_path + "\\" + subfolder
            else:
                subfolder_full_path = self.graphs_path
            if not os.path.exists(subfolder_full_path):
                os.makedirs(subfolder_full_path)
            try:
                plt.figure(0)
                plt.title(plot_name)
                plt.xlabel("Meters")
                plt.ylabel("Meters above sea level")
                plt.grid(True)
                plt.plot(x_data, y_data, "b-")
                str_file = subfolder_full_path + "\\%s.png" % plot_name
                # Remove Graph first, if it already exists. PyPlot does not overwrite by itself
                if os.path.isfile(str_file):
                    os.remove(str_file)
                plt.savefig(str_file)
                plt.clf()
            except Exception:
                self.logger.log_msg("Problem saving graph: %s"
                                 % str(traceback.format_exc()), "Error")
                self.logger.log_msg("Error saving graph %s: " % str(traceback.format_exc()), "Error")
        except Exception:
            self.logger.log_msg("Problem creating sub folder for graph %s" % str(traceback.format_exc()), "Error")

    def export_as_graph_with_dam_line(self, x_data, y_data, line_segments, line_height,  plot_name, dam_info, subfolder=""):
        # Create subfolder if it not exists already
        try:
            if subfolder != "":
                subfolder_full_path = self.graphs_path + "\\" + subfolder
            else:
                subfolder_full_path = self.graphs_path
            if not os.path.exists(subfolder_full_path):
                os.makedirs(subfolder_full_path)
            try:
                plot_title = plot_name
                if dam_info:
                    length = int(dam_info["best_dam"]["total_length"])
                    height = int(dam_info["best_dam"]["rel_height"])
                    drop_h = int(dam_info["estimated_drop_height"])
                    capacity = int(dam_info["given_capacity"])
                    vol = dam_info["best_dam"]["volume"]
                    area = dam_info["best_dam"]["area"]
                    electric_output = dam_info["best_dam"]["electric_output"]
                    plot_title += " | W: %i | H: %i | DH: %i \n P: %.2f | C: %i | V: %.4f | A: %.4f" \
                                  % (length, height, drop_h, electric_output, capacity, vol, area)
                """
                    W = Width [m]
                    H = Height [m]
                    DH = Drop Height [m]
                    P = Output [MW]
                    V = Volume [km^3]
                    A = Area [km^2]
                    C = Capacity [MW]

                """
                plt.figure(0)
                plt.title(plot_title)
                plt.xlabel("Meters")
                plt.ylabel("Meters above sea level")
                plt.grid(True)
                plt.plot(x_data, y_data, "b-")
                str_file = subfolder_full_path + "\\%s.png" % plot_name
                # Plot every line segment
                for segment in line_segments:
                    plt.plot([segment["x_left"], segment["x_right"]], [line_height, line_height], "g-")
                # Remove Graph first, if it already exists. PyPlot does not overwrite by itself
                if os.path.isfile(str_file):
                    os.remove(str_file)
                plt.savefig(str_file)
                plt.clf()
            except Exception:
                self.logger.log_msg("Problem saving graph: %s"
                                 % str(traceback.format_exc()), "Error")
                self.logger.log_msg("Error saving graph %s: " % str(traceback.format_exc()), "Error")
        except Exception:
            self.logger.log_msg("Problem creating sub folder for graph %s" % str(traceback.format_exc()), "Error")

    def save_csv2(self, data, best_dam):

        with open(self.project_path + "\\"+ self.dam_layer_name + "_height_steps.csv", mode="a+", newline="") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)

            # For every data list
            for d in data["per_height"]:
                str_data = [str(data["dam_id"]), str(int(best_dam["rel_height"])), str(int(d["rel_height"])), str(d["area"]),
                            str(d["volume"]), str(data["estimated_drop_height"]), str(data["given_capacity"]), str(d["abs_height"]),
                            str(d["electric_output"]), str(d["cost"])]
                csv_file.writerow(str_data)

    def save_csv(self, data):
        with open(self.project_path + "\\Human_results.csv", mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            csv_file.writerow(["", "", "", "", "", "", "", "", ""])
            # Write column headers
            csv_file.writerow(["Dam: %i" % data["dam_id"],
                               "Drop Height %s [m]" % str(int(data["estimated_drop_height"])),
                               "Capacity [MW]: %s" % str(int(data["given_capacity"])), "", "", "", "", "", "", "", ""])

            csv_file.writerow(["Abs Height [m]", "Rel Height [m]", "Area [km^2]", "Volume [km^3]",
                               "Electric Output [MW]", "Cost [Mio. $]", "Volume per Cost",
                               "Dam Segments", "Total Length [m]", "Width x Height", "Fragmentation",
                               "Dam Line Interruption"])

            # For every data list
            for d in data["per_height"]:
                str_data = [str(d["abs_height"]).replace(".", ","), str(d["rel_height"]).replace(".", ","),
                            str(d["area"]).replace(".", ","), str(d["volume"]).replace(".", ","),
                            str(d["electric_output"]).replace(".", ","), str(d["cost"]).replace(".", ","),
                            str(d["vol_cost"]).replace(".", ","), str(d["segment_count"]).replace(".", ","),
                            str(d["total_length"]).replace(".", ","),
                            str(d["total_length"] * d["rel_height"]).replace(".", ","),
                            str(d["fragmentation"]).replace(".", ","), str(d["share_above_dam_line"]).replace(".", ",")]
                csv_file.writerow(str_data)

        with open(self.project_path + "\\PC_results.csv", mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)

            # For every data list
            for d in data["per_height"]:
                str_data = [str(data["dam_id"]), str(int(data["estimated_drop_height"])), str(int(data["given_capacity"])),
                            str(d["abs_height"]), str(d["rel_height"]), str(d["area"]), str(d["volume"]),
                            str(d["electric_output"]), str(d["cost"]), str(d["vol_cost"]), str(d["segment_count"]),
                            str(d["total_length"]), str(d["total_length"] * d["rel_height"]),
                            str(d["fragmentation"]), str(d["share_above_dam_line"])]
                csv_file.writerow(str_data)

    def save_selected_dam_csv(self, data):
        with open(self.project_path + "\\Human_selected_dams.csv", mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)

            # For every data list
            str_data = ["Dam", str(data["dam_id"]), "Drop Height [m]", str(int(data["estimated_drop_height"])),
                        "Capacity [MW]", str(int(data["given_capacity"])),
                        "Abs Height [m]", str(data["abs_height"]).replace(".", ","),
                        "Rel Height [m]", str(data["rel_height"]).replace(".", ","),
                        "Area [km^2]", str(data["area"]).replace(".", ","),
                        "Volume [km^3]", str(data["volume"]).replace(".", ","),
                        "Electric Output [MW]", str(data["electric_output"]).replace(".", ","),
                        "Cost [Mio. $]", str(data["cost"]).replace(".", ","),
                        "Volume per Cost", str(data["vol_cost"]).replace(".", ","),
                        "Dam Segments", str(data["segment_count"]).replace(".", ","),
                        "Total Length [m]", str(data["total_length"]).replace(".", ","),
                        "Width x Height", str(data["total_length"] * data["rel_height"]).replace(".", ","),
                        "Fragmentation", str(data["fragmentation"]).replace(".", ","),
                        "Dam Line Interruption", str(data["share_above_dam_line"]).replace(".", ","),
                        "Comment", str(data["comment"]),
                        "Decision", str(data["decision"])]

            csv_file.writerow(str_data)

        with open(self.project_path + "\\PC_selected_dams.csv", mode="a+") as csv_file:
            csv_file = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            str_data = [str(data["dam_id"]), str(int(data["estimated_drop_height"])), str(int(data["given_capacity"])),
                        str(data["abs_height"]), str(data["rel_height"]), str(data["area"]), str(data["volume"]),
                        str(data["electric_output"]), str(data["cost"]), str(data["vol_cost"]), str(data["segment_count"]),
                        str(data["total_length"]), str(data["total_length"] * data["rel_height"]),
                        str(data["fragmentation"]), str(data["share_above_dam_line"])]
            csv_file.writerow(str_data)

    def AddField_helper(self, fc, field, type, **kwargs):
        "Helper function for the build in AddField to check for field existence before run"
        #self.logger.log_msg(arcpy.ListFields(fc))
        if not field in [f.name for f in arcpy.ListFields(fc)]:
            arcpy.AddField_management(fc, field, type, **kwargs)



    def get_projections_strings(self, study_area):
        """Get angle and area true projection strings based on the study area"""
        # "Africa", "North Asia", "South Asia", "North America", "South America", "Canada"
        projections = ["", ""]
        if study_area == "Africa":
            projections = ["Africa Lambert Conformal Conic", "Africa Albers Equal Area Conic"]
        elif study_area == "North Asia":
            projections = ["Asia North Lambert Conformal Conic", "Asia North Albers Equal Area Conic"]
        elif study_area == "South Asia":
            projections = ["Asia South Lambert Conformal Conic", "Asia South Albers Equal Area Conic"]
        elif study_area == "North America":
            projections = ["North America Lambert Conformal Conic", "North America Albers Equal Area Conic"]
        elif study_area == "South America":
            projections = ["South America Lambert Conformal Conic", "South America Albers Equal Area Conic"]
        elif study_area == "Canada":
            projections = ["Canada Lambert Conformal Conic", "Canada Albers Equal Area Conic"]

        return projections

    def write_settings_to_log(self):
        self.logger.log_msg("############ SETTINGS USED ############")
        self.logger.log_msg("Dam Points: %s" % self.dam_positions)
        self.logger.log_msg("DEM: %s" % self.dem)
        #self.logger.log_msg("Discharge: %s" % self.discharge)
        self.logger.log_msg("HydroATLAS" + str(self.river_seg))
        self.logger.log_msg("Step Size: %i" % self.height_step)
        self.logger.log_msg("Upriver segments: %s" % self.upriver_segments)
        self.logger.log_msg("Snap Range: %s" % str(self.snap_range))
        self.logger.log_msg("Max Dam Height: %i" % self.max_dam_height)
        self.logger.log_msg("Max Dam Width: %i" % self.max_dam_width)
        self.logger.log_msg("Max Reservoir Area: %i" % self.max_reservoir_area)
        self.logger.log_msg("Max Reservoir Volume: %i" % self.max_reservoir_volume)
        self.logger.log_msg("Max Dam Line Interruption: %f" % self.max_dam_interruption)
        self.logger.log_msg("#######################################")

class Enumerate_river_segments(object):
    def __init__(self, logfile_path = None):
        """Define the tool (tool name is the name of the class)."""
        self.label = "I. Enumerate river segments"
        self.description = "A tool to enumerate HydroATLAS river segments and potentially use a threshold to delete unnecessary segments with low discharge."
        self.canRunInBackground = False
        self.logger = Logger(path=logfile_path)

    def getParameterInfo(self):
        """ Parameter definitions """
        # RiverATLAS river segments
        param_river_segment = arcpy.Parameter(
            displayName = "RiverATLAS river segments",
            name = "in_river_segments",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Input")

        # Threshold for the max. annual discharge
        param_threshold = arcpy.Parameter(
            displayName = "Threshold for max. ann. discharge",
            name = "in_threshold",
            datatype = "GPDouble",
            parameterType = "Optional",
            direction = "Input")

        # Create new FC with deletes features?
        param_delete_below_threshold = arcpy.Parameter(
            displayName = "Create a new Feature Class with the enumeration and deletet features below the threshold?",
            name = "in_delete_below_thresh",
            datatype = "GPBoolean",
            parameterType = "Optional",
            direction = "Input")

        param_river_segment_out = arcpy.Parameter(
            displayName = "Output river network (only if delete is checked)",
            name = "out_river_segments",
            datatype = "DEFeatureClass",
            parameterType = "Optional",
            direction = "Output")


        params = [param_river_segment, param_threshold, param_delete_below_threshold, param_river_segment_out]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        #self.logger.log_msg("Executing 'Enumerate river segments'")
        arcpy.env.overwriteOutput = True # override

        # Step 01/04: Gets the parameters
        #arcpy.AddMessage("Step 01/04: Gets the parameters")
        self.logger.log_msg("Step 01/04: Gets the parameters")
        if not isinstance(parameters[0], str):
            river_segment = parameters[0].valueAsText
            threshold = float(parameters[1].valueAsText)
            delete_below_threshold = parameters[2].value
            river_segment_out = parameters[3].valueAsText
        else:
            river_segment = parameters[0]
            threshold = float(parameters[1])
            delete_below_threshold = parameters[2]
            river_segment_out = parameters[3]
        #df = pd.DataFrame.spatial.from_featureclass(river_segment)

        # Step 02/04: Interates through all river segments (by HYRIV_ID) and numbers them
        #arcpy.AddMessage("Step 02/04: Interates through all river segments (by HYRIV_ID) and numbers them")
        self.logger.log_msg("Step 02/04: Interates through all river segments (by HYRIV_ID) and numbers them")
        Rs_dict_rivid = {} # dictionary for all river segments
        cursor = arcpy.da.SearchCursor(river_segment, ["HYRIV_ID","NEXT_DOWN","dis_m3_pmx"]) # HydroRiverID & ID of the adjacent downstream river segment
        for row in cursor:
            Rs_dict_rivid[row[0]]={"NEXT_DOWN":row[1],"DIS":row[2],"COUNTER":0,"SEQ_NR":0} # adds all river segments by their HydroID with Informations about their downstream element NEXT_DOWN
        counter = 1 # initiates the counter as 1
        for rivid in Rs_dict_rivid: # iterates through the dictionary
            if Rs_dict_rivid[rivid]["COUNTER"] == 0 and Rs_dict_rivid[rivid]["DIS"] > threshold: # if the river segment has no counter jet
                Rs_dict_rivid[rivid]["COUNTER"] = counter # set it to the current counter in the dictionary
                tmp_id = rivid
                tmp_down = Rs_dict_rivid[rivid]["NEXT_DOWN"] # get the downstream element
                counter += 1 # increment the counter
                while not tmp_down == 0: # interate through all downriver river segments until the sea is reached
                    Rs_dict_rivid[tmp_down]["COUNTER"]=counter # and give them incrementing numbers
                    counter += 1
                    tmp_id = tmp_down
                    tmp_down = Rs_dict_rivid[tmp_id]["NEXT_DOWN"]


        # Step 03/04: Sorts all river segments by their number and renumbers them with their SEQ_NR
        #arcpy.AddMessage("Step 03/04: Sorts all river segments by their number and renumbers them with their SEQ_NR")
        self.logger.log_msg("Step 03/04: Sorts all river segments by their number and renumbers them with their SEQ_NR")
        Rs_array = np.zeros((len(Rs_dict_rivid),2)) # empty array with Nx2 Fields, where N is the number of river segments in the dictionary
        array_counter = 0
        for rivid in Rs_dict_rivid:
            Rs_array[array_counter]=[rivid,Rs_dict_rivid[rivid]["COUNTER"]] # stores the informations in the array
            array_counter += 1
        Rs_array_sorted = Rs_array[Rs_array[:,1].argsort()] # sorts the array by each COUNTER
        sequenceNr = 1 # initates the sequenceNr as 1
        for river_seg in Rs_array_sorted:
            if river_seg[1] > 0: # if the segment has enough discharge or is needed between other segments
                Rs_dict_rivid[round(river_seg[0])]["SEQ_NR"]=sequenceNr # sets the SEQ_NR in the dictionary
                sequenceNr += 1
            else:
                Rs_dict_rivid[round(river_seg[0])]["SEQ_NR"]=0 # sets the SEQ_NR in the dictionary


        # Step 04/04: Adds the field and the value of the SEQ_NR to the attribute table
        #arcpy.AddMessage("Step 04/04: Adds the field and the value of the SEQ_NR to the attribute table")
        self.logger.log_msg("Step 04/04: Adds the field and the value of the SEQ_NR to the attribute table")
        if not delete_below_threshold:
            arcpy.AddField_management(river_segment,"SEQ_NR","LONG")
            cursor = arcpy.da.UpdateCursor(river_segment, ["HYRIV_ID","SEQ_NR"])
        else:
            #arcpy.AddMessage("Step 05/04: (optional) Creates the new feature class with river segments above the threshold")
            self.logger.log_msg("Step 05/04: (optional) Creates the new feature class with river segments above the threshold")
            #river_segment_tresh = str(river_segment) + "_thld" + str(int(threshold))
            #arcpy.AddMessage(river_segment)
            arcpy.Copy_management(river_segment,river_segment_out)
            arcpy.AddField_management(river_segment_out,"SEQ_NR","LONG")
            cursor = arcpy.da.UpdateCursor(river_segment_out, ["HYRIV_ID","SEQ_NR"])
        for row in cursor:
            tmp_seq_nr = Rs_dict_rivid[row[0]]["SEQ_NR"]
            if tmp_seq_nr > 0:
                row[1] = tmp_seq_nr
                cursor.updateRow(row)
            elif delete_below_threshold:
                cursor.deleteRow()
        return

class Combine_dem_rasters(object):
    """Define the tool (tool name is the name of the class)."""

    def __init__(self):
        self.label = "Combine DEM rasters"
        self.description = "Gathers all still zip-compressed rasters and combines them to one raster"
        self.canRunInBackground = False
        self.logger = Logger()

        self.unzip_folder_name = "unzip"


    def getParameterInfo(self):
        """Define parameter definitions"""

        # Folder containing the zip-files
        param_zip_folder = arcpy.Parameter(
            displayName="Folder containing the zip-compressed DEMs",
            name="in_zip_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")
        #param_zip_folder.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\dem_con\\af_con_3s_zip_grid"

        # output DEM file
        param_dem = arcpy.Parameter(
            displayName="Output DEM file",
            name="out_dem",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Output")
        #param_dem.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\dem_con\\af_con_3s.tif"


        params = [param_zip_folder, param_dem]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):

        """The source code of the tool."""
        arcpy.env.overwriteOutput = True # override
        zip_folder = str(parameters[0].value)
        out_dem = str(parameters[1].value)
        out_dem_split = out_dem.rsplit("\\",1)
        out_dem_folder = out_dem_split[0]
        out_dem_name = out_dem_split[1]
        self.logger.log_msg(out_dem_folder + " " + out_dem_name)
        zip_list = os.listdir(zip_folder)
        self.logger.log_msg("Unzipping %i DEM files" % len(zip_list))
        unzip_folder = zip_folder + "\\" + self.unzip_folder_name
        counter = 1
        dem_list = []
        for zip_file in zip_list:
            self.logger.log_msg("DEM file %i of %i" % (counter,len(zip_list)))
            counter += 1
            with zipfile.ZipFile(zip_folder + "\\" + zip_file, 'r') as zip_read:
                zip_read.extractall(unzip_folder)
                dem_name = str(zip_file).split(".")
                dem_name = dem_name[0] + "." + dem_name[2]
                dem_path = unzip_folder + "\\" + dem_name
                dem_list.append(dem_path)
        self.logger.log_msg("Combining all DEM files to one file")
        arcpy.MosaicToNewRaster_management(input_rasters=dem_list,output_location=out_dem_folder,
            raster_dataset_name_with_extension=out_dem_name,pixel_type="16_BIT_SIGNED",
            number_of_bands=1,mosaic_method="MEAN")

        """
        arcpy.env.overwriteOutput = True # override
        zip_folder = str(parameters[0].value)
        out_dem = str(parameters[1].value)
        out_dem_split = out_dem.rsplit("\\",1)
        out_dem_folder = out_dem_split[0]
        out_dem_name = out_dem_split[1]

        zip_list = os.listdir(zip_folder)
        self.logger.log_msg("Unzipping %i DEM files" % len(zip_list))
        unzip_folder = zip_folder + "\\" + self.unzip_folder_name
        counter = 1
        dem_list = []
        for zip_file in zip_list:
            self.logger.log_msg("DEM file %i of %i" % (counter,len(zip_list)))
            counter += 1
            with zipfile.ZipFile(zip_folder + "\\" + zip_file, 'r') as zip_read:
                zip_read.extractall(unzip_folder)
                dem_name = str(zip_file).split("_")
                dem_name = dem_name[0] + "_" + dem_name[1]
                dem_path = unzip_folder + "\\" + dem_name + "\\" + dem_name
                dem_list.append(dem_path)

        self.logger.log_msg("Combining all DEM files to one file")
        arcpy.MosaicToNewRaster_management(input_rasters=dem_list,output_location=out_dem_folder,
            raster_dataset_name_with_extension=out_dem_name,pixel_type="16_BIT_SIGNED",
            number_of_bands=1,mosaic_method="MEAN")
        """

class Merge_results(object):
    """Define the tool (tool name is the name of the class)."""

    def __init__(self):
        self.label = "Merge automated reservoir estimation results"
        self.description = "Takes different input GDBs and merges to contents in a new GDB per 'type' "
        self.canRunInBackground = False
        self.logger = Logger()

        self.unzip_folder_name = "unzip"

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Folder containing the zip-files
        param_folders = arcpy.Parameter(
            displayName="GDBs containing the data",
            name="in_folders",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input",
            multiValue=True)
        #param_zip_folder.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\dem_con\\af_con_3s_zip_grid"

        param_suffix = arcpy.Parameter(
            displayName="New feature suffix (e.g. 'eu')",
            name="in_suffix",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        # output DEM file
        param_output_folder = arcpy.Parameter(
            displayName="GDB to save the new features",
            name="new_gdb",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        #param_dem.value = "C:\\Users\\leobu\\Documents\\GIS\\M.Sc.Projekt\\Daten\\dem_con\\af_con_3s.tif"


        params = [param_folders, param_suffix, param_output_folder]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        arcpy.env.overwriteOutput = True # override

        in_folders_read = str(parameters[0].value)
        suffix = str(parameters[1].value)
        out_folder = str(parameters[2].value)
        in_folders = [x.strip() for x in str(in_folders_read).split(";")]

        old_workspace = arcpy.env.workspace

        fhreds = []
        reservoirs = []
        waterlines = []
        waterpolys = []

        human_res = []
        human_select = []
        pc_res = []
        pc_select = []

        for folder in in_folders:
            arcpy.env.workspace = folder
            features = arcpy.ListFeatureClasses()
            self.logger.log_msg("Reading " + str(len(features)) + " features in folder " + str(folder))
            for feature in features:
                if "FHReD" in str(feature):
                    fhreds.append(folder + "\\" + feature)
                elif "reservoir" in str(feature):
                    reservoirs.append(folder + "\\" + feature)
                elif "line" in str(feature):
                    waterlines.append(folder + "\\" + feature)
                elif "watershed_poly" in str(feature):
                    waterpolys.append(folder + "\\" + feature)
            #up_folder = str(folder)[:-10]
            #self.logger.log_msg(up_folder)
            up_folder = folder + "\\.."
            files = os.listdir(up_folder)
            for file in files:
                #self.logger.log_msg(file)
                if "Human_results" in str(file):
                    human_res.append(up_folder + "\\" + file)
                elif "Human_selected" in str(file):
                    human_select.append(up_folder + "\\" + file)
                elif "PC_results" in str(file):
                    pc_res.append(up_folder + "\\" + file)
                elif "PC_selected" in str(file):
                    pc_select.append(up_folder + "\\" + file)

        up_out_folder = out_folder + "\\.."
        out_fhred = out_folder + "\\FHReD"
        out_reservoirs = out_folder + "\\reservoirs_polygons"
        out_waterlines = out_folder + "\\watershed_lines_projected"
        out_waterpolys = out_folder + "\\watershed_polys_projected"
        out_human_res = up_out_folder + "\\Human_results"
        out_human_select = up_out_folder + "\\Human_selected_dams"
        out_pc_res = up_out_folder + "\\PC_results"
        out_pc_select = up_out_folder + "\\PC_selected_dams"

        if suffix:
            out_fhred += "_" + str(suffix)
            out_reservoirs += "_" + str(suffix)
            out_waterlines += "_" + str(suffix)
            out_waterpolys += "_" + str(suffix)
            out_human_res += "_" + str(suffix)
            out_human_select += "_" + str(suffix)
            out_pc_res += "_" + str(suffix)
            out_pc_select += "_" + str(suffix)


        #self.logger.log_msg("Merging all found files")
        arcpy.env.workspace = out_folder
        arcpy.env.extent = "MAXOF"
        if fhreds:
            self.logger.log_msg("Merging all FHReD position features (%i)" % len(fhreds))
            arcpy.Merge_management(fhreds,out_fhred)
            self.delete_duplicates(out_fhred, "DAM_ID", "Capacity__MW_") #"FHReD_area_msm"
        if reservoirs:
            self.logger.log_msg("Merging all FHReD reservoir features (%i)" % len(reservoirs))
            arcpy.Merge_management(reservoirs,out_reservoirs)
            self.delete_duplicates(out_reservoirs, "DAM_ID", "Reservoir_area_msm")
        if waterlines:
            self.logger.log_msg("Merging all FHReD watershed lines (%i)" % len(waterlines))
            arcpy.Merge_management(waterlines,out_waterlines)
        if waterpolys:
            self.logger.log_msg("Merging all FHReD watershed polygons (%i)" % len(waterpolys))
            #arcpy.Intersect_analysis(waterpolys,out_waterpolys)
            arcpy.Merge_management(waterpolys,out_waterpolys)
            self.delete_duplicates(out_waterpolys, "DAM_ID", "Shape_Area")
        if human_res:
            self.logger.log_msg("Concatinating all human readable results (%i)" % len(human_res))
            self.concat_csv(human_res,out_human_res, 0)
        if human_select:
            self.logger.log_msg("Concatinating all human readable selected results (%i)" % len(human_select))
            self.concat_csv(human_select,out_human_select, 0)
        if pc_res:
            self.logger.log_msg("Concatinating all pc readable results (%i)" % len(pc_res))
            self.concat_csv(pc_res,out_pc_res, 1)
        if pc_select:
            self.logger.log_msg("Concatinating all pc readable selected results (%i)" % len(pc_select))
            self.concat_csv(pc_select,out_pc_select, 1)

        arcpy.env.workspace = old_workspace

    def delete_duplicates(self, fc, id_field, area_field):
        cursor = arcpy.da.SearchCursor(fc, ["OBJECTID", id_field, area_field], sql_clause=(None,'ORDER BY ' + str(id_field) + ' ASC') )
        last_dam = [None, None,None]
        to_delete = []
        for row in cursor:
            if row[1] == last_dam[1]:
                self.logger.log_msg("now: " + str(row[1]) + " vs last: " + str(last_dam[1]))
                if row[2] > last_dam[2]:
                    to_delete.append(last_dam[0])
                else:
                    to_delete.append(row[0])
                    continue
            last_dam = [row[0],row[1],row[2]]
        self.logger.log_msg("Deleting " + str(len(to_delete)) +" row(s) in " + str(fc))
        del_cursor = arcpy.da.UpdateCursor(fc, ["OBJECTID"])
        for row in del_cursor:
            if row[0] in to_delete:
                del_cursor.deleteRow()

    def concat_csv(self, in_tables, out_table, head):

        first_table = True
        to_write = []
        for table in in_tables:
            skip_header = head
            with open(table) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=";", quotechar="'")
                for row in csv_reader:
                    if skip_header:
                        skip_header = 0
                        if first_table:
                            to_write.append(row)
                            first_table = False
                    else:
                        to_write.append(row)
        with open(out_table + ".csv", 'ab+') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=";", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            for row in to_write:
                if row:
                    csv_writer.writerow(row)


class Logger(object):
    """Handels the logs to a text file and to the console"""
    def __init__(self, path=None):
        # Logfile

        self.logfile_path = None
        if path:
            try:
                self.logfile_path = path + "/log.txt"
                logfile = open(self.logfile_path, "a+")
                logfile.close()
            except Exception:
                self.log_msg("Error creating logfile: %s"
                                 % str(traceback.format_exc()), "Error")
                return False
        return
    def log_msg(self, message, type="Message"):
        """Logs a message to the console and logfile"""

        # Conditionally show message to user
        message_call = arcpy.AddMessage
        if type == "Warning":
            message_call = arcpy.AddWarning
        elif type == "Error":
            message_call = arcpy.AddError
        message_call(message)

        # Log to logfile if a logfile is given
        if self.logfile_path:
            try:
                logfile = open(self.logfile_path, "a+")
                logfile.write("[%s]: %s \n" % (type, message))
            except Exception as e:
                self.logger.log_msg("[log_message] Error writing to logfile: %s"
                                 % str(traceback.format_exc()), "Error")
            finally:
                logfile.close()

    def check_SA_licence(self):
        """Checks the Spacial Analyst licence"""
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            self.log_msg("Unable to get spatial analyst extension","Error")
            self.log_msg(arcpy.GetMessages(0))
            sys.exit(0)
