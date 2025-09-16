#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 11 18:33:30 2025

@author: bobco-08
"""

import copernicusmarine

copernicusmarine.subset(
  dataset_id="cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m",
  dataset_version="202311",
  variables=["spco2"],
  minimum_longitude=40,
  maximum_longitude=110,
  minimum_latitude=-40,
  maximum_latitude=30,
  start_datetime="2022-01-01T00:00:00",
  end_datetime="2025-09-19T00:00:00",
  coordinates_selection_method="strict-inside",
  netcdf_compression_level=1,
  disable_progress_bar=True,
)