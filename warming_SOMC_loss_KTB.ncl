; Nov 2015
; clm45_CN_pools/clm45_CN_pools.ncl
; Will Wieder

; read & plots  changes in soil C from HWSD, SoilGrids using Crowther & Bradford's emperical estimates
;       (a) 2005 - 2100
; soil temperature data from Dave's permafrost RCN runs, rcp8.5 run with anomoly forcing
; OR
; TSOI from CMIP5 results using CESM1 (clm4.0)
; makes figures 

;*********************************************************************
; reads in spin-up results from control (nfixNPP) & modified (nfixET) runs 
; plots global changes in Nfix & Nlimitation
;*********************************************************************

load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_LIB/ncarg/nclscripts/csm/shea_util.ncl"

begin
;---------------------------------------------------------------
;Read in input variables
;--------------------------------------------------------------
 HWSD_path  = (/"/project/tss/wwieder/soilCN/global_run/input_data/HWSD/surfdata_0.9x1.25_simyr2000_c120906_HWSD_soil.nc"/)
 HWSD_data  = addfile(HWSD_path, "r")
 HWSD_SOC   = HWSD_data->AWT_SOC(0,:,:)   ; Area weighted SOC pools, kgC/m2, 0-30 cm
 HWSD_SOC   = HWSD_SOC * 0.6              ; 60% of SOC (0-30 cm) in top 15 cm (from soilgrids)
 HWSD_depth = 15.                         ; here assumed...
 areaupsc   = HWSD_data->areaupsc
 landfrac   = HWSD_data->landfrac
 grid_area  = areaupsc * landfrac

 HWSD_SOC = where(HWSD_SOC .eq. 0, 1.e36, HWSD_SOC)
 HWSD_SOC@_FillValue = 1.e36
 printVarSummary(HWSD_SOC)
 globHWSD = sum(HWSD_SOC * 1e6 * grid_area) / 1e12 ;PgC globally
 print(globHWSD)
 HWSD_mask = where(ismissing(HWSD_SOC), 0, 1)
 
 SGRD_path  = (/"/project/tss/wwieder/soilCN/SoilGrids/SoilGrids_0.9x1.25.nc"/)
 SGRD_data  = addfile(SGRD_path, "r")
 SGRD_soc   = SGRD_data->OCSTHA_M(0:1,:,:)    ;SOC (0-15cm tonnes ha^-1)
 SGRD_SOC   = dim_sum_n_Wrap(SGRD_soc, 0)
 SGRD_depth = 15.
 SGRD_SOC   = SGRD_SOC / 10.                ;convert kg C m^2
 SGRD_SOC@units = "kgC/m2"

 SGRD_SOC = where(SGRD_SOC .eq. 0, -9999., SGRD_SOC)
 printVarSummary(SGRD_SOC)
 globSGRD = sum(SGRD_SOC * 1e6 * grid_area) / 1e12
 print(globSGRD)

;--------------------------------------------------------------
; For Dave's CLM4.5 using rcp8.5 anamoly forcing...

 path  = (/"/project/tss/wwieder/clm45_CN_pools/"/)
 EXP   = 0                                     ; 0 = CLM45 fully transient, w/ Cru-NCEP
 title = (/"clm45_fullTransient"/)
 model = (/"clm45bgc_hrv_1deg4508_"/)
 exper = (/"rcp85"/)
 foot  = (/".monthly.2006-2299.nc"/)
; fout  = path + title(EXP)
 var   = (/"TSOI"/)
; data  = addfile(path+model+exper(0)+"/"+model+exper(0)+"."+var(0) +foot(0), "r")   
; tsoi  = data->$var(0)$(0:1139,2,:,:)            ; only read data 3rd soil layer through 2100

;--------------------------------------------------------------
; For CESM1.0 (clm4.0) results in CMIP5 archive

 rcp26_path = (/"/project/tss/wwieder/soilCN/global_run/input_data/CMIP5_data/rcp26/tsl_Lmon_CESM1-CAM5_rcp26_r1i1p1_200601-210012.nc"/)
 data26    = addfile(rcp26_path, "r")
 Tsoi26    = data26->tsl(:,0:3,:,:)    ; select first four layers, 0-16.6 cm
 depth     = data26->depth

 rcp45_path = (/"/project/tss/wwieder/soilCN/global_run/input_data/CMIP5_data/rcp45/tsl_Lmon_CESM1-CAM5_rcp45_r1i1p1_200601-210012.nc"/)
 data45    = addfile(rcp45_path, "r")
 Tsoi45    = data45->tsl(:,0:3,:,:)

 rcp60_path = (/"/project/tss/wwieder/soilCN/global_run/input_data/CMIP5_data/rcp60/tsl_Lmon_CESM1-CAM5_rcp60_r1i1p1_200601-210012.nc"/)
 data60    = addfile(rcp60_path, "r")
 Tsoi60    = data60->tsl(:,0:3,:,:)

 rcp85_path = (/"/project/tss/wwieder/soilCN/global_run/input_data/CMIP5_data/rcp85/tsl_Lmon_CESM1-CAM5_rcp85_r1i1p1_200601-210012.nc"/)
 data85    = addfile(rcp85_path, "r")
 Tsoi85    = data85->tsl(:,0:3,:,:)

;---------------------------------------------------------------
; manipulate soil temperature data to calculate depth averaged annual TSOI
;---------------------------------------------------------------

 nlevground = dimsizes(depth)
 bottom     = new(nlevground, double)   ; bottom of each soil layer
 do i = 0, nlevground - 2
   bottom(i) = avg(depth(i:i+1))
 end do

 dz    = new(nlevground, double)   ; thickness of soil layer (m)
 dz(0) = avg(depth(0:1))
 do i  = 1, nlevground-2
   dz(i) = avg(depth(i:i+1)) - sum(dz)
 end do

 frw = sum(dz(0:3))    ; soil thickness (frw) summed for first 4 laysers
 wgt = new(4, double)  ; calculate relative weight for each layer- 

 do i = 0, 3
   wgt(i) = dz(i) / frw
 end do

 tsoi26 = dim_avg_wgt_n_Wrap(Tsoi26, wgt,1,1)   ; calculated depth weighted soil temperature
 tsoi45 = dim_avg_wgt_n_Wrap(Tsoi45, wgt,1,1)  
 tsoi60 = dim_avg_wgt_n_Wrap(Tsoi60, wgt,1,1)   
 tsoi85 = dim_avg_wgt_n_Wrap(Tsoi85, wgt,1,1)   


 TSOI26 = month_to_annual(tsoi26, 1)    ; calculate annula mean soil temperature
 TSOI45 = month_to_annual(tsoi45, 1)
 TSOI60 = month_to_annual(tsoi60, 1)
 TSOI85 = month_to_annual(tsoi85, 1)     ; convert to annual average
 ndims  = dimsizes(TSOI85)
 nyear  = ndims(0)
 nlat   = ndims(1)
 nlon   = ndims(2) 
 TSOI85!0 = "year"
 printVarSummary(TSOI85)
 year  = ispan(2006,2100,1)
 TSOI85&year = year
 lat   = TSOI85&lat
 lon   = TSOI85&lon
 copy_VarCoords(TSOI85, TSOI26)
 copy_VarCoords(TSOI85, TSOI45)
 copy_VarCoords(TSOI85, TSOI60)

 iTSOI26 = dim_avg_n_Wrap(TSOI26(0:9,:,:), 0)
 fTSOI26_2050 = dim_avg_n_Wrap(TSOI26(40:49,:,:), 0)
 fTSOI26_2100 = dim_avg_n_Wrap(TSOI26(nyear-10:nyear-1,:,:), 0)
 dTSOI26_2050 = fTSOI26_2050-iTSOI26
 dTSOI26_2100 = fTSOI26_2100-iTSOI26

 iTSOI45 = dim_avg_n_Wrap(TSOI45(0:9,:,:), 0)
 fTSOI45_2050 = dim_avg_n_Wrap(TSOI45(40:49,:,:), 0)
 fTSOI45_2100 = dim_avg_n_Wrap(TSOI45(nyear-10:nyear-1,:,:), 0)
 dTSOI45_2050 = fTSOI45_2050-iTSOI45
 dTSOI45_2100 = fTSOI45_2100-iTSOI45

 iTSOI60 = dim_avg_n_Wrap(TSOI60(0:9,:,:), 0)
 fTSOI60_2050 = dim_avg_n_Wrap(TSOI60(40:49,:,:), 0)
 fTSOI60_2100 = dim_avg_n_Wrap(TSOI60(nyear-10:nyear-1,:,:), 0)
 dTSOI60_2050 = fTSOI60_2050-iTSOI60
 dTSOI60_2100 = fTSOI60_2100-iTSOI60

 iTSOI85 = dim_avg_n_Wrap(TSOI85(0:9,:,:), 0)
 fTSOI85_2050 = dim_avg_n_Wrap(TSOI85(40:49,:,:), 0)
 fTSOI85_2100 = dim_avg_n_Wrap(TSOI85(nyear-10:nyear-1,:,:), 0)
 dTSOI85_2050 = fTSOI85_2050-iTSOI85
 dTSOI85_2100 = fTSOI85_2100-iTSOI85

 print("mean global dTSOI26 (2050) = "+ avg(dTSOI26_2050)+" (2100) = "+ avg(dTSOI26_2100))
 print("mean global dTSOI45 (2050) = "+ avg(dTSOI45_2050)+" (2100) = "+ avg(dTSOI45_2100))
 print("mean global dTSOI60 (2050) = "+ avg(dTSOI60_2050)+" (2100) = "+ avg(dTSOI60_2100))
 print("mean global dTSOI85 (2050) = "+ avg(dTSOI85_2050)+" (2100) = "+ avg(dTSOI85_2100))

;--------------------------------------------------------------
;calculate running mean soil temperature & change from 2006-2100 mean
 TSOI26_s  = runave_n_Wrap(TSOI26, 5, 1, 0)        ;calculates 10 year running average
 TSOI45_s  = runave_n_Wrap(TSOI45, 5, 1, 0)        ;calculates 10 year running average
 TSOI60_s  = runave_n_Wrap(TSOI60, 5, 1, 0)        ;calculates 10 year running average
 TSOI85_s  = runave_n_Wrap(TSOI85, 5, 1, 0)        ;calculates 10 year running average
 printVarSummary(TSOI85_s)

 TSOI26_i = conform_dims(dimsizes(TSOI26_s), TSOI26_s(0,:,:), (/1,2/))
 TSOI45_i = conform_dims(dimsizes(TSOI45_s), TSOI45_s(0,:,:), (/1,2/))
 TSOI60_i = conform_dims(dimsizes(TSOI60_s), TSOI60_s(0,:,:), (/1,2/))
 TSOI85_i = conform_dims(dimsizes(TSOI85_s), TSOI85_s(0,:,:), (/1,2/))
 printVarSummary(TSOI85_i)

 dTSOI26 = TSOI26_s - TSOI26_i
 dTSOI26 = where(dTSOI26.lt.0, 0, dTSOI26)           ;prevent NANs from being calculated w/ local cooling
 copy_VarMeta(TSOI26_s, dTSOI26)
 dTSOI26@long_name = "change in soil temperature from initial (rcp26)" 
 printVarSummary(dTSOI26)

 dTSOI45 = TSOI45_s - TSOI45_i
 dTSOI45 = where(dTSOI45.lt.0, 0, dTSOI45)           ;prevent NANs from being calculated w/ local cooling
 copy_VarMeta(TSOI45_s, dTSOI45)
 dTSOI45@long_name = "change in soil temperature from initial (rcp45)"
 printVarSummary(dTSOI45)

 dTSOI60 = TSOI60_s - TSOI60_i
 dTSOI60 = where(dTSOI60.lt.0, 0, dTSOI60)           ;prevent NANs from being calculated w/ local cooling
 copy_VarMeta(TSOI60_s, dTSOI60)
 dTSOI60@long_name = "change in soil temperature from initial (rcp60)"
 printVarSummary(dTSOI60)

 dTSOI85 = TSOI85_s - TSOI85_i
 dTSOI85 = where(dTSOI85.lt.0, 0, dTSOI85)           ;prevent NANs from being calculated w/ local cooling
 copy_VarMeta(TSOI85_s, dTSOI85)
 dTSOI85@long_name = "change in soil temperature from initial (rcp85)"

 ; Calculate average weighted change in TSOI85
 sum_area   = sum(grid_area)
 wgt_area   = grid_area / sum(grid_area)
 wgt_areaT  = conform_dims(dimsizes(dTSOI85),  wgt_area , (/1,2/))                 

 Global_dTSOI26 = dim_sum_n((dTSOI26 * wgt_areaT), (/1,2/))
 Global_dTSOI45 = dim_sum_n((dTSOI45 * wgt_areaT), (/1,2/))
 Global_dTSOI60 = dim_sum_n((dTSOI60 * wgt_areaT), (/1,2/))
 Global_dTSOI85 = dim_sum_n((dTSOI85 * wgt_areaT), (/1,2/))
 printVarSummary(Global_dTSOI85)

;--------------------------------------------------------------
 ; emperical losses of C based on Tom's meta-analysis
 ; Here using simple linear regressions 
 ; SOC[,i] <-  SOC[,i-1]  + (a * SOC[,i-1] * dTSOI[i] ) + b 

 ; parameter estimates for 50%, 5%, and 95% respectively 
 a = (/ -0.0369983, -0.0441048, -0.0281913 /) 
 b = (/  0.0009036,  0.0006128,  0.0011552 /)

 HWSD_SOCb = HWSD_SOC * 1e3 / (1e4 * HWSD_depth)  ; convert obser to g/cm3, 0-30 cm depth
 copy_VarMeta(HWSD_SOC, HWSD_SOCb)
 HWSD_SOCb@units = "gC/cm3, 0-30cm"

 SGRD_SOCb = SGRD_SOC * 1e3 / (1e4 * SGRD_depth)  ; convert obser to g/cm3, 0-15 cm depth
 copy_VarMeta(SGRD_SOC, SGRD_SOCb)
 SGRD_SOCb@units = "gC/cm3, 0-"+SGRD_depth+" cm"

 dims        = (/3,nyear,nlat,nlon/)
 x           = (/0.50,0.05,0.95/)
 HWSD_SOC_o  = conform_dims( dims , HWSD_SOC , (/2,3/))
 HWSD_SOC_26 = conform_dims( dims , HWSD_SOCb, (/2,3/))
 HWSD_SOC_45 = conform_dims( dims , HWSD_SOCb, (/2,3/))
 HWSD_SOC_60 = conform_dims( dims , HWSD_SOCb, (/2,3/))
 HWSD_SOC_85 = conform_dims( dims , HWSD_SOCb, (/2,3/))
 SGRD_SOC_o  = conform_dims( dims , SGRD_SOC , (/2,3/))
 SGRD_SOC_26 = conform_dims( dims , SGRD_SOCb, (/2,3/))
 SGRD_SOC_45 = conform_dims( dims , SGRD_SOCb, (/2,3/))
 SGRD_SOC_60 = conform_dims( dims , SGRD_SOCb, (/2,3/))
 SGRD_SOC_85 = conform_dims( dims , SGRD_SOCb, (/2,3/))
 SGRD_SOC_D  = SGRD_SOC_85(0,:,:,:) 
 
 HWSD_SOC_o!0 = "range"
 HWSD_SOC_o!1 = "year"
 HWSD_SOC_o!2 = "lat"
 HWSD_SOC_o!3 = "lon"

 HWSD_SOC_o&range = x
 HWSD_SOC_o&year  = year
 HWSD_SOC_o&lat   = lat
 HWSD_SOC_o&lon   = lon
 copy_VarAtts( HWSD_SOC,  HWSD_SOC_o )
 copy_VarAtts( HWSD_SOCb, HWSD_SOC_26 )
 copy_VarAtts( HWSD_SOCb, HWSD_SOC_45 )
 copy_VarAtts( HWSD_SOCb, HWSD_SOC_60 )
 copy_VarAtts( HWSD_SOCb, HWSD_SOC_85 )
 copy_VarCoords( HWSD_SOC_o,  HWSD_SOC_26 )
 copy_VarCoords( HWSD_SOC_o,  HWSD_SOC_45 )
 copy_VarCoords( HWSD_SOC_o,  HWSD_SOC_60 )
 copy_VarCoords( HWSD_SOC_o,  HWSD_SOC_85 )

 copy_VarAtts( SGRD_SOC,  SGRD_SOC_o )
 copy_VarAtts( SGRD_SOCb, SGRD_SOC_26 )
 copy_VarAtts( SGRD_SOCb, SGRD_SOC_45 )
 copy_VarAtts( SGRD_SOCb, SGRD_SOC_60 )
 copy_VarAtts( SGRD_SOCb, SGRD_SOC_85 )
 copy_VarCoords( HWSD_SOC_o,  SGRD_SOC_o )
 copy_VarCoords( HWSD_SOC_o,  SGRD_SOC_26 )
 copy_VarCoords( HWSD_SOC_o,  SGRD_SOC_45 )
 copy_VarCoords( HWSD_SOC_o,  SGRD_SOC_60 )
 copy_VarCoords( HWSD_SOC_o,  SGRD_SOC_85 )

 printVarSummary(HWSD_SOC_85)
 Treat26 = dTSOI26(0,:,:)
 Treat45 = dTSOI45(0,:,:)
 Treat60 = dTSOI60(0,:,:)
 Treat85 = dTSOI85(0,:,:)
 do j =0,2
   do i = 1,nyear-1                                
     dT26  = (TSOI26_s(i,:,:)-TSOI26_s(i-1,:,:))
     dT45  = (TSOI45_s(i,:,:)-TSOI45_s(i-1,:,:))
     dT60  = (TSOI60_s(i,:,:)-TSOI60_s(i-1,:,:))
     dT85  = (TSOI85_s(i,:,:)-TSOI85_s(i-1,:,:))

     ; SOC[,i] <-  SOC[,i-1]  + (a * SOC[,i-1] + b ) * dT[i]
     HWSD_SOC_26(j,i,:,:) = HWSD_SOC_26(j,i-1,:,:)+(a(j)*HWSD_SOC_26(j,i-1,:,:) +b(j)) * dT26 
     HWSD_SOC_45(j,i,:,:) = HWSD_SOC_45(j,i-1,:,:)+(a(j)*HWSD_SOC_45(j,i-1,:,:) +b(j)) * dT45
     HWSD_SOC_60(j,i,:,:) = HWSD_SOC_60(j,i-1,:,:)+(a(j)*HWSD_SOC_60(j,i-1,:,:) +b(j)) * dT60
     HWSD_SOC_85(j,i,:,:) = HWSD_SOC_85(j,i-1,:,:)+(a(j)*HWSD_SOC_85(j,i-1,:,:) +b(j)) * dT85

     SGRD_SOC_26(j,i,:,:) = SGRD_SOC_26(j,i-1,:,:)+(a(j)*SGRD_SOC_26(j,i-1,:,:) +b(j)) * dT26
     SGRD_SOC_45(j,i,:,:) = SGRD_SOC_45(j,i-1,:,:)+(a(j)*SGRD_SOC_45(j,i-1,:,:) +b(j)) * dT45
     SGRD_SOC_60(j,i,:,:) = SGRD_SOC_60(j,i-1,:,:)+(a(j)*SGRD_SOC_60(j,i-1,:,:) +b(j)) * dT60
     SGRD_SOC_85(j,i,:,:) = SGRD_SOC_85(j,i-1,:,:)+(a(j)*SGRD_SOC_85(j,i-1,:,:) +b(j)) * dT85

   end do
 end do

  HWSD_SOC_26   = HWSD_SOC_26 * HWSD_depth * 1e4 / 1e3           ;convert from g/cm3 to kg/m2
  dHWSD_SOC26   = HWSD_SOC_26 - HWSD_SOC_o
  HWSD_SOC_45   = HWSD_SOC_45 * HWSD_depth * 1e4 / 1e3           ;convert from g/cm3 to kg/m2
  dHWSD_SOC45   = HWSD_SOC_45 - HWSD_SOC_o
  HWSD_SOC_60   = HWSD_SOC_60 * HWSD_depth * 1e4 / 1e3           ;convert from g/cm3 to kg/m2
  dHWSD_SOC60   = HWSD_SOC_60 - HWSD_SOC_o      
  HWSD_SOC_85   = HWSD_SOC_85 * HWSD_depth * 1e4 / 1e3           ;convert from g/cm3 to kg/m2
  dHWSD_SOC85   = HWSD_SOC_85 - HWSD_SOC_o                       ;calculate SOC losses each year
  HWSD_SOC_26@units = HWSD_SOC@units
  HWSD_SOC_45@units = HWSD_SOC@units
  HWSD_SOC_60@units = HWSD_SOC@units
  HWSD_SOC_85@units = HWSD_SOC@units

  SGRD_SOC_26   = SGRD_SOC_26 * SGRD_depth * 1e4 / 1e3
  dSGRD_SOC26   = SGRD_SOC_26 - SGRD_SOC_o
  SGRD_SOC_45   = SGRD_SOC_45 * SGRD_depth * 1e4 / 1e3
  dSGRD_SOC45   = SGRD_SOC_45 - SGRD_SOC_o
  SGRD_SOC_60   = SGRD_SOC_60 * SGRD_depth * 1e4 / 1e3
  dSGRD_SOC60   = SGRD_SOC_60 - SGRD_SOC_o  
  SGRD_SOC_85   = SGRD_SOC_85 * SGRD_depth * 1e4 / 1e3
  dSGRD_SOC85   = SGRD_SOC_85 - SGRD_SOC_o  
  SGRD_SOC_26@units = SGRD_SOC@units
  SGRD_SOC_45@units = SGRD_SOC@units
  SGRD_SOC_60@units = SGRD_SOC@units
  SGRD_SOC_85@units = SGRD_SOC@units
  SGRD_SOC_D    = SGRD_SOC_D * SGRD_depth * 1e4 / 1e3

  dHWSD26_flag    = where(-dHWSD_SOC26.gt.HWSD_SOC_o, 1., 0. )       ;flag losses > initial pools
  dSGRD26_flag    = where(-dSGRD_SOC26.gt.SGRD_SOC_o, 1., 0. )
  dHWSD45_flag    = where(-dHWSD_SOC45.gt.HWSD_SOC_o, 1., 0. )       ;flag losses > initial pools
  dSGRD45_flag    = where(-dSGRD_SOC45.gt.SGRD_SOC_o, 1., 0. )  
  dHWSD60_flag    = where(-dHWSD_SOC60.gt.HWSD_SOC_o, 1., 0. )       ;flag losses > initial pools
  dSGRD60_flag    = where(-dSGRD_SOC60.gt.SGRD_SOC_o, 1., 0. )
  dHWSD85_flag    = where(-dHWSD_SOC85.gt.HWSD_SOC_o, 1., 0. ) 
  dSGRD85_flag    = where(-dSGRD_SOC85.gt.SGRD_SOC_o, 1., 0. ) 

 rel_dHWSD26      = 100. * dHWSD_SOC26 / HWSD_SOC_o
 rel_dSGRD26      = 100. * dSGRD_SOC26 / SGRD_SOC_o
 rel_dHWSD45      = 100. * dHWSD_SOC45 / HWSD_SOC_o
 rel_dSGRD45      = 100. * dSGRD_SOC45 / SGRD_SOC_o
 rel_dHWSD60      = 100. * dHWSD_SOC60 / HWSD_SOC_o
 rel_dSGRD60      = 100. * dSGRD_SOC60 / SGRD_SOC_o
 rel_dHWSD85      = 100. * dHWSD_SOC85 / HWSD_SOC_o
 rel_dSGRD85      = 100. * dSGRD_SOC85 / SGRD_SOC_o

 print(year(44))
 midC = 44      ; data in 2050:
 
 grid_areaC       = conform_dims(dimsizes(HWSD_SOC_85), (grid_area * 1e6), (/2,3/))      ;Convert to m2
 grid_areaD       = conform_dims(dimsizes(SGRD_SOC_D), (grid_area * 1e6), (/1,2/))       ;Convert to m2

 GlobalHWSD26     = dim_sum_n_Wrap((HWSD_SOC_26 * grid_areaC), (/2,3/)) / 1.e12           ;Pg C globally
 GlobalHWSD26_init= GlobalHWSD26(0,0)
 GlobalHWSD26_mean= GlobalHWSD26(0,nyear-1)
 GlobalHWSD26_min = GlobalHWSD26(1,nyear-1)
 GlobalHWSD26_max = GlobalHWSD26(2,nyear-1)
 dGlobalHWSD26    = GlobalHWSD26 - GlobalHWSD26_init
 print("HWSD26 init = "+GlobalHWSD26_init+", mean = "+GlobalHWSD26_mean+ \
             ", min = "+GlobalHWSD26_min +", max = " +GlobalHWSD26_max)

 GlobalHWSD45     = dim_sum_n_Wrap((HWSD_SOC_45 * grid_areaC), (/2,3/)) / 1.e12           ;Pg C globally
 GlobalHWSD45_init= GlobalHWSD45(0,0)
 GlobalHWSD45_mean= GlobalHWSD45(0,nyear-1)
 GlobalHWSD45_min = GlobalHWSD45(1,nyear-1)
 GlobalHWSD45_max = GlobalHWSD45(2,nyear-1)
 dGlobalHWSD45    = GlobalHWSD45 - GlobalHWSD45_init
 print("HWSD45 init = "+GlobalHWSD45_init+", mean = "+GlobalHWSD45_mean+ \
             ", min = "+GlobalHWSD45_min +", max = " +GlobalHWSD45_max)

 GlobalHWSD60     = dim_sum_n_Wrap((HWSD_SOC_60 * grid_areaC), (/2,3/)) / 1.e12           ;Pg C globally
 GlobalHWSD60_init= GlobalHWSD60(0,0)
 GlobalHWSD60_mean= GlobalHWSD60(0,nyear-1)
 GlobalHWSD60_min = GlobalHWSD60(1,nyear-1)
 GlobalHWSD60_max = GlobalHWSD60(2,nyear-1)
 dGlobalHWSD60    = GlobalHWSD60 - GlobalHWSD60_init
 print("HWSD60 init = "+GlobalHWSD60_init+", mean = "+GlobalHWSD60_mean+ \
             ", min = "+GlobalHWSD60_min +", max = " +GlobalHWSD60_max)

 GlobalHWSD85     = dim_sum_n_Wrap((HWSD_SOC_85 * grid_areaC), (/2,3/)) / 1.e12           ;Pg C globally
 GlobalHWSD85_init= GlobalHWSD85(0,0) 
 GlobalHWSD85_mean= GlobalHWSD85(0,nyear-1) 
 GlobalHWSD85_min = GlobalHWSD85(1,nyear-1)
 GlobalHWSD85_max = GlobalHWSD85(2,nyear-1)
 dGlobalHWSD85    = GlobalHWSD85 - GlobalHWSD85_init
 print("HWSD85 init = "+GlobalHWSD85_init+", mean = "+GlobalHWSD85_mean+ \
             ", min = "+GlobalHWSD85_min +", max = " +GlobalHWSD85_max)

 GlobalSGRD26     = dim_sum_n_Wrap((SGRD_SOC_26 * grid_areaC), (/2,3/)) / 1.e12        ;Pg C globally
 GlobalSGRD26_init= GlobalSGRD26(0,0)
 GlobalSGRD26_2050= GlobalSGRD26(:,midC)
 GlobalSGRD26_2100= GlobalSGRD26(:,nyear-1)
 dGlobalSGRD26    = GlobalSGRD26 - GlobalSGRD26_init
 print("SoilGrids 26 init = "+GlobalSGRD26_init)
 print("SoilGrids 26 2050 = "+GlobalSGRD26_2050)
 print("SoilGrids 26 2100 = "+GlobalSGRD26_2100)

 GlobalSGRD45     = dim_sum_n_Wrap((SGRD_SOC_45 * grid_areaC), (/2,3/)) / 1.e12        ;Pg C globally
 GlobalSGRD45_init= GlobalSGRD45(0,0)
 GlobalSGRD45_2050= GlobalSGRD45(:,midC)
 GlobalSGRD45_2100= GlobalSGRD45(:,nyear-1)
 dGlobalSGRD45    = GlobalSGRD45 - GlobalSGRD45_init
 print("SoilGrids 45 init = "+GlobalSGRD45_init)
 print("SoilGrids 45 2050 = "+GlobalSGRD45_2050)
 print("SoilGrids 45 2100 = "+GlobalSGRD45_2100)

 GlobalSGRD60     = dim_sum_n_Wrap((SGRD_SOC_60 * grid_areaC), (/2,3/)) / 1.e12        ;Pg C globally
 GlobalSGRD60_init= GlobalSGRD60(0,0) 
 GlobalSGRD60_2050= GlobalSGRD60(:,midC)
 GlobalSGRD60_2100= GlobalSGRD60(:,nyear-1)
 dGlobalSGRD60    = GlobalSGRD60 - GlobalSGRD60_init 
 print("SoilGrids 60 init = "+GlobalSGRD60_init)
 print("SoilGrids 60 2050 = "+GlobalSGRD60_2050)
 print("SoilGrids 60 2100 = "+GlobalSGRD60_2100)

 GlobalSGRD85     = dim_sum_n_Wrap((SGRD_SOC_85 * grid_areaC), (/2,3/)) / 1.e12        ;Pg C globally
 GlobalSGRD_D     = dim_sum_n_Wrap((SGRD_SOC_D * grid_areaD), (/1,2/)) / 1.e12         
 GlobalSGRD85_init= GlobalSGRD85(0,0) 
 GlobalSGRD85_2050= GlobalSGRD85(:,midC)
 GlobalSGRD85_2100= GlobalSGRD85(:,nyear-1)
 dGlobalSGRD85    = GlobalSGRD85 - GlobalSGRD85_init
 print("SoilGrids 85 init = "+GlobalSGRD85_init)
 print("SoilGrids 85 2050 = "+GlobalSGRD85_2050)
 print("SoilGrids 85 2100 = "+GlobalSGRD85_2100)
 print("SoilGrids final(mean linear) = "+ GlobalSGRD_D(nyear-1))

 copy_VarMeta(HWSD_SOC_o, dHWSD_SOC26  )
 copy_VarMeta(HWSD_SOC_o, dHWSD26_flag )
 copy_VarMeta(HWSD_SOC_o, rel_dHWSD26  )
 dHWSD_SOC26@long_name  = "absolute changes to area weighted SOC, rcp26"
 dHWSD_SOC26@units      = HWSD_SOC@units
 dHWSD26_flag@long_name = "losses > initial pools"
 rel_dHWSD26@long_name  = "relative changes to area weighted SOC, rcp26"
 rel_dHWSD26@units      = "percent"

 copy_VarMeta(HWSD_SOC_o, dHWSD_SOC45  )
 copy_VarMeta(HWSD_SOC_o, dHWSD45_flag )
 copy_VarMeta(HWSD_SOC_o, rel_dHWSD45  )
 dHWSD_SOC45@long_name  = "absolute changes to area weighted SOC, rcp45"
 dHWSD_SOC45@units      = HWSD_SOC@units
 dHWSD45_flag@long_name = "losses > initial pools"
 rel_dHWSD45@long_name  = "relative changes to area weighted SOC, rcp45"
 rel_dHWSD45@units      = "percent"

 copy_VarMeta(HWSD_SOC_o, dHWSD_SOC60  )
 copy_VarMeta(HWSD_SOC_o, dHWSD60_flag )
 copy_VarMeta(HWSD_SOC_o, rel_dHWSD60  )
 dHWSD_SOC60@long_name  = "absolute changes to area weighted SOC, rcp60"
 dHWSD_SOC60@units      = HWSD_SOC@units
 dHWSD60_flag@long_name = "losses > initial pools"
 rel_dHWSD60@long_name  = "relative changes to area weighted SOC, rcp60"
 rel_dHWSD60@units      = "percent"

 copy_VarMeta(HWSD_SOC_o, dHWSD_SOC85  )
 copy_VarMeta(HWSD_SOC_o, dHWSD85_flag )
 copy_VarMeta(HWSD_SOC_o, rel_dHWSD85  )
 dHWSD_SOC85@long_name  = "absolute changes to area weighted SOC, rcp85"
 dHWSD_SOC85@units      = HWSD_SOC@units
 dHWSD85_flag@long_name = "losses > initial pools"
 rel_dHWSD85@long_name  = "relative changes to area weighted SOC, rcp85"
 rel_dHWSD85@units      = "percent"

 copy_VarMeta(SGRD_SOC_o, dSGRD_SOC26  )
 copy_VarMeta(SGRD_SOC_o, dSGRD26_flag )
 copy_VarMeta(SGRD_SOC_o, rel_dSGRD26  )
 dSGRD_SOC26@long_name  = "absolute changes to area weighted SOC, rcp26"
 dSGRD_SOC26@units      = SGRD_SOC@units
 dSGRD26_flag@long_name = "losses > initial pools"
 rel_dSGRD26@long_name  = "relative changes to area weighted SOC, rcp26"
 rel_dSGRD26@units      = "percent"

 copy_VarMeta(SGRD_SOC_o, dSGRD_SOC45  )
 copy_VarMeta(SGRD_SOC_o, dSGRD45_flag )
 copy_VarMeta(SGRD_SOC_o, rel_dSGRD45  )
 dSGRD_SOC45@long_name  = "absolute changes to area weighted SOC, rcp45"
 dSGRD_SOC45@units      = SGRD_SOC@units
 dSGRD45_flag@long_name = "losses > initial pools"
 rel_dSGRD45@long_name  = "relative changes to area weighted SOC, rcp45"
 rel_dSGRD45@units      = "percent"

 copy_VarMeta(SGRD_SOC_o, dSGRD_SOC60  )
 copy_VarMeta(SGRD_SOC_o, dSGRD60_flag )
 copy_VarMeta(SGRD_SOC_o, rel_dSGRD60  )
 dSGRD_SOC60@long_name  = "absolute changes to area weighted SOC, rcp60"
 dSGRD_SOC60@units      = SGRD_SOC@units
 dSGRD60_flag@long_name = "losses > initial pools"
 rel_dSGRD60@long_name  = "relative changes to area weighted SOC, rcp60"
 rel_dSGRD60@units      = "percent"

 copy_VarMeta(SGRD_SOC_o, dSGRD_SOC85  )
 copy_VarMeta(SGRD_SOC_o, dSGRD85_flag )
 copy_VarMeta(SGRD_SOC_o, rel_dSGRD85  )
 dSGRD_SOC85@long_name  = "absolute changes to area weighted SOC, rcp85"
 dSGRD_SOC85@units      = SGRD_SOC@units
 dSGRD85_flag@long_name = "losses > initial pools"
 rel_dSGRD85@long_name  = "relative changes to area weighted SOC, rcp85"
 rel_dSGRD85@units      = "percent"

 ;--------------------------------------------------------------
 ; cacluate carbon-climate feedback gamma = flux(to atmosphere) / Tchange / year
 printVarSummary(GlobalSGRD85)
 printVarSummary(Global_dTSOI85)
 delete([/dT26,dT45,dT60,dT85/]) 
 dT26 = conform_dims(dimsizes(GlobalSGRD26), Global_dTSOI26, 1)
 dT26 = where(dT26 .eq. 0, dT26@_FillValue, dT26)
 dC26  = GlobalSGRD26

 dT45 = conform_dims(dimsizes(GlobalSGRD45), Global_dTSOI45, 1)
 dT45 = where(dT45 .eq. 0, dT45@_FillValue, dT45)
 dC45  = GlobalSGRD45

 dT60 = conform_dims(dimsizes(GlobalSGRD60), Global_dTSOI60, 1)
 dT60 = where(dT60 .eq. 0, dT60@_FillValue, dT60)
 dC60  = GlobalSGRD60  

 dT85 = conform_dims(dimsizes(GlobalSGRD85), Global_dTSOI85, 1)
 dT85 = where(dT85 .eq. 0, dT85@_FillValue, dT85)
 dC85  = GlobalSGRD85

 do i = 1,nyear-1
    dC26(:,i)  = GlobalSGRD26(:,i) - GlobalSGRD26(:,i-1)
    dC45(:,i)  = GlobalSGRD45(:,i) - GlobalSGRD45(:,i-1)
    dC60(:,i)  = GlobalSGRD60(:,i) - GlobalSGRD60(:,i-1)
    dC85(:,i)  = GlobalSGRD85(:,i) - GlobalSGRD85(:,i-1)
 end do

 gSOIL26  = dC26 / dT26
 gSOIL45  = dC45 / dT45
 gSOIL60  = dC60 / dT60
 gSOIL85  = dC85 / dT85

 print(dT85(0,0:20))

 diffSGRD_26_85  = SGRD_SOC_26 - SGRD_SOC_85
 copy_VarMeta(SGRD_SOC_26, diffSGRD_26_85)
 diffSGRD_26_85@long_name = "SGRD_SOC rcp26 - rcp85"

 diffTreat26_85 = Treat26 - Treat85
  
 ;---------------------------------------------------------------
 ; write .nc file of output
 ;---------------------------------------------------------------

 pathout      = (/"/project/tss/wwieder/soilCN/global_run/warming/"/) 
 fout         = pathout + "cmip5_"+title(EXP) + "_SOC_losses_KTB.nc"
 system("/bin/rm -f "+fout )   ; remove any pre-existing file
 ncdf = addfile(fout ,"c")  ; open output netCDF file
 filedimdef(ncdf,"time",-1,True)

 ncdf->TSOI26       = TSOI26
 ncdf->TSOI26_s     = TSOI26_s
 ncdf->dTSOI26      = dTSOI26
 ncdf->Treat26      = Treat26
; ncdf->HWSD_SOC_26  = HWSD_SOC_26
; ncdf->dHWSD_SOC26  = dHWSD_SOC26
; ncdf->dHWSD26_flag = dHWSD26_flag
; ncdf->rel_dHWSD26  = rel_dHWSD26
 ncdf->SGRD_SOC_26  = SGRD_SOC_26
 ncdf->dSGRD_SOC26  = dSGRD_SOC26
; ncdf->dSGRD26_flag = dSGRD26_flag
; ncdf->rel_dSGRD26  = rel_dSGRD26
 ncdf->diffSGRD_26_85 = diffSGRD_26_85 
 ncdf->diffTreat26_85 = diffTreat26_85

 ncdf->TSOI45       = TSOI45
 ncdf->TSOI45_s     = TSOI45_s
 ncdf->dTSOI45      = dTSOI45
 ncdf->Treat45      = Treat45
; ncdf->HWSD_SOC_45  = HWSD_SOC_45
; ncdf->dHWSD_SOC45  = dHWSD_SOC45
; ncdf->dHWSD45_flag = dHWSD45_flag
; ncdf->rel_dHWSD45  = rel_dHWSD45
 ncdf->SGRD_SOC_45  = SGRD_SOC_45
 ncdf->dSGRD_SOC45  = dSGRD_SOC45
; ncdf->dSGRD45_flag = dSGRD45_flag
; ncdf->rel_dSGRD45  = rel_dSGRD45

 ncdf->TSOI60       = TSOI60
 ncdf->TSOI60_s     = TSOI60_s
 ncdf->dTSOI60      = dTSOI60
 ncdf->Treat60      = Treat60
; ncdf->HWSD_SOC_60  = HWSD_SOC_60
; ncdf->dHWSD_SOC60  = dHWSD_SOC60
; ncdf->dHWSD60_flag = dHWSD60_flag
; ncdf->rel_dHWSD60  = rel_dHWSD60
 ncdf->SGRD_SOC_60  = SGRD_SOC_60
 ncdf->dSGRD_SOC60  = dSGRD_SOC60
; ncdf->dSGRD60_flag = dSGRD60_flag
; ncdf->rel_dSGRD60  = rel_dSGRD60
 
 ncdf->TSOI85       = TSOI85
 ncdf->TSOI85_s     = TSOI85_s
 ncdf->dTSOI85      = dTSOI85 
 ncdf->Treat85      = Treat85
; ncdf->HWSD_SOCb    = HWSD_SOCb
; ncdf->HWSD_SOC_85  = HWSD_SOC_85
; ncdf->dHWSD_SOC85  = dHWSD_SOC85 
; ncdf->dHWSD85_flag = dHWSD85_flag
; ncdf->rel_dHWSD85  = rel_dHWSD85
 ncdf->SGRD_SOCb    = SGRD_SOCb
 ncdf->SGRD_SOC_85  = SGRD_SOC_85
 ncdf->dSGRD_SOC85  = dSGRD_SOC85
 ncdf->dSGRD85_flag = dSGRD85_flag
 ncdf->rel_dSGRD85  = rel_dSGRD85

 ncdf->HWSD_SOCb    = HWSD_SOCb
 ncdf->SGRD_SOCb    = SGRD_SOCb

 print("wrote "+fout)

 ;---------------------------------------------------------------
 ; plot SOC losses
 ;---------------------------------------------------------------

  fout = pathout + "SOC_losses_KTB"   
  wks  = gsn_open_wks("ps" , fout); open a X11 or ps file

  plot = new(2,graphic)                             ; create a plot array

  res                   = True                      ; plot mods desired
  res@gsnFrame          = False           ; don't advance fram
  res@gsnDraw           = False           ; don't draw
  res@xyDashPatterns    = (/ 2,0,0,0,0 /)
  res@xyLineThicknesses = (/ 1,5,5,5,5 /)               ; make lines thicker
  res@xyLineColors      = (/ "black","red","yellow","green","blue" /)    ; change line color
  res@trXMinF  = year(1)
  res@trXMaxF  = year(nyear-1)
  res@trYMinF = min(GlobalSGRD85)*0.95
  res@trYMaxF = max(GlobalSGRD85)*1.05
  res@tiYAxisString     = "Global soil C (Pg C, 0-15 cm)"
  res@tiXAxisString     = "Year"
  res@vpXF              = 0.1                 ; make plot bigger
; res@vpYF              = 0.9
  res@vpWidthF          = 0.8
; res@vpHeightF         = 0.8

;-- Legend resources--
  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Left"              ; Change location of
  res@pmLegendParallelPosF   = 0.20                ; move units down
  res@pmLegendOrthogonalPosF = -0.50               ; move units right
  res@lgPerimOn              = False               ; turn off box around

  res@pmLegendWidthF         = 0.18                ; Change width and
  res@pmLegendHeightF        = 0.20                ; height of legend.
  res@lgLabelFontHeightF     = .022                ; change font height
;  res@lgTitleOn              = True                ; turn on legend title
;  res@lgTitleString          = "Example"           ; create legend title
;  res@lgTitleFontHeightF     = .025                ; font of legend title
  res@xyExplicitLegendLabels = (/""," RCP85"," RCP60"," RCP45"," RCP26"/)
;  res@lgItemOrder            = (/4,3,2,0,1/)               ; reverse legend


  ynew       = new(dimsizes(year),typeof(GlobalSGRD85))
  ynew(:)    = 0.
  xplot      = new((/5,nyear/), typeof(nyear))
  yplot      = new((/5,nyear/), typeof(GlobalSGRD85))
  xplot(0,:) = year 
  xplot(1,:) = year  
  xplot(2,:) = year
  xplot(3,:) = year 
  xplot(4,:) = year 

  yplot(0,:) = GlobalSGRD85(0,0)
  yplot(1,:) = GlobalSGRD85(0,:)
  yplot(2,:) = GlobalSGRD60(0,:)
  yplot(3,:) = GlobalSGRD45(0,:)
  yplot(4,:) = GlobalSGRD26(0,:)

  plot(0)  = gsn_csm_xy (wks,xplot,yplot,res) ; create plot

;************************************************
; define a polygon centered with width of 2 sigma
;************************************************
 xp    = new((/2*nyear/), typeof(year))
 xp(0:nyear-1) = year
 xp(nyear: )   = year(::-1)
 yp    = new((/2*nyear/), typeof(GlobalSGRD85))
 yp(0:nyear-1) = (/GlobalSGRD85(1,:)/)
 temp          = (/GlobalSGRD85(2,:)/)
 yp(nyear: )   = temp(::-1)
 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
 gsres@gsFillColor       = "red"                      ; color chosen
 gsres@gsFillOpacityF    =  0.2                        
 dummy85 = gsn_add_polygon (wks,plot(0),xp,yp,gsres)

 yp(0:nyear-1) = (/GlobalSGRD26(1,:)/)
 temp          = (/GlobalSGRD26(2,:)/)
 yp(nyear: )   = temp(::-1)
 gsres@gsFillColor       = "blue"                      ; color chosen
 dummy26 = gsn_add_polygon (wks,plot(0),xp,yp,gsres)     

 yp(0:nyear-1) = (/GlobalSGRD45(1,:)/)
 temp          = (/GlobalSGRD45(2,:)/)
 yp(nyear: )   = temp(::-1)
 gsres@gsFillColor       = "green"                      ; color chosen
 dummy45 = gsn_add_polygon (wks,plot(0),xp,yp,gsres)

 yp(0:nyear-1) = (/GlobalSGRD60(1,:)/)
 temp          = (/GlobalSGRD60(2,:)/)
 yp(nyear: )   = temp(::-1)
 gsres@gsFillColor       = "yellow"                      ; color chosen
 dummy60 = gsn_add_polygon (wks,plot(0),xp,yp,gsres)


 delete(res)

;************************************************
; plot maps of SOC change 
;************************************************
 res                        = True
 res@gsnDraw                = False
 res@gsnFrame               = False
 res@cnSmoothingOn          = False
 res@mpProjection           = "Robinson"
 res@mpOutlineOn            = True
 res@lbOrientation          = "Horizontal"
 res@mpPerimOn              = False
 res@mpGridAndLimbOn        = True         ;draws Limb line around plots, w/o lat-lon lines
 res@mpGridLatSpacingF      = 180
 res@mpGridLonSpacingF      = 180
 res@mpGridLineThicknessF   = 0.
 res@mpGridLineColor        = "transparent"
 res@mpGridMaskMode         = "MaskLand"

;  gsn_define_colormap(wks,"MPL_gist_earth")         ; change colormap 
;  gsn_define_colormap(wks,"MPL_terrain") 
;  gsn_define_colormap(wks,"GMT_seis")
;  res@gsnSpreadColors     = True                   ; use full colormap
;  res@gsnSpreadColorStart = 2                      ; start with last color
;  res@gsnSpreadColorEnd   = 160                    ; end with a beginning color

   gsn_define_colormap(wks,"BlWhRe")  
   res@gsnSpreadColors     = True                   ; use full colormap
   res@gsnSpreadColorEnd   = 65                      ; start with last color
   gsn_reverse_colormap(wks)                        ; reverse colormap

;  colors = (/"white","black","Tan4","Tan3","Tan2","Tan1","moccasin","darkolivegreen3","darkolivegreen4" /)
;  gsn_define_colormap(wks,colors)            ; change colormap

  res@gsnLeftString        = ""   
  res@gsnRightString       = ""  
  res@cnFillOn             = True
  res@cnLinesOn            = False         ; Turn lines off
  res@cnLineLabelsOn       = False         ; Turn labels off
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -3.5
  res@cnMaxLevelValF       =  0.5
  res@cnLevelSpacingF      =  0.5
  res@lbLabelFontHeightF   =  0.025                ; make labels larger
  res@lbTitleOn            = True                  ; turn on title
  res@lbTitleString        = "kg C m~S~-2~N "
  res@lbTitleFontHeightF   = .030                  ; make title smaller
  res@pmLabelBarOrthogonalPosF = .10               ; move whole thing down

  res@vpXF            = 0.1                 ; make plot bigger
  res@vpYF            = 0.9
  res@vpWidthF        = 0.8
  res@vpHeightF       = 0.8
  plot(1)    = gsn_csm_contour_map(wks,dSGRD_SOC85(0,nyear-1,:,:),res)

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
;  resP@gsnMaximize      = True                   ; maximize plots
  resP@gsnPanelFigureStrings= (/"A", "B"/)
  resP@gsnPanelFigureStringsPerimOn = False      ; turns off box around fig labels
  resP@wkBackgroundOpacityF = 0.
  resP@gsnStringFont  = "helvetica-bold"
  resP@amJust         = "TopLeft" 

  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
  frame(wks)
  print("wrote "+fout+".ps")
  delete([/plot, res, resP, wks,fout/])

;************************************************
; plot maps SOC change
;************************************************
  fout  = pathout + "dSOC_KTB"
  wks2  = gsn_open_wks("ps" , fout); open a X11 or ps file
  gsn_define_colormap(wks2,"BlWhRe")
  gsn_reverse_colormap(wks2)                        ; reverse colormap

  plot = new(4,graphic)                             ; create a plot array
  res                        = True
  res@gsnDraw                = False
  res@gsnFrame               = False
  res@gsnLeftString          = ""  
  res@gsnRightString         = ""
  res@cnSmoothingOn          = False
  res@mpProjection           = "Robinson"
  res@mpOutlineOn            = True
  res@lbOrientation          = "Horizontal"
  res@mpPerimOn              = False
  res@mpGridAndLimbOn        = True         ;draws Limb line around plots, w/o lat-lon lines
  res@mpGridLatSpacingF      = 180
  res@mpGridLonSpacingF      = 180
  res@mpGridLineThicknessF   = 0.
  res@mpGridLineColor        = "transparent"
  res@mpGridMaskMode         = "MaskLand"
  res@cnFillOn             = True
  res@gsnSpreadColors      = True             ; use full colormap
  res@cnLinesOn            = False            ; Turn lines off

  res@lbLabelBarOn         = False            ; turn on/off individual cb's
  res@cnLineLabelsOn       = False            ; Turn labels off
  res@gsnSpreadColorEnd    = 65               ; start with last color

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -3.5 
  res@cnMaxLevelValF       =  0.5
  res@cnLevelSpacingF      =  0.5

  plot(0)    = gsn_csm_contour_map(wks2,dSGRD_SOC26(0,94,:,:),res)
  plot(1)    = gsn_csm_contour_map(wks2,dSGRD_SOC45(0,94,:,:),res)
  plot(2)    = gsn_csm_contour_map(wks2,dSGRD_SOC60(0,94,:,:),res)
  plot(3)    = gsn_csm_contour_map(wks2,dSGRD_SOC85(0,94,:,:),res)

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@amJust           = "BottomLeft"
  resP@gsnPanelLabelBar = True                     ; add common colorbar
  resP@gsnPanelFigureStrings= (/" A RCP2.6, 2100", \
                                " B RCP4.5, 2100", \
                                " C RCP6.0, 2100", \
                                " D RCP8.5, 2100" /)

  resP@lbLabelFontHeightF   =  0.025                ; make labels larger
  resP@lbTitleOn            = True                  ; turn on title
  resP@lbTitleString        = "kg C m~S~-2~N "
  resP@lbTitleFontHeightF   = .030                  ; make title smaller

  gsn_panel(wks2,plot,(/2,2/),resP)               ; now draw as one plot
  frame(wks2)
  print("wrote "+fout+".ps")
  delete([/plot, res, resP, wks2, fout/])

;************************************************
; plot maps TEMP change
;************************************************
  fout  = pathout + "dTSOI_KTB"
  wks2  = gsn_open_wks("ps" , fout); open a X11 or ps file
  gsn_define_colormap(wks2,"WhBlGrYeRe")               ; change colormap

  plot = new(8,graphic)                             ; create a plot array
  res                        = True
  res@gsnDraw                = False
  res@gsnFrame               = False
  res@gsnLeftString          = ""
  res@gsnRightString         = ""
  res@cnSmoothingOn          = False
  res@mpProjection           = "Robinson"
  res@mpOutlineOn            = True
  res@lbOrientation          = "Horizontal"
  res@mpPerimOn              = False
  res@mpGridAndLimbOn        = True         ;draws Limb line around plots, w/o lat-lon lines
  res@mpGridLatSpacingF      = 180
  res@mpGridLonSpacingF      = 180
  res@mpGridLineThicknessF   = 0.
  res@mpGridLineColor        = "transparent"
  res@mpGridMaskMode         = "MaskLand"

  res@cnFillOn             = True
  res@gsnSpreadColors      = True              ; use full colormap
  res@cnLinesOn            = False            ; Turn lines off

  res@lbLabelBarOn         = False            ; turn on/off individual cb's
  res@cnLineLabelsOn       = False            ; Turn labels off
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       =  0.
  res@cnMaxLevelValF       =  5. 
  res@cnLevelSpacingF      =  0.5

  plot(0)    = gsn_csm_contour_map(wks2,dTSOI26(44,:,:),res)
  plot(1)    = gsn_csm_contour_map(wks2,dTSOI26(94,:,:),res)
  plot(2)    = gsn_csm_contour_map(wks2,dTSOI45(44,:,:),res)
  plot(3)    = gsn_csm_contour_map(wks2,dTSOI45(94,:,:),res)
  plot(4)    = gsn_csm_contour_map(wks2,dTSOI60(44,:,:),res)
  plot(5)    = gsn_csm_contour_map(wks2,dTSOI60(94,:,:),res)
  plot(6)    = gsn_csm_contour_map(wks2,dTSOI85(44,:,:),res)
  plot(7)    = gsn_csm_contour_map(wks2,dTSOI85(94,:,:),res)

  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@amJust           = "BottomLeft"
  resP@gsnPanelLabelBar = True                ; add common colorbar
  resP@gsnPanelFigureStrings= (/" A RCP2.6, 2050"," B RCP2.6, 2100", \
                                " C RCP4.5, 2050"," D RCP4.5, 2100", \ 
                                " E RCP6.0, 2050"," F RCP6.0, 2100", \ 
                                " G RCP8.5, 2050"," H RCP8.5, 2100" /)

  resP@lbLabelFontHeightF   =  0.025                ; make labels larger
  resP@lbTitleOn            = True                  ; turn on title
  resP@lbTitleString        =  "~S~o~N~C"
  resP@lbTitleFontHeightF   = .030                  ; make title smaller


  gsn_panel(wks2,plot,(/4,2/),resP)               ; now draw as one plot
  frame(wks2)
  delete(res)

  print("wrote "+fout+".ps")
exit




;************************************************
;   plot change in pools & attribution
;************************************************
  wks  = gsn_open_wks("ps"  , fout); open a X11 or ps file

;----------plot zonal mean for each 'source' of N---------------
 res               = True
 res@trXMinF       = -1.                  ; min value on x-axis
 res@trXMaxF       =  1.                  ; max value on x-axis
 res@xyDashPattern = (/0,0,0,0,0,1/)      ; Make curves all solid, except error
 res@xyLineColors     = (/"black","orange","green","brown","blue","red"/)
 res@xyLineThicknessF = 3.

;-- Legend resources--
  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Bottom"            ; Change location of
  res@pmLegendParallelPosF   = 0.22                ; move units right
  res@pmLegendOrthogonalPosF =  -0.4                ; move units down
  res@lgPerimOn              = False               ; turn off box around

  res@pmLegendWidthF         = 0.15                ; Change width and
  res@pmLegendHeightF        = 0.18                ; height of legend.
  res@lgLabelFontHeightF     = .018                ; change font height
;  res@lgTitleOn              = True                ; turn on legend title
;  res@lgTitleString          = "Example"           ; create legend title
;  res@lgTitleFontHeightF     = .025                ; font of legend title
  res@xyExplicitLegendLabels = (/" New N"," Redistribution", " Veg C:N", \
                                 " Soil C:N"," Interaction", " Error" /)     ; explicit labels
  res@lgItemOrder            = (/5,4,3,2,1,0/)               ; reverse legend

  res@tiMainString = "Relative change 2005, "+title(EXP)      
  plot=gsn_csm_xy(wks,zonal_mean(0,:,:),lat,res)       ; plot zonal ave

if (EXP .ne. 2) then
  res@tiMainString = "Relative change 2100, "+title(EXP)  
  plot=gsn_csm_xy(wks,zonal_mean(1,:,:),lat,res)       ; plot zonal ave
end if
  delete([/plot, res/])



;----------plot total change in C & 'source' of N---------------

  plot = new(6,graphic)                             ; create a plot array

; gsn_define_colormap(wks,"WhViBlGrYeOrRe")  ; choose color map
  gsn_define_colormap(wks,"nrl_sirkes")
 
 res                     = True            ; plot mods desired
 res@gsnMaximize         = True
 res@gsnPaperOrientation = "portrait"
 res@gsnDraw             = False           ; don't draw
 res@gsnFrame            = False           ; don't advance frame
 res@cnInfoLabelOn       = False           ; turn off cn info label
 res@cnFillOn            = True            ; turn on color
 res@cnLinesOn           = False           ; turn off contour lines
 res@cnLineLabelsOn      = False           ; turn off contour labels
 res@gsnSpreadColors     = True            ; spread out color table
 res@lbLabelBarOn        = True            ; turn on/off individual cb's
 res@lbOrientation       =   "vertical"    ; vertical label bar
 res@pmLabelBarOrthogonalPosF = -0.01      ; move label bar closer to axis

  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       = -2500.
  res@cnMaxLevelValF       = 2500.
  res@cnLevelSpacingF      = 250.

 res@vpXF             = 0.12                ; default is 0.2 (aspect ratio)
 res@vpYF             = 0.8                 ; default is 0.8
 res@vpHeightF        = 0.4                 ; default is 0.6
 res@vpWidthF         = 0.75                ; default is 0.6

plot(0) = gsn_csm_contour_map_ce(wks,  d_TOTC(0,:,:), res)

  res@cnMinLevelValF       = -2500.
  res@cnMaxLevelValF       = 2500.
  res@cnLevelSpacingF      = 250.
plot(1) = gsn_csm_contour_map_ce(wks,  d_TOTC(1,:,:), res)

 resP                  = True                   ; modify the panel plot
 resP@gsnFrame         = False                  ; don't advance panel plot
 resP@amJust           = "BottomLeft"           
 resP@gsnPanelLabelBar = False               ; add common colorbar
 resP@gsnPanelFigureStrings= (/"(a) 2005-1900 "+title(EXP),"(b) 2100-2005 "+title(EXP)/)
                               
 gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
 frame(wks)
 delete(resP)

 res@lbLabelBarOn        = False           ; turn on/off individual cb's

  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       = -1. 
  res@cnMaxLevelValF       = 1.
  res@cnLevelSpacingF      = 0.2

 plot(0) = gsn_csm_contour_map_ce(wks,  RdCTdN(0,:,:), res)
 plot(1) = gsn_csm_contour_map_ce(wks,  RdCTdf(0,:,:), res)
 plot(2) = gsn_csm_contour_map_ce(wks, RdCTdqV(0,:,:), res)
 plot(3) = gsn_csm_contour_map_ce(wks, RdCTdqS(0,:,:), res)
 plot(4) = gsn_csm_contour_map_ce(wks, Ri_TOT(0,:,:), res)
 plot(5) = gsn_csm_contour_map_ce(wks,  Rerror(0,:,:), res)

 resP                  = True                   ; modify the panel plot
 resP@gsnFrame         = False                  ; don't advance panel plot
 resP@amJust           = "BottomLeft"
 resP@gsnPanelLabelBar = True                ; add common colorbar
 resP@txString = "Relative change 2005, "+title(EXP)
 resP@gsnPanelFigureStrings= (/"(a) Total N" ,"(b) Redistribution","(c) Veg C:N", \
                               "(d) Soil C:N","(e) Interactions"  ,"(f) Error"/)

 gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot
 frame(wks)

if (EXP .ne. 2) then
  plot(0) = gsn_csm_contour_map_ce(wks,  RdCTdN(1,:,:), res)
  plot(1) = gsn_csm_contour_map_ce(wks,  RdCTdf(1,:,:), res)
  plot(2) = gsn_csm_contour_map_ce(wks, RdCTdqV(1,:,:), res)
  plot(3) = gsn_csm_contour_map_ce(wks, RdCTdqS(1,:,:), res)
  plot(4) = gsn_csm_contour_map_ce(wks, Ri_TOT(1,:,:), res)
  plot(5) = gsn_csm_contour_map_ce(wks,  Rerror(1,:,:), res)

  resP@txString = "Relative change 2100, "+title(EXP)
  gsn_panel(wks,plot,(/3,2/),resP)               ; now draw as one plot
  frame(wks)
end if

  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       = -10.
  res@cnMaxLevelValF       = 10.
  res@cnLevelSpacingF      = 2.

 plot(0) = gsn_csm_contour_map_ce(wks,  difTOTC(0,:,:), res)
 plot(1) = gsn_csm_contour_map_ce(wks,  difTOTN(0,:,:), res)

 resP@gsnPanelFigureStrings= (/"(a) dif TOTC" ,"(b) dif TOTN","(c) Veg C:N", \
                               "(d) Soil C:N","(e) Interactions"  ,"(f) Error"/)

 gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
 frame(wks)

  delete([/plot, resP/])

 ;===================================================
 print("wrote " + fout)
 print ("=========================================")
 print ("Finish Time: "+systemfunc("date") )
 print ("=========================================")

 
end
