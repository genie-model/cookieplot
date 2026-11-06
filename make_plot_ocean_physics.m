
% *********************************************************************** %
% *** INITIALIZE ******************************************************** %
% *********************************************************************** %
%
% *** initialize MATLAB environment ************************************* %
% 
clear;
close all;
% add paths
addpath('y:\___git\cookieplot');
%
% *********************************************************************** %

% *********************************************************************** %
% *** plot_timeseries TEST ********************************************** %
% *********************************************************************** %
%
% *** set common parameters ********************************************* %
% 
% NOTE: these parameter values, unless assigned a different value later, 
%       will be used for all plots
% NOTE: before setting up each individual plot, 
%       the contents of structure s are reset by assigning: s=s_cmn; 
%       so options set for earlier plots do not impact subsequent plots
% x-axis (time) range
s_cmn.tmin_kyr=0.0; s_cmn.tmax_kyr=10.0;
% format of the samed plot (or e.g., 'jpg', or 'pdf')
s_cmn.plotformat='png'; 
%
% *** climate (+ age) only ********************************************** %
%
% NOTE: if you are using the same experiment for all plots,
%       you only need specify the experiment name just once
% (1) basic plot
% specify the experiment name, and string for the plot filename (optional)
s=s_cmn; s.exp_1='cookie.C.p_worjh2.r.SPIN'; s.altfilename='_TEST.p_worjh2.r';
plot_timeseries(s);
% (2) adding an optional SSS time-series
s.exp_1='cookie.C.p_worjh2.r.SPIN'; s.altfilename='_TEST.p_worjh2.r.plusSSS';
% dataid_1  ==  optioal first time-series name -- omitting
%               'timeseries_' and '.txt'
% datacol_1 ==  which data column (remember that column 1 is always time)
% axis1min  ==  minimum y-axis value (of timeseries data #1)
% axis1max  ==  maximum y-axis value (of timeseries data #1)
s=s_cmn; s.dataid_1='climate_SSS'; s.datacol_1=2; s.axis1min=34.0; s.axis1max=36.0;
plot_timeseries(s);
%
% *** Cao et al. ******************************************************** %
%
% NOTE: the existence of CO2 (and isotopes) in an experiment 
%       is automatically identified and plotted as the 3rd panel down
% (1) basic plot
s=s_cmn; s.exp_1='cookie.CB.p_worjh2.BASES.caoetal.SPIN'; s.altfilename='_TEST.p_worjh2.BASES.1';
plot_timeseries(s);
% (2) climate (and pCO2) + one additional optional timeseries
s=s_cmn; s.exp_1='cookie.CB.p_worjh2.BASES.caoetal.SPIN'; s.altfilename='_TEST.p_worjh2.BASES.2';
% specify AMOC maximum, which is column #5 in experiments using a standard
% modern configuration (worbe2, worjh2)
s.dataid_1='climate_opsi'; s.datacol_1=5; s.axis1min=0.0; s.axis1max=15.0;
plot_timeseries(s);
% (3) climate (and pCO2) + two additional optional timeseries
s=s_cmn; s.exp_1='cookie.CB.p_worjh2.BASES.caoetal.SPIN'; s.altfilename='_TEST.p_worjh2.BASES.3';
s.dataid_1='climate_opsi'; s.datacol_1=5; s.axis1min=0.0; s.axis1max=15.0; s.axis1str='AMOC (Sv)';
% add global POC export (in addition to AMOC), which for mol yr-1 is col #2
% NOTE: (a) convert from mol yr-1 to PgC yr-1 by setting an axis scaling:
%           axis2scl=1.2e-14
%       (b) specify a sensiblel axis lable rather than the default:
%           axis2str='POC export (PgC yr^{-1})'
s.dataid_2='biop_fexport_POC'; s.datacol_2=2; s.axis2scl=1.2e-14; s.axis2min=4.0; s.axis2max=8.0; s.axis2str='POC export (PgC yr^{-1})';
plot_timeseries(s,'');
%
% *********************************************************************** %

% *********************************************************************** %
% *** plot_ocean_2d TEST ************************************************ %
% *********************************************************************** %
%
% *** set common parameters ********************************************* %
% 
% set common parameters (to be assigned each time to structure s) as before
% NOTE: here, the same experiment is used throughout 2D and 3D plotting
%       and so the value of exptstr1 need only be set once (in s_cmn)
s_cmn.exptstr1='cookie.C.p_worjh2.r.SPIN';
s_cmn.timeslice1=9999.5;
s_cmn.plot.format='png';
%
% *** 2D physics ******************************************************** %
%
% (1) seafloor topography (2D netCDF) / no contours
% fieldstr1     ==  2D netCDF variable name to plot
%                   (setting fieldstr2 will create anomolies)
% data.scale    ==  any scaling of the data 
%                   (default is 1.0 so it need not actually be specified)
% data.min      ==  minimum value on the color scale
% data.max      ==  maximum value on the color scale
% data.intervals == how many color invervals between min and max 
s=s_cmn; s.fieldstr1='2Dgrid_topo'; s.data.scale=1.0; s.data.min=0.0; s.data.max=6000.0; s.data.intervals=24;
plot_ocean_2d(s);
% (2) annual average seaice (2D netCDF) / no contours
s=s_cmn; s.fieldstr1='climate_seaice'; s.data.scale=1.0; s.data.min=0.0; s.data.max=100.0; s.data.intervals=20;
% here, change the colorscale to 'red': plot.colors='red'
% (refer to file 'colorscales.pdf' for options)
% and invert the order of the colors: plot.invertcolors=true
s.plot.colors='red'; s.plot.invertcolors=true;
plot_ocean_2d(s);
% (3) PSI only
% for plotting PSI, enable via: psi.show=true and then set the
% min and max values together with intervals of major and minor contours
% NOTE: this is different from ploting normal variables which use the
%       data.min, data.max, data.intervals parameter settings
s=s_cmn; s.psi.show=true; s.psi.min=-50.0; s.psi.max=50.0; s.psi.major=10; s.psi.minor=5;
plot_ocean_2d(s);
% (4) PSI + SST
% you can overlay PSI on another 2D variable -- to do this,
% (a) specify a variable, e.g., fieldstr1='climate_SST'
% (b) specify data.min, data.max, data.intervals parameter settings
%     (and any data scaling factor)
s=s_cmn; s.fieldstr1='climate_SST'; s.data.min=0.0; s.data.max=25.0; s.data.intervals=25; s.psi.show=true; s.psi.min=-50.0; s.psi.max=50.0; s.psi.major=10; s.psi.minor=5;
plot_ocean_2d(s);
% (5) SST + CONTOURS
% a plot of SST, but now with contours overlaying the block colors
s=s_cmn; s.fieldstr1='climate_SST'; s.data.min=0.0; s.data.max=25.0; s.data.intervals=25;
% turn on overlying contours with: contours.show=true
% and then specify the contour intervals (contours.major, contours.minor)
s.contours.show=true; s.contours.major=5.0; s.contours.minor=1.0;
plot_ocean_2d(s);
%
% *** 3D physics -- horizontal ****************************************** %
%
% (6) ocean surface temperature
% another way to plot SST is temperature distribution at the surface
% (k=16 in a 16-level ocean configuration), i.e., as per 2D plotting
% except now specifying the ocean layer: slice=16 
s=s_cmn; s.fieldstr1='climate_temp'; s.slice=16; s.data.min=0.0; s.data.max=25.0; s.data.intervals=25;
plot_ocean_3d_horz(s);
% (7) ocean surface temperature + CONTOURS
s=s_cmn; s.fieldstr1='climate_temp'; s.slice=16; s.data.min=0.0; s.data.max=25.0; s.data.intervals=25;
% adding overlying contours as per 2D example
s.contours.show=true; s.contours.major=5.0; s.contours.minor=1.0;
plot_ocean_3d_horz(s);
% (8) ocean surface salinity + CONTOURS
s=s_cmn; s.fieldstr1='climate_sal'; s.slice=16; s.data.min=32.0; s.data.max=36.0; s.data.intervals=16;
s.contours.show=true; s.contours.major=0.5; s.contours.minor=0.125;
plot_ocean_3d_horz(s);
% (9) ocean sea-floor temperature
% you can plot more than just horizontal slices -- setting s.slice=-1
% will give you the benthic (immediately above seafloor) surface
% all ohter parameter setting are as per for ocena surface temperature
s=s_cmn; s.fieldstr1='climate_temp'; s.slice=-1; s.data.min=0.0; s.data.max=6.0; s.data.intervals=30;
plot_ocean_3d_horz(s);
% (10) annual average ocean surface currents (vectors) + speed (colors)
% the parameter s.uv.show=true will add curent vectors on top of a plot
% if no variable field is specified, speed will be plotted in color
% (with the vector field on top)
s=s_cmn; s.uv.show=true; s.uv.scale=1.5; s.slice=16; s.data.min=0.0; s.data.max=1.0; s.data.intervals=20;
plot_ocean_3d_horz(s);
% (11) annual average mid-depth currents (vectors) + SST (colors)
% now to try a deeper (than surface) slice (s.slice=6)
% NOTE: the vector length scale is changed: uv.scale=0.5
%       to reflect the slower velocities at depth: 
s=s_cmn; s.uv.show=true; s.uv.scale=1.5; s.slice=6; 
s.fieldstr1='climate_temp'; s.data.min=0.0; s.data.max=6.0; s.data.intervals=30;
plot_ocean_3d_horz(s);
%
% *** 3D physics -- vertical ******************************************** %
%
% <<< WORK IN PROGRESS >>>
% % OPSI (global zonal mean / meridional overturning circulation)
% plot_fields_biogem_3d_i(DUM_EXP,'','ocn_temp','',DUM_T,-1,0,'',1.0,0.0,30.0,30,'','plot_fields_SETTINGS_OPSI',[DUM_NAME '.PHYS_OPSI']);
%
% *********************************************************************** %


% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
close all;
%
% *********************************************************************** %
