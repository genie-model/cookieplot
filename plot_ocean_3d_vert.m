function [OUTPUT] = plot_ocean_3d_vert(s_in,str_in)
arguments
    s_in
    str_in = '';
end
% plot_ocean_3d_vert
%
%   *******************************************************************   %
%   *** biogem 3-D (LON-LAT) DATA PLOTTING ****************************   %
%   *******************************************************************   %
%
%   plot_ocean_3d_vert(s_in,str_in)
%   plots the BIOGEM 3D netCDF data file 'fields_ocean_3d.nc'
%
%   *******************************************************************   %

% *********************************************************************** %
% ***** HISTORY ********************************************************* %
% *********************************************************************** %
%
%   25/10/31    FORKED FROM MUFFINPLOT
%   25/10/31    changed aassumed biogem directory -> results
%               updated netCDF filename
%               *** VERSION 2.00 ******************************************
%   26/09/23    updated for new cookie netCDF variable names
%               REVISED PARAMETER PASSING -- NOW STRUCTURE
%               (for now, the options file is left alone, 
%               as are the original variable names)
%               [move to new format needs to be completed]
%               add option to save main plot as image (bitmap) 
%               instead of (only) PDF format
%               (to improve speed ... hopefully ...)
%               *** VERSION 2.01 ******************************************
%
% *********************************************************************** %
%%

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
% 
% set version!
par_ver = 2.01;
% set function name
str_function = mfilename;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% set current path
str_current_path = pwd;
% path to cookieplot directory from the location of this function
% NOTE: assuming that this function still resides in that directory ...
%       and that cookieplot is in the MATLAB search path
str_cookieplot_path = fileparts(which(str_function));
% add SETTINGS path
par_pathin = [str_cookieplot_path '/' 'SETTINGS'];
if ~(exist(par_pathin,'dir') == 7)
    disp([' * ERROR: Cannot find SETTINGS directory of cookieplot: ' par_pathin]); disp([' ']); return;
else
    addpath(par_pathin);
end
% determine which stupid version of MUTLAB we are using
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);
% initialize internal structure
s = [];
%
% *** ingest parameters settings **************************************** %
% 
% load default plotting options
if ~isfield(s_in,'popt'), s_in.popt='settings_ocean'; end; popt = s_in.popt;
eval(popt); 
% identify missing s_in structure elements and create defaults
% set old format parameters equal to s_in element values
if ~isfield(s_in,'exptstr1'), s_in.exptstr1 = ''; end;  exp_1 = s_in.exptstr1;
if ~isfield(s_in,'exptstr2'), s_in.exptstr2 = ''; end;  exp_2 = s_in.exptstr2;
if ~isfield(s_in,'timeslice1'), s_in.timeslice1 = -1; end;  timesliceid_1 = s_in.timeslice1;
if ~isfield(s_in,'timeslice2'), s_in.timeslice2 = -1; end;  timesliceid_2 = s_in.timeslice2;
if ~isfield(s_in,'fieldstr1'), s_in.fieldstr1 = ''; end;  dataid_1 = s_in.fieldstr1;
if ~isfield(s_in,'fieldstr2'), s_in.fieldstr2 = ''; end;  dataid_2 = s_in.fieldstr2;
if ~isfield(s_in,'slice'), s_in.slice = 99;  end
if ~isfield(s_in,'maskname'), s_in.maskname = ''; end
if ~isfield(s_in,'datastr'), s_in.datastr = ''; end;  overlaydataid = s_in.datastr;
if ~isfield(s_in,'savestr'), s_in.savestr = ''; end; altfilename = s_in.savestr;
%
% *** parse structure-only settings ************************************* %
%
% main plots
if isfield(s_in,'plot')
    if ~isfield(s_in.plot,'showmain'),     s_in.plot.showmain = true; end
    if ~isfield(s_in.plot,'showextra'),    s_in.plot.showextra = false; end
    if ~isfield(s_in.plot,'format'),       s_in.plot.format = 'png'; end
    if ~isfield(s_in.plot,'colors'),       s_in.plot.colors = 'plasma'; end
    if ~isfield(s_in.plot,'invertcolors'), s_in.plot.invertcolors = false; end
    if ~isfield(s_in.plot,'title'),        s_in.plot.title = ''; end
else
    s_in.plot.showmain = true;
    s_in.plot.showextra = false;
    s_in.plot.format = 'png';
    s_in.plot.colors = 'plasma';
    s_in.plot.invertcolors = false;
    s_in.plot.title = '';
end
% grid settings
if isfield(s_in,'grid')
    if ~isfield(s_in.grid,'origin'),   s_in.grid.origin = -180.0; end
    if ~isfield(s_in.grid,'equallat'), s_in.grid.equallat = false; end
else
    s_in.grid.origin = -180.0;
    s_in.grid.equallat = false;
end
% data (model)
if isfield(s_in,'data')
    if ~isfield(s_in.data,'scale'),     s_in.data.scale = 1.0; end
    if ~isfield(s_in.data,'min'),       s_in.data.min = 0.0; end
    if ~isfield(s_in.data,'max'),       s_in.data.max = 0.0; end
    if ~isfield(s_in.data,'intervals'), s_in.data.intervals = 10; end
    if ~isfield(s_in.data,'offset'),    s_in.data.offset = 0.0; end
    if ~isfield(s_in.data,'log10'),     s_in.data.log10 = false; end
    if ~isfield(s_in.data,'save'),      s_in.data.save = false; end
else
    s_in.data.scale = 1.0;
    s_in.data.min = 0.0;
    s_in.data.max = 0.0;
    s_in.data.intervals = 10;
    s_in.data.offset = 0.0;
    s_in.data.log10 = false;
    s_in.data.save = false;
end
% contours
if isfield(s_in,'contours')
    if ~isfield(s_in.contours,'show'),    s_in.contours.show = false; end
    if ~isfield(s_in.contours,'major'),   s_in.contours.major = 1; end
    if ~isfield(s_in.contours,'minor'),   s_in.contours.minor = 1; end
    if ~isfield(s_in.contours,'showlabels'), s_in.contours.showlabels = true; end
    if ~isfield(s_in.contours,'dashneg'), s_in.contours.dashneg = true; end
else
    s_in.contours.show = false;
    s_in.contours.major = 1;
    s_in.contours.minor = 1;
    s_in.contours.showlabels = true;
    s_in.contours.dashneg = true;
end
% overlays
if isfield(s_in,'opsi')
    if ~isfield(s_in.opsi,'show'),  s_in.opsi.show = false; end
    if ~isfield(s_in.opsi,'min'),   s_in.opsi.min = -50.0; end
    if ~isfield(s_in.opsi,'max'),   s_in.opsi.max = 50.0; end
    if ~isfield(s_in.opsi,'major'), s_in.opsi.major = 10.0; end
    if ~isfield(s_in.opsi,'minor'), s_in.opsi.minor = 5.0; end
    if ~isfield(s_in.opsi,'basin'), s_in.opsi.basin = 'g'; end
else
    s_in.opsi.show = false;
    s_in.opsi.min = -50.0;
    s_in.opsi.max = 50.0;
    s_in.opsi.major = 10.0;
    s_in.opsi.minor = 5.0;
    s_in.opsi.basin = 'g';
end
if isfield(s_in,'uv')
    if ~isfield(s_in.uv,'show'),  s_in.uv.show = false; end
    if ~isfield(s_in.uv,'scale'), s_in.uv.scale = 1.0; end
else
    s_in.uv.show = false;
    s_in.uv.scale = 1.0;
end
if isfield(s_in,'overlay')
    if ~isfield(s_in.overlay,'only'), s_in.overlay.only = false; end
else
    s_in.overlay.only = false;
end


%
% *** filter options **************************************************** %
%
% normalize and validate structure-only plotting options
s_in.plot.format = lower(char(s_in.plot.format));
% OPSI-only plotting
% NOTE: create dummy plotting limits for OPSI
if (s_in.opsi.show && isempty(s_in.fieldstr1))
    s_in.fieldstr1='climate_rho';
    dataid_1 = s_in.fieldstr1;
    s_in.slice = 0;
    s_in.overlay.only = true;
    s_in.data.min = s_in.opsi.min;
    s_in.data.max = s_in.opsi.max;
    s_in.data.intervals = floor((s_in.opsi.max - s_in.opsi.min)/s_in.opsi.minor);
end
% set slice == 0 if an averaging mask is applied
if (~isempty(s_in.maskname))
    s_in.slice = 0;
end
%
% *** backwards compatability ******************************************* %
% 
% data point scaling
if ~exist('data_scalepoints','var'), data_scalepoints = 'n'; end
% data saving
% if ~exist('data_save','var'),        data_save = 'y'; end % save (mode) data?
if ~exist('data_saveall','var'),     data_saveall = 'n'; end
if ~exist('data_saveallinfo','var'), data_saveallinfo = 'n'; end
if ~exist('data_output_old','var'),  data_output_old = 'y'; end % return STATM
% extracting min / max / range from seasonal data
if ~exist('data_minmax','var'),      data_minmax  = ''; end
if ~exist('data_nseas','var'),       data_nseas   = 0; end
% model-data
if ~exist('data_seafloor','var'),    data_seafloor = 'n'; end
if ~exist('data_ijk_near','var'),    data_ijk_near = 'n'; end
% plotting
if ~exist('contour_hlt2','var'),     contour_hlt2 = contour_hlt; end
if ~exist('contour_hltval2','var'),  contour_hltval2 = contour_hltval; end
if ~exist('plot_psi','var'),         plot_psi = 'n'; end
% paths
if ~exist('par_pathin','var'),   par_pathin   = 'cgenie_output'; end
if ~exist('par_pathlib','var'),  par_pathlib  = 'source'; end
if ~exist('par_pathout','var'),  par_pathout  = 'PLOTS'; end
if ~exist('par_pathdata','var'), par_pathdata = 'DATA'; end
if ~exist('par_pathmask','var'), par_pathmask = 'MASKS'; end
if ~exist('par_pathexam','var'), par_pathexam = 'EXAMPLES'; end
% plotting panel options
if ~exist('plot_profile','var'), plot_profile = 'y'; end % PLOT PROFILE
if ~exist('plot_zonal','var'),   plot_zonal   = 'y'; end % PLOT ZONAL
if ~exist('plot_hist','var'),    plot_hist    = 'n'; end % PLOT HISTOGRAM
if ~exist('plot_histc_SETTINGS','var'), plot_histc_SETTINGS = 'plot_histc_SETTINGS'; end % histc plotting settings
%
% *** initialize parameters ********************************************* %
% 
% set axes
lat_min = -090;
lat_max = +090;
D_min   = 0000;
D_max   = 5000;
% null data value
par_data_null = 9.9E19;
%
par_rEarth = 6371000.0;
%
% *** DEFINE COLORS ***************************************************** %
%
% define contonental color
color_g = [0.75 0.75 0.75];
% define no-data color
color_b = [0.00 0.00 0.00];
%
% *** SCALING *********************************************************** %
% 
% set default data scaling
if s_in.data.scale == 0.0
    s_in.data.scale = 1.0;
    s_in.data.min   = 0.0;
    s_in.data.max   = 0.0;
    s_in.data.intervals = 10;
end
if strcmp(data_scalepoints,'n')
    datapoint_scale = 1.0;
else
    datapoint_scale = s_in.data.scale;
end
% set global flag if no alt plotting scale is set
% NOTE: catch possibility of one axis being set, but the other @ default
%       (min and max with indetical values)
if ((plot_lat_min == plot_lat_max) && (plot_D_min == plot_D_max))
    plot_global = true;
    plot_xy_scaling = 1.0;
    plot_lat_min = lat_min;
    plot_lat_max = lat_max;
    plot_D_min = D_min;
    plot_D_max = D_max;
else
    plot_global = false;
    if (plot_lat_min == plot_lat_max)
        plot_lat_min = lat_min;
        plot_lat_max = lat_max;
    end
    if (plot_D_min == plot_D_max)
        plot_D_min = D_min;
        plot_D_max = D_max;
    end
    plot_xy_scaling = ((plot_D_max - plot_D_min)/(D_max - D_min)) / ((plot_lat_max - plot_lat_min)/(lat_max - lat_min));
end
%
% *** OPTIONAL PLOTTING SCALE ******************************************* %
%
if ~isempty(contour_file)
    % load data
    contour_data = load(contour_file,'-ascii');
    % adjust if necessary so that contour_data(1) is the lowest value
    if (contour_data(1) > contour_data(end)), contour_data = flipud(contour_data); end
    % adjust number of contours
    % NOTE: remember that s_in.data.intervals is the number of intervals, not the number of contours (which  is s_in.data.intervals+1)
    s_in.data.intervals = length(contour_data) - 1;
    % set max,min limits
    s_in.data.min = contour_data(1);
    s_in.data.max = contour_data(end);
end
%
% *** SET PATHS & DIRECTORIES ******************************************* %
% 
% find current path
str_current_path = pwd;
% find where str_function lives ...
% remove its name (+ '.m' extension) from the returned path ...
str_function_path = which(str_function);
str_function_path = str_function_path(1:end-length(str_function)-3);
% check source code directory and add search path
if ~(exist([str_function_path '/' par_pathlib],'dir') == 7)
    disp([' * ERROR: Cannot find source directory']);
    disp([' ']);
    return;
else
    addpath([str_function_path '/' par_pathlib]);
end
% check masks directory and add search path
if (exist([str_current_path '/' par_pathmask],'dir') == 7)
    addpath([str_current_path '/' par_pathmask]);
elseif (exist([str_function_path '/' par_pathmask],'dir') == 7)
    addpath([str_function_path '/' par_pathmask]);
else
    disp([' * ERROR: Cannot find MASKS directory -- was it moved ... ?']);
    disp([' ']);
    return;
end
% set input path
par_pathin = [str_current_path '/' par_pathin];
if ~(exist(par_pathin,'dir') == 7)
    disp([' * ERROR: Cannot find experiment results directory']);
    disp([' ']);
    return;
end
% set output path
par_pathout = [str_current_path '/' par_pathout];
if ~(exist(par_pathout,'dir') == 7), mkdir(par_pathout);  end
% check/add data path
if ~(exist([str_current_path '/' par_pathdata],'dir') == 7)
    mkdir([str_current_path '/' par_pathdata]); 
end
addpath([str_current_path '/' par_pathdata]);
% now make make str_function text-friendly
str_function = strrep(str_function,'_','-');
%
% *** FILTER OPTIONS **************************************************** %
% 
% there will be no site labels if data is averaged ...
if (data_ijk_mean == 'y'), data_sitelabel = 'n'; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** INITIALIZE ARRAYS ************************************************* %
% *********************************************************************** %
%
xm = [];
ym = [];
zm = [];
lonm = [];
lone = [];
lonw = [];
latm = [];
latn = [];
lats = [];
laym = [];
layt = [];
layb = [];
rawdata=[];
data_0=[];
data_1=[];
data_2=[];
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA & LOAD (OPTIONAL) MASK FILE ********************** %
% *********************************************************************** %
%
% open netCDF file -- test for 'experiment format' or not
% NOTE: other format is indicated by '.nc' extension as experiment ID
if strcmp(exp_1(end-2:end),'.nc')
    ncid_1=netcdf.open(exp_1,'nowrite');
    loc_flag_unpack = false;
else
    % test for experiment
    data_dir = [par_pathin '/' exp_1];
    if (exist(data_dir, 'dir') == 0)
        disp(['ERROR: Experiment cannot be found.']);
        if (exist(par_pathin, 'dir') == 0)
            disp(['INFO: Path: ' par_pathin ' cannot be found.']);
        else
            disp(['INFO: Path: ' par_pathin ' exists.']);
            disp(['INFO: Experiment name: ' exp_1 ' cannot be found.']);
        end
        if (exist([data_dir '.tar.gz'],'file'))
            disp(['INFO: Archive: ' [data_dir '.tar.gz'] ' exists.']);
            disp([' * UN-PACKING ...']);
            untar([data_dir '.tar.gz'],par_pathin);
            loc_flag_unpack = true;
        else
            return;
        end
    else
        loc_flag_unpack = false;
    end
    ncid_1=netcdf.open([par_pathin '/' exp_1 '/results/fields_ocean_3d.nc'],'nowrite');
end
% read netCDf information
[~,nvars,~,~] = netcdf.inq(ncid_1); % [ndims,nvars,ngatts,unlimdimid]
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid data
varid  = netcdf.inqVarID(ncid_1,'2Dgrid_level');
grid_k1 = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_k1 = grid_k1';
% calculate grid dimensions
varid  = netcdf.inqDimID(ncid_1,'lat');
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen] 
jmax = dimlen;
varid  = netcdf.inqDimID(ncid_1,'lon');
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen] 
imax = dimlen;
varid  = netcdf.inqDimID(ncid_1,'zt');
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen] 
kmax = dimlen;
% load remaining grid information
varid  = netcdf.inqVarID(ncid_1,'zt');
grid_zt = netcdf.getVar(ncid_1,varid);
grid_zt = flipud(grid_zt);
varid  = netcdf.inqVarID(ncid_1,'zt_edges');
grid_zt_edges = netcdf.getVar(ncid_1,varid);
grid_zt_edges = flipud(grid_zt_edges);
% determine grid limits
zt_min = min(grid_zt_edges);
zt_max = max(grid_zt_edges);
% set data limts in k space (based on plotting depth limits)
loc_k = find(grid_zt > plot_D_min);
data_kmax = max(loc_k);
loc_k = find(grid_zt < plot_D_max);
data_kmin = min(loc_k);
% calculate topography
for i = 1:imax,
    for j = 1:jmax,
        if grid_k1(j,i) <= kmax
            topo(j,i) = -grid_zt_edges(grid_k1(j,i));
        else
            topo(j,i) = 0.0;
        end
    end
end
% load and calculate remaining grid information
varid  = netcdf.inqVarID(ncid_1,'lat');
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
grid_lon = netcdf.getVar(ncid_1,varid);
[latm laym] = meshgrid(grid_lat,-grid_zt);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid);
[lats layb] = meshgrid(grid_lat_edges(1:jmax),-grid_zt_edges(1:kmax));
[latn layt] = meshgrid(grid_lat_edges(2:jmax+1),-grid_zt_edges(2:kmax+1));
% create area grid
% NOTE: based on:
%       gi_area(:,:) = 2.0*pi*(par_rEarth^2)*(sin((pi/180.0)*gi_latn) - sin((pi/180.0)*gi_lats)).*((gi_lone-gi_lonw)/360.0);
grid_area = NaN(jmax,imax);
for i = 1:imax,
    for j = 1:jmax,
        if grid_k1(j,i) > kmax
            grid_area(j,i) = 0.0;
        else
            grid_area(j,i) = 2.0*pi*(par_rEarth^2)*(sin((pi/180.0)*grid_lat_edges(j+1)) - sin((pi/180.0)*grid_lat_edges(j))).*((grid_lon_edges(i+1)-grid_lon_edges(i))/360.0);
        end
    end
end
% create cell volume array; also depth
% NOTE: assume equal area grid, normaalized area
% NOTE: multiple by 1.0 becasue ... (?) (to ensure correct format?)
data_V = zeros(kmax,jmax,imax);
data_D = zeros(kmax,jmax,imax);
for k = 1:kmax,
    data_V(k,:,:) = grid_area(:,:)*(grid_zt_edges(k) - grid_zt_edges(k+1));
    data_D(k,:,:) = grid_zt(k);
end
%
grid_lon_origin = grid_lon_edges(1);
% test for mask 'type'
% load mask data or create mask
% NOTE: mask single location coordinate is (i,j)
%       but written to the array as (j,i)
% NOTE: when loading from file,
%       flip in j-direction to make consistent with netCDF grid
if (~isempty(s_in.maskname))
    if isnumeric(s_in.maskname)
        mask = zeros(jmax,imax);
        mask(s_in.maskname(2),s_in.maskname(1)) = 1.0;
        s_in.maskname = ['i', num2str(s_in.maskname(1)), 'j', num2str(s_in.maskname(2))];
    elseif ischar(s_in.maskname)
        maskfile = s_in.maskname;
        mask = load([maskfile],'-ascii');
        mask = flipdim(mask,1);
    else
        disp([' ']);
        error('*WARNING*: Unknown mask parameter type (must be character array or vector location) ... ')
    end
end
%
if (~isempty(s_in.maskname))
    topo = mask.*topo;
elseif ((s_in.slice == 0) || (s_in.slice == -1))
    mask = zeros(jmax,imax);
    mask(:,:) = 1.0;
    topo = mask.*topo;
else
    if ((s_in.slice > imax) || (s_in.slice < -1))
        disp([' ']);
        error('*WARNING*: Value of s.slice out-of-range (-1 to imax): ENDING ... ');
    end
    mask = zeros(jmax,imax);
    mask(:,s_in.slice) = 1.0;
    topo = mask.*topo;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET PRIMARY GRIDDED DATASET *************************************** %
% *********************************************************************** %
%
% check that the year exists
varid  = netcdf.inqVarID(ncid_1,'time');
timeslices = netcdf.getVar(ncid_1,varid);
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen]
clear time;
while exist('time','var') == 0
    for n = 1:dimlen,
        if double(int32(100*timeslices(n)))/100 == timesliceid_1
            time = timesliceid_1;
            tid_1 = n;
        end
    end
    if exist('time','var') == 0
        disp('   > WARNING: Year #1 must be one of the following;');
        format long g;
        double(int32(100*timeslices(:)))/100
        format;
        timesliceid_1 = input('   > Time-slice year: ');
    end
end
% check that the variable name exists
varid_1 = [];
while isempty(varid_1)
    for n = 0:nvars-1,
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
        if strcmp(varname,dataid_1)
            varid_1 = n;
        end
    end
    if isempty(varid_1)
        disp('   > WARNING: Variable #1 name must be one of the following;');
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
            varname
        end
        dataid_1 = input('   > Variable name: ','s');
    end
end
% load data
% flip array around diagonal to give (j,i) array orientation
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid_1);
netcdfdata = netcdf.getVar(ncid_1,varid_1);
if length(dimids) == 4
    rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax,tid_1);
    switch data_minmax,
        case {'min'}
            for j = 1:jmax,
                for i = 1:imax,
                    for k = 1:kmax,
                        rawdata(i,j,k) = min(netcdfdata(i,j,k,tid:tid+data_nseas-1));
                    end
                end
            end
        case {'max','minmax'}
            for j = 1:jmax,
                for i = 1:imax,
                    for k = 1:kmax,
                        rawdata(i,j,k) = max(netcdfdata(i,j,k,tid:tid+data_nseas-1));
                    end
                end
            end
    end
    for n = 1:kmax,
        data_1(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
    end
elseif length(dimids) == 3
    rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax);
    for n = 1:kmax,
        data_1(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
    end
elseif length(dimids) == 2
    rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax);
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
else
    data_1(:,:,:) = NaN;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET ALTERNATIVE GRIDDED DATASET *********************************** %
% *********************************************************************** %
%
% *** ALT EXPERIMENT **************************************************** %
%
% open new netCDF file if necessary, else reuse 1st netCDF dataset
if ~isempty(exp_2)
    % open netCDF file -- test for 'experiment format' or not
    % NOTE: other format is indicated by '.nc' extension as experiment ID
    if strcmp(exp_2(end-2:end),'.nc'),
        ncid_2=netcdf.open(exp_2,'nowrite');
    else
        ncid_2=netcdf.open([par_pathin '/' exp_2 '/results/fields_ocean_3d.nc'],'nowrite');
    end
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
    % load data
    varid_2 = netcdf.inqVarID(ncid_2,dataid_2);
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid_2);
else
    ncid_2 = ncid_1;
end
%
% *** ALT TIME-SLICE **************************************************** %
%
if timesliceid_2 >= 0.0
    % check that the year exists
    varid  = netcdf.inqVarID(ncid_2,'time');
    timeslices = netcdf.getVar(ncid_2,varid);
    [dimname, dimlen] = netcdf.inqDim(ncid_2,varid);
    clear time;
    while exist('time','var') == 0
        for n = 1:dimlen,
            if double(int32(100*timeslices(n)))/100 == timesliceid_2
                time = timesliceid_2;
                tid_2 = n;
            end
        end
        if exist('time','var') == 0
            disp('   > WARNING: Year #2 must be one of the following;');
            format long g;
            double(int32(100*timeslices(:)))/100
            format;
            timesliceid_2 = input('   > Time-slice year: ');
        end
    end
else
    tid_2 = tid_1;
end
%
% *** ALT DATA FIELD **************************************************** %
%
if ~isempty(dataid_2)
    % check that the variable name exists
    varid_2 = [];
    while isempty(varid_2)
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
            if strcmp(varname,dataid_2)
                varid_2 = n;
            end
        end
        if isempty(varid_2)
            disp('   > WARNING: Variable #2 name must be one of the following;');
            for n = 0:nvars-1,
                [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
                varname
            end
            dataid_2 = input('   > Variable name: ','s');
        end
    end
else
    varid_2 = varid_1;
end
%
% *** SET DATA ********************************************************** %
%
if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2)),
    % load data
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid_2);
    netcdfdata = netcdf.getVar(ncid_2,varid_2);
    if length(dimids) == 4
        rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax,tid_2);
        switch data_minmax,
            case {'min','minmax'}
                for j = 1:jmax,
                    for i = 1:imax,
                        for k = 1:kmax,
                            rawdata(i,j,k) = min(netcdfdata(i,j,k,tid_2:tid_2+data_nseas-1));
                        end
                    end
                end
            case 'max'
                for j = 1:jmax,
                    for i = 1:imax,
                        for k = 1:kmax,
                            rawdata(i,j,k) = max(netcdfdata(i,j,k,tid_2:tid_2+data_nseas-1));
                        end
                    end
                end
        end
        for n = 1:kmax,
            data_2(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 3
        rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax);
        for n = 1:kmax,
            data_2(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax);
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_2(:,:,:) = NaN;
    end
    data_anomoly = 'y';
else
    data_2(:,:,:) = 0.0;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** STREAMFUNCTION **************************************************** %
% *********************************************************************** %
%
% open streamfunction data (if selected)
%
if (s_in.opsi.show)
        %
    if (plot_opsi_calc == 'n')
        %
        % *** USE PRE-SAVED 2D DATA ********************************************* %
        %
        % open netCDF file
        ncid_0=netcdf.open([par_pathin '/' exp_1 '/results/fields_ocean_2d.nc'],'nowrite');
        % read netCDf information
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_0);
        % load grid information
        varid  = netcdf.inqVarID(ncid_0,'zt');
        [dimname, dimlen] = netcdf.inqDim(ncid_0,varid);
        zmax = dimlen;
        [latm ztm] = meshgrid(grid_lat,-grid_zt);
        varid  = netcdf.inqVarID(ncid_0,'zt_edges');
        opsigrid_zt_edges = flipud(netcdf.getVar(ncid_0,varid));
        varid  = netcdf.inqVarID(ncid_0,'lat_edges');
        opsigrid_lat_edges = netcdf.getVar(ncid_0,varid);
        [opsigrid_lat opsigrid_zt] = meshgrid(opsigrid_lat_edges(1:jmax+1),-opsigrid_zt_edges(1:zmax+1));
        % set variable name
        switch (s_in.opsi.basin)
            case {'g'}
                varid  = netcdf.inqVarID(ncid_0,'climate_opsi');
            case 'a'
                varid  = netcdf.inqVarID(ncid_0,'climate_opsia');
            case 'p'
                varid  = netcdf.inqVarID(ncid_0,'climate_opsip');
            otherwise
                disp('Unknown opsi definition.')
        end
        % open data
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_0,varid);
        rawdata = netcdf.getVar(ncid_0,varid);
        if length(dimids) == 3
            rawdata = rawdata(1:jmax+1,1:zmax+1,tid_1);
            data_0 = rawdata(1:jmax+1,1:zmax+1)';
        elseif length(dimids) == 2
            rawdata = rawdata(1:jmax+1,1:zmax+1);
            data_0 = rawdata(1:jmax+1,1:zmax+1)';
        else
            data_0 = NaN*data_0;
        end
        % invert data to match grid
        data_0=flipud(data_0);
    else
        %
        % *** CALCULATE OVERTURNING ********************************************* %
        %
        % (1a) setup basic grid
        zmax = kmax;
        [latm ztm] = meshgrid(grid_lat,-grid_zt);
        opsigrid_zt_edges = grid_zt_edges;
        opsigrid_lat_edges = grid_lat_edges;
        [opsigrid_lat opsigrid_zt] = meshgrid(opsigrid_lat_edges(1:jmax+1),-opsigrid_zt_edges(1:zmax+1));
        % (1b) create additional grid details
        loc_cv = zeros(jmax);
        %%%[loc_sv, loc_s, loc_dz, loc_dza] = get_grid_genie_full(imax,jmax,kmax);
        [loc_sv, loc_s, loc_dz, loc_dza] = make_genie_grid(imax,jmax,kmax,zt_max,0.0,plot_equallat,0.0);
        for j = 1:jmax,
            loc_cv(j) = sqrt(1 - loc_sv(j)*loc_sv(j));
        end
        % (2) load velocity data
        varid  = netcdf.inqVarID(ncid_1,'climate_v');
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
        rawdata = netcdf.getVar(ncid_1,varid);
        if length(dimids) == 4
            rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax,tid_1);
            for n = 1:kmax,
                data_v(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
            end
        elseif length(dimids) == 3
            rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax);
            for n = 1:kmax,
                data_v(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
            end
        elseif length(dimids) == 2
            rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
            data_v(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
        else
            data_v(:,:,:) = NaN;
        end
        % scale velocity
        % % NOTE: saved velocity field has been scaled by usc (0.05);
        % %       this factor has to be undone for the MOC calculation ...
        % data_v = data_v/0.05;
        % create streamfunction arrays
        loc_ou = zeros(jmax,kmax);
        loc_opsi = zeros(jmax+1,kmax+1);
        % calculate streamfunction
        % NOTE: compared to FORTRAN code [below]:
        %       k loop in reverse order because vertical array is
        %       up-side-down (and k+1 needed rather than k-1)
        %       also: shifted in j (=> differences in 't' and 'c' grids)
        for j = 1:jmax,
            for k = kmax:-1:2,
                for i=1:imax,
                    if ((data_v(k,j,i) < 0.999E19) && (mask(j,i) ~= 0)),
                        loc_ou(j,k) = loc_ou(j,k) + loc_cv(j)*data_v(k,j,i)*(2.0*pi/imax);
                    end
                end
                loc_opsi(j+1,k) = loc_opsi(j+1,k+1) + loc_dz(k)*loc_ou(j,k); 
            end
            if (sum(mask(j,:)) == 0), loc_opsi(j+1,:) = NaN; end
        end
        % scale streamfunction and re-orientate array
        % goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6
        % 5e3, 0.05, 6.37E+06
        % NOTE: use zt_max rather than a fixed 5.0E3 m scale depth value
        loc_opsi = (zt_max*0.05*6.37E+06*1.0E-6)*loc_opsi;
        data_0=loc_opsi';
    end
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
%
% create an output filestring for data and plot saving
if (~isempty(s_in.maskname))
    if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2))
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '_MINUS_', exp_2, '.', 'y', num2str(timesliceid_2), '.', dataid_2, '.', s_in.maskname];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', s_in.maskname];
    end
else
    if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2))
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '_MINUS_', exp_2, '.', 'y', num2str(timesliceid_2), '.', dataid_2, '.', 'i', num2str(s_in.slice)];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', 'i', num2str(s_in.slice)];
    end
end
if ~isempty(overlaydataid)
    filename = [filename, '_VS_', overlaydataid];
    if (data_anomoly == 'y')
        filename = [filename, '.ANOM'];
    end
    if (s_in.overlay.only)
        filename = [filename, '.DO'];
    end
end
if (~isempty(altfilename)), filename = altfilename; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FILTER & PROCESS RAW DATA ***************************************** %
% *********************************************************************** %
%
% *** INITALIZE ********************************************************* %
%
xm = latm;
ym = laym;
if ( (data_anomoly == 'y') || (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2)) ),
    data = data_1 - data_2;
else
    data = data_1;
    data_2 = data_1;
end
data = data + s_in.data.offset;
if (s_in.opsi.show)
    opsidata = data_0;
end
% define initial array sizes
zm = zeros(kmax,jmax);
zm_count = zm;
zl = zeros(kmax,1);
zl_V = zl;
zz = zeros(kmax,jmax);
zz_V = zz;
%
z = 0.0;
z_V = 0.0;
n = 0;
%
% *** PROCESS MAIN DATASET ********************************************** %
%
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            % create special case of just benthic and surface data
            if (s_in.slice == -1)
                if ((k ~= grid_k1(j,i)) && (k ~=kmax))
                    data(k,j,i) = NaN;
                end
            end
            % (normal)
            if (topo(j,i) < -grid_zt(k))
                if ((data(k,j,i) > -0.999E19) && (data(k,j,i) < 0.999E19) &&  ~isnan(data(k,j,i)))
                    zm(k,j) = zm(k,j) + data(k,j,i);
                    zm_count(k,j) = zm_count(k,j) + 1.0;
                    zl(k) = zl(k) + data_V(k,j,i)*data(k,j,i);
                    zl_V(k) = zl_V(k) + data_V(k,j,i);
                    zz(k,j) = zz(k,j) + data_V(k,j,i)*data(k,j,i);
                    zz_V(k,j) = zz_V(k,j) + data_V(k,j,i);
                    z = z + data_V(k,j,i)*data(k,j,i); 
                    z_V = z_V + data_V(k,j,i);
                else
                    data(k,j,i) = NaN;
                    data_1(k,j,i) = NaN;
                    data_2(k,j,i) = NaN;
                    data_D(k,j,i) = NaN;
                end
            else
                data(k,j,i) = NaN;
                data_1(k,j,i) = NaN;
                data_2(k,j,i) = NaN;
                data_D(k,j,i) = NaN;
            end
        end
        if (zm_count(k,j) > 0.0)
            if (s_in.data.log10)
                if (zm(k,j) > 0.0)
                    zm(k,j) = log10(zm(k,j)/s_in.data.scale/zm_count(k,j));
                else
                    zm(k,j) = NaN;
                end
            else
                zm(k,j) = zm(k,j)/s_in.data.scale/zm_count(k,j);
                if contour_noneg == 'y'
                    if (zm(k,j) < 0.0)
                        zm(k,j) = 0.0;
                    end
                end
            end
        else
            zm(k,j) = NaN;
        end
        if ~isnan(zm(k,j)), n = n + 1; end
        if (zz_V(k,j) > 0.0)
            zz(k,j) = zz(k,j)/s_in.data.scale/zz_V(k,j);
        else
            zz(k,j) = NaN;
        end
    end
    if (zl_V(k) > 0.0)
        if (s_in.data.log10)
            if (zl(k) > 0.0)
                zl(k) = log10(zl(k)/s_in.data.scale/zl_V(k));
            else
                zl(k) = NaN;
            end
        else
            zl(k) = zl(k)/s_in.data.scale/zl_V(k);
        end
    else
        zl(k) = NaN;
    end
end % end k loop
nmax = n;
% copy zm before it gets transformed ...
overlaydata_zm(:,:) = zm(:,:);
% set topography uniform in i-direction and equal to deepest point found
for j = 1:jmax
    for i = 2:imax
        if (topo(j,i) < topo(j,i-1))
            topo(j,1:i-1) = topo(j,i);
        else
            topo(j,i) = topo(j,i-1);
        end
    end
end
% save data!!!
if (data_saveall == 'y')
    % data 1
    fid = fopen([par_pathout '/' filename '_DATA1POINTS.', dataid_1, '.dat'], 'wt');
    fprintf(fid, '%% Model data field #1 values at mask locations');
    fprintf(fid, '\n');
    fprintf(fid, '%% Format: grid lon, grid lat, grid depth, model data field #1 value');
    fprintf(fid, '\n');
    % loop through all data points
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                if (~isnan(data_1(k,j,i)))
                    fprintf(fid, '%8.3f %8.3f %8.3f %8.6e %s \n', grid_lon(i), grid_lat(j), grid_zt(k), data_1(k,j,i), '%');
                end
            end
        end
    end
    fclose(fid);
    % data 2
    if (~isempty(dataid_2))
    fid = fopen([par_pathout '/' filename '_DATA2POINTS.', dataid_2, '.dat'], 'wt');
        fprintf(fid, '%% Model data field !2 values at mask locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: grid lon, grid lat, grid depth, model data field #2 value');
        fprintf(fid, '\n');
        % loop through all data points
        for k = 1:kmax
            for j = 1:jmax
                for i = 1:imax
                    if (~isnan(data_2(k,j,i)))
                        fprintf(fid, '%8.3f %8.3f %8.3f %8.6e %s \n', grid_lon(i), grid_lat(j), grid_zt(k), data_2(k,j,i), '%');
                    end
                end
            end
        end
        fclose(fid);
    end
end
%
% *** PROCESS OVERTURNING STREAMFUNCTION (IF SELECTED) ****************** %
%
if (s_in.opsi.show)
    opsizm = zeros(zmax+1,jmax+1);
    for j = 1:jmax+1
        for z = 1:zmax+1
            if (opsidata(z,j) < -1.0E6) || (opsidata(z,j) > 1.0E30)
                opsizm(z,j) = NaN;
            else
                opsizm(z,j) = opsidata(z,j);
            end
            %if isnan (zm(z,j)), opsizm(z,j) = NaN; end
        end
    end
    for z = 1:zmax
        if isnan (zm(z,1)), opsizm(z,1) = NaN; end
        if isnan (zm(z,jmax)), opsizm(z,jmax+1) = NaN; end
        for j = 2:jmax+1
            if isnan (zm(z,j-1)), opsizm(z,j) = NaN; end
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD (OPTIONAL) OVERLAY DATA ************************************** %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % add warning ...
    if (data_ijk == 'y')
        disp([' ']);
        error('*WARNING*: Currently there is no facility for loading data in GENIE (i,j,k) format (it must be explicit in lon/lat/depth): ENDING ... ');
    end
    % read data(!)
    [file_data] = fun_read_file(overlaydataid);
    % extract data cell array, and vector format
    cdata    = file_data.cdata;
    v_format = file_data.vformat;
    % determine number of rows and columns
    n_size    = size(cdata);
    n_rows    = n_size(1);
    n_columns = length(v_format);
    % parse call array data
    flag_format = false;
    switch n_columns
        case 4
            if (sum(v_format(1:4)) == 0)
                % lon, lat, depth, value
                overlaydata_raw  = cell2mat(cdata(:,1:4));
                % create dummay (blank space) labels
                overlaylabel_raw = (blanks(n_rows))';
                data_shapecol    = 'n';
                % flag for a valid format
                flag_format      = true;
            end
        case 5
            if ((sum(v_format(1:4)) == 0) && (v_format(5) == 1))
                % lon, lat, depth, value, LABEL
                overlaydata_raw  = cell2mat(cdata(:,1:4));
                overlaylabel_raw = char(cdata(:,5));
                data_shapecol    = 'n';
                % flag for a valid format
                flag_format      = true;
            end
        case 8
            if ((sum(v_format(1:4)) == 0) && sum(v_format(5:8)) == 4)
                % lon, lat, depth, value, LABEL, SHAPE, EDGE COLOR, FILL COLOR
                overlaydata_raw   = cell2mat(cdata(:,1:4));
                overlaylabel_raw  = char(cdata(:,5));
                overlaydata_shape = char(cdata(:,6));
                overlaydata_ecol  = char(cdata(:,7));
                overlaydata_fcol  = char(cdata(:,8));
                data_shapecol     = 'y';
                % flag for a valid format
                flag_format       = true;
            end
        otherwise
            % (caught by the default of ~flag_format)
    end
    if (~flag_format)
        disp([' ']);
        disp([' * ERROR: Data format not recognized:']);
        disp(['   Number of data columns found == ' num2str(n_columns)']);
        disp(['   Columns must be: comma/tab/space-separated.']);
        disp([' ']);
        fclose(fid);
        return;
    end
    % filter label
    for n = 1:n_rows,
        overlaylabel_raw(n,:) = strrep(overlaylabel_raw(n,:),'_',' ');
    end
    % determine data size
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % check for incomplete file read
    if (nmax < n_rows),
        disp([' ']);
        disp([' * ERROR: Problem in data format (e.g. string must not contain spaces):']);
        disp(['   Number of data rows read == ' num2str(nmax)]);
        disp(['   Number of lines in file  == ' num2str(n_rows)]);
        disp([' ']);
        return;        
    end
    % check for mixed up lon-lat ... i.e. not (LON, LAT) format
    if ( (min(overlaydata_raw(:,2)) < -90) || (max(overlaydata_raw(:,2)) > 90) ),
        loc_tmp = overlaydata_raw(:,1);
        overlaydata_raw(:,1) = overlaydata_raw(:,2);
        overlaydata_raw(:,2) = loc_tmp;
        disp([' ']);
        disp([' * WARNING: lon-lat is not in (LON, LAT) column order format:']);
        disp(['   Swapping ...'])
        disp([' ']);
    end
    % determine equivalent (i,j,k)
    % NOTE: filter out masked area embodied in array <data>
    %       (this also includes bathymetry)
    % NOTE: for data deeper than the depth of layer 1, k=0 is returned
    %       (this is a recent change ... it used to be k=1)
    % NOTE: format: data(k,j,i)
    % NOTE: format: mask(j,i)
    % NOTE: accept overlaydata_raw(n,4) by default, then filter
    overlaydata_ijk(:,:) = zeros(size(overlaydata_raw(:,:)));
    for n = 1:nmax
        overlaydata_ijk(n,1:2) = calc_find_ij(overlaydata_raw(n,1),overlaydata_raw(n,2),grid_lon_origin,imax,jmax);
        overlaydata_ijk(n,3)   = calc_find_k(overlaydata_raw(n,3),kmax);
        overlaydata_ijk(n,4)   = overlaydata_raw(n,4);
        % exclude:
        % * below seafloor if (data_seafloor == 'n') (otherwise adjust k)
        % * outside masked area
        % * outside of specified depth range
        % NOTE: mark
        if (overlaydata_ijk(n,3) < grid_k1(overlaydata_ijk(n,2),overlaydata_ijk(n,1)))
            if (data_seafloor == 'y')
                overlaydata_ijk(n,3) = grid_k1(overlaydata_ijk(n,2),overlaydata_ijk(n,1));
            else
                overlaydata_raw(n,4) = NaN;
                overlaydata_ijk(n,4) = NaN;
            end
        end
        if (~mask(overlaydata_ijk(n,2),overlaydata_ijk(n,1)))
            overlaydata_raw(n,4) = NaN;
            overlaydata_ijk(n,4) = NaN;
        end
        if ( (overlaydata_ijk(n,3) < data_kmin) || (overlaydata_ijk(n,3) > data_kmax) )
            overlaydata_raw(n,4) = NaN;
            overlaydata_ijk(n,4) = NaN;
        end
    end
    % remove data in land cells
    if (data_land == 'n')
        for n = 1:nmax
            if ( ~isnan(overlaydata_ijk(n,4)) )
                if ( isnan(overlaydata_zm(overlaydata_ijk(n,3),overlaydata_ijk(n,2))) )
                    overlaydata_raw(n,4) = NaN;
                    overlaydata_ijk(n,4) = NaN;
                end
            end
        end
    end
    % remove filtered data
    overlaylabel_raw(isnan(overlaydata_raw(:,4)),:) = [];
    overlaydata_raw(isnan(overlaydata_raw(:,4)),:)  = [];
    overlaydata_ijk(isnan(overlaydata_ijk(:,4)),:)  = [];
    % update value of nmax
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % change: +vs depths to -vs for plotting
    overlaydata_raw(find(overlaydata_raw(:,3)>0.0),3) = -1.0*overlaydata_raw(find(overlaydata_raw(:,3)>0.0),3);
    % BLAH
    overlaylabel(:,:) = overlaylabel_raw(:,:);
    overlaydata(:,:)  = overlaydata_raw(:,:);
    % copy data
    overlaydata_ijk_old(:,:) = overlaydata_ijk(:,:);
    % grid (and average per cell) data if requested
    % NOTE: data vector length is re-calculated and the value of nmax reset
    if (data_ijk_mean == 'y')
        overlaydata_ijk(:,:) = [];
        overlaydata(:,:)     = [];
        overlaylabel(:,:)    = [];
        overlaylabel         = 'n/a';
        m=0;
        for k = data_kmin:data_kmax
            for j = 1:jmax
                if (~isnan(overlaydata_zm(k,j)))
                    samecell_locations = find((int32(overlaydata_ijk_old(:,3))==k)&(int32(overlaydata_ijk_old(:,2))==j));
                    samecell_n = size(samecell_locations);
                    if (samecell_n(1) > 0)
                        m=m+1;
                        samecell_mean = mean(overlaydata_ijk_old(samecell_locations,4));
                        overlaydata_ijk(m,:) = [NaN j k samecell_mean];
                        %%%overlaydata_ijk(m,:) = [overlaydata_ijk_old(m,1) j k samecell_mean];
                        overlaydata(m,1)     = grid_lon_origin + 360.0*(overlaydata_ijk(m,1) - 0.5)/imax;
                        overlaydata(m,2)     = 180.0*asin(2.0*(overlaydata_ijk(m,2) - 0.5)/jmax - 1.0)/pi;
                        overlaydata(m,3)     = double(laym(k,j));
                        overlaydata(m,4)     = samecell_mean;
                        overlaylabel         = [overlaylabel; 'n/a'];
                    end
                end
            end
        end
        overlaylabel(end,:) = [];
        nmax=m;
    end
    % scale overlay data
    overlaydata(:,4) = overlaydata(:,4)/datapoint_scale;
    % calculate profile
    overlaydata_k(:) = NaN(data_kmax-data_kmin+1,1);
    for k = data_kmin:data_kmax
        samecell_locations = find(int32(overlaydata_ijk_old(:,3))==k);
        samecell_n = size(samecell_locations);
        if (samecell_n(1) > 0)
            samecell_mean = mean(overlaydata_ijk_old(samecell_locations,4));
            overlaydata_k(k) = samecell_mean;
        end
    end   
end
% create and print mask
if ( isempty(dataid_2) && ~isempty(overlaydataid))
    loc_mask = zeros(jmax,imax);
    if (data_ijk_mean == 'y')
        loc_mask = mask;
    else
       for n = 1:nmax,
           loc_mask(overlaydata_ijk(n,2),overlaydata_ijk(n,1)) = 1.0;
       end
    end
    if (s_in.data.save), fprint_2D(loc_mask,[par_pathout '/' filename '.DATAMASK.' str_date '.dat'],'%4.1f','%4.1f',false,false); end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** DATA PROCESSING AND STATS ***************************************** %
% *********************************************************************** %
%
% calculate stats needed for Taylor Diagram (and plot it!)
%
% *** 3D (GRIDDED) DATA ************************************************* %
%
if (~isempty(dataid_2))
    % transform data sets in vectors
    % NOTE: data_1 is format (k,j,i)
    if ((s_in.slice > 0) && isempty(s_in.maskname))
        data_vector_1 = reshape(data_1(:,:,s_in.slice),jmax*kmax,1)/s_in.data.scale;
        data_vector_2 = reshape(data_2(:,:,s_in.slice),jmax*kmax,1)/s_in.data.scale;
        data_vector_D = reshape(data_D(:,:,s_in.slice),jmax*kmax,1);   
    else
        data_vector_1 = reshape(data_1(:,:,:),imax*jmax*kmax,1)/s_in.data.scale;
        data_vector_2 = reshape(data_2(:,:,:),imax*jmax*kmax,1)/s_in.data.scale;
        data_vector_D = reshape(data_D(:,:,:),imax*jmax*kmax,1);
    end
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    if isempty(overlaydataid), nmax = length(data_vector_2); end
    % calculate stats
    % NOTE: STATM = allstats(Cr,Cf)
    % 	    STATM(1,:) => Mean
    % 	    STATM(2,:) => Standard Deviation (scaled by N)
    % 	    STATM(3,:) => Centered Root Mean Square Difference (scaled by N)
    % 	    STATM(4,:) => Correlation
    %       STATM(5,:) => N
    if (data_stats=='y')
        STATM = calc_allstats(data_vector_1,data_vector_2);
        if (s_in.plot.showextra)
            % plot Taylor diagram
            taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
            exportgraphics(gcf,[par_pathout '/' filename '.TaylorDiagram.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
            %%%% plot Target diagram
            %%%targetdiag_vargout = plot_target(STATM(7,1:2),STATM(8,1:2),'r',1.0,[],[]);
            %%%print('-depsc2', [filename, '_TargetDiagram.', str_date, '.eps']);
        end
    end
else
    % set default vector_1
    % NOTE: if discrete data is present, this will be replaced
    % transform data sets in vectors
    % NOTE: data_1 is format (k,j,i)
    if ((s_in.slice > 0) && isempty(s_in.maskname))
        data_vector_1 = reshape(data_1(:,:,s_in.slice),jmax*kmax,1)/s_in.data.scale;
        data_vector_D = reshape(data_D(:,:,s_in.slice),jmax*kmax,1);   
    else
        data_vector_1 = reshape(data_1(:,:,:),imax*jmax*kmax,1)/s_in.data.scale;
        data_vector_D = reshape(data_D(:,:,:),imax*jmax*kmax,1);
    end
    % (don't) transpose vectors
%     data_vector_1 = data_vector_1;
%     data_vector_D = data_vector_D;
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
end
%
% *** DISCRETE DATA ***************************************************** %
%
% NOTE: no scale transformation has been appplied
%       to either gridded or % overlay data
% NOTE: valid only for data on a single depth level
% NOTE: data_vector_2 == model, data_vector_1 == data
if (~isempty(overlaydataid))
    % set overlay data vector
    data_vector_1 = [];
    data_vector_1 = overlaydata(:,4);
    % disable stats if necessary
    if ( length(data_vector_1) < 2 )
        data_stats = 'n';
    elseif ( (max(data_vector_1) - min(data_vector_1)) == 0.0 )
        data_stats = 'n';
    end
    % populate the gridded dataset vector with values corresponding to
    % the overlay data locations
    % NOTE: !!! overlaydata_zm is (k,j) !!!
    % NOTE: <overlaydata_zm> (copied form <zm>) already has scaling applied
    % NOTE: overlaydata_zm values are already potentially averaged, 
    %       hence use data (as per original code)
    %       (why the changce ... ???)
    % NOTE: re-orientate data_vector_2 to match data_vector_1
    % NOTE: format: data(k,j,i)
    data_vector_2 = [];
    data_vector_D = [];
    data_vector_k = [];
    for n = 1:nmax
        if (data_ijk_mean == 'y')
            data_vector_2(n) = overlaydata_zm(int32(overlaydata_ijk(n,3)),int32(overlaydata_ijk(n,2)));
            data_vector_D(n) = -grid_zt(int32(overlaydata_ijk(n,3)));
            data_vector_k(n) = int32(overlaydata_ijk(n,3));
        else
            data_vector_2(n) = data(int32(overlaydata_ijk(n,3)),int32(overlaydata_ijk(n,2)),int32(overlaydata_ijk(n,1)))/s_in.data.scale;
            data_vector_D(n) = data_D(int32(overlaydata_ijk(n,3)),int32(overlaydata_ijk(n,2)),int32(overlaydata_ijk(n,1)));
            data_vector_k(n) = int32(overlaydata_ijk(n,3));
        end
    end
    data_vector_2 = data_vector_2';
    data_vector_D = data_vector_D';
    data_vector_k = data_vector_k';
    % filter data
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    % calculate stats
    if ( (data_stats == 'y') && ~s_in.overlay.only )
        STATM = calc_allstats(data_vector_1,data_vector_2);
        if (s_in.plot.showextra)
            % plot Taylor diagram
            taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
            print('-depsc2', [par_pathout '/' filename, '_TaylorDiagram.', str_date, '.eps']);
            %%%% plot Target diagram
            %%%targetdiag_vargout = plot_target(STATM(7,1:2),STATM(8,1:2),'r',1.0,[],[]);
            %%%print('-depsc2', [filename, '_TargetDiagram.', str_date, '.eps']);
        end
    end
end
%
% *** SAVE STATS DATA *************************************************** %
%
% STATM(1,:) = MEAN;
% STATM(2,:) = STD;
% STATM(3,:) = RMSD;
% STATM(4,:) = CORRELATIONS;
% STATM(5,:) = N;
% STATM(6,:) = TOTAL RMSD;
% STATM(7,:) = STANDARD DEVIATION NORMALISED RMSD;
% STATM(8,:) = NORMALISED BIAS;
% STATM(9,:) = R2;
% STATM(10,:) = M;
if (data_stats == 'y' && data_output_old == 'y')
    if (~isempty(dataid_2) || (~isempty(overlaydataid) && ~s_in.overlay.only) )
        fid = fopen([par_pathout '/' filename '_STATS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '\n');
        fprintf(fid, '=== STATS SUMMARY ===');
        fprintf(fid, '\n');
        fprintf(fid, 'Number of data points, N                           : %8.6e \n', STATM(5,2));
        fprintf(fid, '\n');
        fprintf(fid, 'Mean                                               : %8.6e \n', STATM(1,2));
        fprintf(fid, 'R2                                                 : %8.6e \n', STATM(9,2));
        fprintf(fid, '\n');
        fprintf(fid, 'Standard Deviation (scaled by N)                   : %8.6e \n', STATM(2,2));
        fprintf(fid, 'Root Mean Square Difference (scaled by N)          : %8.6e \n', STATM(3,2));
        fprintf(fid, 'TOTAL Root Mean Square Difference                  : %8.6e \n', STATM(6,2));
        fprintf(fid, 'Correlation                                        : %8.6e \n', STATM(4,2));
        fprintf(fid, 'STANDARD DEVIATION NORMALISED RMSD                 : %8.6e \n', STATM(7,2));
        fprintf(fid, 'NORMALISED BIAS                                    : %8.6e \n', STATM(8,2));
        fprintf(fid, 'M-score                                            : %8.6e \n', STATM(10,2));
        fclose(fid);
    end
end
%
% *** SAVE EQUIVALENT MODEL DATA **************************************** %
%
% save model data at the data locations
if (s_in.data.save)
    if (~isempty(overlaydataid) && ~s_in.overlay.only )
        fid = fopen([par_pathout '/' filename '_MODELDATAPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at data locations');
        fprintf(fid, '\n');
        % print header
        if (data_ijk == 'y')
            fprintf(fid, '%% Format: i, j, k // model lon, model lat, model depth, model value // [data location is (i,j,k)], data value, (data label)');
        elseif (data_ijk_mean == 'y')
            fprintf(fid, '%% Format: i, j, k // model lon, model lat, model depth, model value // [re-gridded data location same as model], re-gridded data value');
        else
            fprintf(fid, '%% Format: i, j, k // model lon, model lat, model depth, model value // data lon, data lat, data depth, data value, (data label)');
        end
        fprintf(fid, '\n');
        % loop through all data points
        for n = 1:nmax
            % set local values
            % NOTE: data_vector_1 == overlaydata(:,4);
            % NOTE: overlaydata(n,2) is already corrected for plot_equallat
            % NOTE: overlaydata(:,3) is negative (for plotting purposes)
            loc_i = int16(overlaydata_ijk(n,1));
            loc_j = int16(overlaydata_ijk(n,2));
            loc_k = int16(overlaydata_ijk(n,3));
            loc_lon_model = grid_lon_origin + 360.0*(overlaydata_ijk(n,1) - 0.5)/imax;
            loc_lon_data  = overlaydata(n,1);
            if (plot_equallat == 'y')
                loc_lat_model = 180.0*((overlaydata_ijk(n,2) - 0.5)/jmax - 0.5);
            else
                loc_lat_model = 180.0*asin(2.0*(overlaydata_ijk(n,2) - 0.5)/jmax - 1.0)/pi;
            end
            loc_lat_data    = overlaydata(n,2);
            loc_depth_model = data_vector_D(n);
            loc_depth_data  = -overlaydata(n,3);
            loc_value_model = data_vector_2(n);
            loc_value_data  = overlaydata(n,4);
            loc_label       = overlaylabel(n,:);
            % print data
            if (data_ijk == 'y')
                fprintf(fid, '%d %d %d %8.3f %8.3f %8.3f %8.6e %8.6e %s \n', loc_i, loc_j, loc_k, loc_lon_model, loc_lat_model, loc_depth_model, loc_value_model, loc_value_data, loc_label);
            elseif (data_ijk_mean == 'y')
                fprintf(fid, '%d %d %d %8.3f %8.3f %8.3f %8.6e %8.6e \n',    loc_i, loc_j, loc_k, loc_lon_model, loc_lat_model, loc_depth_model, loc_value_model, loc_value_data);
            else
                fprintf(fid, '%d %d %d %8.3f %8.3f %8.3f %8.6e %8.3f %8.3f %8.3f %8.6e %s \n', loc_i, loc_j, loc_k, loc_lon_model, loc_lat_model, loc_depth_model, loc_value_model, loc_lon_data, loc_lat_data, loc_depth_data, loc_value_data, loc_label);
            end
        end
        fclose(fid);
    elseif (data_saveall == 'y')
        fid = fopen([par_pathout '/' filename '_ALLMODELPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at mask locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: i, j, k, model lon, model lat, model depth, model value');
        fprintf(fid, '\n');
        for k = 1:kmax
            for j = 1:jmax
                loc_k = k;
                loc_j = j;
                loc_lat = xm(k,j);
                loc_depth = laym(k,j);
                loc_value = zm(k,j);
                fprintf(fid, '%2d %2d %8.3f %8.3f %8.6e %s \n', loc_j, loc_k, loc_lat, -loc_depth, loc_value, '%');
            end
        end
        fclose(fid);
    else
        fid = fopen([par_pathout '/' filename '_MODELPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at plotted locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: i, j, k, model lon, model lat, model depth, model value');
        fprintf(fid, '\n');
        for k = data_kmin:data_kmax
            for j = 1:jmax
                loc_k = k;
                loc_j = j;
                loc_lat = xm(k,j);
                loc_depth = laym(k,j);
                loc_value = zm(k,j);
                fprintf(fid, '%2d %2d %8.3f %8.3f %8.6e %s \n', loc_j, loc_k, loc_lat, -loc_depth, loc_value, '%');
            end
        end
        fclose(fid);
    end
end
%
% *** CREATE DATA VECTOR FOR HISTOGRAM ********************************** %
%
if s_in.slice > 0
    data_vector = reshape(data(:,s_in.slice,:),jmax*kmax,1);
    data_vector_V = reshape(data_V(:,s_in.slice,:),jmax*kmax,1);
else
    data_vector = reshape(data(:,:,:),imax*jmax*kmax,1);
    data_vector_V = reshape(data_V(:,:,:),imax*jmax*kmax,1);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** ANOMOLY PLOTTING DATA ADJUSTMENTS ********************************* %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % calculate model-data anomoly
    % NOTE: data_vector_2 == model, data_vector_1 == data
    if (data_anomoly == 'y')
        overlaydata(:,4) = data_vector_2(:) - data_vector_1(:);
    end
end
% redefine model grid location values so as to all appear white
if (s_in.overlay.only)
    zm = zeros(size(zm(:,:)));
    zm(find(zm(:,:) == 0)) = NaN;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT MAIN FIGURE ************************************************** %
% *********************************************************************** %
%
if (s_in.plot.showmain)
    %
    % *** CONFIGURE AND CREATE PLOTTING WINDOW ************************** %
    %
    % create figure
    % NOTE: explicitly specify renderer is using useless recent version
    scrsz = get(0,'ScreenSize');
    hfig = figure('Position',[((1 - plot_dscrsz)/2)*plot_dscrsz*scrsz(3) (1 - plot_dscrsz)*plot_dscrsz*scrsz(4) plot_dscrsz*scrsz(3) plot_dscrsz*scrsz(4)]);
    if (par_mutlab > 2015), hfig.Renderer='Painters'; end    
    clf;
    % define plotting regions
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.10 0.10 0.60 0.80]);
    fh(3) = axes('Position',[0.80 0.27 0.20 0.46],'Visible','off');
    % define colormap
    cmap = make_cmap(s_in.plot.colors,s_in.data.intervals+2);
    if (s_in.plot.invertcolors), cmap = flipdim(cmap,1); end,
    colormap(cmap);
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' v.', num2str(par_ver), ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    % *** SET PLOT SCALE ************************************************ %
    %
    % if min > max, then reverse min and max
    if s_in.data.min > s_in.data.max
        loc_minmax = [s_in.data.max s_in.data.min];
        s_in.data.min = loc_minmax(1);
        s_in.data.max = loc_minmax(2);
    end
    % set minimum contour value
    if (s_in.data.min == s_in.data.max)
        s_in.data.min = min(max(zm));
        s_in.data.max = max(max(zm));
    end
    % ensure min and max are not identical ...
    if (s_in.data.min == s_in.data.max)
        s_in.data.min = min(-1,s_in.data.min);
        s_in.data.max = max(1,s_in.data.max);
    end
    %
    % *** CREATE MAIN PLOT ********************************************** %
    %
    set(gcf,'CurrentAxes',fh(2));
    hold on;
    % set color and lat/lon axes and labels
    caxis([s_in.data.min-(s_in.data.max-s_in.data.min)/s_in.data.intervals s_in.data.max]);
    set(gca,'PlotBoxAspectRatio',[1.0 plot_xy_scaling*0.5 1.0]);
    if plot_global
        axis([lat_min lat_max -D_max -D_min]);
        set(gca,'XLabel',text('String','Latitude','FontSize',15),'XTick',[-90 -60 -30 0 30 60 90]);
        set(gca,'YLabel',text('String','Depth (km)','FontSize',15),'YTick',[-D_max:1000:-D_min],'YTickLabel',{'5';'4';'3';'2';'1';'0'});
    else
        axis([plot_lat_min plot_lat_max -plot_D_max -plot_D_min]);
        set(gca,'XLabel',text('String','Latitude','FontSize',15),'XTick',[plot_lat_min plot_lat_max]);
        set(gca,'YLabel',text('String','Depth (km)','FontSize',15),'YTick',[-plot_D_max -plot_D_min],'YTickLabel',{num2str(-plot_D_max);num2str(-plot_D_min)});
    end
    set(gca,'TickDir','out');
    if ~isempty(s_in.plot.title)
        title(s_in.plot.title,'FontSize',18);
    else
        if (~isempty(s_in.maskname))
            title(['Data ID: ',strrep(dataid_1,'_','-'),' / mask = ', strrep(s_in.maskname,'_','-')],'FontSize',12);
        else
            title(['Data ID: ',strrep(dataid_1,'_','-'),' / i = ', num2str(s_in.slice)],'FontSize',12);
        end
    end
    % assign dummy s_in.slice value
    if ((s_in.slice == 0) || (s_in.slice == -1)), s_in.slice = 1; end
    % draw filled rectangles
    for j = jmax:-1:1
        for k = 1:kmax
            if topo(j,s_in.slice) > ym(k,j)
                h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],color_g);
                set(h,'EdgeColor',color_g);
            else
                if (isnan(zm(k,j)))
                    if (~s_in.overlay.only)
                        h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],[1 1 1]);
                        set(h,'EdgeColor',[1 1 1]);
                    end
                else
                    if isempty(contour_file)
                        col = 1 + round(0.5+s_in.data.intervals*(zm(k,j)-s_in.data.min)/(s_in.data.max-s_in.data.min));
                        if col < 1, col = 1; end
                        if col > s_in.data.intervals+2, col = s_in.data.intervals+2; end
                    else
                        col = 1;
                        for c=1:s_in.data.intervals+1
                            if (zm(k,j) > contour_data(c)), col = c+1; end
                        end
                    end
                    if (~s_in.overlay.only)
                        h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],cmap(col,:));
                        set(h,'EdgeColor',cmap(col,:));
                    end
                end
            end
        end
    end
    %
    % *** PLOT CONTINENTAL OUTLINE ************************************** %
    %
    for k = 1:kmax
        for j = 1:jmax-1
            if topo(j,s_in.slice) > ym(k,j)
                if topo(j+1,s_in.slice) <= ym(k,j+1)
                    h = plot([latn(k,j) latn(k,j)],[layb(k,j) layt(k,j)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
        for j = 2:jmax
            if topo(j,s_in.slice) > ym(k,j)
                if topo(j-1,s_in.slice) <= ym(k,j-1)
                    h = plot([lats(k,j) lats(k,j)],[layb(k,j) layt(k,j)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
    end
    for j = 1:jmax
        for k = 2:kmax
            if topo(j,s_in.slice) < ym(k,j)
                if topo(j,s_in.slice) > ym(k-1,j)
                    h = plot([lats(k,j) latn(k,j)],[layb(k,j) layb(k,j)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
    end
    %
    % *** OVERLAY CONTOURS ********************************************** %
    %
    if ( s_in.contours.showlabels && ~s_in.overlay.only && ~s_in.opsi.show )
        if ( (s_in.data.min < 0.0) && (s_in.data.max > 0.0) && s_in.contours.dashneg )
            v = [0.0:s_in.contours.minor:s_in.data.max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.25);
            v = [0.0:s_in.contours.major:s_in.data.max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.5);
            if (s_in.contours.showlabels && ~s_in.data.log10)
                clabel(C,h);
            end
            v = [s_in.data.min:s_in.contours.minor:0.0];
            [C,h] = contour(xm,ym,zm,v,'k--');
            set(h,'LineWidth',0.25);
            v = [s_in.data.min:s_in.contours.major:0.0];
            [C,h] = contour(xm,ym,zm,v,'k--');
            set(h,'LineWidth',0.5);
            if (s_in.contours.showlabels && ~s_in.data.log10)
                clabel(C,h);
            end
        else
            v = [s_in.data.min:s_in.contours.minor:s_in.data.max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.25);
            v = [s_in.data.min:s_in.contours.major:s_in.data.max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.5);
            if (s_in.contours.showlabels && ~s_in.data.log10)
                clabel(C,h);
            end
        end
        % additional highlight contours
        loc_lim = 2.0*(abs(s_in.data.min) + abs(s_in.data.max));
        if (contour_hlt == 'y'),
            v = [-loc_lim+contour_hltval:loc_lim:loc_lim+contour_hltval];
            [C,h] = contour(xm,ym,zm,v,'w-');
            set(h,'LineWidth',1.0);
        end
        if (contour_hlt2 == 'y'),
            v = [-loc_lim+contour_hltval2:loc_lim:loc_lim+contour_hltval2];
            [C,h] = contour(xm,ym,zm,v,'w--');
            set(h,'LineWidth',1.0);
        end
    end
    %
    % *** OVERLAY CONTOURS -- OVERTURNING STREAMFUNCTION (IF SELECTED) ** %
    %
    if (s_in.opsi.show)
        if (s_in.overlay.only)
            s_in.contours.minor = s_in.opsi.min;
            s_in.contours.minor = s_in.opsi.max;
            s_in.data.intervals = (s_in.opsi.max - s_in.opsi.min)/s_in.opsi.minor;
            % re-define colormap
            cmap = make_cmap('anom',s_in.data.intervals+2);
            if (s_in.plot.invertcolors), cmap = flipdim(cmap,1); end,
            colormap(cmap);
            caxis([s_in.data.min-(s_in.data.max-s_in.data.min)/s_in.data.intervals s_in.data.max]);
            % contour
            v = [s_in.opsi.min:s_in.opsi.minor:s_in.opsi.max];
            [C,h] = contourf(opsigrid_lat,opsigrid_zt,opsizm,v);
            set(h,'LineColor','none')
        end
        v = [0.0:s_in.opsi.minor:s_in.opsi.max];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k-');
        set(h,'LineWidth',0.25);
        v = [0.0:s_in.opsi.major:s_in.opsi.max];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k');
        set(h,'LineWidth',0.5);
        if (s_in.contours.showlabels), clabel(C,h); end
        v = [s_in.opsi.min:s_in.opsi.minor:0.0];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k-.');
        set(h,'LineWidth',0.25);
        v = [s_in.opsi.min:s_in.opsi.major:0.0];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k-.');
        set(h,'LineWidth',0.5);
        if (s_in.contours.showlabels), clabel(C,h); end
    end
    %
    % *** OVERLAY DATA ************************************************** %
    %
    if ~isempty(overlaydataid)
        % set uniform marker shape and color
        if (data_shapecol == 'n'),
            for n = 1:nmax,
                overlaydata_shape(n) = 'o';
                overlaydata_ecol(n) = data_sitecolor;
                if ( s_in.overlay.only && (data_siteonly == 'y') )
                    overlaydata_fcol(n) = data_sitecolor;
                else
                    overlaydata_fcol(n) = '-';                    
                end
            end
        end
        % plot overlay data
        if (data_siteonly == 'n'),
            scatter(overlaydata(:,2),overlaydata(:,3),4,overlaydata(:,4),overlaydata_shape(n),'Filled','LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
        else
            if (overlaydata_fcol(n) == '-'),
                scatter(overlaydata(:,2),overlaydata(:,3),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
            else
                scatter(overlaydata(:,2),overlaydata(:,3),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n),'MarkerFaceColor',overlaydata_fcol(n));
            end
        end
    end
    %
    % *** PLOT BORDER *************************************************** %
    %
    h = plot([lat_min lat_max],[-D_max -D_max],'k-');
    set(h,'LineWidth',1.0);
    h = plot([lat_min lat_max],[-D_min -D_min],'k-');
    set(h,'LineWidth',1.0);
    h = plot([lat_min lat_min],[-D_max -D_min],'k-');
    set(h,'LineWidth',1.0);
    h = plot([lat_max lat_max],[-D_max -D_min],'k-');
    set(h,'LineWidth',1.0);
    %
    hold off;
    %
    % *** CREATE COLOR BAR ********************************************** %
    %
    if (~(s_in.overlay.only && (data_siteonly == 'y')))
        %
        set(gcf,'CurrentAxes',fh(3));
        hold on;
        %
        set(gca,'XTick',[],'YTick',[]);
        axis([0 1 0 s_in.data.intervals+2]);
        % draw and label color bar rectangles
        % draw and label start triangle
        c = 1;
        h = fill([0.1 0.2 0.3],[c c-1.0 c],cmap(c,:));
        if isempty(contour_file),
            str = [num2str(s_in.data.min + (c-1)*(s_in.data.max-s_in.data.min)/s_in.data.intervals)];
        else
            str = num2str(contour_data(c));
        end
        textsize = 6+round(80/s_in.data.intervals);
        if textsize > 12, textsize = 12; end
        text(0.40,c,str,'FontName','Arial','FontSize',textsize);
        set(h,'LineWidth',0.5);
        set(h,'EdgeColor','k');
        % draw and label bars
        for c = 2:s_in.data.intervals+1,
            h = fill([0.1 0.1 0.3 0.3],[c-1.0 c c c-1.0],cmap(c,:));
            if isempty(contour_file),
                str = [num2str(s_in.data.min + (c-1)*(s_in.data.max-s_in.data.min)/s_in.data.intervals)];
            else
                str = num2str(contour_data(c));
            end
            text(0.40,c,str,'FontName','Arial','FontSize',textsize);
            set(h,'LineWidth',0.5);
            set(h,'EdgeColor','k');
        end
        % draw end triangle
        c = s_in.data.intervals+2;
        h = fill([0.1 0.2 0.3],[c-1.0 c c-1.0],cmap(c,:));
        set(h,'LineWidth',0.5);
        set(h,'EdgeColor','k');
        %
        hold off;
        %
    end
    %
    % *** PRINT PLOT **************************************************** %
    %
    set(gcf,'CurrentAxes',fh(1));
    switch s_in.plotformat
        case {'pdf'}
            exportgraphics(gcf,[par_pathout '/' filename '.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
        case {'png','jpg','tif'}
            exportgraphics(gcf,[par_pathout '/' filename  '.' str_date '.' s_in.plotformat],'ContentType','image');
        otherwise
            % NO SAVED PLOT!
    end
    %
    % *** SAVE DATA ***************************************************** %
    %
    if (s_in.data.save), fprint_1Dn_d([data_vector_1 data_vector_D],[par_pathout '/' filename '.iDATA.', str_date, '.res']); end
    %
    % ******************************************************************* %
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SECONDARY FIGURES ************************************************* %
% *********************************************************************** %
%
if (s_in.plot.showextra)
    %
    % % *** SET PLOT SCALE ************************************************ %
    % %
    % % set minimum contour value
    % if exist('con_min','var') == 0
    %     con_min = min(min(zz));
    % end
    % % set maximum contour value
    % if exist('con_max','var') == 0
    %     con_max = max(max(zz));
    % end
    % % ensure min and max are not identical ...
    % if con_min == con_max
    %     if con_max == 0.0
    %         con_max = 1.0;
    %     else
    %         con_min = (1.00/1.01)*con_min;
    %         con_max = (1.01)*con_max;
    %     end
    % end
    % % if min > max, then reverse min and max
    % if con_min > con_max
    %     con_min_TEMP = con_min;
    %     con_max_TEMP = con_max;
    %     con_min = con_max_TEMP;
    %     con_max = con_min_TEMP;
    % end
    %
    % *** PLOT FIGURE (profile) ***************************************** %
    %
    % NOTE: dashed line (zl(:)) is the proifle of ALL model data,
    %       not just model cell that have corresponding data
    %      (i.e. it will not necessarily go through the blue triangles
    %       depending on obserevd data coverage)
    if ( ~s_in.overlay.only && (plot_profile == 'y'))
        %
        figure
        hold on;
        if ~isempty(overlaydataid)
            % plot mean model profile (blue dashed line)
            plot(zl(:),-grid_zt(:),'b--');
            set(h,'LineWidth',2.0);
            % plot data + corresponding model values
            % NOTE: for overlay data:
            %       data_vector_2 == model, data_vector_1 == data
            %       for 2 netCDF datasets:
            %       data_vector_2 == data2, data_vector_1 == data1
            scatter(data_vector_1(:),-grid_zt(data_vector_k(:)),10,'mo');
            scatter(data_vector_2(:),-grid_zt(data_vector_k(:)),10,'c^');
            %
            data_vector_12_mean = NaN(kmax,3);
            for k = data_kmin:data_kmax
                samelayer_locations = find((int32(data_vector_k(:))==k));
                samecell_n = size(samelayer_locations);
                if (samecell_n(1) > 0)
                    data_vector_12_mean(k,1) = -grid_zt(k);
                    samelayer_mean = mean(data_vector_1(samelayer_locations));
                    data_vector_12_mean(k,2) = samelayer_mean;
                    scatter(samelayer_mean,-grid_zt(k),50,'ro','filled');
                    samelayer_mean = mean(data_vector_2(samelayer_locations));
                    data_vector_12_mean(k,3) = samelayer_mean;
                    scatter(samelayer_mean,-grid_zt(k),50,'b^','filled');
                else
                    data_vector_12_mean(k,1) = -grid_zt(k);
                end
            end
        else
            %
            plot(zl(:),-grid_zt(:));
            set(h,'LineWidth',1.0);
            scatter(zl(:),-grid_zt(:),25,'r');
        end
        axis([s_in.data.min s_in.data.max -plot_D_max -plot_D_min]);
        xlabel(strrep(dataid_1,'_','-'));
        ylabel('Elevation (m)');
        if ~isempty(s_in.plot.title)
            title(s_in.plot.title,'FontSize',18);
        else
            if (~isempty(s_in.maskname))
                title(['Data ID: ',strrep(dataid_1,'_','-'),' / i = ', strrep(s_in.maskname,'_','-')],'FontSize',12);
            else
                title(['Data ID: ',strrep(dataid_1,'_','-'),' / i = ', num2str(s_in.slice)],'FontSize',12);
            end
        end
        exportgraphics(gcf,[par_pathout '/' filename '.PROFILE.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
        %
    end
    %
    % *** SAVE DATA (profile) ******************************************* %
    %
    % model mean profile
    if ( s_in.data.save && ~s_in.overlay.only && (plot_profile == 'y'))
        fprint_1D2_d([flipud(grid_zt(:)) flipud(zl(:))],[par_pathout '/' filename '.PROFILE.', str_date, '.res']); 
    end
    % data and modal @ data, mean profile
    if ( s_in.data.save && (plot_profile == 'y') && ~isempty(overlaydataid) )
        fid = fopen([par_pathout '/' filename '.PROFILE.MODELDATAMEANS.', str_date, '.res'], 'wt');
        fprintf(fid, '%% Mean data values + mean model values at data locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: k, depth (m), data mean, model mean, ocean volumn (m3) represented by that layer');
        fprintf(fid, '\n');
        for k = kmax:-1:1
            fprintf(fid, '%2d %8.3f %8.6e %8.6e %8.6e %s \n', k, -data_vector_12_mean(k,1), data_vector_12_mean(k,2), data_vector_12_mean(k,3), zl_V(k), '%');
        end
        fclose(fid);
    end
    %
    % *** PLOT FIGURE (surface zonal mean) ****************************** %
    %
    if ( ~s_in.overlay.only && (plot_zonal == 'y'))
        %
        figure
        plot(grid_lat,zz(kmax,:));
        hold on;
        scatter(grid_lat,zz(kmax,:),25,'r');
        axis([-90.0 90.0 s_in.data.min s_in.data.max ]);
        xlabel('Latitude');
        ylabel(strrep(dataid_1,'_','-'));
        exportgraphics(gcf,[par_pathout '/' filename '.ZONAL.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
        %
    end
    %
    % *** SAVE DATA (surface zonal mean) ******************************** %
    %
    if ( s_in.data.save && ~s_in.overlay.only && (plot_zonal == 'y'))
        fprint_1Dn_d([flipud(grid_lat) rot90(zz(kmax,:),1)],[par_pathout '/' filename '.ZONAL.', str_date, '.res']);
    end
    %
    % *** PLOT FIGURE (cross-plot) ************************************** %
    %
    % NOTE: data_vector_1 == DATA  (or model field 1)
    % NOTE: data_vector_2 == MODEL (or model field 2)
    %       => swap data2 and data1 so model (data1) is on y-axis
    %          as per for overlay data
    %          (assuming data2 are regridded observations)
    if ( ~isempty(dataid_2) || ~isempty(overlaydataid) )
        %
        if ~isempty(dataid_2)
            loc_x_data = reshape(data_2(data_kmin:data_kmax,:,:),[],1);
            loc_y_data = reshape(data_1(data_kmin:data_kmax,:,:),[],1);
            loc_D_data = reshape(data_D(data_kmin:data_kmax,:,:),[],1);
            loc_x_label = [strrep(dataid_2,'_','-')];
            loc_y_label = [strrep(dataid_1,'_','-')];
            loc_D_label = ['Depth (m)'];
        elseif ~isempty(overlaydataid)
            loc_x_data = data_vector_1;
            loc_y_data = data_vector_2;
            loc_D_data = -overlaydata(:,3);
            loc_x_label = [strrep(overlaydataid,'_','-')];
            loc_y_label = [strrep(dataid_1,'_','-')];
            loc_D_label = ['Depth (m)'];
        end
        % plot with and without depth coding
        % NOTE: test for insufficient data for scaling the plot
        % NOTE: function range has been moved ...
        if ((max(loc_x_data)-min(loc_x_data)) > 0.0)
            plot_crossplotc(loc_x_data,loc_y_data,[],loc_x_label,loc_y_label,'',POPT,[par_pathout '/' filename '.CROSSPLOT']);
            plot_crossplotc(loc_x_data,loc_y_data,loc_D_data,loc_x_label,loc_y_label,loc_D_label,POPT,[par_pathout '/' filename '.CROSSPLOTD']);
        end
        %
    end
    %
    % *** SAVE DATA (cross-plot relationships) ************************** %
    %
    if ( s_in.data.save && (~isempty(dataid_2) || (~isempty(overlaydataid) && ~s_in.overlay.only)) )
        fprint_1Dn_d([loc_x_data loc_y_data loc_D_data],[par_pathout '/' filename '.CROSSPLOT.', str_date, '.res']);
    end
    %
    % *** PLOT FIGURE (histogram) *************************************** %
    %
    if (plot_hist == 'y')
        % single histogram
        loc_bins1 = [s_in.data.min:(s_in.data.max-s_in.data.min)/s_in.data.intervals:s_in.data.max];
        str_name = [par_pathout '/' filename '.HIST1'];
        plot_histc_2d(data_vector_1,loc_bins1,strrep(dataid_1,'_','-'),[],[],'',plot_histc_SETTINGS,str_name);
        % double histogram
        str_name = [par_pathout '/' filename '.HIST2'];
        if (~isempty(dataid_2))
            loc_min = min(data_vector_2);
            loc_max = max(data_vector_2);
            loc_bins2 = [loc_min:(loc_max-loc_min)/10:loc_max];
            plot_histc_2d(data_vector_1,loc_bins1,strrep(dataid_1,'_','-'),data_vector_2,loc_bins2,[strrep(dataid_2,'_','-')],plot_histc_SETTINGS,[str_name 'D']);
        else
            loc_bins2 = fliplr(grid_zt_edges');
            loc_bins2(find(loc_bins2 < plot_D_min)) = [];
            loc_bins2(find(loc_bins2 > plot_D_max)) = [];
            plot_histc_2d(data_vector_1,loc_bins1,strrep(dataid_1,'_','-'),data_vector_D,loc_bins2,'Depth (m)',plot_histc_SETTINGS,str_name);
        end
    end
    %
    % ******************************************************************* %
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FUNCTION RETURN *************************************************** %
% *********************************************************************** %
%
% NOTE: 
%       STATM(1,2) = MEAN;
%       STATM(2,2) = STD;
%       STATM(3,2) = RMSD;
%       STATM(4,2) = CORRELATIONS;
%       STATM(5,2) = N;
%       STATM(6,2) = TOTAL RMSD;
%       STATM(7,2) = STANDARD DEVIATION NORMALISED RMSD;
%       STATM(8,2) = NORMALISED BIAS;
%       STATM(9,2) = R2;
%       STATM(10,2) = M;
%
% NOTE2:
% xds = datastats(x) returns statistics for the column vector x 
% to the structure xds. Fields in xds are:
% num       -- The number of data values
% max       -- The maximum data value
% min       -- The minimum data value
% mean      -- The mean value of the data
% median    -- The median value of the data
% range     -- The range of the data
% std       -- The standard deviation of the data
%
if (data_output_old == 'y')
    % return diagnostics
    DIAG = [z/z_V 1027.649*z];
    %
    if exist('STATM')
        OUTPUT = [STATM; DIAG];
    else
        OUTPUT = [DIAG];
    end
else
    % NOTE: use data_vector_1 which is the full grid values
    %       when there is no data
    % NOTE: remove NaNs first (also from depth vector)
    output = [];
    if (license('test','Curve Fitting Toolbox'))
        if (~isempty(overlaydataid) && (~s_in.overlay.only || (data_anomoly == 'y')))
            % basic data stats and those of corresponding model locations
            data_vector_1(find(isnan(data_vector_1))) = [];
            output.data = datastats(reshape(data_vector_1,[],1));
            output.data.sum  = sum(data_vector_1); % add sum
            data_vector_2(find(isnan(data_vector_2))) = [];
            output.model = datastats(reshape(data_vector_2,[],1));
            output.model.sum  = sum(data_vector_2); % add sum
        else
            % basic stats
            data_vector_1(find(isnan(data_vector_1))) = [];
            output.model = datastats(reshape(data_vector_1,[],1));
            output.model.sum  = sum(data_vector_1); % add sum
            % add inventory/global volumn-weighted mean
            output.model.inventory = 1027.649*z;
            output.model.mean      = z/z_V;
            % add old min,max
            output.model.min   = min(reshape(zm,[],1));
            output.model.max   = max(reshape(zm,[],1));
            % add MOC properties
            if (s_in.opsi.show)
                loc_opsi = opsizm(find(opsigrid_zt<(-plot_D_min)));
                output.moc_min = min(loc_opsi);
                output.moc_max = max(loc_opsi);
            end
        end
    else
        if (~isempty(overlaydataid) && (~s_in.overlay.only || (data_anomoly == 'y')))
            % basic data stats and those of corresponding model locations
            data_vector_1(find(isnan(data_vector_1))) = [];
            output.data.n    = length(data_vector_1);
            output.data.sum  = sum(data_vector_1);
            output.data.mean = mean(data_vector_1);
            output.data.min  = min(data_vector_1);
            output.data.max  = max(data_vector_1);
            data_vector_2(find(isnan(data_vector_2))) = [];
            output.model.n    = length(data_vector_2);
            output.model.sum  = sum(data_vector_2);
            output.model.mean = mean(data_vector_2);
            output.model.min  = min(data_vector_2);
            output.model.max  = max(data_vector_2);
        else
            output.model.n    = sum(~isnan(reshape(zm,[],1)));
            output.model.mean = mean(reshape(zm,[],1),1,'omitnan');
            output.model.min  = min(reshape(zm,[],1),1,'omitnan');
            output.model.max  = max(reshape(zm,[],1),1,'omitnan');
        end
    end
    % add profile stats (if selected)
    % NOTE: remove log10 transformation
    if (plot_profile == 'y')
        if (s_in.data.log10)
            zl(:) = 10.0.^zl(:);
        end
        output.profile.n    = length(zl);
        output.profile.mean = mean(zl,'omitnan');
        output.profile.min  = min(zl,[],'omitnan');
        output.profile.max  = max(zl,[],'omitnan');
        output.profile.minD = grid_zt(find(zl == output.profile.min));
        output.profile.maxD = grid_zt(find(zl == output.profile.max));
        output.profile.sur  = zl(end);
    end
    % add model-data/model stats
    if exist('STATM')
        output.statm          = STATM;
        output.statm_mean     = STATM(1,2);
        output.statm_std      = STATM(2,2);
        output.statm_rmsd     = STATM(3,2);
        output.statm_corr     = STATM(4,2);
        output.statm_n        = STATM(5,2);
        output.statm_tot_rmsd = STATM(6,2);
        output.statm_sdn_rmsd = STATM(7,2);
        output.statm_bias     = STATM(8,2);
        output.statm_r2       = STATM(9,2);
        output.statm_m        = STATM(10,2);
    end
    % set returned data
    OUTPUT = output;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% close netCDF files
netcdf.close(ncid_1);
if ~isempty(exp_2)
    netcdf.close(ncid_2);
end
if ( s_in.opsi.show && plot_opsi_calc == 'n')
    netcdf.close(ncid_0);
end
%
% *** clean up ****************************************************** %
%
% (optional) remove unpacked dir
if loc_flag_unpack
    disp(['    REMOVE DIR']);
    rmdir([data_dir],'s');
end
%
% *********************************************************************** %
