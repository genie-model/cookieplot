function [] = plot_timeseries(s_in,str_in)
arguments
    s_in
    str_in = '';
end
% plot_timeseries
%
%   ***********************************************************************
%
%   ***********************************************************************
%   *** PLOT TIME-SERIES OUTPUT *******************************************
%   ***********************************************************************
%
%   plot_timeseries(s_in,str_in)
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   26/10/21    FORKED FROM plot_timeseries_biogem
%               *** VERSION 2.00 ******************************************
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% set version!
par_ver = 2.00;
% set function name
str_function = mfilename;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% set paths
str_current_path = pwd; % find current path
s_tmp = what('cookieplot'); % find path to cookieplot directory
str_cookieplot_path = s_tmp.path;
% add path to default settings
par_pathin = [str_cookieplot_path '/' 'SETTINGS'];
if ~(exist(par_pathin,'dir') == 7)
    disp([' * ERROR: Cannot find SETTINGS directory of cookieplot: ' par_pathin]); disp([' ']); return;
else
    addpath(par_pathin);
end
% determine whcih stupid version of MUTLAB we are using
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);
% initialize internal structure
s = [];
% define output naming ocnvention
s.str_resdir = 'results';
s.str_tsbg   = 'timeseries';
s.str_tssg   = '';
s.str_tsorb  = '';
%
% *** ingest original parameters as structure *************************** %
%
% load default plotting options -- OLD non-structure format
if ~isfield(s_in,'popt'), s_in.popt='settings_timeseries'; end; popt = s_in.popt;
eval(popt);
% load in optional settings -- NEW structure format
if ~isempty(str_in)
    if (exist([str_in],'file') || exist([str_in '.n'],'file'))
        eval(str_in);
    else
        disp([' * ERROR: Cannot find structure settings file: ' str_in]); disp([' ']); return;
    end
end
% identify missing s_in structure elements and create defaults
% set old format parameters equal to s_in element values
if ~isfield(s_in,'exp_1'),       s_in.exp_1 = ''; end;       expid1 = s_in.exp_1;
if ~isfield(s_in,'exp_2'),       s_in.exp_2 = ''; end;       expid2 = s_in.exp_2;
if ~isfield(s_in,'dataid_1'),    s_in.dataid_1 = ''; end;    data1_name = s_in.dataid_1;
if ~isfield(s_in,'dataid_2'),    s_in.dataid_2 = ''; end;    data2_name = s_in.dataid_2;
if ~isfield(s_in,'datacol_1'),   s_in.datacol_1 = NaN; end;  data1_n = s_in.datacol_1;
if ~isfield(s_in,'datacol_2'),   s_in.datacol_2 = NaN; end;  data2_n = s_in.datacol_2;
if ~isfield(s_in,'axistmin'),    s_in.axistmin = 0.0; end;   axis_tmin = s_in.axistmin;
if ~isfield(s_in,'axistmax'),    s_in.axistmax = 0.0; end;   axis_tmax = s_in.axistmax;
if ~isfield(s_in,'altfilename'), s_in.altfilename = ''; end; altfilename = s_in.altfilename;
% identify changes to default settings
if ~isfield(s_in,'axis1min'),   s_in.axis1min = 0.0; end; axis_data1_min = s_in.axis1min;
if ~isfield(s_in,'axis1max'),   s_in.axis1max = 0.0; end; axis_data1_max = s_in.axis1max;
if ~isfield(s_in,'axis1scl'),   s_in.axis1scl = 1.0; end; axis_data1_scale = s_in.axis1scl;
if ~isfield(s_in,'axis1str'),   s_in.axis1str = ''; end;  plot_data1_title = s_in.axis1str;
if ~isfield(s_in,'axis2min'),   s_in.axis2min = 0.0; end; axis_data2_min = s_in.axis2min;
if ~isfield(s_in,'axis2max'),   s_in.axis2max = 0.0; end; axis_data2_max = s_in.axis2max;
if ~isfield(s_in,'axis2scl'),   s_in.axis2scl = 1.0; end; axis_data2_scale = s_in.axis2scl;
if ~isfield(s_in,'axis2str'),   s_in.axis2str = '';  end; plot_data2_title = s_in.axis2str;
%
% *** parse structure settings ****************************************** %
%
% identify missing s_in structure elements and create defaults
% identify missing s loaded structure elements (and use s_in default)
if ~isfield(s_in,'plotformat'), s_in.plotformat = 'png'; end;     if ~isfield(s,'plotformat'), s.plotformat = s_in.plotformat; end
if ~isfield(s_in,'dispwarnings'), s_in.dispwarnings = false; end; if ~isfield(s,'dispwarnings'), s.dispwarnings = s_in.dispwarnings; end
if ~isfield(s_in,'log10time'), s_in.log10time = false; end;       if ~isfield(s,'log10time'), s.log10time = s_in.log10time; end
%
% *** backwards compatability ******************************************* %
%
% time scale offset parameter
if ~exist('plot_toffset','var'), plot_toffset = 0.0; end
% time scale direction parameter
if ~exist('plot_tdir','var'),    plot_tdir = 1.0; end
% paths
if ~exist('par_pathin','var'),   par_pathin   = 'cgenie_output'; end
if ~exist('par_pathlib','var'),  par_pathlib  = 'source'; end
if ~exist('par_pathout','var'),  par_pathout  = 'PLOTS'; end
if ~exist('par_pathdata','var'), par_pathdata = 'DATA'; end
if ~exist('par_pathmask','var'), par_pathmask = 'MASKS'; end
if ~exist('par_pathexam','var'), par_pathexam = 'EXAMPLES'; end
% data saving
if ~exist('opt_save_diagnostics','var'), opt_save_diagnostics = false; end
%
% *** MISC ************************************************************** %
%
% define grey color
color_g = [0.75 0.75 0.75];
% flags
data_orb1 = false;
data_orb2 = false;
%
% *** SET PATHS & DIRECTORIES ******************************************* %
%
% find where str_function lives ...
% remove its name (+ '.m' extension) from the returned path ...
str_function_path = which(str_function);
str_function_path = str_function_path(1:end-length(str_function)-3);
% check source code directory and add search path
if ~(exist([str_function_path '/' par_pathlib],'dir') == 7)
    disp([' * ERROR: Cannot find source directory']); disp([' ']); return;
else
    addpath([str_function_path '/' par_pathlib]);
end
% set input path
par_pathin = [str_current_path '/' par_pathin];
if ~(exist(par_pathin,'dir') == 7)
    disp([' * ERROR: Cannot find experiment results directory']); disp([' ']); return;
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
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD DATA ********************************************************* %
% *********************************************************************** %
%
% test for experiment
data_dir = [par_pathin '/' expid1];
if (exist(data_dir, 'dir') == 0)
    if (exist(par_pathin, 'dir') == 0)
        disp(['INFO: Path: ' par_pathin ' cannot be found.']);
    else
        disp(['INFO: Path: ' par_pathin ' exists.']);
        disp(['INFO: Experiment directory: ' expid1 ' cannot be found.']);
    end
    if (exist([data_dir '.tar.gz'],'file'))
        disp(['INFO: Archive: ' [data_dir '.tar.gz'] ' exists.']);
        disp([' * UN-PACKING ...']);
        untar([data_dir '.tar.gz'],par_pathin);
        loc_flag_unpack = true;
    else
        disp(['ERROR: Experiment cannot be found.']); disp([' ']); return;
    end
else
    loc_flag_unpack = false;
end
% load basic data
data_Tatm   = load([par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_climate_SAT.txt']);
data_seaice = load([par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_climate_seaice.txt']);
data_opsi   = load([par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_climate_opsi.txt']);
% test for age-enabled experiment
data_file = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_climate_ventilation_age.txt'];
if (exist(data_file, 'file') ~= 2)
    if s.dispwarnings, disp(['WARNING: No ventilation age output.']); end
    data_age(:,1)   = data_Tatm(:,1);
    data_age(:,2)   = NaN*data_age(:,1);
    flag_noage = true;
else
    data_age   = load([par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_climate_ventilation_age.txt']);
    flag_noage = false;
end
% test for an CO2-enabled experiment
data_file = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_atm_pCO2.txt'];
if (exist(data_file, 'file') ~= 2)
    if s.dispwarnings, disp(['WARNING: The experiment needs to be CO2-enabled in order to fully use this plotting function.']); end
    data_pCO2(:,1)   = data_Tatm(:,1);
    data_pCO2(:,2:3) = rand(length(data_Tatm(:,1)),2);
    flag_noCO2 = true;
else
    data_pCO2  = load([par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_atm_pCO2.txt'],'ascii');
    flag_noCO2 = false;
end
% test for an isotope-enabled experiment
data_file = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_atm_pCO2_13C.txt'];
if (exist(data_file, 'file') ~= 2)
    if s.dispwarnings, disp(['WARNING: The experiment needs to be 13C-enabled in order to fully use this plotting function.']); end
    data_pCO2_13C(:,1)   = data_Tatm(:,1);
    data_pCO2_13C(:,2:3) = rand(length(data_Tatm(:,1)),2);
    flag_no13CO2 = true;
else
    data_pCO2_13C = load([par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_atm_pCO2_13C.txt'],'ascii');
    flag_no13CO2  = false;
    flag_noCO2    = false;
end
% set time units
switch plot_tunits
    case 'kyr'
        data_pCO2(:,1)     = data_pCO2(:,1)*1.0E-3;
        data_pCO2_13C(:,1) = data_pCO2_13C(:,1)*1.0E-3;
        data_Tatm(:,1)     = data_Tatm(:,1)*1.0E-3;
        data_seaice(:,1)   = data_seaice(:,1)*1.0E-3;
    case 'Myr'
        data_pCO2(:,1)     = data_pCO2(:,1)*1.0E-6;
        data_pCO2_13C(:,1) = data_pCO2_13C(:,1)*1.0E-6;
        data_Tatm(:,1)     = data_Tatm(:,1)*1.0E-6;
        data_seaice(:,1)   = data_seaice(:,1)*1.0E-6;
end
% load optional data file #1
if ~isempty(data1_name)
    data1_file    = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_' data1_name '.txt'];
    data1_fileALT = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tssg '_' data1_name '.txt'];
    data1_fileORB = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsorb '_' data1_name '.txt'];
    if (exist(data1_file, 'file') == 2)
        data1 = load(data1_file,'ascii');
        [m,n] = size(data1);
        if (data1_n > n || data1_n < 2)
            disp(['ERROR: You must chose a column number between 2 and ' n ' in BIOGEM time-series file: ', data1_file]); disp([' ']); return;
        end
        % if 2nd experiment => load file and create anomoly
        % NOTE: interpolate 2nd dataset onto grid of 1st dataset if
        %       they are of different lengths (assuming this reflects different time points)
        if ~isempty(expid2)
            loc_data1_file = [par_pathin '/' expid2 '/' s.str_resdir '/' s.str_tsbg '_' data1_name '.txt'];
            if (exist(loc_data1_file, 'file') == 2)
                loc_data1 = load(loc_data1_file,'ascii');
                if (length(loc_data1(:,1)) ~= length(data1(:,1))),
                    % truncate the longer dataset
                    loc_tmax    = data1(end,1);
                    loc_tmaxloc = loc_data1(end,1);
                    if (loc_tmax > loc_tmaxloc),
                        data1(find(data1(:,1) > loc_tmaxloc),:) = [];
                    else
                        loc_data1(find(loc_data1(:,1) > loc_tmax),:) = [];
                    end
                    % interpolate
                    loc_data = interp1(loc_data1(:,1),loc_data1(:,2:end),data1(:,1),'cubic');
                    loc_data1 = [data1(:,1), loc_data];
                end
                % create anomoly
                % NOTE: EXP1 - EXP2
                data1(:,2:end) = data1(:,2:end) - loc_data1(:,2:end);
            else
                disp(['ERROR: BIOGEM time-series file ', loc_data1_file, ' does not exist.']); disp([' ']); return;
            end
        end
    elseif (exist(data1_fileALT, 'file') == 2)
        data1 = load(data1_fileALT,'ascii');
        [m,n] = size(data1);
        % NOTE: convert kyr to yr
        data1(:,1) = 1000.0*data1(:,1);
    elseif (exist(data1_fileORB, 'file') == 2)
        data1 = load(data1_fileORB,'ascii');
        [m,n] = size(data1);
        % set flag for 'orbit' output data
        data_orb1 = true;
    else
        disp(['ERROR: BIOGEM time-series file ', data1_file, ' does not exist.']);
        disp(['(or SEDGEM sedcorenv file ', data1_fileALT, ' does not exist)']);
        disp(['(or BIOGEM orbital data file ', data1_fileORB, ' does not exist)']);
        disp([' ']); return;
    end
    data1(:,data1_n) = axis_data1_scale*data1(:,data1_n);
    axis_data1_min = axis_data1_min;
    axis_data1_max = axis_data1_max;
    switch plot_tunits
        case 'kyr'
            data1(:,1) = data1(:,1)*1.0E-3;
        case 'Myr'
            data1(:,1) = data1(:,1)*1.0E-6;
    end
    if ~isempty(overlaydata1_file),
        if (exist(overlaydata1_file, 'file') == 2)
            overlaydata1 = load(overlaydata1_file,'ascii');
            overlaydata1(:,1) = overlaydata1(:,1) + plot_toffset;
        else
            disp(['ERROR: overlay data file ', overlaydata1_file, ' does not exist.']); disp([' ']); return;
        end
    end
else
    data1 = [];
    overlaydata1_file = [];
    data2_name = ''; % disable 2nd optional dataset
end
% load optional data file #2
if ~isempty(data2_name)
    data2_file = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsbg '_' data2_name '.txt'];
    data2_fileALT = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tssg '_' data2_name '.txt'];
    data2_fileORB = [par_pathin '/' expid1 '/' s.str_resdir '/' s.str_tsorb '_' data2_name '.txt'];
    if (exist(data2_file, 'file') == 2)
        data2 = load(data2_file,'ascii');
        [m,n] = size(data2);
        if (data2_n > n || data2_n < 2)
            disp(['ERROR: You must chose a column number between 2 and ' n ' in BIOGEM time-series file: ', data2_file]); disp([' ']); return;
        end
        % if 2nd experiment => load file, interpolate, and create anomoly
        if ~isempty(expid2)
            % load file
            loc_data2_file = [par_pathin '/' expid2 '/' s.str_resdir '/' s.str_tsbg '_' data2_name '.txt'];
            if (exist(loc_data2_file, 'file') == 2)
                loc_data2 = load(loc_data2_file,'ascii');
                if (length(loc_data2(:,1)) ~= length(data2(:,1))),
                    % truncate the longer dataset
                    loc_tmax    = data2(end,1);
                    loc_tmaxloc = loc_data2(end,1);
                    if (loc_tmax > loc_tmaxloc),
                        data2(find(data2(:,1) > loc_tmaxloc),:) = [];
                    else
                        loc_data2(find(loc_data2(:,1) > loc_tmax),:) = [];
                    end
                    % interpolate
                    loc_data = interp1(loc_data2(:,1),loc_data2(:,2:end),data2(:,1),'cubic');
                    loc_data2 = [data2(:,1), loc_data];
                end
                % create anomoly
                % NOTE: EXP1 - EXP2
                data2(:,2:end) = data2(:,2:end) - loc_data2(:,2:end);
            else
                disp(['ERROR: BIOGEM time-series file ', loc_data2_file, ' does not exist.']); disp([' ']); return;
            end
        end
    elseif (exist(data2_fileALT, 'file') == 2)
        data2 = load(data2_fileALT,'ascii');
        [m,n] = size(data2);
        % NOTE: convert kyr to yr
        data2(:,1) = 1000.0*data2(:,1);
    elseif (exist(data2_fileORB, 'file') == 2)
        data2 = load(data2_fileORB,'ascii');
        [m,n] = size(data2);
        % set flag for 'orbit' output data
        data_orb2 = true;
    else
        disp(['ERROR: BIOGEM time-series file ', data2_file, ' does not exist.']);
        disp(['(or SEDGEM sedcorenv file ', data2_fileALT, ' does not exist)']);
        disp(['(or BIOGEM orbital data file ', data2_fileORB, ' does not exist)']);
        disp([' ']); return;
    end
    data2(:,data2_n) = axis_data2_scale*data2(:,data2_n);
    axis_data2_min = axis_data2_min;
    axis_data2_max = axis_data2_max;
    switch plot_tunits
        case 'kyr'
            data2(:,1) = data2(:,1)*1.0E-3;
        case 'Myr'
            data2(:,1) = data2(:,1)*1.0E-6;
    end
    if ~isempty(overlaydata2_file),
        if (exist(overlaydata2_file, 'file') == 2)
            overlaydata2 = load(overlaydata2_file,'ascii');
            overlaydata2(:,1) = overlaydata2(:,1) + plot_toffset;
        else
            disp(['ERROR: overlay data file ', overlaydata2_file, ' does not exist.']); disp([' ']); return;
        end
    end
else
    data2 = [];
    overlaydata2_file = [];
end
%
% *** PROCESS DATA ****************************************************** %
%
% set length justin case ...
n_data = length(data_Tatm(:,1));
% set time-axis
% NOTE: correct to nearest integer
if (axis_tmin == axis_tmax)
    axis_tmin = int32(data_Tatm(1,1)-0.001);
    axis_tmax = int32(data_Tatm(end,1)+0.001);
end
if s.log10time
    axis_tmin = log10(max(1.0,axis_tmin));
    axis_tmax = log10(axis_tmax);
end
% set plot limits in terms of 'n'
% NOTE: double() to deal with min, max time limits inputted as integers
loc_data = data_Tatm;
loc_terr = min(abs(loc_data(:,1) - double(axis_tmin)));
axis_n_min = find(abs(loc_data(:,1) - double(axis_tmin)) == loc_terr);
loc_terr = min(abs(loc_data(:,1) - double(axis_tmax)));
axis_n_max = find(abs(loc_data(:,1) - double(axis_tmax)) == loc_terr);
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT DATA ********************************************************* %
% *********************************************************************** %
%
% *** DEFINE FIGURE ***************************************************** %
%
% set plot width scaling
if (opt_fatplot)
    loc_dx = plot_fatdx;
else
    loc_dx = 1.0;
end
% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[((1.0-plot_dscrsz)/2)*plot_dscrsz*scrsz(3) ((1.0-plot_dscrsz)/2)*plot_dscrsz*scrsz(4) 0.50*plot_dscrsz*scrsz(3) 0.90*plot_dscrsz*scrsz(4)]);
if (par_mutlab > 2015), hfig.Renderer='Painters'; end
clf;
% define plotting regions
fh(1) = axes('Position',[0 0 1 1],'Visible','on','XTick',[],'YTick',[],'Box','on');
fh(2) = axes('Position',[0.00 0.00 1.00 0.05],'Visible','off');
if (~flag_noCO2)
    fh(5) = axes('Position',[0.10 0.42 loc_dx*0.70 0.14],'XTickLabel',[]);
    if (~isempty(data1)), fh(4) = axes('Position',[0.10 0.26 loc_dx*0.70 0.14],'XTickLabel',[]); end
    if (~isempty(data2)), fh(3) = axes('Position',[0.10 0.10 loc_dx*0.70 0.14],'XTickLabel',[]); end
else
    if (~isempty(data1)), fh(5) = axes('Position',[0.10 0.42 loc_dx*0.70 0.14],'XTickLabel',[]); end
    if (~isempty(data2)), fh(4) = axes('Position',[0.10 0.26 loc_dx*0.70 0.14],'XTickLabel',[]); end
end
fh(6) = axes('Position',[0.10 0.58 loc_dx*0.70 0.14],'XTickLabel',[]);
fh(7) = axes('Position',[0.10 0.74 loc_dx*0.70 0.14],'XTickLabel',[]);
fh(8) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
fh(9) = axes('Position',[0.00 0.90 loc_dx*0.70 0.10],'Visible','off');
% date-stamp plot
set(gcf,'CurrentAxes',fh(1));
text(0.95,0.50,[str_function, ' : ', strrep(expid1,'_',' '), ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
% add title
set(gcf,'CurrentAxes',fh(9));
if ~isempty(plot_title)
    str_title = plot_title;
    str_title_sze = 15;
elseif ~isempty(altfilename)
    str_title = strrep(altfilename,'_',' ');
    str_title_sze = 15;
else
    str_title = ['Time-series: ',strrep(expid1,'_',' ')];
    str_title_sze = 12;
end
text(0.15,0.25,str_title,'FontName','Arial','FontSize',str_title_sze);
% set units
switch plot_tunits
    case 'kyr'
        str_tunits = 'kyr';
    case 'Myr'
        str_tunits = 'Myr';
    otherwise
        str_tunits = 'yr';
end
%
% *** PLOT PANEL #1 -- SURFACE CLIMATE ********************************** %
% 
% set up plot
set(gcf,'CurrentAxes',fh(7));
hold on;
% add time-scale offset/inversion
data_Tatm(:,1)   = plot_toffset + plot_tdir*data_Tatm(:,1);
data_seaice(:,1) = plot_toffset + plot_tdir*data_seaice(:,1);
%
% create plot on first y-axis
% auto-scale axes if limits not set
if (axis_Tatmmin == axis_Tatmmax)
    axis_Tatmmin = min(data_Tatm(:,2));
    axis_Tatmmax = max(data_Tatm(:,2));
end
% if all the values are the same ...
if (axis_Tatmmin == axis_Tatmmax)
    axis_Tatmmin = min(0.0,axis_Tatmmin);
    axis_Tatmmax = max(1.0,axis_Tatmmax);
end
% plot 1st y-axis data
if s.log10time
    hl1 = line(log10(data_Tatm(:,1)),data_Tatm(:,2),'Color','r','LineWidth',1.0);
else
    hl1 = line(data_Tatm(:,1),data_Tatm(:,2),'Color','r','LineWidth',1.0);
end
% set x-axis, 1st y-axis axes
ax1 = gca;
set(ax1,'XColor','k','TickDir','out','XTickLabel','');
set(ax1,'XLim',[axis_tmin axis_tmax]);
set(ax1,'YColor','r','TickDir','out');
set(ax1,'YLim',[axis_Tatmmin axis_Tatmmax]);
set(ax1,'YLabel',text('String',['T_{atm} (' char(176) 'C)'],'FontSize',10));
% create 2nd y-axis
if (axis_icemin == axis_icemax)
    axis_icemin = min(data_seaice(:,2));
    axis_icemax = max(data_seaice(:,2));
end
if (axis_icemin == axis_icemax)
    axis_icemin = 0.0;
    axis_icemax = 100.0;
end
% set 2nd axes
ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','b');
% plot 2nd y-axis data
if s.log10time
    hl2 = line(log10(data_seaice(:,1)),data_seaice(:,2),'Color','b','LineWidth',1.0,'Parent',ax2);
else
    hl2 = line(data_seaice(:,1),data_seaice(:,2),'Color','b','LineWidth',1.0,'Parent',ax2);
end
set(ax2,'YLim',[axis_icemin axis_icemax]);
set(ax2,'YLabel',text('String',['Seaice (%)'],'FontSize',10));
% create x-axis
set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabel','');
set(ax2,'XLim',[axis_tmin axis_tmax]);
%
% *** PLOT PANEL #2 -- DEEP OCEAN PROPERTIES **************************** %
%
% set up plot
set(gcf,'CurrentAxes',fh(6));
hold on;
% add time-scale offset
data_opsi(:,1) = plot_toffset + plot_tdir*data_opsi(:,1);
data_age(:,1)  = plot_toffset + plot_tdir*data_age(:,1);
%
% 1st y-axis
% auto-scale axes if limits not set
if (axis_opsimin == axis_opsimax)
    axis_opsimin = min( min(abs(data_opsi(:,2))),min(abs(data_opsi(:,3))) );
    axis_opsimax = max( max(abs(data_opsi(:,2))),max(abs(data_opsi(:,3))) );
end
% if all the values are still the same ...
if (axis_opsimin == axis_opsimax)
    axis_opsimin = min(0.0,axis_opsimin);
    axis_opsimax = max(1.0,axis_opsimax);
end
% plot 1st y-axis data
if s.log10time
    hl1 = line(log10(data_opsi(:,1)),abs(data_opsi(:,2)),'Color','r','LineWidth',1.0,'LineStyle','--');
    hl1 = line(log10(data_opsi(:,1)),abs(data_opsi(:,3)),'Color','r','LineWidth',1.0,'LineStyle','-');
else
    hl1 = line(data_opsi(:,1),abs(data_opsi(:,2)),'Color','r','LineWidth',1.0,'LineStyle','--');
    hl1 = line(data_opsi(:,1),abs(data_opsi(:,3)),'Color','r','LineWidth',1.0,'LineStyle','-');
end
% set x-axis, 1st y-axis axes
ax1 = gca;
set(ax1,'YColor','r','TickDir','out');
set(ax1,'YLim',[axis_opsimin axis_opsimax]);
set(ax1,'YLabel',text('String',['OPSI (Sv)'],'FontSize',10));
if (flag_noage)
    if (flag_noCO2 && isempty(data1))
        set(ax1,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabelMode','auto');
        set(ax1,'XLim',[axis_tmin axis_tmax]);
        if s.log10time
            set(ax1,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
        else
            set(ax1,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
        end
    else
        set(ax1,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabel','');
        set(ax1,'XLim',[axis_tmin axis_tmax]);
    end
end
%
% 2nd y-axis
if (~flag_noage)
    % auto-scale axes if limits not set
    if (axis_agemin == axis_agemax)
        axis_agemin = min(data_age(:,2));
        axis_agemax = max(data_age(:,2));
    end
    % if all the values are still the same ...
    if (axis_agemin == axis_agemax)
        axis_agemin = min(0.0,axis_agemin);
        axis_agemax = max(1.0,axis_agemax);
    end
    % set 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','b');
    % plot 2nd y-axis data
    hl2 = line(data_age(:,1),data_age(:,2),'Color','b','LineWidth',1.0,'Parent',ax2);
    if s.log10time
        hl1 = line(log10(data_age(:,1)),data_age(:,2),'Color','b','LineWidth',1.0,'Parent',ax2);
    else
        hl1 = line(data_age(:,1),data_age(:,2),'Color','b','LineWidth',1.0,'Parent',ax2);
    end
    % set 2nd y-axis axes
    set(ax2,'YLim',[axis_agemin axis_agemax]);
    set(ax2,'YLabel',text('String',['Benthic age (yrs)'],'FontSize',10));
    %
    % create x-axis (depending on number of panels)
    if (flag_noCO2 && isempty(data1))
        set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabelMode','auto');
        set(ax2,'XLim',[axis_tmin axis_tmax]);
        if s.log10time
            set(ax2,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
        else
            set(ax2,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
        end
    else
        set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabel','');
        set(ax2,'XLim',[axis_tmin axis_tmax]);
    end
end
%
% *** PLOT PANEL #3 -- CO2 [if selected] ******************************** %
%
if (~flag_noCO2)
    %
    % set up plot
    set(gcf,'CurrentAxes',fh(5));
    hold on;
    % add time-scale offset
    data_pCO2(:,1)     = plot_toffset + plot_tdir*data_pCO2(:,1);
    data_pCO2_13C(:,1) = plot_toffset + plot_tdir*data_pCO2_13C(:,1);
    %
    % 1st y-axis
    % auto-scale axes if limits not set
    if (axis_pCO2min == axis_pCO2max)
        axis_pCO2min = 1.0E6*min(data_pCO2(:,3));
        axis_pCO2max = 1.0E6*max(data_pCO2(:,3));
    end
    % if all the values are still the same ...
    if (axis_pCO2min == axis_pCO2max)
        axis_pCO2min = min(0.0,axis_pCO2min);
        axis_pCO2max = max(1.0,axis_pCO2max);
    end
    % plot 1st y-axis data
    if s.log10time
        hl1 = line(log10(data_pCO2(:,1)),1.0E6*data_pCO2(:,3),'Color','r','LineWidth',1.0);
    else
        hl1 = line(data_pCO2(:,1),1.0E6*data_pCO2(:,3),'Color','r','LineWidth',1.0);
    end
    % set x-axis, 1st y-axis axes
    ax1 = gca;
    set(ax1,'YColor','r','TickDir','out');
    set(ax1,'YLim',[axis_pCO2min, axis_pCO2max]);
    set(ax1,'YLabel',text('String',['pCO_{2} (\muatm)'],'FontSize',10));
    if (flag_no13CO2)
        if (isempty(data1))
            set(ax1,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabelMode','auto');
            set(ax1,'XLim',[axis_tmin axis_tmax]);
            if s.log10time
                set(ax1,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
            else
                set(ax1,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
            end
        else
            set(ax1,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabel','');
            set(ax1,'XLim',[axis_tmin axis_tmax]);
        end
    end
    %
    % 2nd y-axis
    if (~flag_no13CO2)
        % 2nd y-axis
        % auto-scale axes if limits not set
        if (axis_d13Cmin == axis_d13Cmax)
            axis_d13Cmin = min(data_pCO2_13C(:,3));
            axis_d13Cmax = max(data_pCO2_13C(:,3));
        end
        % if all the values are still the same ...
        if (axis_d13Cmin == axis_d13Cmax)
            axis_d13Cmin = min(0.0,axis_d13Cmin);
            axis_d13Cmax = max(1.0,axis_d13Cmax);
        end
        % create 2nd y-axis
        ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','b');
        % plot 2nd y-axis data
        if s.log10time
            hl2 = line(log10(data_pCO2_13C(:,1)),data_pCO2_13C(:,3),'Color','b','LineWidth',1.0,'Parent',ax2);
        else
            hl2 = line(data_pCO2_13C(:,1),data_pCO2_13C(:,3),'Color','b','LineWidth',1.0,'Parent',ax2);
        end
        % set 2nd y-axis axes
        set(ax2,'YLim',[axis_d13Cmin, axis_d13Cmax]);
        set(ax2,'YLabel',text('String',['\delta^{13}C (' char(8240) ')'],'FontSize',10));
    end
    %
    % create x-axis (depending on number of panels)
    if (isempty(data1))
        set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabelMode','auto');
        set(ax2,'XLim',[axis_tmin axis_tmax]);
        if s.log10time
            set(ax2,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
        else
            set(ax2,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
        end
    else
        set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabel','');
        set(ax2,'XLim',[axis_tmin axis_tmax]);
    end
end
%
% *** PLOT PANEL #4 (OPTIONAL DATA #1) ********************************** %
%
if (~isempty(data1))    
    % add time-scale offset
    data1(:,1) = plot_toffset + plot_tdir*data1(:,1);
    %
    if (flag_noCO2)
        set(gcf,'CurrentAxes',fh(5));
    else
        set(gcf,'CurrentAxes',fh(4));
    end
    hold on;
    %
    % 1st y-axis
    % auto-scale axes if limits not set
    if (axis_data1_min >= axis_data1_max)
        axis_data1_min = min(data1(axis_n_min:axis_n_max,data1_n));
        axis_data1_max = max(data1(axis_n_min:axis_n_max,data1_n));
    end
    if (axis_data1_min == axis_data1_max)
        axis_data1_min = min(0.0,axis_data1_min);
        axis_data1_max = max(1.0,axis_data1_max);
    end
    % plot data
    if s.log10time
        hl1 = line(log10(data1(:,1)),data1(:,data1_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(log10(data1(:,1)),data1(:,data1_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata1_file)), hp1 = scatter(log10(overlaydata1(:,1)),overlaydata1(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    else
        hl1 = line(data1(:,1),data1(:,data1_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(data1(:,1),data1(:,data1_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata1_file)), hp1 = scatter(overlaydata1(:,1),overlaydata1(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    end
end
% AXES
if (~isempty(data1))
    ax1 = gca;
    % create and label y-axis
    set(ax1,'YColor','k','TickDir','out');
    if (axis_data1_max > axis_data1_min), set(ax1,'YLim',[axis_data1_min axis_data1_max]); end
    if (~isempty(data1))
        if (isempty(plot_data1_title))
            plot_data1_title = data1_name;
            plot_data1_title(find(plot_data1_title(:)=='_')) = '.';
            set(ax1,'YLabel',text('String',[plot_data1_title],'FontSize',10));
        else
            set(ax1,'YLabel',text('String',[plot_data1_title],'FontSize',12));
        end
    else
        set(ax1,'YLabel',text('String','Emissions (PgC yr^{-1})','FontSize',10));
    end
    % create and label x-axis depending on number of panels
    if (isempty(data2))
        set(ax1,'XColor','k','TickDir','out','XTickLabelMode','auto');
        set(ax1,'XLim',[axis_tmin axis_tmax]);
        if s.log10time
            set(ax1,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
        else
            set(ax1,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
        end
    else
        set(ax1,'XColor','k','TickDir','out','XTickLabel','');
        set(ax1,'XLim',[axis_tmin axis_tmax]);
    end
end
%
% *** PLOT PANEL #5 (OPTIONAL DATA #2) ********************************** %
%
if (~isempty(data2))
    % add time-scale offset/inversion
    data2(:,1) = plot_toffset + plot_tdir*data2(:,1);
    %
    if (flag_noCO2)
        set(gcf,'CurrentAxes',fh(4));
    else
        set(gcf,'CurrentAxes',fh(3));
    end
    hold on;
    % set axes limits
    % NOTE: auto-scale only within the requested plotting limits
    if (axis_data2_min >= axis_data2_max)
        axis_data2_min = min(data2(axis_n_min:axis_n_max,data2_n));
        axis_data2_max = max(data2(axis_n_min:axis_n_max,data2_n));
    end
    if (axis_data2_min == axis_data2_max)
        axis_data2_min = (1.0 - sign(axis_data2_min)/1000.0)*axis_data2_min;
        axis_data2_max = (1.0 + sign(axis_data2_max)/1000.0)*axis_data2_max;
    end
    % plot data
    if s.log10time
        hl1 = line(log10(data2(:,1)),data2(:,data2_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(log10(data2(:,1)),data2(:,data2_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata2_file)), hp1 = scatter(log10(overlaydata2(:,1)),overlaydata2(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    else
        hl1 = line(data2(:,1),data2(:,data2_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(data2(:,1),data2(:,data2_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata2_file)), hp1 = scatter(overlaydata2(:,1),overlaydata2(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    end
end
% AXES
if (~isempty(data2))
    ax1 = gca;
    set(ax1,'XLim',[axis_tmin axis_tmax]);
    set(ax1,'XColor','k','XTick',[]);
    set(ax1,'YLim',[axis_data2_min axis_data2_max]);
    set(ax1,'YColor','k','YTick',[]);
    % create 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','k');
    % create and label y-axis
    set(ax2,'YColor','k','TickDir','out','YTickLabelMode','auto');
    if (axis_data2_max > axis_data2_min), set(ax2,'YLim',[axis_data2_min axis_data2_max]); end
    if (~isempty(data2))
        if (isempty(plot_data2_title))
            plot_data2_title = data2_name;
            plot_data2_title(find(plot_data2_title(:)=='_')) = '.';
            set(ax2,'YLabel',text('String',[plot_data2_title],'FontSize',10));
        else
            set(ax2,'YLabel',text('String',[plot_data2_title],'FontSize',12));
        end
    else
        set(ax2,'YLabel',text('String','\Sigma(emissions) (PgC)','FontSize',10));
    end
    % create and label x-axis
    set(ax2,'XColor','k','TickDir','out','XTickLabelMode','auto');
    set(ax2,'XLim',[axis_tmin axis_tmax]);
    if s.log10time
        set(ax2,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
    else
        set(ax2,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
    end
end
%
% *** PRINT PLOT ******************************************************** %
%
% set filename
str_filename = ['timeseries.' expid1];
if (~isempty(altfilename)), str_filename = altfilename; end
if (~isempty(expid2)), str_filename = [expid1 '_MINUS_' expid2]; end
str_filename = [str_filename '.' str_date];
% save plot
set(gcf,'CurrentAxes',fh(1));
switch s_in.plotformat
    case {'pdf'}
        exportgraphics(gcf,[par_pathout '/' str_filename '.pdf'],'BackgroundColor','none','ContentType','vector');
    case {'png','jpg','tif'}
        exportgraphics(gcf,[par_pathout '/' str_filename '.' s_in.plotformat],'ContentType','image');
    otherwise
        % NO SAVED PLOT!
        %%%exportgraphics(gcf,[par_pathout '/' str_filename '.' s_in.plotformat],'ContentType','image','Width', 11.5, 'Height', 8.0, 'Units', 'inches', 'PreserveAspectRatio', 'on');
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** EXTRACT DIAGNOSTIC  DATA ****************************************** %
% *********************************************************************** %
%
if (opt_save_diagnostics)
    %
    % *** \/\/\/\/ ********************************************************** %
    %
    loc_str_data = struct('name', {}, 'time', {}, 'value', {});
    %
    loc_data = data_Tatm;
    loc_terr = min(abs(loc_data(:,1) - axis_tmin));
    loc_n_min = find(abs(loc_data(:,1) - axis_tmin) == loc_terr);
    loc_terr = min(abs(loc_data(:,1) - axis_tmax));
    loc_n_max = find(abs(loc_data(:,1) - axis_tmax) == loc_terr);
    %
    % *** DATA #1 *********************************************************** %
    %
    loc_str = 'pCO2';
    loc_data = data_pCO2;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,3));
    loc_data_min_t = loc_data(find(loc_data(:,3) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,3));
    loc_data_max_t = loc_data(find(loc_data(:,3) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {1}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {1}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {1}, 'value', loc_data(loc_n_min,3));
    loc_str_data = setfield(loc_str_data, {2}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {2}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {2}, 'value', loc_data(loc_n_max,3));
    loc_str_data = setfield(loc_str_data, {3}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {3}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {3}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {4}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {4}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {4}, 'value', loc_data_max);
    %
    % *** DATA #2 *********************************************************** %
    %
    loc_str = 'pCO2_13C';
    loc_data = data_pCO2_13C;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,3));
    loc_data_min_t = loc_data(find(loc_data(:,3) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,3));
    loc_data_max_t = loc_data(find(loc_data(:,3) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {5}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {5}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {5}, 'value', loc_data(loc_n_min,3));
    loc_str_data = setfield(loc_str_data, {6}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {6}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {6}, 'value', loc_data(loc_n_max,3));
    loc_str_data = setfield(loc_str_data, {7}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {7}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {7}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {8}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {8}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {8}, 'value', loc_data_max);
    %
    % *** DATA #3 *********************************************************** %
    %
    loc_str = 'Tatm';
    loc_data = data_Tatm;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,2));
    loc_data_min_t = loc_data(find(loc_data(:,2) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,2));
    loc_data_max_t = loc_data(find(loc_data(:,2) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {9}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {9}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {9}, 'value', loc_data(loc_n_min,2));
    loc_str_data = setfield(loc_str_data, {10}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {10}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {10}, 'value', loc_data(loc_n_max,2));
    loc_str_data = setfield(loc_str_data, {11}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {11}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {11}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {12}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {12}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {12}, 'value', loc_data_max);
    %
    % *** DATA #4 *********************************************************** %
    %
    loc_str = 'seaice';
    loc_data = data_seaice;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,3));
    loc_data_min_t = loc_data(find(loc_data(:,3) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,3));
    loc_data_max_t = loc_data(find(loc_data(:,3) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {13}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {13}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {13}, 'value', loc_data(loc_n_min,3));
    loc_str_data = setfield(loc_str_data, {14}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {14}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {14}, 'value', loc_data(loc_n_max,3));
    loc_str_data = setfield(loc_str_data, {15}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {15}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {15}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {16}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {16}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {16}, 'value', loc_data_max);
    %
    % *** DATA #5 (OPTIONAL DATA #1) **************************************** %
    %
    if (~isempty(data1) && ~opt_invanalysis)
        loc_str = strrep(plot_data1_title,' ','_');
        loc_data = data1;
        loc_data_min = min(loc_data(loc_n_min:loc_n_max,data1_n));
        loc_data_min_t = loc_data(find(loc_data(:,data1_n) == loc_data_min),1);
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max,data1_n));
        loc_data_max_t = loc_data(find(loc_data(:,data1_n) == loc_data_max),1);
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {17}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {17}, 'time', loc_data(loc_n_min,1));
        loc_str_data = setfield(loc_str_data, {17}, 'value', loc_data(loc_n_min,data1_n));
        loc_str_data = setfield(loc_str_data, {18}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {18}, 'time', loc_data(loc_n_max,1));
        loc_str_data = setfield(loc_str_data, {18}, 'value', loc_data(loc_n_max,data1_n));
        loc_str_data = setfield(loc_str_data, {19}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {19}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {19}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {20}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {20}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {20}, 'value', loc_data_max);
    end
    %
    % *** DATA #6 (OPTIONAL DATA #2) **************************************** %
    %
    if (~isempty(data2) && ~opt_invanalysis)
        loc_str = strrep(plot_data2_title,' ','_');
        loc_data = data2;
        loc_data_min = min(loc_data(loc_n_min:loc_n_max,data2_n));
        loc_data_min_t = loc_data(find(loc_data(:,data2_n) == loc_data_min),1);
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max,data2_n));
        loc_data_max_t = loc_data(find(loc_data(:,data2_n) == loc_data_max),1);
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {21}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {21}, 'time', loc_data(loc_n_min,1));
        loc_str_data = setfield(loc_str_data, {21}, 'value', loc_data(loc_n_min,data2_n));
        loc_str_data = setfield(loc_str_data, {22}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {22}, 'time', loc_data(loc_n_max,1));
        loc_str_data = setfield(loc_str_data, {22}, 'value', loc_data(loc_n_max,data2_n));
        loc_str_data = setfield(loc_str_data, {23}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {23}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {23}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {24}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {24}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {24}, 'value', loc_data_max);
    end
    %
    % *** SAVE DATA ********************************************************* %
    %
    if (par_mutlab >= 2014)
        loc_table = struct2table(loc_str_data);
        if ~isempty(altfilename)
            writetable(loc_table,[par_pathout '/' altfilename '.txt'],'Delimiter',' ');
        else
            writetable(loc_table,[par_pathout '/' str_filename '.txt'],'Delimiter',' ');
        end
    end
    %
    % *** /\/\/\/\ ********************************************************** %
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
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
