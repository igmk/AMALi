%%function critterpath

addpath('/atm_meas/polar_5_6/amali/processing/nadir/cloud/Amali_subroutines');  


% addpath('/lidar4/lidar/amali2006/programme/matlab/astrid');        %eigene Programme
% addpath('/lidar4/lidar/amali2006/programme/matlab/xsignal');
% addpath('/lidar3/ozonuser/critter/inversion/');
% addpath('/lidar3/ozonuser/critter/inversion/matlab/');
% addpath('/lidar4/lidar/amali2006/programme/matlab/xcontour');
% addpath('/lidar4/lidar/amali2006/programme/matlab/wurstel');
% addpath('/lidar4/lidar/amali2006/programme/matlab/wavelet/wavelet_mw');
% addpath('/lidar4/lidar/amali2006/programme/matlab/wavelet');
% addpath('/lidar4/lidar/amali2006/programme/matlab/watervapor');
% addpath('/lidar4/lidar/amali2006/programme/matlab/trajecto');
% addpath('/lidar4/lidar/amali2006/programme/matlab/tools');
% addpath('/lidar4/lidar/amali2006/programme/matlab/tmp/nasaames');
% addpath('/lidar4/lidar/amali2006/programme/matlab/tmp');
% addpath('/lidar4/lidar/amali2006/programme/matlab/stextfun');
% addpath('/lidar4/lidar/amali2006/programme/matlab/statistics');
% addpath('/lidar4/lidar/amali2006/programme/matlab/so3anl');
% addpath('/lidar4/lidar/amali2006/programme/matlab/simulation');
% addpath('/lidar4/lidar/amali2006/programme/matlab/sav');
% addpath('/lidar4/lidar/amali2006/programme/matlab/rayoptic_tropo/current');
% addpath('/lidar4/lidar/amali2006/programme/matlab/rayoptic_tropo');
% addpath('/lidar4/lidar/amali2006/programme/matlab/prog_rs/earlinet');
% addpath('/lidar4/lidar/amali2006/programme/matlab/prog_rs/extinction');
% addpath('/lidar4/lidar/amali2006/programme/matlab/prog_rs');
% addpath('/lidar4/lidar/amali2006/programme/matlab/prog_mg');
% addpath('/lidar4/lidar/amali2006/programme/matlab/plot/auswertung_astar');
% addpath('/lidar4/lidar/amali2006/programme/matlab/plot/filter');
% addpath('/lidar4/lidar/amali2006/programme/matlab/plot/mpl');
% addpath('/lidar4/lidar/amali2006/programme/matlab/plot/georg');
% addpath('/lidar4/lidar/amali2006/programme/matlab/plot');
% addpath('/lidar4/lidar/amali2006/programme/matlab/physics');
% addpath('/lidar4/lidar/amali2006/programme/matlab/pez/pez31/docs');
% addpath('/lidar4/lidar/amali2006/programme/matlab/pez/pez31');
% addpath('/lidar4/lidar/amali2006/programme/matlab/pez');
% addpath('/lidar4/lidar/amali2006/programme/matlab/optics/sav');
% addpath('/lidar4/lidar/amali2006/programme/matlab/optics');
% addpath('/lidar4/lidar/amali2006/programme/matlab/netcdf');
% addpath('/lidar4/lidar/amali2006/programme/matlab/nasaames');
% addpath('/lidar4/lidar/amali2006/programme/matlab/myprograms/jens');
% addpath('/lidar4/lidar/amali2006/programme/matlab/myprograms');
% addpath('/lidar4/lidar/amali2006/programme/matlab/os');
% addpath('/lidar4/lidar/amali2006/programme/matlab/models');
% addpath('/lidar4/lidar/amali2006/programme/matlab/mlo-files');
% addpath('/lidar4/lidar/amali2006/programme/matlab/mathews');
% addpath('/lidar4/lidar/amali2006/programme/matlab/logo');
% addpath('/lidar4/lidar/amali2006/programme/matlab/lidarcalc');
% addpath('/lidar4/lidar/amali2006/programme/matlab/lidar');
% addpath('/lidar4/lidar/amali2006/programme/matlab/hexapod');
% addpath('/lidar4/lidar/amali2006/programme/matlab/guitools');
% addpath('/lidar4/lidar/amali2006/programme/matlab/graf_mg');
% addpath('/lidar4/lidar/amali2006/programme/matlab/graf/sav');
% addpath('/lidar4/lidar/amali2006/programme/matlab/graf');
% addpath('/lidar4/lidar/amali2006/programme/matlab/gps');
% addpath('/lidar4/lidar/amali2006/programme/matlab/games');
% addpath('/lidar4/lidar/amali2006/programme/matlab/fitfun');
% addpath('/lidar4/lidar/amali2006/programme/matlab/filter');
% addpath('/lidar4/lidar/amali2006/programme/matlab/evaluate/v2');
% addpath('/lidar4/lidar/amali2006/programme/matlab/evaluate');
% addpath('/lidar4/lidar/amali2006/programme/matlab/earlinet/case_2');
% addpath('/lidar4/lidar/amali2006/programme/matlab/earlinet');
% addpath('/lidar4/lidar/amali2006/programme/matlab/doc');
% addpath('/lidar4/lidar/amali2006/programme/matlab/digitalfilter');
% addpath('/lidar4/lidar/amali2006/programme/matlab/datetime');
% addpath('/lidar4/lidar/amali2006/programme/matlab/data/cira');
% addpath('/lidar4/lidar/amali2006/programme/matlab/data');
% addpath('/lidar4/lidar/amali2006/programme/matlab/constell');
% addpath('/lidar4/lidar/amali2006/programme/matlab/compatib');
% addpath('/lidar4/lidar/amali2006/programme/matlab/bat/quick');
% addpath('/lidar4/lidar/amali2006/programme/matlab/bat/nya_astar');
% addpath('/lidar4/lidar/amali2006/programme/matlab/bat');
% addpath('/lidar4/lidar/amali2006/programme/matlab/averages');
% addpath('/home/critter/matlab2/mein/matlab');  %%Ende
% %addpath('/lidar4/lidar/amali2006/programme/matlab/aenderungen'); %%Diese mu/3 das letzte file sein!




