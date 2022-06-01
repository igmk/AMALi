function readamali_lukas(year,month,day)
%--------------------------------------------------------------------------
% read all amali files of given date <cdate> and save it in a matlab format
% version by Lukas 19 August 2011  based on readamali_pam2011
% created 040310 Holger Deckelmann
%
% Assuming 1700 bins and 6 channels in preallocation and in
% amalireadraw_lukas

% NOCH AENDERN: letztes file rausschmeissen, weil nicht sicher ist, dass
% der header sich waehrend des files nicht geaendert hat!
%--------------------------------------------------------------------------

fprintf('Reading lschmidtpath.\n')

%% Path settings

datadir=['/home/mech/projects/ac3/amali/preprocessor/data/' day '/' ]
matfilepath=['/home/mech/projects/ac3/amali/preprocessor/mat/' year '/'];

logfilename=[year month day '_log.txt'];
fid = fopen(logfilename, 'w');

%% Parameter Settings

datasets = 120;                                                        % Number of timesteps included in one file. (usually 120!)


%% Read AMALi files

% Get directory information, throw away . and ..
files=dir( [ datadir '*' ] )

if files(1).name =='.', files(1:end-1)=files(2:end); files=files(1:end-1); end
if files(1).name =='..', files(1:end-1)=files(2:end); files=files(1:end-1); end
[ cntfiles d ] = size(files);

hlp=datestr(now);
fprintf ('START: %s\n', hlp);
fprintf (fid, 'START: %s\r\n', hlp);
fprintf ('SETTINGS: number of datasets = %d\n', datasets);
fprintf (fid, 'SETTINGS: number of datasets = %d\r\n', datasets);
fprintf ('COMMENTS: ""\n');
fprintf (fid, 'COMMENTS: ""\r\n');


dm1=0;
i=1;
k=1;
while  (i<=cntfiles)
    % Initialize new file path.
    apath=[datadir files(i).name ];
    
    % First file of measurement?
    if (k==1)
        matfilename=[ matfilepath files(i).name(2:12) '.mat' ];
        fprintf ('-----------------------------------------------------------\n');
        fprintf (fid, '-----------------------------------------------------------\r\n');
        fprintf ('new matfile %s.mat\n\n', files(i).name(2:12));
        fprintf (fid, 'new matfile %s.mat\r\n\r\n', files(i).name(2:12));
        %      oldheader=amalireadheader(path);
        
        % preallocate variables
        alldata  = zeros(6,120*cntfiles,1700);
        allinfo  = zeros(4,120*cntfiles);
        shotinfo = zeros(6,120*cntfiles);
        
    end;
    fprintf (' now file %d / %d : %s\n',i,cntfiles,files(i).name);
    fprintf (fid, ' now file %d / %d : %s\r\n',i,cntfiles,files(i).name);
    
    %   % Detect any changes in header.
    %   header=amalireadheader(path); % this is the header of the next file!
    %   ok=1;
    %   lastcorrupt=0;
    %   if ( length(header) ~= length(oldheader)),
    %       ok=0;
    %       lastcorrupt=1;
    %       fprintf('header length changed!\n');
    %       fprintf(fid, 'header length changed!\r\n');
    %   elseif ( strcmp(header, oldheader) == 0 )
    %       ok=0;
    %       lastcorrupt=1;
    %       fprintf('header content changed!\n');
    %       fprintf(fid, 'header content changed!\r\n');
    %   end;
    ok = 1;
    % No header changes from i-1 to i? Read file i into variables.
    if ( ok==1 )
        [ headinfo trackinfo tm shots shots6 cntbins dbin rawdata errcountshot ] = amalireadraw_lukas(apath);
        fprintf(' sampled file %s\n',files(i).name);
        fprintf(fid, ' sampled file %s\r\n',files(i).name);
        fprintf(' counter "shots not 14" (532p 532s 355):  %3i %3i %3i %3i %3i %3i\n\n',errcountshot(:)); %shotfehler fuer jeden kanal
        fprintf(fid, ' counter "shots not 14" (532p 532s 355):  %3i %3i %3i %3i %3i %3i\r\n\r\n',errcountshot(:));
        
        % store all interesting data:
        [d1 dm2 d2]=  size(rawdata); % d1 Anzahl Kanaele dm2 Zeitschritte d2 Hoehenbins
        % time index is the 2nd. concatenate different sizes of time slices
        % into one variable alldata...
        alldata(1:d1,dm1+1:dm2+dm1,1:d2)=rawdata(1:d1,1:dm2,1:d2);
        
        % allinfo contains only information about last (6th) channel!
        allinfo(1,dm1+1:dm1+dm2) = tm;            % time???
        allinfo(2,dm1+1:dm1+dm2) = shots;
        allinfo(3,dm1+1:dm1+dm2) = cntbins;       % bins
        allinfo(4,dm1+1:dm1+dm2) = dbin;          % resolution [m]
        
        shotinfo(1:d1,dm1+1:dm1+dm2) = shots6;
        
        % Datafile is not complete? Indicate save and start new mat file.
        if ( dm2<datasets),
            ok=0;
            fprintf ('file %d incomplete (%d timesteps), measurement interrupted\n', i, dm2);
            fprintf (fid, 'file %d incomplete (%d timesteps), measurement interrupted\r\n', i, dm2);
        end
        
        dm1=dm2+dm1;  % set consequent index dm1,dm2
        % alldata is an array where the measured times go from 1 to dm2 over
        % all corresponding files
        k=k+1;
        i=i+1;
    end
    
    
    % End of this measurement? Any data to save? Save mat-file, delete vars. New sequence of measurements begins.
    if ( ok==0 && k>1),
        % Is the last file corrupt (header change)? Dump it.
        lastcorrupt = 0;
        if ( lastcorrupt == 1),
            fprintf ('file %s dumped (may be corrupt)\n', files(i-1).name);
            fprintf (fid, 'file %s\r dumped (may be corrupt)\n', files(i-1).name);
            if ( k>2 ),                                                             % More than 2 files read? Otherwise there is no data to save (first file is dumped, second read again).
                fprintf ('last file included: %s\n', files(i-2).name);
                fprintf (fid, 'last file included: %s\r\n', files(i-2).name);
            else
                fprintf ('last file included: %s\n', 'none');
                fprintf (fid, 'last file included: %s\r\n', 'none');
            end
            dm1=dm1-dm2;
        else
            fprintf ('last file included: %s\n', files(i-1).name);
            fprintf (fid, 'last file included: %s\r\n', files(i-1).name);
        end
        
        %fprintf('d1: %i  dm1:%i  d2: %i\n', d1, dm1, d2)
        
        % shrink variables (due to preallocation). assuming that last file is
        % representative for all.
        alldata   = alldata(1:d1,1:dm1,1:d2);
        allinfo   = allinfo(:,1:dm1);
        shotinfo  = shotinfo(1:d1,1:dm1);
        
        fprintf ('saving %s\n', matfilename);
        fprintf (fid, 'saving %s\r\n', matfilename);
        mkdir  matfilepath
        save (matfilename ,'headinfo', 'trackinfo', 'allinfo', 'alldata', 'shotinfo')
        % headinfo and trackinfo are valid for whole matfile, because
        % headers have been compared before
        clear alldata;
        clear allinfo;
        clear shotinfo;
        dm1=0;
        k=1;          % Indicate first round of new measurement.
    end
    
end

%% Create .mat file

% Any data to save?
if ( k > 1)
    % shrink variables (due to preallocation). assuming that last file is
    % representative for all.
    alldata   = alldata(1:d1,1:dm1,1:d2);
    allinfo   = allinfo(:,1:dm1);
    shotinfo  = shotinfo(1:d1,1:dm1);
    
    fprintf ('saving %s\n', matfilename);
    fprintf (fid, 'saving %s\r\n', matfilename);
    save (matfilename ,'headinfo', 'trackinfo', 'allinfo', 'alldata', 'shotinfo')
end
fprintf ('END: %s\n', datestr(now));
fprintf (fid, 'END: %s\r\n', datestr(now));
fclose(fid);

return
