%-----------------------------------------------------------
% read the complete information from amali binary files
% created 10.03.04 Holger Deckelmann
%
% added just +1 at k=ii+1:head.ntracks   19.04.11 Roland
% changed errcountshot: now it counts occurrence of shots6(jj,ii)~=14, for each
% channel separately. Removed break. Aug 2011 Lukas
%-----------------------------------------------------------

function [head, trackinfo,tm,shots,shots6,cntbins,dbin, rawdata, errcountshot]= amalireadraw_lukas(apath);

% open file with little endian format
[ fid errmsg ] = fopen(apath,'r','ieee-le' );

if ( fid < 0 )
    display ( [ 'ERROR opening file ' apath ' : ' errmsg ]);
    return
end

% read first 2 lines incl. blanks and CRLF

header=fgetl(fid);
[head.op header]=strtok(header);
[head.hw header]=strtok(header);
[head.sw header]=strtok(header);

header=fgetl(fid);
[x header]=strtok(header);
head.res=str2num(x);
[x header]=strtok(header);
head.qsdelay=str2num(x);
[x header]=strtok(header);
head.rate=str2num(x);
[x header]=strtok(header);
head.pretrigger=str2num(x);
[x header]=strtok(header);
head.angle=str2num(x);
[x header]=strtok(header);
head.ntracks=str2num(x);

for i=1:head.ntracks,
    header=fgetl(fid);
    [x header]=strtok(header);
    trackinfo(i).wvl=str2num(x);
    [trackinfo(i).pol header]=strtok(header);
    [x header]=strtok(header);
    trackinfo(i).uhigh=str2num(x);
    [x header]=strtok(header);
    trackinfo(i).ue=str2num(x);
    [trackinfo(i).tr header]=strtok(header);
end

lastNotValid=0;

% preallocate variables (saves time)
tm          =zeros(1,120);
shots       =zeros(1,120);
shots6      =zeros(head.ntracks,120);
cntbins     =zeros(1,120);
dbin        =zeros(1,120);
rawdata     =zeros(head.ntracks,120,1700);          % !!!Assuming 1700 bins
errcountshot =  zeros(1,6);
jj=0;
%read data lines
while ( ~feof(fid) )
    jj=jj+1;
    ii=0;
    for i=1:head.ntracks,
        ii=ii+1;
        % read record header 29 bytes
        binhead= fread (fid,31,'uchar');
        
        % End of file? Last (=this) record is incomplete. Shrink dimension by one, exit loop.
        if ( feof(fid))
            % lastNotValid = 1;
            % fprintf ('eof');
            jj=jj-1;
            break;
        end;
        
        header(1:29)=char(binhead(3:31));
        tm(jj)=str2num(header(1:13));
        shots(jj)=str2num(header(15:16));
        shots6(ii,jj)=str2num(header(15:16));
        cntbins(jj)=str2num(header(18:21));
        dbin(jj)=str2num(header(23:27));
        
        fmt=sprintf('%d*ushort',cntbins(jj));
        rawdata(ii,jj,1:cntbins(jj))=fread(fid,cntbins(jj),fmt);
        
        
        % Vergleich der Shotzahl mit vorheriger
        %if ( jj>1)
        %if ( shots(jj)~=shots(jj-1))
        % for k=ii+1:head.ntracks  % get rid of all tracks of this time
        % dummy= fread (fid,31,'uchar');
        % dummy=fread(fid,cntbins(jj),fmt);
        % end
        % disp( header(1:13) )
        % jj=jj-1   ;
        % lastNotValid=1;
        % errcountshot = errcountshot + 1;
        % break;
        % else
        % lastNotValid=0;
        % end
        % else
        % end;
        if (shots6(ii,jj)~=14)
            errcountshot(ii) = errcountshot(ii) + 1;
        else
        end;
        
        
        % plausibility tests
        %     if ( cntbins(jj) ~= 1700)                  % 1000 auf 1500 geaendert, AR 13.07.06
        %       for k=ii:head.ntracks  % get rid of all tracks of this time
        %          dummy= fread (fid,31,'uchar');
        %          dummy=fread(fid,cntbins(jj),fmt);
        %        end
        %        jj=jj-1;
        %        break;
        %     end; % dataline is not valid
        % read additional CRLF
    end;
    
    
end;
% Shrink variables to actual size
tm      =tm(1:jj);
shots   =shots(1:jj);
shots6  =shots6(:,1:jj);
cntbins =cntbins(1:jj);
dbin    =dbin(1:jj);
rawdata =rawdata(:,1:jj,:);



%if ( lastNotValid )   %%  Dimension verkuerzen, da letzter record ungueltig
%     rawtmp=rawdata(:,1:jj,:);
%     rawdata=rawtmp;
%end


% disp(['latest header: ' header(1:13)] )
fclose(fid);
return;

