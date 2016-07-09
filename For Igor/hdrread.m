function info = hdrread(filename)
% function info = hdrread(filename)
%
% loads header information from *.hdr file saved by Yokogawa DL series
% Oscilloscopes. The structure "info" contains the header information
% with name convention as in Yokogawa ASCII header file format
% description in Appendix 3 of the DL 1640 user manual
%
% optional input arguments:
% Group, Trace, Block: specifies which trace should be read (omitted values
%  are set to 1)
% datapoints, step,startind:
%  read 'datapoints' values in steps of 'step', starting at position 'startind'
%  from the wvf file. If datapoints is omitted, all data are read, if step 
%  is omitted, step=1. If startind omitted, startind=1
%
% Reading of *.wvf files written by other than the above Oscilloscopes may
% result in errors, because I could not test with those wvf files.
%
% Author:
% Erik Benkler
% Physikalisch-Technische Bundesanstalt
% Section 4.53: Optical Femtosecond Metrology
% Bundesallee 100
% D-38116 Braunschweig
% Germany
% Erik.Benkler a t ptb.de
%
% The implementation is based on Yokogawa ASCII header file format
% description in Appendix 3 of the DL 1640 user manual and on the
% Yokogawa Technical Information TI 7000-21 E 1998 1 3rd edition:
% "Understanding the structure of Binary data (xxxxxxxx.WVF) file created
% by the DL, AR series"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% current state of the project and change history:
%
% Version 1.0, May 21, 2007: first released version tested with DL 1600 scopes
% requires function stringparse.m by E. Benkler
%
% Version 1.1, July 31, 2008:
% (a) adapted to work with DL 7440 scopes (writes empty lines at different positions of header file)
% (b) requires functions stringparse.m and fileparse.m by E. Benkler
%
% Version 1.2, September 24, 2008
% (a) added compatibility with file format written by DL750 scopes
% (b) moved stringparse and fileparse functions to hdrread.m (old
% stringparse.m and fileparse.m should be deleted from Matlab path)
%
% Version 1.3, January 26, 2009
% (a) removed a bug (missing initialization of variable "rep" when using
%     step parameter (thanks to Peter Noronha for pointing out this bug)
% (b) improved performance when using step argument by preallocating data
%     array (thanks to Peter Noronha for pointing out this improvement)
% (c) added support for logic signals (raw binary format) written by DL750
%     model only.
% (d) changed the way the header file is parsed. Only if an entry really
%     exists and makes sense (e.g., a question mark as the entry
%     VPlusOverData does not make sense), the corresponding info structure
%     element is generated.
% (e) print out warning concerning wrong values only if output of these
%     values is requested by presence of the corresponding output argument
%     of wvfread
% (f) info structure added as optional output argument of wvfread
% (g) The $PrivateInfo is not read any more, because this contains too much
%     instrument-specific stuff which is probably not really necessary.
%     Please inform the author if you need reading of the $PrivateInfo
% (h) Error message if there are multiple blocks in the wvf file
% (i) Replaced "eval" by using dynamic structure entries
%
% Version 1.4, April 30, 2009
% (a) improved file name checking (re-submitted to FileExchange)
%
% Version 1.5, June 04, 2009: Major changes!
% zip-file contents: hdrread.m, wvfread.m, wvfreadb.m, wvfreadall.m,
% wvf__logicsignals_example.m
% (a) Added function wvfreadb.m (wvfread with block>1 support)
%     This required modification of the input parameter structure, so, in
%     order to maintain downward compatibility of wvfread.m,  I
%     wrote the new function wvfreadb.m which takes the additional input
%     parameter "Block".
%     wvfread.m is now a wrapper function with its old in- and
%     output structures which calls wvfreadb.m and from here on, all
%     changes are maintained in wvfreadb.m
% (b) Added support for accessing any block in a trace and group:
%     wvfreadb.m can now read any specified block in a trace and group,
%     whereas wvfread.m remains limited to reading of block 1 only.
% (c) Added support for wvf files written in the "Block" format (untested,
%     please report to the author if you can provide a test wvf file in the
%     block format): both wvfread.m and wvfreadb.m now should work with wvf
%     files written in the block format (but wvfread.m is restricted to
%     block number 1, whereas wvfreadb.m can access any specified block)
% (d) Removed a bug in hdrread.m: Entry "TraceName" is now read to
%     tracenummax, not only to Group1.trace1.tracenum(thanks to Frederic
%     Sun and Ronan Danioux for pointing out this bug)
% (e) Removed function "assign"
% (f) Added the wrapper function wvfreadall.m to read all or selected
%     data from a wvf file in a call to the single function wvfreadall.m
% (g) Added sample file "wvf_logicsignals_example.m" explaining access of
%     logic data traces written by the DL7XX scopes.
%
% Version 1.6, 08 December 2009,
% zip-file contents: hdrread.m, wvfread.m, wvfreadb.m, wvfreadall.m,
% wvf__logicsignals_example.m
% (a) Adapted hdrread.m to correctly read files written by newer (verified
%     with model version 6.22) DL750 scopes with muultiple blocks in the
%     "trace" format (some older versions generate as many date entries as
%     blocks exist, the newer models at least sometimes seem to generate
%     only one entry). Thanks to Mathias Graf for reporting this issue.
%
% Version 1.7, 22 March 2011,
%  (a) Removed bug related to optional input parameters datapoints, step
%      and startin. datapoints now is the number of points requested to be
%      returned. Added a warning when requested number
%      is larger than possible due to number of points in the trace.
%      Check ranges of input parameters startind, step, datapoints, and set
%      them to 1 if <1. Added warning if this is done.
%      Thanks to Matthias Hunstig, Univ. Paderborn, for pointing out these
%      issues.
%  (b) Removed bug in wvfread wrapper function regarding varargout
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% beginning of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%checking of file name etc.
if nargin==0
    filename='';
end

if isempty(filename)
    [fname,pname]=uigetfile({'*.hdr', 'Yokogawa Waveform Header Files (*.hdr)'},'Choose Yokogawa waveform header file');
    filename=[pname fname];
    if fname==0 %user pressed Cancel
        info=[];
        return
    end
end

[pname,fname] = fileparts(filename);
if isempty(pname)
    pname='.';
end

if exist([pname,'\', fname,'.hdr'],'file')~=2
    error(['Invalid file name: ' filename]);
end

filename=[pname,'\',fname,'.hdr'];

[fid,message]=fopen(filename);

if fid==-1
    error(message);
end

% //YOKOGAWA ASCII FILE FORMAT
%
% $PublicInfo
fileparse(fid,'$PublicInfo');

% FormatVersion 1.11
info.FormatVersion= sscanf(native2unicode(fileparse(fid,'FormatVersion')), 'FormatVersion %s');      % Version # of header file format

% Model DL1600
info.Model= sscanf(native2unicode(fileparse(fid,'Model')), 'Model %s'); % model name

% Endian Big
info.Endian= sscanf(native2unicode(fileparse(fid,'Endian')), 'Endian %s');    % Endian-Mode (Big/Ltl)

% DataFormat TRACE
info.DataFormat= sscanf(native2unicode(fileparse(fid,'DataFormat')),  'DataFormat %s');  % saving format (Trace/Block) of binary waveform data file
% Trace : grouped in blocks, each block for a single waveform
% Block : grouped in blocks, each block for a certain time interval

% GroupNumber 1
info.GroupNumber= sscanf(native2unicode(fileparse(fid,'GroupNumber')),  'GroupNumber %d'); % number of groups in the file

% TraceTotalNumber 10
info.TraceTotalNumber= sscanf(native2unicode(fileparse(fid,'TraceTotalNumber')),  'TraceTotalNumber %d');   % number of traces in file

% DataOffset 0
info.DataOffset= sscanf(native2unicode(fileparse(fid,'DataOffset')),  'DataOffset %d');  % Start position in binary waveform file

for GroupNum=1:info.GroupNumber

    %
    % $Group#
    fileparse(fid,['$Group' num2str(GroupNum)]);

    % TraceNumber 4
    info.(['Group' num2str(GroupNum)]).TraceNumber=sscanf(native2unicode(fileparse(fid,'TraceNumber')), 'TraceNumber %d');

    % BlockNumber 1
    info.(['Group' num2str(GroupNum)]).BlockNumber= sscanf(native2unicode(fileparse(fid,'BlockNumber')),  'BlockNumber %d');
    traceNumMax = info.(['Group' num2str(GroupNum)]).TraceNumber; %maximum trace number in this group

    % TraceName CH1 CH2 CH3 CH4
    tok='TraceName';
    dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax); %#ok
    for traceNum=1:traceNumMax
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy{1}(traceNum);
    end

    % BlockSize 1002 1002 1002 1002
    tok='BlockSize';
    dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%d', traceNumMax); %#ok
    for traceNum=1:traceNumMax
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy{1}(traceNum);
    end

    % VResolution 1.5625000E+00 1.5625000E+00 1.5625000E+00 1.5625000E+00
    tok='VResolution';
    dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%f', traceNumMax); %#ok
    for traceNum=1:traceNumMax
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy{1}(traceNum);
    end

    % VOffset 0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00
    tok='VOffset';
    dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%f', traceNumMax); %#ok
    for traceNum=1:traceNumMax
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy{1}(traceNum);
    end

    % VDataType IS2 IS2 IS2 IS2
    % ISn : n-byte signed integer
    % IUn : n-byte unsigned integer
    % FSn : n-byte signed real
    % FUn : n-byte unsigned real
    % Bm : m-byte logical data
    tok='VDataType';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, rowpos]=regexp(tokline,'[IFB]\w*\d\s','match','start'); %get entries and their positions. These positions are the starting positions for all lines of the different channels (the VDataType entry should exist for all channels, regardless of their type (i.e., even for logic data)
    if isempty(dummy)
        fclose(fid);
        error([fname '.hdr contains an unknown binary format or VDataType entry not found.']);
    end
    for traceNum=1:traceNumMax          %for all traces in the group
        switch dummy{traceNum}(1)    %contains entry 'IS2' or similar for the selected trace
            case {'I'}                  %integer numbers
                bytenum=sscanf(dummy{traceNum}(3:end), '%d'); %#ok
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.ByteNum=bytenum;
                switch dummy{traceNum}(2)
                    case {'S'}
                        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.BinFormat= ['*int' num2str(8*bytenum)];
                    case {'U'}
                        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.BinFormat= ['*uint' num2str(8*bytenum)];
                end
            case {'F'}                  %real numbers
                bytenum=sscanf(dummy{traceNum}(3:end), '%d'); %#ok
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.ByteNum=bytenum;
                switch dummy{traceNum}(2)
                    case {'S'}
                        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.BinFormat= ['*float' num2str(8*bytenum)];
                    case {'U'}
                        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.BinFormat= ['*float' num2str(8*bytenum)];
                        warning('HDRread:TypeConversion','Type conversion from unsigned float to float may result in errors. Please report this to the author of wvfread.m');
                end
            case {'B'}
                bytenum=sscanf(dummy{traceNum}(2:end), '%d');
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.ByteNum=bytenum;
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).VDataType.BinFormat= ['*uint' num2str(bytenum)];
                if not((bytenum==16) & strcmp(info.Model,'DL750'))
                    fclose(fid);
                    error('Raw binary format (used for logic signals) is so far only implemented for DL750 scopes in wvfread.m');
                end
        end
    end

    % VUnit V V A V
    tok='VUnit';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %#ok %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy(posind);
        end
    end

    % VPlusOverData 32767 32767 32767 32767
    tok='VPlusOverData';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            content=sscanf(dummy{1,posind},'%d');
            if isnumeric(content)
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= content;
            end
        end
    end

    % VMinusOverData -32769 -32769 -32769 -32769
    tok='VMinusOverData';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            content=sscanf(dummy{1,posind},'%d');
            if isnumeric(content)
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= content;
            end
        end
    end

    % VIllegalData 32767 32767 32767 32767
    tok='VIllegalData';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            content=sscanf(dummy{1,posind},'%d');
            if isnumeric(content)
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= content;
            end
        end
    end

    % VMaxData 32766 32766 32766 32766
    tok='VMaxData';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            content=sscanf(dummy{1,posind},'%d');
            if isnumeric(content)
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= content;
            end
        end
    end

    % VMinData -32768 -32768 -32768 -32768
    tok='VMinData';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            content=sscanf(dummy{1,posind},'%d');
            if isnumeric(content)
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= content;
            end
        end
    end

    % HResolution 5.0000000E-09 5.0000000E-09 5.0000000E-09 5.0000000E-09
    tok='HResolution';
    dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%f', traceNumMax); %#ok
    for traceNum=1:traceNumMax
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy{1}(traceNum);
    end

    % HOffset -2.5000000E-06 -2.5000000E-06 -2.5000000E-06 -2.5000000E-06
    tok='HOffset';
    dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%f', traceNumMax); %#ok
    for traceNum=1:traceNumMax
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy{1}(traceNum);
    end

    % HUnit s s s s
    tok='HUnit';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    [dummy, pos]=regexp(tokline,'\S\w*','match','start'); %#ok %get entries and their positions.
    %dummy=textscan(stringparse(native2unicode(fileparse(fid,tok)), tok), '%s', traceNumMax);
    for traceNum=1:traceNumMax
        posind=find(pos==rowpos(traceNum)); %Zuordnung gefundene Position zu Spalte (aus rowpos)
        if ~isempty(posind) % if entry exists for the current trace
            info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok)= dummy(posind);
        end
    end


    BlockNum=info.(['Group' num2str(GroupNum)]).BlockNumber; % number of blocks in this group

    % Date 2001/07/25 2001/07/25 2001/07/25 2001/07/25
    tok='Date';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    dummy=regexp(tokline,['^' tok '(\d*)[ \t]+'],'tokens','once'); %#ok %get entries and their positions.
    if (isempty(dummy{1})) %if only one block in the file (only Date entry, no further DateN entries
        tok1 = tok;   % Date token is just "Date", either because only 1 block or because in newer DL750 (at least from model version 6.22 on) only a single Date entry may be written in the header file.
    else
        tok1 = [tok '1']; % Date token of first block is "Date1"
    end
    dummy=textscan(stringparse(tokline, tok1), '%s', traceNumMax); %#ok %search the date entries of the first block (tok1= either Date or Date1)
    for traceNum=1:traceNumMax % and save them to the entries in the info structure.
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok1)= dummy{1}(traceNum);
    end
    if ~strcmp(tok,'Date') %if as many date entries as blocks exist,
        for dateNum=2:BlockNum % do the same for them:
            dummy=textscan(stringparse(native2unicode(fileparse(fid,[tok num2str(dateNum)])), [tok num2str(dateNum)]), '%s', traceNumMax); %#ok search the date entries of the n-th block
            for traceNum=1:traceNumMax %and save them to the entries in the info structure.
                info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).([tok num2str(dateNum)])= dummy{1}(traceNum);
            end
        end
    end

    % Time 01:45:00 01:45:00 01:45:00 01:45:00
    tok='Time';
    tokline=native2unicode(fileparse(fid,tok));                 % get content of line containing the entry tok
    dummy=regexp(tokline,['^' tok '(\d*)[ \t]+'],'tokens','once'); %#ok %get entries and their positions.
    if (isempty(dummy{1})) %if only one block in the file (only Time entry, no further TimeN entries
        tok1 = tok;   % Time token of first block is just "Time", because only 1 block
    else
        tok1 = [tok '1']; % Time token of first block is "Time1"
    end
    dummy=textscan(stringparse(tokline, tok1), '%s', traceNumMax); %#ok %search the date entries of the first block (tok1= either Time or Time1)
    for traceNum=1:traceNumMax % and save them to the entries in the info structure.
        info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).(tok1)= dummy{1}(traceNum);
    end
    for timeNum=2:BlockNum % if more than 1 blocks / Time entries exist, do the same for them:
        dummy=textscan(stringparse(native2unicode(fileparse(fid,[tok num2str(timeNum)])), [tok num2str(timeNum)]), '%s', traceNumMax); %#ok %search the Time entries of the n-th block
        for traceNum=1:traceNumMax %and save them to the entries in the info structure.
            info.(['Group' num2str(GroupNum)]).(['Trace' num2str(traceNum)]).([tok num2str(timeNum)])= dummy{1}(traceNum);
        end
    end
end
fclose(fid);

function tline=fileparse(fid, tok)
% reads remaining lines from file fid (already opened for reading) until tok string
% is found in the line or EOF is reached. If tok is found, the line
% containing tok is returned. If EOF is reached, an error message is
% produced.

tline='';
while isempty(regexp(tline,['^\' tok '\w*'],'once')) % read lines from file until entry is found or EOF is reached
    tline= fgetl(fid);
    if tline==-1  % EOF reached
        fclose(fid);
        error(['Could not find "' tok '" entry in the file'])
    end
end

function [remainder, pos]=stringparse(str, tok)
%find the next occurrence of a token string (tok) within a longer string (str) and
%output the remaining string without the token

pos=min(strfind(str,tok));
remainder= str(pos+length(tok):end);