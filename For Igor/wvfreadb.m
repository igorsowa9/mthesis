function [y,t, ind_over, ind_under, ind_illegal, info] = wvfreadb(filename, Group, Trace, Block, datapoints, step, startind)
% function [y,t, ind_over, ind_under, ind_illegal, info] = wvfreadb(filename, Group, Trace, Block, datapoints, step, startind)
%
% loads YT-waveform data from *.wvf file saved by Yokogawa DL 750, DL 1600,
% DL716 series and DL 7440 Oscilloscopes into the variables y (y data) and t
%(time data).  The matrices ind_over and ind_under, illegal
% contain the indices of overranged data points outside the upper / lower
% limit of the TDS AD converter and indices of illegal values.
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
    [fname,pname]=uigetfile({'*.wvf', 'Yokogawa Waveform Files (*.wvf)'},'Choose Yokogawa waveform file');
    filename=[pname fname];
    if fname==0 %user pressed Cancel
        y=[];
        t=[];
        ind_over=[];
        ind_under=[];
        ind_illegal=[];
        info=[];
        return
    end
end


[pname,fname] = fileparts(filename);
if isempty(pname)
    pname='.';
end
if exist([pname,'\',fname,'.wvf'],'file')~=2
    error(['Invalid file name: ' pname,'\',fname,'.wvf']);
end
if exist([pname,'\',fname,'.hdr'],'file')~=2
    error(['Invalid file name: ' pname,'\',fname,'.hdr']);
end

info=hdrread([pname,'\',fname,'.hdr']);

filename=[pname,'\',fname,'.wvf'];

[fid,message]=fopen(filename);

if fid==-1
    error(message);
end

if ((nargin<4) || isempty(Block))
    Block=1;
end
if ((nargin<3) || isempty(Trace))
    Trace=1;
end
if ((nargin<2) || isempty(Group))
    Group=1;
end

if nargin<6
    step=1;
end

if nargin<7
    startind=1;
end

if step<1 || (round(step)~=step)
    step=1; 
    warning('WVFread:stepPosint',['"step" input parameter must be a positive integer.\nSetting step=1.']);
end
if startind<1 || (round(startind)~=startind)
    startind=1; 
    warning('WVFread:startindPosint',['"startind" input parameter must be a positive integer.\nSetting startind=1.']);
end

%range checking for input parameters Group and Trace
if Group<1 || (round(Group)~=Group)
    fclose(fid);
    [pname,fname,ext] = fileparts(filename); %#ok
    error('WVFread:GroupPosInt',['Input parameter "Group" must be a positive integer.']);
end
if info.GroupNumber >= Group
    tracenum= info.(['Group' num2str(Group)]).TraceNumber;
    blocknum= info.(['Group' num2str(Group)]).BlockNumber;
    if tracenum < Trace  %#ok
        fclose(fid);
        [pname,fname,ext] = fileparts(filename); %#ok
        error('WVFread:TraceTooLarge',['Input parameter "Trace" is too large. Maximum number of traces in \ngroup' num2str(Group) ' in ' fname ext ' is ' num2str(tracenum) '.']);
    end
    if tracenum<1 || (round(tracenum)~=tracenum)
        fclose(fid);
        [pname,fname,ext] = fileparts(filename); %#ok
        error('WVFread:TracePosInt',['Input parameter "Trace" must be a positive integer.']);
    end
    if blocknum < Block  %#ok
        fclose(fid);
        [pname,fname,ext] = fileparts(filename); %#ok
        error('WVFread:BlockTooLarge',['Input parameter "Block" is too large. Maximum number of block in \ngroup' num2str(Group) ' in ' fname ext ' is ' num2str(blocknum) '.']);
    end
    if blocknum<1 || (round(blocknum)~=blocknum)
        fclose(fid);
        [pname,fname,ext] = fileparts(filename); %#ok
        error('WVFread:BlockPosInt',['Input parameter "Block" must be a positive integer.']);
    end
else
    fclose(fid);
    [pname,fname,ext] = fileparts(filename); %#ok
    error('WVFread:GroupTooLarge',['Input parameter "Group" is too large. Maximum number of groups \nin ' fname ext ' is ' num2str(info.GroupNumber) '.']);
end

% Endian-Mode of binary data
% Big : big endian, Motorola 68000
% Ltl : little endian, Intel 86
if strcmpi(info.Endian,'Big');
    byteorder='b';
else
    byteorder='l';
end


offset=info.DataOffset;         % Offset at beginning of binary file

for groupind=1:(Group-1)        % first for all groups before the selected group...
    tracenum= info.(['Group' num2str(groupind)]).TraceNumber;
    for traceind=1:tracenum     %...for all traces and blocks in these groups...
        %...add the offsets of all traces...
        offset= offset + info.(['Group' num2str(groupind)]).(['Trace' num2str(traceind)]).BlockSize * info.(['Group' num2str(groupind)]).BlockNumber * info.(['Group' num2str(groupind)]).(['Trace' num2str(traceind)]).VDataType.ByteNum;
    end
end

ByteNum=info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VDataType.ByteNum;

% File Format Trace or Block saves traces and blocks in different
% sequences:
% Trace : Group1,Trace1,Block1, Group1,Trace1,Block2, ..., Group1,Trace1,BlockN, ..., Group1,Trace2,Block1, ..., Group2,Trace1,Block1, ..., GroupN,TraceN,BlockN (each block for a specific waveform)
% Block : Group1,Trace1,Block1, Group1,Trace2,Block1, ..., Group1,TraceN,Block1, ..., Group1,Trace1,Block2, ..., Group2,Trace1,Block1, ..., GroupN,TraceN,BlockN (each block for a specific time interval).
if strcmpi(info.DataFormat,'Trace')
    %trace format:
    for traceind=1:(Trace-1)  %... add offsets of the traces before the selected trace in the selected group
        offset= offset + info.(['Group' num2str(Group)]).(['Trace' num2str(traceind)]).BlockSize * info.(['Group' num2str(Group)]).BlockNumber * info.(['Group' num2str(Group)]).(['Trace' num2str(traceind)]).VDataType.ByteNum;
    end
    %... and finally add offsets of the blocks in the selected trace in the selected group before the selected block
    offset= offset + info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).BlockSize * (Block-1) * ByteNum;
else
    %block format:
    for blockind=1:(Block-1) %... add offsets of the blocks before the selected block in the selected group
        for traceind=1:info.(['Group' num2str(Group)]).TraceNumber;
            offset= offset + info.(['Group' num2str(Group)]).(['Trace' num2str(traceind)]).BlockSize * info.(['Group' num2str(Group)]).(['Trace' num2str(traceind)]).VDataType.ByteNum;
        end
    end
    %... and finally add offsets of the traces in the selected block in the selected group before the selected trace
    for traceind=1:(Trace-1)
        offset= offset + info.(['Group' num2str(Group)]).(['Trace' num2str(traceind)]).BlockSize * info.(['Group' num2str(Group)]).(['Trace' num2str(traceind)]).VDataType.ByteNum;
    end
end

offset=offset+(startind-1)*ByteNum;
fseek(fid, double(offset),'bof');    % go to start byte

nop_all= double(info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).BlockSize); %number of data points stored in the selected trace in the file
nop=nop_all-startind+1;
if nargin>=5
    if datapoints<1 || (round(datapoints)~=datapoints)
        datapoints=floor(nop/step); % set to maximum number of data points which can be securely read from the file, using startind and step parameters
        warning('WVFread:datapointsPosInt',['"datapoints" input parameter must be a positive integer.\nSetting datapoints= ' num2str(datapoints) '.']);
    end
    nop = floor(nop/step); %maximum number of data points which can be securely read from the file, using startind and step parameters
    if datapoints > nop %if more datapoints are requested than provided by the trace
        warning('WVFread:inconsistent_params',['The requested combination of input parameters \n' ...
            'datapoints, step and startind would require at least ' num2str(datapoints*step+startind) ' data points in \n'...
            'trace ' num2str(Trace) ' of group ' num2str(Group) '. ' ...
            'The actual number of data points in the trace \nis only ' num2str(nop_all) '. ' ...
            'The number of data points returned by wvfread is thus \n' ...
            'only ' num2str(nop) ' instead of ' num2str(datapoints) '.']);
    else
        nop=datapoints; %if the number of data points in the trace is sufficient, return requested number of datapoints
    end
end

BinFormat= info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VDataType.BinFormat;

values=double(fread(fid,nop,BinFormat, ByteNum*(step-1),byteorder));%#ok %read data values
fclose(fid);

y= info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VOffset + info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VResolution * values;    %scale data values to obtain in correct units
t= info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).HOffset + info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).HResolution * (startind+(1:step:(nop*step))'-1);    %H-Auflösung x (Daten-Nr. – 1) + H-Offset %%% MH

% %handling over- and underranged values
if nargout >=3
    if exist(['info.Group' num2str(Group) '.Trace' num2str(Trace) '.VPlusOverData'],'var')
        ind_over= find(values==info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VPlusOverData);%find indices of values that are larger than the AD measurement range (upper limit)
        %print warning if there are wrong values because they are lying outside
        %the AD converter digitization window:
        if length(ind_over) %#ok
            warning('WVFread:overvals',[int2str(length(ind_over)), ' over range value(s) in file ' filename]); %#ok
        end
    else
        warning('WVFread:over_values',['No VPlusOverData entry for trace ' num2str(Trace) ' in group ' num2str(Group) '.\nIt might happen that there are such values,\nbut the ind_over output variable of wvfread.m is empty.\nPlease check manually. This warning can be switched off by Matlab command:\nwarning off WVFread:over_values' ])
        ind_over=[];
    end
end

if nargout >=4
    if exist(['info.Group' num2str(Group) '.Trace' num2str(Trace) '.VMinusOverData'],'var')
        ind_under= find(values==info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VMinusOverData); %find indices of values that are larger than the AD measurement range (lower limit)
        %print warning if there are wrong values because they are lying outside
        %the AD converter digitization window:
        if length(ind_under) %#ok
            warning('WVFread:undervals',[int2str(length(ind_under)), ' under range value(s) in file ' filename]); %#ok
        end
    else
        warning('WVFread:under_values',['No VMinusOverData entry for trace ' num2str(Trace) ' in group ' num2str(Group) '.\nIt might happen that there are such values,\nbut the ind_under output variable of wvfread.m is empty.\nPlease check manually. This warning can be switched off by Matlab command:\nwarning off WVFread:under_values' ])
        ind_under=[];
    end
end

if nargout >=5
    if exist(['info.Group' num2str(Group) '.Trace' num2str(Trace) '.VIllegalData'],'var')
        ind_illegal= find(values==info.(['Group' num2str(Group)]).(['Trace' num2str(Trace)]).VIllegalData); %find indices of values that are larger than the AD measurement range (illegal values)
        %print warning if there are wrong values because they are lying outside
        %the AD converter digitization window:
        if length(ind_illegal) %#ok
            warning('WVFread:illegalvals',[int2str(length(ind_illegal)), ' illegal value(s) in file ' filename]); %#ok
        end
    else
        warning('WVFread:illegal_values',['No VIllegalData entry for trace ' num2str(Trace) ' in group ' num2str(Group) '.\nIt might happen that there are such values,\nbut the ind_illegal output variable of wvfread.m is empty.\nPlease check manually. This warning can be switched off by Matlab command:\nwarning off WVFread:illegal_values' ])
        ind_illegal=[];
    end
end