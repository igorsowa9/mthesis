function varargout = wvfread(varargin)
% function [y,t, ind_over, ind_under, ind_illegal, info] = wvfread(filename, Group, Trace, datapoints, step, startind)
%
% loads YT-waveform data from *.wvf file saved by Yokogawa DL 750, DL 1600,
% DL716 series and DL 7440 Oscilloscopes into the variables y (y data) and t
%(time data).  The matrices ind_over and ind_under, illegal
% contain the indices of overranged data points outside the upper / lower
% limit of the TDS AD converter and indices of illegal values.
%
% optional input arguments:
% Group, Trace: specifies which trace should be read (omitted values
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

outparams='y';
if nargout>1
    outparams=[outparams ',t'];
    if nargout>2
        outparams=[outparams ',ind_over'];
        if nargout>3
            outparams=[outparams ',ind_under'];
            if nargout>4
                outparams=[outparams ',ind_illegal'];
                if nargout>5
                    outparams=[outparams ',info'];
                end
            end
        end
    end
end

if nargin<=3
    eval(['[' outparams ']= wvfreadb(varargin{1:nargin});']);
else
    eval(['[' outparams ']= wvfreadb(varargin{1:3},1,varargin{4:nargin});']);
end
varargout{1}=y;
if nargout>1
    varargout{2}=t;
    if nargout>2
        varargout{3}=ind_over;
        if nargout>3
            varargout{4}=ind_under;
            if nargout>4
                varargout{5}=ind_illegal;
                if nargout>5
                    varargout{6}=info;
                end
            end
        end
    end
end



