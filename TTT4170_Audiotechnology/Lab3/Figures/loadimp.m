function [h, Fs, Format, Comment] = loadimp(filename)
%
%LOADIMP  Load WinMLS measurement file.
%   Supported file extensions:
%    .WMB binary file containing impulse response and header
%	  .WMT ASCII (text) file containing impulse response and header
%    (Please use wavread.m to read wav-files)
%
%	  If the file can not be read, an error message is given and the
%	  variables are returned as empty, that is <[]>;
%
%   [h]=LOADIMP(filename) loads a .WMB or .WMT format file specified by "filename", 
%       returning the sampled data in variable "h". The extension 
%       must be included in the filename. 
%
%   [h,Fs]=LOADIMP(filename) loads a .WMB or .WMT format file specified by  
%       "filename", returning the sampled data in variable "h" and the 
%       sample rate in variable "Fs". 
%   
%   [h,Fs,Format]=LOADIMP(filename) loads a .WMB or .WMT format file specified by  
%       "filename", returning the sampled data in variable "h", the 
%       sample rate in variable "Fs" and format information in the variable 
%       "Format". The format information is returned as a 6 element
%       vector with the following order:
%
%               Format(1)       Sequence Order
%               Format(2)       Number of averages
%               Format(3)       Channel
%               Format(4)       Maximum level recorded during the measurement (percent)
%               Format(5)       # of bits used for record
%               Format(6)       # of bits used for play
%
%       For version 3 files, more fields are present in the header. To change the 
%       arguments returned in the Format variable, selecting additional fields, 
%       change the line assigning to the variable in 'loadimp.m'. The line is marked (*).
%
%	  If the file cannot be read, an error message is given and h is returned as -1;
% 	  This function uses the file function Valstr.m
%
% Copyright Morset Sound Development 1998-2007.
% update 21112001 - removed reading of extra header
% update 13062003 - support for file format version 4, problems
% update 23062003 - fixed problems with support for ver. 4
% Contact: www.winmls.com
%
% Example: [h1, Fs, Format, Comment]=loadimp('h1.wmb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KP08012001: Updated for version 3 headers
% LM03072001: Made it work for Matlab ver. 4

VersionNumber=4;  % Reads header version 3
VersionNumber_1=1;	% Old files are version number 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize return variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = []; Fs = NaN; Format = []; Comment = [];		% Fs is scalar and cannot be []


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(filename, 'r');
if isunix == 1,	
	fid=fopen(filename, 'r','ieee-le');
else
	fid=fopen(filename, 'r');
end

if ( fid == -1 ) 
  disp(' ');  disp('Error. Could not open a file with the the specified filename.')
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND EXTENSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(filename);
if l<3 % check that length of string is larger than three.
  disp(' ');  disp('Error. Wrong extension')
  disp('Either .WMB or .WMT must be used as file extensions,')
  fclose(fid);
  return
end 

ext=upper(filename(l-2:l));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ BINARY DATA FROM FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ext,'WMB')
  % Read format identifier
  id=setstr(fread(fid, [1,4],'uchar'));
  if (id~='WMLS') 
    disp(' '); disp('Error. The file is not a WinMLS binary file');
    fclose(fid);
    return;
  end;

  % Read header data, check on header version first.
  Version=fread(fid,1, 'ulong');
  %if (Version ~= VersionNumber_1 & Version ~= VersionNumber)
  if (Version<1 & Version>VersionNumber)
    	disp(' '); disp('Error. Illegal version of file.');
    	fclose(fid);
    	return;
  end;  
  if (Version == VersionNumber_1) 
   % Read version 1 header
  	AvgNo=fread(fid,1,'ulong');
  	SeqOrder=fread(fid,1,'ulong');
    Length = 2^SeqOrder-1;
  	Fs=fread(fid,1,'ulong');
  	RecBitsPerSample=fread(fid,1,'ulong');
  	PlayBitsPerSample=fread(fid,1,'ulong');
 	MaxRecLevel=fread(fid,1,'ulong');
 	Channel=fread(fid,1,'ulong');
 	Format=[SeqOrder AvgNo Channel MaxRecLevel RecBitsPerSample PlayBitsPerSample];
    
   % read comment
	len=fread(fid,1,'ulong');
	Comment=(setstr(fread(fid, [1,len+1],'uchar')))';
  elseif (Version >= 3) 
   % Read version 3 header
         Title = setstr(fread(fid,80,'uchar'));			% Measurement title
         Comment = (setstr(fread(fid,60,'char')))';		% Channel comment
         DateOfMeas = setstr(fread(fid,20,'uchar'));	% Date of measurement
         d_new_header_exists = fread(fid,1,'ulong');	% For internal use
         Channel = fread(fid,1,'ulong');					% Channel number in measurement
         NumberOfChannels = fread(fid,1,'ulong');		% Total number of channels in this measurement
         Fs = fread(fid,1,'float');							% Sampling frequency
         Length = fread(fid,1,'ulong');					% Length of measurement
         UseFeedbackLoop = fread(fid,1,'uchar');		% Flag: Is feedback loop used?
         Mixer = fread(fid,1,'uchar');						% Flag: Is WinMLS mixer used?
         
         % Input settings
         InputDevice = setstr(fread(fid,30,'uchar'));			% Name of input device
         MaxRecLevel = fread(fid,1,'float');						% Maximum recording level
         RecBitsPerSample = fread(fid,1,'ulong');				% Number of bits per sample
         MixerInputVolume = fread(fid,1,'ulong');				% Mixer input volume
         MixerInputIsCalibrated = fread(fid,1,'uchar');		% Flag: Set to 1 if mixer input is calibrated
         dummy = fread(fid,3,'uchar'); % SKIP 3 bytes of padding
         MixerInputVol_db = fread(fid,1,'float');				% Mixer input volume in decibels.
         HardwareInputIsCalibrated = fread(fid,1,'uchar');	% Flag: Set to 1 if hardware input is calibrated
         InputUnitLabel = setstr(fread(fid,11,'uchar'));		% Input Unit Label
         InputConversionFactordB = fread(fid,1,'float');		% Input Conversion Factor
         
         % Output settings
         OutputDevice = setstr(fread(fid,30,'uchar'));		% Name of output device
         dummy = fread(fid,2,'uchar');	% SKIP 2 bytes of padding
         PlayBitsPerSample = fread(fid,1,'ulong');				% Number of bits per sample
         MixerOutputMasterVolume = fread(fid,1,'ulong');		% Mixer output master volume
         MixerOutputVolume = fread(fid,1,'ulong');				% Mixer output volume
         MixerOutputIsCalibrated = fread(fid,1,'ulong');		% Flag: Set to 1 if mixer output is calibrated
         MixerOutputVol_db = fread(fid,1,'float');				% Mixer input volume in decibels.
         HardwareOutputIsCalibrated = fread(fid,1,'uchar');	% Flag: Set to 1 if hardware output is calibrated
         OutputUnitLabel = setstr(fread(fid,11,'uchar'));	% Output Unit Label
         OutputConversionFactordB = fread(fid,1,'float');	% Output Conversion Factor
         
         SpeedOfSound = fread(fid,1,'float');					% Speed of sound in meters per second
         MeasurementMode = setstr(fread(fid,257,'uchar'));	% Measurement Mode 
         MeasurementSystemCorrection = fread(fid,1,'uchar');% Flag: Is measurement system correction used?
         PreEmphasis = fread(fid,1,'uchar');						% Flag: Is preemphasis used?        
         DeEmphasis = fread(fid,1,'uchar');						% Flag: Is deemphasis used?
         EmphasisFileName = setstr(fread(fid,256,'uchar'));	% Name of emphasis file
         SeqOrder = fread(fid,1,'ulong');							% Sequence order
         AvgNo = fread(fid,1,'ulong');								% Number of averages
         NumberOfPreMLS = fread(fid,1,'ulong');					% Number of Pre-MLS
         TypeOfMLS = setstr(fread(fid,11,'uchar'));			% Type of MLS
         
 		   if (Version == 3) 
             dummy = fread(fid,1,'uchar');		% Skip one byte of padding
         end;
          
         % (*) Return selected information. Change this line to alter the informating returned in the format array.
         Format=[SeqOrder AvgNo Channel MaxRecLevel RecBitsPerSample PlayBitsPerSample];
          
         % read trailing informating, skip extra header
         % len_e=fread(fid,1,'ulong');								% Length of extra header
         % extra_header = fread(fid,len_e,'uchar');				% Skip extra header 
  end;
  
  % Read version 4 data
  if (Version >= 4) 
     dummy = fread(fid,6,'uchar');						% Skip one byte of padding
     MeasID = setstr(fread(fid,16,'uchar')); 					% MeasID
     ExtraHeaderSize = fread(fid,1,'ulong'); 					% not used
     IsTransferfunction = fread(fid,1,'uchar');						% Flag: Is this a transfer function or a recording (scope mode)        
     SpecifyImpLen = fread(fid,1,'uchar');						% Flag: If transfer function: is impulse response truncated. If scope mode: is mode transient or stationary        
     Resample = fread(fid,1,'uchar');						% Flag: Is the data resampled (used for low frequency measurements        
     SweepLenInSec = fread(fid,1,'double');						% Sweep length in seconds
     StartFreq = fread(fid,1,'double');						% Sweep start frequency
     EndFreq = fread(fid,1,'double');						% Sweep end frequency
     HPfiltered = fread(fid,1,'uchar');						% Flag: HighPass-filter used to avoid noise blowup for sine-sweep measurement
     if (Version == 4) 
        dummy = fread(fid,7,'uchar');		%7 Skip 7 bytes of padding
     end;
  end;
  
  % read impulse response
  h=zeros(1,Length); % initialize data to make it faster
  h=fread(fid,Length,'float');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ ASCII DATA FROM FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(ext,'WMT')

  % Read first line to check if it is a WinMLS Ascii file
  string=fgetl(fid);
  if strcmp(string, '#WinMLS datafile')~=1
    disp(' '); disp('Error. The specified file either corrupt or not a WinMLS Ascii file')
    fclose(fid);
    return
  end
  
  % Ignore all lines starting with '#'
  Version=str2num(valstr(fid));
  if (Version>VersionNumber) 
    disp('Illegal version of file');
    fclose(fid);
    return;
  end;

  AvgNo=str2num(valstr(fid));
  SeqOrder=str2num(valstr(fid));
  Fs=str2num(valstr(fid));
  RecBitsPerSample=str2num(valstr(fid));
  PlayBitsPerSample=str2num(valstr(fid));
  MaxRecLevel=str2num(valstr(fid));
  Channel=str2num(valstr(fid));

  % read comment
  Comment=valstr(fid);
  Comment=Comment(2:length(Comment)-1);

  Format=[SeqOrder AvgNo Channel MaxRecLevel RecBitsPerSample PlayBitsPerSample];

  %skip all lines beginning with '#'
  string=fgetl(fid);
  while findstr(string, '#')==1,
    string=fgetl(fid);
  end;

  % read impulse response
  %imp_l=2^SeqOrder-1;
  %h=zeros(1,imp_l); % initialize data
  %h(1:imp_l)=fscanf(fid,'%f',inf);
  h=fscanf(fid,'%f',inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE FILE EXTENSION IS NOT SUPPORTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  disp(' '); disp('Error. The file extension you specified is not supported.')
  disp('Either .WMB or .WMT must be used as file extensions,')
  disp('Type <help loadimp> for more information')
end

fclose(fid);
