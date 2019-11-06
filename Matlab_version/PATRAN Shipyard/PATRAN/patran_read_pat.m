function pat=patran_read_pat(file)

% ========================================================================
% SYNTAX:
% patran_read_pat(file.pat)
%
% Example:
% patran_read_pat('vessel.pat')
%
% Description:
% This script is used by patran_ops to read a patran file.
% Not for stand alone usage.
%
% Input: 
% Patran file.
%
% Output: 
% Structure with patran file information.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================

% Open patran file
fid=fopen(file,'r');
if fid==-1;disp([path file ' is invalid. Error in patran_read_pat.']);return;end

% Get general information
line=fgetl(fid);line=fgetl(fid);
pat.nbody=str2double(line(13:14));
for i=1:pat.nbody
    eval(['pat.bpan(' num2str(i) ')=str2double(line(15+(' num2str(i) '-1)*4:18+(' num2str(i) '-1)*4));'])
end
line=fgetl(fid);
pat.ncrd=str2double(line(27:34));
pat.npan=str2double(line(35:42));
pat.info=fgetl(fid);

% Get coordinates
for i=1:pat.ncrd,
    line=fgetl(fid);line=fgetl(fid);
    pat.crd(i,1:3)=strread(line,'%16f',3);
    line=fgetl(fid);
end

% Get panel nodes
for i=1:pat.npan,
    line=fgetl(fid);line=fgetl(fid);line=fgetl(fid);
    pat.pan(i,1:4)=strread(line,'%8f',4);
end

fclose(fid);

























