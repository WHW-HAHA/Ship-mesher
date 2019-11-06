function pat=patran_offset(file,path,offset)

% ========================================================================
% SYNTAX:
% patran_relocate
%
% Example:
% patran_relocate
%
% Description:
% This script is used by patran_ops to give a patran file an offset.
% Not for stand alone usage.
%
% Input:
% 1) Write option: True from gui, false from function calls.
%
% 2) Patran structure (from function calls)
%
% OR:
% 2) Patran file with varargin: 'path\file' (or by gui: Optional)
% 3) Offset values from inifile: [x y z] (or by gui: Optional)
%
% Output:
% Relocated patran file.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================

% Check if used directly or by inifile structure
if true %Direct use

    pat=patran_read_pat([path file]);cd(path);
elseif nargin==2
    pat=varargin{1};
    prompt={'Offset x[m]:','Offset y[m]:','Offset z[m]:'};
    offset=str2double(inputdlg(prompt,'Offset values',1,{'0.0','0.0','0.0'}));
    if isempty(offset);disp('Action cancelled, no offset given in patran_offset');return;end
else
    pat=varargin{1};
    offset=varargin{2};
end

% Relocate patran file
for i=1:3;pat.crd(:,i)=pat.crd(:,i)+offset(i);end
write=1;
if write==1
    filename = evalin('base','plan.filepath');
    [pathstr, name] = fileparts(filename); 
    pat.crd = round(pat.crd,8); %avoid floating point errors
    patran_write(pat,path,[name '.pat']);
    cd ..
end

return













