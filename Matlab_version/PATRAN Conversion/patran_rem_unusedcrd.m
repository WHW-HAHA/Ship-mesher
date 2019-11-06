function patran_rem_unusedcrd(varargin)

% ========================================================================
% SYNTAX:
% patran_rem_unusedcrd(varargin)
%
% Example:
% patran_rem_unusedcrd(pat.pat,path,file)
%
% Description:
% This script is used by patran_ops to remove unused coordinates from a
% patran file.
% Not for stand alone usage.
%
% Input:
% None from gui
%
% OR from patran_rem_wl/patran_wl_zero:
% 1) patran structure
% 2) path
% 3) file
%
% Output: 
% Structure with patran file information -> patran_write is called.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================

% Read patran file
if nargin==0
    % Select patranfile
    [file,path,ind]=uigetfile('*.pat','Select patranfile');
    % Check cancel button
    if ind==0;disp('Action cancelled, no bodies given in patran_rem_unusedcrd');return;end
    % Read patran file
    pat=patran_read_pat([path file]);cd(path);
else
    pat=varargin{1};
    path=varargin{2};
    file=varargin{3};
end

% Find unused coordinates
nuc=0;
for i=1:pat.ncrd
    if isempty(find(pat.pan==i,1))
        nuc=nuc+1;
        uc(nuc)=i;
    end
end

% Delete unused coordinates and order panels
if exist('uc','var')
    uc=sort(uc,'descend');
    for i=1:length(uc)
        pat.pan(find(pat.pan>uc(i)))=pat.pan(find(pat.pan>uc(i)))-1;
    end
    pat.crd(uc,:)=[];
    pat.ncrd=size(pat.crd,1);
end

if nuc==0 && nargin==0
    disp(['No unused coordinates found in ' path file '.']);
    disp('No patran file written.');
else
    if nargin==0;
        disp(['Unused coordinates removed from ' path file '.']);
        file=[file(1:end-4) '_remuc.pat'];
    end
    patran_write(pat,path,file);
end