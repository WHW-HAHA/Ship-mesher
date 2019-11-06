function patran_flip_normals(varargin)

% ========================================================================
% SYNTAX:
% patran_flip_normals(varargin)
%
% Example:
% patran_mirror(pat)
%
% Description:
% This script is used by patran_ops to flip the normals of a patran file.
% Not for stand alone usage.
%
% Input:
% None from gui
%
% OR from patran_mirror:
% 1) patran structure
% 2) path
% 3) file
%
% Output:
% Patran file structure with flipped normals -> patran_write is called.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================

if nargin==0
    % Select patranfile
    [file,path,ind]=uigetfile('*.pat','Select patranfile');
    % Check cancel button
    if ind==0;disp('Action cancelled, no bodies given in patran_scale');return;end
    % Read patran file
    pat=patran_read_pat([path file]);cd(path);
else
    pat=varargin{1};
    path=varargin{2};
    file=varargin{3};
end

% Flip normals
pat.pan=pat.pan(:,[2 1 4 3]);
disp(['Normal panel directions ' file path ' flipped.']);

% Write patran file
if nargin==0
    patran_write(pat,path,[file(1:end-4) '_nor.pat']);
else
    patran_write(pat,path,file);
end