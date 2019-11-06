function patran_mirror(file,path)

% ========================================================================
% SYNTAX:
% patran_mirror
%
% Example:
% patran_mirror
%
% Description:
% This script is used by patran_ops to mirror a patran file.
% Not for stand alone usage.
%
% Input:
% None
%
% Output:
% Mirrored patran file structure, patran_flip_normals -> patran_write is
% called.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================
% Mirror plane options
planes={'xz','xy','yz'};
plcol=[2 3 1];
% Read patran file
pat=patran_read_pat([path file]);cd(path);
% Define mirror plane
mirr='xz';
plane=find(strcmp(mirr,planes));

% Mirror patran file
pat.crd(:,plcol(plane))=-pat.crd(:,plcol(plane));
disp([path file ' mirrored']);

% Flip normals and write patran file
patran_flip_normals(pat,path,[file(1:end-4) '_mir.pat']);


