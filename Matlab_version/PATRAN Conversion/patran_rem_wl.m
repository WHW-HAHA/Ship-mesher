function patran_rem_wl

% ========================================================================
% SYNTAX:
% patran_remove_waterline
%
% Example:
% patran_remove_waterline
%
% Description:
% This script is used by patran_ops to remove the waterline panels from a
% patran file.
% Not for stand alone usage.
%
% Input: 
% None.
%
% Output: 
% Structure with patran file information -> patran_write is called.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================

% Select patranfile
[file,path,ind]=uigetfile('*.pat','Select patranfile');
% Check cancel button
if ind==0;disp('Action cancelled, no bodies given in patran_rem_wl');return;end
% Read patran file
pat=patran_read_pat([path file]);cd(path);

% Find waterline panels
nwl=0;b=1;wl(1)=0;
for i=1:pat.npan
    % Check for waterline panels
    if sum(abs(pat.crd(pat.pan(i,:),3)))==0
        nwl=nwl+1;
        panwl(nwl)=i;
    end
    % Check number of bodypanels left
    if i==sum(pat.bpan(1:b))
        if nwl~=0
            wl(b+1)=length(panwl)-sum(wl(1:b));
            bpan(b)=pat.bpan(b)-wl(b+1);
        else
            disp(['No waterline elementes found in body ' num2str(b) ' in ' path file '.']);
        end
        b=b+1;
    end
end
if nwl==0;disp(['No waterline elementes found in ' path file '.']);return;end
pat.pan(panwl,:)=[];
pat.bpan=bpan;
pat.npan=size(pat.pan,1);

disp(['Waterline panels removed from ' path file '.']);
patran_rem_unusedcrd(pat,path,[file(1:end-4) '_remwl.pat']);
