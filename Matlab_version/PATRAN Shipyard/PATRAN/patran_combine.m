function patran_combine(file1,file2,path)

% ========================================================================
% SYNTAX:
% patran_combine
%
% Example:
% patran_combine
%
% Description:
% This script is used by patran_ops to combine patran structures into one
% file. Not for stand alone usage.
%
% Input:
% From gui, OR from function calls:
% 
% 1) Ini structure
% 2) Multi patran structure (from function calls)
%
% Output: 
% Patran file.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================

% Check input data, call by gui or call by funciton
if 1 % Called by gui, ask for input
    % Nr. of bodies
    ini.GENERAL.NBODY=2;
    
    % Multi or single body
    mb='Single body';
    if strcmp(mb,'Multibody');ini.GENERAL.MULTIBODY=1;
    elseif strcmp(mb,'Single body');ini.GENERAL.MULTIBODY=0;
    else disp('Action cancelled, no options given in patran_combine');return;end
    clear mb;
    
    % Choose and read patran files
        pat=patran_read_pat([path file1]);cd(path);
        multi.pat1 = pat;
        clear pat;
        
        pat=patran_read_pat([path file2]);cd(path);
        multi.pat2 = pat;
        clear pat;
        
    % Specify path and name to save patran file to
    ini.path=path;
    ini.file='combined.pat';
    if ~strcmp(ini.path(end),'\') && ~strcmp(ini.path(end),'/');ini.path(end+1)='\';end
    
else % Called by function
    ini=varargin{1};
    multi=varargin{2};
end

% Combine patran files
if ini.GENERAL.MULTIBODY==1
    % Initialize with first body
    pat.nbody=ini.GENERAL.NBODY;
    pat.bpan=multi.pat1.bpan;
    pat.ncrd=multi.pat1.ncrd;
    pat.npan=multi.pat1.npan;
    pat.crd=multi.pat1.crd;
    pat.pan=multi.pat1.pan;
    % Add other bodies
    for i=2:ini.GENERAL.NBODY
        pat.bpan(i)=eval(['multi.pat' num2str(i) '.bpan;']);
        pat.npan=sum(pat.bpan);
        pat.pan(end+1:pat.npan,:)=eval(['multi.pat' num2str(i) '.pan+pat.ncrd;']);
        pat.ncrd=eval(['multi.pat' num2str(i) '.ncrd;'])+pat.ncrd;
        pat.crd(end+1:pat.ncrd,:)=eval(['multi.pat' num2str(i) '.crd;']);
    end
else
    % Initialize with first body
    pat.nbody=1;
    pat.bpan=multi.pat1.bpan;
    pat.ncrd=multi.pat1.ncrd;
    pat.npan=multi.pat1.npan;
    pat.crd=multi.pat1.crd;
    pat.pan=multi.pat1.pan;
    % Add other bodies
    for i=2:ini.GENERAL.NBODY
        pat.bpan=eval(['multi.pat' num2str(i) '.bpan;'])+pat.bpan;
        pat.npan=pat.bpan;
        pat.pan(end+1:pat.npan,:)=eval(['multi.pat' num2str(i) '.pan+pat.ncrd;']);
        pat.ncrd=eval(['multi.pat' num2str(i) '.ncrd;'])+pat.ncrd;
        pat.crd(end+1:pat.ncrd,:)=eval(['multi.pat' num2str(i) '.crd;']);
    end
end

patran_write(pat,ini.path,'combined.pat');




