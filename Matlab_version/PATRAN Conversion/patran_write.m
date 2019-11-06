function patran_write(pat,path,file)

% ========================================================================
% SYNTAX:
% patran_write
%
% Example:
% patran_write(pat,path,file)
%
% Description:
% This script is used by patran_ops to write patranfiles. Not for stand
% alone usage.
%
% Input: 
% Patran structure, path and file to write to.
%
% Output: 
% Patran file.
%
% Revisions
% 1.0   :   K.Hoefakker, March 2011, part of rewriting MATPAT
%
%=========================================================================
% whos
% global h0
% figure(h0);

% Patran information
[n d]=weekday(date);
c=clock; 

try
    pat.info = evalin('base','info');
catch
    pat.info=[d ' ' strrep(date,'-',' ') ' ' strrep(num2str(round(c(4:6))),'  ',':')];
    pat.info=strrep(pat.info,'::',':0');
end

% Open patran file
fpat=fopen([path file],'w');
if fpat==-1;disp([path file ' could not be opened in patran_write.']);return;end

% Write body information to patran file
fprintf(fpat,'25       0       0       1       0       0       0       0       0\r\n');
string='ISYM=0 nbod=%2.0f';
for i=1:pat.nbody;string=[string '%4.0f'];end
string=[string '\r\n'];
vector=[pat.nbody pat.bpan];
fprintf(fpat,string,vector);

% Write panel/coordinate information to patran file
fprintf(fpat,'26       0       0       1%8.0f%8.0f       0       0       0\r\n',pat.ncrd,pat.npan);
fprintf(fpat,[pat.info '\r\n']);

% Print all coord data
for i=1:pat.ncrd
    fprintf(fpat,' 1%8.0f       0       2       0       0       0       0       0\r\n',i);
    a=sprintf('%17.9E%17.9E%17.9E',pat.crd(i,1:3));
    a=strrep(a,'E+0','E+');a=strrep(a,'E-0','E-');
    fprintf(fpat,'%s\r\n',a);
    fprintf(fpat,'1G       6       0       0  000000\r\n');
end

% Print all panel data
for i=1:pat.npan
    fprintf(fpat,' 2%8.0f       4       2       0       0       0       0       0\r\n',i-1);
    fprintf(fpat,'       4       0       0       0 0.000000000E+00 0.000000000E+00 0.000000000E+00\r\n');
    fprintf(fpat,'%8.0f%8.0f%8.0f%8.0f\r\n',pat.pan(i,1:4));
end

% Print last line
fprintf(fpat,'99       0       0       1       0       0       0       0       0\r\n');

% Close patran output file
fclose(fpat);
cd(path);
disp(['Patran file ' path file ' written.']);
if strcmp(file(end-6:end),'nor.pat')
    patran_show_normals(pat);
else
    patran_open(pat);
end
