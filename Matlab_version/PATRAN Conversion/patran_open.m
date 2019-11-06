function patran_open(varargin)

% ========================================================================
% SYNTAX:
% patran_open
%
% Example:
% patran_open
%
% Description:
% This script is used by matpat to open and show patran structures. Not for 
% stand alone usage.
%
% Input:
% None
%
% Output:
% Patran figure.
%
% Revisions
% 1.0   :   K.Hoefakker, April 2011, part of rewriting MATPAT
%
%=========================================================================

global h32 h33 h34 file path

if nargin == 0
    % Select patranfile
    [file,path,ind]=uigetfile('*.pat','Select patranfile','C:\Users\Torben\Documents\MATLAB\RHDHV\lijnenplannen\');
    % Check cancel button
    if ind==0;disp('Action cancelled, no bodies given in patran_open');return;end
    % Read patran file
    pat=patran_read_pat([path file]);cd(path);
else
    pat=varargin{1};
end

assignin('base','pat',pat)


% Show patran file
% Initialize axis
subplot('position',[0.09 0.25 0.71 0.71]);cla;reset(gca);hold on;
colors=[[64 104 225]/255 0 0 1;...
     [64 104 225]/255 0 0 0;...
    0.3 0.3 0.3 0 0 0];

% Plot bodies
iend=0;
for i=1:pat.nbody
    % Determine panels belonging to body i
    istart=iend+1;
    iend=sum(pat.bpan(1:i));
    % Draw body, color switching red-black-blue
    patch('Vertices',pat.crd,...
        'Faces',pat.pan(istart:iend,:),...
        'FaceColor',colors(mod(i,3)+1,1:3),...
        'EdgeColor',colors(mod(i,3)+1,4:6));
    % Get coordinates for text 'Body i'
    txy=[mean(pat.crd(pat.pan(istart:iend,:),1)),...
        mean(pat.crd(pat.pan(istart:iend,:),2))];
    % Plot text 'Body i'
    text(txy(1),txy(2),2.0,['Body ' num2str(i)],'HorizontalAlignment','center');
    hold on;
end

% Check if called by patran_show_normals
if nargin == 2
    norm=varargin{2};
    for i=1:size(norm,1)
        plot3([norm(i,4) norm(i,4)+norm(i,1)],...
            [norm(i,5) norm(i,5)+norm(i,2)],...
            [norm(i,6) norm(i,6)+norm(i,3)],...
            '-k');
    end
end

% Figure settings
xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]');
view(-45,30);
axis equal;
rotate3d on;
grid on;
% Enable figure storing options
% set(h32,'Enable','on');
% set(h33,'Enable','on');
% set(h34,'Enable','on');

