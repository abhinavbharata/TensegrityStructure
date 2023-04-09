clc;clear;close all; 
% Specify material properties
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties:'linear_elastic'£¬ 'multielastic'£¬ 'plastic'
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)
% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;    % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3
substep=100;            % load steps
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of a tensegrity tower.
R=10; h=500; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle

%for i=1:p               % nodal coordinate matrix N
  %   N(:,i)=R*[cos(2*pi*(i-1)/p),sin(2*pi*(i-1)/p),0];
 %end
%for i=p+1:2*p
 %   N(:,i)=[R*cos(2*pi*(i-1)/p+beta*pi/180),R*sin(2*pi*(i-1)/p+beta*pi/180),h];
%end

N=[0 -216.51 -216.51 125 -125 0; 0 125 -125 216.51 216.51 0; 0 0 0 500 500 500];


% Manually specify connectivity indices.
C_b_in = [1 5;2 6;3 4];   % This is indicating the bar connection
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);
% Manually specify connectivity indices.
C_s_in = [4 5;5 6;6 4;1 4;2 5;3 6];  % This is indicating the string connection
% Convert the above matrices into full connectivity matrices.
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
title('Single Layer Prism');

%% Boundary constraints
pinned_X=(1:3)'; pinned_Y=(1:3)'; pinned_Z=(1:3)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% Group information
%generate group index
gr={(7:9)};     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
%% self-stress design
%Calculate equilibrium matrix and member length 
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);
%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;
%prestress design
index_gp=[1,2];                 % number of groups with designed force
fd=-1e5*ones(2,1);              % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[9.4,10.0],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[1.0,1.6],r_s);
R3Ddata.Nradius=ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Single Layer Prism',R3Ddata);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=1:2;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg,10);