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
R=10; h=53.5; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle

N=[110.55251051 0 -110 136.58888333 -60 -60; 0 190.52558883 0 113.49806587 227 0; 0 0 0 515 515 515];
 


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
%[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
% Plot the structure with radius
%R3Ddata.Bradius=interp1([max(radius)],[9.5,10.0],r_b);
%R3Ddata.Sradius=interp1([max(radius)],[1.0,1.6],r_s);
%R3Ddata.Nradius=ones(nn,1);
%tenseg_plot(N,C_b,C_s,[],[],[],'Single Layer Prism',R3Ddata);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=1:2;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg,10);

%% external force, forced motion of nodes, shrink of strings
ind_w=4*3-2;w=9*1000;
ind_dnb=3*(7:9)'; dnb0=5*ones(3,1);
ind_dl0=1; dl0=-0.3;
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);
% calculate external force and
% ind_w=[];w=[];
% ind_dnb=1*3-2; dnb0=0.5;
% ind_dl0=[]; dl0=[];
% ind_w=[3*[4:6]'-2];w=[500*ones(3,1)];
%ind_dnb=4*3-2; dnb0=0.5;

%% input file of ANSYS
ansys_input_gp_00(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),ind_w,w,ind_dnb,dnb0,fullfile(savePath,'Test_tower_ansys'));
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% rest length of members
data.substep=substep;    % substep
% nonlinear analysis
data_out=static_solver2(data);        %solve equilibrium using mNewton method

% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method
t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;
%% plot member force
tenseg_plot_result(1:substep,t_t([1;7;13],:),{'bar','horizontal string','vertical string'},{'Load step','Force (N)'},fullfile(savePath,'plot_member_force.png'),saveimg);
%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'plot_coordinate.png'),saveimg);
%% Plot final configuration
tenseg_plot_catenary(reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],R3Ddata,l0_t(index_s,end));

%% save output data
if savedata==1
    save(fullfile(savePath,['tower_static_',material{1},'.mat']));
end
%% make video of the dynamic
name=fullfile(savePath,['tower_',material{1},'_slack_',num2str(material{2})]); %File name for saving
tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,[],[],[],min(substep,50),name,savevideo,material{2});
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);
%% linearized dynaimcs
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
