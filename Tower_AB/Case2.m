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
%tenseg_plot(N,C_b,C_s);
%title('Single Layer Prism');
%% Boundary constraints
pinned_X=(1:3)'; pinned_Y=(1:3)'; pinned_Z=(1:3)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%% Group information
%generate group index
gr={(7:9)};     % number of elements in one group
index_gp=[1,2]; 
Gp=tenseg_str_gp(gr,C);    %generate group matrix
fd=-1e5*ones(2,1);              % force in bar is given as -1000

H=N*C';                     % element's direction matrix
l=sqrt(diag(H'*H));         % elements' length
l_gp=pinv(Gp)*l;            % elements' length in group

I=eye(size(Gp,2));
e_d=I(:,index_gp);        % e_d is the matrix to select group of member with designed force
l_d=e_d'*l_gp;            % length of top center circular strings
qd=fd./l_d;



Cell_H=mat2cell(H,3,ones(1,size(H,2)));          % transfer matrix H into a cell: Cell_H

A_1a=Ia'*kron(C',eye(3))*blkdiag(Cell_H{:});     % equilibrium matrix
A_1ag=A_1a*Gp;                                   % equilibrium matrix in group constraints

[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);
%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

z=(e_d'*V2)\(qd-e_d'*pinv(A_1ag)*w0a);   %self-stress coefficient

q1_gp=pinv(A_1ag)*w0a;
q1=Gp*q1_gp;
q2_gp=V2*z;
q2=Gp*q2_gp;
q_gp=q1_gp+q2_gp;               % force density in group
q=q1+q2;                        % force density
t=diag(l)*q;                    % force vector
ne=numel(t);        %number of elements

% cross sectional of bars
index_b=find(t<0);              % index of bar in compression
% cross sectional of strings
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,~,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);

%A_s=t(index_s)/sigmas/c_s;	% area of string
%r_s=sqrt(A_s/pi);               % radius of string

% cross sectional of all members
I3=eye(ne);
Ind_b=I3(:,index_b);            % index matrix for bar
Ind_s=I3(:,index_s);            % index matrix for string
A=[Ind_b,Ind_s]*[A_b;A_s];      % cross sectional area
%mass=rho.*A.*l0; 

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=1:2;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg,10);




