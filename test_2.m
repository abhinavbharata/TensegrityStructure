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

%-----------------------------------------------------------------------------------------------------
%% Group information
%generate group index
gr={(7:9)};     % number of elements in one group
index_gp=[1,2]; 
Gp=tenseg_str_gp(gr,C);    %generate group matrix
I=eye(size(Gp,2));
e_d=I(:,index_gp);        % e_d is the matrix to select group of member with designed force
l_d=e_d'*l_gp;            % length of top center circular strings
qd=fd./l_d;
z=(e_d'*V2)\(qd-e_d'*pinv(A_1ag)*w0a);   %self-stress coefficient

q1_gp=pinv(A_1ag)*w0a;
q1=Gp*q1_gp;
q2_gp=V2*z;
q2=Gp*q2_gp;
q_gp=q1_gp+q2_gp;               % force density in group
q=q1+q2;                        % force density
t=diag(l)*q;                    % force vector
ne=numel(t);        %number of elements

% cross sectional of strings
A_s=t(index_s)/sigmas/c_s;	% area of string
r_s=sqrt(A_s/pi);               % radius of string

% cross sectional of all members
I3=eye(ne);
Ind_b=I3(:,index_b);            % index matrix for bar
Ind_s=I3(:,index_s);            % index matrix for string
A=[Ind_b,Ind_s]*[A_b;A_s];      % cross sectional area
mass=rho.*A.*l0; 

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

