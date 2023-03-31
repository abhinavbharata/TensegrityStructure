

%% N C of the structure
% Manually specify node positions of a tensegrity tower.
R=10; h=53.5; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle

%N=[110.55251051 0 -110 136.58888333 -60 -60; 0 190.52558883 0 113.49806587 227 0; 0 0 0 515 515 515];

  for i=1               % nodal coordinate matrix N
     N(:,i)=[0,0,0];
 end
 for i=3               % nodal coordinate matrix N
     N(:,i)=R*[12;22.1;0];
 end
 for i=2               % nodal coordinate matrix N
     N(:,i)=R*[(cos(2*pi*(i-1)/4.8)*100)-0.6819,sin(2*pi*(i-2)/p),0];
 end
 for i=4               % nodal coordinate matrix N
     N(:,i)=R*[-2.76;7.8;50];
 end
for i=5               % nodal coordinate matrix N
     N(:,i)=R*[18.7;-5.35;50];
 end
 for i=6               % nodal coordinate matrix N
     N(:,i)=R*[18;20;50];
 end


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
grid on
