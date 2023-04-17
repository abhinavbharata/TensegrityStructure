

%% N C of the structure
% Manually specify node positions of a tensegrity tower.
R=10; h=53.5; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle

N=[0 250 12 -27.6 185 185;0 0 221.3 78 200 -53.5; 0 0 0 500 500 500];


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
