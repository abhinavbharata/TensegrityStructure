%% N C of the structure
% Manually specify node positions of a tensegrity tower.
R=10; h=30; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle
 
for i=1:p               % nodal coordinate matrix N
     N(:,i)=[sec(2*pi*(i-1)/0.6382),sin(2*pi*(i-1)/p),0];
 end

