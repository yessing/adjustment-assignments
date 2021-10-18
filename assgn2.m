% adjustment assignment 2
% @author Yang Yexing yessingyang@outlook.com
% 
% @brief
% An example on GPS net-adjustment use
% a) classic method with A4 fixed
% b) rank-deficient free net
% and compare
clear,clc
%% baseline input
% n=3*6=18
b14 = [4528.9347,-1700.3398,4325.6152].';
db14 = formDbb([67600,-90516,-62859,168100,107230,90000]);
b24 = [3912.9908,  319.5758,1241.0104].';
db24 = formDbb([57600,-77617,-52109,144400,89046,72900]);
b34 = [5276.2233, 1818.0175,-163.2124].';
db34 = formDbb([57600,-77595,-54089,144400,92365,78400]);
b21 = [-615.9439, 2019.9151,-3084.6049].';
db21 = formDbb([72900,-98502,-67366,184900,-115488,96100]);
b31 = [747.2886,  3518.3567,-4488.8277].';
db31 = formDbb([67600,-90779,-62899,168100,106849,90000]);
b32 = [1363.2325, 1498.4417,-1404.2228].';
db32 = formDbb([422500,-300539,-278659,532900,332169,291600]);
lmaj0 = [b14;b24;b34;b21;b31;b32];



%% x ini val input
% x1 y1 z1 ... x4 y4 z4
i3 = eye(3);
o3 = zeros(3,3);

x00 = [-14063.675,93526.066,7431.515,...
    -14679.619,95545.981,4346.91,...
    -13316.387,97044.422,2942.687,...
    -18592.569,95226.358,3105.897].';
x40 = [-18592.569,95226.358,3105.897].';
mb0 = [-i3, o3, o3, i3;...
       o3,-i3, o3, i3;...
       o3, o3,-i3, i3;...
       i3,-i3, o3, o3;...
       i3, o3,-i3, o3;...
       o3, i3,-i3, o3];
mb0 = -mb0;
%   
sig0 = 0.001;
qv = blkdiag(db14,db24,db34,db21,db31,db32)./(sig0*sig0);
pv = qv^(-1);
% lmin = Lmaj-mb*x0;


%% a) n=6*3=18,t=3*3=9
x0 = x00(1:9);
mb = mb0(:,1:9);
lmaj = [x40;x40;x40;zeros(9,1)]+lmaj0;
lmin = lmaj - mb*x0;
xhata = (mb.'*pv*mb )^(-1)*mb.'*pv*lmin;

%% b) 
x0 = x00;
mb = mb0;
lmin = lmaj0 - mb*x0;
mn = mb.' * pv * mb;
xhat = pinv(mn) * mb.' * pv * lmin;

function Dbb = formDbb(arr)
    Dbb = [arr(1),arr(2),arr(3);...
           arr(2),arr(4),arr(5);...
           arr(3),arr(5),arr(6)];
    Dbb = 1e-10*Dbb;
end


