%function [list_ib1, list_ib2, nlist_ib1, nlist_ib2]= geom_check
% clear
clc
close all

%----------------------------------------
% Input values
%----------------------------------------

xx = 0256;
% xx = 0016;
% xx = 0128;
% zz = 0256;
zz = 0128;
% zz = 0008;
yy = 0177; %nb of points in the channel
% yy = 0025; %nb of points in the channel

%nb of solid points per tile
dsty = 5;
dstx = 7;
% dsty = 3;
% dstx = 3;
shift_x = 2;
shift_z = 2;

%----------------------------------------
% Reading grid
%----------------------------------------
file = sprintf('output_uv_separate/y_grid_canopy_big.dat');
% file = sprintf('output_uv_separate/y_grid_canopy_small.dat');
fid = fopen(file,'r');    
Ngal = fread(fid,4*5,'int');
dny = fread(fid,1,'int'); %nb of points in canopy
yu = fread(fid,yy+2*dny+1,'real*8');
yv = fread(fid,yy+2*dny,'real*8');
dthdyu=fread(fid,yy+2*dny+1,'real*8');
dthdyv=fread(fid,yy+2*dny,'real*8');
dyu2i=fread(fid,3*(yy+2*dny+1),'real*8');
dyv2i= fread(fid,3*(yy+2*dny),'real*8');
fclose(fid);

Ngal = reshape(Ngal,4,5);
dyu2i = reshape(dyu2i,yy+2*dny+1,3);
dyv2i=  reshape(dyv2i,yy+2*dny,3);
dyu2 = diag(dyu2i(2:end,1),-1) + diag(dyu2i(:,2)) + diag(dyu2i(1:end-1,3),1);

%----------------------------------------
% Reading geometry
%----------------------------------------
%Filename to read
file = sprintf('output_uv_separate/boundary_%4.4dx%4.4dx%4.4d.dat',xx,zz,yy)
fid = fopen(file,'r');    
Lx = fread(fid,1,'real*8');
Ly = fread(fid,1,'real*8');
Lz = fread(fid,1,'real*8');
     
Ngal = fread(fid,4*5,'int');
Ngal = reshape(Ngal, [4 5]);
nlist_ib_s_bot_v = fread(fid,1,'int');
nlist_ib_f_bot_v = fread(fid,1,'int');
nlist_ib_s_bot_u = fread(fid,1,'int');
nlist_ib_f_bot_u = fread(fid,1,'int');

nyu11  = fread(fid,1,'int');
nyu21  = fread(fid,1,'int');
nyu12  = fread(fid,1,'int');
nyu22  = fread(fid,1,'int');

nyv11  = fread(fid,1,'int');
nyv21  = fread(fid,1,'int');
nyv12  = fread(fid,1,'int');
nyv22  = fread(fid,1,'int');    

list_ib_s_bot_v     = fread(fid,3*nlist_ib_s_bot_v,'int');
list_ib_s_bot_v = reshape(list_ib_s_bot_v,[3 nlist_ib_s_bot_v]);
list_ib_s_top_v     = fread(fid,3*nlist_ib_s_bot_v,'int');
list_ib_s_top_v = reshape(list_ib_s_top_v,[3 nlist_ib_s_bot_v]);

list_ib_s_bot_u     = fread(fid,3*nlist_ib_s_bot_u,'int');
list_ib_s_bot_u = reshape(list_ib_s_bot_u,[3 nlist_ib_s_bot_u]);
list_ib_s_top_u     = fread(fid,3*nlist_ib_s_bot_u,'int');
list_ib_s_top_u = reshape(list_ib_s_top_u,[3 nlist_ib_s_bot_u]);

list_ib_f_bot_v     = fread(fid,9*nlist_ib_f_bot_v,'int');
list_ib_f_bot_v = reshape(list_ib_f_bot_v,[9 nlist_ib_f_bot_v]);
list_ib_f_top_v     = fread(fid,9*nlist_ib_f_bot_v,'int');
list_ib_f_top_v = reshape(list_ib_f_top_v,[9 nlist_ib_f_bot_v]);

list_ib_f_bot_u     = fread(fid,9*nlist_ib_f_bot_u,'int');
list_ib_f_bot_u = reshape(list_ib_f_bot_u,[9 nlist_ib_f_bot_u]);
list_ib_f_top_u     = fread(fid,9*nlist_ib_f_bot_u,'int');
list_ib_f_top_u = reshape(list_ib_f_top_u,[9 nlist_ib_f_bot_u]);

list_ib_f_w_bot_v     = fread(fid,2*nlist_ib_f_bot_v,'real*8');
list_ib_f_w_bot_v = reshape(list_ib_f_w_bot_v,[2 nlist_ib_f_bot_v]);
list_ib_f_w_top_v     = fread(fid,2*nlist_ib_f_bot_v,'real*8');
list_ib_f_w_top_v = reshape(list_ib_f_w_top_v,[2 nlist_ib_f_bot_v]);

list_ib_f_w_bot_u     = fread(fid,2*nlist_ib_f_bot_u,'real*8');
list_ib_f_w_bot_u = reshape(list_ib_f_w_bot_u,[2 nlist_ib_f_bot_u]);
list_ib_f_w_top_u     = fread(fid,2*nlist_ib_f_bot_u,'real*8');
list_ib_f_w_top_u = reshape(list_ib_f_w_top_u,[2 nlist_ib_f_bot_u]);

list_ib_s_xz     = fread(fid,2*dstx*dstx,'int');
list_ib_s_xz = reshape(list_ib_s_xz,[2 dstx*dstx]);

list_ib_f_xz     = fread(fid,10*dstx*dstx,'real*8');
list_ib_f_xz = reshape(list_ib_f_xz,[10 dstx*dstx]);

fclose(fid)

%save('geom_canopy.mat','list_ib1','list_ib2','nlist_ib1','nlist_ib2')

% --------------------------------------------
% Plot one element
% --------------------------------------------

%x-z plane points one stem
R = (dstx-1)/2;
z_circ = linspace(1,dstx,1000);
x_circ_top = (R^2-(z_circ-(R+1)).^2).^0.5;
x_circ_bot = -x_circ_top +R+1;
x_circ_top = x_circ_top +R+1;
%Weigthing
% x_from_weights = (list_ib_f_xz(1,:)-list_ib_f_xz(3,:).*list_ib_f_xz(8,:) ...
%                 -list_ib_f_xz(5,:).*list_ib_f_xz(9,:))./list_ib_f_xz(7,:);
%             
% z_from_weights = (list_ib_f_xz(2,:)-list_ib_f_xz(4,:).*list_ib_f_xz(8,:) ...
%                 -list_ib_f_xz(6,:).*list_ib_f_xz(9,:))./list_ib_f_xz(7,:);
figure; hold on
plot(list_ib_s_xz(1,:),list_ib_s_xz(2,:),'k.','MarkerSize',30)
plot(list_ib_f_xz(1,:),list_ib_f_xz(2,:),'r.','MarkerSize',30)
plot(list_ib_f_xz(3,:),list_ib_f_xz(4,:),'bx','MarkerSize',15)
plot(list_ib_f_xz(5,:),list_ib_f_xz(6,:),'bx','MarkerSize',15)
plot(list_ib_f_xz(9,:),list_ib_f_xz(10,:),'g.','MarkerSize',20)
plot(x_circ_bot,z_circ,'g')
plot(x_circ_top,z_circ,'g')
% plot(x_from_weights,z_from_weights,'kx','MarkerSize',10)
title('XZ Plane One Stem');ylabel('z');xlabel('x');
axis equal
hold off;


%---------------------------------------------
% Bottom substrate
%---------------------------------------------

%x-z plane points 
%vgrid
% Weightings
x_from_weights = (list_ib_f_bot_v(1,:)-list_ib_f_bot_v(4,:).*list_ib_f_w_bot_v(1,:) ...
                -list_ib_f_bot_v(7,:).*list_ib_f_w_bot_v(2,:))./...
                (1-list_ib_f_w_bot_v(1,:)-list_ib_f_w_bot_v(2,:));      
z_from_weights = (list_ib_f_bot_v(2,:)-list_ib_f_bot_v(5,:).*list_ib_f_w_bot_v(1,:) ...
                -list_ib_f_bot_v(8,:).*list_ib_f_w_bot_v(2,:))./...
                (1-list_ib_f_w_bot_v(1,:)-list_ib_f_w_bot_v(2,:));
figure; hold on
plot(list_ib_s_bot_v(1,:),list_ib_s_bot_v(2,:),'k.','MarkerSize',10)
plot(list_ib_f_bot_v(1,:),list_ib_f_bot_v(2,:),'r.','MarkerSize',10)
plot(list_ib_f_bot_v(4,:),list_ib_f_bot_v(5,:),'b.','MarkerSize',10)
plot(list_ib_f_bot_v(7,:),list_ib_f_bot_v(8,:),'b.','MarkerSize',10)
% plot(x_circ_bot+shift_x,z_circ+shift_z,'g')
% plot(x_circ_top+shift_x,z_circ+shift_z,'g')
% plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

% plot(list_ib_s_bot_v_check(:,1),list_ib_s_bot_v_check(:,2),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,1),list_ib_f_bot_check_v(:,2),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,4),list_ib_f_bot_check_v(:,5),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,7),list_ib_f_bot_check_v(:,8),'bx','MarkerSize',10)
title('XZ Plane Vgrid Bottom');ylabel('z');xlabel('x');
axis equal
hold off;

%ugrid
% Weightings
x_from_weights = (list_ib_f_bot_u(1,:)-list_ib_f_bot_u(4,:).*list_ib_f_w_bot_u(1,:) ...
                -list_ib_f_bot_u(7,:).*list_ib_f_w_bot_u(2,:))./...
                (1-list_ib_f_w_bot_u(1,:)-list_ib_f_w_bot_u(1,:));      
z_from_weights = (list_ib_f_bot_u(2,:)-list_ib_f_bot_u(5,:).*list_ib_f_w_bot_u(1,:) ...
                -list_ib_f_bot_u(8,:).*list_ib_f_w_bot_u(2,:))./...
                (1-list_ib_f_w_bot_u(1,:)-list_ib_f_w_bot_u(1,:));
figure; hold on
plot(list_ib_s_bot_u(1,:),list_ib_s_bot_u(2,:),'k.','MarkerSize',10)
plot(list_ib_f_bot_u(1,:),list_ib_f_bot_u(2,:),'r.','MarkerSize',10)
plot(list_ib_f_bot_u(4,:),list_ib_f_bot_u(5,:),'b.','MarkerSize',10)
plot(list_ib_f_bot_u(7,:),list_ib_f_bot_u(8,:),'b.','MarkerSize',10)
% plot(x_circ_bot+shift_x,z_circ+shift_z,'g')
% plot(x_circ_top+shift_x,z_circ+shift_z,'g')
% plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

% plot(list_ib_s_bot_u_check(:,1),list_ib_s_bot_u_check(:,2),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,1),list_ib_f_bot_check_u(:,2),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,4),list_ib_f_bot_check_u(:,5),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,7),list_ib_f_bot_check_u(:,8),'bx','MarkerSize',10)
title('XZ Plane Ugrid Bottom');ylabel('z');xlabel('x');
axis equal
hold off;

%x-y plane points 
%vgrid
figure; hold on
plot(list_ib_s_bot_v(1,:),list_ib_s_bot_v(3,:),'k.','MarkerSize',10)
plot(list_ib_f_bot_v(1,:),list_ib_f_bot_v(3,:),'r.','MarkerSize',10)
plot(list_ib_f_bot_v(4,:),list_ib_f_bot_v(6,:),'b.','MarkerSize',10)
plot(list_ib_f_bot_v(7,:),list_ib_f_bot_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_bot_v_check(:,1),list_ib_s_bot_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,1),list_ib_f_bot_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,4),list_ib_f_bot_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,7),list_ib_f_bot_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_s_bot_u(1,:),list_ib_s_bot_u(3,:)-0.5,'kx','MarkerSize',10)
plot(list_ib_f_bot_u(1,:),list_ib_f_bot_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_f_bot_u(4,:),list_ib_f_bot_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_f_bot_u(7,:),list_ib_f_bot_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_bot_u_check(:,1),list_ib_s_bot_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,1),list_ib_f_bot_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,4),list_ib_f_bot_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,7),list_ib_f_bot_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('XY Plane Bottom');ylabel('y');xlabel('x');
hold off;

%z-y plane points 
%vgrid
figure; hold on
plot(list_ib_s_bot_v(2,:),list_ib_s_bot_v(3,:),'k.','MarkerSize',10)
plot(list_ib_f_bot_v(2,:),list_ib_f_bot_v(3,:),'r.','MarkerSize',10)
plot(list_ib_f_bot_v(5,:),list_ib_f_bot_v(6,:),'b.','MarkerSize',10)
plot(list_ib_f_bot_v(8,:),list_ib_f_bot_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_bot_v_check(:,2),list_ib_s_bot_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,2),list_ib_f_bot_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,5),list_ib_f_bot_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,8),list_ib_f_bot_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_s_bot_u(2,:),list_ib_s_bot_u(3,:)-0.5,'kx','MarkerSize',10)
plot(list_ib_f_bot_u(2,:),list_ib_f_bot_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_f_bot_u(5,:),list_ib_f_bot_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_f_bot_u(8,:),list_ib_f_bot_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_bot_u_check(:,2),list_ib_s_bot_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,2),list_ib_f_bot_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,5),list_ib_f_bot_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,8),list_ib_f_bot_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('ZY Plane Bottom');ylabel('y');xlabel('z');
hold off;

%---------------------------------------------
% Top substrate
%---------------------------------------------

%x-z plane points 
%vgrid
% Weightings
x_from_weights = (list_ib_f_top_v(1,:)-list_ib_f_top_v(4,:).*list_ib_f_w_top_v(1,:) ...
                -list_ib_f_top_v(7,:).*list_ib_f_w_top_v(2,:))./...
                (1-list_ib_f_w_top_v(1,:)-list_ib_f_w_top_v(2,:));      
z_from_weights = (list_ib_f_top_v(2,:)-list_ib_f_top_v(5,:).*list_ib_f_w_top_v(1,:) ...
                -list_ib_f_top_v(8,:).*list_ib_f_w_top_v(2,:))./...
                (1-list_ib_f_w_top_v(1,:)-list_ib_f_w_top_v(2,:));  
figure; hold on
plot(list_ib_s_top_v(1,:),list_ib_s_top_v(2,:),'k.','MarkerSize',10)
plot(list_ib_f_top_v(1,:),list_ib_f_top_v(2,:),'r.','MarkerSize',10)
plot(list_ib_f_top_v(4,:),list_ib_f_top_v(5,:),'b.','MarkerSize',10)
plot(list_ib_f_top_v(7,:),list_ib_f_top_v(8,:),'b.','MarkerSize',10)
% plot(x_circ_bot+shift_x,z_circ+shift_z,'g')
% plot(x_circ_top+shift_x,z_circ+shift_z,'g')
% plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

% plot(list_ib_s_top_v_check(:,1),list_ib_s_top_v_check(:,2),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,1),list_ib_f_top_check_v(:,2),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,4),list_ib_f_top_check_v(:,5),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,7),list_ib_f_top_check_v(:,8),'bx','MarkerSize',10)
title('XZ Plane Vgrid Top');ylabel('z');xlabel('x');
axis equal
hold off;

%ugrid
% Weightings
x_from_weights = (list_ib_f_top_u(1,:)-list_ib_f_top_u(4,:).*list_ib_f_w_top_u(1,:) ...
                -list_ib_f_top_u(7,:).*list_ib_f_w_top_u(2,:))./...
                (1-list_ib_f_w_top_u(1,:)-list_ib_f_w_top_u(2,:));      
z_from_weights = (list_ib_f_top_u(2,:)-list_ib_f_top_u(5,:).*list_ib_f_w_top_u(1,:) ...
                -list_ib_f_top_u(8,:).*list_ib_f_w_top_u(2,:))./...
                (1-list_ib_f_w_top_u(1,:)-list_ib_f_w_top_u(2,:)); 
figure; hold on
plot(list_ib_s_top_u(1,:),list_ib_s_top_u(2,:),'k.','MarkerSize',10)
plot(list_ib_f_top_u(1,:),list_ib_f_top_u(2,:),'r.','MarkerSize',10)
plot(list_ib_f_top_u(4,:),list_ib_f_top_u(5,:),'b.','MarkerSize',10)
plot(list_ib_f_top_u(7,:),list_ib_f_top_u(8,:),'b.','MarkerSize',10)
% plot(x_circ_bot+shift_x,z_circ+shift_z,'g')
% plot(x_circ_top+shift_x,z_circ+shift_z,'g')
% plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

% plot(list_ib_s_top_u_check(:,1),list_ib_s_top_u_check(:,2),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,1),list_ib_f_top_check_u(:,2),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,4),list_ib_f_top_check_u(:,5),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,7),list_ib_f_top_check_u(:,8),'bx','MarkerSize',10)
title('XZ Plane Ugrid Top');ylabel('z');xlabel('x');
axis equal
hold off;

%x-y plane points 
%vgrid
figure; hold on
plot(list_ib_s_top_v(1,:),list_ib_s_top_v(3,:),'k.','MarkerSize',10)
plot(list_ib_f_top_v(1,:),list_ib_f_top_v(3,:),'r.','MarkerSize',10)
plot(list_ib_f_top_v(4,:),list_ib_f_top_v(6,:),'b.','MarkerSize',10)
plot(list_ib_f_top_v(7,:),list_ib_f_top_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_top_v_check(:,1),list_ib_s_top_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,1),list_ib_f_top_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,4),list_ib_f_top_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,7),list_ib_f_top_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_s_top_u(1,:),list_ib_s_top_u(3,:)-0.5,'kx','MarkerSize',10)
plot(list_ib_f_top_u(1,:),list_ib_f_top_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_f_top_u(4,:),list_ib_f_top_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_f_top_u(7,:),list_ib_f_top_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_top_u_check(:,1),list_ib_s_top_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,1),list_ib_f_top_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,4),list_ib_f_top_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,7),list_ib_f_top_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('XY Plane Top');ylabel('y');xlabel('x');
hold off;

%z-y plane points 
%vgrid
figure; hold on
plot(list_ib_s_top_v(2,:),list_ib_s_top_v(3,:),'k.','MarkerSize',10)
plot(list_ib_f_top_v(2,:),list_ib_f_top_v(3,:),'r.','MarkerSize',10)
plot(list_ib_f_top_v(5,:),list_ib_f_top_v(6,:),'b.','MarkerSize',10)
plot(list_ib_f_top_v(8,:),list_ib_f_top_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_top_v_check(:,2),list_ib_s_top_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,2),list_ib_f_top_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,5),list_ib_f_top_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,8),list_ib_f_top_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_s_top_u(2,:),list_ib_s_top_u(3,:)-0.5,'kx','MarkerSize',10)
plot(list_ib_f_top_u(2,:),list_ib_f_top_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_f_top_u(5,:),list_ib_f_top_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_f_top_u(8,:),list_ib_f_top_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_top_u_check(:,2),list_ib_s_top_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,2),list_ib_f_top_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,5),list_ib_f_top_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,8),list_ib_f_top_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('ZY Plane Top');ylabel('y');xlabel('z');
hold off;

autoArrangeFigures
