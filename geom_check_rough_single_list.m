clear all;
close all;

%----------------------------------------
% Input values
%----------------------------------------

% --- Big Box ---
% xx = 0256;
% zz = 0128;
% yy = 0177; 

% --- Big Big Box ---
xx = 0512;
zz = 0256;
yy = 0177; 

% xx = 0016;
% zz = 0008;
% yy = 0025; 

% size of diameter in grid points 
% dstx = 11;
dstx = 23;
% amount first cylinder is shifted from (0,0)
shift_x = 2;
shift_z = 2;

%----------------------------------------
% Reading grid
%----------------------------------------
% file = sprintf('boundary_K18/y_grid_canopy.dat');
% file = sprintf('boundary_force_inside_even_diameter/y_grid_rough.dat');
% % file = sprintf('output_Lap_coef_fixed_RHS/y_grid_canopy.dat');
% fid = fopen(file,'r');    
% Ngal = fread(fid,4*5,'int');
% dny = fread(fid,1,'int'); %nb of points in canopy
% yu = fread(fid,yy+2*dny+1,'real*8');
% yv = fread(fid,yy+2*dny,'real*8');
% dthdyu=fread(fid,yy+2*dny+1,'real*8');
% dthdyv=fread(fid,yy+2*dny,'real*8');
% dyu2i=fread(fid,3*(yy+2*dny+1),'real*8');
% dyv2i= fread(fid,3*(yy+2*dny),'real*8');
% fclose(fid);
% 
% Ngal = reshape(Ngal,4,5);
% dyu2i = reshape(dyu2i,yy+2*dny+1,3);
% dyv2i=  reshape(dyv2i,yy+2*dny,3);
% dyu2 = diag(dyu2i(2:end,1),-1) + diag(dyu2i(:,2)) + diag(dyu2i(1:end-1,3),1);
% 
% figure;
% subplot(2,1,1); hold on;
% endd = dny+4; Re_t = 185;
% plot(Re_t*yv(2:endd+1),Re_t*(yv(3:endd+2)-yv(1:endd))/2, 'gx');
% plot(Re_t*yu(2:endd+1),Re_t*(yu(3:endd+2)-yu(1:endd))/2, 'mx');
% plot(Re_t*yv(1:endd+1),Re_t./dthdyv(1:endd+1),'go')
% plot(Re_t*yu(1:endd+1),Re_t./dthdyu(1:endd+1),'mo');
% hold off;
% title('Bot'); ylabel('(dy/dj)^+'); xlabel('y^+'); legend('vgrid','ugrid','Location','southeast');
% subplot(2,1,2); hold on;
% sta = length(yu)-endd;
% plot(Re_t*yv(sta-1:end-1),Re_t*(yv(sta:end)-yv(sta-2:end-2))/2, 'gx');
% plot(Re_t*yu(sta:end-1),Re_t*(yu(sta+1:end)-yu(sta-1:end-2))/2, 'mx');
% plot(Re_t*yv(sta-1:end),Re_t./dthdyv(sta-1:end),'go')
% plot(Re_t*yu(sta:end),Re_t./dthdyu(sta:end),'mo');
% hold off;
% title('Top'); ylabel('(dy/dj)^+'); xlabel('y^+'); legend('vgrid','ugrid');
% 
% figure; hold on; 
% plot(-dny:(length(yv)-dny-1), yv', '-gx');
% plot(-dny-0.5:(length(yu)-dny-1.5), yu', 'mx');
% hold off;

%----------------------------------------
% Reading geometry
%----------------------------------------
%Filename to read
file = sprintf('boundary_force_inside_even_diameter/boundary_%4.4dx%4.4dx%4.4d.dat',xx,zz,yy)
% file = sprintf('boundary_force_inside_odd_diameter/boundary_%4.4dx%4.4dx%4.4d.dat',xx,zz,yy)
% file = sprintf('boundary_K18/boundary_%4.4dx%4.4dx%4.4d.dat',xx,zz,yy)
fid = fopen(file,'r');    
Lx = fread(fid,1,'real*8');
Ly = fread(fid,1,'real*8');
Lz = fread(fid,1,'real*8');
     
Ngal = fread(fid,4*5,'int');
Ngal = reshape(Ngal, [4 5]);
nlist_ib_bot_v = fread(fid,1,'int');
nlist_ib_bot_u = fread(fid,1,'int');

nyu11  = fread(fid,1,'int');
nyu21  = fread(fid,1,'int');
nyu12  = fread(fid,1,'int');
nyu22  = fread(fid,1,'int');

nyv11  = fread(fid,1,'int');
nyv21  = fread(fid,1,'int');
nyv12  = fread(fid,1,'int');
nyv22  = fread(fid,1,'int');    

list_ib_bot_v     = fread(fid,9*nlist_ib_bot_v,'int');
list_ib_bot_v = reshape(list_ib_bot_v,[9 nlist_ib_bot_v]);
list_ib_top_v     = fread(fid,9*nlist_ib_bot_v,'int');
list_ib_top_v = reshape(list_ib_top_v,[9 nlist_ib_bot_v]);

list_ib_bot_u     = fread(fid,9*nlist_ib_bot_u,'int');
list_ib_bot_u = reshape(list_ib_bot_u,[9 nlist_ib_bot_u]);
list_ib_top_u     = fread(fid,9*nlist_ib_bot_u,'int');
list_ib_top_u = reshape(list_ib_top_u,[9 nlist_ib_bot_u]);

list_ib_w_bot_v     = fread(fid,3*nlist_ib_bot_v,'real*8');
list_ib_w_bot_v = reshape(list_ib_w_bot_v,[3 nlist_ib_bot_v]);
list_ib_w_top_v     = fread(fid,3*nlist_ib_bot_v,'real*8');
list_ib_w_top_v = reshape(list_ib_w_top_v,[3 nlist_ib_bot_v]);

list_ib_w_bot_u     = fread(fid,3*nlist_ib_bot_u,'real*8');
list_ib_w_bot_u = reshape(list_ib_w_bot_u,[3 nlist_ib_bot_u]);
list_ib_w_top_u     = fread(fid,3*nlist_ib_bot_u,'real*8');
list_ib_w_top_u = reshape(list_ib_w_top_u,[3 nlist_ib_bot_u]);

list_ib_xz     = fread(fid,6*dstx*dstx,'int');
list_ib_xz = reshape(list_ib_xz,[6 dstx*dstx]);

list_ib_w_xz     = fread(fid,5*dstx*dstx,'real*8');
list_ib_w_xz = reshape(list_ib_w_xz,[5 dstx*dstx]);

fclose(fid);

%save('geom_canopy.mat','list_ib1','list_ib2','nlist_ib1','nlist_ib2')

% --------------------------------------------
% Plot one element
% --------------------------------------------

%x-z plane points one stem
R = (dstx-2)/2;
c = 1.5;
%Weigthing
x_from_weights = (list_ib_xz(1,:)-list_ib_xz(3,:).*list_ib_w_xz(1,:) ...
                -list_ib_xz(5,:).*list_ib_w_xz(2,:))./...
                (1-list_ib_w_xz(1,:)-list_ib_w_xz(2,:));      
z_from_weights = (list_ib_xz(2,:)-list_ib_xz(4,:).*list_ib_w_xz(1,:) ...
                -list_ib_xz(6,:).*list_ib_w_xz(2,:))./...
                (1-list_ib_w_xz(1,:)-list_ib_w_xz(2,:));
figure; hold on
plot(list_ib_xz(1,:),list_ib_xz(2,:),'r.','MarkerSize',30)
plot(list_ib_xz(3,:),list_ib_xz(4,:),'b.','MarkerSize',30)
plot(list_ib_xz(5,:),list_ib_xz(6,:),'b.','MarkerSize',30)
plot(list_ib_w_xz(3,:),list_ib_w_xz(4,:),'k.','MarkerSize',20)
plot(x_from_weights,z_from_weights,'g.','MarkerSize',15)
rectangle('Position',[c, c, R*2,R*2],'Curvature',[1,1],'LineStyle','-.','LineWidth',1);
title('XZ Plane One Stem');ylabel('z');xlabel('x');
axis equal
hold off;


%---------------------------------------------
% Bottom substrate
%---------------------------------------------

%x-z plane points 
%vgrid
% Weightings
x_from_weights = (list_ib_bot_v(1,:)-list_ib_bot_v(4,:).*list_ib_w_bot_v(2,:) ...
                -list_ib_bot_v(7,:).*list_ib_w_bot_v(3,:))./...
                (1-list_ib_w_bot_v(2,:)-list_ib_w_bot_v(3,:));      
z_from_weights = (list_ib_bot_v(2,:)-list_ib_bot_v(5,:).*list_ib_w_bot_v(2,:) ...
                -list_ib_bot_v(8,:).*list_ib_w_bot_v(3,:))./...
                (1-list_ib_w_bot_v(2,:)-list_ib_w_bot_v(3,:));
figure; hold on
plot(list_ib_bot_v(1,:),list_ib_bot_v(2,:),'r.','MarkerSize',10)
plot(list_ib_bot_v(4,:),list_ib_bot_v(5,:),'b.','MarkerSize',10)
plot(list_ib_bot_v(7,:),list_ib_bot_v(8,:),'b.','MarkerSize',10)
plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

% plot(list_ib_s_bot_v_check(:,1),list_ib_s_bot_v_check(:,2),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,1),list_ib_f_bot_check_v(:,2),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,4),list_ib_f_bot_check_v(:,5),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,7),list_ib_f_bot_check_v(:,8),'bx','MarkerSize',10)
title('XZ Plane Vgrid Bottom');ylabel('z');xlabel('x');
axis equal
hold off;

%ugrid
% Weightings
x_from_weights = (list_ib_bot_u(1,:)-list_ib_bot_u(4,:).*list_ib_w_bot_u(2,:) ...
                -list_ib_bot_u(7,:).*list_ib_w_bot_u(3,:))./...
                (1-list_ib_w_bot_u(2,:)-list_ib_w_bot_u(3,:));      
z_from_weights = (list_ib_bot_u(2,:)-list_ib_bot_u(5,:).*list_ib_w_bot_u(2,:) ...
                -list_ib_bot_u(8,:).*list_ib_w_bot_u(3,:))./...
                (1-list_ib_w_bot_u(2,:)-list_ib_w_bot_u(3,:));
figure; hold on
plot(list_ib_bot_u(1,:),list_ib_bot_u(2,:),'r.','MarkerSize',10)
plot(list_ib_bot_u(4,:),list_ib_bot_u(5,:),'b.','MarkerSize',10)
plot(list_ib_bot_u(7,:),list_ib_bot_u(8,:),'b.','MarkerSize',10)
plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

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
plot(list_ib_bot_v(1,:),list_ib_bot_v(3,:),'r.','MarkerSize',10)
plot(list_ib_bot_v(4,:),list_ib_bot_v(6,:),'b.','MarkerSize',10)
plot(list_ib_bot_v(7,:),list_ib_bot_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_bot_v_check(:,1),list_ib_s_bot_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,1),list_ib_f_bot_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,4),list_ib_f_bot_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,7),list_ib_f_bot_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_bot_u(1,:),list_ib_bot_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_bot_u(4,:),list_ib_bot_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_bot_u(7,:),list_ib_bot_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_bot_u_check(:,1),list_ib_s_bot_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,1),list_ib_f_bot_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,4),list_ib_f_bot_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_bot_check_u(:,7),list_ib_f_bot_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('XY Plane Bottom');ylabel('y');xlabel('x');
hold off;

%z-y plane points 
%vgrid
figure; hold on
plot(list_ib_bot_v(2,:),list_ib_bot_v(3,:),'r.','MarkerSize',10)
plot(list_ib_bot_v(5,:),list_ib_bot_v(6,:),'b.','MarkerSize',10)
plot(list_ib_bot_v(8,:),list_ib_bot_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_bot_v_check(:,2),list_ib_s_bot_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,2),list_ib_f_bot_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,5),list_ib_f_bot_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_bot_check_v(:,8),list_ib_f_bot_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_bot_u(2,:),list_ib_bot_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_bot_u(5,:),list_ib_bot_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_bot_u(8,:),list_ib_bot_u(9,:)-0.5,'bx','MarkerSize',10)

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
x_from_weights = (list_ib_top_v(1,:)-list_ib_top_v(4,:).*list_ib_w_top_v(2,:) ...
                -list_ib_top_v(7,:).*list_ib_w_top_v(3,:))./...
                (1-list_ib_w_top_v(2,:)-list_ib_w_top_v(3,:));      
z_from_weights = (list_ib_top_v(2,:)-list_ib_top_v(5,:).*list_ib_w_top_v(2,:) ...
                -list_ib_top_v(8,:).*list_ib_w_top_v(3,:))./...
                (1-list_ib_w_top_v(2,:)-list_ib_w_top_v(3,:));  
figure; hold on
plot(list_ib_top_v(1,:),list_ib_top_v(2,:),'r.','MarkerSize',10)
plot(list_ib_top_v(4,:),list_ib_top_v(5,:),'b.','MarkerSize',10)
plot(list_ib_top_v(7,:),list_ib_top_v(8,:),'b.','MarkerSize',10)
plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

% plot(list_ib_s_top_v_check(:,1),list_ib_s_top_v_check(:,2),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,1),list_ib_f_top_check_v(:,2),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,4),list_ib_f_top_check_v(:,5),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,7),list_ib_f_top_check_v(:,8),'bx','MarkerSize',10)
title('XZ Plane Vgrid Top');ylabel('z');xlabel('x');
axis equal
hold off;

%ugrid
% Weightings
x_from_weights = (list_ib_top_u(1,:)-list_ib_top_u(4,:).*list_ib_w_top_u(2,:) ...
                -list_ib_top_u(7,:).*list_ib_w_top_u(3,:))./...
                (1-list_ib_w_top_u(2,:)-list_ib_w_top_u(3,:));      
z_from_weights = (list_ib_top_u(2,:)-list_ib_top_u(5,:).*list_ib_w_top_u(2,:) ...
                -list_ib_top_u(8,:).*list_ib_w_top_u(3,:))./...
                (1-list_ib_w_top_u(2,:)-list_ib_w_top_u(3,:)); 
figure; hold on
plot(list_ib_top_u(1,:),list_ib_top_u(2,:),'r.','MarkerSize',10)
plot(list_ib_top_u(4,:),list_ib_top_u(5,:),'b.','MarkerSize',10)
plot(list_ib_top_u(7,:),list_ib_top_u(8,:),'b.','MarkerSize',10)
plot(x_from_weights,z_from_weights,'g.','MarkerSize',5)

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
plot(list_ib_top_v(1,:),list_ib_top_v(3,:),'r.','MarkerSize',10)
plot(list_ib_top_v(4,:),list_ib_top_v(6,:),'b.','MarkerSize',10)
plot(list_ib_top_v(7,:),list_ib_top_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_top_v_check(:,1),list_ib_s_top_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,1),list_ib_f_top_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,4),list_ib_f_top_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,7),list_ib_f_top_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_top_u(1,:),list_ib_top_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_top_u(4,:),list_ib_top_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_top_u(7,:),list_ib_top_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_top_u_check(:,1),list_ib_s_top_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,1),list_ib_f_top_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,4),list_ib_f_top_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,7),list_ib_f_top_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('XY Plane Top');ylabel('y');xlabel('x');
hold off;

%z-y plane points 
%vgrid
figure; hold on
plot(list_ib_top_v(2,:),list_ib_top_v(3,:),'r.','MarkerSize',10)
plot(list_ib_top_v(5,:),list_ib_top_v(6,:),'b.','MarkerSize',10)
plot(list_ib_top_v(8,:),list_ib_top_v(9,:),'b.','MarkerSize',10)

% plot(list_ib_s_top_v_check(:,2),list_ib_s_top_v_check(:,3),'kx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,2),list_ib_f_top_check_v(:,3),'rx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,5),list_ib_f_top_check_v(:,6),'bx','MarkerSize',10)
% plot(list_ib_f_top_check_v(:,8),list_ib_f_top_check_v(:,9),'bx','MarkerSize',10)
%ugrid
plot(list_ib_top_u(2,:),list_ib_top_u(3,:)-0.5,'rx','MarkerSize',10)
plot(list_ib_top_u(5,:),list_ib_top_u(6,:)-0.5,'bx','MarkerSize',10)
plot(list_ib_top_u(8,:),list_ib_top_u(9,:)-0.5,'bx','MarkerSize',10)

% plot(list_ib_s_top_u_check(:,2),list_ib_s_top_u_check(:,3)-0.5,'k.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,2),list_ib_f_top_check_u(:,3)-0.5,'r.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,5),list_ib_f_top_check_u(:,6)-0.5,'b.','MarkerSize',10)
% plot(list_ib_f_top_check_u(:,8),list_ib_f_top_check_u(:,9)-0.5,'b.','MarkerSize',10)
title('ZY Plane Top');ylabel('y');xlabel('z');
hold off;

autoArrangeFigures

%---------------------------------------------
% Check Lap coef
%--------------------------------------------- 

% vgrid bot
for ii = 1:length(list_ib_w_bot_v)
    if list_ib_bot_v(3,ii)==nyv21
        if list_ib_w_bot_v(1,ii) ~= 1
            disp('ERROR1');
        end 
    else
        if list_ib_w_bot_v(1,ii) == 0
            if list_ib_w_bot_v(2,ii) ~= 0
                disp('ERROR2');
            end 
        elseif list_ib_w_bot_v(1,ii) == 1
            if list_ib_w_bot_v(2,ii) == 0
                disp('ERROR3');
            end
        else
            disp('ERROR4');
        end
    end
end  

% ugrid bot
for ii = 1:length(list_ib_w_bot_u)
    if list_ib_bot_u(3,ii)==nyu21
        if list_ib_w_bot_u(1,ii) ~= 1
            disp('ERROR1');
        end 
    else
        if list_ib_w_bot_u(1,ii) == 0
            if list_ib_w_bot_u(2,ii) ~= 0
                disp('ERROR2');
            end 
        elseif list_ib_w_bot_u(1,ii) == 1
            if list_ib_w_bot_u(2,ii) == 0
                disp('ERROR3');
            end
        else
            disp('ERROR4');
        end
    end
end

% vgrid top
for ii = 1:length(list_ib_w_top_v)
    if list_ib_top_v(3,ii)==nyv12
        if list_ib_w_top_v(1,ii) ~= 1
            disp('ERROR1');
        end 
    else
        if list_ib_w_top_v(1,ii) == 0
            if list_ib_w_top_v(2,ii) ~= 0
                disp('ERROR2');
            end 
        elseif list_ib_w_top_v(1,ii) == 1
            if list_ib_w_top_v(2,ii) == 0
                disp('ERROR3');
            end
        else
            disp('ERROR4');
        end
    end
end  

% ugrid top
for ii = 1:length(list_ib_w_top_u)
    if list_ib_top_u(3,ii)==nyu12
        if list_ib_w_top_u(1,ii) ~= 1
            disp('ERROR1');
        end 
    else
        if list_ib_w_top_u(1,ii) == 0
            if list_ib_w_top_u(2,ii) ~= 0
                disp('ERROR2');
            end 
        elseif list_ib_w_top_u(1,ii) == 1
            if list_ib_w_top_u(2,ii) == 0
                disp('ERROR3');
            end
        else
            disp('ERROR4');
        end
    end
end  
    


disp('End of Lap coef check ')

