clear all;

% --- Time Stamp Matrix ---
timematrix=[...
% 00000;... 
% 00399;... % fixed missing post (jamon )
% 06296; ... % CSD3 MPI_fixed_no_force 
% 03192; ... % CSD3 solid points only
% 03202; ... % CSD3 no overlap
% 00200; ... % CSD3 no overlap correct
% 00242; ... % CSD3 no overlap correct
% 00305; ... % CSD3 no overlap correct
% 01252; ... % CSD3 no overlap correct
% 01247; ... % CSD3 no overlap correct
% 00033; ... % CSD3 no neg u_f
% 04202; ... % CSD3 solid points only Akshaths trick
% 03673; ... % CSD3 no force u
% 00634; ... % CSD3 uv separate with forcing 
03763; ... % CSD3 uv separate without forcing
];

% --- Output Folder ---
% output = 'output_solid_points_only';
% output = 'output_no_overlap_after_correct';
% output = 'output_no_neg_u_f';
% output = 'output_solid_only_Akshaths_trick';
% output = 'output_no_force_u';
% output = 'output_debug_IB2';
output = 'output_uv_separate';

% --- Big Box Dimensions ---
xx = 0256;
zz = 0128;
yy = 0177;

% --- Small Box Dimensions ---
% xx = 0016;
% zz = 0008;
% yy = 0025;

% --- Height j of the slice ---
jlvl=3;

% --- Delete this? --- 
Re_t = 180;
Re = 3250;
utau = Re_t/Re; 
alpha=1/(2*2*pi);
beta=1/(2*2*pi);

for tt = 1:length(timematrix) 
    time = timematrix(tt)
    file = sprintf('%s/u1_phys_%4.4dx%4.4dx%4.4d_t%5.5d.dat',output,xx,zz,yy,time);
    %file

    %% Reading
    fid = fopen(file,'r');
    % a = fread(fid,[1,1],'int');
    h       = fread(fid,1,'*uint32');
    t       = fread(fid,1,'real*8');
    Re      = fread(fid,1,'real*8');
    alp     = fread(fid,1,'real*8');
    bet     = fread(fid,1,'real*8');
    mpgx    = fread(fid,1,'real*8');
    nband   = fread(fid,1,'int');
    iter    = fread(fid,1,'int');
    dummit  = fread(fid,90,'int');
    N       = fread(fid,4*5,'int');
    N       = reshape(N,4,5);
        fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
    yu       = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');
    dthetai = fread(fid,1,'real*8');
        fread(fid,1,'real*8'); %skip
    dthdyu   = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');


         %fread(fid,1,'int'); %skip
    for index=1:N(4,1+1)-N(4,0+1)+1
        j1(index)        = fread(fid,1,'int');
        one1(index)      = fread(fid,1,'int');
        nx1(index)       = fread(fid,1,'int');
        nz1(index)       = fread(fid,1,'int');
        yj1(index)       = fread(fid,1,'real*8');
        u11(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                          fread(fid,1,'real*8');
    end
    for index=1:N(4,2+1)-N(4,1+1)+1
        j2(index)        = fread(fid,1,'int');
        one2(index)      = fread(fid,1,'int');
        nx2(index)       = fread(fid,1,'int');
        nz2(index)       = fread(fid,1,'int');
        yj2(index)       = fread(fid,1,'real*8');
        u12(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                          fread(fid,1,'real*8');
    end
    for index=1:N(4,3+1)-N(4,2+1)+2-2+2
        j3(index)        = fread(fid,1,'int');
        one3(index)      = fread(fid,1,'int');
        nx3(index)       = fread(fid,1,'int');
        nz3(index)       = fread(fid,1,'int');
        yj3(index)       = fread(fid,1,'real*8');
        u13(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                          fread(fid,1,'real*8');
    end
    fclose(fid);

    u1 = [u11(:,1);u12(:,1);u13(:,1)];

    file = sprintf('%s/u2_phys_%4.4dx%4.4dx%4.4d_t%5.5d.dat',output,xx,zz,yy,time);
    %file

    %% Reading
    fid = fopen(file,'r');
    % a = fread(fid,[1,1],'int');
    h       = fread(fid,1,'*uint32');
    t       = fread(fid,1,'real*8');
    Re      = fread(fid,1,'real*8');
    alp     = fread(fid,1,'real*8');
    bet     = fread(fid,1,'real*8');
    mpgx    = fread(fid,1,'real*8');
    nband   = fread(fid,1,'int');
    iter    = fread(fid,1,'int');
    dummit  = fread(fid,90,'int');
    N       = fread(fid,4*5,'int');
    N       = reshape(N,4,5);
        fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
    yv       = fread(fid,N(3,nband+1)-N(3,0+1)+2,'real*8');
    dthetai = fread(fid,1,'real*8');
        fread(fid,1,'real*8'); %skip
    dthdyv   = fread(fid,N(3,nband+1)-N(3,0+1)+2,'real*8');


         %fread(fid,1,'int'); %skip
    for index=1:N(3,1+1)-N(3,0+1)+1
        j1(index)        = fread(fid,1,'int');
        one1(index)      = fread(fid,1,'int');
        nx1(index)       = fread(fid,1,'int');
        nz1(index)       = fread(fid,1,'int');
        yj1(index)       = fread(fid,1,'real*8');
        u21(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                          fread(fid,1,'real*8');
    end
    for index=1:N(3,2+1)-N(3,1+1)+1
        j2(index)        = fread(fid,1,'int');
        one2(index)      = fread(fid,1,'int');
        nx2(index)       = fread(fid,1,'int');
        nz2(index)       = fread(fid,1,'int');
        yj2(index)       = fread(fid,1,'real*8');
        u22(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                          fread(fid,1,'real*8');
    end
    for index=1:N(3,3+1)-N(3,2+1)+2-2
        j3(index)        = fread(fid,1,'int');
        one3(index)      = fread(fid,1,'int');
        nx3(index)       = fread(fid,1,'int');
        nz3(index)       = fread(fid,1,'int');
        yj3(index)       = fread(fid,1,'real*8');
        u23(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                          fread(fid,1,'real*8');
    end
    fclose(fid);

    u2 = [u21(:,1);u22(:,1);u23(:,1)];

    file = sprintf('%s/u3_phys_%4.4dx%4.4dx%4.4d_t%5.5d.dat',output,xx,zz,yy,time);
    %file

    %% Reading
    fid = fopen(file,'r');
    % a = fread(fid,[1,1],'int');
    h       = fread(fid,1,'*uint32');
    t       = fread(fid,1,'real*8');
    Re      = fread(fid,1,'real*8');
    alp     = fread(fid,1,'real*8');
    bet     = fread(fid,1,'real*8');
    mpgx    = fread(fid,1,'real*8');
    nband   = fread(fid,1,'int');
    iter    = fread(fid,1,'int');
    dummit  = fread(fid,90,'int');
    N       = fread(fid,4*5,'int');
    N       = reshape(N,4,5);
        fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
    yu       = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');
    dthetai = fread(fid,1,'real*8');
        fread(fid,1,'real*8'); %skip
    dthdyu   = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');


         %fread(fid,1,'int'); %skip
    for index=1:N(4,1+1)-N(4,0+1)+1
        j1(index)        = fread(fid,1,'int');
        one1(index)      = fread(fid,1,'int');
        nx1(index)       = fread(fid,1,'int');
        nz1(index)       = fread(fid,1,'int');
        yj1(index)       = fread(fid,1,'real*8');
        u31(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                          fread(fid,1,'real*8');
    end
    for index=1:N(4,2+1)-N(4,1+1)+1
        j2(index)        = fread(fid,1,'int');
        one2(index)      = fread(fid,1,'int');
        nx2(index)       = fread(fid,1,'int');
        nz2(index)       = fread(fid,1,'int');
        yj2(index)       = fread(fid,1,'real*8');
        u32(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                          fread(fid,1,'real*8');
    end
    for index=1:N(4,3+1)-N(4,2+1)+2-2
        j3(index)        = fread(fid,1,'int');
        one3(index)      = fread(fid,1,'int');
        nx3(index)       = fread(fid,1,'int');
        nz3(index)       = fread(fid,1,'int');
        yj3(index)       = fread(fid,1,'real*8');
        u33(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                          fread(fid,1,'real*8');
    end
    fclose(fid);

    u3 = [u31(:,1);u32(:,1);u33(:,1)];




%% U

    usurf=0.5*(u11(jlvl,:)+u13(end-jlvl,:));
    usurf=u11(jlvl,:);
%     usurf=u13(end-jlvl,:);
    usurf=reshape(usurf,xx+2,zz);
    usurf=usurf(1:xx,:);
    usurf = transpose(usurf);
    %usurf(usurf<0) = 0;
    %usurf = usurf./max(usurf(:));
    figure(1)
    clf
    xl=Re_t*[1:xx]/(xx)*2*pi/alp;
    zl=Re_t*[1:zz]/(zz)*2*pi/bet;
%     c1 = contourf(zl,xl,usurf,15);
%     c1 = contourf(zl,xl,[usurf(18:24,17:24),usurf(18:24,1:16);usurf(1:17,17:24),usurf(1:17,1:16)],15);
    
%     midptx = solidptsx/2 + 1;
%     midptz = solidptsz/2 + 1;
%     shiftx = xx/2 + 1 - midptx;
%     shiftz = zz/2 + 1 - midptz;
%     
%     uplot = [usurf(:,end-shiftz+1:end),usurf(:,1:end-shiftz+1)];
%     uplot = [uplot(end-shiftx+1:end,:);uplot(1:end-shiftx+1,:)];
    uplot = usurf;
        
    %uplot = [usurf(18:24,17:24),usurf(18:24,1:17);usurf(1:18,17:24),usurf(1:18,1:17)]
    contour_levels = -1:1:1;
    c1 = contourf(xl,zl,uplot/utau,'linestyle','none');
%     c1 = contourf(zl,xl,uplot/utau,contour_levels,'linestyle','none');
%     c1 = contour(zl,xl,uplot,contour_levels,'k','linewidth',1);
    hold on;
    load uv_sep.mat
    plot((Re_t/(xx*1.5)*2*pi/alp)*list_ib_s_bot_u(1,:),(Re_t/(zz*1.5)*2*pi/bet)*list_ib_s_bot_u(2,:),'w.','MarkerSize',4) 
    hold off;
    set(gca,'Fontn','Times','FontSize',20,'linewidth',1);
    box on
%     print -deps2c ../../OneDrive/'My Documents'/'Fourth Year Project'/'Final Report'/XZContourU.eps
    colorbar
    %shading flat
    shading interp
    axis equal; ylabel('z^+');xlabel('x^+');
    %view(0,90)
    colormap('parula')
%     colormap('pink')
    
    %% V 
    vsurf = u21(jlvl,:);
    vsurf=reshape(vsurf,xx+2,zz);
    vsurf=vsurf(1:xx,:);
    vsurf = transpose(vsurf);
    figure(2)
    clf
    xl=Re_t*[0:xx-1]/(xx)*2*pi/alp;
    zl=Re_t*[0:zz-1]/(zz)*2*pi/bet;
    
    vplot = vsurf;
    contour_levels = -1:0.05:1;
    c1 = contourf(xl,zl,vplot/utau,'linestyle','none');     
%     c1 = contourf(zl,xl,uplot,contour_levels,'linestyle','none');
% %     c1 = contour(zl,xl,uplot,contour_levels,'k');
    hold on;
    plot((Re_t/(xx*1.5)*2*pi/alp)*list_ib_s_bot_u(1,:),(Re_t/(zz*1.5)*2*pi/bet)*list_ib_s_bot_u(2,:),'w.','MarkerSize',4) 
    hold off;
    colorbar
    %shading flat
    shading interp
    axis equal; ylabel('z^+');xlabel('x^+');
    %view(0,90)
    colormap('parula')
%% W
    wsurf=0.5*(u31(jlvl,:)+u33(end-jlvl,:));
    wsurf = u31(jlvl,:);
%     wsurf = u33(end-jlvl,:);
    wsurf=reshape(wsurf,xx+2,zz);
    wsurf=wsurf(1:xx,:);
    wsurf = transpose(wsurf);
    figure(3)
    clf
    xl=Re_t*[0:xx-1]/(xx)*2*pi/alp;
    zl=Re_t*[0:zz-1]/(zz)*2*pi/bet;
    
%     midptx = solidptsx/2 + 1;
%     midptz = solidptsz/2 + 1;
%     shiftx = xx/2 + 1 - midptx;
%     shiftz = zz/2 + 1 - midptz;
%     
%     wplot = [wsurf(:,end-shiftz+1:end),wsurf(:,1:end-shiftz+1)];
%     wplot = [wplot(end-shiftx+1:end,:);wplot(1:end-shiftx+1,:)];
    wplot = wsurf;
    contour_levels = -5:5:5;
    c1 = contourf(xl,zl,wplot/utau,'linestyle','none');
%     c1 = contourf(zl,xl,wplot/utau,contour_levels,'linestyle','none');
%     c1 = contour(zl,xl,uplot,contour_levels,'k');
    hold on;
    plot((Re_t/(xx*1.5)*2*pi/alp)*list_ib_s_bot_u(1,:),(Re_t/(zz*1.5)*2*pi/bet)*list_ib_s_bot_u(2,:),'w.','MarkerSize',4) 
    hold off;   
    colorbar
    %shading flat
    shading interp
    axis equal; ylabel('z^+');xlabel('x^+');
    %view(0,90)
    colormap('parula')

%%    
    %load ../Code_170112_Compsl12_ref/new_postpro/unify_stats/UmC
    %load ../Code_170112_Compsl24/new_postpro/unify_stats/UmC
    %load ../Code_170112_Compsl35/new_postpro/unify_stats/UmC
    %load ../Code_170112_Compsl47/new_postpro/unify_stats/UmC
    
%     hold on
%     UmC_s2 = UmC(:,:,1);
%     UmC_s = UmC(:,:,1);
%     UmC_s(UmC_s<0) = 0;
%     UmC_s = UmC_s./max(UmC_s(:));
%     contour(zl,xl,UmC_s,15,'w--','LineWidth',2)
%     
%     uPL1 = u11(1,:,:);
%     uPL1=reshape(uPL1,xx+2,zz);
%     uPL1 = uPL1(1:xx,:);
%     uPL2 = u11(2,:,:);
%     uPL2=reshape(uPL2,xx+2,zz);
%     uPL2 = uPL2(1:xx,:);
%     
%     
%     Re_DNS = mean(mean(UmC_s2))*(2*pi/1/48)*3250
%     Re_lam = mean(mean(uPL1))*(2*pi/alp)*Re
end