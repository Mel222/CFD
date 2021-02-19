clear all;

timematrix=[...

% 00218;...
% 00399;... % fixed missing post (jamon )
% 06296; ... % CSD3 MPI_fixed_no_force 
% 03192; ... % CSD3 solid points only
03202; ... % CSD3 no overlap
];

output = 'output_no_overlap_after';

xx = 0256;
% xx = 0128;
zz = 0128;
yy = 0177;

%solidptsx = xx/3;
%zval = solidptsx/2 + 1; % middle of solid 
% zval = solidptsx*2 + 1; % middle of air?  
ylim = -100;
zval = 5;
Re_t = 180;
Re = 3250;
utau = Re_t/Re;
alpha=1/(2*2*pi);
beta=1/(2*2*pi);

jlvl=1;

Re = 1;

ntx = 1;
ntz = 1;

pptx = xx;
pptz = zz;


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

%     usurf = zeros(80,24);
%     for i = 1:6 % P driven
    for i = 1:length(u11(:,1))
        usurf1=u11(i,:);
        usurf1=reshape(usurf1,xx+2,zz);
        usurf1=usurf1(1:xx,zval);
%         usurf1(usurf1<0) = 0;
        usurf(i,1:xx) = usurf1;
    end
% %     for i = 2:67
%     % Need to fix this as mid band has dif Ngal 
%     for i = 2:length(u12(:,1)) % P driven    
%         usurf2=u12(i,:);
%         usurf2=reshape(usurf2,xx+2,zz);
%         usurf2=usurf2(1:xx,zval);
% %         usurf2(usurf2<0) = 0;
%         usurf(i-1+length(u11(:,1)),1:xx) = usurf2;
%     end
%     for i = 2:7
    for i = 2:length(u13(:,1))
        usurf3=u13(i,:);
        usurf3=reshape(usurf3,xx+2,zz);
        usurf3=usurf3(1:xx,zval);
%         usurf3(usurf3<0) = 0;
        usurf(i-1+length(u11(:,1))+length(u12(:,1)),1:xx) = usurf3;
    end
    
    usurf = usurf(1:(yy+1),:);
%     usurf = usurf./max(usurf(:));

%    midptx = solidptsx/2 + 1;
%    shiftx = xx/2 + 1 - midptx;
%    uplot = [usurf(:,end-shiftx+1:end),usurf(:,1:end-shiftx+1)];
    uplot = usurf;
    figure(1)
    clf
    xl=Re_t*[0:xx-1]/(xx)*2*pi/alp; 
    contour_levels = [0:2:200];
    c1 = contourf(xl,Re_t*yu(yu<ylim/Re_t),uplot(yu<ylim/Re_t,:)/utau,contour_levels,'linestyle','none');
%     c1 = contour(xl,yu(1:22),uplot(1:22,:),contour_levels,'k');
    set(gca,'Fontn','Times','FontSize',20,'linewidth',1);
    box on
    colorbar
    %shading flat
    shading interp
    axis equal;xlabel('x^+');ylabel('y^+'); 
    %view(0,90)
    colormap('parula')
    
    %% V 
%     vsurf = zeros(80,24);
    for i = 1:length(u21(:,1))
        vsurf1=u21(i,:);
        vsurf1=reshape(vsurf1,xx+2,zz);
        vsurf1=vsurf1(1:xx,zval);
        vsurf(i,1:xx) = vsurf1;
    end
    % Need to fix this as mid band has dif Ngal 
%     for i = 2:length(u22(:,1))   
%         vsurf2=u22(i,:);
%         vsurf2=reshape(vsurf2,xx+2,zz);
%         vsurf2=vsurf2(1:xx,zval);
%         vsurf(i-1+length(u21(:,1)),1:xx) = vsurf2;
%     end
    for i = 2:length(u23(:,1))
        vsurf3=u23(i,:);
        vsurf3=reshape(vsurf3,xx+2,zz);
        vsurf3=vsurf3(1:xx,zval);
        vsurf(i-1+length(u21(:,1))+length(u22(:,1)),1:xx) = vsurf3;
    end
    
    vsurf = vsurf(1:(yy-1),:);

%     midptx = solidptsx/2 + 1;
%     shiftx = xx/2 + 1 - midptx;
%     vplot = [vsurf(:,end-shiftx+1:end),vsurf(:,1:end-shiftx+1)];
    vplot = vsurf;
    figure(2)
    clf
    xl=Re_t*[0:xx-1]/(xx)*2*pi/alp; 
    contour_levels = [-1:0.005:1];
    contour_levels = -1:1:1;
%      c1 = contourf(xl,Re_t*yv(yv<ylim/Re_t),vplot(yv<ylim/Re_t,:)/utau,'linestyle','none');
    c1 = contourf(xl,Re_t*yv(yv<ylim/Re_t),vplot(yv<ylim/Re_t,:)/utau,contour_levels,'linestyle','none');
%     c1 = contour(xl,yu(1:22),uplot(1:22,:),contour_levels,'k');
    set(gca,'Fontn','Times','FontSize',20,'linewidth',1);
    box on
    colorbar
    %shading flat
    shading interp
    axis equal;xlabel('x^+');ylabel('y^+'); 
    %view(0,90)
    colormap('parula')
        
    
    
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