function [psf_trans_site] = get_psf_transcription_site(max_int_pos, img_fish)
% get_psf_transcription_site: Estimates the width of the psf of the
% transcrtipion site, using 2D Gaussian psf and least square fitting.
% SYNOPSIS:
%  [psf_trans_site] = get_psf_transcription_site(max_int_pos, img_fish)
% 
% PARAMETERS:
%     max_int_pos: Position transcription site
% 
%     img_fish: smFISH data 
% 
% OUTPUTS:
%   psf_trans_site: Estimated psf of transcription site
 

% Get ROI around trancription site
img_focus = img_fish(:,:,max_int_pos(3));
boxsize = 3;
box_trans = img_focus(round(max_int_pos(2))+1-boxsize:round(max_int_pos(2))+1+boxsize, round(max_int_pos(1))+1-boxsize:round(max_int_pos(1))+1+boxsize);

%% ---------User Input---------------------
MdataSize = boxsize*2 + 1; % Size of nxn data matrix
% parameters are: [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
x0 = [max(box_trans, [], 'all'),0,1,0,1,0]; %Inital guess parameters
xs = x0;% [2,2.2,7,3.4,4.5,+0.02*2*pi]; %centroid parameters
noise = 0; % noise in % of centroid peak value (x(1))
InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'

%% ---Generate centroid to be fitted--------------------------------------
[X,Y] = meshgrid(-boxsize:boxsize);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
[Xhr,Yhr] = meshgrid(linspace(-MdataSize/2,MdataSize/2,300)); % generate high res grid for plot
xdatahr = zeros(300,300,2);
xdatahr(:,:,1) = Xhr;
xdatahr(:,:,2) = Yhr;

Z = box_trans - min(box_trans, [], 'all');
%% --- Fit---------------------
options=optimset('MaxFunEvals',100000, 'TolFun', 1e-10, 'maxiter', 5000 );
x0 =x0(1:5);

lb = [0,-MdataSize/2,0,-MdataSize/2,0];
ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2];
[x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub, options);
x(6) = 0;

%% ---------Plot data-------------
figure(1)
C = del2(Z);
mesh(X,Y,Z,C) %plot data
hold on
surface(Xhr,Yhr,D2GaussFunctionRot(x,xdatahr),'EdgeColor','none') %plot fit
axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5 -noise noise+x(1)])
alpha(0.2)  
hold off

hf2 = figure(2);
set(hf2, 'Position', [20 20 950 900])
alpha(0)
subplot(4,4, [5,6,7,9,10,11,13,14,15])
imagesc(X(1,:),Y(:,1)',Z)
set(gca,'YDir','reverse')
colormap('parula')
string1 = ['       Amplitude','    X-Coordinate', '    X-Width','    Y-Coordinate','    Y-Width','     Angle'];
string3 = ['Fit      ',num2str(x(1), '% 100.3f'),'             ',num2str(x(2), '% 100.3f'),'         ',num2str(x(3), '% 100.3f'),'         ',num2str(x(4), '% 100.3f'),'        ',num2str(x(5), '% 100.3f'),'     ',num2str(x(6), '% 100.3f')];


% generate points along horizontal axis
m = -tan(x(6));% Point slope formula
b = (-m*x(2) + x(4));
xvh = -MdataSize/2:MdataSize/2;
yvh = xvh*m + b;
hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
% generate points along vertical axis
mrot = -m;
brot = (mrot*x(4) - x(2));
yvv = -MdataSize/2:MdataSize/2;
xvv = yvv*mrot - brot;
vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);
hold on
plot([xvh(1) xvh(size(xvh))] ,[yvh(1) yvh(size(yvh))],':k', 'linewidth', 3) 
plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],':k','linewidth', 3) 
xlabel('Position [pixel]','FontSize',12)
ylabel('Position [pixel]','FontSize',12)
hold off
axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5])

ymin = - noise * x(1);
ymax = x(1)*(1+noise);
xdatafit = linspace(-MdataSize/2-0.5,MdataSize/2+0.5,300);
hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
vdatafit = x(1)*exp(-(xdatafit-x(4)).^2/(2*x(5)^2));
subplot(4,4, [1:3])
xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
plot(xposh+0.5,hPoints,'kx',xdatafit,hdatafit,':k' , 'MarkerSize',10, 'linewidth', 2)
%xlabel('Position [pixel]')
ylabel('Intensity [# photons]','FontSize',12)
axis([-MdataSize/2-0.5 MdataSize/2+0.5 ymin*1.1 ymax*1.1])
subplot(4,4,[8,12,16])
xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
plot(vPoints,xposv+0.5,'kx',vdatafit,xdatafit,':k',  'MarkerSize',10, 'linewidth', 2)
%ylabel('Position [pixel]')
xlabel('Intensity [# photons]','FontSize',12)
axis([ymin*1.1 ymax*1.1 -MdataSize/2-0.5 MdataSize/2+0.5])
set(gca,'YDir','reverse')
figure(gcf) % bring current figure to front

psf_trans_site = sqrt(x(3)^2+x(5)^2);
end
