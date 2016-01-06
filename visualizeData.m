
folder = 'log/';
%format = '-djpeg';
format = '-deps';
step1 = 2;
step2 = 2;
scale = 1.0;

fileID = fopen([folder 'parameters.txt'],'r');
Re = fscanf(fileID,'%f', 1);
Xmax = fscanf(fileID,'%f', 1);
rkap = fscanf(fileID,'%f', 1);
tors = fscanf(fileID,'%f', 1);
epsr = fscanf(fileID,'%f', 1);
Im = fscanf(fileID,'%f', 1);
Jm = fscanf(fileID,'%f', 1);
Km = fscanf(fileID,'%f', 1);
fclose(fileID); 

fileID = fopen([folder 'rth.txt'],'r');
if fileID == -1,
    rth = linspace(0,1,Jm);
else
    rth = fscanf(fileID,'%f', Jm);
    fclose(fileID);
end

fileID = fopen([folder 'p.txt'],'r');
p_ = fscanf(fileID,'%f', Im*Km*Jm);
p = reshape(p_,Km,Jm,Im);
fclose(fileID);

fileID = fopen([folder 'u.txt'],'r');
u_ = fscanf(fileID,'%f', (Im+1)*Jm*Km);
uu = reshape(u_,Km,Jm,Im+1);
u = 0.5*(uu(:,:,1:end-1) + uu(:,:,2:end));
fclose(fileID);

fileID = fopen([folder 'v.txt'],'r');
v_ = fscanf(fileID,'%f', Im*(Jm+1)*Km);
vv = reshape(v_,Km,Jm+1,Im);
v = 0.5*(vv(:,1:end-1,:) + vv(:,2:end,:));
v(:,1,:) = vv(:,2,:);
fclose(fileID);

fileID = fopen([folder 'w.txt'],'r');
w_ = fscanf(fileID,'%f', Im*Jm*(Km+1));
ww = reshape(w_,Km+1,Jm,Im);
w = 0.5*(ww(1:end-1,:,:) + ww(2:end,:,:));
fclose(fileID); 

U = squeeze(u(:,:,1));
V = squeeze(v(:,:,1));
W = squeeze(w(:,:,1));


fprintf('Data: max(global flow) = %f, max(secondary flow) = %f, max(secondary flow)/max(global flow) = %f\n', max(max(sqrt(U.^2 + V.^2 + W.^2))), max(max(sqrt(V.^2 + W.^2))), max(max(sqrt(V.^2 + W.^2)))/max(max(sqrt(U.^2 + V.^2 + W.^2))));


hr = 1./Jm;
ht = 2*pi./Km;
%[R TH] = meshgrid(0.5*hr:hr:(1-0.5*hr), 0.5*ht:ht:(2*pi-0.5*ht) );
[R TH] = meshgrid(rth, 0.5*ht:ht:(2*pi-0.5*ht) );
[X Y] = pol2cart(TH, R);

vel_th = W;
vel_r = V;
[v_pol, w_pol] = pol2cartLoc(vel_th,vel_r, TH);

% figure;
% subplot(1,2,1);
% I = v_pol;
% I_norm = (I-min(min(I)))./max(max(I-min(min(I))));
% warp(X, Y, zeros(size(X)), I_norm);
% view(2); axis square;
% 
% subplot(1,2,2);
% I = w_pol;
% I_norm = (I-min(min(I)))./max(max(I-min(min(I))));
% warp(X, Y, zeros(size(X)), I_norm);
% view(2); axis square;


fig = figure(1);  
quiver(X(1:step1:end,1:step2:end), Y(1:step1:end,1:step2:end), v_pol(1:step1:end,1:step2:end), w_pol(1:step1:end,1:step2:end), scale); 
axis square;% axis off;
print(fig,[folder 'secondaryFlowVector'],format);
% hold on;
% surf(X,Y,sqrt(v_pol.^2 + w_pol.^2),'FaceAlpha',0.6);
% shading flat;

%figure;
% lgmax = (max(max(sqrt(squeeze(u(:,:,1)).^2 + v_pol.^2 + w_pol.^2))));
% lgmin = (min(min(sqrt(squeeze(u(:,:,1)).^2 + v_pol.^2 + w_pol.^2))));
% lgsp = linspace(lgmin,lgmax,20).^3/lgmax^2;
% contour(X,Y,sqrt(v_pol.^2 + w_pol.^2),lgsp);
% colorbar;

fig = figure(2);
secondaryFlow = sqrt(V.^2 + W.^2);
contour([X' X(1,:)']', [Y' Y(1,:)']', [secondaryFlow' secondaryFlow(1,:)']');
axis square;% axis off;
colorbar;
print(fig,[folder 'secondaryFlowContour'],format);

fig = figure(3);
globalFlow = sqrt(U.^2 + V.^2 + W.^2);
contour([X' X(1,:)']', [Y' Y(1,:)']', [globalFlow' globalFlow(1,:)']');
axis square;% axis off;
colorbar;
print(fig,[folder 'globalFlowContour'],format);

%fig = figure(4);
% startx = -1:0.1:1;
% starty = zeros(size(startx));
% %v_pol vel_r
% %w_pol vel_th
% stream = stream2(R,TH,v_pol,w_pol,startx,starty);
% for k=1:max(size(startx))
%     sreamR = stream{k}(:,1);
%     sreamTH = stream{k}(:,2);
%     [streamX streamY] = pol2cart(sreamTH, sreamR);
%     plot(streamX,streamY);
%     hold on;
% end
% hold off;
% 
% figure(5);
% P = squeeze(p(:,:,1));
% contour([X' X(1,:)']', [Y' Y(1,:)']', [P' P(1,:)']');
% axis square;
% colorbar;
% 
% figure(6);
% Vel = sqrt(squeeze(u(:,:,1)).^2 + squeeze(v(:,:,1)).^2 + squeeze(w(:,:,1)).^2);
% surf([X' X(1,:)']', [Y' Y(1,:)']', [Vel' Vel(1,:)']' );
% axis square;
