close
clear all
clc

%% READ MODEL INFORMATION
fid = fopen('model_info.txt','r');
i = 1;
while ~feof(fid)
  b = fscanf(fid,'%c',1);
  if b==','
    a(i) = fscanf(fid,'%f',1);
    i = i +1;
  end
end
fclose(fid);

nNode = a(1);
nElem = a(2);
nDOF = a(3);
nDOF_srf = a(4);
dim_x = a(5);
dim_y = a(6);
h = a(7);
nNodeX = a(8)-2; % left and right nodes are fixed, and will not be displayed
nNodeY = a(9)-1; % bottom node is fixed, and will not be displayed
nElemX = a(10);
nElemY = a(11);

%% READ MODEL
node2xy = importdata('node2xy.txt');
x = node2xy(:,1);
y = node2xy(:,2);
DOFx2node = importdata('DOFx2node.txt')+1; % original data is in C index

X = zeros(nNodeY,nNodeX);
Y = X;
for i = 1 : nNodeX
    for j = 1 : nNodeY
        X(j,i) = x(nNodeY*(i-1)+j);
        Y(j,i) = y(nNodeY*(i-1)+j);
    end
end

%% READ SIMULATION INFORMATION
fid = fopen('simulation_info.txt','r');
i = 1;
while ~feof(fid)
  b = fscanf(fid,'%c',1);
  if b==','
    a(i) = fscanf(fid,'%f',1);
    i = i +1;
  end
end
fclose(fid);

nTstep = a(1);
dt = a(2);

%% READ DATA
fid = fopen('ux_reg.dat','rb');
ux_reg = fread(fid,[nDOF/2,nTstep],'double');
fclose(fid);
ux_max = max(max(abs(ux_reg)));


%% PLOT SNAPSHOTS
status = exist('snapshots','dir');
if status==0
    mkdir('snapshots');
end
figure('Position',[263,164,560,420])
colormap jet
Z = zeros(size(X));
for k = 1 : nTstep
    for i = 1 : nNodeX
        for j = 1 : nNodeY
            Z(j,i) = ux_reg(nNodeY*(i-1)+j,k);
        end
    end
    if mod(k,5)==0
        fprintf('time step %4i\n',k);
    end
    pcolor(X,Y,Z);
    title(sprintf('t=%7.5f/%7.5f',dt*(k-1),(nTstep-1)*dt))
    xlabel('x[m]')
    ylabel('y[m]')
    shading flat
    caxis([-ux_max,ux_max])
    axis image
    axis([h,(nElemX-1)*h,h,nElemY*h])
    print(sprintf('snapshots/step%04i',k),'-dpng')
end

% %% PLOT
% [X,Y] = meshgrid(x,y);
% ux = ux_reg(:,1);
% UX = griddata(x,y,ux,X,Y);
% pcolor(X,Y,UX);shading flat
% colormap jet
% caxis([-ux_max,ux_max]);