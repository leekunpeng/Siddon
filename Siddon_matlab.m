%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%创建人：李昆鹏 leekunpeng@hotmail.com
%日  期：2021年3月31 最后修改：2021.4.1
%修改人：
%日  期：
%功能：<Fast calculation of the exact radiological path  for a
%threedimensiona-original> 论文实现，仅一个角度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;
clear all;
%%
% 参数设置
% 模体
load Shepp-Logan_3D_256.mat
% load img3.mat
%%
Phantom = I;
% get phantom size
[Nx, Ny, Nz] = size(Phantom);
% SOD、SDD、探测器大小、探测器间隔、模体体素间隔
Detector_Size = 400;
Proj_ = zeros(Detector_Size, Detector_Size);    % single proj result

SOD = 800;
Source_Coordinate = [0,floor(Detector_Size/2),floor(Detector_Size/2)];  %(x,y,z)
dx = 1;	dy = 1; dz = 1; % spacing
X1 = Source_Coordinate(1);  Y1 = Source_Coordinate(2);  Z1 = Source_Coordinate(3);  % point 1

% X，Y，Z网格位置
range_Nx = [0:1:Nx-1];  range_Ny = [0:1:Ny-1];  range_Nz = [0:1:Nz-1];
Xplane = zeros(1,Nx);   Yplane = zeros(1,Ny);   Zplane = zeros(1,Nz);

% 模体几何位置
% Xplane1 Yplane1 Zplane1 -> bottom left corner
Xplane(1) = floor(SOD/2);	Yplane(1) = floor(Detector_Size-Ny)/2;   Zplane(1) = floor(Detector_Size-Nz)/2;
% X_i,Y_j,Z_k   e.q(3)
Xplane = Xplane(1) + range_Nx*dx;
Yplane = Yplane(1) + range_Ny*dy;
Zplane = Zplane(1) + range_Nz*dz;

%%
for i = 1:1:Detector_Size
    for j = 1:1:Detector_Size
        % point 2 coordinate
        X2 = SOD; Y2 = i; Z2 = j; % bottom left corner
        
        % 计算αmin和αmax   e.q(4) 注意X2==X1，Y2==Y1，Z2==Z1
        if X2~=X1
            Alpha_x(1)=(Xplane(1)-X1)/(X2-X1);
            Alpha_x(Nx)=(Xplane(Nx)-X1)/(X2-X1);
        end
        
        if Y2~=Y1
            Alpha_y(1)=(Yplane(1)-Y1)/(Y2-Y1);
            Alpha_y(Ny)=(Yplane(Ny)-Y1)/(Y2-Y1);
        end
        
        if Z2~=Z1
            Alpha_z(1)=(Zplane(1)-Z1)/(Z2-Z1);
            Alpha_z(Nz)=(Zplane(Nz)-Z1)/(Z2-Z1);
        end
        
        % quantities e.q (5)
        Alpha_Min = max([0,min(Alpha_x(1),Alpha_x(Nx)), min(Alpha_y(1),Alpha_y(Ny)), min(Alpha_z(1),Alpha_z(Nz))]);
        Alpha_Max = min([1, max(Alpha_x(1),Alpha_x(Nx)), max(Alpha_y(1),Alpha_y(Ny)), max(Alpha_z(1),Alpha_z(Nz))]);
        
        % if alpha_max less than or equal to alpha_min, then the ray does
        % not intersect the CT array
        if Alpha_Max>Alpha_Min
            % 计算交线位置，{iMin，iMax}，{jMin，jMax}，{kMin，kMax}
            % e.q (6)
            if (X2-X1)>=0
                iMin = Nx - floor((Xplane(Nx) - Alpha_Min*(X2-X1) - X1)/dx);
                iMax = 1 + ceil((X1 + Alpha_Max*(X2-X1) - Xplane(1))/dx);
            elseif (X2-X1)<=0
                iMin = Nx - floor((Xplane(Nx) - Alpha_Max*(X2-X1) - X1)/dx);
                iMax = 1 + ceil((X1 + Alpha_Min*(X2-X1) - Xplane(1))/dx);
            end
            
            if (Y2-Y1)>=0
                jMin = Ny - floor((Yplane(Ny) - Alpha_Min*(Y2-Y1) - Y1)/dy);
                jMax = 1 + ceil((Y1 + Alpha_Max*(Y2-Y1) - Yplane(1))/dy);
            elseif (Y2-Y1)<=0
                jMin = Ny - floor((Yplane(Ny) - Alpha_Max*(Y2-Y1) - Y1)/dy);
                jMax = 1 + ceil((Y1 + Alpha_Min*(Y2-Y1) - Yplane(1))/dy);
            end
            
            if (Z2-Z1)>=0
                kMin = Nz - floor((Zplane(Nz) - Alpha_Min*(Z2-Z1) - Z1)/dz);
                kMax = 1 + ceil((Z1 + Alpha_Max*(Z2-Z1) - Zplane(1))/dz);
            elseif (Z2-Z1)<=0
                kMin = Nz - floor((Zplane(Nz) - Alpha_Max*(Z2-Z1) - Z1)/dz);
                kMax = 1 + ceil((Z1 + Alpha_Min*(Z2-Z1) - Zplane(1))/dz);
            end
            
            % for given range index {iMin，iMax}，{jMin，jMax} and {kMin，kMax}
            % the sets parametric values {aplha_x}{alpha_y}{alpha_z}
            % corresponding to the intersections of the ray with these
            % planes if (X2-X1)>0: {alpha_x} = min -> max; if (X2-X1)<0
            % {alpha_x} = max -> min;
            % 计算αx，αy，αz
            % 计算αx
            AlphaX=zeros(1, round(iMax-iMin+0.5));
            if (X2-X1)>0
                for ix=1:(iMax-iMin+1)
                    AlphaX(ix)=(Xplane(round(ix+iMin-1))-X1)/(X2-X1);
                end
            end
            if (X2-X1)<0
                for ix=1:(iMin-iMax+1)
                    AlphaX(ix)=(Xplane(round(iMax+1-ix))-X1)/(X2-X1);
                end
            end
            
            % 计算αy
            AlphaY=zeros(1, round(jMax-jMin+0.5));
            if (Y2-Y1)>0
                for iy=1:(jMax-jMin+1)
                    AlphaY(iy)=(Yplane(round(iy+jMin-1))-Y1)/(Y2-Y1);
                end
            end
            if (Y2-Y1)<0
                for iy=1:(jMax-jMin+1)
                    AlphaY(iy)=(Yplane(round(jMax+1-iy))-Y1)/(Y2-Y1);
                end
            end
            
            % 计算αz
            AlphaZ=zeros(1, round(kMax-kMin+0.5));
            if (Z2-Z1)>0
                for iz=1:(kMax-kMin+1)
                    AlphaZ(iz)=(Zplane(round(iz+kMin-1))-Z1)/(Z2-Z1);
                end
            end
            if (Z2-Z1)<0
                for iz=1:(kMax-kMin+1)
                    AlphaZ(iz)=(Zplane(round(kMax+1-iz))-Z1)/(Z2-Z1);
                end
            end
            
            %merge AlphaX,AlphaY AlphaZ并按照升序合并成一个Alpha数组
            %Alpha里面还应该包含alphamin，和alphamax
            Alpha = zeros(1,(iMax-iMin+1)+(jMax-jMin+1)+(kMax-kMin+1)+2);   % e.q (9)
            Alpha = [Alpha_Min, AlphaX, AlphaY, AlphaZ, Alpha_Max];    % e.q (8)
            Alpha = sort(Alpha);	%升序排列
            Alen = size(Alpha,2);
            
            d12 = sqrt((X2-X1).*(X2-X1)+(Y2-Y1).*(Y2-Y1)+(Z2-Z1).*(Z2-Z1));    % eq.11
            % 计算交线长度
            d = 0;
            len = zeros(1,Alen);
            for m = 2:Alen
                % e.q (12) (13)
                AlphaMid = (Alpha(m)+Alpha(m-1))/2;
                i_m = 1 + floor((X1 + AlphaMid*(X2-X1) - Xplane(1))/dx);
                j_m = 1 + floor((Y1 + AlphaMid*(Y2-Y1) - Yplane(1))/dy);
                k_m = 1 + floor((Z1 + AlphaMid*(Z2-Z1) - Zplane(1))/dz);
                i_m = max(1,i_m); i_m = min(i_m,Nx);
                j_m = max(1,j_m); j_m = min(j_m,Ny);
                k_m = max(1,k_m); k_m = min(k_m,Nz);
                d = d + d12*(Alpha(m)-Alpha(m-1))*Phantom(i_m,j_m,k_m);
            end
            Proj_(i,j) = d;
        end
    end
end

%%
figure;imshow(Proj_,[]);imcontrast;