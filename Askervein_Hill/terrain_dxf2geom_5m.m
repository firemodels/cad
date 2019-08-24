% M. Vanella
% Aug 2019
% terrain_dxf2geom.m
%
% Read Askervein Vector Map in DXF format and make a terrain geometry
% with 5m horizontal resolution, write to askervein.geom ASCII file.

close all
clear all
clc

%% Parameters:
IAXIS = 1; JAXIS = 2; KAXIS = 3; MDIM = 3;
NOD1  = 1; NOD2  = 2; NOD3  = 3; MNOD = 3;

basedir='./';
dxf_file='Askervein_vertices.dxf';

basedirout='./';
fileout ='Askervein.geom';

GEOM_ID='terrain';
SURF_ID='terrain';

grdflg='stretched'; % 'uniform'


%% Load dxf:
disp(['1. Load Vector vertices from DXF file: ' basedir dxf_file ' ..'])

[c_Line,c_Poly,c_Cir,c_Arc,c_Poi] = f_LectDxf([basedir dxf_file]);
np=length(c_Poi);
XYZ=zeros(np,3);
% Define points:
for ip=1:np
    XYZ(ip,IAXIS:KAXIS) = c_Poi{ip,1};
end

figure
plot3(XYZ(:,IAXIS),XYZ(:,JAXIS),XYZ(:,KAXIS),'.b','MarkerSize',11)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
box on

disp(['   Done.'])

%% Location of HT, CP, RS points, From Ask83.pdf report page 7:
HT_x = 75383;
HT_y = 23737;

CP_x = 75678;
CP_y = 23465;

RS_x = 74300;
RS_y = 20980;

%% Make X,Y mesh, Interpolate Z and build VERTS and FACES arrays:
disp(' ')
disp(['2. Build VERTS and FACES arrays ..'])
xmin=73300;
xmax=79300;
ymin=20000;
ymax=26000;
dx  =5;
dy  =5;
if(strcmp(grdflg,'uniform'))
   npx=(xmax-xmin)/dx+1;
   npy=(ymax-ymin)/dy+1;
   x  = linspace(xmin,xmax,npx);
   y  = linspace(ymin,ymax,npy);
elseif(strcmp(grdflg,'stretched'))
   % X direction:
   xs =           HT_x; % Mesh block is defined from xl=xs+Lx(1) to xf=xs+Lx(2)
   Lx(1)= -(HT_x-xmin);
   Lx(2)=  (xmax-HT_x);

   % Low side number m of cells with small uniform size DXS:
   m(1)  =   200;
   DXS(1)=   -dx;
   m(2)  =   220;
   DXS(2)=-DXS(1);

   % Total number of cells:
   Nx(1) = 250;
   Nx(2) = 280;

   [npx,x]=mesh_stretch(xs,Lx,m,DXS,Nx);

   % Y direction:
   ys =           HT_y; % Mesh block is defined from xl=xs+Lx(1) to xf=xs+Lx(2)
   Ly(1)= -(HT_y-ymin);
   Ly(2)=  (ymax-HT_y);

   % Low side number m of cells with small uniform size DXS:
   m(1)  =   300;
   DYS(1)=   -dy;
   m(2)  =   180;
   DYS(2)=-DYS(1);

   % Total number of cells:
   Ny(1) = 360;
   Ny(2) = 230;

   [npy,y]=mesh_stretch(ys,Ly,m,DYS,Ny);

end

nverts = npx*npy;
VERTS=zeros(nverts,KAXIS);
ip=0;
for j=1:npy
    for i=1:npx
        ip = ip+1;
        VERTS(ip,IAXIS:JAXIS) = [x(i) y(j)];
    end
end
vq = griddata(XYZ(:,IAXIS),XYZ(:,JAXIS),XYZ(:,KAXIS),...
              VERTS(:,IAXIS),VERTS(:,JAXIS),'natural');
VERTS(:,KAXIS) = vq;

nfaces = 2*(npx-1)*(npy-1);
FACES=zeros(nfaces,NOD3+1);
SURF_IND = 1;
ifc=0;
for j=1:npy-1
    for i=1:npx-1
        ip1 = (j-1)*npx + i;
        ip2 = (j-1)*npx + i + 1;
        ip3 = (j-0)*npx + i;
        ip4 = (j-0)*npx + i + 1;
        ifc = ifc + 1;
        FACES(ifc,1:NOD3+1) = [ip1 ip2 ip4 SURF_IND];
        ifc = ifc + 1;
        FACES(ifc,1:NOD3+1) = [ip1 ip4 ip3 SURF_IND];
    end
end

disp( '   Done.')
disp(['   Number of points on Regular Grid: ' num2str(nverts)])
disp(['   Number of  faces on Regular Grid: ' num2str(nfaces)])

%% DEVC locations, added by RJM ...............
F = scatteredInterpolant(VERTS(:,IAXIS)-HT_x,VERTS(:,JAXIS)-HT_y,VERTS(:,KAXIS),'natural')
M = importdata('Askervein_devc_loc.csv',',',1)
xdevc = M.data(:,4)-HT_x;
ydevc = M.data(:,5)-(HT_y+8e5);
zdevc_agl = M.data(:,3);
zdevc = F(xdevc,ydevc)+zdevc_agl;
fid = fopen('Askervein_devc_loc.fds','wt');
for i=1:length(xdevc)
  if zdevc(i)>0
    % devc = ['&DEVC ID=''',M.textdata{i+1,1},''', XYZ=',num2str(xdevc(i)),',',num2str(ydevc(i)),',',num2str(zdevc(i)),', QUANTITY=''VELOCITY''/'];
    fprintf( fid, '&DEVC ID=''%s'', XYZ=%2f,%2f,%2f, QUANTITY=''VELOCITY''/\n', M.textdata{i+1,1}, xdevc(i),ydevc(i),zdevc(i) );
  end
end
fclose(fid);
% return
%% end RJM ....................................

figure
hold on
trimesh(FACES(:,NOD1:NOD3),VERTS(:,IAXIS)-HT_x,VERTS(:,JAXIS)-HT_y,VERTS(:,KAXIS))
plot3(HT_x-HT_x,HT_y-HT_y,130,'.k','MarkerSize',13)
text(HT_x-HT_x+20,HT_y-HT_y+20,150,'HT','FontSize',16)
plot3(CP_x-HT_x,CP_y-HT_y,130,'.k','MarkerSize',14)
text(CP_x-HT_x,CP_y-HT_y,150,'CP','FontSize',16)
plot3(RS_x-HT_x,RS_y-HT_y,70,'.k','MarkerSize',14)
text(RS_x-HT_x,RS_y-HT_y,80,'RS','FontSize',16)

colorbar
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
box on

%% Write geom file:
disp(' ')
disp(['3. Writing geom file: ' basedirout fileout ' ..'])
tstart=tic;
[fid]=fopen([basedirout fileout],'w');

% Write &GEOM namelist
[wid]=fprintf(fid,'&GEOM ID=''%s'', SURF_ID=''%s'', IS_TERRAIN=T, \n',GEOM_ID,SURF_ID);

% Vertices:
[wid]=fprintf(fid,'VERTS=\n');
for ip=1:nverts
    [wid]=fprintf(fid,' %14.8f,   %14.8f,   %14.8f,\n',...
          [VERTS(ip,IAXIS)-HT_x VERTS(ip,JAXIS)-HT_y VERTS(ip,KAXIS)]);
end

% Tetrahedra:
% No volume elements.

% Surface Triangles:
[wid]=fprintf(fid,'FACES=\n');
for ifc=1:nfaces
    [wid]=fprintf(fid,' %6d,   %6d,   %6d, %6d\n',FACES(ifc,NOD1:NOD3+1));
end

[wid]=fprintf(fid,'/ \n');

% Close geom file
fclose(fid);

disp(['   Done. Time taken: ' num2str(toc(tstart)) ' sec.'])

disp('4. Askervein GEOM parameters: ')
disp(['   Position X,Y point HT : ' num2str(HT_x-HT_x,'%8.2f') ', ' num2str(HT_y-HT_y,'%8.2f') ' m.'])
disp(['   Position X,Y point CP : ' num2str(CP_x-HT_x,'%8.2f') ', ' num2str(CP_y-HT_y,'%8.2f') ' m.'])
disp(['   Position X,Y point RS : ' num2str(RS_x-HT_x,'%8.2f') ', ' num2str(RS_y-HT_y,'%8.2f') ' m.'])
disp(['   Terrain XMIN, and XMAX: ' num2str(xmin-HT_x,'%8.2f') ', ' num2str(xmax-HT_x,'%8.2f') ' m.'])
disp(['   Terrain YMIN, and YMAX: ' num2str(ymin-HT_y,'%8.2f') ', ' num2str(ymax-HT_y,'%8.2f') ' m.'])
disp(['   Terrain ZMIN, and ZMAX: ' num2str(min(VERTS(:,KAXIS)),'%8.2f') ', ' num2str(max(VERTS(:,KAXIS)),'%8.2f') ' m.'])

[valx,iloc]=min(abs(x-HT_x));
[valy,jloc]=min(abs(y-HT_y));
ivert = (jloc-1)*npx + iloc;

disp(['   Terrain Height at HT  : ' num2str(VERTS(ivert,KAXIS),'%8.2f') ' m.'])

%% Now suggest some meshes:
Lx = xmax-xmin;
Ly = ymax-ymin;

Target_dx = 4; % either 4, 8, 16 m.
Target_dy = Target_dx;

if(Target_dx == 16)
   N_mesh_x = 15;
elseif(Target_dx == 8)
   N_mesh_x = 15;
elseif(Target_dx == 4)
   N_mesh_x = 30;
end
N_mesh_y = 12;

disp(' ')
disp('Suggested Meshes:')
disp(['Single mesh line, DX,DY= ' num2str(Target_dx) ' m, ' num2str(N_mesh_x*N_mesh_y) ' meshes:'])

DX   = Lx / N_mesh_x;
X_lo = xmin-HT_x; X_hi=X_lo+DX;

DY   = Ly / N_mesh_y;
Y_lo = ymin-HT_y; Y_hi=Y_lo+DY;

disp(' ')
disp(['&MULT ID=''Mult_' num2str(Target_dx) 'm'', DX=' num2str(DX,'%8.2f') ...
       ', DY=' num2str(DY,'%8.2f') ', I_UPPER=' num2str(N_mesh_x-1) ...
       ', J_UPPER=' num2str(N_mesh_y-1) '/'])
if(Target_dx == 16)
disp(['&MESH IJK=' num2str(ceil(DX/Target_dx)) ',' ...
                   num2str(ceil(DY/Target_dy)) ',38, ' ...
      'XB=' num2str(X_lo,'%8.2f') ',' num2str(X_hi,'%8.2f') ...
      ','   num2str(Y_lo,'%8.2f') ',' num2str(Y_hi,'%8.2f') ...
      ',-20.25,979.75, MULT_ID=''Mult_16m'', TRNZ_ID=''MY TRANSFORM''/'])
elseif(Target_dx == 8)
disp(['&MESH IJK=' num2str(ceil(DX/Target_dx)) ',' ...
                   num2str(ceil(DY/Target_dy)) ',58, ' ...
      'XB=' num2str(X_lo,'%8.2f') ',' num2str(X_hi,'%8.2f') ...
      ','   num2str(Y_lo,'%8.2f') ',' num2str(Y_hi,'%8.2f') ...
      ',-20.25,979.75, MULT_ID=''Mult_8m'', TRNZ_ID=''MY TRANSFORM''/'])
elseif(Target_dx == 4)
disp(['&MESH IJK=' num2str(ceil(DX/Target_dx)) ',' ...
                   num2str(ceil(DY/Target_dy)) ',90, ' ...
      'XB=' num2str(X_lo,'%8.2f') ',' num2str(X_hi,'%8.2f') ...
      ','   num2str(Y_lo,'%8.2f') ',' num2str(Y_hi,'%8.2f') ...
      ',-20.25,979.75, MULT_ID=''Mult_4m'', TRNZ_ID=''MY TRANSFORM''/'])
end
disp(' ')


return