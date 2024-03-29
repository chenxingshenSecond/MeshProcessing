clc
clear
close all
%% Points cloud :
% I = textscan ('skull_cc.obj');   % fread
filename = 'horse.txt' ; 
fid = fopen(  filename  );

x = [] ; %% 
y = [] ; %% 
z = [] ; %% 

vp_ind = [];  vn_ind = []; 
while( ~feof(fid) )
    SS = fgets(fid) ;
    if( 1 )
       %% scanf(SS,'v %f//%f %f//%f %f//%f',a,b,c,d,e,f)
        C = textscan(SS,  '%d %f %f %f' ) ;  %  
        x = [ x , C{2} ] ;
        y = [ y , C{3} ] ;
        z = [ z , C{4} ] ;
    end
    %     // S = fscanf(fid,'%s');
end

fclose(fid);

%% 
filename = 'horse759.ele' ; 
fid = fopen(  filename  );
t = [] ; 
while( ~feof(fid) )
    SS = fgets(fid) ;
    if( 1 )
        C = textscan(SS,  '%d %d %d %d %d' ) ;  %  
        t = [ t ; [ C{2} C{3}  C{4}  C{5}  ] ];
    end
    %     // S = fscanf(fid,'%s');
end
t = t + 1; 
fclose(fid);



p = [x' y' z'] ; 
%% 
figure(1);
hold on
axis equal
title('Points Cloud','fontsize',14)
data = p ; camproj('perspective') ; 
plot3( data(:,1),data(:,2),data(:,3),'g.') %%  

% t = delaunay(x,y,z);
%% plot of the oyput triangulation
figure(4)
hold on
title('Output Triangulation','fontsize',14)
light ; camproj('perspective') ; 
lighting phong
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'c' , 'edgecolor' , 'b' ) ;  %plot della superficie trattata
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'y' , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'b' , 'edgecolor' , 'b' ) ;  %plot della superficie trattata
trisurf( t(:,1:3) , data(:,1) , data(:,2),data(:,3) , 'facecolor' , [0.8,0.8,0.8] , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
trisurf( t(:,2:4) , data(:,1) , data(:,2),data(:,3) , 'facecolor' , [0.8,0.0,0.0] , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
trisurf( t(:,[1 2 4]) , data(:,1) , data(:,2),data(:,3) , 'facecolor' , [0.0,0.8,0.0] , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
trisurf( t(:,[1 3 4]) , data(:,1) , data(:,2),data(:,3) , 'facecolor' , [0.0,0.0,0.8] , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
axis equal ;
axis vis3d ;
daspect([1 1 1]) ;
axis on; camproj('perspective') ; 
xlabel('x')
ylabel('y')
%%
% %% plot of the oyput triangulation
% fid = fopen([ filename ,num2str(300),'.obj'] , 'w' );  % % 
% fprintf(fid , '## Generated by Matlab \n \n'  ) ;
% fprintf(fid , [ '##  time = ' , datestr(now,0),'  \n \n']  ) ;
% 
% fprintf(fid,'vt %f %f\n', [0,0] );
% fprintf(fid,'v  %f %f %f\n' , point_cloud ); 
% fprintf(fid,'vn %f %f %f\n', vertex_normal );
% fprintf(fid,'f %d/1/%d %d/1/%d %d/1/%d\n',  faces_index );
% % fprintf(fid,'f %d %d %d\n',  faces_index );
% fclose(fid);
% 