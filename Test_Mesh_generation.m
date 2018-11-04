clc
clear
close all
%% Points cloud :
clear all 
% fid = fopen('Liver_out.obj');
fid = fopen('OKOKK1.obj');
addpath( genpath(cd) ) ;
iter  = 0 ;  
iter2 = 0 ;  
%% 
while( ~feof(fid) )
    SS = fgets(fid) ;
    if( SS(1)=='v' &&  SS(2)==' ')
        % scanf(SS,'v %f//%f %f//%f %f//%f',a,b,c,d,e,f)
        iter = iter+1 ;
        C = textscan(SS,  'v %f %f %f' ) ;  % 可以分析这一行数据 保存在C这个cell中
        x(iter) = C{1} ;
        y(iter) = C{2} ;
        z(iter) = C{3} ;
    elseif ( SS(1)=='v' &&  SS(2)=='n')
        iter2 = iter+1 ;
        C = textscan(SS,  'vn %f %f %f' ) ;  % 可以分析这一行数据 保存在C这个cell中
        nx(iter2) = C{1} ;
        ny(iter2) = C{2} ;
        nz(iter2) = C{3} ;
    end
    %     // S = fscanf(fid,'%s');
end

fclose(fid);

p =  [ x' y' z'];
p(  isnan( x'+y'+z') ,   :  )  = [] ;  

figure ; plot(x , y ,'*') ; axis vis3d ; axis equal ; %% 
t = delaunay( x, y ) ; 

decenter =  p - repmat ( mean(p) , [ length(p) , 1  ] )  ; 
decenter_pq =  decenter'*decenter ; 
[u ,v ,d ] = svd( decenter_pq );   

xyz1 = decenter * u  ;
xyz2 = decenter * d  ;  

figure ; plot( xyz1(:,1) , xyz1(:,2) , '*' )
figure ; plot( xyz2(:,1) , xyz2(:,2) , '*' )


%% plot of the current point cloud
figure(1);
hold on
axis equal
title('Points Cloud','fontsize',14)
data = p ;
plot3(data(:,1),data(:,2),data(:,3),'g.')
%% Run  program
% [t] = MyCrust( data ) ; % 
% t=delaunay(x,y,z);
%% plot of the oyput triangulation
figure 
hold on
title('Output Triangulation','fontsize',14)
axis equal ;  
%% 
t_back = t ; 
t(t>length(data) )  = length(data) ;  %   
% t = Mesh_reconfig_Fcn( t , 1 ) ;%%  
trisurf ( t , data(:,1) , data(:,2),data(:,3),'facecolor','c','edgecolor','b' )  % plot della superficie trat tata    
%% %%  制作法向量 


vertex_normal_flag   = 1;
vertex_normal = zeros( 3, length( data ) ) ;
% faces_normal  = zeros( 3, length( t  ) ) ;
first_vertex  = p( t(:,1) ,: );
second_vertex = p( t(:,2) ,: );
third_vertex  = p( t(:,3) ,: );
line_v1_v2    = first_vertex  - second_vertex ;
line_v2_v3    = second_vertex - third_vertex  ;
clear first_vertex third_vertex second_vertex  ;
faces_normal = cross( line_v1_v2 , line_v2_v3 , 2 ) ; % 列的维度进行
faces_normal = faces_normal./repmat(sqrt(sum(faces_normal.^2,2)) , 1 ,size(faces_normal, 2 ) )  ;
if(vertex_normal_flag)
    for i=1:length( data )
        [a,b] = find( t == i) ;
        vertex_normal(:,i) =  mean( faces_normal(a,:) )';
    end
    vertex_normal = vertex_normal./repmat(sqrt(sum(vertex_normal.^2,1)) ,size(vertex_normal, 1 ) , 1  )  ;
end
%%
point_cloud =  p' ;   % 
face_triangle_index   = t' ;  %  
normal_triangle_index = t' ;  % 
faces_index = [  face_triangle_index(1,:) ;normal_triangle_index(1,:);... 
                 face_triangle_index(2,:) ;normal_triangle_index(2,:);...
                 face_triangle_index(3,:) ;normal_triangle_index(3,:)   ] ; % 
 
fid = fopen([ 'Liver_',num2str(301),'.obj'] , 'w' );  %
fprintf(fid,'vt %f %f\n', [0,0] );
fprintf(fid,'v  %f %f %f\n' , point_cloud ); 
fprintf(fid,'vn %f %f %f\n', vertex_normal );
fprintf(fid,'f %d/1/%d %d/1/%d %d/1/%d\n',  faces_index );
% fprintf(fid,'f %d %d %d\n',  faces_index );
fclose(fid);
%% 
