clc
clear
close all
%% Points cloud :
% I = textscan ('skull_cc.obj');   % fread

fid = fopen('obj\skull_cc.obj');
iter  = 0 ;  
iter2 = 0 ;  
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

p=[x' y' z'];
figure(1);
hold on
axis equal
title('Points Cloud','fontsize',14)
data = p ;
plot3(data(:,1),data(:,2),data(:,3),'g.')
%% Run  program
[t] = MyCrust ( data ); 
%% Mesh_reconfig 
t_original = t ;  
% load t_tmp ;
t_0 = t(1,:) ; 
t_1 = t(1,:) ; 
t_2 = t(2:end,:) ; 
iter = 1 ; 
threshold = 260 ; 
clear  t ; 
while( ~isempty( t_2 ) && iter<= threshold) 
    disp(iter) ;
    index_all        = [] ;
    shuffle_flag_all = [] ; 
    t_2_save = t_2 ; % 保存版本 
    for i = 1 : size(t_1,1) % 顺序不绝对 
        %   match_flag1 = find( ( t_2(:,1)==t_1(i,1) ) + ( t_2(:,2)==t_1(i,2) ) + ( t_2(:,3)==t_1(i,3) ) ==2 ) ;  %
        %   match_flag2 = find( ( t_2(:,1)==t_1(i,1) ) + ( t_2(:,3)==t_1(i,2) ) + ( t_2(:,2)==t_1(i,3) ) ==2 ) ;  %
        %   match_flag3 = find( ( t_2(:,2)==t_1(i,1) ) + ( t_2(:,3)==t_1(i,2) ) + ( t_2(:,1)==t_1(i,3) ) ==2 ) ;  %
        %   match_flag4 = find( ( t_2(:,2)==t_1(i,1) ) + ( t_2(:,1)==t_1(i,2) ) + ( t_2(:,3)==t_1(i,3) ) ==2 ) ;  %
        %   match_flag5 = find( ( t_2(:,3)==t_1(i,1) ) + ( t_2(:,1)==t_1(i,2) ) + ( t_2(:,2)==t_1(i,3) ) ==2 ) ;  %
        %   match_flag6 = find( ( t_2(:,3)==t_1(i,1) ) + ( t_2(:,2)==t_1(i,2) ) + ( t_2(:,1)==t_1(i,3) ) ==2 ) ;  %
        
        match_flag1_1  = ( t_2(:,1)==t_1(i,1) ) ;  
        match_flag1_2  = ( t_2(:,2)==t_1(i,1) ) ;  
        match_flag1_3  = ( t_2(:,3)==t_1(i,1) ) ;  
        
        match_flag2_1  = ( t_2(:,1)==t_1(i,2) ) ;  
        match_flag2_2  = ( t_2(:,2)==t_1(i,2) ) ;  
        match_flag2_3  = ( t_2(:,3)==t_1(i,2) ) ;  
        
        match_flag3_1  = ( t_2(:,1)==t_1(i,3) ) ;  
        match_flag3_2  = ( t_2(:,2)==t_1(i,3) ) ;  
        match_flag3_3  = ( t_2(:,3)==t_1(i,3) ) ;  
        
        % index_matrix  =  (match_flag1==2) | (match_flag2==2) | (match_flag3==2) | (match_flag4==2) | (match_flag5==2) | (match_flag6==2)  ;  
        % clear match_flag1  match_flag2 match_flag3   match_flag4   match_flag5 match_flag6 ;
        index =  find(  abs(match_flag1_1 + match_flag1_2 + match_flag1_3 +  match_flag2_1 + match_flag2_2 + match_flag2_3  +  match_flag3_1 + match_flag3_2 + match_flag3_3 - 2 )< 0.1 ) ; 
        % ccc 
        for jj = 1:length( index )
            line1 = [ match_flag1_1( index(jj) ) ,  match_flag1_2( index(jj) ) ,  match_flag1_3( index(jj) ) ;...
                      match_flag2_1( index(jj) ) ,  match_flag2_2( index(jj) ) ,  match_flag2_3( index(jj) ) ;...
                      match_flag3_1( index(jj) ) ,  match_flag3_2( index(jj) ) ,  match_flag3_3( index(jj) ) ] ; 
            zero_line      =  find( abs( sum(line1,2)-0 )<0.1  )     ; 
            
            none_zero_line =  setdiff([1,2,3],zero_line ) ; 
            
            match1   =  find( abs( line1( none_zero_line(1) , :  ) - 1 )<0.1 ) ; %
            match2   =  find( abs( line1( none_zero_line(2) , :  ) - 1 )<0.1 ) ; %
            % 
            clock_wise_flag1 = none_zero_line(2) - none_zero_line(1) ; % 
            if( clock_wise_flag1==-2 || clock_wise_flag1==1 )
                clock_wise_flag1 =  1 ; 
            else
                clock_wise_flag1 = -1 ; 
            end
            % 
            clock_wise_flag2  = match2  - match1  ;
            if( abs(clock_wise_flag2+2)<0.1 || abs(clock_wise_flag2-1)<0.1 )
                clock_wise_flag2 =  1 ; 
            else
                clock_wise_flag2 = -1 ; 
            end
            % line1( zero_line , : ) = [1 , 1 , 1 ] - line1( none_zero_line(1) , :  ) - line1( none_zero_line(2) , :  )  ; 
            
            if ( clock_wise_flag2 * clock_wise_flag1 < 0 )   
                shuffle_flag_all = [ shuffle_flag_all , 1 ] ; % 
            else
                shuffle_flag_all = [ shuffle_flag_all , -1 ] ;
            end
        end
        t_2(index,:) =  0  ; % 防止再次被处理 
        %  unique ([ match_flag1 ; match_flag2; match_flag3 ; match_flag4 ; match_flag5 ; match_flag6] ) ;%  find(index_matrix)   ; 
        index_all =   [ index_all   ; index  ] ; 
        % index_all = unique(index_all)  ;
    end 
    
    t_2 = t_2_save  ; %  恢复t2
    %% 总的处理一次就可以 % 
    if( ~isempty(index_all) && ~isempty( t_2 )  )  
        for jj = 1 : numel( index_all ) %  reshuttle 
            % disp('Ok') ;
            if( shuffle_flag_all(jj) < 0  )  % || sum( index_all(jj)==match_flag4 )  || sum( index_all(jj)==match_flag6 ) 
                tmp  = t_2( index_all(jj) , 1 )  ;
                t_2( index_all(jj) , 1 ) = t_2( index_all(jj) , 2)  ; %
                t_2( index_all(jj) , 2 ) = tmp  ; %
            end
        end
        t_0 = [ t_0 ; t_2(index_all,:) ] ;
        t_1 =  t_2(index_all,:)  ;
        t_2 (index_all ,:) = [] ;
    elseif isempty( index_all )
        break ; 
    end
    iter = iter + 1; 
end

t = t_0 ; % 重新覆盖原有的位置 % 

% t=delaunay(x',y',z');
%% plot of the oyput triangulation
figure(4)
hold on
title('Output Triangulation','fontsize',14)
light
lighting phong
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'c' , 'edgecolor' , 'b' ) ;  %plot della superficie trattata
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'y' , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'b' , 'edgecolor' , 'b' ) ;  %plot della superficie trattata
trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , [0.8,0.8,0.8] , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
axis equal ;
axis vis3d ;
daspect([1 1 1]) ;
axis on;
xlabel('x')
ylabel('y')
%% 
figure(3)
hold on
title('Output Triangulation','fontsize',14)
light
lighting phong ; 
% trisurf( t_2 , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'c' , 'edgecolor' , 'b' ) ;  %plot della superficie trattata
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'y' , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
% trisurf( t , data(:,1) , data(:,2),data(:,3) , 'facecolor' , 'b' , 'edgecolor' , 'b' ) ;  %plot della superficie trattata
trisurf( t_2 , data(:,1) , data(:,2),data(:,3) , 'facecolor' , [0.8,0.8,0.8] , 'edgecolor' , 'none' ) ;  %plot della superficie trattata
axis equal ;
axis vis3d ;
daspect([1 1 1]) ;
axis on;
xlabel('x')
ylabel('y')


%% 法向量计算
vertex_normal_flag = 1;
vertex_normal = zeros( 3, length( data ) ) ;
% faces_normal  = zeros( 3, length( t  ) ) ;
first_vertex  = p( t(:,1) ,: );
second_vertex = p( t(:,2) ,: );
third_vertex  = p( t(:,3) ,: );
line_v1_v2    = first_vertex  - second_vertex ;
line_v2_v3    = second_vertex - third_vertex  ;
clear first_vertex third_vertex second_vertex  ;
faces_normal = cross( line_v1_v2 , line_v2_v3 , 2 ) ; % 列的维度进行  % 
faces_normal = faces_normal./repmat(sqrt(sum(faces_normal.^2,2)) , 1 ,size(faces_normal, 2 ) )  ; % 

if(vertex_normal_flag)  %%
    for i = 1 : length( data )    %%
        [a , b] = find( t == i) ; %%
        vertex_normal(:,i) =  mean( faces_normal(a,:) )' ; %%
    end %%
    vertex_normal = vertex_normal ./ repmat(sqrt(sum(vertex_normal.^2,1)) ,size(vertex_normal, 1 ) , 1  )  ;
end

%{ 　sprint 意思是string printf ； fprintf 是file printf  %}
point_cloud = [x;y;z] ; 
face_triangle_index   = t' ; 
normal_triangle_index = t' ; 
normal_orignal        = [ nx ;ny ;nz ] ;

faces_index = [  face_triangle_index(1,:) ;normal_triangle_index(1,:);... 
                 face_triangle_index(2,:) ;normal_triangle_index(2,:);...
                 face_triangle_index(3,:) ;normal_triangle_index(3,:)] ; % ; t // 
% 
%
fid = fopen([ 'skull_cc_with_normal_',num2str(300),'.obj'] , 'w' );  %
fprintf(fid,'vt %f %f\n', [0,0] );
fprintf(fid,'v  %f %f %f\n' , point_cloud ); 
fprintf(fid,'vn %f %f %f\n', vertex_normal );
fprintf(fid,'f %d/1/%d %d/1/%d %d/1/%d\n',  faces_index );
% fprintf(fid,'f %d %d %d\n',  faces_index );
fclose(fid);
% 