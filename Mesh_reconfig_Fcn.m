%   网格重新绘制 % % 
%   主要是解决网格链接顺序不一致的问题 % 
function remesh = Mesh_reconfig_Fcn(mesh_input,flag)

t = mesh_input ;  
% load t_tmp ;
if flag< 0.1
    tmp = t(1,1) ;
    t(1,1) = t(1,2) ; 
    t(1,2) = tmp; 
end

t_0 = t(1,:) ; 
t_1 = t(1,:) ; 
t_2 = t(2:end,:) ; 
iter = 1 ; 
% threshold = 3 ; && iter<threshold 
clear  t ; 
while( ~isempty( t_2 ) ) 
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
        index =  find(  match_flag1_1 + match_flag1_2 + match_flag1_3 +  match_flag2_1 + match_flag2_2 + match_flag2_3  +  match_flag3_1 + match_flag3_2 + match_flag3_3  == 2 ) ; 
        
        
        for jj = 1:length( index )
            line1 = [ match_flag1_1( index(jj) ) ,  match_flag1_2( index(jj) ) ,  match_flag1_3( index(jj) ) ;...
                      match_flag2_1( index(jj) ) ,  match_flag2_2( index(jj) ) ,  match_flag2_3( index(jj) ) ;...
                      match_flag3_1( index(jj) ) ,  match_flag3_2( index(jj) ) ,  match_flag3_3( index(jj) ) ] ; 
            zero_line      =  find( sum(line1,2)==0 )     ; 
            
            none_zero_line =  setdiff([1,2,3],zero_line ) ; 
            
            match1   =  find( line1( none_zero_line(1) , :  ) == 1 ) ; 
            match2   =  find( line1( none_zero_line(2) , :  ) == 1 ) ; 
            % 
            clock_wise_flag1 = none_zero_line(2) - none_zero_line(1) ; % 
            if( clock_wise_flag1==-2 || clock_wise_flag1==1 )
                clock_wise_flag1 =  1 ; 
            else
                clock_wise_flag1 = -1 ; 
            end
            % 
            clock_wise_flag2  = match2  - match1  ;
            if( clock_wise_flag2==-2 || clock_wise_flag2==1 )
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
        index_all = ( [ index_all   ; index  ]  ) ; 
        % index_all = unique(index_all)  ;
    end
    
    t_2 = t_2_save ; %  恢复t2
    
    save t_1 ;
    save t_2 ; 
    
    %% 总的处理一次就可以 % 
    if( ~isempty(index_all) && ~isempty( t_2 )  )  
        for jj = 1 : numel( index_all ) %  reshuttle 
            disp('Ok') ;
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

remesh = t_0 ; % 重新覆盖原有的位置 %  
end