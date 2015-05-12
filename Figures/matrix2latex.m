function matrix2latex( filename, A, L1, L2, varargin)
    %formatSpec = 'X is %4.2f meters or %8.3f mm\n';
    
    fileID = fopen( filename,'w');
    
    %% header
    fprintf( fileID, '\\begin{tabular}{');
    fprintf( fileID, 'r| ');
    for i=1:size(A,2)
        fprintf( fileID, 'r ');
    end 
    fprintf( fileID, '}\n');
    for i=1:size(A,2)
        if( iscell(L1))
            fprintf( fileID, '&%s ', L1{i})
        else
            fprintf( fileID, '&%s ', num2str(L1(i)))
        end
    end
    fprintf( fileID, '\\\\\n\\hline \n');
    
    for i=1:size(A,1)
        if( iscell(L2))
            fprintf( fileID, L2{i});
        else
            fprintf( fileID, num2str(L2(i)));
        end
        
        if ~isempty(varargin)
            for j=1:size(A,2)
                if A(i,j) == 1
                    fprintf( fileID, '&++ ');
                else
                    if A(i,j) >= 0.1
                        fprintf( fileID, '&+ ');
                    else
                        if A(i,j) >= 0.05
                            fprintf( fileID, '&- ');
                        else
                            fprintf( fileID, '&- - ');
                        end
                    end
                end 
            end
        else
            fprintf( fileID, '&%6.3f ', A(i,:))
        end
        
        fprintf( fileID, '\\\\\n');
    end
    
    %% footer
    fprintf( fileID, '\\end{tabular}\n');
    
    fclose(fileID);
end