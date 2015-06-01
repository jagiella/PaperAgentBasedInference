X = [1:3; 4:6; 7:9];
X'
X(:)

w1=[];
w2=[];
for i=1:3
    for j=1:3
        if( i+1 <= 3)
            w1(end+1) = X(i,j) - 1;% - X(i,j);
        end
        if( j+1 <= 3)
            w2(end+1) = X(i,j) - 1;% - X(i,j);
        end
    end
end
w1, w2
