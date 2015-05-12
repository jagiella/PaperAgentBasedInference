function [substr, common1, common2] = diffSubStr( varargin)
%size(varargin{:})
%varargin
l = minStrLen( varargin{:});
i = 1;
while(i <= l && equalAtIndex( i, varargin{:}))
    i = i + 1;
end

j=0;
while( equalAtLastIndex(j, varargin{:}))
    j = j + 1;
end

for k=1:nargin
    substr{k} = varargin{k}(i:end-j);
end

common1 = varargin{k}(1:i-1);
common2 = varargin{k}(end-j+1:end);

end

function l = minStrLen( varargin)
    l = length( varargin{1});
    for i=2:nargin
        l = min(l, length( varargin{i}));
    end
end

function eq = equalAtIndex( idx, varargin)
    letter = varargin{1}(idx);
    for i=2:nargin-1
        varargin{i}(idx)
        if( letter ~= varargin{i}(idx))
            eq = false;
            return 
        end
    end
    eq = true;
end

function eq = equalAtLastIndex( idx, varargin)
%size(varargin)
%idx,
%disp(varargin{1}),
%varargin{1}(end-idx)
    letter = varargin{1}(end-idx);
    for i=2:nargin-1
        %varargin{i}(end-idx)
        if( letter ~= varargin{i}(end-idx))
            eq = false;
            return 
        end
    end
    eq = true;
end