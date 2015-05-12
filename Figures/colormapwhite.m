function cmap = colormapwhite( n, varargin)

c1 = [1 0 0];
c2 = [0 0 1];
if nargin >= 2
    c1 = varargin{1};
end
if nargin >= 3
    c2 = varargin{2};
end

cmap=[];

a=floor(n/2);
b=a+1;

for i=1:a % excludes white
    alpha = (i-1)/a;
    cmap(i,:) = alpha * [1 1 1] + (1-alpha) * c1;
end
for i=b:n % includes white
    alpha = (i-b)/(n-b);
    cmap(i,:) = (1-alpha) * [1 1 1] + alpha * c2;
end

end