function [y p] = MMD_DTW(x)
% MMD-DTW for multiple multidimensional signals
% x -> cell array of multiple n-D signals.
% p -> the optimal path.
% y -> matrix of the normalized signals in direction of column.
%
% See: P. SANGUANSAT - "Multiple Multidimensional Sequence Alignment Using
% Generalized Dynamic Time Warping" - WSEAS Transactions on Mathematics
%
% NOTE: The computed distance distance, has no meaning for comparing the
% multiple sequences. Each sequence must be previously rearranged as the
% obtained path before calculating the distance, which can computed by any
% classical method such as Euclidean distance.

if ~exist('x','var')
    x{1} = [1 2 3;1 2 3];
    x{2} = [1 3 3 4;1 3 3 4 ];
    x{3} = [1 2 2 3 4;1 2 2 3 4];

    show = 1;
end

if ~exist('show','var')
    show = 0;
end

N = length(x);
% Initial e
sizee = zeros(1,N);
for i=1:N
    sizee(i) = size(x{i},2)+2;
end
e = realmax*ones(sizee);
e(1) = 0;

% Find accumulated distance
str = '';
for q=1:N
    str = [str sprintf('fori%d=2:size(x{%d},2)+1\n',q,q)];
end

str = [str 'cost = dist(['];
for q=1:N
    str = [str sprintf('x{%d}(:,i%d-1) ',q,q)];
end
str = [str ']);'];
str = [str 'block = e('];

for q=1:N
    str = [str sprintf('i%d-1:i%d,',q,q)];
end
str = str(1:end-1);
str = [str ');mine = min(block(:));'];
str = [str 'e('];

for q=1:N
    str = [str sprintf('i%d,',q)];
end
str = str(1:end-1);
str = [str ') = cost + mine;'];

for q=1:N
    str = [str sprintf('\nend\n')];
end
str = [str 'e = e('];

for q=1:N
    str = [str sprintf('2:size(x{%d},2)+2,',q)];
end
str = str(1:end-1);
str = [str ');'];
eval(str);

% Find optimal path
p = ones(1,N);
p = minP(p,e,x);
if show
    if N==3
        plot3(p(:,1),p(:,2),p(:,3));
        xlabel(sprintf('S_1'));
        ylabel(sprintf('S_2'));
        zlabel(sprintf('S_3'));
        xlim([min(p(:,1)) max(p(:,1))]);
        ylim([min(p(:,2)) max(p(:,2))]);
        zlim([min(p(:,3)) max(p(:,3))]);
        grid on
    else
        if N==2;
            imshow(e(1:size(x{1},2),1:size(x{2},2)),[]);
            figure
        end
        combos = combntns(1:N,2);
        mn = factor(size(combos,1));
        mn1 = floor(length(mn)/2);
        m = prod(mn(1:mn1));
        n = prod(mn(mn1+1:end));
        for k=1:size(combos,1)
            subplot(m,n,k);
            plot(p(:,combos(k,1)),p(:,combos(k,2)));
            xlabel(sprintf('S_%d',combos(k,1)));
            ylabel(sprintf('S_%d',combos(k,2)));
        end
    end
end

% Normalized signals
for i=1:N
    y{i} = x{i}(:,p(:,i));
end

function p = minP(p,e,x)
% Find minimum e
% Last step
laststep = p(size(p,1),:);
N = length(laststep) ;

% Find minimum e
str = 'block = e(';
for q=1:N
    str = [str,sprintf('%d:%d,',laststep(q), ...
    laststep(q)+1)];
end
str = str(1:end-1);
str = [str,');'];
eval(str)
block(1) = realmax;
[mine idxe] = min(block(:));

str = '[';
for q=1:N
    str = [str,sprintf('move5(%d) ',q)];
end
str = [str , ']=ind2sub(size(block),idxe);'];
eval(str);
step = laststep+(move5-1);
p = [p;step];

str = 'fin = sum(step == [';
for q=1:N
    str = [str,sprintf('size(x{%d},2) ',q)];
end
str = [str , ']);'];
eval(str);
if fin == N
    return
end

% recursive
p = minP(p,e,x);


function d = dist(x)
% Find distance
d = 0;
for i=1:size(x,2)
    for j=i+1:size(x,2)
        d = d + sum(abs(x(:,i)-x(:,j))); % L1
    end
end


