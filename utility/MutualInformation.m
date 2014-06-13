function [ val ] = MutualInformation( I1, I2 )
%MFunc compute mutual information
val = H1Func( I1 ) + H1Func( I2 ) - H2Func( I1, I2 );


end

function [ val ] = H1Func( I1 )
%H1 function, compute the entropy of one image
h = size(I1, 1);
w = size(I1, 2);
n = histc(reshape(I1, 1, h*w),[0:255]);
val = 0;
for i = 1:256
    if n(i) ~= 0
        p = n(i)/(h*w);
        val = val - p*log2(p);
    end
end

end

function [ val ] = H2Func( I1, I2 )
%H1 function, compute the entropy of one image
h = size(I1, 1);
w = size(I1, 2);
nmap = zeros( 256, 256 );
for i = 1:h
    for j = 1:w
        p1 = I1(i,j);
        p2 = I2(i,j);
        nmap(p1+1, p2+1) = nmap(p1+1, p2+1) + 1;
    end
end
val = 0;
for i = 1:256
    for j = 1:256
        if nmap(i, j) ~= 0
            p = nmap(i,j) / (h*w);
            val = val - p*log2(p);
        end
    end
end


end
