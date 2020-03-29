function q = c2q( c )
%Q2C Summary of this function goes here
%   Detailed explanation goes here
[h,w] = size(c);
% if w==1
 H=h/2;
 a=c(1:H,:);
 b=-conj(c(H+1:h,:));
 q=dc(a,b);
% end
% if w~=1
%    disp([' error in q2c']);
% end

end

