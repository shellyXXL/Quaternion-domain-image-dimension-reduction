function c = realM2c_q( realM )
%REALM2C Summary of this function goes here
%   Detailed explanation goes here
   [h,w]=size(realM);
   
   a= complex(realM(:,1),realM(:,2));
   b= complex(realM(:,3),realM(:,4));
   %b=-conj(b);
   c=[a;b];
end

