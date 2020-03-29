function realMat = c2realM_q( c )

  [h,w]=size(c);
  realMat=zeros(h/2,4);
  realMat(:,1)=real(c(1:h/2,1));
  realMat(:,2)=imag(c(1:h/2,1));
  realMat(:,3)=real(c(h/2+1:end,1));
  realMat(:,4)=imag(c(h/2+1:end,1));

end