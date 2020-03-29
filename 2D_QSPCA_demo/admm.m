function [ B] = admm( X,Gram, A, B,lamda, ro, stops, maxSteps_admm, stopNum)
%ADMM Summary of this function goes here
%   Detailed explanation goes here
     %%%%%% defalut paras
     
     A=adjoint(A); 
     A=A(:,1);
     B=adjoint(B);
     B=B(:,1);
     X=adjoint(X);
     Gram=adjoint(Gram);
     if isempty(Gram)
         Gram=X'*X; 
      end
     [n, p] = size(X);
     a=complex(zeros(p,1));
     U=complex(zeros(p/2,4));

     
     step = 0; % current algorithm iteration number
     converged = false;
     
   while ~converged && step < maxSteps_admm
    step = step + 1;    

     %%%%%%%%%%%% step 1 B
    B=(Gram+(lamda+ro)*eye(p))\(Gram*A+ro*realM2c(U)-a);
     
     %%%%%%%%%%%% step 2  U
     U_old = U;     
     [Uh, Uw]=size(U);
     temp0=B+(1/ro)*a;
     temp=c2realM(temp0);
     U=zeros(size(temp));      

     
   %  fix nonzero number of U ,step2
      Unorm=sqrt(sum(temp.*temp,2));
      if stopNum < p 
        sortedUnorm = sort(Unorm, 'descend');
        UnormT=sortedUnorm(floor(stopNum) + 1);
      end
     gamma=UnormT;
     for i=1: Uh
         thresh=norm(temp(i,:),2);
         if thresh>gamma
              for j=1:4
                 t_temp=gamma*abs(temp(i,j))/thresh;
                 if temp(i,j)>0 
                   U(i,j)=temp(i,j)-t_temp;
                 else 
                   U(i,j)=temp(i,j)+t_temp;
                end
             end
         else
             U(i,:)=0;
         end 
     end
     
     %%%%%%%%%%%% step 3  a
     a=a+ro*(B-realM2c(U));


criterion1 = norm(c2q(B)-c2q(realM2c(U)),2);
 criterion2 = ro*norm(abs(U_old-U),'fro');
 if criterion1>10*criterion2
     ro=2*ro;
 elseif criterion2>10*criterion1
     ro=0.5*ro;
 else
     ro=ro;
 end

  converged = ((criterion1 < stops)&&(criterion2 < stops));
   
   end
   
   B=realM2c(U);
   B=c2q(B);
 


end




