%Shrinkage and Thresholding Operator Routine 
%It works on a real number:
%  S_lambda(qi) =   | qi-lambda     qi>lambda
%                   | qi+lambda     qi<-lambda  
%                   | 0             qi \in [-lambda, lambda]
%-------------------------------------------------
function STO = sto(qi,lambda)
STO=0;
    if qi>lambda
        STO=qi-lambda;   
    elseif qi < -lambda
        STO=qi+lambda; 
    elseif (abs(qi)<=lambda)
        STO=0;
    end
end

