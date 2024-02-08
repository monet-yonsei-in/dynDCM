function A=mnet_dcm2hz(A)
% A=mnet_dcm2hz(A)
% Convert DCM diagnoal terms to be real Hz
% INPUT:
% OUTPUT:
% Hae-Jeong Park @ Yonsei University, Aug. 15, 2017

if size(A,1)~=size(A,2)
    A=-0.5*exp(A);
    return;
end

for s=1:size(A,3)
    for j=1:size(A,1)
        A(j,j,s)=-0.5*exp(A(j,j,s));
    end
end