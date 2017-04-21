function [ F ] = assembleRHS( U, h, R )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F = abs(U(3:end-2)).^2.*U(3:end-2)...
  - (-U(1:end-4) + 16*U(2:end-3) -30*U(3:end-2) + 16*U(4:end-1) - U(5:end))/(12*h^2)...
  - (U(1:end-4) - 8*U(2:end-3) + 8*U(4:end-1) - U(5:end))*4./R/12*h;

end

