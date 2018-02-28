function [ deriv1, deriv2 ] = find_derivative(h, array, deriv_method )
%This function uses central or forward difference to approximate the
%first and second derivatives of an array
% Approximate difference to find first and second derivative of LVP
if strcmp(deriv_method,'forward')
    deriv1=diff(array)/h;
    deriv2=diff(deriv1)/h;
elseif strcmp(deriv_method,'central')
    deriv1=zeros(1,length(array));
    deriv2=zeros(1,length(array));
    deriv1(1)=(array(2)-array(1))/h;                    %forward difference for first element
    deriv1(end)=(array(end)-array(end-1))/h;            %backward difference for last element
   %central difference for 1st and second order
    for i=2:(length(array)-1)
       deriv1(i)=(array(i+1)-array(i-1))/(2*h);
       deriv2(i)=(array(i+1)-2*array(i)+array(i-1))/h^2;
    end
    deriv2(1)=(deriv1(2)-deriv1(1))/h;                  %forward diff
    deriv2(end)=(deriv1(end)-deriv1(end-1))/h;          %backwards diff
end


end

