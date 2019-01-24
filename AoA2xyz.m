function [mxTargets] = AoA2xyz(rvR, rvAlph)

% Build xyz coordinates for each Rx
mxTargets = zeros(3,size(rvR,2));
mxTargets(:,1) = getXYZ(rvR(1), rvAlph([1 4]));
mxTargets(:,2) = getXYZ(rvR(2), rvAlph([1 2]));
mxTargets(:,3) = getXYZ(rvR(3), rvAlph([3 2]));
mxTargets(:,4) = getXYZ(rvR(4), rvAlph([3 4]));
end

function [cvOut] = getXYZ(scR, rvAlph)
% alph(1) is for the x-axis
% alph(2) is for the y-axis
if scR==0
    cvOut=NaN(3,1);
else
    cvOut(1:2,:) = scR*sin(rvAlph');
%     r^2-x^2-y^2=z^2
    cvOut(3) = real(sqrt(scR^2-cvOut(1)^2-cvOut(2)^2));
end
end

