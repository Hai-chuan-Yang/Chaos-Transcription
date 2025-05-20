function [popu] = BoundaryDetection(popu,lu)
[NP, D] = size(popu);
xl = repmat(lu(1,:), NP, 1);
xu = repmat(lu(2,:), NP, 1);
ranX=rand(NP, D);
%check the lower bound
pos1 = popu < xl;
popu(pos1) = xl(pos1) + ranX(pos1) .* (xu(pos1) - xl(pos1));
%check the upper bound
pos2 = popu > xu;
popu(pos2) = xl(pos2) + ranX(pos2) .* (xu(pos2) - xl(pos2));
end