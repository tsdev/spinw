function result = sw_issymop(Op)
% function determines whether the matrix is a symmetry operator
%
% result = SW_ISSYMOP(Op)
%
% Op        Symemtry operators with rotation and translation. Dimensions
%           are [3 4 nOp].
%

if size(Op,1) == 3 && size(Op,2) == 4
    result = true;
else
    result = false;
end

end