function [outtensor,outlegdimensions] = tensordotWM( tensor1, ind1, legdimensions1, tensor2, ind2,legdimensions2 )
%TENSORDOT Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ERROR CHECK!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(ind1)~=length(ind2)
    error('The number of contracted indices must be equal for the two tensors!')
end
for i = 1:length(ind1)
    if (ind1(i) > length(legdimensions1)) || (ind2(i) > length(legdimensions2))
        error('The tensors have less legs, problem with the indices ind1 or ind2!')
    end
    if legdimensions1(ind1(i)) ~= legdimensions2(ind2(i))
        error('The tensor leg dimensions must be compatible!')
    end
    if sum(ind1(i) == ind1) ~= 1
        error('Elements of ind1 must be different. One cannot trace a leg with itself!')
    end
    if sum(ind2(i) == ind2) ~= 1
        error('Elements of ind1 must be different. One cannot trace a leg with itself!')
    end
end

%Calculations of Matlab leg numbers. Matlab cuts the dim=1 legs if they are
%in position > 2, therefore the code must be careful!

matlab_legnumber1 = length(legdimensions1);
i = length(legdimensions1);
doagain = true;
while doagain && ( i > 0)
    if legdimensions1(i) == 1
        matlab_legnumber1 = matlab_legnumber1 - 1;
    else
        doagain = false;
    end
    i = i - 1;
end
if matlab_legnumber1 < 2
    matlab_legnumber1 = 2;
end
if matlab_legnumber1 ~= length(size(tensor1))
    error('The legdimensions1 and tensor1 are not compatible')
end
size1 = size(tensor1);
for j = 1:matlab_legnumber1
    if (j==2) && length(legdimensions1) == 1
        if size1(j) ~= 1
            error('The legdimensions1 and tensor1 are not compatible')
        end
    elseif legdimensions1(j) ~= size1(j)
        error('The legdimensions1 and tensor1 are not compatible')
    end
end

matlab_legnumber2 = length(legdimensions2);
i = length(legdimensions2);
doagain = true;
while doagain && ( i > 0)
    if legdimensions2(i) == 1
        matlab_legnumber2 = matlab_legnumber2 - 1;
    else
        doagain = false;
    end
    i = i - 1;
end
if matlab_legnumber2 < 2
    matlab_legnumber2 = 2;
end
if matlab_legnumber2 ~= length(size(tensor2))
    error('The legdimensions2 and tensor2 are not compatible')
end
size2 = size(tensor2);
for j = 1:matlab_legnumber2
    if (j==2) && length(legdimensions2) == 1
        if size2(j) ~= 1
            error('The legdimensions2 and tensor2 are not compatible')
        end
    elseif legdimensions2(j) ~= size2(j)
        error('The legdimensions2 and tensor2 are not compatible')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determination of appropriate leg permutations, and the results leg
%dimensions
legperm1 = [];
remaining1 = 1:length(legdimensions1);
for i = 1:length(ind1)
    legperm1 = [legperm1,ind1(i)];
    remaining1 = remaining1(remaining1~=ind1(i));
end
legperm1 = [remaining1,legperm1];

legperm2 = [];
remaining2 = 1:length(legdimensions2);
for i = 1:length(ind2)
    legperm2 = [legperm2,ind2(i)];
    remaining2 = remaining2(remaining2~=ind2(i));
end
legperm2 = [legperm2,remaining2];

%convert leg permutations to MatLab leg permutations
matlab_perm1 = legperm1(legperm1<=matlab_legnumber1);
matlab_perm2 = legperm2(legperm2<=matlab_legnumber2);
if length(matlab_perm1) == 1
    matlab_perm1 = [2,1];
end
if length(matlab_perm2) == 1
    matlab_perm2 = [1,2];
end

%explicit calculation
legdimensions1_permuted = legdimensions1(legperm1);
legdimensions2_permuted = legdimensions2(legperm2);

outlegdimensions = [legdimensions1_permuted(1:(length(legdimensions1)-length(ind1))),...
    legdimensions2_permuted((length(ind2)+1):end)];

tensor1matrixform = reshape(permute(tensor1,matlab_perm1),...
    prod(legdimensions1_permuted(1:(length(legdimensions1)-length(ind1)))),...
    prod(legdimensions1_permuted((length(legdimensions1)-length(ind1)+1):end)));
tensor2matrixform = reshape(permute(tensor2,matlab_perm2),...
    prod(legdimensions2_permuted(1:length(ind2))),...
    prod(legdimensions2_permuted((length(ind2)+1):end)));
if length(outlegdimensions) == 0
    outtensor = reshape(tensor1matrixform*tensor2matrixform,[1,1]);
elseif length(outlegdimensions) == 1
    outtensor = reshape(tensor1matrixform*tensor2matrixform,[outlegdimensions,1]);
else
    outtensor = reshape(tensor1matrixform*tensor2matrixform,outlegdimensions);
end;


end

