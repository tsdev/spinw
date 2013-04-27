function cGraph = sw_fsub(conn, ~)
% cGraph = SW_FSUB(conn, nExt) simple graph vertex coloring, determines non
% connected sublattices for Monte-Carlo calculation.
%
% Input:
%
% conn          Contains edge indices which are connected
%               conn(1,idx)-->conn(2,idx), dimensions are [2 nConn].
% nExt          Size of the magnetic unit cell in units of cells.
%
% Output:
%
% cGraph        Vector, that assigns every magnetic moment to a sublattice.
%
% See also SW.ANNEAL.
%

nEdge  = max(max(conn));
cGraph = zeros(1,nEdge);

for ii = 1:nEdge
    idx = [conn(2,conn(1,:)==ii) conn(1,conn(2,:)==ii)];
    
    cNeigh = cGraph(idx);
    cIdx = 1;
    while any(cNeigh==cIdx)
        cIdx = cIdx + 1;
    end
    cGraph(ii) = cIdx;
end

end