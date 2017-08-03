function cGraph = sw_fsub(conn, ~)
% simple graph vertex coloring
%
% cGraph = SW_FSUB(conn, nExt)
%
% It creates a simple graph vertex coloring, determines non connected
% sublattices for Monte-Carlo calculation.
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
% See also SPINW.ANNEAL.
%

if nargin == 0
    help sw_fsub;
    return
end

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