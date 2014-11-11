function fprintf0(fid,varargin)
% same as fprintf, expect if fid is zero, will not produce output.

if fid~=0
    fprintf(fid,varargin{:});
end

end