function x = iindex(x,varargin)
% inline indexing using parenthesis or braces

switch varargin{1}
    case '()'
        x = x(varargin{2:end});
    case '{}'
        x = x{varargin{2:end}};
    otherwise
end

end