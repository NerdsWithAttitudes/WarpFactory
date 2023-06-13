function  surfq(varargin)
%surfs the squeeze of the array
if length(varargin)>=3 && length(varargin{1})>3 && length(varargin{2})>3 && length(varargin{3})>3 && isa(varargin{2},'double') && isa(varargin{3},'double')
    surf(varargin{1},varargin{2},squeeze(varargin{3}),varargin{4:end})
    colormap(redblue(squeeze(varargin{3})))
else
    surf(squeeze(varargin{1})',varargin{2:end})
    colormap(redblue(squeeze(varargin{1})))
end
end

