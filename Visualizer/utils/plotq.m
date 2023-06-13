function  plotq(varargin)
%plots the squeeze of the array
if length(varargin)>=2 && length(varargin{1})>3 && length(varargin{2})>3
    plot(squeeze(varargin{1}),squeeze(varargin{2}),varargin{3:end})
else
    plot(squeeze(varargin{1}),varargin{2:end})
end
end

