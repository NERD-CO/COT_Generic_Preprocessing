function savestruct = make_save_struct(varargin)

for argi = 1:nargin
    savestruct.(inputname(argi)) = varargin{argi};
end