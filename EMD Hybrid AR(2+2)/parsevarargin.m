function varstruct = parsevarargin(v, templatestruct)
%% Input Arguments 
%v: varargin, variable length input argument list
% templatestruct = structure with fieldnames and default values
%% Output Arguments
% varstruct: structure with its field values set to the input ones
  
  if ~exist('templatestruct') %#ok<EXIST>
    [~, vars] = parseparms(v);
    fnames = vars(1:2:end); % names of the variables (strings)
    vals = vars(2:2:end); % values of the variables
    varstruct = [];
    for k = 1:length(fnames) % set the value of the variable
      varstruct.(fnames{k})= vals{k};
    end
  else
    varstruct = templatestruct(1);
    F = fieldnames(templatestruct);
    J = 1;
    while J <= (length(v)-1)
      if ischar(v{J})        
        id  = find(strcmpi(F, v(J)));
        if ~isempty(id)
          varstruct.(F{id(1)})= v{J+1}; 
          J = J + 1;
        end
      end
      J = J + 1;
    end
  end
    
    
    
