function files_names = wildfiles(varargin);

%WILDFILES converts a wildcard expression to a list of files in a cell array.  
%
% 
% NOTE:     * Tested only in (Windows) but not in (UNIX).
%           * If no files found, the function will return (0).
%
% Examples:
% 
%   wildfiles *.*
%   D = wildfiles('m*.*');
%   D = wildfiles('m*.*','s*.m');
%
%   Copyright 2004 Fahad Al Mahmood
%   Version: 1.0 $  $Date: 1-May-2004


k=1;
for n=1:nargin
    [pathstr,name,ext]=fileparts(varargin{n});
    if isempty(name) & ~isempty(ext)
        name = '*';
    end
    if strmatch(ext,'.','exact')
        ext = '.*';
    end
    full_name = [pathstr filesep name ext];
    if isempty(pathstr)
        full_name = [name ext];
    end
    if isdir([pathstr filesep name])
        full_name = [pathstr filesep name filesep '*.*'];
    end
    files = dir(full_name);
    
    % Assigning file names to a cell array after making sure they do not
    % contain directories.
    for m=1:length(files)
        if ~files(m).isdir
            files_names{k} = files(m).name;
            k=k+1;
        end
    end
    
    % Flag (0) is case no file is found.
    if ~exist('files_names','var')
        files_names = [];
    end
    files_names = files_names';
end