function colmap = cmasher(name,varargin)
if isunix 
setenv('PATH', getenv('PATH')+":/usr/local/bin")
end 
[~,syslist] = system('cmr cmlist');
clist = textscan(syslist,'%s');
ntotcolors = length(clist{:});
colorlist = [];
for i = 1:ntotcolors
    colorlist= [colorlist; string(clist{1,1}{i})];
end

valid_color = ~isempty(find(name==colorlist,1));

if valid_color

    if nargin==2
        ncols = varargin{1};
        command = strjoin(strcat(['cmr cmcolors cmr.',string(name),' ',num2str(int8(ncols))]),'');
    else
        command = strjoin(strcat(['cmr rgbtable cmr.',string(name)]),'');
    end

    [~,sysresult] = system(command);
    colmap = cell2mat(textscan(sysresult,'%f%f%f'));

else
    colmap= [];
    disp('not a valid color name, valid colors are:')
    disp(colorlist)
end