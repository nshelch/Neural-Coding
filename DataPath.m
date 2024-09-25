function path = DataPath()

pcID = getenv('COMPUTERNAME');

if strcmpi(pcID, 'OBA-PC-01')
    path = 'S:\UserFolders\CharlesGreenspon\BCI_TactileCoding';
elseif strcmpi(pcID, 'DESKTOP-LEG2SE6')
    path = 'U:\UserFolders\CharlesGreenspon\BCI_TactileCoding';
else
    path = 'U:\UserFolders\CharlesGreenspon\BCI_TactileCoding';
end

end