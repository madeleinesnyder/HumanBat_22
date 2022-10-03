function figCount = saveFig(directory,name,id,go)

%=== Find all open figures
figHandles = findall(0,'Type','figure');
figHandles = flip(figHandles);

for i = 1:numel(figHandles)
    if go
        saveas(figHandles(i),[directory, '/', name '_figure' num2str(id+i-1) '.png']);   close(figHandles(i));
    else
        figCount = id;
    end
end

figCount = id+numel(figHandles);

end

