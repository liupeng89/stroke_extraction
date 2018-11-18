% Get the connected components of character from BW images based on the 
% labeling connected components algorithm.
function components = Components_split(grayscale, conn)
    
    % change black character to white color;
    grayscale = ~grayscale;
    [L N] = bwlabel(grayscale, conn);
    
    components = {};
    for i =1: N
        bg = L ~= i;
        components{end+1} = bg;
    end
end