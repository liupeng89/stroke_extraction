img_path = 'D:\matlab_projects\stroke_extraction\kai_images';
files = dir(img_path);

file_names = {};
for i=1: length(files)
    if contains(files(i).name, '.jpg')
        file_names{end+1} = files(i).name;
    end
end

LEN = length(file_names);
disp(LEN);

for i=1: LEN
    img_url = strcat(img_path, '\', file_names(i));
    
    strokes = Stroke_extraction(img_url);
end