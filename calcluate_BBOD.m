% Tuning parameters: nbin, epsilon, co_max, co_mean

% Calculate BBOD of each pixel
index = find(IMG > 0.5);
[index_rows, index_cols] = ind2sub(size(IMG), index);

num = length(index_rows);
nbin = 30; % theta = 180 / nbin
epsilon = 1e-2;
bbod = zeros(nRows, nColumns, nbin);

orientation = ones(nRows, nColumns) * (-10);
degree = zeros(num, 1);

tic 
pix2dis = ones(nbin, 1);
for ii = 1:nbin
    if ii ~= nbin / 2
        pix2dis(ii) = 1 / abs(cos(ii / nbin * pi));
    end
end

for k = 1:num
    
    id_row = index_rows(k);
    id_col = index_cols(k);
    
    count = zeros(nbin, 1);
    
    for i = 1:2*nbin
        theta = i * pi / nbin;
        
        new_row = id_row; new_col = id_col;
        dist = 0;

        while IMG(new_row, new_col) >= 1
            if i > nbin
                count(i - nbin) = count(i - nbin) + 1;
                dist = dist + 1;
            else
                count(i) = count(i) + 1;
                dist = dist + 1;
            end
            
            toltheta = 0.5 / 180 * pi;
            
            if theta < pi/4 || theta >= pi*7/4
                new_row = round(-tan(theta) * dist + id_row);
                new_col = new_col + 1;
            elseif theta < pi*5/4 && theta > pi*3/4
                new_row = round(tan(theta) * dist + id_row);
                new_col = new_col - 1;
            elseif theta >= pi/4 && theta <= pi*3/4
                new_row = new_row - 1;
                new_col = round(cot(theta) * dist + id_col);
            else
                new_row = new_row + 1;
                new_col = round(-cot(theta) * dist + id_col);
            end
            
            % boundary control
            if new_row > nRows || new_row < 1 || new_col > nColumns || new_col < 1
                break;
            end    
        end
    end
    
    thetas = (pi/nbin: pi/nbin: pi) * 180 / pi;
    
    count = smooth(count, 2);
    count(count < max(1/3*max(count), mean(count))) = 0;
    
    [peaks2, loc2] = findpeaks([count; count]);
    if isempty(loc2)
        continue;
    end
    
    if length(loc2) >= 2
        dist_min = round(45 / (360 / nbin));
        for ii = 1: length(loc2)-1
            if ii < length(loc2)
                dist_ = abs(loc2(ii) - loc2(ii+1));
                if dist_ < dist_min
                    peaks2(ii) = [];
                    loc2(ii) = [];
                end
            end
        end
    end
    
    for ii = 1: length(loc2)
        if loc2(ii) <= nbin
            bbod(id_row, id_col, loc2(ii)) = peaks2(ii);
        else
            bbod(id_row, id_col, loc2(ii)-nbin) = peaks2(ii);
        end
    end
    
    [peak peakInd] = max(peaks2);
    if loc2(peakInd) <= nbin
        orientation(id_row, id_col) = loc2(peakInd) / nbin * 180;
    else
        orientation(id_row, id_col) = (loc2(peakInd) - nbin) / nbin * 180;
    end
end

toc
bbod = cat(3, bbod, bbod);
figure; imagesc(orientation); colorbar; hold on;