%Tuning parameters: nbin, co_max, co_mean
clear;clc;close all;

% img_path = 'D:\DATA\kai_images\0_4327.jpg';
img_path = 'D:\DATA\kai_images\0_0071.jpg';

% read characters into a matrix CH
CH=imread(img_path);
% CH = rgb2gray(CH);

% grayscale to bitmap
bit_thresh = 127;
CH = CH > bit_thresh;

CH=double(CH);
nline=size(CH,1);npix=size(CH,2);
indx=find(CH<0.5);indx0=find(CH>0.5);
CH(indx)=1;CH(indx0)=0;
figure
imagesc(CH)
colorbar

% PBOD of each pixel
index=find(CH>0.5); 
[indx1,indx2]=ind2sub(size(CH),index);
num=length(indx1);
nbin=30;
degree=zeros(num,1);


pix2dis=ones(nbin,1);
for ii=1:nbin
    if ii~=nbin/2
        pix2dis(ii)=1/abs(cos(ii/nbin*pi));
    end
end


tic
% figure
% h=waitbar(0,'Be patient!')
for k=1:num
%     waitbar(k/num)
    line=indx1(k);
    pixel=indx2(k);
%     line=109;
%     pixel=59;
    
    count=zeros(nbin,1);
    for i=1:nbin
        theta=i*2*pi/nbin;
        newline=line;newpixel=pixel;
        newlinesc=newline;newpixelsc=newpixel;
        dist=0;
%          if theta==162/180*pi
        while CH(newline,newpixel)>=1
             count(i)=count(i)+1; 
             dist=dist+1;
            % find the next pixel
            toltheta=1/180*pi;
           
            if theta<pi/4 || theta>pi*7/4
                newline=round(-tan(theta)*dist+line);
                newpixel=newpixel+1;
            elseif theta<pi*5/4 && theta>pi*3/4
                newline=round(tan(theta)*dist+line);
                newpixel=newpixel-1;
            elseif theta>=pi/4&&theta<=pi*3/4
                newline=newline-1;
                newpixel=round(cot(theta)*dist+pixel);
            else
                newline=newline+1;
                newpixel=round(-cot(theta)*dist+pixel);
            end
           
            % boundary control
            if newline>nline || newline<1 || newpixel>npix || newpixel<1
                break;
            end
      
        end
%          end
    end
    thetas=(2*pi/nbin:2*pi/nbin:2*pi)*180/pi;
    % determine the degree of the pixel
%     plot(thetas,count)
    %We are counting distance, not just pixel.
    %count=count.*pix2dis;
    
    count(count<max(1/3*max(count),mean(count)))=0;
%     hold on
%     plot(thetas,count)
%     hold off
%     
    [peaks, loc]=findpeaks(count);
    numfalsepeak=0;
    
    if length(loc)>=2
        dist_=zeros(length(loc)-1,1);
        for ii=1:length(loc)-1
            dist_(ii)=abs(loc(ii)-loc(ii+1));
        end
        distmin=round(45/(360/nbin));
        falsepeak=find(dist_<distmin);
        numfalsepeak=length(falsepeak);
    end
    numpeaks=length(findpeaks(count))-numfalsepeak;
    
    [peaks2, loc2]=findpeaks([count;count]);
    numfalsepeak=0;
    if length(loc2)>=2
        dist_=zeros(length(loc2)-1,1);
        for ii=1:length(loc2)-1
            dist_(ii)=abs(loc2(ii)-loc2(ii+1));
        end
        distmin=round(45/(360/nbin));
        falsepeak=find(dist_<distmin);
        numfalsepeak=length(falsepeak);
    end
    numpeaks2=length(findpeaks([count;count]))-numfalsepeak;
    if numpeaks2>2*numpeaks
        degree(k)=numpeaks+1;
    else
        degree(k)=numpeaks;
    end
end
toc

edgeindx=sub2ind(size(CH),indx1(degree==1),indx2(degree==1));
bodyindx=sub2ind(size(CH),indx1(degree==2),indx2(degree==2));
singindx=sub2ind(size(CH),indx1(degree>=3),indx2(degree>=3));
CH(edgeindx)=2;
CH(singindx)=3;

figure; imagesc(CH); colorbar

%Tao Jia 11/17/16, 2nd code
%v2: separate by 45 degrees
%Tuning parameters: nbin, epsilon, co_max, co_mean

% % BBOD of each pixel
index=find(CH>0.5); 
[indx1,indx2]=ind2sub(size(CH),index);
num=length(indx1);
nbin=40;%dtheta=180/nbin
epsilon=1e-2;
bbod=zeros(nline, npix, nbin);
orientation=ones(nline, npix)*(-10);
degree=zeros(num,1);
tic
pix2dis=ones(nbin,1);
for ii=1:nbin
    if ii~=nbin/2
        pix2dis(ii)=1/abs(cos(ii/nbin*pi));
    end
end

% h=waitbar(0,'Be patient!')
for k=1:num
    line=indx1(k);
    pixel=indx2(k);
    
    count=zeros(nbin,1);
    for i=1:2*nbin
        theta=i*pi/nbin;
        newline=line;newpixel=pixel;
        dist=0;
         if theta==126/180*pi
             theta;
         end
        while CH(newline,newpixel)>=1
            if i>nbin
                count(i-nbin)=count(i-nbin)+1;
                dist=dist+1;
            else
                count(i)=count(i)+1; 
                dist=dist+1;
            end
            % find the next pixel
            toltheta=0.5/180*pi;
           
            if theta<pi/4 || theta>pi*7/4
                newline=round(-tan(theta)*dist+line);
                newpixel=newpixel+1;
            elseif theta<pi*5/4 && theta>pi*3/4
                newline=round(tan(theta)*dist+line);
                newpixel=newpixel-1;
            elseif theta>=pi/4&&theta<=pi*3/4
                newline=newline-1;
                newpixel=round(cot(theta)*dist+pixel);
            else
                newline=newline+1;
                newpixel=round(-cot(theta)*dist+pixel);
            end
           
            % boundary control
            if newline>nline || newline<1 || newpixel>npix || newpixel<1 
                break;
            end
            if theta>pi/2-epsilon && theta<pi/2+epsilon
                theta;
            end
        end
    end
    thetas=(pi/nbin:pi/nbin:pi)*180/pi;
   
    count=smooth(count,2);
    count(count<max(1/3*max(count),mean(count)))=0;
    
    [peaks2, loc2]=findpeaks([count;count]);
    
    if isempty(loc2)
        continue;
    end
        
    if length(loc2)>=2
        distmin=round(45/(360/nbin));
        for ii=1:length(loc2)-1
            if ii<length(loc2)
                dist_=abs(loc2(ii)-loc2(ii+1));
                if dist_<distmin
                    peaks2(ii)=[];
                    loc2(ii)=[];
                end
            end
        end
    end
    
    for ii=1:length(loc2)
        
        if loc2(ii)<=nbin
            bbod(line, pixel, loc2(ii))=peaks2(ii);
        else
            bbod(line, pixel, loc2(ii)-nbin)=peaks2(ii);
        end
    end
    
    [peak, peakInd]=max(peaks2);
    if loc2(peakInd)<=nbin%*9/8
        orientation(line, pixel)=loc2(peakInd)/nbin*180;
    else
        orientation(line, pixel)=(loc2(peakInd)-nbin)/nbin*180;
    end
end
toc
bbod=cat(3,bbod,bbod);
figure
imagesc(orientation)
colorbar
hold on

%Tao Jia 11/17/16, 3rd code
%Tuning parameters: thres, thres_simil, thres_nbin

tic

%Separate the strokes into connected components. Can change "18" to "26"
%for better connectivity.
[nline, npix, nbintimes2]=size(bbod);
nbin=nbintimes2/2;
bbod_knn=bbod;
bbod_knn(bbod_knn<0)=0;
sepa = bwconncomp(bbod_knn,18);
cells = sepa.PixelIdxList;
strokeCells = [];
NCells = max(size(cells));

%Too short a component is probably noise
thres = 40; %for Yingbi
%thres=80;%for LiuGongquan
%thres=120;%for ChuSuiliang
ii = 1;
while ii<=max(size(cells))
    
    if max(size(cell2mat(cells(ii))))<=thres
        cells(ii)=[];
    else

        ii=ii+1;
    end
end

%Convert from cell to individual stroke matrices
default=-10;
strokeMat = ones(max(size(cells)),nline, npix)*default;
pixInStroke=zeros(1,max(size(cells)));
for ii=1:max(size(cells))
%     if ii==24
%         ii;
%     end
    [indx1,indx2, indx3]=ind2sub(size(bbod),cell2mat(cells(ii)));
    pixInStroke(ii)=length(indx1);
        for jj=1:pixInStroke(ii)
            strokeMat(ii,indx1(jj),indx2(jj))=indx3(jj);
        end
%     strokeMat(ii,:,:)=full(sparse(indx1, indx2, indx3, nline, npix));
%     if ii==24
%         strokeMat1=zeros(1,nline,npix);
%         for jj=1:pixInStroke
%             strokeMat1(1,indx1(jj),indx2(jj))=indx3(jj);
%         end
%     end
    ii;
    %strokeMat=full(strokeMat);
end
strokeMat=full(strokeMat);

%As the bbod is doubled to make the "perioditicity", there will be
%redundant strokes.
[npixs,pixRank]=sort(pixInStroke);
fakeStroke=zeros(1,length(pixRank));
simil=zeros(length(pixInStroke));
difNbin=zeros(length(pixInStroke));
thres_simil=0.2;
thres_nbin=1;

% h=waitbar(0,'Be patient!');
for ii=1:length(pixInStroke)
%     waitbar(ii/length(pixInStroke))
    indexii=find(squeeze(strokeMat(pixRank(ii),:,:))>=0);
    for ind=1:length(indexii)
        [indx1,indx2]=ind2sub(size(squeeze(strokeMat(pixRank(ii),:,:))),indexii(ind));
            
        for jj=ii+1:length(pixInStroke)
            if strokeMat(pixRank(ii),indx1,indx2)~=default&&...
               strokeMat(pixRank(jj),indx1,indx2)~=default&&...
               (abs(strokeMat(pixRank(jj),indx1,indx2)-strokeMat(pixRank(ii),indx1,indx2))<=nbin+thres_nbin&&...
               abs(strokeMat(pixRank(jj),indx1,indx2)-strokeMat(pixRank(ii),indx1,indx2))>=nbin-thres_nbin)%||...
               %abs(strokeMat(pixRank(jj),indx1,indx2)-strokeMat(pixRank(ii),indx1,indx2))<=thres_nbin)
%             if (strokeMat(pixRank(jj),indx1,indx2)==strokeMat(pixRank(ii),indx1,indx2)+nbin||...
%                 strokeMat(pixRank(jj),indx1,indx2)==strokeMat(pixRank(ii),indx1,indx2)-nbin) &&...
%                 strokeMat(pixRank(ii),indx1,indx2)~=0
                simil(pixRank(ii),pixRank(jj))=simil(pixRank(ii),pixRank(jj))+1;
                difNbin(pixRank(ii),pixRank(jj))=strokeMat(pixRank(jj),indx1,indx2)...
                    -strokeMat(pixRank(ii),indx1,indx2);
                %simil(pixRank(jj),pixRank(ii))=simil(pixRank(jj),pixRank(ii))+1;
            end
        end
    end
end
for ii=1:length(pixInStroke)
    if ii==11 || ii==20
        ii;
    end
    [countSimii,simii]=max(simil(ii,:));
    if countSimii/pixInStroke(ii)>thres_simil
        fakeStroke(ii)=1;
        
        %Add the rest pixels of the fake stroke
        [rlii,rpii,rbii]=find(squeeze(strokeMat(ii,:,:)~=default));
        for jj=1:length(rlii)
            if strokeMat(simii,rlii,rpii)==default
                strokeMat(simii,rlii,rpii)=rbii+difNbin(ii,simii);
            end
        end
    end
end
fakeStrokeIndex=find(fakeStroke);
strokeMat(fakeStrokeIndex,:,:)=[];
g=figure;
g.OuterPosition=[1200 500 500 500];
noStroke=size(strokeMat,1);
sqrtnoStroke=ceil(sqrt(noStroke));
for ii=1:noStroke
    subplot(sqrtnoStroke,sqrtnoStroke,ii);
    imagesc(squeeze(strokeMat(ii,:,:)/nbin*180));
    title(['Stroke', num2str(ii)]);
end

toc

%Tao Jia 11/17/16, 4th code
%Tuning parameters: singRange, thresSing, epsilonPara, epsilonSimi

EDGE=2;
BODY=1;
SING=3;
NONE=0;
NEARSING=4;
CHStrokeBin=CH;

%If some point is close enough to SING, it is probably not good to be a
%stroke junction
singRange=3;
thresSing=0.1;
nline=size(CH,1);npix=size(CH,2);
nopixel=length(find(CHStrokeBin>=1));
index=find(CHStrokeBin==BODY|CHStrokeBin==EDGE); 
[indx1,indx2]=ind2sub(size(CHStrokeBin),index);
num=length(indx1);
for ii=1:length(index)
    pixNb=0;
    singNb=0;
    for ll=-singRange:singRange
        for pp=-singRange:singRange
            if (indx1(ii)+ll)>0 && (indx1(ii)+ll)<=nline...
                    && (indx2(ii)+pp)>0 && (indx2(ii)+pp)<=npix
                pixNb=pixNb+1;
                if CHStrokeBin(indx1(ii)+ll,indx2(ii)+pp)==SING
                    singNb=singNb+1;
                end
            end
        end
    end
    if singNb/pixNb>thresSing
        CHStrokeBin(indx1(ii),indx2(ii))=NEARSING;
    end
end
g=figure;
g.OuterPosition=[0 200 300 300];
imagesc(CHStrokeBin)
colorbar

%If two strokes share a non-SING pixel where they have very different 
%orientations, their similarity +=1
noStroke=size(strokeMat,1);
similarity=zeros(noStroke);
npixs=zeros(noStroke,1);
epsilonPara=180/nbin*1;%in degrees
for ii=1:noStroke
    strokeMatii=squeeze(strokeMat(ii,:,:));
    [lineii,pixii,binii]=find(strokeMatii>0);
    indii=sub2ind(size(strokeMatii),lineii,pixii);
    npixs(ii)=length(lineii);
    for jj=ii+1:noStroke
        strokeMatjj=squeeze(strokeMat(jj,:,:));
        [linejj,pixjj,binjj]=find(strokeMatjj>0);
        indjj=sub2ind(size(strokeMatjj),linejj,pixjj);
        npixs(jj)=length(linejj);
        for pii=1:npixs(ii)
            sameloc=find(indjj==indii(pii));
                
            %The condition: 
            %1 can find this pixel in stroke jj &&
            %(2 this pixel is BODY or EDGE ||
            %3 the two strokes are paralell at this pixel)
            if ~isempty(sameloc)&&...
                (CHStrokeBin(lineii(pii),pixii(pii))==BODY ||...
                CHStrokeBin(lineii(pii),pixii(pii))==EDGE ||...
                mod(abs(strokeMatii(lineii(pii),pixii(pii))-strokeMatjj(lineii(pii),pixii(pii))),nbin)...
                <epsilonPara/180*nbin)
                if ii==8&&jj==9
                    ii;
                end
                similarity(ii,jj)=similarity(ii,jj)+1;
            end
        end
    end
end

%If two strokes have high similarity, they are connected
epsilon_simi=5e-2;
connected=zeros(noStroke);
for ii=1:noStroke
    for jj=ii+1:noStroke
        %if similarity(ii,jj)/npixs(ii)>epsilon_simi||similarity(ii,jj)/npixs(jj)>epsilon_simi
        if similarity(ii,jj)>=nopixel*epsilon_simi/noStroke
            connected(ii,jj)=1;
            connected(jj,ii)=1;
        end
    end
end

%If two strokes are connected, they share the same label, and are combined
%into a same stroke in newStrokeMat
noRealStroke=noStroke;
labelStroke=1:noStroke;
gStroke=graph(connected);
binLabel=conncomp(gStroke);
g=figure;
g.OuterPosition=[0 500 500 500];
noNewStroke=max(binLabel);
sqrtnoNewStroke=ceil(sqrt(noNewStroke));
newStrokeMat=zeros(noNewStroke,size(strokeMat,2),size(strokeMat,3));

for ii=1:noNewStroke
    for jj=1:noStroke
        if binLabel(jj)==ii
            %newStrokeMat(ii,:,:)=strokeMat(jj,:,:);
            newStrokeMat(ii,:,:)=newStrokeMat(ii,:,:)+strokeMat(jj,:,:);
        end
    end
    subplot(sqrtnoNewStroke,sqrtnoNewStroke,ii);
    imagesc(squeeze(newStrokeMat(ii,:,:)/nbin*180));
    title(['BinnedStroke', num2str(ii)]);
    drawnow
end

