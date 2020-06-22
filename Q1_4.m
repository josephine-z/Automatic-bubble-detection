clear all; close all; clc; clear memory;


%% INVARIANT MOMENT EQUAITONS
xcg=    @(I,xg,yg)      sum(xg(:).*I(:))./sum(I(:));
ycg=    @(I,xg,yg)      sum(yg(:).*I(:))./sum(I(:));
mu=     @(I,p,q,xg,yg)  sum(((xg(:)-xcg(I,xg,yg)).^p).*((yg(:)-ycg(I,xg,yg)).^q).*I(:));
n=      @(I,p,q,xg,yg)	mu(I,p,q,xg,yg)./(mu(I,0,0,xg,yg).^(1+0.5*(p+q)));
M1=     @(BW,xg,yg)     n(BW,2,0,xg,yg)+n(BW,0,2,xg,yg);
M2=     @(BW,xg,yg)     (n(BW,2,0,xg,yg)-n(BW,0,2,xg,yg)).^2 + 4*n(BW,1,1,xg,yg).^2;
M3=     @(BW,xg,yg)     (n(BW,3,0,xg,yg)-3*n(BW,1,2,xg,yg)).^2 + (3*n(BW,2,1,xg,yg)-n(BW,0,3,xg,yg)).^2;
M4=     @(BW,xg,yg)     (n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).^2 + (n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).^2;
M5=     @(BW,xg,yg)     (n(BW,3,0,xg,yg)-3*n(BW,1,2,xg,yg)).*(n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)) + (n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).^2 -3*(n(BW,2,1,xg,yg)-n(BW,0,3,xg,yg)).^2 + (3*n(BW,2,1,xg,yg)-n(BW,0,3,xg,yg)).*(n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).*(3*(n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).^2 - (n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).^2);
M6=     @(BW,xg,yg)     (n(BW,2,0,xg,yg)-n(BW,0,2,xg,yg)).*((n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).^2-(n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).^2) + 4*n(BW,1,1,xg,yg).*(n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).*(n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg));
M7=     @(BW,xg,yg)     (3*n(BW,2,1,xg,yg)-n(BW,0,3,xg,yg)).*(n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).*((n(BW,3,0,xg,yg)+n(BW,1,2,xg,yg)).^2-3*(n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).^2) +(3*n(BW,1,2,xg,yg)-n(BW,3,0,xg,yg)).*(n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).*(3*(n(BW,1,2,xg,yg)+n(BW,3,0,xg,yg)).^2-(n(BW,2,1,xg,yg)+n(BW,0,3,xg,yg)).^2);


%% MAKE A CIRCLE & PERFORM IMAGE MOMENTS
numpix=255;                                     %image size
numpix=1+2*floor(numpix/2);                     %make sure numpix is odd
Im=zeros(numpix,numpix);                        %initialise image
xc=0.5*(1+size(Im,2));                          %x centre of circle
yc=0.5*(1+size(Im,1));                          %y centre of circle
[xgc,ygc]=meshgrid([1:numpix]);                 %create grid with x and y going from 1 to numpix
R=100;											%radius of circle
F=(xgc-xc).^2+(ygc-yc).^2-R^2;	                %calculate equation of circle less radius of circle. This finds all points within the circle
Im(F<0)=1;                                      %make pionts inside circle equal 1
cxg=repmat(1:numpix,1,length(Im));              %x points in Im
cyg=repmat(1:numpix,1,length(Im));              %y points in Im
cM1=M1(Im,cxg,cyg);                             %begin: invariant moment calculations
cM2=M2(Im,cxg,cyg);
cM3=M3(Im,cxg,cyg);
cM4=M4(Im,cxg,cyg);
cM5=M5(Im,cxg,cyg);
cM6=M6(Im,cxg,cyg);
cM7=M7(Im,cxg,cyg);                             %end: invariant moment calculations


%% READ IMAGE & DETECT BOUNDARY
Imrgb=imread('images/Image_Q1.jpg');            %read rgb image
Imgray=rgb2gray(Imrgb);                         %convert to grayscale
level=graythresh(Imgray);                       %calculate the threshold
BW=im2bw(Imgray,level);                         %convert the image to a black and white image with intensities above the threshold equal to 1
BW=edge(BW,'canny',[],4);                       %edge detection on binary threshold image
SE=strel('disk',16);                            %define structural element
BWc=imcomplement(imclose(BW,SE));               %close contours within structural element & inverse colours to make blobs in image have a value of 1 over a background of 0


%% FIND BLOBS/BUBBLE IN IMAGE & THEIR MOMENTS
byg=repmat(transpose(1:size(BWc,1)),1,size(BWc,2));                           %y points in BWc
bxg=repmat(1:size(BWc,2),size(BWc,1),1);                                      %x points in BWc
[L,num]=bwlabel(BWc,4);                                                       %label all blobs in image
for i=1:num                                                                   %from 1 to the number of blobs in BWc:
    [row,col]=find(L==i);                                                     %find rows and columns of L = counter
    k=convhull(col,row);                                                      %find line surrounding points in blob
    mask(:,:,i)=poly2mask(col(k),row(k),size(Imrgb,1),size(Imrgb,2));         %new image of filled in blob
    bM1(i)=M1(mask(:,:,i),bxg,byg);                                           %begin: invariant moment calculations
    bM2(i)=M2(mask(:,:,i),bxg,byg);
    bM3(i)=M3(mask(:,:,i),bxg,byg);
    bM4(i)=M4(mask(:,:,i),bxg,byg);
    bM5(i)=M5(mask(:,:,i),bxg,byg);
    bM6(i)=M6(mask(:,:,i),bxg,byg);
    bM7(i)=M7(mask(:,:,i),bxg,byg);                                           %end: invariant moment calculations
end
[c1,ind1]=min(abs(cM1-bM1));                                                  %begin: find value and index of blob with closest moment to circle
[c2,ind2]=min(abs(cM2-bM2));
[c3,ind3]=min(abs(cM3-bM3));
[c4,ind4]=min(abs(cM4-bM4));
[c5,ind5]=min(abs(cM5-bM5));
[c6,ind6]=min(abs(cM6-bM6));
[c7,ind7]=min(abs(cM7-bM7));                                                  %end: find value and index of blob with closest moment to circle
fi=mode(horzcat(ind1,ind2,ind3,ind4,ind5,ind6,ind7));                         %find index of blob that is most like a circle (completed using mode average)


%% DRAW CIRCLES AROUND BUBBLE
[centres, radii]=imfindcircles(mask(:,:,fi),[100 length(mask(:,:,fi))],'Sensitivity',0.9);     %find radius and centre of outer circle of bubble
d=find(L==fi);                                                                                 %find index of bubble
bIm=zeros(size(BWc));                                                                          %initialise new image
bIm(d)=1;                                                                                      %image of bubble
zoom=radii-25;                                                                                 %amount to zoom in by to find inner circle of bubble
bIm1=imcomplement(bIm(centres(2)-zoom:centres(2)+zoom,centres(1)-zoom:centres(1)+zoom));       %zoomed-in image of inner circle of bubble, and swap 1's and 0's in image
byg1=repmat(transpose(1:size(bIm1,1)),1,size(bIm1,2));                                         %y points in inner circle image            
bxg1=repmat(1:size(bIm1,2),size(bIm1,1),1);                                                    %x points in inner circle image               
[L1,num1]=bwlabel(bIm1,4);                                                                     %label all blobs on image
for i=1:num1                                                                                   
    [row1,col1]=find(L1==i);                                                                   %gives rows and columns of L==i
    k=convhull(col1,row1);                                                                     %find line surrounding points found                                                                     
    mask1(:,:,i)=poly2mask(col1(k),row1(k),size(bIm1,1),size(bIm1,2));                         %image of filled in blob                          
    ibM1(i)=M1(mask1(:,:,i),bxg1,byg1);                                                        %begin: invariant moment calculations     
    ibM2(i)=M2(mask1(:,:,i),bxg1,byg1);
    ibM3(i)=M3(mask1(:,:,i),bxg1,byg1);
    ibM4(i)=M4(mask1(:,:,i),bxg1,byg1);
    ibM5(i)=M5(mask1(:,:,i),bxg1,byg1);
    ibM6(i)=M6(mask1(:,:,i),bxg1,byg1);
    ibM7(i)=M7(mask1(:,:,i),bxg1,byg1);                                                        %end: invariant moment calculations
end
[c1,ind1]=min(abs(cM1-ibM1));                                                                  %begin: find value and index of blob with closest moment to circle
[c2,ind2]=min(abs(cM2-ibM2));
[c3,ind3]=min(abs(cM3-ibM3));
[c4,ind4]=min(abs(cM4-ibM4));
[c5,ind5]=min(abs(cM5-ibM5));
[c6,ind6]=min(abs(cM6-ibM6));
[c7,ind7]=min(abs(cM7-ibM7));                                                                  %end: find value and index of blob with closest moment to circle
fi=mode(horzcat(ind1,ind2,ind3,ind4,ind5,ind6,ind7));                                          %find index of blob that is most like a circle (completed using mode())
[centres1, radii1]=imfindcircles(mask1(:,:,fi),[10 length(mask1(:,:,fi))],'Sensitivity',0.86); %find radius and centre of inner circle
figure(1),hold on;imagesc(flipud(Imrgb));h=viscircles(centres,radii);h1=viscircles(centres1+centres-zoom,radii1); %plot original image with circles around inner and outer radius





