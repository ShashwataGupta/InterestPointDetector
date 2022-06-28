function [ keyPts keyPtsLoc DoG2 ] = getSIFTkeyPoints( imgTemp )

[M N] = size(imgTemp);
imgTemp = double(imgTemp);
imgOrg = imgTemp;

pyramid1 = [];
pyramid2 = [];
pyramid3 = [];

% ---------------------------  Octave 1 ----------------------------------
pyramidLevel = 0;

imgTemp(M : (M + 40), N : (N + 40)) = 0;
imgIn = int8(imgTemp);
% imgTemp = [imgTemp zeros(M,6)];
% imgTemp = [imgTemp ; zeros(6,N)];
clear conSum;
k=(sqrt(2));
for scaleCount = 0 : 3
    
sigmaMask = (k ^ scaleCount) * 1.6;
mask_size = round(3*sigmaMask);
% sigmaMask = (rootTwo ^ (scaleCount + (2 * pyramidLevel))) * 1.6;

%-------------------------- Generating Gaussian Mask ----------------------
for x = -mask_size : mask_size
    
    for y=-mask_size:mask_size
              
        h(x + mask_size+1,y + mask_size+1) = (1 / ((2 * pi) * ((sigmaMask) * (sigmaMask))) ) * exp(-((x * x)+(y * y))/(2 * (sigmaMask * sigmaMask))); %Gaussian Function
        
    end 
    
end

% h = fspecial('gaussian',7,)
%--------------------------------------------------------------------------

%Applying Mask
for i = 1 : M    
    for j = 1 : N
        
        t = imgTemp(i : i + 2*mask_size, j : j + 2*mask_size)' .* h;
        conSum(i, j) = sum(sum(t));
        
    end    
end
% conSum = conv2(imgTemp, h, 'valid');
pyramid1 = [pyramid1 conSum];

end



[M, N]=size(imgOrg);

% ------------------------ Extracting Key Points --------------------------
DoG1 = pyramid1(1 : M, 1 : N) - pyramid1(1 : M, N + 1 : 2 * N);
DoG2 = pyramid1(1 : M, N + 1 : 2 * N) - pyramid1(1 : M, 2 * N + 1 : 3 * N);
DoG3 = pyramid1(1 : M, 2 * N + 1 : 3 * N) - pyramid1(1 : M, 3 * N + 1 : 4 * N);

[M, N] = size(DoG2);
keyPts = [];
keyPtsLoc = [];

for i = 2 : M-1
    
    for j = 2 : N-1
        
        x = DoG1(i - 1 : i + 1, j - 1 : j + 1);
        y = DoG2(i - 1 : i + 1, j - 1 : j + 1);
        y_bkp = y;
        z = DoG3(i - 1 : i + 1, j - 1 : j + 1);
        
        y(1 : 4) = y(1 : 4);
        y(5 : 8) = y(6 : 9);
        
        xMax = max(max(x));
        xMin = min(min(x));
        
        yMax = max(max(y));
        yMin = min(min(y));
        
        zMax = max(max(z));             
        zMin = min(min(z));
        
        gMax = max([xMax yMax zMax]);
        gMin = min([xMin yMin zMin]);
               
        if ((DoG2(i,j) > gMax || DoG2(i,j) < gMin) ) % If global extrema detected, it is potential keypoint   
            
            
            roi_x = (y_bkp(2,3) - y_bkp(2,1)) / 2;
            roi_y = (y_bkp(3,2) - y_bkp(1,2)) / 2;
            roi_s = (z(2,2) - x(2,2)) / 2;
            
            roi_xx = y_bkp(2,3) - 2*y_bkp(2,2) + y_bkp(2,1);
            roi_yy = y_bkp(3,2) - 2*y_bkp(2,2) + y_bkp(1,2);
            roi_ss = z(2,2) - 2*y_bkp(2,2) + x(2,2);
            
            roi_xy = ( (y_bkp(1,1) + y_bkp(3,3)) - (y_bkp(1,3) + y_bkp(3,1)) ) / 4 ;
            roi_xs = ( (x(2,1) + z(2,3)) - (x(2,3) + z(2,1)) ) / 2;
            roi_ys = ( (x(1,2) + z(3,2)) - (x(3,2) + z(1,2)) ) / 2;
            
            TrH = roi_xx + roi_yy;
            DetH = (roi_xx * roi_yy) - (roi_xy * roi_xy);
            
            edge_rat = (TrH ^ 2) / DetH;
            
            if(abs(DoG2(i,j)) > 7.65 && edge_rat < 12.1) % Rejecting edges and low contrast images
                keyPts = [keyPts DoG2(i, j)];
                keyPtsLoc = [keyPtsLoc [i ; j] ];
            end
            
        end
        
    end
    
end

end

