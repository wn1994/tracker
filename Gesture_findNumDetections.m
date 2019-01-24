function [gestureMetrics] = Gesture_findNumDetections(numRangeBins,numDopplerBins,preDetMatrix,detThresh)
 wt=0;maxVal=0;
 maxRangeIdx=0;
 maxDopplerIdx=0;
 numDetections=0;
 numDetectionsNeg=0;
 numDetectionsPos=0;
 dopplerAvg=0;
 rangeAvg=0;
 wtSum=0;
 dopplerAvgNeg=0;
 rangeAvgNeg=0;
 wtSumNeg=0;
 dopplerAvgPos=0;
 rangeAvgPos=0;
 wtSumPos=0;
 valatindx = 0;
%  tempPtr=preDetMatrixPtr+numDopplerBins*NUM_INITIAL_RANGE_BINS_SKIP+NUM_INITIAL_DOPPLER_BINS_SKIP;
%  Computing statistics for positive dopplers
% !!!
NUM_INITIAL_RANGE_BINS_SKIP=1;
NUM_INITIAL_DOPPLER_BINS_SKIP=1;
MAX_RANGE_BIN_USED=numDopplerBins;
MAX_DOPPLER_BIN_USED_ONE_SIDED=numDopplerBins / 2;

    for i=NUM_INITIAL_RANGE_BINS_SKIP:MAX_RANGE_BIN_USED
        % 前半部分是正
        for j=NUM_INITIAL_DOPPLER_BINS_SKIP+1:MAX_DOPPLER_BIN_USED_ONE_SIDED
            valatindx = preDetMatrix(i,j);
                % 如果大于阈值，进行计算
            if valatindx>detThresh
                % 累加符合条件的点
                numDetectionsPos = numDetectionsPos+1;
                wt=valatindx;
                % 对该点的Doppler和range值累加
                dopplerAvgPos=dopplerAvgPos+wt*j;
                rangeAvgPos=rangeAvgPos+(wt+1000)*i;
                % 累加符合计算条件的点的值
                wtSumPos=wtSumPos+wt;
            end
                %记录最大点的位置i,j
            if valatindx>maxVal
                maxVal=valatindx;
                maxRangeIdx=i;
                maxDopplerIdx=j;
            end
        end
    end
  
    for i=NUM_INITIAL_RANGE_BINS_SKIP:MAX_RANGE_BIN_USED
		for j=numDopplerBins-MAX_DOPPLER_BIN_USED_ONE_SIDED:numDopplerBins-NUM_INITIAL_DOPPLER_BINS_SKIP-1
            valatindx = preDetMatrix(i,j);
            if valatindx>detThresh
                numDetectionsNeg=numDetectionsNeg+1;
                wt=valatindx;
                dopplerAvgNeg=dopplerAvgNeg+wt*(j-numDopplerBins);
                rangeAvgNeg=rangeAvgNeg+(wt+1000)*i;
                wtSumNeg=wtSumNeg+wt;
            end
            if valatindx>maxVal
                    maxVal=valatindx;
                    maxRangeIdx=i;
                    maxDopplerIdx=j;
            end
        end
    end
    
    % 检测到的点是正负相加的数目
    numDetections=numDetectionsPos+numDetectionsNeg;
	wtSum=wtSumPos+wtSumNeg;
	dopplerAvg=dopplerAvgPos+dopplerAvgNeg;
	rangeAvg=rangeAvgPos+rangeAvgNeg;
    
    if wtSum~=0
        dopplerAvg=dopplerAvg/wtSum;
		rangeAvg=rangeAvg/wtSum;
    end
    if wtSumPos~=0
        dopplerAvgPos=dopplerAvgPos/wtSumPos;
		rangeAvgPos=rangeAvgPos/wtSumPos;
    end
    if wtSumNeg~=0
        dopplerAvgNeg=dopplerAvgNeg/wtSumNeg;
		rangeAvgNeg=rangeAvgNeg/wtSumNeg;
    end
    
    % storing the range and doppler index corresponding to the maximum in the predetection matrix
	% 存储能量最大的点
	maxIndices(1)=maxRangeIdx;
	maxIndices(2)=maxDopplerIdx;
    
    % 9个信息
    gestureMetrics(1)=11;
	gestureMetrics(2)=numDetections;
	gestureMetrics(3)=fix(dopplerAvg);
	gestureMetrics(4)=fix(rangeAvg); %
	gestureMetrics(5)=numDetectionsPos;  %
	gestureMetrics(6)=fix(dopplerAvgPos);  %
	gestureMetrics(7)=fix(rangeAvgPos);
	gestureMetrics(8)=numDetectionsNeg;  %
	gestureMetrics(9)=fix(dopplerAvgNeg);  %
	gestureMetrics(10)=fix(rangeAvgNeg);
    
end