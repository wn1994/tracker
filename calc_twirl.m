function [out] = calc_twirl(datavec)
GESTURE_PKT_SIZE_BYTES=26;
NUM_PKTS_TO_COLLECT =4;

% global numdetections
% global dopplerave
global rangeave
global angleval
global numdetectionspos
global doppleravepos
global numdetectionsneg
global doppleraveneg
% NN
L_training = 30;

% ��¼����
global xin_data
%% Run
if ~isempty(datavec)
    if datavec(1)==11 % check magic word  ��һλ��11
        %            GESTURE_PKT_SIZE_BYTES/2=13        NUM_PKTS_TO_COLLECT=4
        % �ĸ�оƬ��ÿ��оƬһ�ε����ݴ�СΪ13
        datavec1=reshape(datavec,GESTURE_PKT_SIZE_BYTES/2,NUM_PKTS_TO_COLLECT);
        
        rangeave=[rangeave(1+NUM_PKTS_TO_COLLECT:end) datavec1(4,:) ]; %rangeave
        angleval=[angleval(1+NUM_PKTS_TO_COLLECT:end) datavec1(11,:) ]; %angleval
        numdetectionspos=[numdetectionspos(1+NUM_PKTS_TO_COLLECT:end) datavec1(5,:) ]; %numdetections
        doppleravepos=[doppleravepos(1+NUM_PKTS_TO_COLLECT:end) datavec1(6,:) ]; %dopplerave
        numdetectionsneg=[numdetectionsneg(1+NUM_PKTS_TO_COLLECT:end) datavec1(8,:) ]; %numdetections
        doppleraveneg=[doppleraveneg(1+NUM_PKTS_TO_COLLECT:end) datavec1(9,:) ]; %dopplerave
        % format data for classifier   4 5 6 8 9
        xin=[numdetectionspos([end-L_training+1:end]) numdetectionsneg([end-L_training+1:end]) doppleravepos([end-L_training+1:end]) doppleraveneg([end-L_training+1:end]) rangeave([end-L_training+1:end]) angleval([end-L_training+1:end])];
        xin_data = [xin_data;xin];
        % *****classification*****
        % ͨ��nn�ֳ�1˫����2������3���� �����ַ���
        out = myNeuralNetworkFunction(xin');
    end
end
% save('test_data.mat','records');
% save('test_out.mat','result');