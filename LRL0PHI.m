clear;clc;

garma = 2.2;

element = struct('w',{},'w_mask',{},'Y',{},'M_mean',{},'Y_mean',{},'I_mean',{},G,{},N,{},c,{});
I = zeros(50,49);
U_ = zeros(1,1);

function LRL0PHI=LRL0PHI(I, mask)
    W_ = size(I); % ��þ������״
    W_ = W_(2);  %����Ŀ�ȣ�������
    H_ = size(I); % ��þ������״
    H_ = H_(1);  %����ĸ߶ȣ�������
    I_ = zeros(H_, W_); %Դ�����е�converTo�ǲ��Ƕ��ࣿ��������������ת������Ҫ
    mask_ = zeros(H_, W_);
    
    I_ = ones(H_, W_);
    U_ = zeros(H_, W_);
    U_last_ = zeros(H_, W_);
    M_ = zeros(H_, W_);
    Y_ = zeros(H_, W_);
    
    kernelx_plus_ = [0.0,-1.0,1.0];   % ������״Ϊ��1��3���ľ���
    kernelx_minus_ = [-1.0,1.0,0.0];   % ������״Ϊ��1��3���ľ���
    kernely_plus_ = [0.0;-1.0;1.0];   % ������״Ϊ��3��1���ľ���
    kernely_minus_ = [-1.0;1.0;0.0];   % ������״Ϊ��3��1���ľ���
end

function funs=LRL0PHI;
funs.ffind=@ffind;
funs.g=@g;
funs.setParameters=@setParameters;

%��c++Դ���е�һ��ffind�����ĵ�һ������start��list���͵ĵ��������൱�����飩
function p=ffind(start, end_, value)
    p = start;
    while(p ~= end_)
        if p == value
            break
        ++p;  % ָ��p��һ��ֵ
        end
    end
end


%��c++Դ���еڶ���ffind�����ĵ�һ������start��list���͵ĵ��������൱�����飩
function p=ffind(start, end_, value)
    p = start;
    while(p ~= end_)
        if p == value
            break
        ++p;  % ָ��p��һ��ֵ
        end
    end
end
    

function ans_=g(iter, k, l)
    ans_ = pow((iter/k),garma) * l;
end

function []=setParameters(rho, dt, lambda_l0, lambda_rank, k)
    rho_ = rho;
    dt_ = dt;
    h_ = 1.0;
    alpha_ = 1.0;
    lambda_rank_ = lambda_rank; //10;
    lambda_l0_ = lambda_l0;
    k_ = k;
end

function []=sub_1(k)
    epsilon = 1.0e-04;
    I = struct('int',{},'element',{});  %  map<int, element> I;
    for i=0:H_
        
        
    

        
     
    

    



