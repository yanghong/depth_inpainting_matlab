clear;clc;

garma = 2.2;

element = struct('w',{},'w_mask',{},'Y',{},'M_mean',{},'Y_mean',{},'I_mean',{},G,{},N,{},c,{});
I = zeros(50,49);
U_ = zeros(1,1);

function LRL0PHI=LRL0PHI(I, mask)
    W_ = size(I); % 获得矩阵的形状
    W_ = W_(2);  %矩阵的宽度（列数）
    H_ = size(I); % 获得矩阵的形状
    H_ = H_(1);  %矩阵的高度（行数）
    I_ = zeros(H_, W_); %源程序中的converTo是不是多余？可能是数据类型转化的需要
    mask_ = zeros(H_, W_);
    
    I_ = ones(H_, W_);
    U_ = zeros(H_, W_);
    U_last_ = zeros(H_, W_);
    M_ = zeros(H_, W_);
    Y_ = zeros(H_, W_);
    
    kernelx_plus_ = [0.0,-1.0,1.0];   % 定义形状为（1，3）的矩阵
    kernelx_minus_ = [-1.0,1.0,0.0];   % 定义形状为（1，3）的矩阵
    kernely_plus_ = [0.0;-1.0;1.0];   % 定义形状为（3，1）的矩阵
    kernely_minus_ = [-1.0;1.0;0.0];   % 定义形状为（3，1）的矩阵
end

function funs=LRL0PHI;
funs.ffind=@ffind;
funs.g=@g;
funs.setParameters=@setParameters;

%在c++源码中第一个ffind函数的第一个参数start是list类型的迭代器（相当于数组）
function p=ffind(start, end_, value)
    p = start;
    while(p ~= end_)
        if p == value
            break
        ++p;  % 指向p下一个值
        end
    end
end


%在c++源码中第二个ffind函数的第一个参数start是list类型的迭代器（相当于数组）
function p=ffind(start, end_, value)
    p = start;
    while(p ~= end_)
        if p == value
            break
        ++p;  % 指向p下一个值
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
        
        
    

        
     
    

    



