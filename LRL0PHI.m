% 程序还存在很多问题，主要原因是对matlab编程的不熟悉，不知道有些结构需要怎么实现。
% 所以在程序中需要老师指点的地方，都有注释。
% 下面总结一下大的几个问题：
% 1.一个m文件中多个子函数并存和调用的问题。
% 2.element结构体在matlab中的构造问题
% 3.从上面的第二个问题可以引出matlab中列表（list）类型或者近似list类型的数据结构可以替换它。
% 4.源代码中涉及到的迭代器的问题在很多的for循环中被使用，但是还没有找到好的解决方法或者替代方案。
% 5.代码还有很多地方需要验证和修改的地方。
clear;clc;

garma = 2.2;

% 在源码中element是一个结构体，里面的数据分别是整数w、整数w_mask、浮点类型 Y、浮点类型M_mean、Y_mean、I_mean
% 还有列表（list）类型G、N 和 图（map）类型 c
% 对于matlab中的结构体不是很了解，所以应该是要改的
element = struct('w',{},'w_mask',{},'Y',{},'M_mean',{},'Y_mean',{},'I_mean',{},G,{},N,{},c,{});


% 这里是对各个需要用到的矩阵的初始化变量
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

% 在m文件中的子函数调用定义引用
function funs=LRL0PHI;
funs.ffind=@ffind;
funs.g=@g;
funs.setParameters=@setParameters;
funs.sub_1=@sub_1;
funs.sub_2=@sub_2;
funs.sub_3=@sub_3;
funs.compute=@compute;


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
    lambda_rank_ = lambda_rank; % 10;
    lambda_l0_ = lambda_l0;
    k_ = k;
end

function []=sub_1(k)
    epsilon = 1.0e-04;
    I = struct('int',{},'element',{});  %  map<int, element> I;
    H_ = size(U_);
    H_ = H_(2);
    W_ = size(U_);
    W_ = W_(1);
    for i=0:H_
        for j=0:W_
            index = i*W_+j;
            temp = element;
            temp.G = [temp.G index];
            temp.w = 1;
            temp.w_mask = impixel(i,j,mask_); % 从mask_矩阵中获取
            temp.I_mean = temp.w_mask * impixel(i,j,I_);
            temp.M_mean = impixel(i,j,M_);
            temp.Y_mean = impixel(i,j,Y_);
            temp.Y = (temp.w*rho_*(temp.M_mean - temp.Y_mean) + 2*temp.w_mask*temp.I_mean)/(temp.w*rho_ + 2*temp.w_mask);
            
            if i ~= 0
                temp.N(end + 1) = (i-1)*W_+j;   % temp.N.push_back((i-1)*U_.cols+j); 源码中N是list类型（数组）
                temp.c(end + 1) = [(i-1)*W_+j,1];   %  temp.c.insert(pair<int,int>((i-1)*U_.cols+j,1)); 源码中c是map<int,int>类型。
            end
            
            if i ~= H_-1
                temp.N(end + 1) = (i+1)*W_+j;
                temp.c(end + 1) = [(i+1)*W_+j,1];
            end
            
            if j ~= 0
                temp.N(end + 1) = (i)*W_+j-1;
                temp.c(end + 1) = [(i)*W_+j-1,1];
            end
            
            if j ~= W_-1
                temp.N(end + 1) = (i)*W_+j+1;
                temp.c(end + 1) = [(i)*W_+j+1,1];
            end
            % 在结构体I中插入index和temp，下面的语句肯定是不行的
            I(end + 1) = [index, temp];  % I.insert(pair<int,element>(index,temp));  
        end
    end
        % map<int, element>::iterator i;
        % for(i = I.begin(); i != I.end(); ++i)
        %   for(list<int>::iterator j = i->second.G.begin(); j != i->second.G.end(); ++j)
           %   U_.at<float>(*j/U_.cols, *j % U_.cols) = i->second.Y;
        beta = 0;
        iter = 0;
        while(1)
            % 没有迭代器的话，怎么指向开始和结尾，怎么在循环中递增 ？
            it = I.begin();% map<int,element>::iterator it = I.begin();
            while(it ~= I.end())
                i = it.first;  % i = it->first;
                 % for(list<int>::iterator j = I[i].N.begin(); j != I[i].N.end();)
                 % 这个for 也是需要改的  
                 for(j = I(i).N.begin():I(i).N.end())
                    value1 = 0;
                    value2 = 0;
                    value3 = 0;
                    temp1 = I(i);
                    temp2 = I(j);
                    value1 = (temp1.w * temp1.w_mask * rho_ * pow(temp1.I_mean - temp1.M_mean + temp1.Y_mean,2) / (temp1.w * rho_ + 2 * temp1.w_mask))+(temp2.w * temp2.w_mask * rho_ * pow(temp2.I_mean - temp2.M_mean + temp2.Y_mean, 2) / (temp2.w * rho_ + 2 * temp2.w_mask))+(beta * ffind(I[i].c.begin(), I[i].c.end(), *j)->second);
                    X = 0;
                    X = (temp1.w*rho_*(temp1.M_mean-temp1.Y_mean) + 2*temp1.w_mask*temp1.I_mean+temp2.w*rho_*(temp2.M_mean-temp2.Y_mean) + 2*temp2.w_mask*temp2.I_mean)/(temp1.w*rho_ + 2*temp1.w_mask + temp2.w*rho_ + 2*temp2.w_mask);
                    value2 = (temp1.w*rho_/2.0 * pow(X - temp1.M_mean + temp1.Y_mean, 2))+(temp1.w_mask * pow(X - temp1.I_mean, 2))+(temp2.w*rho_/2.0f * pow(X - temp2.M_mean + temp2.Y_mean, 2))+(temp2.w_mask * pow(X - temp2.I_mean, 2));
                    Xi = 0; 
                    Xj = 0;
                    if(temp2.Y < temp1.Y)
                    
                        Xi = (temp1.w*rho_*(temp1.M_mean-temp1.Y_mean) + 2*temp1.w_mask*temp1.I_mean+temp2.w*rho_*(temp2.M_mean-temp2.Y_mean+1.0) + 2*temp2.w_mask*(temp2.I_mean+1.0f))/(temp1.w*rho_ + 2*temp1.w_mask + temp2.w*rho_ + 2*temp2.w_mask);
                        Xj = Xi - 1.0;
                    
                    elseif(temp2.Y >= temp1.Y)
                    
                        Xi = (temp1.w*rho_*(temp1.M_mean-temp1.Y_mean) + 2*temp1.w_mask*temp1.I_mean+temp2.w*rho_*(temp2.M_mean-temp2.Y_mean-1.0) + 2*temp2.w_mask*(temp2.I_mean-1.0f))/(temp1.w*rho_ + 2*temp1.w_mask + temp2.w*rho_ + 2*temp2.w_mask);
                        Xj = Xi + 1.0;
                     end
                        
                    
                    value3 = (temp1.w*rho_/2.0 * pow(Xi - temp1.M_mean + temp1.Y_mean, 2))+(temp1.w_mask * pow(Xi - temp1.I_mean, 2))+(temp2.w*rho_/2.0 * pow(Xj - temp2.M_mean + temp2.Y_mean, 2))+(temp2.w_mask * pow(Xj - temp2.I_mean, 2)) + k_*lambda_l0_;
                    if(value1 <= epsilon)
                        value1 = 0;
                    end
                        
                    if(value2 <= epsilon)
                        value2 = 0;
                    end
                        
                    if(value3 <= epsilon)
                        value3 = 0;
                    end
                    % value1 = value1 <= epsilon ? 0 : value1;
                    % value2 = value2 <= epsilon ? 0 : value2;
                    % value3 = value3 <= epsilon ? 0 : value3;

                   %  if(value3 < (value1 < value2 ? value1 : value2))  % value3 小于value1和value2中的最小值，相当于论文中的fA <= fB and fA <= fC
                    if(value1 < value2)
                        if(value3 < value1)
                            temp1.Y = Xi;
                            temp2.Y = Xj;
                            j = j + 1;
                        end
                    elseif(value3 < value2)
                            temp1.Y = Xi;
                            temp2.Y = Xj;
                            j = j + 1;
                    
                    % value2小于 value1和value3 中的最小值，相当于论文中的fB <= fA and fB <= fC
                    elseif(value1 < value3 || value2<=value3)
                            temp_value = j; % int temp_value = *j;
                            t1 = I(i).G;
                            t2 = I(j).G;
                            I(i).Y = X;
                            I(i).Y_mean = (I(i).Y_mean*I(i).w + I(j).Y_mean*I(j).w)/(I(i).w + I(j).w);
                            I(i).M_mean = (I(i).M_mean*I(i).w + I(j).M_mean*I(j).w)/(I(i).w + I(j).w);
                            
                            if(I(i).w_mask == 0 && I(j).w_mask == 0)
                                I(i).I_mean = 0;
                            else
                                I(i).I_mean = (I(i).I_mean*I(i).w_mask + I(j).I_mean*I(j).w_mask)/(I(i).w_mask + I(j).w_mask);
                            end
                            I(i).w = I(i).w + I(j).w;
                            I(i).w_mask = I(i).w_mask + I(j).w_mask;
                            I(i).G.merge(I(j).G);   % I[i].G.merge(I[*j].G);
                            I(i).c.erase(temp_value);
                            % 列表（）list的删除方法erase(),这里需要找一个适合的删除方法
                            j = I(i).N.erase(j);
                            % 数组迭代器，在matlab中有没有可以使用的迭代器？
                            % 如果没有的话，怎么改造
                            list<int>::iterator k;   
                            % for(k = I[temp_value].N.begin(); k != I[temp_value].N.end(); ++k)
                            % 这个for循环也是需要改的
                            for(k = I(temp_value).N.begin():I(temp_value).N.end())                          
                                if(k == i)
                                    continue;
                                end
                                if(ffind(I(i).N.begin(), I(i).N.end(), k)~=I(i).N.end())
                                % ffind函数的修改，需要符合两种情况
                                % 这里的  ->second应该是指迭代器下一个，需要确认
                                    ffind(I(i).c.begin(),I(i).c.end(),k)->second +=
                                    ffind(I(temp_value).c.begin(),I(temp_value).c.end(),k)->second;
                                    ffind(I(k).c.begin(),I(k).c.end(),i)->second =
                                    ffind(I(i).c.begin(),I(i).c.end(),k)->second;
                             
                                else
                                    I(i).N.push_back(k);
                                    temp = I(i).N;
                                    I(k).N.push_back(i);
                                    % 这里和前面一样，insert需要改
                                    I(i).c.insert(pair<int,int>(k,ffind(I(temp_value).c.begin(),I(temp_value).c.end(),k)->second));
                                    I(k).c.insert(pair<int,int>(i,ffind(I(temp_value).c.begin(),I(temp_value).c.end(),k)->second));
                                end
                                % 这里删除I中的方法，c是map类型的erase和N是list类型的remove
                                I(k).c.erase(temp_value);
                                I(k).N.remove(temp_value);
                                k = k+1;
                            
                            I.erase(temp_value);
                            % 这里指向的是it迭代器的下一个，next()函数不知道可不可以
                            it = it.next();
                            break;
                            end
                    else
                        j = j+1;
                    end
                 end
               if(it == I.end())
                   break;
               end
                   it = it+1;
            end
            beta = g(++iter, K, lambda_l0_);
            if(beta > lambda_l0_)
                break;
            end
        end
        % 这里的for循环也要改
        for(map<int, element>::iterator i = I.begin(); i != I.end(); ++i)
            for(list<int>::iterator j = i->second.G.begin(); j != i->second.G.end(); ++j)
                U_.at<float>(*j/U_.cols, *j % U_.cols) = i->second.Y;
                    temp = zeros(H_, W_);
        U_.convertTo(temp, CV_8UC1);
        string name;
        name = 'ans_'+to_string(iters1)+'.png';
        imwrite('./result/'+name, temp);
            end
            
    
            
function []=sub_2()
    A = U_ + Y_;
    mask = 255*ones(H_, W_);
    lambda = rho_ / 2.0 / lambda_rank_ / alpha_;
    
    M_ = TNNR(At, mask, 9, 9, lambda);  % 这里需要引入TNNR函数
end
    
function []=sub_3()
    Y_ = Y_ + (U_ - M_);    
end

function U_=compute(int K, int max_iter, string path, Mat &original, string path1)
    % 这里的for循环需要改
  for(int iter = 0; iter < max_iter; ++iter)
    {
      disp('=====================');
      disp('Iter = ', (iter + 1));
      iters1 = iter;
      sub_1(K);
      disp('Piecewise L0 done'); %%  分段 L0  Done
      sub_2();
      disp('LR done');
      sub_3();
        
        
      % 这一句只是关系的转化和矩阵单位的转化
      U_.convertTo(output, CV_8UC1);
      string newpath = path + to_string(iter+1);

      imwrite(newpath + '.png', output);
      disp('result in ', newpath + '.png');
      mask = 255.0 * mask_;
      % 这里需要引入PSNR函数
      psnr = PSNR(original, U_, mask);

      string filePath = path1 + '/';
      filePath = filePath + to_string(iter + 1);
      filePath = filePath + '.txt';
      % 下面是对结果的保存和展示
      fstream file1(filePath, ios::out);
      % 将psnr结果保存在file1文件中
      file1 << psnr << endl;
      disp('PSNR = '+psnr);
      file1.close();
      disp('result psnr in ' + filePath);

      if( iter >= 1 && norm(U_, U_last_) / norm(U_last_) < 1e-03 )
        break;
      end
      disp('relative error = '+norm(U_, U_last_) / norm(U_last_));
      U_.copyTo(U_last_)
  end;
     
    

    



