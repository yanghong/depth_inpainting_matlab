% ���򻹴��ںܶ����⣬��Ҫԭ���Ƕ�matlab��̵Ĳ���Ϥ����֪����Щ�ṹ��Ҫ��ôʵ�֡�
% �����ڳ�������Ҫ��ʦָ��ĵط�������ע�͡�
% �����ܽ�һ�´�ļ������⣺
% 1.һ��m�ļ��ж���Ӻ�������͵��õ����⡣
% 2.element�ṹ����matlab�еĹ�������
% 3.������ĵڶ��������������matlab���б�list�����ͻ��߽���list���͵����ݽṹ�����滻����
% 4.Դ�������漰���ĵ������������ںܶ��forѭ���б�ʹ�ã����ǻ�û���ҵ��õĽ�������������������
% 5.���뻹�кܶ�ط���Ҫ��֤���޸ĵĵط���
clear;clc;

garma = 2.2;

% ��Դ����element��һ���ṹ�壬��������ݷֱ�������w������w_mask���������� Y����������M_mean��Y_mean��I_mean
% �����б�list������G��N �� ͼ��map������ c
% ����matlab�еĽṹ�岻�Ǻ��˽⣬����Ӧ����Ҫ�ĵ�
element = struct('w',{},'w_mask',{},'Y',{},'M_mean',{},'Y_mean',{},'I_mean',{},G,{},N,{},c,{});


% �����ǶԸ�����Ҫ�õ��ľ���ĳ�ʼ������
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

% ��m�ļ��е��Ӻ������ö�������
function funs=LRL0PHI;
funs.ffind=@ffind;
funs.g=@g;
funs.setParameters=@setParameters;
funs.sub_1=@sub_1;
funs.sub_2=@sub_2;
funs.sub_3=@sub_3;
funs.compute=@compute;


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
            temp.w_mask = impixel(i,j,mask_); % ��mask_�����л�ȡ
            temp.I_mean = temp.w_mask * impixel(i,j,I_);
            temp.M_mean = impixel(i,j,M_);
            temp.Y_mean = impixel(i,j,Y_);
            temp.Y = (temp.w*rho_*(temp.M_mean - temp.Y_mean) + 2*temp.w_mask*temp.I_mean)/(temp.w*rho_ + 2*temp.w_mask);
            
            if i ~= 0
                temp.N(end + 1) = (i-1)*W_+j;   % temp.N.push_back((i-1)*U_.cols+j); Դ����N��list���ͣ����飩
                temp.c(end + 1) = [(i-1)*W_+j,1];   %  temp.c.insert(pair<int,int>((i-1)*U_.cols+j,1)); Դ����c��map<int,int>���͡�
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
            % �ڽṹ��I�в���index��temp����������϶��ǲ��е�
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
            % û�е������Ļ�����ôָ��ʼ�ͽ�β����ô��ѭ���е��� ��
            it = I.begin();% map<int,element>::iterator it = I.begin();
            while(it ~= I.end())
                i = it.first;  % i = it->first;
                 % for(list<int>::iterator j = I[i].N.begin(); j != I[i].N.end();)
                 % ���for Ҳ����Ҫ�ĵ�  
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

                   %  if(value3 < (value1 < value2 ? value1 : value2))  % value3 С��value1��value2�е���Сֵ���൱�������е�fA <= fB and fA <= fC
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
                    
                    % value2С�� value1��value3 �е���Сֵ���൱�������е�fB <= fA and fB <= fC
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
                            % �б���list��ɾ������erase(),������Ҫ��һ���ʺϵ�ɾ������
                            j = I(i).N.erase(j);
                            % �������������matlab����û�п���ʹ�õĵ�������
                            % ���û�еĻ�����ô����
                            list<int>::iterator k;   
                            % for(k = I[temp_value].N.begin(); k != I[temp_value].N.end(); ++k)
                            % ���forѭ��Ҳ����Ҫ�ĵ�
                            for(k = I(temp_value).N.begin():I(temp_value).N.end())                          
                                if(k == i)
                                    continue;
                                end
                                if(ffind(I(i).N.begin(), I(i).N.end(), k)~=I(i).N.end())
                                % ffind�������޸ģ���Ҫ�����������
                                % �����  ->secondӦ����ָ��������һ������Ҫȷ��
                                    ffind(I(i).c.begin(),I(i).c.end(),k)->second +=
                                    ffind(I(temp_value).c.begin(),I(temp_value).c.end(),k)->second;
                                    ffind(I(k).c.begin(),I(k).c.end(),i)->second =
                                    ffind(I(i).c.begin(),I(i).c.end(),k)->second;
                             
                                else
                                    I(i).N.push_back(k);
                                    temp = I(i).N;
                                    I(k).N.push_back(i);
                                    % �����ǰ��һ����insert��Ҫ��
                                    I(i).c.insert(pair<int,int>(k,ffind(I(temp_value).c.begin(),I(temp_value).c.end(),k)->second));
                                    I(k).c.insert(pair<int,int>(i,ffind(I(temp_value).c.begin(),I(temp_value).c.end(),k)->second));
                                end
                                % ����ɾ��I�еķ�����c��map���͵�erase��N��list���͵�remove
                                I(k).c.erase(temp_value);
                                I(k).N.remove(temp_value);
                                k = k+1;
                            
                            I.erase(temp_value);
                            % ����ָ�����it����������һ����next()������֪���ɲ�����
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
        % �����forѭ��ҲҪ��
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
    
    M_ = TNNR(At, mask, 9, 9, lambda);  % ������Ҫ����TNNR����
end
    
function []=sub_3()
    Y_ = Y_ + (U_ - M_);    
end

function U_=compute(int K, int max_iter, string path, Mat &original, string path1)
    % �����forѭ����Ҫ��
  for(int iter = 0; iter < max_iter; ++iter)
    {
      disp('=====================');
      disp('Iter = ', (iter + 1));
      iters1 = iter;
      sub_1(K);
      disp('Piecewise L0 done'); %%  �ֶ� L0  Done
      sub_2();
      disp('LR done');
      sub_3();
        
        
      % ��һ��ֻ�ǹ�ϵ��ת���;���λ��ת��
      U_.convertTo(output, CV_8UC1);
      string newpath = path + to_string(iter+1);

      imwrite(newpath + '.png', output);
      disp('result in ', newpath + '.png');
      mask = 255.0 * mask_;
      % ������Ҫ����PSNR����
      psnr = PSNR(original, U_, mask);

      string filePath = path1 + '/';
      filePath = filePath + to_string(iter + 1);
      filePath = filePath + '.txt';
      % �����ǶԽ���ı����չʾ
      fstream file1(filePath, ios::out);
      % ��psnr���������file1�ļ���
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
     
    

    



