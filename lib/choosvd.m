function y = choosvd( n, d)   % n�ǵ�ǰ���ε�λ��  dΪ10
% %  if  choosvd(p,sv) ==1
if n <= 100 
    if d / n <= 0.02      % ��С10/100 = 0.1     0
        y = 1;
    else
        y = 0;
    end
elseif n <= 200                     
    if d / n <= 0.06      % ��С10/200 = 0.05     167����  
        y = 1;
    else
        y = 0;
    end
elseif n <= 300           % ��С10/300 = 0.033    0
    if d / n <= 0.26
        y = 1;
    else
        y = 0;
    end
elseif n <= 400          % ��С10/400 = 0.025     0  
    if d / n <= 0.28
        y = 1;
    else
        y = 0;
    end
elseif n <= 500         % ��С10/500 = 0.02       30����
    if d / n <= 0.34
        y = 1;
    else
        y = 0;
    end
else
    if d / n <= 0.38
        y = 1;
    else
        y = 0;
    end
end




