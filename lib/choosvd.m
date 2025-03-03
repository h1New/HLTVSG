function y = choosvd( n, d)   % n是当前波段的位置  d为10
% %  if  choosvd(p,sv) ==1
if n <= 100 
    if d / n <= 0.02      % 最小10/100 = 0.1     0
        y = 1;
    else
        y = 0;
    end
elseif n <= 200                     
    if d / n <= 0.06      % 最小10/200 = 0.05     167波段  
        y = 1;
    else
        y = 0;
    end
elseif n <= 300           % 最小10/300 = 0.033    0
    if d / n <= 0.26
        y = 1;
    else
        y = 0;
    end
elseif n <= 400          % 最小10/400 = 0.025     0  
    if d / n <= 0.28
        y = 1;
    else
        y = 0;
    end
elseif n <= 500         % 最小10/500 = 0.02       30波段
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




