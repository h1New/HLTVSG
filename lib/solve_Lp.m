function x = solve_Lp(y, lambda, p)
    % Modified by Dr. Weisheng Dong for 3D tensor inputs
    J     = 2;  % Number of iterations
    
        % Check if lambda is a constant or a matrix
    if numel(lambda) == 1
        % If lambda is a scalar, use it as a constant
        lambda = lambda * ones(size(y));  % Convert lambda to a matrix of the same size as y
    end
    
    
    tau   = (2*lambda.*(1-p)).^(1/(2-p)) + p*lambda.*(2*(1-p)*lambda).^((p-1)/(2-p));
    x     = zeros(size(y));  % Initialize output tensor
    [M, N, P] = size(y);     % Get the size of the input tensor y

    for i = 1:M
        for j = 1:N
            for k = 1:P
                % For each element (i,j,k), get the corresponding y and lambda values
                if abs(y(i,j,k)) > tau(i,j,k)
                    y0 = y(i,j,k);
                    lambda0 = lambda(i,j,k);
                    t = abs(y0);
                    for iteration = 1:J
                        t = abs(y0) - p*lambda0*(t)^(p-1);  % Apply the formula
                    end
                    x(i,j,k) = sign(y0) * t;  % Assign the result to the output tensor
                end
            end
        end
    end
end

