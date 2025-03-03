%% solve the problem 
%  --------------------------------------------------------------------------------------------------------
%  ALM algorithm: solve the following problem
%  Original Problem: 
%     min_X, B, S, U, V  sum_i=1^2 |D_i(U)|_{q_i}^{q_i} + lambda * |B|_1 + beta * sum_b=1^P |S_b|_*
%     s.t. ||Y - X - S - B||_F^2 <= xi, rank(S_b) <= r_s
%          D3(X) = U x_3 V, V' * V = I
%     where,D is difference operator
%  Auxiliary Problem: 
%     min_X, B, S, U, V, Z_i  sum_i=1^2 ||Z_i||_{q_i}^{q_i} + lambda * ||B||_1 + beta * sum_b=1^P ||S_b||_*
%     s.t. ||Y - X - S - B||_F^2 <= xi, rank(S_b) <= r_s
%     D3(X) = U x_3 V, V' * V = I
%     D_i(U) = Z_i, i=1,2
%  The lagrange function is :                                          
%     L(X, B, S, U, V, Z_i, Gamma_i) 
%        = sum_i=1^2 ||Z_i||_{q_i}^{q_i} + lambda * ||B||_1 + beta * sum_b=1^P ||S_b||_* 
%            + mu/2 * sum_i=1^2 ||D_i(U) - Z_i||_F^2 + sum_i=1^2 <Gamma_i, D_i(U) - Z_i> 
%            + mu/2 * ||D3(X) - U x_3 V||_F^2 + <Gamma_3, D3(X) - U x_3 V>
%            + mu/2 * ||Y - X - B - S||_F^2 + <Gamma_4, Y - X - B - S> 
%     s.t. rank(S_b) <= r_s, V' * V = I
%  --------------------------------------------------------------------------------------------------------
function [DenoisedHSI] = HLTVSG(Nhsi, beta, lambda, r)
%% Input Varaible
%   Nhsi   --- Polluted HSI data
%   beta   --- trade-off parameter of stripe noise
%   lambda --- trade-off parameter of sparse noise
%   r      --- rank of needed-to-be recovered Tensor along mode-3
%% Output Variable
%   DenoisedHSI --- Denoised HSI
%% Initializing alm variables
tol     = 1e-5;
maxIter = 50;
rho     = 1.3;
max_mu  = 1e6;
mu      = 0.1;
%% initialize
[M,N,p] = size(Nhsi);
sizeD = size(Nhsi);
Y       =  Nhsi;
sv = 10;
r2 = 1;
eps =1e-4;
%% FFT setting
h       = M;
w       = N;
d       = r;
sizeU   = [h,w,d];
%% 
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z   = ( abs(psf2otf([+1, -1], [w,p,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
determX =  Eny_z;
determ  = Eny_x + Eny_y;
%% Initializing optimization variables
DX      = diff_z(Nhsi,sizeD);
Temp_DX = reshape(DX,M*N,p);
[~,~,v] = svd(Temp_DX,'econ');
V       = v(:,1:r);
U       = double(ttm(tensor(DX),V',3));
B       = zeros(M,N,p);  %Sparse
S       = zeros(M,N,p);  %Stripe
M1      = zeros(M,N,r);
M2      = zeros(M,N,r);
M3      = zeros(M,N,p);
M4      = zeros(M,N,p);
%% main loop
iter    = 0;
while iter<maxIter
    iter      = iter + 1;  
    if iter == 1
        preX  = Y;
    else
        preX  = X;
    end
    %% -Update Z1,Z2 
    TempZ1   = diff_x(U,sizeU)+M1/mu;
    Z1       = solve_Lp(TempZ1, 1/mu, 1/2); 
    TempZ2   = diff_y(U,sizeU)+M2/mu;
    Z2       = solve_Lp(TempZ2, 1/mu,1/2);
    %% -Update X
    diffT_pX = diff_zT(double(ttm(tensor(U),V,3))-M3/mu,sizeD);
    tempX    = Y-S-B+M4/mu;
    numer1X  = diffT_pX +tempX;
    X        = real( ifftn( fftn(numer1X) ./ (determX + 1 + eps) ) );
    %% -Update U
    diffT_p  = diff_xT(Z1-M1/mu,sizeU)+diff_yT(Z2-M2/mu,sizeU);
    temp1 = diff_z(X,sizeD);
    temp     = double(ttm(tensor(temp1 + M3/mu),V',3));
    numer1   = diffT_p +temp;
    U        = real( ifftn( fftn(numer1) ./ (determ + 1 + eps) ) );
    %% -Update V
    Temp1    = Unfold(temp1+M3/mu,sizeD,3);
    Temp2    = Unfold(U,[M,N,r],3)';
    T        = Temp1*Temp2;
    [u,~,v]  = svd(T,'econ');
    V        = u*v';
    %% -Update S
    TempS    = Y - X - B + M4/mu;
    for i =1:p
        temp2 = TempS(:,:,i);
        if  choosvd(p,sv) ==1
            [U1, sigma, V1] = lansvd(temp2, sv, 'L');
        else
            [U1,sigma,V1] = svd(temp2,'econ');
        end
        sigma = diag(sigma);
        svp = min(length(find(sigma>beta/(mu))),r2);
        if svp<sv
            sv = min(svp + 1, p);
        else
            sv = min(svp + round(0.05*p), p);
        end
        S(:,:,i) = U1(:, 1:svp) * diag(sigma(1:svp) - beta/(mu)) * V1(:, 1:svp)';
    end
    %% -Update B
    B    = prox_l1(Y-X-S+M4/mu,lambda/mu);
    %% -Update Multiplier
    leq1 = diff_x(U,sizeU)-Z1;
    leq2 = diff_y(U,sizeU)-Z2;
    leq3 = temp1 -double(ttm(tensor(U),V,3));
    leq4 = Y - X - S - B;
    stopC = norm(X(:)-preX(:),'fro') / norm(preX(:),'fro');
    if stopC<tol
        break;
    else
        M1 = M1 + mu*leq1;
        M2 = M2 + mu*leq2;
        M3 = M3 + mu*leq3;
        M4 = M4 + mu*leq4;
        mu = min(max_mu,mu*rho); 
    end
end
%% Denoised image
DenoisedHSI = X;
if p >= 150
    for i=1:p
    DenoisedHSI(:,:,i)=DenoisedHSI(:,:,i)+ mean(mean(S(:,:,i)));
    end
end
end