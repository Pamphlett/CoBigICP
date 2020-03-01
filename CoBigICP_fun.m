function [ R, T, Flag ] = CoBigICP_fun(Md, Mo, MovData, RefData, MovInfo, RefInfo, Tf0, DistThr, AngThr, sigma_times)
if nargin == 0
    clc; close all;
    RefData = readCloudCsv('Hokuyo_0.csv' , 0.9)';
    MovData = readCloudCsv('Hokuyo_1.csv' , 0.9)';
    showInitVal = 1;
    if showInitVal
        figure;
        hold on;
        grid on;
        pcshow(RefData(1:3, :)', 'g');
        pcshow(MovData(1:3, :)', 'b');
        title('Initial Value');
    end
    if length(RefData) >= length(MovData)
        Tf0 = [ eye(3) mean(RefData(:, 1:length(MovData))-MovData, 2); 0 0 0 1];
    else
        Tf0 = [ eye(3) mean(RefData-MovData(:, 1:length(RefData)), 2); 0 0 0 1];
    end
    DistThr = Inf;
    AngThr = cosd(30);
    MovInfo = CalNormalsFun( MovData );
    RefInfo = CalNormalsFun( RefData );
    Md = createns(RefData');
    Mo = createns(MovData');
end

Flag = 1;
sigma_itr = 31.62;
DecayPram = 0.97;%1.03;
R0 = Tf0(1:3, 1:3); T0 = Tf0(1:3, end);
R = R0; T = T0;
dR = eye(3);
dT = zeros(3, 1);
U0 = cat(2, RefInfo(:).U );
U1 = cat(2, MovInfo(:).U );
Norm_Ref = cat( 2, RefInfo(:).normal );
Norm_Mov = cat( 2, MovInfo(:).normal );
epsilon = 1e-3;
MaxIter = 50;
JArray = [];
for nIter = 1 : 1 : MaxIter
    %%%%%%%%%%
    T_curr = [R T; 0 0 0 1];
    T_curr_inv = inv(T_curr);
    R_inv = T_curr_inv(1:3, 1:3); T_inv = T_curr_inv(1:3, end);
    AftData = Loc2Glo( MovData, R', T );   % apply transformation to move points.
    refAftData = Loc2Glo( RefData, R_inv', T_inv );   % apply transformation to ref points.
    [NNIdx, DD] = knnsearch( Md, AftData' ); % establish correspondence.
    % bidirectional correspondence
    [NNIdx_inverse, ~] = knnsearch( Mo, refAftData' );
    bi_eff = find(vecnorm(AftData - AftData(:, NNIdx_inverse(NNIdx))) < 0.01);
    Angle = sum( Norm_Ref(:, NNIdx) .* (R * Norm_Mov) );
    EffIdx_sim = find( abs(Angle) > AngThr & DD' < DistThr );
    EffIdx = intersect(EffIdx_sim, bi_eff);
    MovIdx = EffIdx;
    RefIdx = NNIdx(EffIdx);
    
    tmpAft = AftData(:, MovIdx);
    tmpRef = RefData(:, RefIdx);
    
    %%%%%%%%%% after transformation, we need align AftData to RefData.
    %%%%%%%%%%% obtain rotation and translation via solving quadratic programming.
    USE_NOTATION = 0;
    if USE_NOTATION
        A = zeros(6, 6);
        b = zeros(6, 1);
        tmpDiff = tmpAft - tmpRef;
        J = 0.0;
        for i = 1 : 1 : size(tmpAft, 2)
            id0 = RefIdx(i);
            id1 = MovIdx(i);
            omega = MovInfo(id1).omega + R*RefInfo(id0).omega*R';
            H = [-SkewFun(tmpAft(:, i)) eye(3)];
            v = tmpDiff(:, i);
            deltaA = H'*omega*H;
            deltab = H'*omega*v;
            deltaJ = v'*omega*v;
            wi = - 1/ (sigma_itr^2) * exp(- deltaJ/(2*sigma_itr^2));
            A = A + wi * deltaA;
            b = b + wi * deltab;
            J = J + deltaJ;
        end
        x = -A\b;
        tmp = 0.5 * x(1:3);
        quat = [ 1 - norm(tmp) tmp'];
        dR = quat2rotm(quat);
        % dR = expm(SkewFun(x(1:3)));
        dT = x(4:end);
        R = dR * R;
        T = dR*T + dT;
    else
        [H, b, J] = CalHbCobig_Gabor(tmpAft,  tmpAft - tmpRef, R, RefInfo(RefIdx), MovInfo(MovIdx), sigma_itr);
        dx = -pinv(H)*b;
        dR = expm(SkewFun(dx(1:3)));
        dT = dx(4:6);
        R = R*dR;
        T = T + dT;
        bTest = 1;
    end
    Err = max( norm(dR - eye(3)), norm(dT) );
    JArray(end+1) = J/length(MovIdx);
    if length(JArray) >= 2
        if abs((JArray(end) - JArray(end-1))/ JArray(end-1) ) <= 1e-4
            break;
        end
    end
    str = sprintf('nIter = %02d/%02d, J = %.6f', nIter, MaxIter, JArray(end));
    disp(str);
    bTest = 1;
    sigma_itr = sigma_itr * DecayPram;
end
IS_SHOW = 1;
if IS_SHOW
    figure;
    hold on;
    grid on;
    plot(JArray, 'bx-');
    title('Convergence Curve');
    % Res = CalRes(MovData);
    figure;
    hold on;
    grid on;
    pcshow(RefData(1:3, :)', 'g');
    pcshow(AftData(1:3, :)', 'b');
    title('Result of CoBigICP');
end
Tf_est = [R T; 0 0 0 1];
% Tf_gt
bTest = 1;
end