clear
count = 0;
RT2 = 100000.;
VT2 = -6000.;
BETA = 500.;
RT2H = 100025.;
VT2H = -6150.;
BETAH = 800.;
ORDER = 3;
TS = .05;
TF = 30.;
Q33 = 0./TF;
T = 0.;
S = 0.;
H = .001;
SIGNOISE = sqrt(500.);
PHI = zeros(ORDER,ORDER);
P = zeros(ORDER,ORDER);
Q = zeros(ORDER,ORDER);
IDN = eye(ORDER);
P(1,1) = SIGNOISE*SIGNOISE;
P(2,2) = 20000.;
P(3,3) = 300.^2;
HMAT = zeros(1,ORDER);
HMAT(1,1) = 1.;

while RT2 >= 0.
    RT2OLD = RT2;
    VT2OLD = VT2;
    STEP = 1;
    FLAG = 0;
    while STEP <= 1
        if FLAG == 1
            RT2 = RT2+H*RT2D;
            VT2 = VT2+H*VT2D;
            T = T+H;
            STEP = 2;
        end
        RT2D = VT2;
        VT2D = .0034*32.2*VT2*VT2*exp(-RT2/22000.)/(2.*BETA)-32.2;
        FLAG = 1;
    end
    FLAG = 0;
    % 
    % true target simulation 
    RT2 = .5*(RT2OLD+RT2+H*RT2D);
    VT2 = .5*(VT2OLD+VT2+H*VT2D);
    
    
    S = S+H;
    if S >= (TS-.00001)
        S = 0.;
        
        % preset 
        RHOH = .0034*exp(-RT2H/22000.);
        F21 = -32.2*RHOH*VT2H*VT2H/(2.*22000.*BETAH);
        F22 = RHOH*32.2*VT2H/BETAH;
        F23 = -RHOH*32.2*VT2H*VT2H/(2.*BETAH*BETAH);
        PHI(1,1) = 1.;
        PHI(1,2) = TS;
        PHI(2,1) = F21*TS;
        PHI(2,2) = 1.+F22*TS;
        PHI(2,3) = F23*TS;
        PHI(3,3) = 1.;
        Q(2,2) = F23*F23*Q33*TS*TS*TS/3.;
        Q(2,3) = F23*Q33*TS*TS/2.;
        Q(3,2) = Q(2,3);
        Q(3,3) = Q33*TS;
        
        % 
        % covariance matrix 
        %%
        % predict:
        M = PHI*P*PHI'+Q;
        
        % correct 
        HMHT = HMAT*M*HMAT';                    % riccati 1
        HMHTR = HMHT(1,1)+SIGNOISE*SIGNOISE;    % riccati 2
        HMHTRINV = 1./HMHTR;                    % riccati 3
        MHT = M*HMAT';                          % riccati 4
        for I = 1:ORDER
            % 
            % kalman gains 
            %% 
            GAIN(I,1) = MHT(I,1)*HMHTRINV;      % riccati 5 and final. 
        end
        P = (IDN-GAIN*HMAT)*M;
        
        
        XNOISE = SIGNOISE * randn; %  .86217;% gaussc7(SIGNOISE);
        
        %
        % state equations \ predict?  
        %%
        RT2DB = VT2H;
        VT2DB = .0034*32.2*VT2H*VT2H*exp(-RT2H/22000.)/(2.*BETAH)-32.2;
        
        RES = RT2+XNOISE-(RT2H+RT2DB*TS);
        %
        % integration with kalman gains 
        %%        
        RT2H = RT2H+RT2DB*TS+GAIN(1,1)*RES;
        VT2H = VT2H+VT2DB*TS+GAIN(2,1)*RES;
        BETAH = BETAH+GAIN(3,1)*RES;
        
        
        ERRY = RT2-RT2H;
        SP11 = sqrt(P(1,1));
        ERRV = VT2-VT2H;
        SP22 = sqrt(P(2,2));
        ERRBETA = BETA-BETAH;
        SP33 = sqrt(P(3,3));
        RT2K = RT2/1000.;
        count = count+1;
        ArrayT(count) = T;
        ArrayRT2(count) = RT2;
        ArrayRT2H(count) = RT2H;
        ArrayRT2K(count) = RT2K;
        ArrayVT2(count) = VT2;
        ArrayVT2H(count) = VT2H;
        ArrayBETA(count) = BETA;
        ArrayBETAH(count) = BETAH;
        ArrayERRBETA(count) = ERRBETA;
        ArraySP33(count) = SP33;
    end
end


figure
plot(ArrayRT2K,ArrayBETA,ArrayRT2K,ArrayBETAH),grid
xlabel('Altitude (Kft)')
ylabel('Ballistic Coefficient (Lb/Ft^2)')
axis([0 100 0 1000])
figure
plot(ArrayRT2K,ArraySP33,ArrayRT2K,-ArraySP33,ArrayRT2K,...
ArrayERRBETA),grid
xlabel('Altitude (Kft)')
ylabel('Error in Ballistic Coefficient (Lb/Ft^2)')
clc
output = [ArrayT',ArrayRT2K',ArrayRT2',ArrayRT2H',ArrayVT2',...
ArrayVT2H',ArrayBETA',ArrayBETAH'];
save datfil.txt output /ascii
disp 'simulation finished'
