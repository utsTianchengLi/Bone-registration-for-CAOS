clear all;
close all;


feaEM = importdata("real_hip_fea.csv");
keyEM= feaEM.data(4:13,4:6);  % for calculate target error
feaEM = feaEM.data(1:3,4:6);  % for get initial value


Contour_Partial_all = importdata("hip1.csv"); % for registration
Contour_Partial_all = Contour_Partial_all.data(:,4:6);
Contour_Partial = Contour_Partial_all';
nz=[1:1:length(Contour_Partial)];
Contour_Partial = Contour_Partial(:,nz);
Partial_new=Contour_Partial;

load('Gsdf_hip.mat');
load('mesh_hip.mat');

n=length(Contour_Partial);

figure
% plot3(Contour_FS(1,:),Contour_FS(2,:),Contour_FS(3,:),'g.','Markersize',5);
Draw(ver, tri, [0.8 0.8 1.0],[],1)
hold on;
plot3(Contour_Partial(1,:),Contour_Partial(2,:),Contour_Partial(3,:),'g.','Markersize',5);
hold on;
title('Input');
hold on;


xyz_bound=[min(min(min(XX))); max(max(max(XX))); min(min(min(YY))); max(max(max(YY))); min(min(min(ZZ))); max(max(max(ZZ)))];
DTgrid = griddedInterpolant(XX,YY,ZZ,DT);
GXgrid = griddedInterpolant(XX,YY,ZZ,FX);
GYgrid = griddedInterpolant(XX,YY,ZZ,FY);
GZgrid = griddedInterpolant(XX,YY,ZZ,FZ);

GXgrid_bound = griddedInterpolant(XX,YY,ZZ,FX,'nearest');
GYgrid_bound = griddedInterpolant(XX,YY,ZZ,FY,'nearest');
GZgrid_bound = griddedInterpolant(XX,YY,ZZ,FY,'nearest');

dgrid = [GXgrid(Contour_Partial(1,:),Contour_Partial(2,:),Contour_Partial(3,:));...
    GYgrid(Contour_Partial(1,:),Contour_Partial(2,:),Contour_Partial(3,:));...
    GZgrid(Contour_Partial(1,:),Contour_Partial(2,:),Contour_Partial(3,:))]';
dError = [DTgrid(Contour_Partial(1,:),Contour_Partial(2,:),Contour_Partial(3,:))]';


% 3 bony landmarks for initial guess
Tar_fea=[21.7512,   27.1817,   62.2995;
         20.2518,   2.6655,    29.3579;
         34.1895,   14.4357,  -0.0423];  %hip landmarks from CT

[XX0,Cov_icp]=correspondenceICP(Tar_fea,feaEM');
R_icp=eul2rotm(XX0(1:3)); T_icp=XX0(4:6)';
Contour_initial_all=R_icp*Contour_Partial_all'+T_icp;

figure
Draw(ver, tri, [0.8 0.8 1.0],[],1)
hold on;
plot3(Contour_initial_all(1,:),Contour_initial_all(2,:),Contour_initial_all(3,:),'g.','Markersize',5);
hold on;
title('Initial Guess');
hold on;
% set(gca,'Visible','off');
% set(gcf,'color','w');

%% M estimator - Cauchy Bootstrapper: input (initial estimate, zeros(1,6)), output (bootstrapper’s initial estimate)
% XX0=zeros(1,6);
RR0=eul2rotm(XX0(1:3)); TT0=XX0(4:6)';
X=ones(6,1)*0;   %[dPhi;dx;dy;dz] for lie group
alpha = 1; %  1 for Cauchy, 2 for Geman-McClure, 0.5 is similiar to Huber

id=[1:n]';
for outK=1:10 % iteration for move outlier

    for Mk=1:100 % iteration for M estimator
        r=zeros(n,1);
        w=zeros(n,1);
        F=zeros(n,1);
        J=zeros(n,6);
        B=zeros(6,1);
        H=zeros(6,6);
        for Mi=1:n
            XYZ_i=[Partial_new(1,Mi); Partial_new(2,Mi); Partial_new(3,Mi)];
            XXYYZZ=expm(VectortoSkew(X(1:3)))*RR0*XYZ_i + TT0 + X(4:6);
            r(Mi)=norm(DTgrid(XXYYZZ(1),XXYYZZ(2),XXYYZZ(3)));
            w(Mi)=1/(1+r(Mi)^2)^alpha;
            F(Mi)=w(Mi)*DTgrid(XXYYZZ(1),XXYYZZ(2),XXYYZZ(3));
            dMdP=w(Mi)*[GXgrid(XXYYZZ(1),XXYYZZ(2),XXYYZZ(3)), GYgrid(XXYYZZ(1),XXYYZZ(2),XXYYZZ(3)), GZgrid(XXYYZZ(1),XXYYZZ(2),XXYYZZ(3))];
            dPdX=[-VectortoSkew(RR0*XYZ_i),eye(3)];
            J(Mi,:)=dMdP*dPdX;

            H=H+J(Mi,:)'*J(Mi,:);
            B=B-J(Mi,:)'*F(Mi);
        end

        delta=H\B;
        RR0=expm(VectortoSkew(delta(1:3)))*RR0;
        TT0=TT0+delta(4:6);
        X=zeros(6,1);
        delta_sq = sum(delta.^2);

        if delta_sq <= 1e-8
            break;
        end

    end

    threshold=0.1;
    outlier_id=find(w<threshold);
    id(outlier_id)=[];
    inlinear_id=find(w>threshold);


    if isempty(outlier_id)
        break;
    else
        Partial_new=Partial_new(:,inlinear_id);
        n=length(Partial_new);
    end

end
RR=RR0; TT=TT0;



%% Result
R_result=RR;
t_result=TT;
Result=[rotm2eul(R_result),t_result']
Contour_Result=R_result*Contour_Partial+t_result;
Target_EM=R_result*keyEM'+t_result;
figure
Draw(ver, tri, [0.8 0.8 1.0], [],1)
hold on;
plot3(Contour_Result(1,inlinear_id),Contour_Result(2,inlinear_id),Contour_Result(3,inlinear_id),'g.','Markersize',7);
hold on;
title('Result');
hold on;

Error_sdf=Dist_2PTs(ver,Contour_Result(:,inlinear_id)')

hip_keypoints=[-6.0453,   -6.5948,    79.8407;
               -51.1726,   17.6529,   71.8320;
               -61.7729,   36.5040,  -28.1030;
               -24.7480,   41.2975,  -39.9093;
                28.9991,   37.0443,  -31.9439;
                58.4133,   50.6013,  -16.3706;
                69.4652,  -57.3617,   7.9928;
                62.7382,  -10.4568,   19.0905;
                29.4459,   1.6819,    43.8377;
                27.6139,   1.9227,    14.2663]; %hip Target landmarks from CT
key_error= mean(vecnorm(Target_EM-hip_keypoints'))

