clear all;
close all;

load('z010_Tibia_L_mesh.mat');
load('z010_Tibia_L_partial.mat');
load('z010_Tibia_L_gSDF.mat');

ver=mesh.vertices;
cen=mean(ver);
ver=ver-cen;
tri=mesh.faces;

% Contour_Partial=Contour_Partial';
n=length(Contour_Partial);

figure
% plot3(Contour_FS(1,:),Contour_FS(2,:),Contour_FS(3,:),'g.','Markersize',5);
Draw(ver, tri, [0.8 0.8 1.0],[],1)
hold on;
plot3(Contour_Partial(1,:),Contour_Partial(2,:),Contour_Partial(3,:),'g.','Markersize',5);
hold on;
title('Initial');
hold on;
% set(gca,'Visible','off');
% set(gcf,'color','w');



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


%% M estimator - Cauchy Bootstrapper: input (initial estimate, zeros(1,6)), output (bootstrapper’s initial estimate)

XX0=zeros(1,6);
tic
RR0=eul2rotm(XX0(1:3)); TT0=XX0(4:6)';
X=ones(6,1)*0;   %[dPhi;dx;dy;dz] for lie group
alpha = 0.001; %  1 for Cauchy, 2 for Geman-McClure

id=[1:n]';

for outK=1:100 % iteration for move outlier


    for Mk=1:100 % iteration for M estimator
        r=zeros(n,1);
        w=zeros(n,1);
        F=zeros(n,1);
        J=zeros(n,6);
        B=zeros(6,1);
        H=zeros(6,6);
        for Mi=1:n
            XYZ_i=[Contour_Partial(1,Mi); Contour_Partial(2,Mi); Contour_Partial(3,Mi)];
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

        R_result=RR0;
        t_result=TT0;
        Contour_Result=R_result*Contour_Partial+t_result;
        % figure
        % % plot3(Contour_FS(1,:),Contour_FS(2,:),Contour_FS(3,:),'g.','Markersize',8);
        % Draw(ver, tri, [0.8 0.8 1.0],[],1)
        % hold on;
        % plot3(Contour_Result(1,:),Contour_Result(2,:),Contour_Result(3,:),'g.','Markersize',5);
        % hold on;
        % set(gca,'Visible','off');
        % set(gcf,'color','w');


        if delta_sq <= 1e-8
            break;
        end

    end

    threshold=0;
    outlier_id=find(w<threshold);
    id(outlier_id)=[];
    inlinear_id=find(w>threshold);


    if isempty(outlier_id)
        break;
    else
        Contour_Partial=Contour_Partial(:,inlinear_id);
        n=length(Contour_Partial);
    end

end
inlinear_id=id;
RR=RR0; TT=TT0;

toc

figure
% plot3(Contour_FS(1,:),Contour_FS(2,:),Contour_FS(3,:),'g.','Markersize',8);
Draw(ver, tri, [0.8 0.8 1.0],[],1)
hold on;
plot3(Contour_Result(1,:),Contour_Result(2,:),Contour_Result(3,:),'g.','Markersize',5);
hold on;
title('Result');
hold on;
% set(gca,'Visible','off');
% set(gcf,'color','w');

%% Result
T_gt = [35;-25;-50];
theta_gt=[0.1; 0.15; -0.05];
GT=[theta_gt;T_gt]';
R_result=RR;
t_result=TT;
Result=[rotm2eul(R_result),t_result'];


Contour_Result=R_result*Contour_Partial+t_result;
sdf=[Result-GT];
sdf_error=[ mean(abs(sdf(1:3)))/pi*180, mean(abs(sdf(4:6))),Dist_2PTs(ver,Contour_Result(:,inlinear_id)')]