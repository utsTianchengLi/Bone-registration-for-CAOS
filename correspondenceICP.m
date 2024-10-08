function [Result,Cov] = correspondenceICP (Target, Observation)


n=length(Observation);

% figure
% plot3(Target(1,:),Target(2,:),Target(3,:),'g*','Markersize',10);
% hold on;
% plot3(Observation(1,:),Observation(2,:),Observation(3,:),'ro','Markersize',10);
% hold on;
% title('Input');
% hold on;

e_sum=[];
delta_sum=[];
X=ones(6,1)*0;   %[dPhi;dx;dy;dz]
XX=ones(6,1)*0;  %[Alpha,Beta,Gamma,x,y,z]
RR=eul2rotm(XX(1:3)');
TT=XX(4:6);

Source_initial=RR*Observation+TT;
% figure
% plot3(Target(1,:),Target(2,:),Target(3,:),'g*','Markersize',10);
% hold on;
% plot3(Source_initial(1,:),Source_initial(2,:),Source_initial(3,:),'ro','Markersize',10);
% hold on;
% title('Initial Value'); xlabel('X');  ylabel('Y');  zlabel('Z');   
% hold on;


for k=1:100 %Iteration
    dMdP_all=[];
    J=zeros(3*n,6);
    F=zeros(3*n,1);
    H=zeros(6,6);
    B=zeros(6,1);

    for i=1:n
        Obs_i=[Observation(1,i); Observation(2,i); Observation(3,i)];
        Tar_i=[Target(1,i); Target(2,i); Target(3,i)];

        F(3*i-2:3*i)=expm(VectortoSkew(X(1:3)))*RR*Obs_i+TT + X(4:6) - Tar_i;
        J(3*i-2:3*i,:)=[-VectortoSkew(RR*Obs_i),eye(3)];
       
        H=H+J(3*i-2:3*i,:)'*J(3*i-2:3*i,:);
        B=B-J(3*i-2:3*i,:)'*F(3*i-2:3*i);

    end

    error=sum(abs(F));
    delta=H\B;

    X=delta; 
    RR=expm(VectortoSkew(X(1:3)))*RR;
    TT=TT+X(4:6);

    X=ones(6,1)*0;

    delta_sq = sum(delta.^2);
    e_sum=[e_sum;error];

    if delta_sq <= 1e-8
        break;
    end

end
R_result=RR;
T_result=TT;

%% Result
% GT=[theta_gt,T_gt']
Result=[rotm2eul(R_result),T_result'];
Cov=inv(H);
% Cov=diag(diag(inv(H)));

Obs_Result=R_result*Observation+T_result;

% figure
% plot3(Target(1,:),Target(2,:),Target(3,:),'g*','Markersize',10);
% hold on;
% plot3(Obs_Result(1,:),Obs_Result(2,:),Obs_Result(3,:),'ro','Markersize',10);
% hold on;
% title('Result');
% hold on;


end