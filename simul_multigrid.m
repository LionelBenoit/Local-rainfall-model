function[SimulatedRainTS]=simul_multigrid(Coord_Euler,CondiRainTS_gaussian,m,nb_ep, step_t,nb_max_cond)

%!!CondiRainTS_gaussian.t have to start at 0.

M_res_sim=zeros(length(Coord_Euler),nb_ep);

%generate multigrid simulation path
Path_to_sim=[1;nb_ep];
To_sim=2:nb_ep-1;
To_sim=To_sim';
iter=1;
while ~isempty(To_sim) 
    V_path_temp=[];
    sim_path_sort=sort(Path_to_sim);
    for i=2:length(sim_path_sort)
        
        ind_min=sim_path_sort(i-1)+1;
        ind_max=sim_path_sort(i)-1;
        if ind_min>ind_max
            Path_to_sim=[Path_to_sim;To_sim];
            To_sim=[];
            break;
        end
        bmin=find(To_sim==ind_min);
        bmax=find(To_sim==ind_max);
        temp=floor(median((To_sim(bmin:bmax))));
        rem=find(To_sim==temp);
        To_sim=[To_sim(1:rem-1);To_sim(rem+1:end)];
        V_path_temp=[V_path_temp;temp];
    end
    if ~isempty(V_path_temp)
        Path_to_sim=[Path_to_sim;V_path_temp];
    end
    iter=iter+1;
end

Simulated_pts=[];

for i=1:length(Path_to_sim)
    display(strcat('Simulcondi: random path=',num2str(i),'/',num2str(length(Path_to_sim))))
    %display(strcat('pt_num=',num2str(i)))
    %search conditioning points
    V_dist_t=abs(Simulated_pts-Path_to_sim(i));
    V_dist_tempo_sort=sort(V_dist_t);
    ind_sort=min(nb_max_cond,length(V_dist_t));
    %Path_to_sim(i)
    if ind_sort>0
        dist_max=V_dist_tempo_sort(ind_sort);
        %Simulated_pts
        V_ind_cond=find(V_dist_t<=dist_max);%
        %Simulated_pts(V_ind_cond)
        pts_condi=[];
        Z_condi=[];
        %conditioning to previous simulations
        for j=1:length(V_ind_cond)
            for k=1:length(Coord_Euler(:,1))
                temp=[Coord_Euler(k,1)-Simulated_pts(V_ind_cond(j))*step_t*m(10)*cos(m(11)*pi/180), Coord_Euler(k,2)-Simulated_pts(V_ind_cond(j))*step_t*m(10)*sin(m(11)*pi/180), Simulated_pts(V_ind_cond(j))*step_t];
                pts_condi=[pts_condi; temp];
                Z_condi=[Z_condi; M_res_sim(k,Simulated_pts(V_ind_cond(j)))];
            end
        end
        
        tmin=min(Simulated_pts(V_ind_cond))*step_t;
        tmax=max(Simulated_pts(V_ind_cond))*step_t;
        
        nb_harddata=0;
        for j=1:length(CondiRainTS_gaussian)

            my_idx=find((CondiRainTS_gaussian(j).t > (Path_to_sim(i)-10)*step_t) & (CondiRainTS_gaussian(j).t < (Path_to_sim(i)+10)*step_t));
            ind_min=min(my_idx);
            ind_max=max(my_idx);
            
            temp=[CondiRainTS_gaussian(j).X-CondiRainTS_gaussian(j).t(ind_min:ind_max)*m(10)*cos(m(11)*pi/180), CondiRainTS_gaussian(j).Y-CondiRainTS_gaussian(j).t(ind_min:ind_max)*m(10)*sin(m(11)*pi/180), CondiRainTS_gaussian(j).t(ind_min:ind_max)];
            pts_condi=[pts_condi; temp];
            Z_condi=[Z_condi; CondiRainTS_gaussian(j).RainRate(ind_min:ind_max)];
            nb_harddata=nb_harddata+length(CondiRainTS_gaussian(j).RainRate(ind_min:ind_max));
            
        end
    else
        
        pts_condi=[];
        Z_condi=[];
        nb_harddata=0;
        for j=1:length(CondiRainTS_gaussian)

            my_idx=find((CondiRainTS_gaussian(j).t > (Path_to_sim(i)-10)*step_t) & (CondiRainTS_gaussian(j).t < (Path_to_sim(i)+10)*step_t));
            ind_min=min(my_idx);
            ind_max=max(my_idx);
            temp=[CondiRainTS_gaussian(j).X-CondiRainTS_gaussian(j).t(ind_min:ind_max)*m(10)*cos(m(11)*pi/180), CondiRainTS_gaussian(j).Y-CondiRainTS_gaussian(j).t(ind_min:ind_max)*m(10)*sin(m(11)*pi/180), CondiRainTS_gaussian(j).t(ind_min:ind_max)];
            pts_condi=[pts_condi; temp];
            Z_condi=[Z_condi; CondiRainTS_gaussian(j).RainRate(ind_min:ind_max)];
            nb_harddata=nb_harddata+length(CondiRainTS_gaussian(j).RainRate(ind_min:ind_max));
            
        end
    end
    
    pts_target=[];

    for k=1:length(Coord_Euler(:,1))
        temp=[Coord_Euler(k,1)-Path_to_sim(i)*step_t*m(10)*cos(m(11)*pi/180), Coord_Euler(k,2)-Path_to_sim(i)*step_t*m(10)*sin(m(11)*pi/180), Path_to_sim(i)*step_t];
        pts_target=[pts_target;temp];
    end
    
    
    %covariance matrices
    if ~isempty(pts_condi)
        [Sigma_target]=M_cov_target(pts_target,m);
        [Sigma_condi]=M_cov_obs(pts_condi,nb_harddata,m);
        inv_Sigma_condi=inv(Sigma_condi);
        [Sigma_cross_cov]=M_cross_cov(pts_target,pts_condi,nb_harddata,m);
    else
        Z_condi=0;
        Sigma_cross_cov=0;
        inv_Sigma_condi=0;
        [Sigma_target]=M_cov_target(pts_target,m);
    end
    
    %simulation
    [Z_target]=core_simul(Z_condi, Sigma_cross_cov, Sigma_target, inv_Sigma_condi);

    for k=1:length(Coord_Euler(:,1))
        M_res_sim(k,Path_to_sim(i))=Z_target(k);
    end
    
    Simulated_pts=[Simulated_pts; Path_to_sim(i)];
end

%create final structure
SimulatedRainTS=struct();
for i=1:length(Coord_Euler(:,1))
    SimulatedRainTS(i).X=Coord_Euler(i,1);
    SimulatedRainTS(i).Y=Coord_Euler(i,2);
    SimulatedRainTS(i).t=zeros(length(Path_to_sim),1);
    SimulatedRainTS(i).RainRate=zeros(length(Path_to_sim),1);
    for j=1:nb_ep
        SimulatedRainTS(i).t(j)=j*step_t;
        SimulatedRainTS(i).RainRate(j)=M_res_sim(i,j);
    end
end

%------subfunctions-------
function[Z_target]=core_simul(Z_obs, Sigma_cross_cov, Sigma_target, inv_Sigma_obs)
    mu=Sigma_cross_cov'*inv_Sigma_obs*Z_obs;
    Sigma=Sigma_target-Sigma_cross_cov'*inv_Sigma_obs*Sigma_cross_cov;
    
    L=chol(Sigma,'lower');
    V_IID=randn(length(Sigma_target(:,1)),1);
    Z_target=mu+L*V_IID;
end

function[Sigma_target]=M_cov_target(pts_target,m)

    V_X=pts_target(:,1);
    V_Y=pts_target(:,2);
    V_t=pts_target(:,3);

    MVX1=repmat(V_X,1,length(V_X));
    MVX2=repmat(V_X',length(V_X),1);
    clear V_X
    MVY1=repmat(V_Y,1,length(V_Y));
    MVY2=repmat(V_Y',length(V_Y),1);
    clear V_Y
    
    M_ds=sqrt((MVX2-MVX1).^2+(MVY2-MVY1).^2);
    clear MVX1
    clear MVX2
    clear MVY1
    clear MVY2
    
    MVt1=repmat(V_t,1,length(V_t));
    MVt2=repmat(V_t',length(V_t),1);
    clear V_t

    M_dt=abs(MVt2-MVt1);
    clear MVt1
    clear MVt2

    to=1;
    c=m(1)^(-2*m(2));
    a=m(3)^(-2*m(4));
    Elem=a.*M_dt.^(2*m(4))+1;
    Sigma=1./(Elem.^to).*exp(-c.*(M_ds.^(2.*m(2)))./(Elem.^(m(5).*m(2))));
        
    Sigma_target=(1-m(6)^2)*Sigma; %no measurement noise between target points
end

function[Sigma_obs]=M_cov_obs(pts_obs,nb_harddata,m)

    V_X=pts_obs(:,1);
    V_Y=pts_obs(:,2);
    V_t=pts_obs(:,3);

    MVX1=repmat(V_X,1,length(V_X));
    MVX2=repmat(V_X',length(V_X),1);
    clear V_X
    MVY1=repmat(V_Y,1,length(V_Y));
    MVY2=repmat(V_Y',length(V_Y),1);
    clear V_Y
    
    M_ds=sqrt((MVX2-MVX1).^2+(MVY2-MVY1).^2);
    clear MVX1
    clear MVX2
    clear MVY1
    clear MVY2
    
    MVt1=repmat(V_t,1,length(V_t));
    MVt2=repmat(V_t',length(V_t),1);
    clear V_t

    M_dt=abs(MVt2-MVt1);
    clear MVt1
    clear MVt2

    to=1;
    c=m(1)^(-2*m(2));
    a=m(3)^(-2*m(4));
    Elem=a.*M_dt.^(2*m(4))+1;
    Sigma=1./(Elem.^to).*exp(-c.*(M_ds.^(2.*m(2)))./(Elem.^(m(5).*m(2))));
        
    %measurement noise only for hard data (and not internal conditioning)
    [si,~]=size(M_ds);
    Sigma_obs=(1-m(6)^2)*Sigma;
    M_noise=[[zeros(si-nb_harddata,si-nb_harddata) zeros(si-nb_harddata,nb_harddata)];[zeros(nb_harddata,si-nb_harddata) eye(nb_harddata)*(m(6)^2)]];
    Sigma_obs=Sigma_obs+M_noise;
end

function[Sigma_cross]=M_cross_cov(pts_target,pts_obs,nb_harddata,m)

    V_X_t=pts_target(:,1);
    V_Y_t=pts_target(:,2);
    V_X_o=pts_obs(:,1);
    V_Y_o=pts_obs(:,2);

    MVX1=repmat(V_X_t',length(V_X_o),1);
    MVY1=repmat(V_Y_t',length(V_Y_o),1);
    MVX2=repmat(V_X_o,1,length(V_X_t));
    MVY2=repmat(V_Y_o,1,length(V_Y_t));
    clear V_X_t;
    clear V_Y_t;
    clear V_X_o;
    clear V_Y_o;


    M_ds=sqrt((MVX2-MVX1).^2+(MVY2-MVY1).^2);
    clear MVX1;
    clear MVY1;
    clear MVX2;
    clear MVY2;

    V_time_t=pts_target(:,3);
    V_time_o=pts_obs(:,3);

    MVt1=repmat(V_time_t',length(V_time_o),1);
    MVt2=repmat(V_time_o,1,length(V_time_t));
    clear V_time_t;
    clear V_time_o;

    M_dt=abs(MVt2-MVt1);
    clear MVt1;
    clear MVt2;

    to=1;
    c=m(1)^(-2*m(2));
    a=m(3)^(-2*m(4));
    Elem=a.*M_dt.^(2*m(4))+1;
    Sigma_cross=1./(Elem.^to).*exp(-c.*(M_ds.^(2.*m(2)))./(Elem.^(m(5).*m(2))));
    Sigma_cross=(1-m(6)^2)*Sigma_cross;
    %measurement noise only for hard data (and not internal conditioning)
    [si,~]=size(M_ds);
    M_noise=(M_ds==0).*(M_dt==0).*[ones(si-nb_harddata,length(pts_target(:,1)))*(m(6)^2); zeros(nb_harddata,length(pts_target(:,1)))];
    Sigma_cross=Sigma_cross+M_noise;
end

end