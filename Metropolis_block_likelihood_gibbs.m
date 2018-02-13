function [V_tot, V_sample]=Metropolis_block_likelihood_gibbs(m0,V_step_walk,M_prior,MeasuredRainTS,size_block,nb_ep_warmup,nb_sample,step, nb_iter_gibbs)


%-----------------------------Prepare data--------------------------
display('Prepare data');

Coord=[];
Rm=[];

nb_ep_mes=length(MeasuredRainTS(1).RainRate);
nb_blocks=floor(nb_ep_mes/size_block);

size_block_coord=length(MeasuredRainTS)*size_block;

ind_Z=0;

for b=1:nb_blocks
    for i=1:size_block
        for j=1:length(MeasuredRainTS)
            ind_Z=ind_Z+1;
            Coord=[Coord;[MeasuredRainTS(j).X, MeasuredRainTS(j).Y, MeasuredRainTS(j).t((b-1)*size_block+i)]];
            Rm=[Rm;MeasuredRainTS(j).RainRate((b-1)*size_block+i)];
        end
    end
end

Z_k=Rm;
for i=1:length(Z_k)
    if Rm(i)==0
        Z_k(i)=abs(-2-m0(7))*rand-2;
    end
end

%---------------------------------Initialization-----------------------------------
display('Initialization');
[Sigma_k,inv_Sigma_k,Z_k]=apply_model(m0, Coord, size_block_coord, Rm, Z_k);
[Z_k]=gibbs_update_Z0(m0, Z_k, Rm, inv_Sigma_k, size_block_coord,500);


%--------------------Metropolis inversion (Mosegaard & Tarantola, 1995)---------------
display('Sampling');
V_tot=zeros(nb_ep_warmup+nb_sample,length(m0)+2+1);
V_sample=zeros(nb_sample,length(m0)+length(Z_k));
nb_accept=0;
alpha_covariance=0;
alpha_anamorphosis=0;

q0=ones(length(m0),1);

for i=1: nb_ep_warmup + nb_sample*step
    
    L0=likelihood_tot(Sigma_k, inv_Sigma_k, Z_k, Rm, m0, size_block_coord);
    
    [m1,q1]=random_walk(m0,q0,V_step_walk,M_prior);

    [Sigma_p,inv_Sigma_p,Z_p]=apply_model(m1, Coord, size_block_coord, Rm, Z_k);
    L1=likelihood_tot(Sigma_p, inv_Sigma_p, Z_p, Rm, m1, size_block_coord);
    
    alpha=exp(L1-L0)*prod(q1)/prod(q0);
    
    if alpha > 1 || rand > 1-alpha
                
        m0=m1;
        q0=q1;
        Sigma_k=Sigma_p;
        inv_Sigma_k=inv_Sigma_p;
        Z_k=Z_p;
                
        [Z_k]=gibbs_update_Z0(m0, Z_k, Rm ,inv_Sigma_k, size_block_coord ,nb_iter_gibbs);
                
        nb_accept=nb_accept+1;
       
    end
            
    V_tot(i,:)=[m0, alpha_anamorphosis, alpha_covariance, L1];

    if mod(i-nb_ep_warmup,step)==0 && i > nb_ep_warmup
        V_sample(floor((i-nb_ep_warmup)/step),:)=[m0, Z_k'];
    end
    
    display(strcat('sampling:',num2str(i),'/',num2str(+nb_ep_warmup + nb_sample*step),' Acceptance rate:',num2str(nb_accept/i)));
    
end
display(strcat('Acceptance rate:',num2str(nb_accept/(nb_ep_warmup + nb_sample*step))));

%-------------------------------Likelihood functions---------------------------------
function [sum_log_likelihood]=likelihood_tot(Sigma, inv_Sigma, Z, Rm, m, size_block_coord)

    nb=length(Z)/size_block_coord;
    sum_log_likelihood_structure=0;
    sum_log_likelihood_anamorphose=0;
    
    for i1=1:nb
        Rmb=Rm((i1-1)*size_block_coord+1:i1*size_block_coord);
        Zb=Z((i1-1)*size_block_coord+1:i1*size_block_coord);
        
        NI=sum(Rmb>0);
        log_likelihood_structure=-0.5*logdet(Sigma)-0.5*Zb'*inv_Sigma*Zb-NI*log(2*pi);
        sum_log_likelihood_structure=sum_log_likelihood_structure+log_likelihood_structure;
        
        log_likelihood_anamorphose=0;
        nb_pos_data=0;
        for i2=1:size_block_coord
            if Rmb(i2)>0.5 
                L=m(8)*m(9)*Rmb(i2)^(m(9)-1);
                log_likelihood_anamorphose=log_likelihood_anamorphose+log(abs(L));
                nb_pos_data=nb_pos_data+1;
            end
        end
       
        sum_log_likelihood_anamorphose=sum_log_likelihood_anamorphose+log_likelihood_anamorphose;
    end

    sum_log_likelihood=sum_log_likelihood_structure+sum_log_likelihood_anamorphose;
    if isinf(sum_log_likelihood)||isnan(sum_log_likelihood)
        display('Undefined likelihood!!')
    end
    
end



%----------------------------------Gibbs sampler---------------------------------
function [Z]=gibbs_update_Z0(m, Z, Rm, inv_Sigma,size_block_coord, nb_iter_gibbs)

    
    for i_gibbs=1:nb_iter_gibbs
        
        %initial values
        Z0=Z(1:size_block_coord);
        Rmb=Rm(1:size_block_coord);
        for i2=1:size_block_coord
            if Rmb(i2)==0
                [ mu_gibbs, sigma_gibbs ] = conditional_normal_sim(inv_Sigma,i2,Z0);
                sim_value=mu_gibbs+randn*sigma_gibbs;
                it=1;   
                while sim_value>m(7)
                    sim_value=mu_gibbs+randn*sigma_gibbs;
                    it=it+1;
                    if it>1000
                        sim_value=Z0(i2);
                        break
                    end
                end
                Z0(i2)=sim_value;
            end
        end
        Z(1:size_block_coord)=Z0;
        
        %final values
        Zf=Z(end-size_block_coord+1:end);
        Rmb=Rm(end-size_block_coord+1:end);
        for i2=1:size_block_coord
            if Rmb(i2)==0
                [ mu_gibbs, sigma_gibbs ] = conditional_normal_sim(inv_Sigma,i2,Zf);
                sim_value=mu_gibbs+randn*sigma_gibbs;
                it=1;   
                while sim_value>m(7)
                    sim_value=mu_gibbs+randn*sigma_gibbs;
                    it=it+1;
                    if it>1000
                        sim_value=Zf(i2);
                        break
                    end
                end
                Zf(i2)=sim_value;
            end
        end
        Z(end-size_block_coord+1:end)=Zf; 
        
        
        %all non extramal values
        for ic=size_block_coord+1:length(Z)-size_block_coord-1
            if Rm(ic)==0
                ind_ini_Zc=ic-floor(size_block_coord/2);
                Zc=Z(ind_ini_Zc:ind_ini_Zc+size_block_coord-1);
                [ mu_gibbs, sigma_gibbs ] = conditional_normal_sim(inv_Sigma,floor(size_block_coord/2)+1,Zc);
                sim_value=mu_gibbs+randn*sigma_gibbs;
                it=1;   
                while sim_value>m(7) || sim_value<-3
                    sim_value=mu_gibbs+randn*sigma_gibbs;
                    it=it+1;
                    if it>1000
                        sim_value=Z(ic);
                    break
                    end
                end
            Z(ic)=sim_value;
            end 
        end
        
    end
    
end

%------------------------------Forward model => for likelihood and gibbs computation------------------

function[Sigma,inv_Sigma,Z]=apply_model(m, Coord, size_block_coord, Rm, Z)
    %marginal distribution
    for jj=1:length(Rm)
        if Rm(jj)>0
            Z(jj)=m(8)*( Rm(jj)^m(9) )+m(7);
        end
    end
    %covariance structure
    my_coord=Coord(1:size_block_coord,:);
    if m(10)==0 && m(11)==0
        Coord_Lagrang=[my_coord(:,1), my_coord(:,2), my_coord(:,3)]; 
    else
        Coord_Lagrang=[my_coord(:,1)-my_coord(:,3)*m(10)*cos(m(11)*pi/180), my_coord(:,2)-my_coord(:,3)*m(10)*sin(m(11)*pi/180), my_coord(:,3)]; 
    end
    
    
    V_X=Coord_Lagrang(:,1);
    V_Y=Coord_Lagrang(:,2);
    V_t=Coord_Lagrang(:,3);

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
        
    Sigma=(1-m(6)^2)*Sigma+eye(size_block_coord)*(m(6)^2);

    inv_Sigma=inv(Sigma);
        
end



%--------------------------------------random walk--------------------------------

function [m1,q1]=random_walk(m0,q0,V_step_walk,M_prior)
    m1=m0; 
    q1=q0;
    
    %spatial range
    step_walk=V_step_walk(1);
    B1=max(M_prior(1,1),m0(1)-step_walk);
    B2=min(M_prior(1,2),m0(1)+step_walk);
    m1(1)=rand*(B2-B1)+B1;
    q1(1)=1/(B2-B1);
    
    %space regularity parameter (smoothness)
    step_walk=V_step_walk(2);
    B1=max(M_prior(2,1),m0(2)-step_walk);
    B2=min(M_prior(2,2),m0(2)+step_walk);
    m1(2)=rand*(B2-B1)+B1;
    q1(2)=1/(B2-B1);
    
    %temporal range
    step_walk=V_step_walk(3);
    B1=max(M_prior(3,1),m0(3)-step_walk);
    B2=min(M_prior(3,2),m0(3)+step_walk);
    m1(3)=rand*(B2-B1)+B1;
    q1(3)=1/(B2-B1);
    
    %time regularity parameter (smoothness)
    step_walk=V_step_walk(4);
    B1=max(M_prior(4,1),m0(4)-step_walk);
    B2=min(M_prior(4,2),m0(4)+step_walk);
    m1(4)=rand*(B2-B1)+B1;
    q1(4)=1/(B2-B1);
    
    %space-time interaction parameter
    step_walk=V_step_walk(5);
    B1=max(M_prior(5,1),m0(5)-step_walk);
    B2=min(M_prior(5,2),m0(5)+step_walk);
    m1(5)=rand*(B2-B1)+B1;
    q1(5)=1/(B2-B1);
    
    %noise
    step_walk=V_step_walk(6);
    B1=max(M_prior(6,1),m0(6)-step_walk);
    B2=min(M_prior(6,2),m0(6)+step_walk);
    m1(6)=rand*(B2-B1)+B1;
    q1(6)=1/(B2-B1);
    
    %advection velocity
    step_walk=V_step_walk(10);
    B1=max(M_prior(10,1),m0(10)-step_walk);
    B2=min(M_prior(10,2),m0(10)+step_walk);
    m1(10)=rand*(B2-B1)+B1;
    q1(10)=1/(B2-B1);
    
    %advection direction
    step_walk=V_step_walk(11);
    B1=max(M_prior(11,1),m0(11)-step_walk);
    B2=min(M_prior(11,2),m0(11)+step_walk);
    m1(11)=rand*(B2-B1)+B1;
    q1(11)=1/(B2-B1);
    
    %a1
    step_walk=V_step_walk(8);
    B1=max(M_prior(8,1),m0(8)-step_walk);
    B2=min(M_prior(8,2),m0(8)+step_walk);
    m1(8)=rand*(B2-B1)+B1;
    q1(8)=1/(B2-B1);
    
    %a2
    step_walk=V_step_walk(9);
    B1=max(M_prior(9,1),m0(9)-step_walk);
    B2=min(M_prior(9,2),m0(9)+step_walk);
    m1(9)=rand*(B2-B1)+B1;
    q1(9)=1/(B2-B1);
    
end


end