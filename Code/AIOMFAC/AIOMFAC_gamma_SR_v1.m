%--UNIFAC--
function [Ln_gamma_is_SR] = AIOMFAC_gamma_SR_v1(molecules,x_i,molecule_group_flag,molecule_group_stoich,T)
%A) ----------IMPROVING COMPUTATIONAL EFFICIENCY-------------
%**************DEFINING PERSISTENT VARIABLES***************************
persistent q_i
persistent r_i
persistent AFCData_Main_save
persistent s_i_n
persistent Q_k_i
persistent u_i_m_n
%***********************************************************************
molecule_group_flag=sparse(molecule_group_flag);
molecule_group_stoich=sparse(molecule_group_stoich);
%****************PRE ALLOCATING MATRICES******************
Ln_gamma_i_C=sparse(1,molecules);%----Combinatorial----
Ln_gamma_i_R=sparse(1,molecules);%----Residual---------
%----Final activity coefficients-------
Ln_gamma_is_SR=sparse(1,molecules);
%gamma_i=sparse(1,molecules);
%-----------REDUCING SIZE OF STOICHIOMETRIC MATRICES------------
test=unique(molecule_group_flag);
test2=nonzeros(test);
[r_unique,c_unique]=size(test2);
max_group_num=r_unique;
%now we want to re-arrange the
%molecule_group_flag file into a matrix with same rows as molecules and
%same columns as the number of functional groups within all molecules.
molecule_group_stoich_new=sparse(molecules,max_group_num); %the new matrix to be populated
group_flag_array_new(1:max_group_num)=1;
size_matrix=size(molecule_group_stoich);
step=1;
for i=1:max_group_num
    group_flag_array_new(step)=test2(i,1);
    clear ind
    %[r,c]=find(molecule_group_flag==test2(i,1));
    ind=find(molecule_group_flag==test2(i,1));
    [rtest,ctest]=ind2sub(size_matrix,ind);
    if (isempty(ind)==0)
    %molecule_group_stoich_new([r],i)=molecule_group_stoich([r],[c]);
    molecule_group_stoich_new(rtest,i)=molecule_group_stoich(ind);
    end
    step=step+1;
end
clear test test2 c_unique r_unique molecule_group_stoich size_matrix
%-----------------------------------------------------------------


%*****************************************************************
%---------------PERSISTENT CALCULATIONS-----------------------------
if isempty(q_i)
    q_i=sparse(1,molecules);
    w_i=sparse(1,molecules);
    ss_n=sparse(1,max_group_num);
    AFCData_ri_save = AFC_datari_v1();
    AFCData_qi_save = AFC_dataqi_v1();
    AFCData_Main_save = AFC_datamain_v1();
    AFC_SR_save = AFC_SR_interact_params();
    %Q_k=sparse(molecules,max_group_num);
    %R_k=sparse(molecules,max_group_num);
    %vk_i=sparse(molecules,max_group_num);
    Q_k_i=sparse(molecules,max_group_num);
    R_k_i=sparse(molecules,max_group_num);
    q_i=sparse(molecules);
    r_i=sparse(molecules);
    %-----------COMBINATORIAL CONTRIBUTION-------------------------
    %h = waitbar(0,'Creating persistent molecular matrix in unifac...');
    for mol=1:molecules
    %waitbar(mol/molecules)
    %order the functional groups in some way
    %[r,c]=size(molecule_group_flag(mol,:));
    for i=1:max_group_num
        %here look at each group identifier then pull out
        %the group number and stoichiometry seperately
        %dosnt matter if a group is repeated as the interactions
        %between each are set to zero anyway (which you can do manually)
        if (molecule_group_stoich_new(mol,i)~=0)
        %%Q_k(mol,i)=AFCData_qi_save(group_flag_array_new(i));
        %%R_k(mol,i)=AFCData_ri_save(group_flag_array_new(i));
        %%vk_i(mol,i)=molecule_group_stoich_new(mol,i);
        %mol_func_groups_rec(i)=molecule_group_flag(mol,i)     
        Q_k_i(mol,i)=AFCData_qi_save(group_flag_array_new(i))*molecule_group_stoich_new(mol,i);
        R_k_i(mol,i)=AFCData_ri_save(group_flag_array_new(i))*molecule_group_stoich_new(mol,i);
        end
    end
    q_i(mol)=sum(Q_k_i(mol,:));
    r_i(mol)=sum(R_k_i(mol,:));
    end
    %close(h)
    clear R_k_i R_k AFCData_ri_save AFCData_qi_save
    %--------------RESIDUAL CONTRIBUTION---------------------------
    s_i_n=sparse(molecules,max_group_num);
    t_i_m_n=zeros(molecules,max_group_num,max_group_num); %cannot work with 3D sparse matrices
    u_i_m_n=zeros(molecules,max_group_num,max_group_num); %cannot work with 3D sparse matrices
    %h = waitbar(0,'Creating first persistent residual matrix in unifac...');
    for mol=1:molecules
    %waitbar(mol/molecules)
    for i=1:max_group_num
        for j=1:max_group_num
        a_m_n_1=AFC_SR_save(AFCData_Main_save(group_flag_array_new(i)),AFCData_Main_save(group_flag_array_new(j)));
        a_m_n_2=AFC_SR_save(AFCData_Main_save(group_flag_array_new(j)),AFCData_Main_save(group_flag_array_new(i)));
        %-------------------------
        alpha_m_n_1=a_m_n_1/T;
        alpha_m_n_2=a_m_n_2/T;
        tau_m_n_1=exp(-1*alpha_m_n_1);
        tau_m_n_2=exp(-1*alpha_m_n_2);
        test=Q_k_i(mol,i);
        t_i_m_n(mol,i,j)=test.*tau_m_n_1;
        u_i_m_n(mol,i,j)=test.*tau_m_n_2;
        %------------------------
        end
    end
    end
    clear AFCData_Main_save AFC_SR_save
    %close(h)
    %h = waitbar(0,'Creating second persistent residual matrix in unifac...');
    for mol=1:molecules
   % waitbar(mol/molecules)
        for m=1:max_group_num
        s_i_n(mol,m)=sum(t_i_m_n(mol,:,m));   %sum over all functional groups in molecule i            
        end
    end
    %close(h)
    clear t_i_m_n alpha_m_n_1 alpha_m_n_2 tau_m_n_1 tau_m_n_2
end

%----------------------------------------------------------------
%-----------FULL CALCULATION------------------
%--------------combinatorial------------------
Q=sum(x_i.*q_i);
%for mol=1:molecules
    w_i=q_i/Q;
%end
omega_i=log(w_i);

R=sum(x_i.*r_i);
%for mol=1:molecules
    row_i=r_i/R;
%end
P_i=log(row_i);
%row_i=cell(molecules);
%--combine the two--
%%for mol=1:molecules
    delta_i=row_i*(5*Q-1);
    cross_i=5*q_i;
%%end
%--calculate combinatorial part--
%for mol=1:molecules
temp_const1(1,1:molecules)=1;
temp_const2(1,1:molecules)=1;
Ln_gamma_i_C(1,1:molecules)=temp_const1+delta_i+P_i+cross_i.*(omega_i-P_i-temp_const2);
%end
clear Q R temp_const1 temp_const2 delta_i P_i cross_i 
%----------------residual----------------------
%a_m_n=sparse(max_group_num,max_group_num);
temp=sparse(molecules);
uu_m_n=sparse(max_group_num,max_group_num);
for i=1:max_group_num
    for j=1:max_group_num
        temp=u_i_m_n(1:molecules,i,j)';
        uu_m_n(i,j)=sum(x_i.*temp);  %sum over all molecules         
    end
end
clear temp

temp=sparse(molecules);
ss_n=sparse(max_group_num);
xx_weird_i_n=sparse(molecules,max_group_num);
for m=1:max_group_num
     temp=s_i_n(:,m)';
     ss_n(m)=sum(x_i.*temp);  %sum over all molecules%
     xx_weird_i_n(1:molecules,m)=log(s_i_n(1:molecules,m)/ss_n(m));   
end
clear temp
%--calculate the residual part--

%h = waitbar(0,'Calculating temporary residual contribution in unifac...');
for mol=1:molecules
   % waitbar(mol/molecules)
    summation1=0;
    for n=1:max_group_num
        summation2=0;
        brac1=0;
        brac2=0;
        for m=1:max_group_num
            brac1=brac1+(u_i_m_n(mol,m,n))/s_i_n(mol,m);
            brac2=brac2+(uu_m_n(m,n))/ss_n(m);             
        end
        summation2=xx_weird_i_n(mol,n)-omega_i(mol)+brac1-brac2;
        summation1=summation1+Q_k_i(mol,n)*summation2;
    end
Ln_gamma_i_R(mol)=summation1;
end
%close(h)
clear uu_m_n ss_n xx_weird_i_n omega_i brac1 brac2 summation1 summation2

%-------------------------------
%--final activity coefficient---
Ln_gamma_is_SR=Ln_gamma_i_R+Ln_gamma_i_C;
%gamma_i=exp(Ln_gamma_i);