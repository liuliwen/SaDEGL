function bestval = SaDEGL(fname, fun, D, Lbound, Ubound, NP, itermax, runindex, fid);

f_name=sprintf('para/f_%2d_%2d.txt',fun,runindex);
fid_f=fopen(f_name,'w');

cr_name=sprintf('para/cr_%2d_%2d.txt',fun,runindex);
fid_cr=fopen(cr_name,'w');

pop_name=sprintf('para/pop_%2d_%2d.txt',fun,runindex);
fid_pop=fopen(pop_name,'w');

%%初始化种群
pop = Lbound + rand(NP, D) .* (Ubound-Lbound);
if fun == 7
    XRmin = 0*ones(NP,D); 
    XRmax = 600*ones(NP,D);
    pop = XRmin + rand(NP, D) .* (XRmax-XRmin);
    XRmin = [];
    XRmax = [];
end
if fun == 25
    XRmin = -2*ones(NP,D); 
    XRmax = 5*ones(NP,D); 
    pop = XRmin + rand(NP, D) .* (XRmax-XRmin);
    XRmin = [];
    XRmax = [];
end
vals=feval(fname,pop,fun);%计算适应度。
[bestVal,ibest]=min(vals);%当前最优值与其下标

f_rec=zeros(NP,1);
iter=1;
fprintf(fid,'iter:%d. 最优值为：%e\n',iter,bestVal);

k=0.3;
cr_i=ones(NP,D)*0.6;
% w_i=zeros(NP,1);
fm1=0.8;
fm2=0.2;
% f=1;
Wm=0.8;
f_rec=zeros(NP,1);
cr_all=zeros(NP,D);
upop=zeros(NP,D);
deltaF=zeros(NP,1);
p=0.1;
rot=0:1:NP-1;
while iter<itermax   
    ind = randperm(4);              % index pointer array
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a4  = a3(rt+1);
    rt = rem(rot+ind(2),NP);
    a5  = a4(rt+1);
    U=rand(NP,D);
    jrand=randi(D,NP,1);
    
    %变异操作    
    f1_i=normrnd(fm1,0.5,NP,1);
    f1_i=abs(f1_i);
    f2_i=normrnd(fm2,0.3,NP,1);
    f2_i=abs(f2_i);
    vpop=pop+repmat(f1_i,1,D).*(pop(a2,:)-pop(a3,:))+repmat(f2_i,1,D).*(repmat(pop(ibest,:),NP,1)-pop)+repmat(f2_i,1,D).*(pop(a4,:)-pop(a5,:));
    %越界处理
    if fun~=7
        index=find(vpop<Lbound);
        vpop(index)=(Lbound(index)+pop(index))/2;
        index=find(vpop>Ubound);
        vpop(index)=(Ubound(index)+pop(index))/2;
    end
    %交叉操作
    upop=pop;
    upop(U<cr_i)=vpop(U<cr_i);
    temp=(1:NP)';
    upop([temp,jrand])=vpop([temp,jrand]);
    cr_all(U<cr_i)=1;%记录交叉时选自变异向量的列
    cr_all([temp,jrand])=1;
    cr_suc=cr_all;%记录选择后第j列选自变异向量且成功进入下一代的个体
    tempvals = feval(fname, upop, fun);
    %选择
    for i=1:NP       
       if(tempvals(i)<=vals(i))
           if(tempvals(i)<vals(i))
               deltaF(i)=vals(i)-tempvals(i);
               f_rec(i)=1;               
           end
           pop(i,:)=upop(i,:);
           vals(i)=tempvals(i);                   
       else
           cr_suc(i,:)=0; 
       end
    end
    %记录参数信息
    fprintf(fid_f,'\niter:%4d,fm1:%f,fm2:%f.numbesrs success to next generate:%3d.\n\n',iter,fm1,fm2,sum(f_rec));
    fprintf(fid_cr,'\niter:%4d,.numbesrs success to next generate:%3d.\n\n',iter,sum(f_rec));
%   更新Fm
    if sum(f_rec)~=0       
        fm1=0.8*fm1+(0.2)*sum(f_rec.*f1_i.*deltaF.*f1_i./(f1_i+f2_i))/sum(f_rec.*deltaF.*f1_i./(f1_i+f2_i)); 
        fm2=0.8*fm2+(0.2)*sum(f_rec.*f2_i.*deltaF.*f2_i./(f1_i+f2_i))/sum(f_rec.*deltaF.*f2_i./(f1_i+f2_i)); 
        if fm1<0.2
           fm1=0.2;
       end
       if fm2>1
           fm2=1;
       end
       if fm2<0.2
           fm2=0.2;
       end
       if fm1>1
           fm1=1;
       end
       f_rec=zeros(NP,1);
    end
    %更新cr_i
        allSum=sum(cr_all,1);
        sucSum=sum(cr_suc,1);
        cr_all=zeros(NP,D);
        if sucSum~=0
%             cr_i=cr_i+repmat(0.1*(cr_i_sum-cr_i_sum_mean)/(cr_i_sum_max-cr_i_sum_min),NP,1);%这里写错了，应该是下面这种
%             cr_i=cr_i+repmat(0.05*(cr_i_sum*(2/(cr_i_sum_max-cr_i_sum_min))-(cr_i_sum_min+cr_i_sum_max)/(cr_i_sum_max-cr_i_sum_min)),NP,1);
            rate=sucSum./allSum;
            meanRate=mean(rate);
            cr_i=cr_i+repmat(rate-meanRate,NP,1);
            cr_i(cr_i>=0.9)=0.7;
            cr_i(cr_i<0.2)=0.4;
            cr_i(cr_i(:,allSum>5)<0.3)=0.5;
        end
    
    
     for i=1:NP
         fprintf(fid_f,'%f ',f1_i(i));
     end
     for i=1:D
         fprintf(fid_cr,'%f ',cr_i(1,i));
     end
    %记录参数信息
    fprintf(fid_f,'\n');
    fprintf(fid_cr,'\n');
    fprintf(fid_pop,'\n iter:%d\n',iter);
    for i=1:NP
        fprintf(fid_pop,'适应度为：%e---',vals(i));
        for j=1:D
            fprintf(fid_pop,'%f ',pop(i,j));
        end
        fprintf(fid_pop,'\n');
    end
   
    [bestVal,ibest]=min(vals);
    fprintf(fid,'fun:%4d, iter:%d, 最优值为：%e.最优值id：%d.\n',fun,iter,bestVal,ibest);
    iter=iter+1;
    fprintf('fun:%d,iter：%d,最优值为：%e.最优值id：%d.\n',fun,iter,bestVal,ibest);
end  
fclose(fid_f);
fclose(fid_cr);
bestval=bestVal;



