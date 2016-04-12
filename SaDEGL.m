function bestval = SaDEGL(fname, fun, D, Lbound, Ubound, NP, itermax, runindex, fid);

w_name=sprintf('para/w_%2d_%2d.txt',fun,runindex);
fid_w=fopen(w_name,'w');

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
cr_i=ones(NP,D)*0.8;
w_i=zeros(NP,1);
fm1=0.5;
fm2=1;
% f=1;
Wm=0.7;
w_rec=zeros(NP,1);
cr_rec=zeros(NP,D);
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
%     pbest_index=randi(floor(D*p),NP,1);
    %变异操作
%     [sortedVal,sortedIndex]=sort(vals,'descend');
    
    w_i=normrnd(Wm,0.2,NP*5,1);
    w_i=abs(w_i);        
    tempindex=find(w_i<1);
    w_i=w_i(tempindex(1:NP));
    
    w=1-sigmf(10*iter/itermax,[10,0.5])*0.2;
%     f1=1-sigmf(10*iter/itermax,[10,0.5])*0.5;
%     f2=1+sigmf(10*iter/itermax,[10,0.5])*0.5;
%     w=1-sigmf(10*iter/itermax,[10,0.5])*0.8;
    f1_i=normrnd(fm1,0.3,NP,1);
    f1_i=abs(f1_i);
    f2_i=normrnd(fm2,0.3,NP,1);
    f2_i=abs(f2_i);
    vpop1=pop+repmat(f1_i,1,D).*(pop(a2,:)-pop(a3,:));
    vpop2=pop+repmat(f2_i,1,D).*(repmat(pop(ibest,:),NP,1)-pop)+repmat(f2_i,1,D).*(pop(a4,:)-pop(a5,:));
%      f=1;
%     vpop1=pop+f1*(pop(a2,:)-pop(a3,:))+f1*(pop(a4,:)-pop(a5,:));
%     vpop2=pop+f2*(repmat(pop(ibest,:),NP,1)-pop)+f2*(repmat(pop(ibest,:),NP,1)-pop(a1,:));
%     vpop=repmat(w_i,1,D).*vpop1+repmat((1-w_i),1,D).*vpop2;
    vpop=w*vpop1+(1-w)*vpop2;
    %越界处理
    if fun~=7
%         tempLbound=repmat(Lbound,[NP,1]);
        index=find(vpop<Lbound);
        vpop(index)=(Lbound(index)+pop(index))/2;
%         tempUbound=repmat(Ubound,[NP,1]);
        index=find(vpop>Ubound);
        vpop(index)=(Ubound(index)+pop(index))/2;
    end
    %交叉操作
    upop=pop;
    upop(U<cr_i)=vpop(U<cr_i);
    temp=(1:NP)';
    upop([temp,jrand])=vpop([temp,jrand]);
    cr_rec(U<cr_i)=1;
    cr_rec([temp,jrand])=1;
    
    tempvals = feval(fname, upop, fun);
    for i=1:NP
        
       if(tempvals(i)<=vals(i))
           if(tempvals(i)<vals(i))
               deltaF(i)=vals(i)-tempvals(i);
               w_rec(i)=1;
           end
           pop(i,:)=upop(i,:);
           vals(i)=tempvals(i);                   
       else
           cr_rec(i,:)=0; 
       end
    end
    
%     fprintf(fid_w,'\niter:%4d,Wm:%f,f:%f.numbesrs success to next generate:%3d.\n\n',iter,Wm,f,sum(w_rec));
    fprintf(fid_w,'\niter:%4d,Wm:%f.numbesrs success to next generate:%3d.\n\n',iter,Wm,sum(w_rec));
    fprintf(fid_f,'\niter:%4d,fm1:%f,fm2:%f.numbesrs success to next generate:%3d.\n\n',iter,fm1,fm2,sum(w_rec));
    fprintf(fid_cr,'\niter:%4d,.numbesrs success to next generate:%3d.\n\n',iter,sum(w_rec));
%     更新Wm
    if sum(w_rec)~=0
%      Wm=sum(w_rec.*w_i.*deltaF)/sum(w_rec.*deltaF);
       Wm=0.8*Wm+(0.2)*(sum(w_rec.*w_i)/sum(w_rec)); 
       fm1=0.8*fm1+(0.2)*(sum(w_rec.*f1_i)/sum(w_rec)); 
       fm2=0.8*fm2+(0.2)*(sum(w_rec.*f2_i)/sum(w_rec)); 
       if fm1<0.3
           fm1=0.3;
       end
       if fm2>1.5
           fm2=1.5;
       end
       if fm2<0.3
           fm2=0.3;
       end
       if fm1>1.5
           fm1=1.5;
       end
       if Wm<0.3
           Wm=0.3;
       end
       w_rec=zeros(NP,1);
    end
    %更新cr_i
    if rem(iter,5)==0 & iter~=0
        cr_i_sum=sum(cr_rec,1);
        cr_i_sum_min=min(cr_i_sum);
        cr_i_sum_max=max(cr_i_sum);
        cr_i_sum_mean=mean(cr_i_sum);
        if cr_i_sum_max-cr_i_sum_min~=0
            cr_i=cr_i+repmat(0.1*(cr_i_sum-cr_i_sum_mean)/(cr_i_sum_max-cr_i_sum_min),NP,1);
            cr_i(cr_i>0.9)=0.7;
            cr_i(cr_i<0.2)=0.4;
            cr_i(cr_i(:,cr_i_sum>5)<0.3)=0.5;
        end
    end
    
    
     for i=1:NP
         fprintf(fid_w,'%f ',w_i(i));
         fprintf(fid_f,'%f ',f1_i(i));
     end
     for i=1:D
         fprintf(fid_cr,'%f ',cr_i(1,i));
     end
     fprintf(fid_w,'\n');
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
fclose(fid_w);
fclose(fid_f);
fclose(fid_cr);
bestval=bestVal;



