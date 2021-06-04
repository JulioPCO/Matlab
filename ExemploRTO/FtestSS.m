function [ISS] = FtestSS(keepx)
%Steady State identification


%creating array for steadystate
steadyarray=[];
steadyarray=[steadyarray; keepx(1,:)]; %X1
steadyarray=[steadyarray; keepx(2,:)]; %X2
steadyarray=[steadyarray; keepx(3,:)]; %X3
steadyarray=[steadyarray; keepx(4,:)]; %X4
steadyarray=steadyarray';
size(steadyarray);

%Inicializing
j=0;
Xf=[];
vf=[];
df=[];
issteady=[];

%parameter from Cao-Rhineheart filter

lambda1=[1     0.650   0.650   0.450];
lambda2=[1      0.35    0.35  0.75];
lambda3=[1       0.3     0.3    0.45];

% Cao-Rhinheart Method based on R test
for i=1:size(steadyarray,2)
            k=2;
            
    for j=1:size(steadyarray,1)
        if j==1
           Xf(1,i)= mean(steadyarray(end-20:end,i)); %if the first u0 change is below the size array the program show bad results
           vf(1,i)= var(steadyarray(end-20:end,i))+0.0001;
           df(1,i)= 2*var(steadyarray(end-20:end,i))+0.0001;
        else
            Xf(j,i)= lambda1(k) .* steadyarray(j,i) + (1-lambda1(k)) .* Xf(j-1,i);
            vf(j,i)= lambda2(k) .* (steadyarray(j,i) - Xf(j-1,i)).^2 + (1-lambda2(k)) .* vf(j-1,i);
            df(j,i)= lambda3(k) .* (steadyarray(j,i) - steadyarray(j-1,i)).^2 + (1-lambda3(k)) .* df(j-1,i);
        end
        
        R(j,i)= (2 - lambda1(k)) * (vf(j,i)./df(j,i));
        Rc=[1.9 1.9 1.9 1.9];
        Rcrit=ones(size(R,1),1)*Rc;
        issteady(j,i) = R(j,i)<Rcrit(j,i);
        
    end
end

RP=(prod(R'))';
if size(steadyarray,2)>1
    ISS=(prod(issteady'))';
else
    ISS=issteady;
end

figure(7)
for i=1:size(steadyarray,2)   
    subplot(2,5,i)  
    plot(R(3:end,i));
    hold on
    plot(Rcrit(3:end,i),'-.');
end

figure(9)
plot(ISS,'+')

count=0;
for i=1:length(ISS)
   if ISS(i) == 0
       count = count +1;
   end
end
end