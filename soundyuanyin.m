function zzg
[aa,bb]=audioread('D:matlab2020a\bin\a长元音.m4a'); 
%读取文件
%[aa,bb]=audioread('D:\matlab2020a\bin\eyin2.M4A');
aa=aa(:,1); 
ma=max(aa);
aa=aa/ma;
sound(aa,bb)
cc=aa(abs(aa)>0.0017);
cc=cc(22000:32000); %截取部分音频  %40000点开始较短 a取22000:32000 u取30000
%到40000  o取18000：28000 e取10000：20000 i取30000：40000
scc=size(cc)        %取10000以内的点
sound(cc,bb)
figure(9)
 plot(aa)
hold on
% plot(cc)
hold off
 sFF=size(cc);
 xx=0:1:sFF(1)-1;
 ym=xx;
%  %----------------------------------------------------分析小波空间
%  % 初始化变量
% square_coe3 = 0;
% square_coe4 = 0;
% square_coe5 = 0;
% square_coe6 = 0;
% square_coe7 = 0;
% square_coe8 = 0;
% for qq=3:8
% [TT,coe]=signaldecompositionsixspline(cc,qq,xx);
% kk=TT/2^qq;
% figure(qq)
% stem(kk,coe)
% %---------------创建树状图
% if(qq==3)
%     coe3=coe;
%     square_coe3 = sum(coe3.^2);
% end
% if(qq==4)
%     coe4=coe;
%     square_coe4 = sum(coe4.^2);
% end
% if(qq==5)
%     coe5=coe;
%     square_coe5 = sum(coe5.^2);
% end
% if(qq==6)
%     coe6=coe;
%     square_coe6 = sum(coe6.^2);
% end
% if(qq==7)
%     coe7=coe;
%     square_coe7 = sum(coe7.^2);
% end
% if(qq==8)
%     coe8=coe;
%     square_coe8 = sum(coe8.^2);
% end
% % 创建一些示例数据
% x = 3:8;
% y = [square_coe3,square_coe7, square_coe5, square_coe6, square_coe4,square_coe8];
% 
% % 使用 bar 函数创建柱状图
% bar(x, y);
% % 添加标题和标签
% xlabel('Wavelet Space');
% ylabel('Wavelet Energy');
% %---------------------------------------------------柱状图结束
% 
% if(qq>4)
%  ym(qq-4,:)=recovwav(xx,kk,coe,-qq);
% end
% if(qq==4)
%     ym1=recovwav(xx,kk,coe,-qq);
% end;
% qq
% end
% %------------------------------部分结束
% %plot(xx,ym)

[ym1,yout1]=SoFoProc(cc,7,xx,1.2,50);
[ym2,yout2]=SoFoProc(cc,6,xx,1.2,50);
[ym3,yout3]=SoFoProc(cc,5,xx,1.2,50);
ym=ym1+ym2+ym3; 
figure(10)
plot(xx,yout3) 
title('Original Signal') % 设置标题为 Original Signal

figure(11)
plot(xx,ym) 
title('Synthetic Signal') % 设置标题为 Original Signal
% figure(10)
% plot(xx,yout3,'r-')
% hold on 
% plot(xx,ym,'b-')
% hold off
yout=yout1+yout2+yout3;
    sound(yout,bb) 
pause(2) 
sound(ym,bb)
pause(2)
sound(yout,bb) 
pause(2) 
sound(ym,bb)
pause(2) 
sound(yout,bb)
pause(2)
sound(ym,bb)
pause(2)
sound(yout,bb)
pause(2) 
sound(ym,bb)
pause(2)
sound(yout,bb)
pause(2)
sound(ym,bb)
%------------------------------------------------------------------------
% qq=5;
% sFF=size(cc);
% xx=0:1:sFF(1)-1;
% [TT,coe1]=signaldecompositionsixspline(cc,qq,xx);
% kk=TT/2^qq;
% %---------------------------------------------------------
% sdd=size(dd);
% xd=0:1:sdd(1)-1;
% [TT2,coe2]=signaldecompositionsixspline(dd,qq,xd);
% kk2=TT2/2^qq;
%-------------------------------------

return
function [ym,yout]=SoFoProc(yy,qq,xx,thrHold,time)
[TT,coe]=signaldecompositionsixspline(yy,qq,xx);
kk=TT/2^qq;
%-----------------------------
xa=0:2^(qq-1):xx(end);
xa=xa(2:2:end);
yx=recovwav(xa,kk,coe,-qq);
ym=recovwav(xx,kk,coe,-qq);
%--------------------------------------------
ww=-pi:0.001:pi;
yx=yx';
re1=fourierzzg(ww,yx,1);
energy = abs(re1).^2;
figure(qq)
% 计算频谱的能量
plot(ww,energy)
%------------------------------------------
syx=size(yx);
NN=syx(2);
nn=1:NN;
% % -------------------------------------------
re1b=abs(re1);
mnre1=mean(re1b);
re1ref=re1.*(re1b>thrHold*mnre1);
rerealA=real(re1ref);
reimA=imag(re1ref);
[rexx,rereal]=diffseq(ww,rerealA);
figure(qq+10)
% stem(rexx,rereal,'b') %这一部分挪到后面去  %要插值则挪到后面
[rexx,reim]=diffseq(ww,reimA);   %取极值在这里插入值  
% %%%%____________自己插值  
% %%第一部分为qq=5
% if qq == 5 
%    %reim = zeros(size(reim)); %%%%虚部全为0
%    %rereal = zeros(size(rereal)); %%%%实部全为0
%    %%设置实部
%    rereal(:) = 0;
%    putf = [-3, -2.2, -0.34];   %设置放置频谱数值
%    for a = 1:3
%    p1 = (putf(a)-(-pi))/0.001;   %计算相应位置
%    dp1 = floor(p1); %向下取整
%    if a == 1 
%    for b = 1:5                            %设置加几组正负相见的数值
%    randomInteger = randi([2, 4]);  %加一组范围之间的随机数
%    rereal(dp1) = 1*randomInteger;  %乘以随机数
% %    figure(b+15)
% %    stem(rexx,rereal,'b')%debug
%    ff=dp1+5;
%    % randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
%    rereal(ff)=-rereal(dp1);   %不加随机数
%    % rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
%    dp1 = dp1+7;
% %    figure(b+30)
% %    stem(rexx,rereal,'b') %debug
%    end
%    elseif a == 2
%     for b = 1:2                             %设置加几组正负相见的数值
%    randomInteger = randi([1, 1]);  %加一组范围内之间的随机数
%    rereal(dp1) = 1*randomInteger;   %乘以随机数
%    % randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
%    % rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
%    rereal(dp1+10)=-rereal(dp1);   %不加随机数
%    dp1 = dp1+b*10;
%     end
%    else a == 3 
%     for b = 1:6                             %设置加几组正负相见的数值
% randomInteger = randi([2, 10]);  %加一组0-5之间的随机数
% rereal(dp1) = 1*randomInteger/10;   %乘以随机数
% %可能是这个原因 去掉加的随机数
% %randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
% %rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
% rereal(dp1+10)=-rereal(dp1);
% dp1 = dp1+48;
%     end
% end
%    end
% subArray = rereal(1:3142);
% reversedArray = flip(subArray);
% rereal = [subArray,reversedArray];
% reim = circshift(rereal, 10)*0.8;
% elseif qq == 6
% %reim = zeros(size(reim)); %%%%虚部全为0
%    %rereal = zeros(size(rereal)); %%%%实部全为0
%    %%设置实部
%    rereal(:) = 0;
%    putf = [-2.8];   %设置放置频谱数值
%    for a = 1:1
%    p1 = (putf(a)-(-pi))/0.001;   %计算相应位置
%    dp1 = floor(p1); %向下取整
%    if a == 1 
%    for b = 1:11                            %设置加几组正负相见的数值
%    randomInteger = randi([2, 4]);  %加一组范围之间的随机数
%    rereal(dp1) = 1*randomInteger;   %乘以随机数
%    % randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
%    rereal(dp1+5)=-rereal(dp1);   %不加随机数
%    % rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
%    dp1 = dp1+45;
%    end
%    elseif a == 2
%     for b = 1:2                             %设置加几组正负相见的数值
%    randomInteger = randi([1, 2]);  %加一组范围内之间的随机数
%    rereal(dp1) = 1*randomInteger;   %乘以随机数
%    % randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
%    % rereal(dp1+20)=-rereal(dp1)+ randomNumber;     %设置相同相反的数值在附近
%    rereal(dp1+10)=-rereal(dp1);   %不加随机数
%    dp1 = dp1+b*20;
%     end
%    else a == 3 
%     for b = 1:1                             %设置加几组正负相见的数值
% randomInteger = randi([0, 5]);  %加一组0-5之间的随机数
% rereal(dp1) = 1*randomInteger;   %乘以随机数
% %可能是这个原因 去掉加的随机数
% randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
% rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
% rereal(dp1+20)=-rereal(dp1);
% dp1 = dp1+b*50;
%     end
% end
%    end
% subArray = rereal(1:3142);
% reversedArray = flip(subArray);
% rereal = [subArray,reversedArray];
% reim = circshift(rereal, 10)*0.8;
% else qq == 7
%     %reim = zeros(size(reim)); %%%%虚部全为0
%    %rereal = zeros(size(rereal)); %%%%实部全为0
%    %%设置实部
%    rereal(:) = 0;
%    putf = [-2.7];   %设置放置频谱数值
%    for a = 1:1
%    p1 = (putf(a)-(-pi))/0.001;   %计算相应位置
%    dp1 = floor(p1); %向下取整
%    if a == 1 
%    for b = 1:5                            %设置加几组正负相见的数值
%    randomInteger = randi([1, 10]);  %加一组范围之间的随机数
%    rereal(dp1) = 1+randomInteger/10;   %乘以随机数
%    % randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
%    rereal(dp1+5)=-rereal(dp1);   %不加随机数
%    % rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
%    dp1 = dp1+100;
%    end
%    elseif a == 2
%     for b = 1:2                             %设置加几组正负相见的数值
%    randomInteger = randi([1, 2]);  %加一组范围内之间的随机数
%    rereal(dp1) = 1*randomInteger;   %乘以随机数
%    % randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
%    % rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
%    rereal(dp1+10)=-rereal(dp1);   %不加随机数
%    dp1 = dp1+b*20;
%     end
%    else a == 3 
%     for b = 1:1                             %设置加几组正负相见的数值
% randomInteger = randi([0, 5]);  %加一组0-5之间的随机数
% rereal(dp1) = 1*randomInteger;   %乘以随机数
% %可能是这个原因 去掉加的随机数
% randomNumber = -1 + 2 * rand();  %设置加一组随机数使得并不为0.测试一下情况
% rereal(dp1+20)=-rereal(dp1)+randomNumber;     %设置相同相反的数值在附近
% rereal(dp1+20)=-rereal(dp1);
% dp1 = dp1+b*50;
%     end
% end
%    end
% subArray = rereal(1:3142);
% reversedArray = flip(subArray);
% rereal = [subArray,reversedArray];
% reim = circshift(rereal, 10)*0.8;
% end
   
re1ref=rereal+i*reim;
sizere=size(rereal);
randre=rand(sizere);
randre=(randre)*2;
re1ref=re1ref.*randre;  %%去掉随机数

figure(qq+10)
stem(ww,rereal,'b')

hold on
stem(ww,reim,'r')
hold off
LreB=(rereal>0);
LreS=(rereal<0);
rereal=LreB-LreS;
LimB=(reim>0);
LimS=(reim<0);
reim=LimB-LimS;
reim=abs(reimA).*reim;
rereal=abs(rerealA).*rereal;
re1ref=rereal+i*reim;
re1ref=re1ref*time;
%-----------------------------------

ro=ifourierzzg(ww,re1ref,nn);
ro=real(ro);
yro=intwav6(-qq,nn-1,xx);%特征样本还原成小波空间信号。
yout=yro*ro';%特征样本还原成小波空间信号。
return
function [rexx,reyy]=diffseq(xx,yy)
y1=yy(1:end-1);
y2=yy(2:end);
x1=xx(1:end-1);
x2=xx(2:end);
diffleft=(y1-y2)./(x1-x2);
diffright=(y2-y1)./(x2-x1);
diffleft=[diffleft,diffright(end)];
diffright=[diffleft(1),diffright];
muldiff=diffleft.*diffright;
Judge=(muldiff<=0);
rexx=xx.*Judge;
reyy=yy.*Judge;
return
function re=fourierzzg(ww,yy,T)%ww,yy are row vectors, T is a real number
syy=size(yy);
NN=syy(2);
nn=1:NN;
nn=nn';
nw=kron(nn,ww);
exnw=exp(-i*nw*T);
re=yy*exnw;
return
function ro=ifourierzzg(ww,yy,kk)%kk and yy are row vectors
dw=ww(2)-ww(1);
Num=(ww(end)-ww(1))/(2*pi);

kk=kk/Num;
kw=kron(kk,ww');
ekw=exp(i*kw);
ro=yy*ekw*dw/(ww(end)-ww(1));
return
function [TT,coe]=signaldecompositionsixspline(FF,qq,xx)

xxsize=size(xx)
FFsize=size(FF)
%----------------
 [KK,cc]=wavcoeobtain(xx,FF,qq);
 TT=KK*2^qq;
 coe=cc;

return
function yx=recovwav(xx,kk,coe,jj)
AA=wavelet6(jj,kk,xx);

yx=AA*coe;
return;
function yx=recovery(xx,yy1,qq)
dx=xx(2)-xx(1);
JJ=log2(1/dx);
[Nq,coe]=scalingcoeobtain(xx,yy1,qq);
%coe=coe';
yy=new6(JJ-qq,Nq,xx);
yx=yy*coe;
return;
function [Nq,coe]=scalingcoeobtain(x1,yy1,qq)
[ka,lun]=coefficientlunsixspline;
lun=lun.*2^(-(qq)/2);
T1=x1(2)-x1(1);
jj1=log2(1/T1)-qq;
begin1=x1(1);
ending1=x1(end);
kk1=ceil(2^jj1*begin1)+3:floor(2^jj1*ending1)-3;
AA1=new6(jj1,kk1,x1);
AA1=conv2(lun,AA1');
AA1=AA1';
FAA1=((AA1'*AA1)^(-1))*AA1';
Nq=kk1;
yy1=conv2(lun,yy1');
yy1=yy1';
coe=FAA1*yy1;
return;


%-------------------------------------------------------------
function [KK,cc]=wavcoeobtain(x1,yy,qq)%I0 begining number of samples, I1 ending number of samples,qq space difference between VJ and Vj0, qq=J-j0.dx sampling interval.
[ka,lun]=coefficientlunsixspline;
lun=lun.*2^(-(qq)/2);
begin1=x1(1);
ending1=x1(end);
dx1=x1(2)-x1(1);
JJ1=-log2(dx1);
jj1=JJ1-qq;
KK0=ceil((2.^jj1)*begin1+5);
KK1=floor((2.^jj1)*ending1-6);
kk1=KK0:KK1;
%-----------------
AA1=wavelet6(jj1,kk1,x1);
AA1=conv2(lun,AA1');
AA1=AA1';
FAA1=((AA1'*AA1)^(-1))*AA1';
yy=conv2(lun,yy');
KK=kk1;
size(FAA1)
size(yy)
cc=FAA1*yy';

return
function [kk,coe]=coefficientlunsixspline
coe=[0.000001067401872  -0.000001809559269   0.000003049418566  -0.000005134645802   0.000008667324169  -0.000014699552483];
coe=[coe,0.000025086665871  -0.000043129935732   0.000074752906052  -0.000130661830868   0.000230312800281  -0.000409156750070];
coe=[coe,0.000731692760992  -0.001314363043593   0.002363694999505  -0.004233427595015   0.007489059137909  -0.012906551660607];
coe=[coe,0.021150688424150  -0.031686223098915   0.041398957596618   0.954537045432857   0.041398957596618  -0.031686223098915];
coe=[coe,0.021150688424150  -0.012906551660607   0.007489059137909  -0.004233427595015   0.002363694999505  -0.001314363043593];
coe=[coe,0.000731692760992  -0.000409156750070   0.000230312800281  -0.000130661830868   0.000074752906052  -0.000043129935732];
coe=[coe,0.000025086665871  -0.000014699552483   0.000008667324169  -0.000005134645802   0.000003049418566  -0.000001809559269  0.000001067401872];
kk=-21:21;
return
function yy=new6(JJ,k,mm)%xx行k列
mm=mm';
smm=ones(size(mm));
sk=ones(size(k));
mmk=kron(mm,sk);
kmm=kron(k,smm);
xx=2^JJ*mmk-kmm;
LL1=(xx>3);
LL2=(xx<=3&xx>2);
LL3=(xx<=2&xx>1);
LL4=(xx<=1&xx>0);
LL5=(xx<=0&xx>-1);
LL6=(xx<=-1&xx>-2);
LL7=(xx<=-2&xx>-3);
LL8=(xx<=-3);
yy=LL1*0;
yy=yy+LL2.*(-(xx-3).^5./(2^3*15));
yy=yy+LL3.*(1/3*xx.^5-3*xx.^4+10*xx.^3-14*xx.^2+5*xx+17/5)/2^3;
yy=yy+LL4.*(-1/3*xx.^5+xx.^4-2*xx.^2+11/5)/2^2;
yy=yy+LL5.*(1/3*xx.^5+xx.^4-2*xx.^2+11/5)/2^2;
yy=yy+LL6.*(-1/3*xx.^5-3*xx.^4-10*xx.^3-14*xx.^2-5*xx+17/5)/2^3;
yy=yy+LL7.*((xx+3).^5./(2^3*15));
yy=yy+LL8*0;
return;
function yy=wavelet6(JJ,kk,mm)%xx行k列
fa7=-1/2554675200;fa6=1021/1277337600;fa5=-1249/19353600;
fa4=314987/255467520;fa3=-6322333/638668800;fa2=6127141/141926400;fa1=-74131711/638668800;
a0=10504567/51093504;a1=-21112517/85155840;a2=10504567/51093504;a3=-74131711/638668800;a4=6127141/141926400;
a5=-6322333/638668800;a6=314987/255467520;a7=-1249/19353600;a8=1021/1277337600;a9=-1/2554675200;
JJ=JJ+1;kk=2*kk;
yy1=new6(JJ,kk-7,mm);
yy2=new6(JJ,kk-6,mm);
yy3=new6(JJ,kk-5,mm);
yy4=new6(JJ,kk-4,mm);
yy5=new6(JJ,kk-3,mm);
yy6=new6(JJ,kk-2,mm);
yy7=new6(JJ,kk-1,mm);
yy8=new6(JJ,kk,mm);
yy9=new6(JJ,kk+1,mm);
yy10=new6(JJ,kk+2,mm);
yy11=new6(JJ,kk+3,mm);
yy12=new6(JJ,kk+4,mm);
yy13=new6(JJ,kk+5,mm);
yy14=new6(JJ,kk+6,mm);
yy15=new6(JJ,kk+7,mm);
yy16=new6(JJ,kk+8,mm);
yy17=new6(JJ,kk+9,mm);
yy=fa7*yy1+fa6*yy2+fa5*yy3+fa4*yy4+fa3*yy5+fa2*yy6+fa1*yy7+a0*yy8+a1*yy9+a2*yy10+a3*yy11+a4*yy12+a5*yy13+a6*yy14;
yy=yy+a7*yy15+a8*yy16+a9*yy17;
return;
%--------------------------------------------------------
function yy=intwav6(jj,kk,xx)%xx rows,kk columns
rr=[-0.024748803086600,  -0.004228831280609,   0.032735075614982,   0.066640824859928,  -0.207375249801789,   0.148604245264020];
rr=[rr,   -0.223595393669259,   1.332229288004729,  -3.333259781367311,   4.404035771664184,  -3.333316947049802,   1.332304643395008];
rr=[rr,   -0.223655938434532,   0.148647351384730,  -0.207410980100174,   0.066672539604844,   0.032708774573182,  -0.004207132565580];
mm=-8:1:9;
MO=ones(size(mm));
KO=ones(size(kk));
KK=kron(kk,MO');
MM=kron(KO,mm');
KM=MM+2*KK;
LL=0;
for ii=mm+9
   aa=KM(ii,:);
    bb=rr(ii)*new6(jj+1,aa,xx);
    LL=LL+bb;
end
yy=LL;
return;

