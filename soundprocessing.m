function zzg
% % % 原音频
% [aa,bb]=audioread('D:\MATLAB\bin\TTL.m4a');
% %提取有效音频信号
% aa=aa(:,1);
% % sound(aa,bb)
% ma=max(aa);
% aa=aa/ma;
% cc=aa(176700:197100);
% dd=aa(284500:311400);
% cc=cc/max(cc);
% cc=cc(abs(cc)>0.06);
% dd=dd/max(dd);
% dd=dd(abs(dd)>0.06);
% % sound(cc,bb);
% %sound(dd,bb);
% % ------------------------------------------------------------------------

%这部分为自己录的o的音频信号 补0法
[aa,bb]=audioread('D:\matlab2020a\bin\iii.M4A'); %连续音频信号1
aa=aa(:,1);%取一个声道的音频
% plot(aa);%绘制原声道
ma=max(aa);
aa=aa/ma;
% % cc=aa(40000:55000); %截取男生有效音频
% % % cc=aa(100500:115500); %试试女生的
% % dd=aa(100500:115500); %截取女生有效音频  o部分
%%%%----------------------------------------------
%cc=aa(35000:55000);
%dd=aa(93000:113000);    %e部分
%%%%%-------------------------------------------
cc=aa(40000:55000);   %取第一部分  因为两部分幅值貌似不太一样
dd=aa(70000:85000); %取第二部分    %i部分
%%%%-----------------------------
cc=cc/max(cc);
%%cc=cc(abs(cc)>0.06); %%改为小于0.06的赋值为0
[m,n]=size(cc); %得到矩阵的行号和列号
for i=1:m;
for j=1:n;
if(abs(cc(i,j))<0.06);     %设定阈值
cc(i,j)=0;;
end
end
end
dd=dd/max(dd);
%%dd=dd(abs(dd)>0.06); %%改为小于0.06的赋值为0
[m,n]=size(dd); %得到矩阵的行号和列号
for i=1:m;
for j=1:n;
if(abs(dd(i,j))<0.06);     %设定阈值
dd(i,j)=0;;
end
end
end
% sound(cc,bb)
 %sound(dd,bb)
%  %%% 部分结束

% % % 不补0试试
% [aa,bb]=audioread('D:\MATLAB\bin\e2.M4A');
% %提取有效音频信号
% aa=aa(:,1);
% % sound(aa,bb)
% ma=max(aa);
% aa=aa/ma;
% cc=aa(35000:55000);
% dd=aa(93000:113000);
% cc=cc/max(cc);
% cc=cc(abs(cc)>0.06);
% dd=dd/max(dd);
% dd=dd(abs(dd)>0.06);
% % sound(cc,bb);
% %sound(dd,bb);
% % ------------------------------------------------------------------------

%小波变换提取系数
qq=5;
sFF=size(cc);
xx=0:1:sFF(1)-1;
% cc=wavelet6(0,0,xx);
[TT,coe1]=signaldecompositionsixspline(cc,qq,xx);
kk=TT/2^qq;
sdd=size(dd);
 xd=0:1:sdd(1)-1;

[TT2,coe2]=signaldecompositionsixspline(dd,qq,xd);
kk2=TT2/2^qq;

%%%%%%部分结束

% sscoe1=size(coe1);
% sscoe2=size(coe2);
% figure(1)
% stem(TT,coe1)
% hold on
% stem(TT2,coe2,'r')
% hold off
%-------------------------------------



%%%取奇次列特征样本
xa=0:2^4:xx(end);
xa=xa(2:2:end);
yx=recovwav(xa,kk,coe1,-qq);
yx2=recovwav(xa,kk2,coe2,-qq);
%%部分结束

% figure(2)
% stem(xa,yx,'r*')
% hold on
% stem(xa,yx2,'b*')
% hold off
%-----------------------------------


% syx=size(yx);
% kx=0:syx(1)-1;
% yy=intwav6(-qq,kx,xx);
% yy=yy*yx;
% % sound(yy,bb)
ym=recovwav(xx,kk,coe1,-qq);
% figure(4)
% plot(xx,yy)
% hold on
% plot(xx,ym)
% hold off
%--------------------------------------------------------------

%%傅里叶变换分析频谱
ww=-pi:0.0005:pi;
yx=yx';
 yx2=yx2';
re1=fourierzzg(ww,yx,1);
syx=size(yx);
NN=syx(2);
nn=1:NN;
re1b=abs(re1);
mnre1=1.2*mean(re1b); %改变平均值倍数
are1=re1b.*(re1b>mnre1); %raw
% are1=re1b.*(re1b>20); %粗暴去峰值
re1ref=re1.*(re1b>mnre1);%raw
% re1ref=re1.*(re1b>20);%粗暴去峰值
% re1ref=0; %测试

% %%去除第一个尖峰
% startIndex0 = 1187; %%定义数据起始量
% endIndex0 = 1218;   %%定义数组结束
% re1ref(startIndex0:endIndex0) = 0;
% startIndex0 = 1187; %%定义数据起始量
% endIndex0 = 1218;   %%定义数组结束
% are1(startIndex0:endIndex0) = 0;
% %去除对称部分
% startIndex1 = 11350; %%定义数据起始量
% endIndex1 = 11381;   %%定义数组结束
% re1ref(startIndex1:endIndex1) = 0;
% startIndex1 = 11350; %%定义数据起始量
% endIndex1 = 11381;   %%定义数组结束
% are1(startIndex1:endIndex1) = 0;
% %%部分完毕
% 
% 
% 
% %去除第二个尖峰
% startIndex2 = 1297; %%定义数据起始量
% endIndex2 = 1395;   %%定义数组结束
% re1ref(startIndex2:endIndex2) = 0;
% are1(startIndex2:endIndex2) = 0;
% %去除对称部分
% startIndex3 = 11174; %%定义数据起始量
% endIndex3 = 11271;   %%定义数组结束
% re1ref(startIndex3:endIndex3) = 0;
% are1(startIndex3:endIndex3) = 0;
% %%部分完毕
% 
% %去除第三个尖峰
% startIndex4 = 1454; %%定义数据起始量
% endIndex4 = 1527;   %%定义数组结束
% re1ref(startIndex4:endIndex4) = re1ref(1,1527);
% are1(startIndex4:endIndex4) = are1(1,1527);
% %去除对称部分
% startIndex5 = 11041; %%定义数据起始量
% endIndex5 = 11114;   %%定义数组结束
% re1ref(startIndex5:endIndex5) = re1ref(1,11041);
% are1(startIndex5:endIndex5) = are1(1,11041);
% %%部分完毕
% 
% %去除第五个尖峰（跳过最大的）
% startIndex6 = 1640; %%定义数据起始量
% endIndex6 = 1723;   %%定义数组结束
% re1ref(startIndex6:endIndex6) = re1ref(1,1640);
% are1(startIndex6:endIndex6) = are1(1,1640);
% %去除对称部分
% startIndex7 = 10845; %%定义数据起始量
% endIndex7 = 10928;   %%定义数组结束
% re1ref(startIndex7:endIndex7) = re1ref(1,10928);
% are1(startIndex7:endIndex7) = are1(1,10928);
% %%部分完毕
% 
% %去除尖峰最大
% startIndex8 = 1536; %%定义数据起始量
% endIndex8 = 1640;   %%定义数组结束
% re1ref(startIndex8:endIndex8) = re1ref(1,1640);
% are1(startIndex8:endIndex8) = are1(1,1640);
% %去除对称部分
% startIndex9 = 10928; %%定义数据起始量
% endIndex9 = 11033;   %%定义数组结束
% re1ref(startIndex9:endIndex9) = re1ref(1,10928);
% are1(startIndex9:endIndex9) = are1(1,10928);
% %%部分完毕
% 
% %%去除末尾
% startIndex0 = 4216; %%定义数据起始量
% endIndex0 = 4357;   %%定义数组结束
% re1ref(startIndex0:endIndex0) = 0;
% are1(startIndex0:endIndex0) = 0;
% %去除对称部分
% startIndex1 = 8212; %%定义数据起始量
% endIndex1 = 8353;   %%定义数组结束
% re1ref(startIndex1:endIndex1) = 0;
% are1(startIndex1:endIndex1) = 0;
% %%部分完毕


 ro=ifourierzzg(ww,re1ref,nn);  %还原到特征样本空间
 ro=real(ro);
 yro=intwav6(-qq,nn,xx);%特征样本还原成小波空间信号。
yout=yro*ro';%特征样本还原成小波空间信号。
 figure(5)
 stem(nn,ro,'r*')
 hold on
 stem(nn,(yx),'b*')
 hold off
 figure(6)
plot(xx,ym,'r-')
hold on
plot(xx,yout,'b-')
hold off
 sound(cc,bb);
 pause(5)
 sound(yout*80,bb)
%------------------------------------------------
re2=fourierzzg(ww,yx2,1);
re2=abs(re2);
mnre2=mean(re2);
are2=re2.*(re2>mnre2);
% 
% size(re2)
% size(re1)
figure(3)
plot(ww,abs(re1),'r-')
hold on
plot(ww,are1,'k-')
%plot(ww,are2,'b-')
hold off
% hold on
% plot(ww,re2,'k-')
% hold off
% --------------------------------
% y1=recovwav(xx,kk,coe1,-qq);
% y2=recovwav(xx,kk2,coe2,-qq);
% sound(y1,bb);
% pause(5)
% ii=100000
% sound(y2,bb);


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

