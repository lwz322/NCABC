function y=eobj(x,fun_no)


if fun_no==1
        for i=1:size(x,1)
        if x(i)<15
            y(i)=(160/15)*(15-x(i));
        else
            y(i)=40*(x(i)-15);
        end
    end
elseif fun_no==2
    for i=1:size(x,1)
        if x(i)<10
            y(i)=16*x(i);
        elseif x(i)<15
            y(i)=32*(15-x(i));
        else
            y(i)=40*(x(i)-15);
        end
    end
elseif fun_no==3
    for i=1:size(x,1)
        if x(i)<2.5
            y(i)=80*(2.5-x(i));
        elseif x(i)<5
            y(i)=64*(x(i)-2.5);
        elseif x(i)<7.5
            y(i)=64*(7.5-x(i));
        elseif x(i)<12.5
            y(i)=28*(x(i)-7.5);
        elseif x(i)<17.5
            y(i)=28*(17.5-x(i));
        elseif x(i)<22.5
            y(i)=32*(x(i)-17.5);
        elseif x(i)<27.5
            y(i)=32*(27.5-x(i));
        else 
            y(i)=80*(x(i)-27.5);
        end
    end 
elseif fun_no==4
      y=(sin(5*pi*x)).^6;
elseif fun_no==5
    y=exp(-2*log(2)*((x-0.1)/0.8).^2).*(sin(5*pi*x)).^6;
elseif fun_no==6
     y=(sin(5*pi*(x.^0.75-0.05))).^6;
elseif fun_no==7
    y=exp(-2*log(2)*((x-0.08)/0.854).^2).*(sin(5*pi*(x.^(3/4)-0.05))).^6;
elseif fun_no==8
     y=200-(x(:,1).*x(:,1)+x(:,2)-11).^2-(x(:,1)+x(:,2).*x(:,2)-7).^2;
elseif fun_no==9
    factor1=(4-2.1*(x(:,1).^2)+(x(:,1).^4)/3).*(x(:,1).^2)+x(:,1).*x(:,2);
    factor2=(-4+4*(x(:,2).^2)).*(x(:,2).^2);
    y=-4*(factor1+factor2);
elseif fun_no==10    
temp=zeros(size(x,1),1);
for i=0:24,
    a=(mod(i,5)-2)*16;b=(floor(i/5)-2)*16;
    temp=temp+1./(1+i+(x(:,1)-a).^6+(x(:,2)-b).^6);
end
y=500-1./(0.002+temp);

elseif fun_no==11
    a=size(x,1);
s1=zeros(a,1);
s2=zeros(a,1);

for i=1:5
    s1=s1+i*cos((i+1).*x(:,1)+i);
    s2=s2+i*cos((i+1).*x(:,2)+i);
end
y=-s1.*s2;

elseif fun_no==12
a=size(x,1);
s1=zeros(a,1);
s2=zeros(a,1);
s3=zeros(a,1);

for i=1:5
    s1=s1+i*cos((i+1).*x(:,1)+i);
    s2=s2+i*cos((i+1).*x(:,2)+i);
    s3=s3+i*cos((i+1).*x(:,3)+i);
end
y=-s1.*s2.*s3;

elseif fun_no==13
   a=size(x,1);
s1=zeros(a,1);
s2=zeros(a,1);
s3=zeros(a,1);
s4=zeros(a,1);

for i=1:5
    s1=s1+i*cos((i+1).*x(:,1)+i);
    s2=s2+i*cos((i+1).*x(:,2)+i);
    s3=s3+i*cos((i+1).*x(:,3)+i);
    s4=s4+i*cos((i+1).*x(:,4)+i);
end
y=-s1.*s2.*s3.*s4; 

elseif fun_no==14
    a=log(x);
    y=sin(10*a);
    
elseif fun_no==15
    a1=log(x(:,1));
    a2=log(x(:,2));
    s1=sin(10*a1);
    s2=sin(10*a2);
    y=0.5*(s1+s2);
elseif fun_no==16
    a1=log(x(:,1));
    a2=log(x(:,2));
    a3=log(x(:,3));
    s1=sin(10*a1);
    s2=sin(10*a2);
    s3=sin(10*a3);
    y=(s1+s2+s3)/3;
end















