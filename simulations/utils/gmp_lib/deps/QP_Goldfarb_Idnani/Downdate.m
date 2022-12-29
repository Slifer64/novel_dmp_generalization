function [R,J,q]=Downdate(R,J,q,pos)
m=size(R,2);
    tmp=R(:,m);
    R(:,pos)=[];R=[R,zeros(size(tmp))];
    [Q,R]=qr(R);
    J=J*Q;
    q=q-1;