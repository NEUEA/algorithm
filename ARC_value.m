function v=ARC_value(point)
    global B_Angle;
    global pro;
    m=pro.od;
    for i=1:m-1
        up(i)=B_Angle(i,point(m+i)+1);
        down(i)=B_Angle(i,point(m+i));
    end
    w=max(sum(up-point(2:m)),sum(point(2:m)-down));
%     v=(1+w/(sum(up-down)/(m-1)))*point(1);
    v=(2+w/sum(up-down))*point(1);
end