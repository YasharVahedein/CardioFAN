function GEOMETRY(LL,ALF,N,NODE_CONNECT,ARZ1,ARZ4,NNODE,threshold,threshold1)
%**************
% Copyright (c) Rochester Institute of Technology (RIT) 2018 <Yashar Seyed Vahedein, Alexander Liberson>
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% The user is recommended to reference the first released publication based on this code:
% 
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%**************
Xin(1:N)=0;
Yin(1:N)=0;
for iv=1:N
    for i=1:NNODE
        if NODE_CONNECT(i,3)==iv
            Xend(iv)=Xin(iv)+LL(iv)*cos(ALF(iv));
            Yend(iv)=Yin(iv)+LL(iv)*sin(ALF(iv));
            if Xin(iv)==0 && Yin(iv)==0
                redo(iv)=1;
            else
                redo(iv)=0;
            end
        elseif NODE_CONNECT(i,4)==iv
            Xend(iv)=Xin(iv)+LL(iv)*cos(ALF(iv));
            Yend(iv)=Yin(iv)+LL(iv)*sin(ALF(iv));
            if Xin(iv)==0 && Yin(iv)==0
                redo(iv)=1;
            else
                redo(iv)=0;
            end
        end
    end
    for i=1:NNODE
        if NODE_CONNECT(i,2)==iv && NODE_CONNECT(i,3)~=0
            Xin(NODE_CONNECT(i,3))=Xend(iv);
            Yin(NODE_CONNECT(i,3))=Yend(iv);
            if (NODE_CONNECT(i,4)~=0)
                Xin(NODE_CONNECT(i,4))=Xend(iv);
                Yin(NODE_CONNECT(i,4))=Yend(iv);
            end
        end
    end
end

%make sure that geometry is correct
for iv=1:N
    if redo(iv)==1
        Xend(iv)=Xend(iv)+Xin(iv);
        Yend(iv)=Yend(iv)+Yin(iv);
    end
end



figure (2)
counter=0;
for iv=1:N
    if ARZ1(iv)>threshold && abs(Xend(iv)-Xin(iv))>0.02
        plot([Yin(iv) Yend(iv)],[-Xin(iv) -Xend(iv)],'color','k','LineWidth',4);
        yy=Xin(iv)+(Xend(iv)-Xin(iv))/2;
        xx=Yin(iv)+(Yend(iv)-Yin(iv))/2;
        counter=counter+1;
        text(xx,-yy,num2str(iv),'VerticalAlignment','top')
        hold on
    elseif ARZ1(iv)>threshold && abs(Xend(iv)-Xin(iv))<=0.02
        plot([Yin(iv) Yend(iv)],[-Xin(iv) -Xend(iv)],'color','k','LineWidth',4);
        yy=Xin(iv)+(Xend(iv)-Xin(iv))/2;
        xx=Yin(iv)+(Yend(iv)-Yin(iv))/2;
        counter=counter+1;
        text(xx,-yy,num2str(iv),'HorizontalAlignment','Left')
        hold on
    end
    if ARZ1(iv)>threshold1 && ARZ1(iv)<=threshold && abs(Xend(iv)-Xin(iv))>0.02
        plot([Yin(iv) Yend(iv)],[-Xin(iv) -Xend(iv)],'color','k','LineWidth',2);
        yy=Xin(iv)+(Xend(iv)-Xin(iv))/2;
        xx=Yin(iv)+(Yend(iv)-Yin(iv))/2;
        counter=counter+1;
        text(xx,-yy,num2str(iv),'VerticalAlignment','top')
        hold on
    elseif ARZ1(iv)>threshold1 && ARZ1(iv)<=threshold && abs(Xend(iv)-Xin(iv))<=0.02
        plot([Yin(iv) Yend(iv)],[-Xin(iv) -Xend(iv)],'color','k','LineWidth',2);
        yy=Xin(iv)+(Xend(iv)-Xin(iv))/2;
        xx=Yin(iv)+(Yend(iv)-Yin(iv))/2;
        counter=counter+1;
        text(xx,-yy,num2str(iv),'HorizontalAlignment','Left')
        hold on
    end
    if ARZ1(iv)<=threshold1 && abs(Xend(iv)-Xin(iv))>0.02
        plot([Yin(iv) Yend(iv)],[-Xin(iv) -Xend(iv)],'color','k','LineWidth',1);
        yy=Xin(iv)+(Xend(iv)-Xin(iv))/2;
        xx=Yin(iv)+(Yend(iv)-Yin(iv))/2;
        counter=counter+1;
        text(xx,-yy,num2str(iv),'VerticalAlignment','top')
        hold on
    elseif ARZ1(iv)<=threshold1 && abs(Xend(iv)-Xin(iv))<=0.02
        plot([Yin(iv) Yend(iv)],[-Xin(iv) -Xend(iv)],'color','k','LineWidth',1);
        yy=Xin(iv)+(Xend(iv)-Xin(iv))/2;
        xx=Yin(iv)+(Yend(iv)-Yin(iv))/2;
        text(xx,-yy,num2str(iv),'HorizontalAlignment','Left')
        hold on
    end
    Casenum = int2str(N);
    title([Casenum, ' Vessel Model Geometry'])
    hold on
end
grid on

end
