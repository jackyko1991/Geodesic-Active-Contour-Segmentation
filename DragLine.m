function h = DragLine(f,aH,h)   
%鼠标改变直线方向，
    global xy_start
    global xdata_start
    global ydata_start
    set(h,'ButtonDownFcn',@startDragFcn)
    set(f,'WindowButtonUpFcn',@stopDragFcn);
    function startDragFcn(varargin)
        set(f,'WindowButtonMotionFcn',@draggingFcn);
        xy_start=get(aH,'CurrentPoint');
        xdata_start=get(h,'xdata');
        ydata_start=get(h,'ydata');
    end
    function draggingFcn(varargin)
        xy_current=get(aH,'CurrentPoint');
%         dx=xy_current(1)-xy_start(1)
%         dy=xy_current(3)-xy_start(3)
        xdata_start(2) = xy_current(1);
        ydata_start(2) = xy_current(3);
%         set(h,'xdata',xdata_start+dx,'ydata',ydata_start+dy);
        set(h,'xdata',xdata_start,'ydata',ydata_start);
    end
    function stopDragFcn(varargin)
        set(f,'WindowButtonMotionFcn','')
    end
end
