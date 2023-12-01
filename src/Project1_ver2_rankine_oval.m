% The following code was contributed by QIN Qijun and LIN Ju. 
% It is intended for use only in the course AAE3001 Fundamentals of Aerodynamics (Sem 1 2023).

syms x y;

V_inf = 10;
lamda = 20;

a = sqrt(1+20/(10*pi));
b = 1;
% c = sqrt(a^2 - b^2);
% c = 0.646268;
c = 0.6396466943575;
e = b/a;


% oval = x^2/a^2 + y^2/c^2 - 1;
oval = 10*y + 10/pi *(atan(y/(x+1)) - atan(y/(x-1))) - 0;

psi(x, y) = V_inf*y + (lamda/(2*pi))*(atan(y/(x+b)) - atan(y/(x-b)));

u(x, y) = diff(psi(x, y), y);
v(x, y) = -diff(psi(x, y), x);

V(x, y) = sqrt(u^2 + v^2);

C_p(x, y) = 1 - (V/V_inf)^2;


syms theta;

x_(theta) = a*cos(theta);
y_(theta) = c*sin(theta);

u_(theta) = u(x_, y_);
v_(theta) = v(x_, y_);

V_(theta) = sqrt(u_^2 + v_^2);

C_p_(theta) = 1 - (V_/V_inf)^2;


%--------------------------------------------------------------------------
% Draw Rankine Oval

end_at_stag = true;

x_plot_oval = [-a];
y_plot_oval = [0];

precision_improve = 1;

precision_xplot = 0.0001;
precision_yplot = precision_xplot;
% x_plot = [-2];


x_plot = [-a];
y_plot = [precision_yplot];

x_i = -a;
y_i = precision_yplot;

while x_i <= -a/2
    
    streamLine_diff_inverse = streamLine_diff_inverse_func(x_i, y_i);
    x_i = x_i + streamLine_diff_inverse * precision_yplot;
    y_i = y_i + precision_yplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

end

while -a/2 < x_i && x_i <= a/2
    
    streamLine_diff = streamLine_diff_func(x_i, y_i);
    x_i = x_i + precision_xplot;
    y_i = y_i + streamLine_diff * precision_xplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

end

precision_xplot = precision_xplot / precision_improve;
precision_yplot = precision_xplot;

if(end_at_stag)
    while a/2 < x_i && x_i < a && y_i > 0
    
    streamLine_diff_inverse = streamLine_diff_inverse_func(x_i, y_i);
    x_i = x_i - streamLine_diff_inverse * precision_yplot;
    y_i = y_i - precision_yplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

    end
else
    while a/2 < x_i && y_i > 0
    
    streamLine_diff_inverse = streamLine_diff_inverse_func(x_i, y_i);
    x_i = x_i - streamLine_diff_inverse * precision_yplot;
    y_i = y_i - precision_yplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

    end
end

precision_xplot = precision_xplot * precision_improve;
precision_yplot = precision_xplot;

if(end_at_stag && x_i >= a)
    while(y_i > 0)
        x_i = x_i;
        y_i = y_i - precision_yplot;
        x_plot = [x_plot, x_i];
        y_plot = [y_plot, y_i];
    end
end

x_plot = [x_plot, a];
y_plot = [y_plot, 0];

x_plot_oval = [x_plot_oval, x_plot];
y_plot_oval = [y_plot_oval, y_plot];
x_plot = [];
y_plot = [];


x_plot = [-a];
y_plot = [-precision_yplot];

x_i = -a;
y_i = -precision_yplot;

while x_i <= -a/2
    
    streamLine_diff_inverse = streamLine_diff_inverse_func(x_i, y_i);
    x_i = x_i - streamLine_diff_inverse * precision_yplot;
    y_i = y_i - precision_yplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

end

while -a/2 < x_i && x_i <= a/2
    
    streamLine_diff = streamLine_diff_func(x_i, y_i);
    x_i = x_i + precision_xplot;
    y_i = y_i + streamLine_diff * precision_xplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

end

precision_xplot = precision_xplot / precision_improve;
precision_yplot = precision_xplot;

if(end_at_stag)
    while a/2 < x_i && x_i < a && y_i < 0
    
    streamLine_diff_inverse = streamLine_diff_inverse_func(x_i, y_i);
    x_i = x_i + streamLine_diff_inverse * precision_yplot;
    y_i = y_i + precision_yplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

    end
else
    while a/2 < x_i && y_i < 0
    
    streamLine_diff_inverse = streamLine_diff_inverse_func(x_i, y_i);
    x_i = x_i + streamLine_diff_inverse * precision_yplot;
    y_i = y_i + precision_yplot;
    x_plot = [x_plot, x_i];
    y_plot = [y_plot, y_i];

    end
end

precision_xplot = precision_xplot * precision_improve;
precision_yplot = precision_xplot;

if(end_at_stag && x_i >= a)
    while(y_i < 0)
        x_i = x_i;
        y_i = y_i + precision_yplot;
        x_plot = [x_plot, x_i];
        y_plot = [y_plot, y_i];
    end
end

x_plot = flip(x_plot);
y_plot = flip(y_plot);

x_plot_oval = [x_plot_oval, x_plot];
y_plot_oval = [y_plot_oval, y_plot];


plot(x_plot_oval, y_plot_oval, 'b', 'LineWidth',1.5);

%--------------------------------------------------------------------------



title('Rankine Oval');
xlabel('x');
ylabel('y');
% xlim([x_plot(1), x_plot(end)]);
% ylim([x_plot(1), x_plot(end)]);
axis equal;
axis([-2 2 -2 2]);

legend('Rankine Oval');

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% xlim([-2.5, 2.5]);
% ylim([-2.5, 2.5]);
% zlim([-1, 1]);
% zticks([-1:0.25:1])
% zticks('manual');
grid on
% shift_axis_to_origin( gca );


function streamLine_diff = streamLine_diff_func(x, y)
    streamLine_diff = v_func(x, y) / u_func(x, y);
end

function streamLine_diff_inverse = streamLine_diff_inverse_func(x, y)
    streamLine_diff_inverse = u_func(x, y) / v_func(x, y);
end

function u_f = u_func(x, y)
    % u_f = 1791925356007081/(562949953421312*(y^2/(x + 1)^2 + 1)*(x + 1)) - 1791925356007081/(562949953421312*(y^2/(x - 1)^2 + 1)*(x - 1)) + 10;
    % u_f = vpa( (10 + (10/pi)*((x+1)/((x+1)^2 + y^2) - (x-1)/((x-1)^2 + y^2))) , 10);
    u_f = (10 + (10/pi)*((x+1)/((x+1)^2 + y^2) - (x-1)/((x-1)^2 + y^2)));
end

function v_f = v_func(x, y)
    % v_f = (1791925356007081*y)/(562949953421312*(y^2/(x + 1)^2 + 1)*(x + 1)^2) - (1791925356007081*y)/(562949953421312*(y^2/(x - 1)^2 + 1)*(x - 1)^2);
    % v_f = vpa( (10/pi)*(y/((x + 1)^2 + y^2) - y/((x-1)^2 + y^2)) , 10);
    v_f = (10/pi)*(y/((x + 1)^2 + y^2) - y/((x-1)^2 + y^2));
end