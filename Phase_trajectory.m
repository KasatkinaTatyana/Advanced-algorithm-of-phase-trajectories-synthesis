clear all
close all
clc

global constr
% 1) y' < constr

Y_0 = [0.1 0.2 0.5];
Y_end = [2.5 0.5 0.5];

N = 1000;

% —троитс€ фазова€ крива€ вида
% Psi = c0 + c1*(y - y0) + c2*(y - y0)^2 + c3*(y - y0)^3 + d*(y - y0)^2*(y - yend)^2,
% соедин€юща€
% начальное положение Y_0 и конечное положение Y_end
delta = Y_end(1) - Y_0(1);

c0 = Y_0(2);
c1 = Y_0(3)/Y_0(2);

Ma = [delta^2     delta^3;
      2*delta     3*delta^2];
  
Mb = [Y_end(2) - c1*delta - c0;
      Y_end(3)/Y_end(2) - c1];
Mc = inv(Ma)*Mb;

c2 = Mc(1);
c3 = Mc(2);

Expr = c0 + c1*delta + c2*delta^2 + c3*delta^3;

Expr_1 = c1 + 2*c2*delta + 3*c3*delta^2;

d = 0;

figure();
hold on; grid on;

dy = (Y_end(1) - Y_0(1)) / N;
y=Y_0(1):dy:Y_end(1);

Psi = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1)).^2 + c3*(y - Y_0(1)).^3;

plot(y,Psi);

constr = 0.9;

N_shift = 1;

cond = false;

d_min = -100;
d_max = 100;
hd = 0.01;

% while not(cond)
%     flag = 0;
%     % »щем границы интервала, на котором нарушаетс€ ограничение
%     for i=1:(N+1)
%         if ((flag == 0)&&(Psi(i) >= constr))
%             y_left(1) = Y_0(1) + (i-1 - N_shift)*dy;  % отступили N_shift шагов назад от критической точки
%             y_left(2) = Psi(i - N_shift);
%             y_left(3) = (c1 + 2*c2*(y_left(1) - Y_0(1)) + 3*c3*(y_left(1) - Y_0(1))^2)*y_left(2);
%             flag = 1;
%         end
%         if ((flag==1)&&(Psi(i) < constr))
%             y_right(1) = Y_0(1) + (i-1 + N_shift)*dy;
%             y_right(2) = Psi(i + N_shift);
%             y_right(3) = (c1 + 2*c2*(y_right(1) - Y_0(1)) + 3*c3*(y_right(1) - Y_0(1))^2)*y_right(2);
%             break;
%         end
%     end
%     
%     for d=d_min:hd:d_max
%         cond = IsCurveExist(y_left, y_right, dy, d);
%         if (cond == true)
%             break;
%         end
%     end
%     
%     if not(cond)
%         N_shift = N_shift + 1;
%     end
% end


% -------------∆елание добитьс€ более плавного поведени€ кривой------------
% -------------ƒл€ этого эмпирически введено ограничение y_left(1) <= 0.4)
N_shift = 1;

while (true)
    flag = 0;
    for i=1:(N+1)
        if ((flag == 0)&&(Psi(i) >= constr))
            y_left(1) = Y_0(1) + (i-1 - N_shift)*dy;  % отступили N_shift шагов назад от критической точки
            y_left(2) = Psi(i - N_shift);
            y_left(3) = (c1 + 2*c2*(y_left(1) - Y_0(1)) + 3*c3*(y_left(1) - Y_0(1))^2)*y_left(2);
            flag = 1;
        end
        if ((flag==1)&&(Psi(i) < constr))
            y_right(1) = Y_0(1) + (i-1 + N_shift)*dy;
            y_right(2) = Psi(i + N_shift);
            y_right(3) = (c1 + 2*c2*(y_right(1) - Y_0(1)) + 3*c3*(y_right(1) - Y_0(1))^2)*y_right(2);
            break;
        end
    end
    if (y_left(1) <= 0.4)
        break;
    end
    N_shift = N_shift + 1;
end

for d=d_min:hd:d_max
    cond = IsCurveExist(y_left, y_right, dy, d);
    if (cond == true)
        break;
    end
end

if (cond)
    % Ќаходим интервал по d, на котором ограничение не нарушаетс€
    d_lims = d_interval(y_left, y_right, dy, d);
    % —троим корректирующий отрезок при максимальном d
    replace_part = curve_synthesis(y_left, y_right, dy, d_lims(2));
    x = y_left(1):dy:y_right(1);
    plot(x,replace_part,'r');
else
    disp('јлгоритм завершил свою работу');
end




