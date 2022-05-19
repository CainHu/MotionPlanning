% 读取数据
% uiopen('H:\CLionProjects\MotionPlanning\demo5\matlab_code\data.csv',1);

[m, n] = size(data);

num_steps = 100;
u_buf = (0:num_steps-1)/num_steps;

% 这里的vel, acc, jerk, snap均为以u为单位的
% 以时间为单位需要分别乘上相应的Δt的幂
% 但这里只考察连续性，所以忽略了此系数
pos_buf = zeros([m*num_steps, 1]);
vel_buf = zeros([m*num_steps, 1]);
acc_buf = zeros([m*num_steps, 1]);
jerk_buf = zeros([m*num_steps, 1]);
snap_buf = zeros([m*num_steps, 1]);
time_buf = (1:m*num_steps)'/num_steps;

for i = 1:m
    for j = 1:num_steps
        [pos, vel, acc, jerk, snap] = bezier5(u_buf(j), data{i,:});
        pos_buf(j + (i - 1) * num_steps, 1) = pos;
        vel_buf(j + (i - 1) * num_steps, 1) = vel;
        acc_buf(j + (i - 1) * num_steps, 1) = acc;
        jerk_buf(j + (i - 1) * num_steps, 1) = jerk;
        snap_buf(j + (i - 1) * num_steps, 1) = snap;
    end
end

figure(1)
plot(time_buf, pos_buf)
title('position')
grid on;

figure(2)
plot(time_buf, vel_buf)
title('velocity')
grid on;

figure(3)
plot(time_buf, acc_buf)
title('acceleration')
grid on;

figure(4)
plot(time_buf, jerk_buf)
title('jerk')
grid on;

figure(5)
plot(time_buf, snap_buf)
title('snap')
grid on;


