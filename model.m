clc; clear;

%As velocity increases, the rate of bonds breaking over time increases
Vi = 5; %cm/s
Vo = 15; %cm/s
deltaP = 2; % dPa
L = 10; %cm
u = 0.045; %g/(cm*s)
z = 3; %cm
Vol_t = 50; %cm^3

dur = 1000; %ms
RBC = 0.93;%ratio of RBC to total volume
P = 5.2; %kPa
R = [100];

a = 0.9998; %constant for RBC decay rate during compresion
b = 0.001; %constant deteriming the effect of Pressure on decay rate

t = [0:1:dur];

Vz = abs(deltaP/(2*u*L)*z^2 + ((Vo-Vi)/L-deltaP/(2*u))*z + Vi);

dBdt = -Vz/9000;
RBC = 0.93;

Bi = Vol_t*(1-RBC); %cm^3: initial fibrinogen volume
B = [Bi];

for i = 1:dur

    dRdt = R(end)*(a*exp(b*(5.2-P))-1);

    if B(end) > 0.35*Bi
        B = [B B(end)+dBdt+dRdt];

    elseif B(end) <= 0
        B = [B B(end)];

    else
        B = [B B(end)+dRdt];
    end
end

B1 = B;

RBC = 0.9;
Bi2 = Vol_t*(1-RBC); %cm^3: initial fibrinogen volume
B = [Bi2];

for i = 1:dur

    dRdt = R(end)*(a*exp(b*(5.2-P))-1);

    if B(end) > 0.35*Bi2
        B = [B B(end)+dBdt+dRdt];
S
    elseif B(end) <= 0
        B = [B B(end)];

    else
        B = [B B(end)+dRdt];
    end
end

B2 = B;

RBC = 0.8;
Bi3 = Vol_t*(1-RBC); %cm^3: initial fibrinogen volume
B = [Bi3];

for i = 1:dur

    dRdt = R(end)*(a*exp(b*(5.2-P))-1);

    if B(end) > 0.35*Bi3
        B = [B B(end)+dBdt+dRdt];

    elseif B(end) <= 0
        B = [B B(end)];

    else
        B = [B B(end)+dRdt];
    end
end

B3 = B;

plot(t,B1/Bi*100)
hold on
plot(t,B2/Bi2*100)
hold on
plot(t,B3/Bi3*100)
title('Thrombus dissolution over time vs. RBC ratio')
ylabel('Thrombus integrity (%)')
xlabel('time (ms)')
legend('0.93','0.90','0.80')

% B1 = B;
% 
% Vz = 40;
% dBdt = -Vz/9000;
% B = [Bi];
% for i = 1:dur
%     if B(end) > 0.35*Bi
%         B = [B B(end)+dBdt];
%     else
%         B = [B B(end)];
%     end
% end
% 
% B2 = B;
% 
% Vz = 80;
% dBdt = -Vz/9000;
% B = [Bi];
% for i = 1:dur
%     if B(end) > 0.35*Bi
%         B = [B B(end)+dBdt];
%     else
%         B = [B B(end)];
%     end
% end
% 
% B3 = B;
% 
% plot(t,B1)
% hold on
% plot(t,B2)
% hold on
% plot(t,B3)
% title('Fibrinogen bond decay over time at V = 20 cm/s, 40 cm/s, and 80 cm/s')
% xlabel('time (ms)')
% ylabel('Fibrinogen Volume (cm^3)')
% legend('V = 20 cm/s','V = 40 cm/s', 'V = 80 cm/s')

% R1 = R;
% 
% P = 6.2;
% R = [100];
% 
% for i = 1:dur
%     dRdt = R(end)*(a*exp(b*(5.2-P))-1);
%     R = [R R(end)+dRdt];
% end
% 
% R2 = R;
% 
% P = 8.2;
% R = [100];
% 
% for i = 1:dur
%     dRdt = R(end)*(a*exp(b*(5.2-P))-1);
%     R = [R R(end)+dRdt];
% end
% 
% R3 = R;

% plot(t,R1)
% hold on
% plot(t,R2)
% hold on
% plot(t,R3)
% legend('5.2 kPa','6.2 kPa','8.2 kPa')
% title('RBC decay over time at compression P = 5.2kPa, 6.2kPa, and 8.2kPa')
% xlabel('Time (ms)')
% ylabel('Alive RBC (%)')


% Vz = abs(deltaP/(2*u*L)*z^2 + ((Vo-Vi)/L-deltaP/(2*u))*z + Vi);
% plot(deltaP,Vz)
% title('Magnitude of velocity at z = 3 cm, L = 10 cm, Vi = 5 cm/s, and Vo = 15cm/s, as a function of pressure')
% xlabel('Pressure (dPa)')
% ylabel('Velocity (cm/s)')


% plot(z,Vz)
% title('Magnitude of velocity at deltaP = 2 dPa, L = 10 cm, Vi = 5 cm/s, and Vo = 15cm/s, as a function of thrombus position')
% ylabel('Veloctiy (cm/s)')
% xlabel('Thrombus Position (cm)')


