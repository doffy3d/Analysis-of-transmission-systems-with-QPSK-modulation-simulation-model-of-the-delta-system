clc, clear all, close all;

Ebpn = 0:7;

% teorijski
Pb_teor = erfc(sqrt(10.^(Ebpn/10)))/2;

figure, semilogy(Ebpn,Pb_teor,'LineWidth',1.5);
grid on
xlabel('E_b/p_N [dB]');
ylabel('BER');
hold on
% merenja simulacije


Pb_rrc_r1 = [0.0937 0.0700 0.0522 0.0362 0.0236 0.0140 0.0081 0.0043];
Pb_rrc_r2 = [0.0907 0.0679 0.0485 0.0324 0.0199 0.0107 0.0061 0.0027];
Pb_rrc_r3 = [0.0879 0.0646 0.0445 0.0287 0.0164 0.0084 0.0042 0.0015];


Pb_rrc_z_r1 = [0.0602 0.0368 0.0201 0.0095 0.0047 0.0015 5.5000e-04 2.8000e-04];
Pb_rrc_z_r2 = [0.0544 0.0336 0.0182 0.0083 0.0032 0.0012 2.7000e-04 1.5000e-04];
Pb_rrc_z_r3 = [0.0554 0.0319 0.0169 0.0071 0.0023 6.2000e-04 7.0000e-05 0];



Pb_fir_B2 = [0.0993 0.0795 0.0599 0.0441 0.0305 0.0206 0.0135 0.0075];
Pb_fir_B4 = [0.1008 0.0792 0.0596 0.0450 0.0310 0.0213 0.0137 0.0075];
Pb_fir_B8 = [0.0957 0.0743 0.0548 0.0382 0.0265 0.0170 0.0100 0.0054];

Pb_fir_z_B2 = [0.0649 0.0433 0.0259 0.0145 0.0071 0.0033 0.0012 3.7000e-04];
Pb_fir_z_B4 = [0.0671 0.0452 0.0270 0.0142 0.0075 0.0031 0.0012 5.7000e-04];
Pb_fir_z_B8 = [0.0612 0.0395 0.0223 0.0118 0.0055 0.0023 8.5000e-04 2.1000e-04];


semilogy(Ebpn,Pb_fir_B2,Ebpn,Pb_fir_B4,Ebpn,Pb_fir_B8,'LineWidth',1.5);
semilogy(Ebpn,Pb_rrc_r1,Ebpn,Pb_rrc_r2,Ebpn,Pb_rrc_r3,'LineWidth',1.5);
legend('Teoretska','Beta=2','Beta=4','Beta=8','r=0.25','r=0.5','r=1','Location','southwest');


figure, semilogy(Ebpn,Pb_teor,'LineWidth',1.5);
grid on
xlabel('E_b/p_N [dB]');
ylabel('BER');
hold on
semilogy(Ebpn,Pb_fir_z_B2,Ebpn,Pb_fir_z_B4,Ebpn,Pb_fir_z_B8,'LineWidth',1.5);
semilogy(Ebpn,Pb_rrc_z_r1,Ebpn,Pb_rrc_z_r2,Ebpn,Pb_rrc_z_r3,'LineWidth',1.5);
legend('Teoretska','Beta=2','Beta=4','Beta=8','r=0.25','r=0.5','r=1','Location','southwest');