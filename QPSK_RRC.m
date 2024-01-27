clc, clear all, close all;

%     
%% 1. Generisanje informacione sekvence - diskretan izvor bez memorije     
Nsim = 100000;               % broj generisanih bita
x = rand(Nsim,1);           % generisanje informacione binarne sekvence
P = 0.5;                    % verovatnoca pojave jedinice

%davanje vrednosti binarnoj sekvenci
for i=1:Nsim
    if x(i) > P
        x(i) = 1;
    else
        x(i) = 0;
        
    end
end


%% blok za mapiranje QPSK
x = x';
QPSK_signal = [];
for n_brojac_1 = 1:2:length(x)
    if x(1,n_brojac_1)== 0 && x(1,n_brojac_1+1)== 0
        x1 = exp(1i*pi/4);
    elseif x(1,n_brojac_1)== 0 && x(1,n_brojac_1+1)== 1
        x1 = exp(1i*3*pi/4);
    elseif x(1,n_brojac_1)== 1 && x(1,n_brojac_1+1)== 1
        x1 = exp(1i*5*pi/4);
    elseif x(1,n_brojac_1)== 1 && x(1,n_brojac_1+1)== 0
        x1 = exp(1i*7*pi/4);
    end
    QPSK_signal = [QPSK_signal x1];
end




%% 2. obrada signala

V_sim = 2000;                   %protok po simbolu [sim/s]
Nbps = 2;                       %broj bita po simbolu
V_bit = V_sim * Nbps;           %bitski protok [b/s]
Ts = 1 / V_sim;                 %trajanje simbola
M1up = 32;                      %faktor uvecanja ucestanosti                      
Nsps = M1up;                    %broj odbiraka po simbolu
Tb = Ts / Nbps;                 %trajanje bita
fb = 1/Tb;
Todb = Ts/Nsps;                 %trajanje odbirka
fs = 1/Todb;                    %ucestanost odabiranja


r1 = 0.25;                      %roll-off faktor RRC filtra
r2 = 0.5; 
r3 = 1;

% podela signala na I i Q granu
I_grana = real(QPSK_signal);
Q_grana = imag(QPSK_signal);

% povecanje frekvencije odabiranja
I_grana_up = upsample(I_grana,M1up);
Q_grana_up = upsample(Q_grana,M1up);



%% uoblicavanje filtrom
N_rrc = 100;         %red rcc filtra
h_rrc = srrcf(N_rrc, M1up, r1);   % odredjivanje impulsnog odziva filtra
h_rrc2 = srrcf(N_rrc, M1up, r2);
h_rrc3 = srrcf(N_rrc, M1up, r3);

% crtanje impulsnog odziva filtara
osa_rrc = (-N_rrc/2:N_rrc/2);
figure, stem(osa_rrc, h_rrc);
hold on
stem(osa_rrc, h_rrc2);
stem(osa_rrc, h_rrc3);
grid on
xlabel('n');
ylabel('h(n)');
title('Impulsni odziv RRC filtra');
legend('r = 0.25','r = 0.5','r = 1');

% prelazak u spektralni domen
[H_rrc, wosa] = freqz(h_rrc, 1, 1024);
[H_rrc2, wosa] = freqz(h_rrc2, 1, 1024);
[H_rrc3, wosa] = freqz(h_rrc3, 1, 1024);

% crtanje amplitudske karakteristike filtara
figure, plot(wosa/pi, [abs(H_rrc) abs(H_rrc2) abs(H_rrc3)])
ylabel('H(w)');
xlabel('w/pi'), title('Amplitudska karakteristika RRC filtra');
legend('r = 0.25','r = 0.5','r = 1');
grid on
% crtanje pojacanja filtara
figure, plot(wosa/pi, [20*log10(abs(H_rrc)),20*log10(abs(H_rrc2)),20*log10(abs(H_rrc3))]);
xlabel('w/pi'), ylabel('g[dB]'),title('Pojacanje RRC');
legend('r = 0.25','r = 0.5','r = 1');
grid on
% crtanje fazne karakteristike filtara
figure, plot(wosa/pi, [unwrap(angle(H_rrc)),unwrap(angle(H_rrc2)),unwrap(angle(H_rrc3))]);
ylabel('{\phi}')
xlabel('w/pi'), title('Fazna karakteristika RRC');
legend('r = 0.25','r = 0.5','r = 1');
grid on


% filtriranje signala
[I_grana_rrc, zf1] = filter(h_rrc, 1, I_grana_up);
[Q_grana_rrc, zf2] = filter(h_rrc, 1, Q_grana_up);

% poravnavanje signala usled kasnjenja
duzs = length(I_grana_rrc);
duzr = length(zf1);
I_grana_rrc(duzs+1:duzs+duzr) = zf1;
I_grana_rrc = I_grana_rrc * M1up;
Q_grana_rrc(duzs+1:duzs+duzr) = zf2;
Q_grana_rrc = Q_grana_rrc * M1up;


% spektar I grane
 [P0_I, w] = pwelch(I_grana_rrc,bartlett(M1up),0,512,fs);
 figure,plot(w/max(w),10*log10(P0_I),'b-o');
 hold on

% spektar Q grane
[P0_Q, w] = pwelch(Q_grana_rrc,bartlett(M1up),0,512,fs);
plot(w/max(w),10*log10(P0_Q),'r*')
grid on
xlabel('Normalizovana frekvencija');
ylabel('Snaga/Frekvencija [dB/rad]');
legend('I grana','Q grana');

% kreiranje predajnog signala
Predaja = I_grana_rrc + 1i*Q_grana_rrc;


%% ABGS - kanal
Ebpn = 15;           %odnos srednje bitske energije i SGSS suma u dB
SNRdB = Ebpn + 10*log10(Nbps) - 10*log10(Nsps);   % definisanje SNR

% dodavanje suma u signal
y = awgn(Predaja,SNRdB,'measured');


%% prijemni deo
% podela signala u Q i I granu na prijemu
I_grana_prijem = real(y);
Q_grana_prijem = imag(y);

% filtriranje signala na prijemu
[I_grana_prijem_rrc, zf1 ] = filter(h_rrc, 1, I_grana_prijem);
[Q_grana_prijem_rrc, zf2 ] = filter(h_rrc, 1, Q_grana_prijem);

% poravnavanje signala na prijemu
duzs = length(I_grana_prijem_rrc);
duzr = length(zf1);
I_grana_prijem_rrc(duzs+1:duzs+duzr) = zf1;
Q_grana_prijem_rrc(duzs+1:duzs+duzr) = zf2;

% smanjenje frekvencije
I_grana_prijem_down = downsample(I_grana_prijem_rrc,M1up);
Q_grana_prijem_down = downsample(Q_grana_prijem_rrc,M1up);
delay = N_rrc/(M1up*2);
I_grana_prijem_down = I_grana_prijem_down(2*delay+1:end-2*delay);
Q_grana_prijem_down = Q_grana_prijem_down(2*delay+1:end-2*delay);

% formiranje prijemnog signala
Prijemni_QPSK = I_grana_prijem_down + 1i*Q_grana_prijem_down;

% crtanje konstelacionog dijagrama
figure, plot(real(Prijemni_QPSK),imag(Prijemni_QPSK),'x');
axis([-1.5 1.5 -1.5 1.5])
grid on
xlabel('In-Phase');
ylabel('Quadrature');
hold on
plot(real(QPSK_signal),imag(QPSK_signal),'o','MarkerSize',10,'MarkerFaceColor','r');


%% blok za demapiranje QPSK
biti_prijem = [];

for n_brojac_2 = 1:length(Prijemni_QPSK)
    if real(Prijemni_QPSK(n_brojac_2)) >= 0 && imag(Prijemni_QPSK(n_brojac_2)) >= 0
        biti_prijem_pom = [0 0];
    elseif real(Prijemni_QPSK(n_brojac_2)) >= 0 && imag(Prijemni_QPSK(n_brojac_2)) < 0
        biti_prijem_pom = [1 0];
    elseif real(Prijemni_QPSK(n_brojac_2)) < 0 && imag(Prijemni_QPSK(n_brojac_2)) < 0
        biti_prijem_pom = [1 1];
    elseif real(Prijemni_QPSK(n_brojac_2)) < 0 && imag(Prijemni_QPSK(n_brojac_2)) >= 0
        biti_prijem_pom = [0 1];
    end
    biti_prijem = [biti_prijem biti_prijem_pom];
end



[Broj_gresaka,BER] = biterr(x,biti_prijem)


