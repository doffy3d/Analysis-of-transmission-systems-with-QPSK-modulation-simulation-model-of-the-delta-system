clc, clear all, close all;

     
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




% podela signala na I i Q granu
I_grana = real(QPSK_signal);
Q_grana = imag(QPSK_signal);

% povecanje frekvencije odabiranja
I_grana_up = upsample(I_grana,M1up);
Q_grana_up = upsample(Q_grana,M1up);



%% uoblicavanje filtrom
N_fir = 100;         %red rcc filtra
h_fir = fir1(N_fir, 1/M1up, kaiser(N_fir+1,2));
h_fir2 = fir1(N_fir, 1/M1up, kaiser(N_fir+1,4));
h_fir3 = fir1(N_fir, 1/M1up, kaiser(N_fir+1,8));

% crtanje impulsnog odziva filtara
osa_rrc = (-N_fir/2:N_fir/2);
figure, stem(osa_rrc, h_fir);
hold on
stem(osa_rrc, h_fir2);
stem(osa_rrc, h_fir3);
grid on
xlabel('n');
ylabel('h(n)');
title('Impulsni odziv');
legend('Beta = 2','Beta = 4','Beta = 8');

% prelazak u spektralni domen
[H_rrc, wosa] = freqz(h_fir, 1, 1024);
[H_rrc2, wosa] = freqz(h_fir2, 1, 1024);
[H_rrc3, wosa] = freqz(h_fir3, 1, 1024);

% crtanje amplitudske karakteristike filtara
figure, plot(wosa/pi, [abs(H_rrc) abs(H_rrc2) abs(H_rrc3)])
ylabel('H(w)');
xlabel('w/pi'), title('Amplitudska karakteristika');
legend('Beta = 2','Beta = 4','Beta = 8');
grid on
% crtanje pojacanja filtara
figure, plot(wosa/pi, [20*log10(abs(H_rrc)),20*log10(abs(H_rrc2)),20*log10(abs(H_rrc3))]);
xlabel('w/pi'), ylabel('g[dB]'),title('Pojacanje');
legend('Beta = 2','Beta = 4','Beta = 8');
grid on
% crtanje fazne karakteristike filtara
figure, plot(wosa/pi, [unwrap(angle(H_rrc)),unwrap(angle(H_rrc2)),unwrap(angle(H_rrc3))]);
ylabel('{\phi}')
xlabel('w/pi'), title('Fazna karakteristika');
legend('Beta = 2','Beta = 4','Beta = 8');
grid on


% filtriranje signala
[I_grana_fir, zf1] = filter(h_fir2, 1, I_grana_up);
[Q_grana_fir, zf2] = filter(h_fir2, 1, Q_grana_up);

% poravnavanje signala usled kasnjenja
duzs = length(I_grana_fir);
duzr = length(zf1);
I_grana_fir(duzs+1:duzs+duzr) = zf1;
I_grana_fir = I_grana_fir * M1up;
Q_grana_fir(duzs+1:duzs+duzr) = zf2;
Q_grana_fir = Q_grana_fir * M1up;

% spektar I grane
 [P0_I, w] = pwelch(I_grana_fir,bartlett(M1up),0,512,fs);
 figure,plot(w/max(w),10*log10(P0_I),'b-o');
 hold on

% spektar Q grane
[P0_Q, w] = pwelch(Q_grana_fir,bartlett(M1up),0,512,fs);
plot(w/max(w),10*log10(P0_Q),'r*')
grid on
xlabel('Normalizovana frekvencija');
ylabel('Snaga/Frekvencija [dB/rad]');
legend('I grana','Q grana');


% kreiranje predajnog signala
Predaja = I_grana_fir + 1i*Q_grana_fir;


%% ABGS - kanal
Ebpn = 15;           %odnos srednje bitske energije i SGSS suma u dB
SNRdB = Ebpn + 10*log10(Nbps) - 10*log10(Nsps);   % definisanje SNR

% dodavanje suma u signal
y = awgn(Predaja,SNRdB,'measured');


%% prijemni deo
% podela signala u Q i I granu na prijemu
I_grana_prijem = real(y);
Q_grana_prijem = imag(y);


% % spektar I grane
% [P0_I, w] = pwelch(I_grana_prijem,bartlett(M1up),0,2048,fs);
% figure,plot(w/max(w),10*log10(P0_I),'b-o')
% hold on
% 
% % spektar Q grane
% [P0_Q, w] = pwelch(Q_grana_prijem,bartlett(M1up),0,512,fs);
% plot(w/max(w),10*log10(P0_Q),'r*')


% filtriranje signala na prijemu
[I_grana_prijem_fir, zf1 ] = filter(h_fir2, 1, I_grana_prijem);
[Q_grana_prijem_fir, zf2 ] = filter(h_fir2, 1, Q_grana_prijem);

% poravnavanje signala na prijemu
duzs = length(I_grana_prijem_fir);
duzr = length(zf1);
I_grana_prijem_fir(duzs+1:duzs+duzr) = zf1;
Q_grana_prijem_fir(duzs+1:duzs+duzr) = zf2;

% smanjenje frekvencije
I_grana_prijem_down = downsample(I_grana_prijem_fir,M1up);
Q_grana_prijem_down = downsample(Q_grana_prijem_fir,M1up);
delay = N_fir/(M1up*2);
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


