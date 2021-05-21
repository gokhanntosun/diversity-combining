clear
load EGC_SER.mat
load MRC_SER.mat
load SC_SER.mat
load A21_SER.mat
load A22_SER.mat
load MIMO22_SER.mat

SNR_DB=0:2:20;

%grouped by combining method
figure;%selection combining
semilogy(SNR_DB,SC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,SC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,SC_SER(3,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,SC_SER(4,:),"s-","LineWidth",2);hold on;
legend("M=1","M=2","M=3","M=4","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("Selection Combining");
ylim([10^-5 1]);grid on;axis square;

figure;%equal gain combining
semilogy(SNR_DB,EGC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,EGC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,EGC_SER(3,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,EGC_SER(4,:),"s-","LineWidth",2);hold on;
legend("M=1","M=2","M=3","M=4","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("Equal Gain Combining");
ylim([10^-5 1]);grid on;axis square;

figure;%maximal-ratio combining
semilogy(SNR_DB,MRC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MRC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MRC_SER(3,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MRC_SER(4,:),"s-","LineWidth",2);hold on;
legend("M=1","M=2","M=3","M=4","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("Maximal-Ratio Combining");
ylim([10^-5 1]);grid on;axis square;

%grouped by receive antenna count

figure;%two antennas
semilogy(SNR_DB,SC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,SC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,EGC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MRC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,A22_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MIMO22_SER(1,:),"s-","LineWidth",2);hold on;
legend("No Diversity","SC","EGC","MRC",...
    "2x2 Alamouti","2x2 MIMO","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("2 Antenna @ RX");
ylim([10^-5 1]);grid on;axis square;

figure;%three antennas
semilogy(SNR_DB,SC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,SC_SER(3,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,EGC_SER(3,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MRC_SER(3,:),"s-","LineWidth",2);hold on;
legend("No Diversity","SC","EGC","MRC","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("3 Antenna @ RX");
ylim([10^-5 1]);grid on;axis square;

figure;%four antennas
semilogy(SNR_DB,SC_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,SC_SER(4,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,EGC_SER(4,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,MRC_SER(4,:),"s-","LineWidth",2);hold on;
legend("No Diversity","SC","EGC","MRC","Location","Southwest");
xlabel("SNR(dB)");ylabel("SER");title("4 Antenna @ RX");
ylim([10^-5 1]);grid on;axis square;

%alamouti plots
%a21 vs. mrc
figure;
semilogy(SNR_DB,MRC_SER(2,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,A21_SER(1,:),"s-","LineWidth",2);hold on;
legend("MRC | M=2","2x1 Alamouti");
xlabel("SNR(dB)");ylabel("SER");title("MRC with 2 Antennas vs. 2x1 Alamouti");
ylim([10^-5 1]);grid on;axis square;

%a21 vs. a22
figure;
semilogy(SNR_DB,A21_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,A22_SER(1,:),"s-","LineWidth",2);hold on;
legend("2x1 Alamouti","2x2 Alamouti");
xlabel("SNR(dB)");ylabel("SER");title("2x1 Alamouti vs. 2x2 Alamouti");
ylim([10^-5 1]);grid on;axis square;

%mimo22 vs. a22
figure;
semilogy(SNR_DB,MIMO22_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,A21_SER(1,:),"s-","LineWidth",2);hold on;
semilogy(SNR_DB,A22_SER(1,:),"s-","LineWidth",2);hold on;
legend("2x2 MIMO","2x1 Alamouti","2x2 Alamouti");
xlabel("SNR(dB)");ylabel("SER");title("MIMO vs. Alamouti");
ylim([10^-5 1]);grid on;axis square;





