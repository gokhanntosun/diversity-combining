clear;close all;
snr_db=0:2:20;%snr values in db
N02=1./(10.^(snr_db./10));%snr values
sigmas=sqrt(N02);%std dev of awgn
fading=sqrt(1);%fading std dev

A=[2 2];%antennas at [rx tx]
SNR_L=length(snr_db);%constants

snr_ber=zeros(1,length(snr_db));%snr vs ber vector
snr_ser=zeros(1,length(snr_db));%snr vs fer vector
fr_lim=3000;fer_lim=100000;%MC limits

M=4;%modulation alphabet size
m=1:M;%symbol index

spf=5000;%symbols per frame
bps=log2(M);%bits per symbol
bpf=bps*spf;%bits per frame

bb=de2bi((m-1)')';%bit book
sym=((exp(1j*2*pi*(m-1)/M))/sqrt(A(2)))';%complex symbol book
%allocate tx power among two antennas by dividing by sqrt2

for s=1:SNR_L
    %variance of awgn noise
    sigma=sigmas(s);
    fr=0;fr_ber=0;fr_ser=0;
    while fr<fr_lim && fr_ser<fer_lim
        %rayleigh fading coeffs
        h_ip=normrnd(0,fading,[A(2),A(1),spf]);
        h_q=normrnd(0,fading,[A(2),A(1),spf]);
        h=((h_ip+1j*h_q)./sqrt(2));
        %complex gaussian noise
        n_ip=normrnd(0,sigma,[A(2),spf]);
        n_q=normrnd(0,sigma,[A(2),spf]);
        n=(n_ip+1j*n_q)./sqrt(2);
        
        %transmitter
        bits=randi([0 1],[1 A(2)*bpf]);%generate information sequence
        bits=reshape(bits,[bps A(2)*spf]);%divide sequence into log2(M) bits
        sym_ind=(bi2de(bits')+1)';%compute symbol indexes
        x_td=zeros(1,A(2)*spf);
        for i=1:A(2)*spf
            %map to symbols, x tilda
            x_td(1,i)=sym(sym_ind(:,i),:);
        end
        %reshape to obtain parallel data streams
        x_td=reshape(x_td,[A(2) spf]);
        y_td=zeros(A(1),spf);
        %transmit
        for i=1:spf
            H=h(:,:,i);
            %singular value decomp. of channel
            [U,S,V]=svd(H);
            %transmit precoding
            x=V*x_td(:,i);
            %channel
            y=awgn(H*x,snr_db(s));%+n(:,i);
            %receiver
            %receiver shaping
            y_td(:,i)=S\(U'*y);%[y1;y2 y3;y4 ...]
        end
        %receiver cont.
        y_td=reshape(y_td,[1,A(1)*spf]);
        y_td=repmat(y_td,[M 1]);
        [~,det_sym_ind]=min(abs(sym-y_td),[],1);%minimum distance for mapping symbols back
        %errors
        sym_err=nnz(sym_ind-det_sym_ind);%symbol errors
        bit_err=nnz(bits-bb(:,det_sym_ind));%calculate bit error
        
        fr_ber=fr_ber+bit_err;%update bit errors
        fr_ser=fr_ser+sym_err;%update symbol errors
        fr=fr+1;%increase frame count
    end%end mc loop
    [snr_db(s) fr fr_ber]
    snr_ber(1,s)=fr_ber/(fr*bpf);
    snr_ser(1,s)=fr_ser/(fr*spf);
end%end snr loop
MIMO22_SER=snr_ser;
save("MIMO22_SER","MIMO22_SER");