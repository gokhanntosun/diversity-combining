clear;close all;
snr_db=0:2:20;%snr values in db
N02=1./(10.^(snr_db./10));%snr values
sigmas=sqrt(N02);%std dev of awgn
fading=sqrt(1);%fading std dev

antennas=[1 2 3 4];%antenna count
S=length(snr_db);%constants
A=length(antennas);

snr_ber=zeros(1,length(snr_db));%snr vs ber vector
snr_ser=zeros(1,length(snr_db));%snr vs fer vector
fr_lim=3000;fer_lim=10000;%MC limits

M=4;%modulation alphabet size
m=1:M;%symbol index

spf=5000;%symbols per frame
bps=log2(M);%bits per symbol
bpf=bps*spf;%bits per frame

bb=de2bi((m-1)')';%bit book
sym=(exp(1j*2*pi*(m-1)/M))';%complex symbol book

for a=1:A
    for s=1:S
        %variance of awgn noise
        sigma=sigmas(s);
        fr=0;fr_ber=0;fr_ser=0;
        while fr<fr_lim && fr_ser<fer_lim
            %transmitter
            bits=randi([0 1],[1 bpf]);%generate information sequence
            bits=reshape(bits,[bps spf]);%divide sequence into log2(M)-bits
            sym_ind=(bi2de(bits')+1)';%compute symbol indexes
            x=zeros(1,spf);%transmit signal
            for i=1:spf
                %map to symbols
                x(1,i)=sym(sym_ind(:,i),:);
            end
            x=repmat(x,[a,1]);%repeat for all antennas
            
            %channel
            %complex gaussian noise
            n_ip=normrnd(0,sigma,[a,spf]);
            n_q=normrnd(0,sigma,[a,spf]);
            n=(n_ip+1j*n_q)./sqrt(2);
            %rayleigh fading coeffs
            h_ip=normrnd(0,fading,[a,spf]);
            h_q=normrnd(0,fading,[a,spf]);
            h=(h_ip+1j*h_q)./sqrt(2);
            %impose channel
            r=h.*x+n;
            
            %receiver
            %combining
            r_comb=sum(r.*conj(h),1);
            %detection
            r_comb=repmat(r_comb,[M 1]);%stack received signal, divide by channel impulse
            distance=abs(sym-r_comb);%calculate distance between received symbols and mod. symbols
            [~,det_sym_ind]=min(distance,[],1);%minimum distance for mapping symbols back
            %errors
            sym_err=nnz(sym_ind-det_sym_ind);%symbol errors
            bit_err=nnz(bits-bb(:,det_sym_ind));%calculate bit error
            
            fr_ber=fr_ber+bit_err;%update bit errors
            fr_ser=fr_ser+sym_err;%update symbol errors
            fr=fr+1;%increase frame count
        end%end mc loop
        [a snr_db(s) fr fr_ber]
        snr_ber(a,s)=fr_ber/(fr*bpf);
        snr_ser(a,s)=fr_ser/(fr*spf);
    end%end snr loop
end%end antenna loop
MRC_SER=snr_ser;
save("MRC_SER","MRC_SER");