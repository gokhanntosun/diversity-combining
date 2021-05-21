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
sym=((exp(1j*2*pi*(m-1)/M))/sqrt(2))';%complex symbol book
%allocate tx power among two antennas by dividing by sqrt2

a=2;
    for s=1:S
        %variance of awgn noise
        sigma=sigmas(s);
        fr=0;fr_ber=0;fr_ser=0;
        while fr<fr_lim && fr_ser<fer_lim
            %transmitter
            bits=randi([0 1],[1 bpf]);%generate information sequence
            bits=reshape(bits,[bps spf]);%divide sequence into log2(M)-bits
            sym_ind=(bi2de(bits')+1)';%compute symbol indexes
            sym_ind=reshape(sym_ind,[a,spf/a]);
            x_norm=zeros(a,spf/a);%transmit signal
            for i=1:spf/a
                %t1
                x_norm(1,i)=sym(sym_ind(1,i),:);%s1
                x_norm(2,i)=sym(sym_ind(2,i),:);%s2
            end
            %t2 symbols
            x_conj=[-1;1].*conj(flipud(x_norm));
            x=zeros(a,spf);x(:,1:2:end)=x_norm;x(:,2:2:end)=x_conj;
            
            %channel
            %complex gaussian noise
            n_ip=normrnd(0,sigma,[a,spf]);
            n_q=normrnd(0,sigma,[a,spf]);
            n=(n_ip+1j*n_q)./sqrt(2);
            %rayleigh fading coeffs
            %2x2 at time t, 3rd dim. represents the time indexes
            h_ip=normrnd(0,fading,[2,2,spf/2]);
            h_q=normrnd(0,fading,[2,2,spf/2]);
            h=(h_ip+1j*h_q)./sqrt(2);
            h=repelem(h,1,1,2);%constant channel for two symbol periods
            r=zeros(a,spf);
            for i=1:spf
                %multiply by fading at each symbol period
                r(:,i)=h(:,:,i)*x(:,i);
            end
            r=r+n;
            
            %receiver
            r=reshape(r,[a+a spf/a]);%[y1; y2; y3; y4]
            h=h(:,:,1:2:end);%remove repeating indexes
            z=zeros(a,spf/a);
            for i=1:spf/a
                z(1,i)=conj(h(1,1,i))*r(1,i)+conj(h(2,1,i))*r(2,i)+...
                    h(1,2,i)*conj(r(3,i))+h(2,2,i)*conj(r(4,i));
                z(2,i)=conj(h(1,2,i))*r(1,i)+conj(h(2,2,i))*r(2,i)-...
                    h(1,1,i)*conj(r(3,i))-h(2,1,i)*conj(r(4,i));
            end
            z=reshape(z,[1 spf]);%return back to vector 
            
            %detection
            z=repmat(z,[M 1]);%stack the received signal
            [~,det_sym_ind]=min(abs(sym-z),[],1);%minimum distance for mapping symbols back
            det_sym_ind=reshape(det_sym_ind,[a spf/a]);
            %errors
            sym_err=nnz(sym_ind-det_sym_ind);%symbol errors
            bit_err=nnz(bits-bb(:,det_sym_ind));%calculate bit error
            
            fr_ber=fr_ber+bit_err;%update bit errors
            fr_ser=fr_ser+sym_err;%update symbol errors
            fr=fr+1;%increase frame count
        end%end mc loop
        [a snr_db(s) fr fr_ber]
        snr_ber(1,s)=fr_ber/(fr*bpf);
        snr_ser(1,s)=fr_ser/(fr*spf);
    end%end snr loop
A22_SER=snr_ser;
save("A22_SER","A22_SER");