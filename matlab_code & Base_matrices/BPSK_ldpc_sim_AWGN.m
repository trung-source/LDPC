
clear
close all
format short g

%EbNodB = 2;
MaxItrs = 8;

%Basegraphnumber_j_z
load base_matrices/NR_1_0_16.txt
B = NR_1_0_16;
[mb,nb] = size(B);
z = 16;

Slen = sum(B(:)~=-1); %So luong cac gia tri khac -1 trong B
min_reg = zeros(max(sum(B ~= -1,2)),z); %Thanh ghi luu tru cho cac gia tri minsum
k = (nb-mb)*z; %So cot bit message trong B
n = nb*z; %So cot bit codeword trong B
Rate = k/n;  %Toc do


dB=[0.5:0.25:3];                % Khoang SNR chay theo dB
FER=zeros(1,length(dB));        % Mang luu tru Frame Error Rate
BER=zeros(1,length(dB));        % Mang luu tru Bit Error Rate

FER_uncoded=zeros(1,length(dB));        % Mang luu tru Frame Error Rate khong ma hoa
BER_uncoded=zeros(1,length(dB));        % Mang luu tru Bit Error Rate khong ma hoa

FER_iter_2=zeros(1,length(dB));        % Mang luu tru Frame Error Rate voi lan lap 2
BER_iter_2=zeros(1,length(dB));        % Mang luu tru Bit Error Rate voi lan lap 2

FER_iter_4=zeros(1,length(dB));        % Mang luu tru Frame Error Rate voi lan lap 4
BER_iter_4=zeros(1,length(dB));        % Mang luu tru Bit Error Rate voi lan lap 4

EbNo=10.^(dB/10);                % Eb/No doi tu dB thanh thap phan



Nblocks = 100;                      % So luong Block Message mo phong
disp('          dB          FER        BER           Nblkerrs    Nbiterrs     Nblocks')
for g=1:length(EbNo) % Vong lap de test trong khoang SNR
    Nbiterrs = 0; Nblkerrs = 0; 
    Nbiterrs_uncoded = 0; Nblkerrs_uncoded = 0; 
    Nbiterrs_iter_2 = 0; Nblkerrs_iter_2 = 0; 
    Nbiterrs_iter_4 = 0; Nblkerrs_iter_4 = 0; 
    sigma = sqrt(1/(2*Rate*EbNo(g)));
    for i = 1: Nblocks
        msg = randi([0 1],1,k); %Tao k-bit message ngau nhien
        %msg = zeros(1,k); %all-zero message
        %cword = zeros(1,n); %all-zero codeword
        
        %Encoding 
        cword = ldpc_encode(B,z,msg);
        cword = cword(1:n);  
        s = 1 - 2 * cword; %Ma hoa BPSK
        %h = 1/sqrt(2)*[randn(1,length(s)) + j*randn(1,length(s))];
        r = s + sigma * randn(1,n); %AWGN channel I
        
        %Soft-decision, Layer decoding
        L = r; %Cac belief nhan duoc
        itr = 0; %so lan lap hien tai
        R = zeros(Slen,z); %Thanh ghi chua du lieu xu ly tren hang
        
        % Thuc hien Layer Decoding
        while itr < MaxItrs
            Ri = 0;
            for lyr = 1:mb %lyr la layer, mac dinh la 1
                ti = 0; %so luong cac gia tri khac -1 tai row=lyr
                for col = find(B(lyr,:) ~= -1)
                       ti = ti + 1;
                       Ri = Ri + 1;
                       %Subtraction (Thuc hien phep tru)
                       L((col-1)*z+1:col*z) = L((col-1)*z+1:col*z)-R(Ri,:);
                       %can chinh lai gia tri hang va luu tru trong min_reg
                       min_reg(ti,:) = mul_sh(L((col-1)*z+1:col*z),B(lyr,col)); 
                end
                %minsum tai min_reg: ti x z
                for i1 = 1:z %min_reg(1:ti,i1)
                    [min1,pos] = min(abs(min_reg(1:ti,i1))); %Cuc tieu thu nhat theo gia tri tuyet doi
                    min2 = min(abs(min_reg([1:pos-1 pos+1:ti],i1))); %Cuc tieu thu 2 theo gia tri tuyet doi
                    S = sign(min_reg(1:ti,i1)); %Luu tru dau cua 1 hang
                    parity = prod(S); %Tich cac dau cua 1 hang
                    min_reg(1:ti,i1) = min1; %Thay the thanh gia tri tuyet doi min1
                    min_reg(pos,i1) = min2; %Thay the cho vi tri cua min1 thanh gia tri tuyet doi min2
                    min_reg(1:ti,i1) = parity*S.*min_reg(1:ti,i1); %Nhan them dau cua hang
                end
                %Thuc hien can chinh lai gia tri cot, addition va luu tru trong R
                Ri = Ri - ti; %reset Ri
                ti = 0;
                for col = find(B(lyr,:) ~= -1)
                        Ri = Ri + 1;
                        ti = ti + 1;
                        %Can chinh gia tri cot
                        R(Ri,:) = mul_sh(min_reg(ti,:),z-B(lyr,col));
                        %Addition (Thuc hien phep cong)
                        L((col-1)*z+1:col*z) = L((col-1)*z+1:col*z)+R(Ri,:);
                end
            end      
            if itr == 2
               msg_iter_2 = L(1:k)<0;
            end
            if itr == 4
               msg_iter_4 = L(1:k)<0;
            end       
            itr = itr + 1;       
        end
        msg_cap = L(1:k) < 0; %Quyet dinh
        uncoded_msg = r(1:k) < 0;
        %Dem loi
        Nerrs_uncoded = sum(msg ~= uncoded_msg);
        Nerrs_iter_2 = sum(msg ~= msg_iter_2);
        Nerrs_iter_4 = sum(msg ~= msg_iter_4);
        Nerrs = sum(msg ~= msg_cap);
        if Nerrs > 0
            Nbiterrs = Nbiterrs + Nerrs;
            Nblkerrs = Nblkerrs + 1;
        end
        if Nerrs_uncoded > 0
            Nbiterrs_uncoded = Nbiterrs_uncoded + Nerrs_uncoded;
            Nblkerrs_uncoded = Nblkerrs_uncoded + 1;
        end
        if Nerrs_iter_2 > 0
            Nbiterrs_iter_2 = Nbiterrs_iter_2 + Nerrs_iter_2;
            Nblkerrs_iter_2 = Nblkerrs_iter_2 + 1;
        end
        if Nerrs_iter_4 > 0
            Nbiterrs_iter_4 = Nbiterrs_iter_4 + Nerrs_iter_4;
            Nblkerrs_iter_4 = Nblkerrs_iter_4 + 1;
        end  
    end
   
    BER(g) = Nbiterrs/k/Nblocks;
    FER(g) = Nblkerrs/Nblocks;
    
    BER_uncoded(g) = Nbiterrs_uncoded/k/Nblocks;
    FER_uncoded(g) = Nblkerrs_uncoded/Nblocks;
    
    BER_iter_2(g) = Nbiterrs_iter_2/k/Nblocks;
    FER_iter_2(g) = Nblkerrs_iter_2/Nblocks;
    
    BER_iter_4(g) = Nbiterrs_iter_4/k/Nblocks;
    FER_iter_4(g) = Nblkerrs_iter_4/Nblocks;

    %disp([EbNodB FER_sim BER_sim Nblkerrs Nbiterrs Nblocks])
    disp([dB(g) FER(g) BER(g) Nblkerrs Nbiterrs Nblocks])
end
	
subplot(2,2,1)
plot(dB,BER,'b-o')
legend('BER MaxIter ');
title('Bit Error Rate')
ylabel('BER')
xlabel('Eb/No (dB)')

subplot(2,2,2)
plot(dB,BER,'b-o')
hold on;
plot(dB,BER_uncoded,'r--o')
plot(dB,BER_iter_2,'k:o')
plot(dB,BER_iter_4,'m-.o')
legend('BER MaxIter ','BER Uncoded','BER iter 2','BER iter 4');
title('Bit Error Rate')
ylabel('BER')
xlabel('Eb/No (dB)')
grid
%figure
hold off;

subplot(2,2,3)
plot(dB,FER,'b-o')
legend('FER MaxIter ');
title('Frame Error Rate')
ylabel('FER')
xlabel('Eb/No (dB)')

subplot(2,2,4)
plot(dB,FER,'b-o')
hold on;
plot(dB,FER_uncoded,'r--o')
plot(dB,FER_iter_2,'k:o')
plot(dB,FER_iter_4,'m-.o')
legend('FER MaxIter ','FER Uncoded','FER iter 2','FER iter 4');
title('Frame Error Rate')
ylabel('FER')
xlabel('Eb/No (dB)')
grid





