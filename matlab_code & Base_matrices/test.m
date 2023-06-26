

clear;
clc;
%Basegraphnumber_j_z
load base_matrices/NR_1_0_16.txt
B = NR_1_0_16;
[mb,nb] = size(B);
z = 16;

k = (nb-mb)*z; %So cot bit message trong B
n = nb*z; %So cot bit codeword trong B
Rate = k/n;  %Toc do

msg = randi([0 1],1,k);
msg = zeros(1,k); %all-zero message
cword = ldpc_encode(B,z,msg);
c = cword;
cword = cword(1:n);


EbNo=2; % Eb/No theo db
sigma = sqrt(1/(2*Rate*EbNo));
s = 1 - 2 * cword; %Ma hoa BPSK
%h = 1/sqrt(2)*[randn(1,length(s)) + j*randn(1,length(s))];
r = s + sigma * randn(1,n); %AWGN channel I


