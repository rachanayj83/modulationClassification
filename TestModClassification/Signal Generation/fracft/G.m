gray = [0; 1];
B=2;
for n = 1:B-1
    
    gray = [[zeros(2^n,1); ones(2^n,1)] , [gray; flipud(gray)]];
    
end


BVec = 0:B-1;

Out=2.^BVec;


A=gray.*Out;
