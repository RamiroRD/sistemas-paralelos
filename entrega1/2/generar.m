pkg load image

arg_list = argv();
N = str2num(arg_list{1});
file = fopen(sprintf('data-%d.bin', N),'w');
A = 1 * rand(N, N);
B = 2 * rand(N, N);
C = 3 * rand(N, N);
D = 4 * rand(N, N);
E = 5 * rand(N, N);
F = 6 * rand(N, N);
L = 7 * tril(rand(N,N));
U = 8 * triu(rand(N,N));


fwrite(file, A', 'double', 'l');
fwrite(file, B', 'double', 'l');
fwrite(file, C', 'double', 'l');
fwrite(file, D', 'double', 'l');
fwrite(file, E', 'double', 'l');
fwrite(file, F', 'double', 'l');
fwrite(file, L', 'double', 'l');
fwrite(file, U', 'double', 'l');


res = mean2(U) * mean2(L) * A * A * C + mean2(B) * ( (L * B * E) + (D * U * F) );

fwrite(file, res', 'double', 'l');
fclose(file);
