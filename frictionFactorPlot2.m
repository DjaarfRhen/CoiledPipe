K = find(-diff(t)>500);
tmp = Dp_raw([K-1; end]);
Dp = reshape(tmp, [11 51]);
imagesc(1:0.1:2, 0:0.1:5, Dp);
set(gca,'YDir','normal')
axis image


K = find(-diff(VarName1)>500);
Dp = VarName3([K-1; end]);
loglog(100:100:2000, 32./[100:100:2000], 'green');
hold on;
loglog(100:100:2000, 8*Dp, 'red');