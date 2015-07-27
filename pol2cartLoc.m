function [vx,vy] = pol2cartLoc(th,r,TH)

vx = r.*cos(TH)-th.*sin(TH);
vy = r.*sin(TH)+th.*cos(TH);

