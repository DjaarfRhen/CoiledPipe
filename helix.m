function [] = helix()

    helixCurvTors(0.2, 0.4);
%    hold on;
%     helixCurvTors(0.01, 0.1);
%     hold on;
%     helixCurvTors(0.01, 0.13);
%     hold on;
%     helixCurvTors(0.01, 1.0);
%     hold on;
%    helixRadPitch(5,100/(2*pi));
%    hold off;
    
end

function [ ] = helixCurvTors( Curv, Tors )

    step = 1e-2;
    s = 0:step:200;
    
    Omega = sqrt(Curv^2 + Tors^2);
    Rad = Curv/Omega^2;
    Pitch = Tors/Omega^2;
    
    x = Rad*cos(Omega*s);
    y = Rad*sin(Omega*s);
    z = Pitch*Omega*s;
        
    length = s(end)*sqrt(Rad^2+Pitch^2);
    
    fprintf('Rad = %f\nPitch = %f\nstep = %f\nlength = %f\n\n', Rad, Pitch, 2*pi*Pitch, length);
    
    plot3(x,y,z,'green');
    axis equal;
    
end

function [ ] = helixRadPitch( Rad, Pitch )

    step = 1e-2;
    s = 0:step:10;
    
    x = Rad*cos(s);
    y = Rad*sin(s);
    z = Pitch*s;
    
    Curv = abs(Rad)/(Rad^2+Pitch^2);
    Tors = Pitch/(Rad^2+Pitch^2);
    length = s(end)*sqrt(Rad^2+Pitch^2);
    
    fprintf('curv = %f\ntors = %f\nlength = %f\n\n', Curv, Tors, length);
    
    plot3(x,y,z,'red');
    axis equal;
    
end

