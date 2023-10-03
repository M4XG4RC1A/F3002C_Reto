
%% getmodes
clc;

Xs = round(linspace(650,1000,10)*1000);
Ys = round(linspace(550,1000,10)*1000);
Ls = ["\lambda", "A_{eff}"];

lenX = length(Xs);
lenY = length(Ys);
lenL = length(Ls);

%Get min 3 mode
minArea = 100;
mins = [100,100,100];
xMin = [0,0,0,0];
yMin = [0,0,0,0];

%Get 
minVert = [100, 0, 0];
minHor = [100, 0, 0];
sminVert = [100, 0, 0];
sminHor = [100, 0, 0];
minVertArea = [100, 0, 0];
minHorArea = [100, 0, 0];
for y = Ys
        for x=Xs
            load("./nModes/Waveguide"+string(x)+"_"+string(y)+"_1580.mat");
            if nmodes >= 3
                for i =1:3
                    if mins(i)>Aeff(i)*1e12
                        mins(i) = Aeff(i)*1e12;
                        xMin(i) = x/1000;
                        yMin(i) = y/1000;
                    end
                end
                if Aeff(4)~=0
                    if mins(3)>Aeff(4)*1e12
                        mins(3) = Aeff(4)*1e12;
                        xMin(3) = x/1000;
                        yMin(3) = y/1000;
                    end
                end
                if x>y
                    if minHor(1)>Aeff(1)*1e12
                        minHor(1) = Aeff(1)*1e12;
                        minHor(2) = x/1000;
                        minHor(3) = y/1000;
                    end
                    if minHorArea(1)>x*y/(1000000*1000000)
                        minHorArea(1) = x*y/(1000000*1000000);
                        minHorArea(2) = x/1000;
                        minHorArea(3) = y/1000;
                    end
                    if sminHor(1)>Aeff(i)*1e12
                        sminHor(1) = Aeff(i)*1e12;
                        sminHor(2) = x/1000;
                        sminHor(3) = y/1000;
                    end
                elseif y>x
                    if minVert(1)>Aeff(2)*1e12
                        minVert(1) = Aeff(2)*1e12;
                        minVert(2) = x/1000;
                        minVert(3) = y/1000;
                    end
                    if minVertArea(1)>x*y/(1000000*1000000)
                        minVertArea(1) = x*y/(1000000*1000000);
                        minVertArea(2) = x/1000;
                        minVertArea(3) = y/1000;
                    end
                    if sminVert(1)>Aeff(i)*1e12
                        sminVert(1) = Aeff(i)*1e12;
                        sminVert(2) = x/1000;
                        sminVert(3) = y/1000;
                    end
                end
                if minArea>x*y/(1000000*1000000)
                    minArea=x*y/(1000000*1000000);
                    xMin(4) = x/1000;
                    yMin(4) = y/1000;
                end
            end
        end
end

%Min for modes
disp("Min effective areas for modes");
for i=1:3
    disp("Min Effective Area Mode"+num2str(i)+": "+num2str(mins(i)));
    disp("     X: "+num2str(xMin(i)));
    disp("     Y: "+num2str(yMin(i)));
end
disp("------------------------");

%Min for area
disp("Min areas");
disp("Min Area: "+num2str(minArea));
disp("     X: "+num2str(xMin(4)));
disp("     Y: "+num2str(yMin(4)));
disp("------------------------");

%Min for Vertical
disp("Min Vertical");
disp("Min Effective Area: "+num2str(sminVert(1)));
disp("     X: "+num2str(sminVert(2)));
disp("     Y: "+num2str(sminVert(3)));
disp("Min Effective Area (TE): "+num2str(minVert(1)));
disp("     X: "+num2str(minVert(2)));
disp("     Y: "+num2str(minVert(3)));
disp("Min Area: "+num2str(minVertArea(1)));
disp("     X: "+num2str(minVertArea(2)));
disp("     Y: "+num2str(minVertArea(3)));
disp("------------------------");

%Min for horizontal
disp("Min Horizontal");
disp("Min Effective Area: "+num2str(sminHor(1)));
disp("     X: "+num2str(sminHor(2)));
disp("     Y: "+num2str(sminHor(3)));
disp("Min Effective Area (TE): "+num2str(minHor(1)));
disp("     X: "+num2str(minHor(2)));
disp("     Y: "+num2str(minHor(3)));
disp("Min Area: "+num2str(minHorArea(1)));
disp("     X: "+num2str(minHorArea(2)));
disp("     Y: "+num2str(minHorArea(3)));
disp("------------------------");