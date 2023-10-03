
%% getmodes
clc;

Xs = linspace(300,1000,5);
Ys = linspace(100,1000,5);
Ls = [510,750,1000,1530,1555,1580];

lenX = length(Xs);
lenY = length(Ys);
lenL = length(Ls);

fileID = fopen("Table.txt",'w');

%Start table
fprintf(fileID,'%s\n',"\begin{table}[H]\centering");
fprintf(fileID,'%s\n',"\resizebox{0.45\textwidth}{!}{");
fprintf(fileID,'%s\n',"  \begin{tabular}{rr|"+convertCharsToStrings(repmat('c',1,lenX))+"|ll}");
fprintf(fileID,'%s\n',"     & &\multicolumn{"+num2str(lenX)+"}{c}{Wg_x}\\");

strLine = "     & ";
for x = Xs
    strLine = strLine+"&"+num2str(x)+"nm";
end
strLine = strLine+"\\ \cline{1-"+num2str(lenX+4)+"} \\[-1em]";
fprintf(fileID,'%s\n',strLine);

firstStrStart = "     \multirow{"+num2str(lenL*lenY)+"}{*}{Wg_y}";
endStrLine = "&\multirow{"+num2str(lenL*lenY)+"}{*}{$\lambda$}";

for y = Ys
    strLine = firstStrStart+"&\multirow{"+num2str(lenL)+"}{*}{"+num2str(y)+"nm}";
    for l = Ls
        for x=Xs
            load("./nModes/Waveguide"+string(x)+"_"+string(y)+"_"+string(l)+".mat");
            strLine = strLine+"&"+num2str(nmodes);
        end
        strLine = strLine+endStrLine+"\\";
        fprintf(fileID,'%s\n',strLine);
        strLine = "&";
        endStrLine="";
    end
    fprintf(fileID,'%s\n',"\cline{2-"+num2str(lenX+3)+"} \\[-1em]");
    firstStrStart="";
end

fprintf(fileID,'%s\n'," \end{tabular}}");
fprintf(fileID,'%s\n'," \caption{Tabla de numero de modos para cada guía simulada a diferentes longitudes de onda} \label{table5_1}");
fprintf(fileID,'%s\n',"\end{table}");

fclose(fileID);