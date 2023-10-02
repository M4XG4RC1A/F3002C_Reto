
%% getmodes
clc;

Xs = linspace(650,1000,10);
Ys = linspace(55,1000,10);
Ls = ["\lambda", "A_{eff}"];

lenX = length(Xs);
lenY = length(Ys);
lenL = length(Ls);

fileID = fopen("Table.txt",'w');

%Start table
fprintf(fileID,'%s\n',"\begin{table}[H]\centering");
fprintf(fileID,'%s\n',"\resizebox{0.45\textwidth}{!}{");
fprintf(fileID,'%s\n',"  \begin{tabular}{rr|"+convertCharsToStrings(repmat('c',1,lenX))+"|l}");
fprintf(fileID,'%s\n',"     & &\multicolumn{"+num2str(lenX)+"}{c}{Wg_x}\\");

strLine = "     & ";
for x = Xs
    strLine = strLine+"&"+num2str(x)+"nm";
end
strLine = strLine+"\\ \cline{1-"+num2str(lenX+4)+"} \\[-1em]";
fprintf(fileID,'%s\n',strLine);

firstStrStart = "     \multirow{"+num2str(lenL*lenY)+"}{*}{Wg_y}";
endStrLine = "";%"&\multirow{"+num2str(lenL*lenY)+"}{*}{$\lambda$}";

nmodes = 0;
for y = Ys
    strLine = firstStrStart+"&\multirow{"+num2str(lenL)+"}{*}{"+num2str(y)+"nm}";
    for l = Ls
        for x=Xs
            load("./nModes/Waveguide"+string(x*1000)+"_"+string(y*1000)+".mat");
            if l == "\lambda"
                strLine = strLine+"&"+num2str(nmodes)+"("+num2str(x*y)+")";
            else
                line = sprintf('%.0f,' , nonzeros(Aeff));
                line  = line(1:end-1);
                strLine = strLine+"&"+line;
            end
        end
        strLine = strLine+"&"+l+"nm"+endStrLine+"\\";
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