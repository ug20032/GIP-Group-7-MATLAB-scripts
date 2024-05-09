% Drag Coefficients
function Cd = TOMIYAMA_DragCoeff(Re,Eo,n)
    % Tomiyama, 1998 "Drag coefficients of single bubbles under normal and micro gravity conditions."
    % https://www.jstage.jst.go.jp/article/jsmeb1993/41/2/41_2_472/_pdf/-char/ja

    Cd(1) = max(min(16/Re*(1 + 0.15*Re^0.687),48/Re),8/3*Eo/(Eo+4)); % pure water
    Cd(2) = max(min(24/Re*(1 + 0.15*Re^0.687),72/Re),8/3*Eo/(Eo+4)); % slightly contaminated
    Cd(3) = max(24/Re*(1 + 0.15*Re^0.687),8/3*Eo/(Eo+4)); % contaminated

    Cd = Cd(n);
end
