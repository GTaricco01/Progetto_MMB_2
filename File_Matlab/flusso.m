function f = flusso(u,v,a,lam,s)

% scelta dello schema numerico

switch s
    case 1
        f = eulero(u,v,a);
    case 2
        f = laxfri(u,v,a,lam);
    case 3
        f = upwind(u,v,a);
    case 4
        f = laxwen(u,v,a,lam);
end

return
