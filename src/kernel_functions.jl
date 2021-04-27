# simple kernel_function

function kappa(x1,x2,y1,y2,rho)
    d = sqrt((x1-y1)^2+(x2-y2)^2);
    return (1+d/rho+d^2/(3*rho^2))*exp(-d/rho)
end
  
function Delta_x_kappa(x1,x2,y1,y2,rho)
    d = sqrt((x1-y1)^2+(x2-y2)^2); 
    return (d^2-2*rho*(rho+d))/(3*rho^4)*exp(-d/rho)
end

function Delta_x_y_kappa(x1,x2,y1,y2,rho)
    d = sqrt((x1-y1)^2+(x2-y2)^2);
    return (8*rho^2+d^2-7*rho*d)/(3*rho^6)*exp(-d/rho)
end