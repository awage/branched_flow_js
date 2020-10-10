function dxV = dxVdefpot (T, kx, ky, x, y, phi, c, Ed)

hbar = 1.0545718e-34;
kB = 1.3806e-23;

comps1=length(kx);
comps2=length(ky);
dxV = zeros(length(x), length(y));
tau = 0;

for i=1:comps1
    
    for j=1:comps2
        
        k = sqrt(kx(i)*kx(i) + ky(j)*ky(j));
        w = c *k; % linear regime where the frequency is proportional to the modulus of vector k
        
        dxVnew = -kx(i)*Ed/c*sqrt(2*hbar*w) / sqrt(exp(hbar*w/ kB /T)-1) * sin(kx(i)*x + ky(j)*y' - w*tau + phi(i,j));
        
        if isnan(dxVnew)==1
            dxVnew = zeros(length(x), length(y)); % to avoid dividing 0 by 0
        end
        
        dxV = dxV +dxVnew';

    end
end



