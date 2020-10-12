% Script Soft disk solver

clear



np=50; % number of points of the calculation
tmax=1e-11; % maximum integration time
h=1e-14; % step size of the integration


% Potential parameters

E=2;
v0=5e5; % initial (Fermi) velocity
me=9.10938356e-31; % electron (effective) mass
T = 10; % temperature
c=2.5e4; % speed of phonons
Ed=2e-7; % potential strength
kxmin=-2.5e8; % phonon wavevectors definition
kxmax=2.5e8;
ncomps=50;
% for(var j=kxmin, var k=0; j<kxmax; j+=(kxmax-kxmin)/(ncomps-1), k++){
%   kx[k]=j; ky[k]=j;
%}
kx=kxmin:(kxmax-kxmin)/(ncomps-1):kxmax;
ky=kx;

% var phi = new Array(ncomps*ncomps);
% for(var j=0; j<ncomps;j++){
%  for(var k=0; k<ncomps; k++){
%    phi[j+ncomps*k]=2*Math.pi()*Math.random();
%}
%}
phi=2*pi*rand(ncomps,ncomps); % random phases for the phonons

% Plotting the potential

npV=200; % resolution of the potential picture
xmin=0;
xmax=2e-6;
% for(var j=xmin, var k=0; j<xmax; j+=(xmax-xmin)/(npV-1), k++){
%   xg[k]=j;
%}
xg=xmin:(xmax-xmin)/(npV-1):xmax;
ymin=-1e-6;
ymax=1e-6;
% for(var j=ymin, var k=0; j<ymax; j+=(ymax-ymin)/(npV-1), k++){
%   yg[k]=j;
%}
yg=ymin:(ymax-ymin)/(npV-1):ymax;

V = Vdefpot (T, kx, ky, xg, yg, phi, c, Ed);

figure()
pcolor(xg,yg,V')
xlabel('x')
ylabel('y')
shading flat
colormap(flipud(hot))
axis square
hold on


% Initial conditions

x0=0;
y0min=-0.2e-6;
y0max=0.2e-6;
% for(var j=y0min, var k=0; j<y0max; j+=(y0max-y0min)/(np-1), k++){
%   y0[k]=j;
%}
y0= y0min:(y0max-y0min)/(np-1):y0max; %initial conditions

beam=0; % to switch between beam and source modes

anglemin=-15;
anglemax=15;
% for(var j=anglemin, var k=0; j<anglemax; j+=(anglemax-anglemin)/(np-1), k++){
%   angle[k]=j;
%}
angle=anglemin:(anglemax-anglemin)/(np-1):anglemax;

if beam==1
    param=y0;
else
    param=angle;
end



for i=1:length(param)

    remaining =  length(param) - i

    tic
    %  var trajx = new Array(round(tmax/h)) %-876*ones(1,round(tmax/h));
    trajx=-876*ones(1,round(tmax/h));
    trajy=-876*ones(1,round(tmax/h));
    trajvx=-876*ones(1,round(tmax/h));
    trajvy=-876*ones(1,round(tmax/h));
    trajE=-876*ones(1,round(tmax/h));
    trajT=-876*ones(1,round(tmax/h));
    trajV=-876*ones(1,round(tmax/h));

    x=zeros(1,round(tmax/h));
    vx=zeros(1,round(tmax/h));
    ax=zeros(1,round(tmax/h));
    y=zeros(1,round(tmax/h));
    vy=zeros(1,round(tmax/h));
    ay=zeros(1,round(tmax/h));


    if beam==1
        x(1)=x0;
        y(1)=y0(i);
        vx(1)=v0;
        vy(1)=0;
    else
        x(1)=x0;
        y(1)=0;
        vx(1)=v0*cosd(angle(i));
        vy(1)=v0*sind(angle(i));
    end

     V0 = Vdefpot (T, kx, ky, x(1), y(1), phi, c, Ed);
     E0 = V0 + 0.5*me*v0*v0;
     v0=sqrt(2/me*(E0-V0)-vy(1)*vy(1));


    t=0;
    n=0;
    escape=0;


    % Start of symplectic integration

    while n<round(tmax/h)-1 && escape==0

        t=t+h;
        n=n+1;

        x(n+1) = x(n) + vx(n)*h + 0.5*ax(n)*h*h;
        y(n+1) = y(n) + vy(n)*h + 0.5*ay(n)*h*h;


        ax(n+1) = -dxVdefpot (T, kx, ky, x(n+1), y(n+1), phi, c, Ed)/me;
        ay(n+1) = -dyVdefpot (T, kx, ky, x(n+1), y(n+1), phi, c, Ed)/me;

        dxVdefpot (T, kx, ky, x(n+1), y(n+1), phi, c, Ed)
        dyVdefpot (T, kx, ky, x(n+1), y(n+1), phi, c, Ed)

        vx(n+1) = vx(n) +0.5 * (ax(n)+ax(n+1))*h;
        vy(n+1) = vy(n) +0.5 * (ay(n)+ay(n+1))*h;

        x(n+1)
        y(n+1)
        ax(n+1)
        ay(n+1)
        vx(n+1)
        vy(n+1)

pause()

        trajx(n)=x(n+1);
        trajy(n)=y(n+1);
        trajvx(n)=vx(n+1);
        trajvy(n)=vy(n+1);
%       trajT(n)=0.5*me*(vx(n+1)*vx(n+1)+vy(n+1)*vy(n+1)) ;
%        trajV(n)= VSoftDiskpoint(alpha, r, x01, y01, x(n+1),y(n+1)) + VSoftDiskpoint(alpha, r, x02, y02, x(n+1),y(n+1)) + VSoftDiskpoint(alpha, r, x03, y03, x(n+1),y(n+1)) + VSoftDiskpoint(alpha, r, x04, y04, x(n+1),y(n+1));
  %      trajE(n)=trajT(n)+trajV(n);

        if x(n+1)>xmax || y(n+1)>ymax || y(n+1)<ymin
            escape=1;
        end

    end % End of symplectic integration


        trajE=trajE(trajE~=0);
        trajT=trajT(trajT~=0);
        trajV=trajV(trajV~=0);

        trajx=trajx(trajx~=-876);
        trajy=trajy(trajy~=-876);

        % Trajectories in the potential picture

        plot(trajx,trajy,'-b','LineWidth',1)
        drawnow
        % Energy conservation picture
        % figure(),plot(trajT,'b'),hold on, plot(trajV,'r'), plot(trajE,'k'), hold off


    toc

end % end of y0
hold off
