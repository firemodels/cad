% M. Vanella
% Aug 2019
% mesh_stretch.m
function [np,xpl]=mesh_stretch(xs,Lx,m,DXS,Nx)

% Define tolerance:
etol_rel = 1.e-10;
max_iter = 1000;

for iside=1:2
   disp(' ')
   if(iside==1)
      disp(['Computing Low side parameter.'])
   else
      disp(['Computing High side parameter.'])
   end
    % Newton iteration:
    %t=cputime;
    an(iside) = 1.;
    for n=1:max_iter
        % Compute f, fp:
        suma = 0;;
        sumap= 0.;
        for i=m(iside)+1:Nx(iside)
            suma = suma  + an(iside)^(i-m(iside));
            sumap= sumap + (i-m(iside))*an(iside)^(i-m(iside)-1);
        end
        f_an = (m(iside)-Lx(iside)/DXS(iside)) + suma;
        fp_an= sumap;

        % New value of a:
        an1 = an(iside) - f_an/fp_an;

        % Convergence test:
        if(abs(an1-an(iside))/an(iside) < etol_rel)
            disp(['Iter ' num2str(n,'%4.4d') ...
                 ', Convergence found, stretch factor a=' num2str(an1) ...
                 ', relative error=' num2str(abs(an1-an(iside))/an(iside)) '.'])
            break
        end
%         if(mod(n,ceil(max_iter/1000)) == 0)
%             disp(['Iter ' num2str(n,'%4.4d') ...
%                  ', Relative error=' num2str(abs(an1-an(iside))/an(iside),'%18.12f') '.'])
%         end
        % Update an:
        an(iside)=an1;
    end
    %disp(['Time taken :' num2str(cputime-t) ' sec.'])
end

% Build xpl:
Nxt = Nx(1)+Nx(2);
Lxt = Lx(2)-Lx(1);

xplm(1) = xs;
% Low side:
iside=1;
for i=1:m(iside)
   xplm(i+1) =xs+i*DXS(iside);
end
for i=m(iside)+1:Nx(iside)
    dx = DXS(iside);
    sdx= 0.;
    for k=m(iside)+1:i
        dx = an(iside)*dx;
        sdx= sdx + dx;
    end
    x = m(iside)*DXS(iside) + sdx;
    xplm(i+1) =xs+x;
end
xplm(Nx(iside)+1) =xs+Lx(1);

% High side:
iside=2;
xplp(1) = xs;
for i=1:m(iside)
   xplp(i+1) =xs+i*DXS(iside);
end
for i=m(iside)+1:Nx(iside)
    dx = DXS(iside);
    sdx= 0.;
    for k=m(iside)+1:i
        dx = an(iside)*dx;
        sdx= sdx + dx;
    end
    x = m(iside)*DXS(iside) + sdx;
    xplp(i+1) =xs+x;
end
xplp(Nx(iside)+1) =xs+Lx(iside);

% Now final xpl, xipl:
xpl =  [xplm(end:-1:1) xplp(2:end)];

np = Nxt + 1;

return

end