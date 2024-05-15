% Autor: Guillermo Fernando Regodón Harkness
% Github: https://github.com/GuillermoRegodon/matlab-class-BlochEqSolver.git

% Clase de Matlab para resolver y mostrar en figuras la solución de las
% ecuaciones de Bloch.

% Los objetos son una estructura de datos más las funciones permitidas.
% En Matlab, significa que tenemos que devolver siempre la estructura y
% guardarla. Si no, la computación se pierde (aunque otros efectos, como 
% plot sí puedan ser observables).

% La clase viene acompañada de un LiveScript de Matlab para ilustrar
% conceptos de respuesta en frecuencia de la magnetización de un material
% sometido a polarización circular

classdef BlochEqSolver
    properties
        tau1
        tau2
        B0
        B1
        w
        M0
        tspan
        Mi
        f
        t
        M
    end
    methods
        function obj = BlochEqSolver(tau1, tau2, B0, B1, w, chi0, tmax, dt, Mi)
            % gamma = 1;
            % mu0 = 1;
            if length(Mi) == 3
                obj.tau1 = tau1;
                obj.tau2 = tau2;
                obj.B0 = B0;
                obj.B1 = B1;
                obj.w = w;
                obj.M0 = chi0/(1+chi0)*B0;
                obj.tspan = 0:dt:tmax;
                obj.Mi = Mi;
                obj.f = @(t,M)obj.bloch(t, M);
                obj = obj.ode();
            end
        end

        function Mp = bloch(obj, t, M)
            Mp = -cross(M, [obj.B1*cos(obj.w*t); obj.B1*sin(obj.w*t); obj.B0]) + ...
                + [-M(1)/obj.tau2; -M(2)/obj.tau2; (obj.M0 - M(3))/obj.tau1];
        end

        function obj = ode(obj)
            [obj.t, obj.M] = ode89(obj.f, obj.tspan, obj.Mi);
        end
        
        function obj = plot(obj, tmax)
            plot(obj.t, obj.M)
            legend({"Mx", "My", "Mz"})
            xlabel("Time")
            if nargin > 1
                xlim([0 tmax])
            end
        end

        function obj = plot3(obj)
            plot3(obj.M(:,1), obj.M(:,2), obj.M(:,3))
            hold on
            plot3(0, 0, 0, "+k");
            quiver3(0, 0, 0, obj.M(end,1), obj.M(end,2), obj.M(end,3), ...
                0, "Color", [48, 112, 183]/255)
            hold off
            axis(1.2*[max(max(obj.M(:, 1:2)))*[-1 1 -1 1] ...
                min([min(obj.M(:,3)) 0]) max([max(obj.M(:,3)) 0])])
            pbaspect([1 1 1])
            grid on
            xlabel("x")
            ylabel("y")
            zlabel("z")
        end

        function obj = fft(obj, wmax)
            dt = obj.t(2)-obj.t(1);
            if nargin < 2
                wmax = pi/dt;
            else
                wmax = min(pi/dt, abs(wmax));
            end
            ft = abs(fft(obj.M));
            Lmax = length(ft(:,1));
            L = floor(Lmax/2);
            wx = linspace(0,2*pi,Lmax+1);
            wx = wx(1:L)/dt;
            plot(wx(1:L), ft(1:L, :))
            axis([0 wmax 0 1.2*max(max(ft))])
            legend({"|fft(Mx)|", "|fft(My)|", "|fft(Mz)|"})
        end

        function obj = plotd(obj, N, inc)
            if nargin < 3
                N = 1;
                inc = 1;
            end
            for i = 0:inc:N*inc-1
                Mt = [obj.M(end-i,1) obj.M(end-i,2)];
                Mt = Mt/sqrt(Mt*Mt');
                quiver(0, 0, Mt(1), Mt(2), ...
                    0, "Color", [48, 112, 183]/255)
                hold on
                quiver(0, 0, cos(obj.w*obj.t(end-i)), sin(obj.w*obj.t(end-i)), ...
                    0, "Color", [201, 92, 46]/255)
            end
            xlabel("Mx, Bx")
            ylabel("My, By")
            legend({"M", "B"})
            axis(1.1*[-1 1 -1 1])
            pbaspect([1 1 1])
            hold off
        end

    end
end