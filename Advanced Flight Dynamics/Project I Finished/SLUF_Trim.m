function [AOA,del_t,theta,Vel_0] = SLUF_Trim(m,g,rho,V,s,Isp,Ca,Ac)
% Solve f(X_trim,U_trim)=0 for the acceleration states, using the existing
% body-axis force and kinetics functions. Angles are in degrees.

% Unknowns: angle of attack and throttle. For wings-level SLUF, theta=AOA,
% beta=phi=psi=0, p=q=r=0, and body velocity is [V*cosd(AOA);0;V*sind(AOA)].
% These constraints give zero lateral/vertical earth velocity. North
% position is allowed to increase: dNorth/dt=V, not zero.

options=optimoptions('fsolve','Display','iter',...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12);

[sol,residual,exitflag]=fsolve(@balance,[5;0.05],options); %input data and starting guesses

assert(exitflag>0 && norm(residual,inf)<1e-8,'SLUF trim did not converge.'); %sanity check

AOA=sol(1); del_t=sol(2);

assert(del_t>=0 && del_t<=1,'Required throttle is outside [0,1].'); %sanity check

theta=AOA; %for SLUF

Vel_0=[V*cosd(AOA);0;V*sind(AOA)];

    function residual=balance(unknowns)
        
        alpha=unknowns(1); throttle=unknowns(2); %AoA in degrees, Throttle fraction
        
        velocity=[V*cosd(alpha);0;V*sind(alpha)];% u = VcosAOA; v =0, and Vsin AOA
        
        C=lift_drag_fun(alpha);
        
        qs=0.5*rho*V^2*s; %Command dynamic pressure times area oftrim condition N
        
        a=speed_of_sound(15); % Same constant temperature as original mode
        % l
        Ft=thrust([throttle,Isp,g,rho,V/a,a,Ca,Ac]);

        F=Force_Total([C(1),C(2),qs,alpha,0,Ft]);
        
        acceleration=ThreeDOF_Kinetics(...
            [m;g;0;0;0;velocity;alpha;0;F(:)]);

        % Lateral acceleration is identically zero for the chosen symmetry.
        residual=acceleration([1,3]); % Require du/dt=0 and dw/dt=0
    end
end