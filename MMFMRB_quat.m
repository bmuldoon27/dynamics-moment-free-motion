function [fig1, fig2, fig3] = MMFMRB_quat(dimensions,IC,T)
    %% MOMENT FREE MOTION OF A RIGID BODY USING QUATERNIONS
    % INPUTS -------------------------------------------------------------
    % dimensions = [a, b ,c] --- Vector of body dimensions
    % IC = [xdot0, x0, omega0, p0, omegadot0, pdot0]; -- initial conditions
    % WHERE:
    % xdot0     = 1 x 3 vector of Initial Center of Mass Velocity
    % x0        = 1 x 3 vector of Initial Center of Mass Position 
    % omega0    = 1 x 3 vector of Initial Angular Velocity Components
    % p0        = 1 x 4 vector of Initial Quaternion Components
    % omegadot0 = 1 x 4 vector of Initial Angular Acceleration Components
    % pdot0     = 1 x 4 vector of Intiial Quaternion Components Rates
    
    % OUTPUTS ------------------------------------------------------------
    % fig1 ---- figure of all solution parameters
    % fig2 ---- figure with a table of all initial conditions
    % fig3 ---- figure of the rigid body motion animation

    %% Set Output Subplot Dimensions
    nrowplot = 3;
    ncolplot= 4;
    %--------------------------------------------------------------------------
    %% Define Parallelipiped Geometry and Constants
    %--------------------------------------------------------------------------
    % Mass
    m = 5;        % kg
    % Gravitational Constant 
    g = 9.81;        % m/s^2
    % Define geometry (INPUTS)
    % Example;
    a = dimensions(1);     % m  -- length along 3 direction dimension 
    b = dimensions(2);     % m  -- length along 2 direction dimension
    c = dimensions(3);     % m  -- length along 1 direction dimension
    lambda1 = (m/12)*(b^2 + c^2); 	% kg-m^2 % LARGEST
    lambda2 = (m/12)*(a^2 + c^2); 	% kg-m^2
    lambda3 = (m/12)*(a^2 + b^2);  	% kg-m^2 % SMALLEST

    %% Set Simulation Parameters
    %--------------------------------------------------------------------------
    dt = 0.005;             % time step
    tf = T;                 % final time
    tsim = [0 : dt : tf]';  % initial mesh size
    tol = 1e-8;             % solver convergence tolerance
    options = odeset('abstol', tol, 'reltol', tol); % set the option changes
    %% Integrate Equations of Motion
    %--------------------------------------------------------------------------
    %IC = [xdot0, x0, omega0, p0, omegadot0,pdot0'];
    [T, Y] = ode45(@EOM, tsim, IC, options, m, g, ...
        lambda1, lambda2, lambda3);
    % Extract the position and orientation solutions:
    xdot= Y(:,1);
    ydot= Y(:,2);
    zdot= Y(:,3);
    x=Y(:,4);
    y=Y(:,5);
    z=Y(:,6);
    omega1 = Y(:,7);
    omega2 = Y(:,8);
    omega3 = Y(:,9);
    e0= Y(:,10);
    e1= Y(:,11);
    e2= Y(:,12);
    e3= Y(:,13);
    % Initialize arrays for corotational basis vectors over time:
    eb1 = zeros(3,length(e0));
    eb2 = zeros(3,length(e0));
    eb3 = zeros(3,length(e0));
    % Unpack the corotational basis vectors over time from the rotation matrix
    % -> The columns of the rotation matrix are the corotational basis vectors:
    for i= 1:length(e0)
        q0 = e0(i);
        q1 = e1(i);
        q2 = e2(i);
        q3 = e3(i);
        % The rotation is parameterized using unit quaternions (Euler-param.)
        B1 = (q0^2-q1^2-q2^2-q3^2)*eye(3);
        B2 = [ 2*q1*q1, 2*q1*q2, 2*q1*q3;
               2*q1*q2, 2*q2*q2, 2*q2*q3;
               2*q1*q3, 2*q2*q3, 2*q3*q3;];
        B3 = [   0     ,  -2*q0*q3,  2*q0*q2;
               2*q0*q3 ,     0    , -2*q0*q1;
               -2*q0*q2, 2*q0*q1  ,    0];
        Qmat = B1+B2+B3;
        % unpack basis vectors (cols of Q)
        eb1(:,i) = Qmat(:,1);
        eb2(:,i) = Qmat(:,2);
        eb3(:,i) = Qmat(:,3);
    end
    %% Initialize Figure of Subplots
    %--------------------------------------------------------------------------
    % make figure fill screen, fullfig - download from Matlab Website
    fig1 = fullfig; 
    set(gcf, 'color', 'w', 'name', 'Solution Parameters Window');
    %% Plot the evolution of the corotational basis by tracing the trajectories
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,[7 8 11 12]);
    [X, Y, Z] = sphere;
    sphere_obj = surf(X, Y, Z);
    set(sphere_obj, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05);
    tit1=title({'Trajectories traced by the corotational basis on the unit sphere:';...
        '$\bf{e}_1 =$ blue, $\bf{e}_2 =$ red, $\bf{e}_3 =$  black, $\bf{r}=$ cyan, $\bf{i}=$ green'});
    set(tit1,'Interpreter','latex','fontsize',15);
    xl1=xlabel('$\bf{E}_1$'); 
    yl1=ylabel('$\bf{E}_2$'); 
    zl1=zlabel('$\bf{E}_3$');
    set(xl1,'Interpreter','latex','fontsize',12);
    set(yl1,'Interpreter','latex','fontsize',12);
    set(zl1,'Interpreter','latex','fontsize',12,'rotation',0);
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    axis equal;
    hold on;
    % Draw the basis vectors 
    eb1v = quiver3(0,0,0, eb1(1,1), eb1(2,1), eb1(3,1), 'b', 'LineWidth', 2,...
        'AutoScale', 'off');
    eb2v = quiver3(0,0,0, eb2(1,1), eb2(2,1), eb2(3,1), 'r', 'LineWidth', 2,...
        'AutoScale', 'off');
    eb3v = quiver3(0,0,0, eb3(1,1), eb3(2,1), eb3(3,1), 'k', 'LineWidth', 2,...
        'AutoScale', 'off');
    % get the axis of rotation components over time
    sin_acos_e0=sin(acos(e0));
    rr1 = e1./sin_acos_e0;
    rr2 = e2./sin_acos_e0;
    rr3 = e3./sin_acos_e0;
    % draw vector for axis of rotation
    rrv = quiver3(0,0,0, rr1(1), rr2(1), rr3(1), 'c', 'LineWidth', 3,...
        'AutoScale', 'off');
    % get the angular velocity vector unit direction over time
    omega=[omega1';
           omega2';
           omega3'];
    % get angular velocity unit direction over time
    omega_axis = zeros(3,length(omega));
    for k = 1:1:length(omega)
        omega_axis(:,k) = omega(:,k)./norm(omega(:,k));
    end
    omegav = quiver3(0,0,0, omega_axis(1,1), omega_axis(2,1), omega_axis(3,1),...
        'g', 'LineWidth', 3,'AutoScale', 'off');
    % Draw their tips' trajectories:
    e1p = line(eb1(1,:), eb1(2,:), eb1(3,:), 'Color', 'b', 'LineWidth', 2);
    e2p = line(eb2(1,:), eb2(2,:), eb2(3,:), 'Color', 'r', 'LineWidth', 2);
    e3p = line(eb3(1,:), eb3(2,:), eb3(3,:), 'Color', 'k', 'LineWidth', 2);
    %rrp = line(rr1(:), rr2(:), rr3(:), 'Color','c', 'LineWidth', 2); % axis of rotation
    omegap = line(omega_axis(1,:), omega_axis(2,:), omega_axis(3,:), 'Color', 'g', 'LineWidth', 3);
    grid on;

    %% Plot e3, e2, e0 Quaternion Components
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,6);
    eparam=plot3(e3,e2,e0,'Linewidth',3); grid on;
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    axis([-1 1 -1 1 -1 1])
    zh2= get(gca,'zlabel');
    ax_eparam = eparam.Parent;   % Important
    set(ax_eparam, 'XTick', -1:0.5:1)
    set(ax_eparam, 'YTick', -1:0.5:1)
    set(ax_eparam, 'ZTick', -1:0.5:1)
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    set(gcf, 'Position',  [100, 500, 300, 300])
    ax_eparam.XTickLabel = {'-1','','$e_3$','','1'};
    ax_eparam.YTickLabel = {'-1','','$e_2$','','1'};
    ax_eparam.ZTickLabel = {'-1','','$e_0$','','1'};
    axis equal;

    %% Plot e1, e2, e0 Quaternion Components
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,10);
    eparam=plot3(e1,e2,e0,'Linewidth',3); grid on;
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    axis([-1 1 -1 1 -1 1])
    zh2= get(gca,'zlabel');
    ax_eparam = eparam.Parent;   % Important
    set(ax_eparam, 'XTick', -1:0.5:1)
    set(ax_eparam, 'YTick', -1:0.5:1)
    set(ax_eparam, 'ZTick', -1:0.5:1)
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    set(gcf, 'Position',  [100, 500, 300, 300])
    ax_eparam.XTickLabel = {'-1','','$e_1$','','1'};
    ax_eparam.YTickLabel = {'-1','','$e_2$','','1'};
    ax_eparam.ZTickLabel = {'-1','','$e_0$','','1'};
    axis equal;

    %% Plot Quaternions over Time
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,[1 2]); grid on;
    plot(T,e0,'Linewidth',2); hold on;
    plot(T,e1,'Linewidth',2); grid on;
    plot(T,e2,'Linewidth',2); 
    plot(T,e3,'Linewidth',2); 
    xl3=xlabel('Time'); 
    yl3=ylabel('Quaternion Value'); 
    set(xl3,'Interpreter','latex','fontsize',15);
    set(yl3,'Interpreter','latex','fontsize',15);
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    leg3=legend('$e_0$','$e_1$','$e_2$','$e_3$','Location','southoutside');
    leg3.NumColumns=4;
    set(leg3,'Interpreter','latex','fontsize',15);


    %% Plot Rotation Tensor Components 
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,5);
    R32=zeros(length(e1));
    R13=zeros(length(e1));
    R21=zeros(length(e1));
    for i = 1:1:length(e1)
        q0=e0(i); q1=e1(i); q2=e2(i); q3=e3(i);
        R1 = (q0^2-q1^2-q2^2-q3^2)*eye(3);
        R2 = [ 2*q1*q1, 2*q1*q2, 2*q1*q3;
               2*q1*q2, 2*q2*q2, 2*q2*q3;
               2*q1*q3, 2*q2*q3, 2*q3*q3;];
        R3 = [   0     ,  -2*q0*q3,  2*q0*q2;
               2*q0*q3 ,     0    , -2*q0*q1;
               -2*q0*q2, 2*q0*q1  ,    0];
        R = R1+R2+R3;
        R32(i) = R(3,2);
        R13(i) = R(1,3);
        R21(i) = R(2,1);
    end

    Rcomp=plot3(R32,R13,R21,'Linewidth',3); 
    zh4= get(gca,'zlabel');
    ax_Rcomp = Rcomp.Parent;   % Important
    set(ax_Rcomp, 'XTick', -1:0.5:1)
    set(ax_Rcomp, 'YTick', -1:0.5:1)
    set(ax_Rcomp, 'ZTick', -1:0.5:1)
    set(gca,'TickLabelInterpreter','latex','fontsize',20)
    set(gcf, 'Position',  [100, 100, 300, 300])
    ax_Rcomp.XTickLabel = {'-1','','$R_{32}$','','1'};
    ax_Rcomp.YTickLabel = {'-1','','$R_{13}$','','1'};
    ax_Rcomp.ZTickLabel = {'-1','','$R_{21}$','','1'};

    % Steiner's Roman Surface (if enforce e_1=0 case for RP^2 representation)
    PSI = linspace(0,2*pi,100); 
    THETA=linspace(pi,-pi,100);
    [psi,theta] = meshgrid(PSI,THETA);
    x1=sin(2*psi).*(sin(theta/2)).^2;
    x2=cos(psi).*sin(theta);
    x3=sin(psi).*sin(theta);
    hold on; grid on;
    mesh(x1,x2,x3,'FaceAlpha',0.05,'EdgeAlpha',0.1)


    %% Plot Angular Velocity Components Over Time
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,[3 4]);
    plot(T,omega1,'Linewidth',2); hold on; grid on;
    plot(T,omega2,'Linewidth',2); hold on; grid on;
    plot(T,omega3,'Linewidth',2); hold on; grid on;
    xl4a=xlabel('Time, $t$'); 
    yl4a=ylabel('$\omega_i$'); 
    set(xl4a,'Interpreter','latex','fontsize',15);
    set(yl4a,'Interpreter','latex','fontsize',15,'rotation',0);
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    leg4=legend('$\omega_1$','$\omega_2$','$\omega_3$','Location',...
                                                               'southoutside');
    leg4.NumColumns=3;
    set(leg4,'Interpreter','latex','fontsize',15);


    %% Plot Energy Over Time
    %--------------------------------------------------------------------------
    subplot(nrowplot,ncolplot,9);
    Ttrans = 1/2*m*(xdot.^2 + ydot.^2 + zdot.^2);
    Trot = 1/2*(lambda1*omega1.^2 + lambda2*omega2.^2 + lambda3*omega3.^2);
    Ttot = Ttrans + Trot;
    U = m*g*z;
    E = Ttot+U;

    plot(T,E/E(1),'Linewidth',2);
    xlE = xlabel('Time');
    ylE = ylabel('$\frac{E}{E_0}$','Rotation',0);
    titleE = title('Energy Over Time');
    set(xlE,'Interpreter','latex','fontsize',15);
    set(ylE,'Interpreter','latex','fontsize',15);
    set(titleE,'Interpreter','latex','fontsize',15);
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    axis equal; grid on;

    %% Plot A Figure with Initial Conditions
    fig2= figure;
    set(gcf, 'color', 'w', 'name', 'Inputs & Initial Conditions');
    dat =  {'a', a;
            'b', b;
            'c', c;
            'xdot_1(t=0)', IC(1); 
            'xdot_2(t=0)', IC(2); 
            'xdot_3(t=0)', IC(3); 
            'x_1(t=0)', IC(4);
            'x_2(t=0)', IC(5);
            'x_3(t=0)', IC(6); 
            'omega_1(t=0)', IC(7);...
            'omega_2(t=0)', IC(8);
            'omega_3(t=0)', IC(9);
            'e_0(t=0)', IC(10);...
            'e_1(t=0)', IC(11);...   
            'e_2(t=0)', IC(12);...
            'e_3(t=0)', IC(13);...
            'omegadot_1(t=0)', IC(14);...
            'omegadot_2(t=0)', IC(15);
            'omegadot_3(t=0)', IC(16);
            'edot_0(t=0)', IC(17);
            'edot_1(t=0)', IC(18);
            'edot_2(t=0)', IC(19);
            'edot_3(t=0)', IC(20)};
    columnname =   {'Parameter', 'Value'};
    columnformat = {'char', 'numeric', 'char'}; 
    t = uitable('Units','normalized','Position',...
                [0.05 0.05 0.755 0.87], 'Data', dat,... 
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'RowName',[]);

    %% Plot Animation of Motion
    %--------------------------------------------------------------------------
    % Set up the animation window:
    fig3=figure;
    set(gcf, 'color', 'w', 'name', 'Animation');

    axis equal;
    bound=max([a,b,c]);
    xlim([min(x)-bound, max(x)+bound]);
    ylim([min(y)-bound, max(y)+bound]);
    zlim([min(z)-bound, max(z)+bound]);
    xl5=xlabel('$x_1$ (m) '); 
    yl5=ylabel('$x_2$ (m) '); 
    zl5=zlabel('$x_3$ ', 'rotation', 0);
    set(xl5,'Interpreter','latex','fontsize',15);
    set(yl5,'Interpreter','latex','fontsize',15);
    set(zl5,'Interpreter','latex','fontsize',15);
    set(gca,'TickLabelInterpreter','latex','fontsize',15)
    view([135 30]);
    grid on;

    % We need to animate 6 planes to form a rectangular prism. Track 8 material 
    % points, chosen to be the vertices of the prism:

    rCM = [x, y, z]';

    r1 = rCM + 0.5*a*eb1 + 0.5*b*eb2 + 0.5*c*eb3;
    r2 = rCM - 0.5*a*eb1 + 0.5*b*eb2 + 0.5*c*eb3;
    r3 = rCM - 0.5*a*eb1 - 0.5*b*eb2 + 0.5*c*eb3;
    r4 = rCM + 0.5*a*eb1 - 0.5*b*eb2 + 0.5*c*eb3;
    r5 = rCM + 0.5*a*eb1 + 0.5*b*eb2 - 0.5*c*eb3;
    r6 = rCM - 0.5*a*eb1 + 0.5*b*eb2 - 0.5*c*eb3;
    r7 = rCM - 0.5*a*eb1 - 0.5*b*eb2 - 0.5*c*eb3;
    r8 = rCM + 0.5*a*eb1 - 0.5*b*eb2 - 0.5*c*eb3;

    % Calculate the vertices locations over time for the 6 planes:

    vertices1_x = [r1(1,:); r2(1,:); r3(1,:); r4(1,:)];
    vertices1_y = [r1(2,:); r2(2,:); r3(2,:); r4(2,:)];
    vertices1_z = [r1(3,:); r2(3,:); r3(3,:); r4(3,:)];
    surf_1 = patch(vertices1_x(:,1), vertices1_y(:,1), vertices1_z(:,1), 'FaceColor', 'c', 'FaceAlpha', 0.5);

    vertices2_x = [r5(1,:); r6(1,:); r7(1,:); r8(1,:)];
    vertices2_y = [r5(2,:); r6(2,:); r7(2,:); r8(2,:)];
    vertices2_z = [r5(3,:); r6(3,:); r7(3,:); r8(3,:)];
    surf_2 = patch(vertices2_x(:,1), vertices2_y(:,1), vertices2_z(:,1), 'FaceColor', 'c', 'FaceAlpha', 0.5);

    vertices3_x = [r3(1,:); r4(1,:); r8(1,:); r7(1,:)];
    vertices3_y = [r3(2,:); r4(2,:); r8(2,:); r7(2,:)];
    vertices3_z = [r3(3,:); r4(3,:); r8(3,:); r7(3,:)];
    surf_3 = patch(vertices3_x(:,1), vertices3_y(:,1), vertices3_z(:,1), 'FaceColor', 'c', 'FaceAlpha', 0.5);

    vertices4_x = [r1(1,:); r2(1,:); r6(1,:); r5(1,:)];
    vertices4_y = [r1(2,:); r2(2,:); r6(2,:); r5(2,:)];
    vertices4_z = [r1(3,:); r2(3,:); r6(3,:); r5(3,:)];
    surf_4 = patch(vertices4_x(:,1), vertices4_y(:,1), vertices4_z(:,1), 'FaceColor', 'c', 'FaceAlpha', 0.5);

    vertices5_x = [r1(1,:); r4(1,:); r8(1,:); r5(1,:)];
    vertices5_y = [r1(2,:); r4(2,:); r8(2,:); r5(2,:)];
    vertices5_z = [r1(3,:); r4(3,:); r8(3,:); r5(3,:)];
    surf_5 = patch(vertices5_x(:,1), vertices5_y(:,1), vertices5_z(:,1), 'FaceColor', 'c', 'FaceAlpha', 0.5);

    vertices6_x = [r2(1,:); r3(1,:); r7(1,:); r6(1,:)];
    vertices6_y = [r2(2,:); r3(2,:); r7(2,:); r6(2,:)];
    vertices6_z = [r2(3,:); r3(3,:); r7(3,:); r6(3,:)];
    surf_6 = patch(vertices6_x(:,1), vertices6_y(:,1), vertices6_z(:,1), 'FaceColor', 'c', 'FaceAlpha', 0.5);

    % Highlight one of the vertices:
    point6 = line(r6(1,1), r6(2,1), r6(3,1), 'marker', 'o', 'markerfacecolor', 'b');

    % Animate the body
    % ----Animation Video Save Options
    % animation = VideoWriter('tossed-book.avi');
    % animation.FrameRate = 1/dt/4;
    % open(animation);

    for jjj = 1:length(eb1)

        set(surf_1, 'xdata', vertices1_x(:,jjj), 'ydata', vertices1_y(:,jjj), 'zdata', vertices1_z(:,jjj));
        set(surf_2, 'xdata', vertices2_x(:,jjj), 'ydata', vertices2_y(:,jjj), 'zdata', vertices2_z(:,jjj));
        set(surf_3, 'xdata', vertices3_x(:,jjj), 'ydata', vertices3_y(:,jjj), 'zdata', vertices3_z(:,jjj));
        set(surf_4, 'xdata', vertices4_x(:,jjj), 'ydata', vertices4_y(:,jjj), 'zdata', vertices4_z(:,jjj));
        set(surf_5, 'xdata', vertices5_x(:,jjj), 'ydata', vertices5_y(:,jjj), 'zdata', vertices5_z(:,jjj));
        set(surf_6, 'xdata', vertices6_x(:,jjj), 'ydata', vertices6_y(:,jjj), 'zdata', vertices6_z(:,jjj));
        set(point6, 'xdata', r6(1,jjj), 'ydata', r6(2,jjj), 'zdata', r6(3,jjj));
        drawnow;
        % writeVideo(animation, getframe(gcf));

    end
    
    %% Create the 3 plot animations for rotations.berkeley.edu website
    % {Plot 1: Angle of Rotation over Time} 
    % {Plot 2: Axis of rotation, and instantaneous ang. vel. dir. over
    % time}
    % {Plot 3: Evolution of the corotational basis vectors over time.}

    phi = 2*acos(e0);
    
    iter0 = 1;
    iter = iter0;

    rwp = figure('Renderer', 'painters', 'Position', [100 100 800 200]);
    imageDirectory='gif_figures_asym';
    % setup sub plot 1
    subplot(1,3,1);
    plot(T(iter0),phi(iter0),'-b','Linewidth',1); hold on;
    xlim([0,max(T)]); ylim([min(phi), max(phi)]);

    while iter <= length(e0) % step through time and plot data at each snap
        % plot on the first figure axes
        subplot(1,3,1);
        plot(T(iter0:iter), phi(iter0:iter),'-b','Linewidth',2); drawnow;
        hold on;
        plot(T(iter), phi(iter),'o','MarkerFaceColor','b');
        xlim([0,max(T)]); ylim([0.875*min(phi), 1.125*max(phi)]); 
        grid on;
        xl1=xlabel('Time');
        yl1=ylabel('$\phi$','rotation',0);
        set(xl1,'Interpreter','Latex','Fontsize',15);
        set(yl1,'Interpreter','Latex','Fontsize',15);
        set(gca,'TickLabelInterpreter','latex','fontsize',12)
        
        % plot on the second figure axes
        sp2=subplot(1,3,2);
        omegav = quiver3(0,0,0, omega_axis(1,iter), omega_axis(2,iter), omega_axis(3,iter),...
        'b', 'LineWidth', 1,'AutoScale', 'off');
        omegap = line(omega_axis(1,iter0:iter), omega_axis(2,iter0:iter), omega_axis(3,iter0:iter), 'Color', 'b', 'LineWidth', 1);
        hold on; grid off;
        rrv = quiver3(0,0,0, rr1(iter), rr2(iter), rr3(iter), 'r', 'LineWidth', 1,...
        'AutoScale', 'off');
        rrp = line(rr1(iter0:iter), rr2(iter0:iter), rr3(iter0:iter), 'Color','r', 'LineWidth', 1); % axis of rotation
        axis([-1 1 -1 1 -1 1]);
        xl2=xlabel('$\bf{E}_1$');
        yl2=ylabel('$\bf{E}_2$');
        zl2=zlabel('$\bf{E}_3$','rotation',0);
        set(xl2,'Interpreter','Latex','Fontsize',12);
        set(yl2,'Interpreter','Latex','Fontsize',12);
        set(zl2,'Interpreter','Latex','Fontsize',12);
        ax_axes= get(gca,'zlabel');
        ax_axes_par = ax_axes.Parent;   % Important
        set(ax_axes_par, 'XTick', -1:0.5:1)
        set(ax_axes_par, 'YTick', -1:0.5:1)
        set(ax_axes_par, 'ZTick', -1:0.5:1)
        set(gca,'TickLabelInterpreter','latex','fontsize',12)

        % plot on third axes
        sp3=subplot(1,3,3);
        % Draw the basis vectors 
        eb1v = quiver3(0,0,0, eb1(1,iter), eb1(2,iter), eb1(3,iter), 'k', 'LineWidth', 1,...
        'AutoScale', 'off'); hold on;
        eb2v = quiver3(0,0,0, eb2(1,iter), eb2(2,iter), eb2(3,iter), 'k', 'LineWidth', 1,...
        'AutoScale', 'off');
        eb3v = quiver3(0,0,0, eb3(1,iter), eb3(2,iter), eb3(3,iter), 'k', 'LineWidth', 1,...
        'AutoScale', 'off');

        % Draw their tips' trajectories:
         e1p = line(eb1(1,iter0:iter), eb1(2,iter0:iter), eb1(3,iter0:iter), 'Color', 'k', 'LineWidth', 1);
         e2p = line(eb2(1,iter0:iter), eb2(2,iter0:iter), eb2(3,iter0:iter), 'Color', 'k', 'LineWidth', 1);
         e3p = line(eb3(1,iter0:iter), eb3(2,iter0:iter), eb3(3,iter0:iter), 'Color', 'k', 'LineWidth', 1);
        grid off; axis([-1 1 -1 1 -1 1]);
        xl3=xlabel('$\bf{E}_1$');
        yl3=ylabel('$\bf{E}_2$');
        zl3=zlabel('$\bf{E}_3$','rotation',0);
        set(xl3,'Interpreter','Latex','Fontsize',12);
        set(yl3,'Interpreter','Latex','Fontsize',12);
        set(zl3,'Interpreter','Latex','Fontsize',12);
        ax_bv= get(gca,'zlabel');
        ax_bv_par = ax_bv.Parent;   % Important
        set(ax_bv_par, 'XTick', -1:0.5:1)
        set(ax_bv_par, 'YTick', -1:0.5:1)
        set(ax_bv_par, 'ZTick', -1:0.5:1)
        set(gca,'TickLabelInterpreter','latex','fontsize',12)
        
  
        % Say Cheese :D
        saveas(rwp, num2str(iter-1, [imageDirectory, '/iter=%f.png']))

        % Clear figs
        for i = 1:1:3; subplot(1,3,i); clf; end % clear each 
        
        % time step iteration increase
        iter = iter+1;
    end

end
