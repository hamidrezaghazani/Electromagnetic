clc
clear
close all;

%By <<Hamidreza Sanati Ghazani>> 

%Tel_id: @hamid_sg

eps0 = 1/(36*pi);



waitfor(helpdlg('Hi! This is a program for calculation Electric field And potentiel!'))

h1 = menu('Choose Type of charges:','continuous','discrete');
opts =struct('WindowStyle','modal','Interpreter','tex','Resize','on');
if h1 == 2
h2 = inputdlg('Enter number of charges:','Charges Num',1,{'2'},opts);
h2 = str2num(h2{:});


h3 = menu('Choose Type of coodinate:','Cartesian','Polar','Sphere');

if h3 == 1
        if h2 == 2
            promt = {'q1(nC)','X','Y','Z';'q2(nC)','X','Y','Z'};
            defans ={'1','0','0','1';'-1','0','0','-1'};
        else
            promt  = cell(h2,4);
            defans = cell(h2,4);
            for i = 1:h2
                promt{i,1} = ['q',num2str(i),'(nC)'];
                promt{i,2} = 'X';
                promt{i,3} = 'Y';
                promt{i,4} = 'Z';
                defans{i,1}= '1';
                defans{i,2}= '2';
                defans{i,3}= '3';
                defans{i,4}= '4';
            end
        end
        answer = inputdlgcol(promt,'Charges',1,defans,opts,4);
        for i = 1:size(answer,1)
            answer1(i) = str2num(answer{i});
        end
        clear answer;
        answer = answer1;
        clear answer1;
        answer = reshape(answer,h2,4);


        q = answer(:,1)';
        x = answer(:,2)';
        y = answer(:,3)';
        z = answer(:,4)';

        answer_p = inputdlgcol({'X','Y','Z'},'Enter P info',1,{'0','10','0'},opts,3);
        for i = 1:size(answer_p,1)
            answer1(i) = str2num(answer_p{i});
        end
        clear answer_p;
        answer_p = answer1;
        clear answer1;
        p = reshape(answer_p,1,3);
elseif h3 == 2
        if h2 == 2
            promt = {'q1(nC)','r','phi','Z';'q2(nC)','r','phi','Z'};
            defans ={'1','0','0','1';'-1','0','0','-1'};
        else
            promt  = cell(h2,4);
            defans = cell(h2,4);
            for i = 1:h2
                promt{i,1} = ['q',num2str(i),'(nC)'];
                promt{i,2} = 'r';
                promt{i,3} = 'phi';
                promt{i,4} = 'Z';
                defans{i,1}= '1';
                defans{i,2}= '2';
                defans{i,3}= '3';
                defans{i,4}= '4';
            end
        end
        answer = inputdlgcol(promt,'Charges',1,defans,opts,4);
        for i = 1:size(answer,1)
            answer1(i) = str2num(answer{i});
        end
        clear answer;
        answer = answer1;
        clear answer1;
        answer = reshape(answer,h2,4);


        q = answer(:,1)';
        x = answer(:,2)';
        y = answer(:,3)';
        z = answer(:,4)';

        answer_p = inputdlgcol({'r','phi','Z'},'Enter P info',1,{'10','pi/2','0'},opts,3);
        for i = 1:size(answer_p,1)
            answer1(i) = str2num(answer_p{i});
        end
        clear answer_p;
        answer_p = answer1;
        clear answer1;
        p = reshape(answer_p,1,3);
        a = [x;y;z];
        for i = 1:size(a,2)
             b(:,i) = [a(1,i)*cos(a(2,i));a(1,i)*sin(a(2,i));a(3,i)];
        end
        b = round(b.*10^4)/(10^4);
        
        x = b(1,:);
        y = b(2,:);
        z = b(3,:);
        p = [p(1)*cos(p(2));p(1)*sin(p(2));p(3)]';
        
        
elseif h3 == 3
        if h2 == 2
            promt = {'q1(nC)','R','theta','phi';'q2(nC)','R','theta','phi'};
            defans ={'1','1','0','0';'-1','1','pi','0'};
        else
            promt  = cell(h2,4);
            defans = cell(h2,4);
            for i = 1:h2
                promt{i,1} = ['q',num2str(i),'(nC)'];
                promt{i,2} = 'R';
                promt{i,3} = 'theta';
                promt{i,4} = 'phi';
                defans{i,1}= '1';
                defans{i,2}= '2';
                defans{i,3}= '3';
                defans{i,4}= '4';
            end
        end
        answer = inputdlgcol(promt,'Charges',1,defans,opts,4);
        for i = 1:size(answer,1)
            answer1(i) = str2num(answer{i});
        end
        clear answer;
        answer = answer1;
        clear answer1;
        answer = reshape(answer,h2,4);


        q = answer(:,1)';
        x = answer(:,2)';
        y = answer(:,3)';
        z = answer(:,4)';

        answer_p = inputdlgcol({'R','theta','phi'},'Enter P info',1,{'10','pi/2','pi/2'},opts,3);
        for i = 1:size(answer_p,1)
            answer1(i) = str2num(answer_p{i});
        end
        clear answer_p;
        answer_p = answer1;
        clear answer1;
        p = reshape(answer_p,1,3);
        a = [x;y;z];
        for i = 1:size(a,2)
             b(:,i) = [a(1,i)*sin(a(2,i))*cos(a(3,i));a(1,i)*sin(a(2,i))*sin(a(3,i));a(1,i)*cos(a(2,i))];
        end
        b = round(b.*10^4)/(10^4);
        
        x = b(1,:);
        y = b(2,:);
        z = b(3,:);
        p = [p(1)*sin(p(2))*cos(p(3));p(1)*sin(p(2))*sin(p(3));p(1)*cos(p(2))]';
        
end

% x = [ 0 , 0 , 0];
% y = [ 0 , 0 , 0];
% z = [ 1 , -1 , 0];
% q = [ 1 , -1 , 8];
% p = [ 0 , 8 , 0 ];

for i = 1:size(q,2)
    r(:,i) = [ p(1)-x(i);
               p(2)-y(i);
               p(3)-z(i)];
end
for i = 1:size(q,2)
    v(i) = q(i)/((4*pi*eps0)*(norm(r(:,i))));
    Ex(i,1) =  q(i)/((4*pi*eps0)*(norm(r(:,i)))^3)*r(1,i);
    Ey(i,1) =  q(i)/((4*pi*eps0)*(norm(r(:,i)))^3)*r(2,i);
    Ez(i,1) =  q(i)/((4*pi*eps0)*(norm(r(:,i)))^3)*r(3,i);
end
E_Total = [ sum(Ex,1) , sum(Ey,1) , sum(Ez,1) ]
V_Total = sum(v)

E_hat = [E_Total(1)/norm(E_Total),E_Total(2)/norm(E_Total),E_Total(3)/norm(E_Total)];
E_hat = E_hat + p;

for i = 1:size(q,2)
    if q(i)<0
        plot3(x(i),y(i),z(i),'o','markersize',13,'markerfacecolor','r')
        hold on;
    else
        plot3(x(i),y(i),z(i),'o','markersize',13,'markerfacecolor','b')
        hold on;
    end
end
hold on;
grid on;
plot3(p(1),p(2),p(3),'s','markersize',13,'markerfacecolor','y')
plot3([p(1),E_hat(1)],[p(2),E_hat(2)],[p(3),E_hat(3)],'linewidth',3)
xlabel('X');ylabel('Y');zlabel('Z');
%legend('q','p','Field');

for i = 1:size(q,2)
    plot3([0,x(i)],[0,y(i)],[0,z(i)],'-g');
    hold on;
end
elseif h1 ==1

        opts =struct('WindowStyle','modal','Interpreter','tex','Resize','on');
        h4 = menu('Type of continuous charges:','A line charge symmetrical in axis z','A ring charge symmetrical in plane xy');
        if h4 == 1
                a = inputdlgcol({'length'},'Linear charge length',1,{'2'},opts);
                a = str2num(a{:});
                a = a/2;
                rho = inputdlgcol({'Rho (nC/m)'},'Linear charge density',1,{'1'},opts);
                rho = str2num(rho{:});
                answer_p = inputdlgcol({'X','Y','Z'},'Enter P info',1,{'0','10','0'},opts,3);
                for i = 1:size(answer_p,1)
                    answer1(i) = str2num(answer_p{i});
                end
                clear answer_p;
                answer_p = answer1;
                clear answer1;
                p = reshape(answer_p,1,3);

                syms x y z L rh eps0

        %         p = [ 0 , 0 , 10 ];
        %         a = 2;
        %         rho = 2;

                v = log((z+L+(x^2+y^2+(z+L)^2)^0.5)/((z-L+(x^2+y^2+(z-L)^2)^0.5)));
                Ex = simplify(diff(v,x));
                Ey = simplify(diff(v,y));
                Ez = simplify(diff(v,z));
                m = rh/(4*pi*eps0);

                disp('V');
                pretty(m*v)
                fprintf('\n\n\n')
                disp('Vp =')

                disp(double(subs(subs(subs(subs(subs(subs(m*v,x,p(1)),y,p(2)),z,p(3)),eps0,1/(36*pi)),rh,rho),L,a)))

                fprintf('\n\n\n')
                disp('Ex');
                pretty(-m*Ex)
                fprintf('\n\n\n')
                disp('Exp =')
                Exp1 = double(subs(subs(subs(subs(subs(subs(-m*Ex,x,p(1)),y,p(2)),z,p(3)),eps0,1/(36*pi)),rh,rho),L,a));
                disp(Exp1)



                fprintf('\n\n\n')
                disp('Ey');
                pretty(-m*Ey)
                fprintf('\n\n\n')
                disp('Eyp =')
                Eyp1 = double(subs(subs(subs(subs(subs(subs(-m*Ey,x,p(1)),y,p(2)),z,p(3)),eps0,1/(36*pi)),rh,rho),L,a));
                disp(Eyp1)




                fprintf('\n\n\n')
                disp('Ez');
                pretty(-m*Ez)
                fprintf('\n\n\n')
                disp('Ezp =')
                Ezp1 = double(subs(subs(subs(subs(subs(subs(-m*Ez,x,p(1)),y,p(2)),z,p(3)),eps0,1/(36*pi)),rh,rho),L,a));
                disp(Ezp1)

                E_Total = [Exp1,Eyp1,Ezp1];
                E_hat = [E_Total(1)/norm(E_Total),E_Total(2)/norm(E_Total),E_Total(3)/norm(E_Total)];
                E_hat = E_hat + p;

                plot3([0,0],[0,0],[-a,a],'-b','linewidth',5)
                hold on;
                plot3(p(1),p(2),p(3),'s','markersize',13,'markerfacecolor','y')
                plot3([p(1),E_hat(1)],[p(2),E_hat(2)],[p(3),E_hat(3)],'linewidth',3)
                xlabel('X');ylabel('Y');zlabel('Z');
                grid on;

        elseif h4 == 2

            waitfor(helpdlg('This program calculate Electric field And potantiel in Z axis P(0,0,z)!!'))
            aa = inputdlgcol({'radius'},'Linear charge radius',1,{'1'},opts);
            aa = str2num(aa{:});
            rho = inputdlgcol({'Rho (nC/m)'},'Linear charge density',1,{'1'},opts);
            rho = str2num(rho{:});
            zz = inputdlgcol({'Z'},'P point height',1,{'5'},opts);
            zz = str2num(zz{:});

            syms a z Rho eps0;
            V = (Rho*2*pi*a)/(4*pi*eps0*(a^2+z^2)^0.5);
            Ez = (Rho*2*pi*a*z)/(4*pi*eps0*(a^2+z^2)^(1.5));
            fprintf('\n\n\n')
            disp('V');
            pretty(V)
            fprintf('\n\n\n')
            disp('V =')
            disp(double(subs(subs(subs(subs(V,a,aa),Rho,rho),z,zz),eps0,1/(36*pi))));




            fprintf('\n\n\n')
            disp('Ez');
            pretty(Ez)
            fprintf('\n\n\n')
            fprintf('\n\n\n')
            disp('Ez =')

            Ez1 = double(subs(subs(subs(subs(Ez,a,aa),Rho,rho),z,zz),eps0,1/(36*pi)));
            disp(Ez1);
            clear x y z;
            x = -aa:0.01:aa;
            y = [sqrt(aa^2-x.^2);-sqrt(aa^2-x.^2)];
            z = zeros(1,size(x,2));
            if rho>0
                plot3(x,y,z,'-b','linewidth',3);
                hold on;
            else
                plot3(x,y,z,'-r','linewidth',3);
                hold on;
            end

            plot3(0,0,zz,'s','markersize',13,'markerfacecolor','y')
            plot3([0,0],[0,0],[zz,(Ez1/abs(Ez1))+zz],'linewidth',5)
            xlabel('X');ylabel('Y');zlabel('Z');
            grid on;


        end
end

