%Limpieza de pantalla
clear all
close all
clc
% Declaración de las variables simbólicas
syms q1(t) l1 t q2(t) l2 q3(t) l3
% Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP = [0 0 0];
% Creamos el vector de coordenadas articulares
Q = [q1 q2 q3];
disp("Coordenadas articulares"); pretty(Q);
% Creamos el vector de velocidades articulares
Qp = diff(Q, t);
disp("Velocidades articulares "); pretty(Qp);
% Número de grado de libertad del robot
GDL = size(RP, 2);
GDL_str = num2str(GDL);

% Articulación 1
% Definición de la posición de las articulaciones con respecto a la
% anterior
P(:,:,1) = [l1*cos(q1); 
            l1*sin(q1); 0]; 
P(:,:,2) = [l2*cos(q2);
            l2*sin(q2);0]; 
P(:,:,3) = [l3*cos(q3);
            l3*sin(q3);0];
% Matriz de rotación de la articulaciones actuales con respecto a la
% anterior
R(:,:,1) = [cos(q1) -sin(q1) 0;
            sin(q1) cos(q1)  0;
            0          0     1];
R(:,:,2) = [cos(q2) -sin(q2) 0;
            sin(q2)  cos(q2) 0;
            0        0       1];
R(:,:,3) = [cos(q3) -sin(q3) 0;
            sin(q3)  cos(q3) 0;
            0        0       1];
% Creamos un vector de ceros
Vector_zeros = zeros(1, 3);

% Inicializamos las matrices de transformación homogéneas locales
A(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_zeros 1]);

% Inicializamos las matrices de transformación homogénea locales
T(:,:,GDL) = simplify([R(:,:,GDL) P(:,:,GDL); Vector_zeros 1]);

% Inicializamos los vectores de posición y rotacion
PO(:,:,GDL) = P(:,:,GDL);
RO(:,:,GDL) = R(:,:,GDL);

for i = 1:GDL
    i_str = num2str(i);
    % Locales 
    disp(strcat("Matriz de transformación local A", i_str));
    A(:,:,i) = simplify([R(:,:,i) P(:,:,i); Vector_zeros 1]);
    pretty(A(:,:,i));
    % Globales
    try
        T(:,:,i) = T(:,:,i-1)*A(:,:,i);
    catch
        T(:,:,i) = A(:,:,i);
    end
    disp(strcat("Matriz de transformación global T", i_str));
    T(:,:,i) = simplify(T(:,:,i));
    pretty(T(:,:,i));
    % Obtenemos la matriz de rotación "RO" y el vector de translación PO de
    % la matriz de transformación homogénea global T(:,:,GDL)
    RO(:,:,i) = T(1:3, 1:3, i);
    PO(:,:,i) = T(1:3, 4, i);
    pretty(RO(:,:,i));
    pretty(PO(:,:,i))
end

% Calculamos el jacobiano lineal de forma diferencial
% Derivadas parciales de x respecto a q1
Jv11 = functionalDerivative(PO(1, 1, GDL), q1);
Jv12 = functionalDerivative(PO(1, 1, GDL), q2);
Jv13 = functionalDerivative(PO(1, 1, GDL), q3);
% Derivadas parciales de y respecto a q1
Jv21 = functionalDerivative(PO(2, 1, GDL), q1);
Jv22 = functionalDerivative(PO(2, 1, GDL), q2);
Jv23 = functionalDerivative(PO(2, 1, GDL), q3);
% Derivadas parciales de z respecto a q1
Jv31 = functionalDerivative(PO(3, 1, GDL), q1);
Jv32 = functionalDerivative(PO(3, 1, GDL), q2);
Jv33 = functionalDerivative(PO(3, 1, GDL), q3);
% Creamos la matriz del jacobiano lineal
jv_d = simplify([Jv11 Jv12 Jv13; Jv21 Jv22 Jv23; Jv31 Jv32 Jv33]);
disp("Jacobiano lineal obtenido de forma diferencial"); pretty(jv_d)

%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);

for k= 1:GDL
    if RP(k)==0 
       %Para las juntas de revolución
        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));
            Jw_a(:,k)=[0,0,1];
         end
     else
%         %Para las juntas prismáticas
        try
            Jv_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)=[0,0,1];
        end
            Jw_a(:,k)=[0,0,0];
     end
 end    
Jv_a = simplify(Jv_a);
Jw_a = simplify(Jw_a);
disp("Jacobiano lineal obtenido de forma analítica"); pretty(Jv_a);
disp("Jacobiano angular obtenido de forma analítica");pretty(Jw_a);
V =  simplify(Jv_a*Qp');
disp("Velocidad lineal obtenida mediante el Jacobiano lineal"); pretty(V);
W = simplify(Jw_a*Qp');
disp("Velocidad angular mediante el Jacobiano angular"); pretty(W);
