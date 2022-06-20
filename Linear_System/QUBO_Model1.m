% Copyright 2022 Jun's research group
% Licensed under the MIT License
% https://www.mit.edu/~amini/LICENSE.md

close all, clear all, clc;

Qubits = 2;
Selection = 2;   % 1 is random, 2 is manual

% qubits initiallization
vec_qb = zeros(Qubits,2^Qubits);   
for qb = 1:2^Qubits-1 %Qubits
    qb_st = qb;
    for k = 1:Qubits
        next = fix(qb_st/2);
        vec_qb(k, qb+1) = rem(qb_st, 2);
        qb_st = next;
    end
end

if Selection == 1
    Matrix_size = 4;
    Range = 100;
    Sub_int = 2^Qubits;     % subrange interval
    Sub_num = 30;           % number of subinterval for each x_i
    mat_M = Range*rand(Matrix_size)-Range/2;                            % A in Ax = b
    vec_x = Sub_int*Sub_num*rand(Matrix_size,1) - Sub_int*Sub_num/2;    % x
    vec_b = mat_M*vec_x;
else
    mat_M = [3, 1 ; -1, 2];
    vec_b = [-1; 5];
    Matrix_size = max(size(mat_M));
end

Matrix = mat_M
Vector_b = vec_b

T0 = zeros(Matrix_size,1);
vec_c = vec_b - mat_M*T0;

QM = zeros(2*Qubits*Matrix_size, 2*Qubits*Matrix_size);
%%% Linear terms %%%
for k = 1:Matrix_size
    for i = 1:Matrix_size
        for l = 1:Qubits
            cef1 = (2^(2*(l-1)))*(mat_M(k,i)^2);
            cef2 = (2^(l))*(mat_M(k,i)*vec_c(k));
            po1 = 2*Qubits*(i-1) + l;
            po2 = 2*Qubits*(i-1) + l + Qubits;
            QM(po1,po1) = QM(po1,po1) + cef1 - cef2;
            QM(po2,po2) = QM(po2,po2) + cef1 + cef2;
        end
    end
end
 
%%% First quadratic term %%% 
for k = 1:Matrix_size
    for i = 1:Matrix_size
        for l1 = 1:Qubits-1
            for l2 = l1+1:Qubits
                qcef = 2^(l1+l2-1)*mat_M(k,i)^2;
                po1 = 2*Qubits*(i-1) + l1;
                po2 = 2*Qubits*(i-1) + l2;
                QM(po1,po2) = QM(po1,po2) + qcef;
                po3 = 2*Qubits*(i-1) + l1 + Qubits;
                po4 = 2*Qubits*(i-1) + l2 + Qubits;
                QM(po3,po4) = QM(po3,po4) + qcef;
            end
        end
    end
end
                
%%% Second quadratic term %%% 
for k = 1:Matrix_size
    for i = 1:Matrix_size-1
        for j  = i+1:Matrix_size
            for l1 = 1:Qubits
                for l2 = 1:Qubits 
                    qcef = 2^(l1+l2-1)*mat_M(k,i)*mat_M(k,j); 
                    po1 = 2*Qubits*(i-1) + l1;
                    po2 = 2*Qubits*(j-1) + l2;
                    QM(po1,po2) = QM(po1,po2) + qcef;
                    po3 = 2*Qubits*(i-1) + l1 + Qubits;
                    po4 = 2*Qubits*(j-1) + l2 + Qubits;
                    QM(po3,po4) = QM(po3,po4) + qcef;
                    po5 = 2*Qubits*(i-1) + l1;
                    po6 = 2*Qubits*(j-1) + l2 + Qubits;
                    QM(po5,po6) = QM(po5,po6) - qcef;
                    po7 = 2*Qubits*(i-1) + l1 + Qubits;
                    po8 = 2*Qubits*(j-1) + l2;
                    QM(po7,po8) = QM(po7,po8) - qcef;
                end
            end
        end
    end
end

%QM
Abs_Min = -dot(vec_c, vec_c);
fprintf('Absolute minimum value from vector product -(b*)*(b) = %f\n\n', Abs_Min)
Annealing = 0;
vec_Anl = [];
  
for k = 0:2^(2*Qubits*Matrix_size)-1
    T = [];
    temp_po = k;
    for Order = 1:2*Matrix_size
        po = rem(temp_po, 2^Qubits);
        next_po = (temp_po - po)/2^Qubits;
        temp_po = next_po;
        T = [T; vec_qb(:,po+1)];
    end
    qubo_val = T.'*QM*T;
    if qubo_val < Annealing
        Annealing = qubo_val;
        vec_Anl = T;
    end
end
fprintf('Absolute minimum value from annealing = %f\n\n', Annealing)
fprintf('Vector for absolute minimum value = \n')
fprintf('q1 q2 q3 q4 q5 q6 q7 q8 ---\n')
fprintf('%d  ', vec_Anl)
sol_x = zeros(Matrix_size,1);
for k = 1:Matrix_size
    for i = 1:Qubits
        sol_x(k) = sol_x(k) + vec_Anl((k-1)*2*Qubits+i)*2^(i-1);
        sol_x(k) = sol_x(k) - vec_Anl((k-1)*2*Qubits+i+Qubits)*2^(i-1);
    end
end
fprintf('\n\nConverted vector from annealing = \n')
fprintf('%f  ', sol_x)

Ax_b = mat_M*sol_x - vec_b;
fprintf('\n\nAx - b = \n')
fprintf('%f  ', Ax_b)
fprintf('\n')

