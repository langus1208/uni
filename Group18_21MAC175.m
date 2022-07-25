r_simplex_combined([7 0 11 -10 0 26],[1 -1 1 0 1 1; 0 1 -1 1 0 3; 1 1 -3 1 1 0; 0 0 1 0 0 1],[6 28 13 7]',['==','<=','>=','>='],'min',100000000000000)


function [z_min, z_max] = r_simplex_combined(c,A,b,V,Con,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Section 1 converts a problem to a maximisation problem and to standard
% form

%Solve a problem of the form: optimise z = c'x subject to Ax = b, x>=0

%%%Before running this function, define M as a large (comparative to other
%values of c and b) positive number e.g

%Example of input: (problem question) r_simplex_combined([7 0 11 -10 0 26],[1 -1 1 0 1 1; 0 1 -1 1 0 3; 1 1 -3 1 1 0; 0 0 1 0 0 1],[6 28 13 7]',['==','<=','>=','>='],'min',100000)

%Input variables:
% c:initial cost coefficients
% A:constraint matrix
% b:solutions/RHS vector
% V:signs of the constraint equations, put in in the order from top to bottom
% Con:maximise or minimise  
V2=V;

if Con == 'max'
    M = -M;
end
%transfer A to canonical form
% n_con:number of constraint equtions
[Row,Column]=size(A);
n_con=length(V)/2;
for y=1:n_con;
    if V(2*y-1:2*y)=='==';
        V(2*y-1:2*y)=['aa'];       
    elseif V(2*y-1:2*y)=='>=';
        V(2*y-1:2*y)=['aa'];
        A=[A,zeros(Row,1)];
        [Row,Column]=size(A);
         A(y,Column)=-1;
         c=[c,0];
    else V(2*y-1:2*y)=='<=';
        V(2*y-1:2*y)=['aa'];
    end
end

[Row,Column]=size(A);
for y=1:n_con;
    if V2(2*y-1:2*y)=='==';
        V2(2*y-1:2*y)=['aa'];
        A=[A,zeros(Row,1)];
        [Row,Column]=size(A);
        A(y,Column)=1;
        c=[c,M];
    elseif V2(2*y-1:2*y)=='>=';
        V2(2*y-1:2*y)=['aa'];
        A=[A,zeros(Row,1)];
        [Row,Column]=size(A);
         A(y,Column)=1;
         c=[c,M];
    else V2(2*y-1:2*y)=='<=';
        V2(2*y-1:2*y)=['aa'];
        A=[A,zeros(Row,1)];
        [Row,Column]=size(A);
        A(y,Column)=1;
        c=[c,0];
    end
end

%change the cost coefficient to maximisation problem for the following code
if Con=='min'
    c=-c;
else Con=='max'
    c=c;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format rat

%generate BFS;
[m,n] = size(A);
bfs = zeros(1,n-m);
bfs = [bfs b'];
A;
b
c;
bfs;

%indices of A corresponding to the Basis matrix
c_b_indices = find(bfs);
%equivalently c_b_indices - indices of c corresponding to Cb'

%indices not corrensponding to the Basis matrix
V_indices = find(ones(1,n)-abs(sign(bfs)));
%generate this for later when the optimal is reached
solution = zeros(n,1);

%Max iteration prevents cycling
iters_max = 500;
for i = 1:iters_max

    %STEP 1
    %Compute inverse of B
    B = (A(:,c_b_indices));
    Binv = inv(B);

    %STEP 2
    %find invB * b
    d = Binv * b;

    %STEP 3
    %Optimality computation
    Opt_com = zeros(1,n);
    %This is Cb'*B^-1*
    Opt_com(:,V_indices) = c(:,c_b_indices)*Binv*A(:,V_indices)-c(:,V_indices);
    
    %STEP 4
    %Choosing entering vector
    %cj is most negative value
    %j is the column corrensponding to the entering vector
    [cj,j] = min(Opt_com);

    %STEP 5
    %if most negative value (cj) is >= 0 then solution is optimal already
    if cj >= 0
        solution(c_b_indices,:) = d;
        d
        z_max = c(:,c_b_indices)*d
        z_min = -z_max
        %isnan here is used to detect division by 0; has to be done at the
        %solution as 0 is an approximation after iterations.
        if (isnan(z_min) == 1)
            fprintf('<strong>System is unbounded</strong>\n')
        end
        %z_max > M shows that the artificial variables has not been forced
        %to 0 - this requires an arbitrarily large M to be defined by the
        %user.
        %And clearly z_min > z_max is an infeasible solution
        if (z_min > z_max)
            fprintf('<strong>System is infeasible</strong>\n')
        elseif (z_max > M)
            fprintf('<strong>System is infeasible</strong>\n')
        end
        return;
    end
    
    %STEP 6
    %Calculating the leaving vector   

    %A(:,j) is Pj, the entering vector
    w = Binv * A(:,j);
    %Finding minimum positive of d(i)/w(i)
    min_com = find(w>=0)';
    [min_val,val_pos] = min(d(min_com) ./ w(min_com));
    k = min_com(val_pos(1));

    %Step 7
    %Swapping the leaving vector out of B with the entering vector
    leaving_col = c_b_indices(k);
    %Swap entering vector into indices of c corresponding Cb'
    c_b_indices(k) = j;
    %Swap leaving vector into vector of indices of A not corresponding to the
    %basis matrix
    V_indices(j == V_indices) = leaving_col;
    
    if i==iters_max
        fprintf("System is degenerate")
    end
end
end
