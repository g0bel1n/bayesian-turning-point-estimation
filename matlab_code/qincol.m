function Qcol = qincol(Q)

% wandelt eine quadrat Matrix der dim dd in eine Spalte der dim dd*(dd+1)/2 um
% wobei nur die obige Dreiecksmatrix spaltenweise hintereinander geschrieben wird

% input: Q .. ddxdd
% output: Qcol .. dd*(dd+1)/2x1

dd = size(Q,2);
hilf = 0;
for i = 1:dd,
   Qcol(hilf+1:hilf+i,1) = Q(1:i,i);
   hilf =hilf + i;
end   