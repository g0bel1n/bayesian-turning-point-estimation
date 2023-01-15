function Q = qinmatr(Qcol)

% fuer symm. Matrizen Umkehrung zu Programm qincol

% wandelt eine Spalte der dim dd*(dd+1)/2 eine quadrat Matrix der dim dd in eine quadrat Matrix der dim dd um
% wobei die Elemente der Spalte in die obige Dreiecksmatrix spaltenweise nebeneinander geschrieben werden
% und die Matrix zu einer symmetr. Matrix ergaenzt wird

% input: Qcol .. [dd*(dd+1)/2]x1
% output: Q .. ddxdd

D = size(Qcol,1);

dd = -.5 + (.25 + 2*D)^.5;

hilf = 0;
for i = 1:dd,
   Q(1:i,i) = Qcol(hilf+1:hilf+i,1);
   Q(i,1:i) = Qcol(hilf+1:hilf+i,1)';
   hilf =hilf + i;
end   