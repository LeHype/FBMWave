function [Value] = HamiltonianPathForSimplex(Simplex,p)
%% Returns the Value of the Simplex at the Position of the Path
% |O---->---->--|
% ||----<----<--|
% ||---->---->--|
% ||----<----<--|
Pathlength = size(Simplex,2)*size(Simplex,1);
p = mod(p,Pathlength);
if p ==0
    p = Pathlength;
end
Row = ceil(p/(size(Simplex,1)-1));
Column = mod(p,(size(Simplex,1)-1));
if Column == 0
    Column = size(Simplex,2);
else 
    Column = Column+1;
end


if mod(Row,2)==0
    Column = size(Simplex,1)-Column+2;
end
if p > (Pathlength -size(Simplex,1));
Column = 1;
Row = size(Simplex,1)-(p-(Pathlength -size(Simplex,1)))+1;
end
% Simplex(Row,Column) = 0;
Value = Simplex(Row,Column);
end