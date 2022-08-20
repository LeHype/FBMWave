function [Value varargout] = HamiltonianPathForSimplex(Simplex,p)
%% Returns the Value of the Simplex at the Position of the Path
% |O---->---->--|
% |xxxxxxxxxxxxx|
% ||----<----<--|
% |xxxxxxxxxxxxx|
% ||---->---->--|
% |xxxxxxxxxxxxx|
% ||----<----<--|
Pathlength = size(Simplex,2)*size(Simplex,1);
p = mod(p,Pathlength);

    
if p ==0
    p = Pathlength;
end
Row = ceil(p/(size(Simplex,1)));
Row = (Row*2)-1;
if mod(p,size(Simplex,1)) ==0
    Row = Row+1;
    disp('1');
end


Column = mod(p,(size(Simplex,1)));
if Column == 0
    [a b Column] = HamiltonianPathForSimplex(Simplex,p-1);
disp('2')
else 
    Column = Column+1;
end


if mod(Row,3)==0 
    disp('3123123')
    [a b Column] = HamiltonianPathForSimplex(Simplex,p-1);
end
if p > (Pathlength -size(Simplex,1));
Column = 1;
Row = size(Simplex,1)-(p-(Pathlength -size(Simplex,1)))+1;
end
if p ==1
    Row = 1;
    Column =2;
end
Simplex(Row,Column) = 0
Value = Simplex(Row,Column);
varargout =cell(2);
varargout{1}=Row;
varargout{2} =Column;
end