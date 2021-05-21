function [F,X,Y,P] = SphereDipPotential(XYZ,Q,D,R,r0,a,b,Dx,Dy,Nxy)

%about user's plane
e1 = a/norm(a); 
e2 = b - dot(e1, b).*e1; 
e2 = e2/norm(e2);%basis vectors in plane
e3 = cross(e1, e2)/norm(cross(e1, e2));
P = [e1, e2, e3]; %transition matrix from global 3d basis to basis in plane

xy = P\(XYZ - r0);% coordinates of charges in 2d basis
P = [e1, e2]; % matrix of mapping from plane to 3d 

X = zeros(Nxy(1), Nxy(2));
Y = zeros(Nxy(1), Nxy(2));
for ii = 1 : Nxy(1) % y-coordinates of points of plane
Y(ii, :) = Dy(1) : (Dy(2) - Dy(1))/(Nxy(1) - 1) : Dy(2); 
end
for jj = 1 : Nxy(2) % x-coordinates of points of plane
X(:, jj) = Dx(1) : (Dx(2) - Dx(1))/(Nxy(2) - 1) : Dx(2); 
end

F = zeros(Nxy(1), Nxy(2)); %potentials of points in plane

% finds distance in 2d basis between point with coordinates (ii, jj) and
% (qq)th ball
function [vec] = plane_vec(ii, jj, qq) 
    vec = [X(ii, jj); Y(ii, jj); 0] - xy(:, qq); 
end

for ii = 1 : Nxy(1)
    for jj = 1 : Nxy(2) %running through all points
        for qq = 1 : size(Q, 1) %running through all charges
            if norm(plane_vec(ii, jj, qq)) < R(qq)
                F(ii, jj) = F(ii, jj) + Q(qq)/R(qq) + dot(D(qq,:), plane_vec(ii, jj, qq))/R(qq)^2; %if point is in ball
            else
                F(ii, jj) = F(ii, jj) + Q(qq)/norm(plane_vec(ii, jj, qq)) + dot(D(qq,:),plane_vec(ii, jj, qq))/norm((plane_vec(ii, jj, qq)))^3;               
            end
        end
    end
end

end

