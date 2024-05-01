function Mat = vectorSkewSymMat(vec)
% Compute skew symmetric matrix representation of a vector cross product 
% [vec x]

Mat = [ 0     -vec(3)   vec(2);
      vec(3)    0    -vec(1);
     -vec(2)   vec(1)    0 ];