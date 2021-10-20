function dist = distance_point_line(x,y,a,b,c)

dist = abs(a*x+b*y+c)/sqrt(a^2 + b^2);

end