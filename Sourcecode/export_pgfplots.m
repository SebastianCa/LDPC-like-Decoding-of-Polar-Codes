function [  ] = export_pgfplots( x, y, filename)

filename=[filename '.dat'];

if (length(x)~=length(y))
    error('x and y do not have the same length!');
end

data = [ x(:) y(:) ];
save(filename, '-ascii','data');

end

