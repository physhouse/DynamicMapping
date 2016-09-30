function pdf2d(X, Y)

xlo = min(X);
xhi = max(X);
ylo = min(Y);
yhi = max(Y);

difx = (xhi - xlo) / 100;
dify = (yhi - ylo) / 100;

x_axis = xlo:difx:xhi;
y_axis = ylo:dify:yhi;

[x_mesh_upper, y_mesh_upper] = meshgrid(x_axis(2:end), y_axis(2:end));
[x_mesh_lower, y_mesh_lower] = meshgrid(x_axis(1:end-1), y_axis(1:end-1));

x_centers = (x_axis(2:end) + x_axis(1:end-1)) / 2;
y_centers = (y_axis(2:end) + y_axis(1:end-1)) / 2;

pdf = mean( bsxfun(@le, X(:), x_mesh_upper(:).') ...
    & bsxfun(@gt, X(:), x_mesh_lower(:).') ...
    & bsxfun(@le, Y(:), y_mesh_upper(:).') ...
    & bsxfun(@gt, Y(:), y_mesh_lower(:).') );

pdf = reshape(pdf, length(x_axis)-1, length(y_axis)-1);

pdf = pdf ./ (y_mesh_upper - y_mesh_lower) ./ (x_mesh_upper - x_mesh_lower);

figure
imagesc(x_centers, y_centers, pdf)
axis([xlo xhi ylo yhi])
%axis equal
colorbar