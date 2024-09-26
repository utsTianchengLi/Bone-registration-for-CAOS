function draw= Draw(V1, f1,color,aaa,light)

patch('Faces', f1, 'Vertices', V1,'FaceColor', color, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.15);
if ~isempty(aaa)
alpha(0.5); 
end
alpha(0.9); 
hold on;
axis equal;
% Add a camera light, and tone down the specular highlighting
if light==1
% camlight('left');
camlight(90,65)
end
material('dull');
% xlim([-100 150]);
view(160,20)
% title('bone model (world frame)');
xlabel('X[mm]');ylabel('Y[mm]');zlabel('Z[mm]');
end