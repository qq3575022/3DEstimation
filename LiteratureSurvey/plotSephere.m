[X,Y,Z] = sphere;

figure
surface(1.2*X,1.2*Y,1.2*Z,'FaceColor', [0.3 0.3 1],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor',[0.2,0.2,0.2])
axis equal

hold on
r = 0.8;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

surf(1.2*X2+0.5,1.2*Y2-0.8,Z2-0.7,'FaceColor', [1 0.3 0.3],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor',[0.2,0.2,0.2])

hold on
r = 0.8;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

surf(1.2*X2-0.5,1.2*Y2-0.8,Z2-0.7,'FaceColor', [0.3 1 0.3],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor',[0.2,0.2,0.2])

hold on
r = 0.8;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

surf(1.2*X2+1.0,1.2*Y2+0.3,Z2-0.7,'FaceColor', [1 0.7 0.3],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor',[0.2,0.2,0.2])


xlabel('X')
ylabel('Y')
zlabel('Z')
camlight