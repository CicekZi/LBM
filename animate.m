%% animate


vorticity = ((circshift(V, [0, -1]) - circshift(V, [0, 1])) ./ (2.*dx)) - ((circshift(U, [-1, 0]) - circshift(U, [1, 0]))./ (2.*dx));

Vorticity= reshape(vorticity, [nx, ny]);

if mod(loop_counter,1)==0

  u = (sqrt(U.^2+V.^2));
  quiver(Vorticity',Vorticity', 'Color', 'r', 'LineWidth',2) 
  % imagesc(Vorticity');

  % imagesc(u');
  axis equal off;
  drawnow
  fprintf(' Iteration = %i, Uresidual = %f, Vresidual = %f \n', loop_counter, Uresidual, Vresidual)
end
