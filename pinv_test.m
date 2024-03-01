for i = 1:50
  A = randn(1000);

  B = inv(A);
end

for i = 1:50
  A = randn(1000);

  B = pinv(A);
end