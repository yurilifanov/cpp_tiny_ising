function [ffpt, jab_Q] = ffpt_and_rate(Q, A, B, t)
  n = size(Q, 1);
  I = setdiff(1:n, [A, B]);

  pi_ = get_stationary_distribution(Q);
  v = pi_;
  v([I, B]) = 0;
  v = v ./ sum(v);

  ffpt = get_ffpt(Q, A, I, v, t);
  jab_Q = get_rate_via_generator(Q, A, B, v);
end
function pi_ = get_stationary_distribution(Q)
  n = size(Q, 1);
  Y = [zeros(1, n), 1];
  M = [Q, ones(n, 1)];
  pi_ = Y / M;
end
function ffpt = get_ffpt(Q, A, I, v, t)
  Q_RR = Q([A, I], [A, I]);
  v_R = v([A, I]);
  ffpt = zeros(size(t));
  for i = 1:length(t)
    ffpt(i) = 1 - sum(v_R * expm(t(i) * Q_RR));
  end
end
function jab = get_rate_via_generator(Q, A, B, v)
  n = size(Q, 1);
  I = setdiff(1:n, [A, B]);
  Q_prime = Q;
  Q_prime(:, A) = 0;
  Q_prime(B, :) = 0;
  Q_prime_RR = Q_prime([A, I], [A, I]);
  Q_prime_RB = Q_prime([A, I], B);
  v_R = v([A, I]);
  jab = sum(v_R * expm(1e10 * Q_prime_RR) * Q_prime_RB);
end
function jab = get_rate_via_transition_matrix(Q, A, B, v)
  n = size(Q, 1);
  I = setdiff(1:n, [A, B]);
  W = eye(n) - Q ./ diag(Q);
end