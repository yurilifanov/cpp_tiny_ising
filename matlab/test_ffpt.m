close all, clear all, fclose all;
kbt = 1.0;
h = 0.0;
n_smp = 1e6;

cd('..');
system(sprintf('./a.out %d %0.6f %0.6f', n_smp, kbt, h));
cd('matlab');

t = linspace(0, 500, 100);
g = get_sparse_generator(kbt, h, '../');
[f, r] = ffpt_and_rate(g, 1, size(g, 1), t);

t_smp = readbin('../fpts.bin', 'double');
[f_smp, x] = ecdf(t_smp);

plot(t, f, 'ok'), hold all,
plot(x, f_smp, '.')

function data = readbin(fname_str, type_str)
  f = fopen(fname_str, 'r');
  data = fread(f, type_str);
  fclose(f);
end