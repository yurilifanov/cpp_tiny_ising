function mat = get_sparse_generator(kbt, h, dir_str)
  i_from = readbin([dir_str, '/cfgs_from.bin'], 'uint64') + 1;
  i_to = readbin([dir_str, '/cfgs_to.bin'], 'uint64') + 1;
  d_ones = readbin([dir_str, '/ones_deltas.bin'], 'int64');
  d_aligned = readbin([dir_str, '/aligned_deltas.bin'], 'int64');

  p_gen_acc = min(1, exp((d_aligned + h * d_ones) / kbt));
  accum = accumarray(i_from, p_gen_acc);
  num_cfg = length(accum);
  num_sites = log2(num_cfg);

  p_gen_acc = [p_gen_acc ./ accum(i_from); -ones(num_cfg, 1)] ./ num_sites;
  i_from = [i_from; (1:num_cfg)'];
  i_to = [i_to; (1:num_cfg)'];

  mat = sparse(i_from, i_to, p_gen_acc);
end
function data = readbin(fname_str, type_str)
  f = fopen(fname_str, 'r');
  data = fread(f, type_str);
  fclose(f);
end


