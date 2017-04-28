import numpy as np
from scipy.sparse import csr_matrix, hstack
from scipy.sparse.linalg import lsqr, expm_multiply

def get_ffpt(Q, A, I, v, t_min, t_max, num_t_pts):
  R = np.hstack((A, I))
  Q_RR = Q[R, :][:, R]
  v_R = v[R]
  return 1. - np.sum(expm_multiply(Q_RR.T, v_R.T, t_min, t_max, num_t_pts), 1)

def get_sparse_generator(kbt, h, r_dir=''):
  i_from = np.fromfile(r_dir + './cfgs_from.bin', dtype=np.uint64)
  i_to = np.fromfile(r_dir + './cfgs_to.bin', dtype=np.uint64)
  d_ones = np.fromfile(r_dir + './ones_deltas.bin', dtype=np.int64)
  d_aligned = np.fromfile(r_dir + './aligned_deltas.bin', dtype=np.int64)

  p_gen_acc = np.minimum(1., np.exp((d_aligned + h * d_ones) / kbt))
  accum = np.bincount(i_from.astype(np.int64), weights=p_gen_acc)
  n_cfg = len(accum)
  n_sites = np.log2(n_cfg)

  p_gen_acc = np.hstack((p_gen_acc / accum[i_from], -np.ones(n_cfg))) * n_sites
  i_from = np.hstack((i_from, np.arange(n_cfg)))
  i_to = np.hstack((i_to.astype(np.int64), np.arange(n_cfg)))
  return csr_matrix((p_gen_acc, (i_from, i_to)), shape=(n_cfg, n_cfg))

def get_stationary_distribution(sparse_Q):
  n = sparse_Q.shape[0]
  Y = np.hstack((np.zeros(n), 1))
  M = hstack((sparse_Q, np.ones((n, 1))))
  ans = lsqr(M.T, Y, atol=1e-15, btol=1e-15)
  print('Residuals: {}'.format(ans[3]))
  return ans[0]

if __name__ == '__main__':
  Q = get_sparse_generator(1, 0)
  v = get_stationary_distribution(Q)
  A = 0
  n_cfg = Q.shape[0]
  B = n_cfg - 1
  I = np.setdiff1d(np.arange(n_cfg), np.hstack((A, B)))
  #ffpt = get_ffpt(Q, A, I, v, 0, 1000, 10)
  ffpt = np.array([
              0.07878145,  
              0.99236894,  
              0.99989866,  
              0.99999865,  
              0.99999998,
              1.        ,  
              1.        ,  
              1.        ,  
              1.        ,  
              1.        ])
