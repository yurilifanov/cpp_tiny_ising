import numpy as np
from plotly import tools
from plotly.offline import plot
from plotly.graph_objs import Scatter, Layout, Figure
from scipy.sparse import csr_matrix, hstack
from scipy.sparse.linalg import lsqr, expm_multiply

def get_e_ffpt(r_dir=''):
  fpts = np.fromfile(r_dir + './fpts.bin', dtype=np.float64)
  h, bins = np.histogram(fpts, bins='auto', density=False)
  h = h / np.sum(h)
  return np.cumsum(h), bins

def get_ffpt(Q, A, I, v, t_min, t_max, num_t_pts):
  R = np.hstack((A, I))
  Q_RR = Q[R, :][:, R]
  v_R = v[R]
  return 1. - np.sum(expm_multiply(Q_RR.T, v_R.T, t_min, t_max, num_t_pts), 1)

def get_sparse_generator(kbt, h, r_dir='', lazy=False):
  i_from = np.fromfile(r_dir + './cfgs_from.bin', dtype=np.uint64)
  i_to = np.fromfile(r_dir + './cfgs_to.bin', dtype=np.uint64)
  d_ones = np.fromfile(r_dir + './ones_deltas.bin', dtype=np.int64)
  d_aligned = np.fromfile(r_dir + './aligned_deltas.bin', dtype=np.int64)

  p_acc = np.minimum(1., np.exp((d_aligned + h * d_ones) / kbt))
  accum = np.bincount(i_from.astype(np.int64), weights=p_acc)
  n_cfg = len(accum)
  N = np.log2(n_cfg)

  if lazy:
    #note that there is no division by N of the acc probabilities,
    #so no need to multiply by N here 
    rates = np.hstack((p_acc, -accum))
  else:
    rates = np.hstack((p_acc / accum[i_from], -np.ones(n_cfg))) * N

  i_from = np.hstack((i_from, np.arange(n_cfg)))
  i_to = np.hstack((i_to.astype(np.int64), np.arange(n_cfg)))
  return csr_matrix((rates, (i_from, i_to)), shape=(n_cfg, n_cfg))

def get_stationary_distribution(sparse_Q):
  n = sparse_Q.shape[0]
  Y = np.hstack((np.zeros(n), 1))
  M = hstack((sparse_Q, np.ones((n, 1))))
  ans = lsqr(M.T, Y, atol=1e-15, btol=1e-15)
  print('Residuals: {}'.format(ans[3]))
  return ans[0]

if __name__ == '__main__':
  Q = get_sparse_generator(1, 0, lazy=False)
  v = get_stationary_distribution(Q)
  A = 0
  n_cfg = Q.shape[0]
  B = n_cfg - 1
  I = np.setdiff1d(np.arange(n_cfg), np.hstack((A, B)))

  v[np.hstack((I, B))] = 0
  v = v / np.sum(v)

  t_min, t_max, num_pts = (10, 40, 10)
  ffpt = get_ffpt(Q, A, I, v, t_min, t_max, num_pts)
  t = np.linspace(t_min, t_max, num_pts)

  e_ffpt, e_t = get_e_ffpt()
  f_trace = Scatter(x=t, y=ffpt, name='Analytical', mode='markers')
  e_trace = Scatter(x=e_t, y=e_ffpt, name='Empirical', mode='line')
  layout = Layout(title='First passage time CDF',
                  xaxis={'title': 't (Survival time)'},
                  yaxis={'title': 'CDF(t)'})
  figure = Figure(data=[f_trace, e_trace], layout=layout)
  plot(figure, filename='ffpt_plot.html')

