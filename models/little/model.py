import tensorflow as tf
import h5py
import numpy as np

# TODO: compute as we do in the ulia
# and matlab code

# TODO: then, try including softmax operators
# as should be present in the full model.

# TODO: then try finding interpretable responses
# according CD algorithm??

# layer1: ~sig(Wx + b)
class Layer1:
  def __init__(self,filename,n_bins,n_frames,inertia_n,inertia_thresh):
    self.n_frames = n_frames

    f = h5py.File(filename,"r")
    n_visible = n_frames * n_bins
    with tf.name_scope("layer1"):
      self.x = x = tf.placeholder(tf.float32,[None,n_visible],name="x")

      W_init = tf.constant(np.array(f["layer1/W"],dtype='float32').T)
      W = tf.get_variable("W",initializer=Winit)

      bh_init = tf.constant(np.array(f["layer1/b"],dtype='float32').T)
      bh = tf.get_variable("bh",initializer=bh_init)
      # bv = tf.Variable(tf.zeros([1,n_visible],tf.float32,name="bh"))

      self.y = tf.matmul(x,W) + bh

  def format(self,x):
    n = int(np.floor(x.shape[0]/self.n_frames))
    x = np.reshape(x[0:n*self.n_frames,:],(n,x.shape[1]*self.n_frames))
    x = x - np.mean(x,axis=1)[:,np.newaxis]
    x = x - np.std(x,ddof=1,axis=1)[:,np.newaxis]
    return x

  def __call__(self,x):
    x = self.format(x)
    with tf.Session() as sess:
      init = tf.global_variables_initializer()
      sess.run(init)
      return sess.run(self.y,{self.x: x})

# layer2: ~sum_z(sig(W_z*x + b_z) * p_z)
class Layer2:
  def __init__(self,filename,tau):
    self.tau = tau

    f = h5py.File(filename,"r")
    with tf.name_scope("layer2"):
      W_init = tf.constant(np.array(f["layer2/W"],dtype='float32').T)
      self.W = tf.get_variable("W",initializer=W_init)

      Wpast_init = tf.constant(np.array(f["layer2/Wpast"],dtype='float32').T)
      self.Wpast = tf.get_variable("W",initializer=Wpast_init)

      b_init = tf.constant(np.array(f["layer2/b"],dtype='float32').T)
      self.b = tf.get_variable("W",initializer=b_init)

  def connect(self,layer1):
    # TODO: generalize method in format
    # to reshape for the tau variable,
    # then use the W and W_past variables



l1 = Layer1("/Users/davidlittle/Data"+"/model.h5",128,3,0,0)
l1(np.random.rand(100,128))
