# A small physics simulation of masses randomly connected by springs, moving around under the influence of Hooke's law.
# https://mastodon.social/@safest_integer/115128042642307128

k = 1 # spring constant
m = 1 # mass
d0 = 1 # rest length
n = 20 # n of particles
# init random positions in nx2 matrix 
xx = matrix(runif(2*n), n, 2)
# init velocities (all zero)
vv = 0 * xx

# random connections, stored in list cc where cc[[i]] is vector of indices to
# which mass i connects
cc = replicate(n, numeric(0))
while (any(sapply(cc, length) < 2)) {
  i = sample(n,1)
  j = sample((1:n)[-i],1)
  cc[[i]] = unique(c(cc[[i]], j))
  cc[[j]] = unique(c(cc[[j]], i))
}

# Plotting function, showing masses as points and springs between them as lines

plot_them = function() {
  par(mar=rep(0,4))
  plot(xx, xlim=c(-1,2), ylim=c(-1,2), pch=16, cex=2, ax=F, ann=F)
  for (i in 1:n) {
    x0 = xx[i,1]
    y0 = xx[i,2]
    x1 = xx[cc[[i]], 1]
    y1 = xx[cc[[i]], 2]
    segments(x0, y0, x1, y1)
  }
  invisible(NULL)
}
plot_them()

# Forces via Hooke's law F[i] = sum -k * (dij - d0) * eij (dij is distance and eij is unit vector from i to j)
# Acceleration via F[i] = m * a[i]

calc_accel = function() {
  aa = matrix(0,n,2)
  for (i in 1:n) {
    # distance vectors between i and neighbors 
    xn = xx[cc[[i]], ]
    vd = t(xx[i,] - t(xn))
    # absolute distances
    dd = sqrt(rowSums(vd^2)) 
    # total force vector on mass i
    ff = colSums(-k * (dd - d0) * vd / dd)
    # accelerate
    aa[i, ] = ff / m
  }
  aa
}

Main simulation loop: Calculate acceleration and update velocity and positions using forward Euler scheme

while(1) {
  aa = calc_accel()
  # euler updates of velocity and position
  vv = vv + .01 * aa
  xx = xx + .01 * vv
  # plot and use flush/hold for smoother animation
  plot_them()
  dev.flush(); dev.hold(); Sys.sleep(.01)
}

