set.seed(0)
mu = 0.01; L = 1; kappa = L/mu
n = 100
D = runif(n); D = 10^D; Dmin = min(D); Dmax = max(D)
D = (D-Dmin) / (Dmax-Dmin)
D = mu + D*(L-mu)
A = diag(D) 
x0 = runif(n, 0, 1)
x_star = rep(0, 100)

f <- function(x) {
  0.5*t(x) %*% A %*% x
}

df <- function(x) {
  A %*% x
}

GradientDescent <- function(x0, x_star, L, f, df, e = 1e-6) {
  iter = 0
  value = f(x0)
  x1 = x0
  while (f(x1) - f(x_star) > e) {
    x1 = x1 - (1/L)*df(x1)
    iter = iter + 1
    value = c(value, f(x1))
  }
  return(list(x0, iter, value))
}

Newton <- function(x0, x, f, df, e = 1e-6) {
  iter = 0
  value = f(x0)
  x1 = x0
  while (f(x1) - f(x_star) > e) {
    x1 = x1 - solve(A)%*%df(x1)
    iter = iter + 1
    value = c(value, f(x1))
  }
  return(list(x0, iter, value))
}

Nesterov <- function(x0, x_star, L, m, f, df, e = 1e-6) {
  iter = 0
  value = f(x0)
  alpha = 1/L
  beta = (sqrt(L) - sqrt(m))/(sqrt(L)+sqrt(m))
  x1 = x0
  while (f(x1) - f(x_star) > e) {
    # Write your updates here:
    y1 = x1 + beta*(x1 - x0)
    x0 = x1
    x1 = y1 - alpha*df(y1)
    iter = iter + 1
    value = c(value, f(x1))
  }
  return(list(x0, iter, value))
}

out_GD_1 <- GradientDescent(x0, x_star, L = 1, f, df)
out_GD_2 <- GradientDescent(x0, x_star, L = 5, f, df)
out_Nest <- Nesterov(x0, x_star, L = 1, m = 0.01, f, df)
out_Newton <- Newton(x0, x_star, f, df)
# Plots:
plot(0, type="n", xlab="number of iterations", ylab="log of function differences", xlim=c(0, 500), ylim=c(-7, 1))
lines(log(out_GD_1[[3]], 10))
lines(log(out_GD_2[[3]], 10), col = 'blue')
lines(log(out_Nest[[3]], 10), col = 'red')
lines(log(out_Newton[[3]], 10), col = 'purple')
legend("topright", legend = c("GD1", "GD2", "Nest", "Newton"), lty = rep(1,4), col = c("black", "blue", "red", "purple"))




