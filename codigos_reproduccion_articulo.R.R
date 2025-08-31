# Instalar paquetes necesarios (si no están instalados)
install.packages(c("HDInterval", "binom", "boot"))

# Cargar librerías
library(HDInterval)
library(binom)
library(boot)

# 1. INTERVALOS PARA PROPORCIÓN BINOMIAL (n=100, k=30, α=0.05)
n <- 100
k <- 30
alpha <- 0.05

# a) Clopper-Pearson (exacto)
ci.binom <- function(n, k, alpha) {
  if (k == 0) {
    p1 <- 0.0
    p2 <- 1 - (alpha/2)**(1/n)
  } else if (k == n) {
    p1 <- (alpha/2)**(1/n)
    p2 <- 1.0
  } else {
    helper <- function(p, k, n, val) {
      return(pbinom(k, n, p) - val)
    }
    r <- uniroot(helper, k = (k-1), n = n, val = 1-alpha/2, interval = c(0,1))
    p1 <- r$root
    r <- uniroot(helper, k = k, n = n, val = alpha/2, interval = c(0,1))
    p2 <- r$root
  }
  return(data.frame(p1 = p1, p2 = p2))
}

clopper_pearson <- ci.binom(n, k, alpha);clopper_pearson

# b) Wilson (aproximación normal)
wilson_interval <- binom.confint(k, n, methods = "wilson");wilson_interval

# c) HPD (bayesiano)
hdi_interval <- hdi(qbeta, 1-alpha, shape1 = (k+1), shape2 = (n-k+1));hdi_interval

# d) Razón de verosimilitud (K=8)
lr.binom <- function(n, k, K) {
  helper <- function(p, n, k, K) {
    return((p**k * (1-p)**(n-k)) / ((k/n)**k * (1-k/n)**(n-k)) - 1/K)
  }
  p1 <- rep(0, length(k))
  p2 <- p1
  if (k == 0) {
    p1 <- 0
  } else {
    r <- uniroot(helper, n = n, k = k, K = K, interval = c(0, k/n))
    p1 <- r$root
  }
  if (k == n) {
    p2 <- 1
  } else {
    r <- uniroot(helper, n = n, k = k, K = K, interval = c(k/n, 1))
    p2 <- r$root
  }
  return(data.frame(p1 = p1, p2 = p2))
}

lr_interval <- lr.binom(n, k, 8);lr_interval

# 2. FIGURA 3: RAZÓN DE VEROSIMILITUD (n=40, k=10)
n_fig <- 40
k_fig <- 10
p_hat <- k_fig/n_fig

# Función de razón de verosimilitud
likelihood_ratio <- function(p) {
  (p^k_fig * (1-p)^(n_fig-k_fig)) / (p_hat^k_fig * (1-p_hat)^(n_fig-k_fig))
}

# Crear secuencia de valores p
p_values <- seq(0.05, 0.45, length.out = 1000)
lr_values <- sapply(p_values, likelihood_ratio)

# Encontrar límites del intervalo (K=8)
roots <- lr.binom(n_fig, k_fig, 8)

# Graficar (versión básica)
plot(p_values, lr_values, type = "l", lwd = 2, col = "blue",
     main = "Figura 3: Razón de Verosimilitud L(p)/L(̂p) para n=40, k=10",
     xlab = "p", ylab = "L(p)/L(̂p)")
abline(h = 1/8, lty = 2, col = "red")
abline(v = roots$p1, lty = 3, col = "green")
abline(v = roots$p2, lty = 3, col = "green")
abline(v = p_hat, lty = 1, col = "black")

# 3. INTERVALOS PARA LA MEDIA (distribución asimétrica)
# Generar datos según f(x)=3x²
set.seed(123)
n_media <- 20
x <- runif(n_media)^(1/3)  # Transformación inversa para f(x)=3x²

# Intervalo t clásico
t_interval <- t.test(x, conf.level = 0.95)$conf.int;t_interval

# Bootstrap BCa
boot_mean <- function(data, indices) {
  return(mean(data[indices]))
}
boot_results <- boot(x, boot_mean, R = 1000);boot_results
bootstrap_ci <- boot.ci(boot_results, type = "bca");bootstrap_ci

