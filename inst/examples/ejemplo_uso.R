# Ejemplo completo de uso del paquete hexeInvMC

# Cargar el paquete
library(hexeInvMC)

# Datos de entrada
x <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
y <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)

# Valor observado
y0  <- 10.5
uy0 <- 0.25

# Selección de distribuciones
dist_x  <- "triangle"
dist_y  <- "norm"
dist_y0 <- "unif"

# Ejecutar regresión inversa con simulación Monte Carlo extendida
# Número de simulaciones ajustable (por defecto 10000)
resultado <- hexeInvMC(
  x, ux, y, uy,
  y0 = y0, uy0 = uy0,
  dist_x = dist_x,
  dist_y = dist_y,
  dist_y0 = dist_y0,
  n_sim = 10000
)

# Imprimir resultado
print(resultado)
