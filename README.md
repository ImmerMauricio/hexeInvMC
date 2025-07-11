
# hexeInvMC 🧪

**Regresión inversa con incertidumbre en ambos ejes, basada en
simulación Monte Carlo extendida.**

Este paquete está diseñado para laboratorios reales que necesitan
estimar concentraciones a partir de señales (`y`), considerando
incertidumbres tanto en `x` como en `y`, con modelos lineales o
cuadráticos, y permitiendo configuraciones flexibles para las
distribuciones de incertidumbre.

------------------------------------------------------------------------

## 📦 Instalación

Puedes instalar la versión más reciente desde GitHub:

``` r
# Instalar 'devtools' si aún no lo tienes
install.packages("devtools")

# Instalar hexeInvMC desde GitHub
devtools::install_github("ImmerMauricio/hexeInvMC")
```

------------------------------------------------------------------------

## 🧪 Ejemplo de uso

``` r
library(hexeInvMC)

# Datos experimentales
x  <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
y  <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)

# Medida instrumental observada
y0  <- 10.5
uy0 <- 0.25

# Llamado a la función principal
hexeInvMC(
  x, ux, y, uy,
  y0 = y0, uy0 = uy0,
  dist_x = "norm",   # Puede ser "norm", "unif" o "triangle"
  dist_y = "norm",
  dist_y0 = "norm",
  n_sim = 10000      # Número de simulaciones Monte Carlo
)
```

------------------------------------------------------------------------

## ⚙️ Parámetros disponibles

| Argumento   | Descripción                                              |
|-------------|----------------------------------------------------------|
| `x`, `ux`   | Valores de concentración y su incertidumbre              |
| `y`, `uy`   | Valores de respuesta y su incertidumbre                  |
| `y0`, `uy0` | Medida instrumental a invertir y su incertidumbre        |
| `dist_x`    | Distribución para `x` (`"norm"`, `"unif"`, `"triangle"`) |
| `dist_y`    | Distribución para `y`                                    |
| `dist_y0`   | Distribución para `y0`                                   |
| `n_sim`     | Número de simulaciones (por defecto `10000`)             |

------------------------------------------------------------------------

## 📚 Más información

El paquete incluye:

- Diagnóstico de modelo (outliers, heterocedasticidad, R²)
- Selección automática del mejor modelo (lineal simple o cuadrático)
- Exclusión del cúbico por razones metrológicas
- Visualización con `ggplot2` con barras de error en ambos ejes
- Estimación inversa con incertidumbre típica combinada (`uc`)

------------------------------------------------------------------------

> **In God we trust; all others must bring data**. *William Edwards
> Deming*.

> Este paquete fue inspirado en una metrología química accesible  
> para todos, buscando un mundo con medidas más confiables,  
> … *y en una Brujita* 🧚🪄
