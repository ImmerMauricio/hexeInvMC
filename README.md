
# hexeInvMC 🧪

**Regresión inversa con incertidumbre en ambos ejes, basada en
simulación Monte Carlo extendida.**

Este paquete está diseñado para laboratorios que necesitan estimar
concentraciones a partir de señales (`y`), considerando incertidumbres
tanto en `x` como en `y`, con modelos de grado 1 o grado 2, y
permitiendo configuraciones flexibles para las distribuciones de
incertidumbre.

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

# Llamado a la función principal (muestra el gráfico principal)
res<-hexeInvMC(
  x, ux, y, uy,
  y0 = y0, uy0 = uy0,
  dist_x = "norm",   # Puede ser "norm", "unif" o "triangle"
  dist_y = "norm",
  dist_y0 = "norm",
  n_sim = 10000      # Número de simulaciones Monte Carlo
)
```

------------------------------------------------------------------------

## 📈 Gráficos adicionales (opcional)

Además del gráfico principal, la función devuelve objetos `ggplot` que
puedes imprimir cuando los necesites. La forma más simple es guardar el
resultado en `res` y luego imprimir el gráfico deseado.

``` r
print(res$plots$x0)
print(res$plots$resid_fitted)
print(res$plots$qq_resid)
print(res$plots$influence)
```

Gráficos disponibles:

- `res$plots$x0`: distribución Monte Carlo de `x0`
- `res$plots$resid_fitted`: residuos vs valores ajustados
- `res$plots$qq_resid`: gráfico Q–Q de residuos
- `res$plots$influence`: influencia (leverage vs residuo estandarizado;
  referencia Cook)

------------------------------------------------------------------------

## ⚙️ Parámetros disponibles

| Argumento   | Descripción                                              |
|-------------|----------------------------------------------------------|
| `x`, `ux`   | Valores de concentración y su incertidumbre típica       |
| `y`, `uy`   | Valores de respuesta y su incertidumbre típica           |
| `y0`, `uy0` | Medida instrumental a invertir y su incertidumbre típica |
| `dist_x`    | Distribución para `x` (`"norm"`, `"unif"`, `"triangle"`) |
| `dist_y`    | Distribución para `y`                                    |
| `dist_y0`   | Distribución para `y0`                                   |
| `n_sim`     | Número de simulaciones (por defecto `10000`)             |

------------------------------------------------------------------------

## 📚 Más información

El paquete incluye:

- Selección automática de modelo (grado 1 o grado 2), con chequeo de
  curvatura tipo Mandel (ISO 8466-1)
- Chequeo de unicidad del modelo cuadrático en el rango (ISO 8466-2)
- Chequeo informativo de consistencia de la dispersión con `u(y)`
  (ISO/TS 28037)
- Diagnóstico: R², outliers (residuo estandarizado), influencia (Cook),
  heterocedasticidad (ncvTest), normalidad (Shapiro–Wilk)
- Visualización con barras de error en ambos ejes y punto invertido
  `x0 ± u(x0)` en rojo
- Estimación inversa `x0` con incertidumbre típica `u(x0)` e intervalo
  central del 95 %

------------------------------------------------------------------------

> **In God we trust; all others must bring data**.  
> *William Edwards Deming*.

> Este paquete fue inspirado en una metrología química accesible  
> para todos, buscando un mundo con medidas más confiables,  
> … *y en una Brujita* 🧚🪄
