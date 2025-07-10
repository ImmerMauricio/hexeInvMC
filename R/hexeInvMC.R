#' hexeInvMC: Regresión inversa con incertidumbre
#'
#' Realiza estimación inversa considerando incertidumbre en ambos ejes mediante simulación Monte Carlo extendida.
#'
#' @param x Vector numérico con valores de x (por ejemplo, concentración)
#' @param ux Vector numérico con incertidumbre estándar asociada a x
#' @param y Vector numérico con valores de y (por ejemplo, absorbancia)
#' @param uy Vector numérico con incertidumbre estándar asociada a y
#' @param y0 Valor de respuesta instrumental observado para estimar su correspondiente x
#' @param uy0 Incertidumbre estándar asociada a y0
#' @param n_sim Número de simulaciones (por defecto 10000)
#' @param dist_x Distribución para simular x: "norm", "unif", "triangle"
#' @param dist_y Distribución para simular y: "norm", "unif", "triangle"
#' @param dist_y0 Distribución para simular y0: "norm", "unif", "triangle"
#'
#' @return Lista con el modelo seleccionado, el valor estimado de x y su incertidumbre combinada uc(x).
#' También imprime diagnóstico y genera una gráfica con ggplot2.
#' @export
#'
#' @examples
#' x <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
#' ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
#' y <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
#' uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)
#' hexeInvMC(x, ux, y, uy, y0 = 10.5, uy0 = 0.25)

hexeInvMC <- function(x, ux, y, uy, y0, uy0,
                      n_sim = 10000,
                      dist_x = "norm",
                      dist_y = "norm",
                      dist_y0 = "norm") {

  # --- Simuladores ---
  get_simulator <- function(distribution) {
    switch(distribution,
           norm     = function(n, mu, sigma) rnorm(n, mu, sigma),
           unif     = function(n, mu, sigma) runif(n, mu - sigma, mu + sigma),
           triangle = function(n, mu, sigma) {
             a <- mu - sigma; b <- mu + sigma; c <- mu
             EnvStats::rtri(n, min = a, max = b, mode = c)
           },
           stop("Distribución no soportada: ", distribution)
    )
  }

  sim_x  <- get_simulator(dist_x)
  sim_y  <- get_simulator(dist_y)
  sim_y0 <- get_simulator(dist_y0)

  # --- Diagnóstico estadístico ---
  check_model_aptitude <- function(modelo, x, y) {
    razones <- c()

    std_res <- rstandard(modelo)
    if (any(abs(std_res) > 2)) {
      razones <- c(razones, "Presencia de outliers moderados.")
    }

    prueba_bp <- try(bptest(modelo), silent = TRUE)
    if (!inherits(prueba_bp, "try-error") && prueba_bp$p.value < 0.05) {
      razones <- c(razones, "Varianza no homogénea detectada.")
    }

    r2 <- summary(modelo)$r.squared
    if (r2 < 0.98) {
      razones <- c(razones, "R² bajo.")
    }

    return(razones)
  }

  # --- Ajuste y selección de modelo ---
  modelo_lineal_simple  <- lm(y ~ x)
  modelo_pol_cuadratico <- lm(y ~ x + I(x^2))
  modelo_pol_cubico     <- lm(y ~ x + I(x^2) + I(x^3))

  f_values <- c(
    lineal_simple    = anova(modelo_lineal_simple)$"F value"[1],
    pol_cuadratico   = anova(modelo_pol_cuadratico)$"F value"[1],
    pol_cubico       = anova(modelo_pol_cubico)$"F value"[1]
  )

  modelo_ajuste_mas_alto <- names(which.max(f_values))
  if (modelo_ajuste_mas_alto == "pol_cubico") {
    f_values <- f_values[names(f_values) != "pol_cubico"]
  }
  modelo_seleccionado <- names(which.max(f_values))

  modelo_usar <- switch(modelo_seleccionado,
                        lineal_simple  = modelo_lineal_simple,
                        pol_cuadratico = modelo_pol_cuadratico)

  diagnostico <- check_model_aptitude(modelo_usar, x, y)

  # --- Monte Carlo extendido ---
  x_estimados <- numeric(n_sim)

  for (i in seq_len(n_sim)) {
    x_sim <- sim_x(length(x), x, ux)
    y_sim <- sim_y(length(y), y, uy)
    y_simulado <- sim_y0(1, y0, uy0)

    if (modelo_seleccionado == "lineal_simple") {
      modelo <- lm(y_sim ~ x_sim)
      b0 <- coef(modelo)[1]; b1 <- coef(modelo)[2]
      se_b0 <- summary(modelo)$coefficients[1, 2]
      se_b1 <- summary(modelo)$coefficients[2, 2]
      b0_sim <- rnorm(1, b0, se_b0)
      b1_sim <- rnorm(1, b1, se_b1)
      x_estimados[i] <- (y_simulado - b0_sim) / b1_sim
    } else {
      x2_sim <- x_sim^2
      modelo <- lm(y_sim ~ x_sim + I(x2_sim))
      b0 <- coef(modelo)[1]; b1 <- coef(modelo)[2]; b2 <- coef(modelo)[3]
      se_b0 <- summary(modelo)$coefficients[1, 2]
      se_b1 <- summary(modelo)$coefficients[2, 2]
      se_b2 <- summary(modelo)$coefficients[3, 2]
      b0_sim <- rnorm(1, b0, se_b0)
      b1_sim <- rnorm(1, b1, se_b1)
      b2_sim <- rnorm(1, b2, se_b2)
      a <- b2_sim; b <- b1_sim; c <- b0_sim - y_simulado
      disc <- b^2 - 4*a*c
      if (disc >= 0) {
        x1 <- (-b + sqrt(disc)) / (2*a)
        x2 <- (-b - sqrt(disc)) / (2*a)
        x_estimados[i] <- ifelse(min(x1, x2) > 0, min(x1, x2), max(x1, x2))
      } else {
        x_estimados[i] <- NA
      }
    }
  }

  valor_medido <- mean(x_estimados, na.rm = TRUE)
  uc_x <- sd(x_estimados, na.rm = TRUE)

  # --- Gráfico final ---
  df_sim <- data.frame(x = x_sim, y = y_sim, ux = ux, uy = uy)
  punto_estimado <- data.frame(x = valor_medido, y = y0, ux = uc_x, uy = uy0)

  g <- ggplot(df_sim, aes(x = x, y = y)) +
    geom_point(size = 2, color = "black", shape = 16) +
    geom_errorbar(aes(ymin = y - uy, ymax = y + uy), width = 0.004, color = "gray60") +
    geom_errorbarh(aes(xmin = x - ux, xmax = x + ux), height = 0.004, color = "gray60") +
    geom_point(data = punto_estimado, aes(x = x, y = y), shape = 21, fill = "red", color = "black", size = 3) +
    geom_errorbar(data = punto_estimado, aes(ymin = y - uy, ymax = y + uy), width = 0.1, color = "red", linewidth = 1) +
    geom_errorbarh(data = punto_estimado, aes(xmin = x - ux, xmax = x + ux), height = 0.1, color = "red", linewidth = 1) +
    labs(
      title = paste0("Modelo de calibración con incertidumbre (", modelo_seleccionado, ")"),
      subtitle = paste0("x estimado = ", round(valor_medido, 5),
                        "   |   uc(x) = ", round(uc_x, 5),
                        "   |   y0 = ", y0,
                        "   |   uc(y0) = ", uy0),
      x = "Concentración simulada (x)",
      y = "Respuesta instrumental simulada (y)"
    )

  coefs <- coef(modelo_usar)
  if (modelo_seleccionado == "lineal_simple") {
    g <- g + stat_function(fun = function(x) coefs[1] + coefs[2]*x,
                           color = "blue", linewidth = 1.2)
  } else {
    g <- g + stat_function(fun = function(x) coefs[1] + coefs[2]*x + coefs[3]*x^2,
                           color = "blue", linewidth = 1.2)
  }

  plot(g)

  # --- Resultado final ---
  cat("\n===========================================\n")
  cat("          [RESUMEN FINAL hexeInvMC]\n")
  cat("===========================================\n")
  cat("Modelo seleccionado     :", modelo_seleccionado, "\n")
  cat("x estimado              :", round(valor_medido, 5), "\n")
  cat("uc(x)                   :", round(uc_x, 5), "\n")
  cat("y0 observado            :", y0, "\n")
  cat("uc(y0)                  :", uy0, "\n")

  if (length(diagnostico) == 0) {
    cat("Diagnóstico del modelo  : APTO para fines metrológicos.\n")
  } else {
    cat("Diagnóstico del modelo  : NO ÓPTIMO.\n")
    cat("Razones:\n")
    for (r in diagnostico) cat(" -", r, "\n")
  }
  cat("===========================================\n")

  invisible(list(
    modelo = modelo_seleccionado,
    x_estimado = valor_medido,
    uc_x = uc_x,
    diagnostico = diagnostico
  ))
}
