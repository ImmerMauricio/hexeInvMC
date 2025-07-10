#' hexeInvMC
#'
#' Regresión inversa con incertidumbre en ambos ejes utilizando simulación Monte Carlo extendida.
#' Selección automática de modelo (lineal simple o cuadrático), diagnóstico y visualización.
#'
#' @param x Vector numérico con los valores de concentración (x).
#' @param ux Vector numérico con las incertidumbres típicas de x.
#' @param y Vector numérico con los valores de respuesta (y).
#' @param uy Vector numérico con las incertidumbres típicas de y.
#' @param y0 Valor observado a invertir.
#' @param uy0 Incertidumbre típica de y0.
#' @param dist_x Distribución para x: "norm", "unif" o "triangle".
#' @param dist_y Distribución para y: "norm", "unif" o "triangle".
#' @param dist_y0 Distribución para y0: "norm", "unif" o "triangle".
#' @param n_sim Número de simulaciones Monte Carlo. Por defecto 10000.
#'
#' @return Un resumen de resultados con x estimado, uc(x), modelo seleccionado y diagnóstico.
#' @examples
#' x <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
#' ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
#' y <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
#' uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)
#' hexeInvMC(x, ux, y, uy, y0 = 10.5, uy0 = 0.25,
#'           dist_x = "norm", dist_y = "norm", dist_y0 = "norm",
#'           n_sim = 10000)
#'
#' @export
hexeInvMC <- function(x, ux, y, uy, y0, uy0,
                      dist_x = "norm", dist_y = "norm", dist_y0 = "norm",
                      n_sim = 10000) {

  # Instala y carga paquetes necesarios
  paquetes <- c("ggplot2", "car", "EnvStats")
  for (pkg in paquetes) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }

  # Función de simulación según distribución
  get_sim <- function(dist) {
    switch(dist,
           norm     = function(n, mu, sigma) rnorm(n, mu, sigma),
           unif     = function(n, mu, sigma) runif(n, mu - sigma, mu + sigma),
           triangle = function(n, mu, sigma) {
             a <- mu - sigma; b <- mu + sigma; c <- mu
             EnvStats::rtri(n, min = a, max = b, mode = c)
           },
           stop("Distribución no soportada: ", dist)
    )
  }

  sim_x <- get_sim(dist_x)
  sim_y <- get_sim(dist_y)
  sim_y0 <- get_sim(dist_y0)

  # Modelos
  modelo_lineal <- lm(y ~ x)
  modelo_cuad <- lm(y ~ x + I(x^2))
  modelo_cub <- lm(y ~ x + I(x^2) + I(x^3))

  f1 <- anova(modelo_lineal)$"F value"[1]
  f2 <- anova(modelo_cuad)$"F value"[1]
  f3 <- anova(modelo_cub)$"F value"[1]

  fval <- c(lineal_simple = f1, pol_cuadratico = f2, pol_cubico = f3)
  best_model <- names(which.max(fval))
  if (best_model == "pol_cubico") fval <- fval[names(fval) != "pol_cubico"]
  modelo_sel <- names(which.max(fval))

  modelo_usar <- switch(modelo_sel,
                        lineal_simple = modelo_lineal,
                        pol_cuadratico = modelo_cuad)

  check_model <- function(modelo, x, y) {
    razones <- c()
    if (any(abs(rstandard(modelo)) > 2)) razones <- c(razones, "Presencia de outliers moderados.")
    prueba_bp <- try(bptest(modelo), silent = TRUE)
    if (!inherits(prueba_bp, "try-error") && prueba_bp$p.value < 0.05)
      razones <- c(razones, "Varianza no homogénea detectada.")
    if (summary(modelo)$r.squared < 0.98) razones <- c(razones, "R² bajo.")
    return(razones)
  }

  diagnostico <- check_model(modelo_usar, x, y)

  # Simulación Monte Carlo extendida
  x_estimados <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    x_sim <- sim_x(length(x), x, ux)
    y_sim <- sim_y(length(y), y, uy)
    y0_sim <- sim_y0(1, y0, uy0)

    if (modelo_sel == "lineal_simple") {
      mod <- lm(y_sim ~ x_sim)
      b0 <- coef(mod)[1]; b1 <- coef(mod)[2]
      se_b0 <- summary(mod)$coefficients[1, 2]
      se_b1 <- summary(mod)$coefficients[2, 2]
      b0_sim <- rnorm(1, b0, se_b0)
      b1_sim <- rnorm(1, b1, se_b1)
      x_estimados[i] <- (y0_sim - b0_sim) / b1_sim
    } else {
      x2 <- x_sim^2
      mod <- lm(y_sim ~ x_sim + I(x2))
      b0 <- coef(mod)[1]; b1 <- coef(mod)[2]; b2 <- coef(mod)[3]
      se_b0 <- summary(mod)$coefficients[1, 2]
      se_b1 <- summary(mod)$coefficients[2, 2]
      se_b2 <- summary(mod)$coefficients[3, 2]
      b0_sim <- rnorm(1, b0, se_b0)
      b1_sim <- rnorm(1, b1, se_b1)
      b2_sim <- rnorm(1, b2, se_b2)
      a <- b2_sim; b <- b1_sim; c <- b0_sim - y0_sim
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

  df_sim <- data.frame(x = x, y = y, ux = ux, uy = uy)
  punto <- data.frame(x = valor_medido, y = y0, ux = uc_x, uy = uy0)

  g <- ggplot(df_sim, aes(x = x, y = y)) +
    geom_point(size = 2, color = "black", shape = 16) +
    geom_errorbar(aes(ymin = y - uy, ymax = y + uy), width = 0.004, color = "gray60") +
    geom_errorbarh(aes(xmin = x - ux, xmax = x + ux), height = 0.004, color = "gray60") +
    geom_point(data = punto, aes(x = x, y = y), shape = 21, fill = "red", color = "black", size = 3) +
    geom_errorbar(data = punto, aes(ymin = y - uy, ymax = y + uy), width = 0.1, color = "red", linewidth = 1) +
    geom_errorbarh(data = punto, aes(xmin = x - ux, xmax = x + ux), height = 0.1, color = "red", linewidth = 1) +
    labs(title = paste0("Modelo de calibración con incertidumbre (", modelo_sel, ")"),
         subtitle = paste0("x estimado = ", round(valor_medido, 5),
                           " | uc(x) = ", round(uc_x, 5),
                           " | y0 = ", y0, " | uc(y0) = ", uy0),
         x = "Concentración (x)", y = "Respuesta (y)")

  coefs <- coef(modelo_usar)
  if (modelo_sel == "lineal_simple") {
    g <- g + stat_function(fun = function(x) coefs[1] + coefs[2]*x, color = "blue", linewidth = 1.2)
  } else {
    g <- g + stat_function(fun = function(x) coefs[1] + coefs[2]*x + coefs[3]*x^2, color = "blue", linewidth = 1.2)
  }

  print(g)

  cat("===========================================\n")
  cat("        [RESUMEN FINAL hexeInvMC]\n")
  cat("===========================================\n")
  cat("Modelo seleccionado     :", modelo_sel, "\n")
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
}
