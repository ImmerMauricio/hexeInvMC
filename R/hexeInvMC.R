#' hexeInvMC
#'
#' Regresión inversa con incertidumbre en ambos ejes mediante simulación Monte Carlo extendida.
#' Selección automática de modelo polinómico de grado 1 o 2 (ambos lineales en los parámetros),
#' diagnóstico y reporte metrológico con un gráfico principal.
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
#' @return (Invisiblemente) un objeto tipo lista con clase "hexeInvMC_result".
#' @export
hexeInvMC <- function(x, ux, y, uy, y0, uy0,
                      dist_x = "norm", dist_y = "norm", dist_y0 = "norm",
                      n_sim = 10000) {

  if (!is.numeric(x) || !is.numeric(ux) || !is.numeric(y) || !is.numeric(uy)) {
    stop("x, ux, y, uy deben ser numéricos.")
  }
  if (length(x) != length(ux) || length(y) != length(uy) || length(x) != length(y)) {
    stop("x, ux, y, uy deben tener la misma longitud.")
  }
  if (length(x) < 3) stop("Se requieren al menos 3 puntos de calibración.")
  if (any(!is.finite(x)) || any(!is.finite(ux)) || any(!is.finite(y)) || any(!is.finite(uy))) {
    stop("x, ux, y, uy no deben contener NA/Inf.")
  }
  if (!is.finite(y0) || !is.finite(uy0)) stop("y0 y uy0 deben ser valores finitos.")
  if (any(ux <= 0) || any(uy <= 0) || uy0 <= 0) stop("ux, uy y uy0 deben ser positivos (> 0).")

  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 1000) stop("n_sim debe ser un entero >= 1000.")
  n_sim <- as.integer(n_sim)

  req_pkgs <- c("ggplot2", "car", "EnvStats")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Faltan paquetes requeridos: ", paste(missing_pkgs, collapse = ", "),
         ". Instálalos y vuelve a intentar.")
  }

  gg_ver <- utils::packageVersion("ggplot2")
  use_orientation <- gg_ver >= "4.0.0"
  geom_xerr <- function(mapping = NULL, data = NULL, ...) {
    if (use_orientation) {
      ggplot2::geom_errorbar(mapping = mapping, data = data, orientation = "y", ...)
    } else {
      ggplot2::geom_errorbarh(mapping = mapping, data = data, ...)
    }
  }

  get_sim <- function(dist) {
    switch(dist,
           norm = function(n, mu, sigma) stats::rnorm(n, mu, sigma),
           unif = function(n, mu, sigma) stats::runif(n, mu - sigma, mu + sigma),
           triangle = function(n, mu, sigma) {
             a <- mu - sigma; b <- mu + sigma; c <- mu
             EnvStats::rtri(n, min = a, max = b, mode = c)
           },
           stop("Distribución no soportada: ", dist, ". Usa: 'norm', 'unif' o 'triangle'.")
    )
  }

  sim_x  <- get_sim(dist_x)
  sim_y  <- get_sim(dist_y)
  sim_y0 <- get_sim(dist_y0)

  n <- length(x)

  modelo_g1 <- stats::lm(y ~ x)
  modelo_g2 <- NULL
  if (n >= 4) {
    modelo_g2 <- stats::lm(y ~ x + I(x^2))
  }

  mandel_Fcalc <- NA_real_
  mandel_Fcrit <- NA_real_
  mandel_p <- NA_real_
  mandel_decision <- "NO_EVALUADO"

  if (!is.null(modelo_g2)) {
    sse1 <- sum(stats::residuals(modelo_g1)^2)
    sse2 <- sum(stats::residuals(modelo_g2)^2)
    ds2 <- max(0, sse1 - sse2)
    sy2_sq <- sse2 / (n - 3)

    if (is.finite(sy2_sq) && sy2_sq > 0) {
      mandel_Fcalc <- ds2 / sy2_sq
      mandel_Fcrit <- stats::qf(0.99, df1 = 1, df2 = n - 3)
      mandel_p <- 1 - stats::pf(mandel_Fcalc, df1 = 1, df2 = n - 3)
      mandel_decision <- if (mandel_Fcalc < mandel_Fcrit) "LINEAL_ADECUADO" else "REQUIERE_CUADRATICO"
    }
  }

  modelo_sel <- if (mandel_decision == "REQUIERE_CUADRATICO") "grado_2" else "grado_1"
  modelo_usar <- if (modelo_sel == "grado_2" && !is.null(modelo_g2)) modelo_g2 else modelo_g1

  x_star <- NA_real_
  no_univoco <- FALSE
  x_min <- min(x); x_max <- max(x)

  if (modelo_sel == "grado_2") {
    co <- stats::coef(modelo_usar)
    b <- as.numeric(co[2])
    c <- as.numeric(co[3])

    if (is.finite(c) && c != 0) {
      x_star <- -b / (2 * c)
      if (is.finite(x_star) && x_star > x_min && x_star < x_max) {
        no_univoco <- TRUE
      }
    }
  }

  chi2_obs <- NA_real_
  chi2_df <- NA_real_
  chi2_crit <- NA_real_
  chi2_p <- NA_real_
  chi2_reject <- NA

  if (modelo_sel == "grado_1") {
    a <- stats::coef(modelo_usar)[1]
    b <- stats::coef(modelo_usar)[2]
    r_w <- (y - (a + b * x)) / uy
    chi2_obs <- sum(r_w^2)
    chi2_df <- n - 2

    if (is.finite(chi2_obs) && chi2_df > 0) {
      chi2_crit <- stats::qchisq(0.95, df = chi2_df)
      chi2_p <- 1 - stats::pchisq(chi2_obs, df = chi2_df)
      chi2_reject <- (chi2_obs > chi2_crit)
    }
  }

  resid <- stats::residuals(modelo_usar)
  fitted <- stats::fitted(modelo_usar)

  r2 <- try(summary(modelo_usar)$r.squared, silent = TRUE)
  if (inherits(r2, "try-error")) r2 <- NA_real_

  coef_tab <- summary(modelo_usar)$coefficients

  rs <- try(stats::rstandard(modelo_usar), silent = TRUE)
  if (!inherits(rs, "try-error")) {
    rstd <- as.numeric(rs)
  } else {
    sdr <- stats::sd(resid)
    rstd <- if (is.finite(sdr) && sdr > 0) as.numeric(resid / sdr) else rep(NA_real_, n)
  }

  cd <- try(stats::cooks.distance(modelo_usar), silent = TRUE)
  cooks <- if (!inherits(cd, "try-error")) as.numeric(cd) else rep(NA_real_, n)

  h <- try(stats::hatvalues(modelo_usar), silent = TRUE)
  hatv <- if (!inherits(h, "try-error")) as.numeric(h) else rep(NA_real_, n)

  out_idx <- which(is.finite(rstd) & abs(rstd) > 2)
  cook_thr <- 4 / n
  cook_idx <- which(is.finite(cooks) & cooks > cook_thr)

  ncv_p <- NA_real_
  ncv_stat <- NA_real_
  ncv_df <- NA_real_
  ncv_obj <- try(car::ncvTest(modelo_usar), silent = TRUE)
  if (!inherits(ncv_obj, "try-error")) {
    if (!is.null(ncv_obj$p)) ncv_p <- as.numeric(ncv_obj$p)
    if (!is.null(ncv_obj$ChiSquare)) ncv_stat <- as.numeric(ncv_obj$ChiSquare)
    if (!is.null(ncv_obj$Df)) ncv_df <- as.numeric(ncv_obj$Df)
  }

  sh_p <- NA_real_
  sh_stat <- NA_real_
  if (length(resid) >= 3 && length(resid) <= 5000) {
    sh_obj <- try(stats::shapiro.test(resid), silent = TRUE)
    if (!inherits(sh_obj, "try-error")) {
      sh_p <- as.numeric(sh_obj$p.value)
      sh_stat <- as.numeric(sh_obj$statistic)
    }
  }

  x0_sim <- rep(NA_real_, n_sim)

  for (i in seq_len(n_sim)) {
    x_sim  <- sim_x(n, x, ux)
    y_sim  <- sim_y(n, y, uy)
    y0_sim <- sim_y0(1, y0, uy0)

    if (modelo_sel == "grado_1") {
      mod <- stats::lm(y_sim ~ x_sim)

      b0 <- stats::coef(mod)[1]
      b1 <- stats::coef(mod)[2]

      se_b0 <- summary(mod)$coefficients[1, 2]
      se_b1 <- summary(mod)$coefficients[2, 2]

      b0_sim <- stats::rnorm(1, b0, se_b0)
      b1_sim <- stats::rnorm(1, b1, se_b1)

      x0_sim[i] <- (y0_sim - b0_sim) / b1_sim

    } else {
      mod <- stats::lm(y_sim ~ x_sim + I(x_sim^2))

      b0 <- stats::coef(mod)[1]
      b1 <- stats::coef(mod)[2]
      b2 <- stats::coef(mod)[3]

      se_b0 <- summary(mod)$coefficients[1, 2]
      se_b1 <- summary(mod)$coefficients[2, 2]
      se_b2 <- summary(mod)$coefficients[3, 2]

      b0_sim <- stats::rnorm(1, b0, se_b0)
      b1_sim <- stats::rnorm(1, b1, se_b1)
      b2_sim <- stats::rnorm(1, b2, se_b2)

      a_q <- b2_sim
      b_q <- b1_sim
      c_q <- b0_sim - y0_sim

      disc <- b_q^2 - 4 * a_q * c_q
      if (is.finite(disc) && disc >= 0 && is.finite(a_q) && a_q != 0) {
        r1 <- (-b_q + sqrt(disc)) / (2 * a_q)
        r2q <- (-b_q - sqrt(disc)) / (2 * a_q)
        x0_sim[i] <- ifelse(min(r1, r2q) > 0, min(r1, r2q), max(r1, r2q))
      } else {
        x0_sim[i] <- NA_real_
      }
    }
  }

  validas <- is.finite(x0_sim)
  x_valid <- x0_sim[validas]

  if (length(x_valid) == 0) {
    x0_hat <- NA_real_
    u_x0 <- NA_real_
    ic95 <- c(NA_real_, NA_real_)
  } else {
    x0_hat <- mean(x_valid)
    u_x0 <- stats::sd(x_valid)
    ic95 <- stats::quantile(x_valid, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  }

  frac_valid <- mean(validas)

  metro <- function(alerta, msg) {
    if (isTRUE(alerta)) paste0("Ojo: ", msg) else msg
  }
  interpret_p <- function(p) {
    if (!is.finite(p)) return("No evaluado / no aplica.")
    if (p < 0.05) return("Evidencia de problema (p < 0,05).")
    "No hay evidencia suficiente de problema (p ≥ 0,05)."
  }

  mandel_txt_stat <- if (!is.finite(mandel_Fcalc) || !is.finite(mandel_Fcrit)) {
    "No evaluado (se requieren ≥ 4 puntos y cálculo estable)."
  } else if (mandel_Fcalc < mandel_Fcrit) {
    "Recta adecuada (Fcalc < Fcrit al 99%)."
  } else {
    "Se justifica grado 2 (Fcalc ≥ Fcrit al 99%)."
  }

  mandel_txt_metro <- if (!is.null(modelo_g2) && is.finite(mandel_Fcalc) && is.finite(mandel_Fcrit) && mandel_Fcalc >= mandel_Fcrit) {
    "Se detecta curvatura; se ajustó un modelo con x^2."
  } else if (!is.null(modelo_g2) && is.finite(mandel_Fcalc) && is.finite(mandel_Fcrit) && mandel_Fcalc < mandel_Fcrit) {
    "No se observa curvatura sistemática en el rango; una recta describe bien la forma."
  } else if (is.null(modelo_g2)) {
    metro(TRUE, "con menos de 4 puntos no se puede evaluar curvatura con Mandel.")
  } else {
    metro(TRUE, "no se pudo evaluar Mandel con confianza.")
  }

  quad_txt_stat <- if (modelo_sel != "grado_2") {
    "No aplica (modelo seleccionado no es cuadrático)."
  } else if (!is.finite(x_star)) {
    "No evaluado (c ≈ 0 o x* no finito)."
  } else if (no_univoco) {
    "No unívoco: x* dentro del rango."
  } else {
    "Unívoco: x* fuera del rango."
  }

  quad_txt_metro <- if (modelo_sel != "grado_2") {
    "Este chequeo solo importa si se usa modelo con x^2."
  } else if (!is.finite(x_star)) {
    metro(TRUE, "la curvatura es muy pequeña/inestable; x* no es interpretable.")
  } else if (no_univoco) {
    metro(TRUE, paste0("la curva tiene un máximo/mínimo dentro del rango (x* = ", round(x_star, 6),
                       "); la inversión y→x deja de ser única."))
  } else {
    "En el rango, la curva se comporta de forma monótona; la inversión y→x es única."
  }

  chi2_txt_stat <- if (modelo_sel != "grado_1") {
    "No aplica (prueba informativa definida para recta)."
  } else if (!is.finite(chi2_obs) || !is.finite(chi2_crit)) {
    "No evaluado."
  } else if (isTRUE(chi2_reject)) {
    "Alerta: chi2_obs > chi2_crit al 95%."
  } else {
    "OK: chi2_obs ≤ chi2_crit al 95%."
  }

  k_infl <- NA_real_
  if (is.finite(chi2_obs) && is.finite(chi2_df) && chi2_df > 0) {
    k_infl <- sqrt(chi2_obs / chi2_df)
  }

  chi2_txt_metro <- if (modelo_sel != "grado_1") {
    "Chequeo informativo para recta."
  } else if (!is.finite(chi2_obs) || !is.finite(chi2_crit)) {
    metro(TRUE, "no se pudo evaluar consistencia con u(y).")
  } else if (isTRUE(chi2_reject)) {
    if (is.finite(k_infl)) {
      metro(TRUE, paste0("la nube está un poco más dispersa de lo que predicen tus u(y). ",
                         "Como guía rápida, k ≈ ", round(k_infl, 3),
                         " haría la dispersión más coherente con u(y) (informativo)."))
    } else {
      metro(TRUE, "la nube está un poco más dispersa de lo que predicen tus u(y).")
    }
  } else {
    "La dispersión observada es coherente con tus u(y)."
  }

  r2_alert <- is.finite(r2) && r2 < 0.98
  r2_txt_metro <- if (!is.finite(r2)) {
    metro(TRUE, "no se pudo calcular R².")
  } else if (r2_alert) {
    metro(TRUE, "R² es bajo según el umbral actual; revisa rango, datos o forma funcional.")
  } else {
    "R² es alto; es buena señal de tendencia global, pero no valida por sí solo u(y)."
  }

  out_alert <- length(out_idx) > 0
  out_txt_metro <- if (out_alert) {
    metro(TRUE, paste0("hay puntos que se salen de la tendencia. Índices {", paste(out_idx, collapse = ", "), "}."))
  } else {
    "No se observan outliers moderados con este criterio."
  }

  cook_alert <- length(cook_idx) > 0
  cook_txt_metro <- if (cook_alert) {
    metro(TRUE, paste0("uno o más puntos dominan el ajuste. Índices {", paste(cook_idx, collapse = ", "), "}."))
  } else {
    "No se observan puntos dominantes según Cook > 4/n."
  }

  ncv_alert <- is.finite(ncv_p) && ncv_p < 0.05
  ncv_txt_metro <- if (!is.finite(ncv_p)) {
    metro(TRUE, "no se pudo evaluar heterocedasticidad.")
  } else if (ncv_alert) {
    metro(TRUE, "la variabilidad parece cambiar con el nivel; podría requerir ponderación.")
  } else {
    "No hay evidencia clara de cambio de variabilidad con el nivel (según ncvTest)."
  }

  sh_alert <- is.finite(sh_p) && sh_p < 0.05
  sh_txt_metro <- if (!is.finite(sh_p)) {
    metro(TRUE, "no se pudo evaluar la forma de los residuos.")
  } else if (sh_alert) {
    metro(TRUE, "los residuos no parecen gaussianos; puede haber colas pesadas o asimetría.")
  } else {
    "Los residuos se ven razonablemente gaussianos (según Shapiro)."
  }

  tests <- data.frame(
    prueba = c(
      "Mandel (ISO 8466-1) Fcalc (99%)",
      "Unicidad cuadrático (ISO 8466-2) x* = -b/(2c)",
      "Chi-cuadrado (ISO/TS 28037; informativo; usa uy)",
      "R² del ajuste",
      "Outliers (|rstandard| > 2)",
      "Influencia (Cook > 4/n)",
      "Heterocedasticidad (car::ncvTest)",
      "Normalidad de residuos (Shapiro–Wilk)"
    ),
    estadistico = c(mandel_Fcalc, x_star, chi2_obs, as.numeric(r2),
                    length(out_idx), length(cook_idx), as.numeric(ncv_stat), as.numeric(sh_stat)),
    gl1 = c(1, NA_real_, chi2_df, NA_real_, NA_real_, NA_real_, as.numeric(ncv_df), NA_real_),
    gl2 = c(ifelse(is.finite(n), n - 3, NA_real_), NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    p_value = c(mandel_p, NA_real_, chi2_p, NA_real_, NA_real_, NA_real_, ncv_p, sh_p),
    umbral = c(mandel_Fcrit, NA_real_, chi2_crit, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )

  tests$interpretacion_estadistica <- c(
    mandel_txt_stat,
    quad_txt_stat,
    chi2_txt_stat,
    if (is.finite(r2)) { if (r2 < 0.98) "Alerta: R² < 0,98 (umbral actual)." else "R² sin alerta (umbral actual)." } else "No evaluado.",
    if (out_alert) "Hay outliers moderados." else "Sin outliers moderados.",
    if (cook_alert) "Hay puntos influyentes." else "Sin puntos influyentes.",
    interpret_p(ncv_p),
    interpret_p(sh_p)
  )

  tests$interpretacion_metrologica <- c(
    mandel_txt_metro,
    quad_txt_metro,
    chi2_txt_metro,
    r2_txt_metro,
    out_txt_metro,
    cook_txt_metro,
    ncv_txt_metro,
    sh_txt_metro
  )

  alerts <- character(0)
  notes <- character(0)

  if (modelo_sel == "grado_2" && no_univoco) alerts <- c(alerts, "ISO 8466-2: inversión no unívoca (x* dentro del rango).")
  if (is.finite(frac_valid) && frac_valid < 0.95) alerts <- c(alerts, "Monte Carlo: fracción válida < 0,95 (inversión inestable).")
  if (modelo_sel == "grado_1" && isTRUE(chi2_reject)) alerts <- c(alerts, "ISO/TS 28037 (informativo): dispersión mayor a la esperada según u(y).")
  if (r2_alert) alerts <- c(alerts, "R² bajo según umbral actual.")
  if (out_alert) alerts <- c(alerts, paste0("Outliers moderados: {", paste(out_idx, collapse = ", "), "}."))
  if (cook_alert) alerts <- c(alerts, paste0("Puntos influyentes: {", paste(cook_idx, collapse = ", "), "}."))
  if (ncv_alert) alerts <- c(alerts, "Heterocedasticidad (ncvTest p < 0,05).")
  if (sh_alert) alerts <- c(alerts, "No-normalidad de residuos (Shapiro p < 0,05).")

  if (!is.null(modelo_g2) && mandel_decision == "LINEAL_ADECUADO") notes <- c(notes, "Mandel: recta adecuada al 99%.")
  if (!is.null(modelo_g2) && mandel_decision == "REQUIERE_CUADRATICO") notes <- c(notes, "Mandel: se aplicó grado 2 por curvatura al 99%.")
  if (modelo_sel == "grado_1" && is.finite(chi2_obs) && is.finite(chi2_crit) && !isTRUE(chi2_reject)) notes <- c(notes, "Chi-cuadrado (informativo): coherente con u(y).")

  status <- "APTO"
  if ((modelo_sel == "grado_2" && no_univoco) || (is.finite(frac_valid) && frac_valid < 0.95)) {
    status <- "NO_APTO"
  } else if (length(alerts) > 0) {
    status <- "APTO_CON_ALERTAS"
  }

  df_cal <- data.frame(x = x, y = y, ux = ux, uy = uy)
  punto <- data.frame(x = x0_hat, y = y0, ux = u_x0, uy = uy0)

  g_main <- ggplot2::ggplot(df_cal, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = 2, color = "black", shape = 16) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = y - uy, ymax = y + uy), width = 0.004, color = "gray60") +
    geom_xerr(ggplot2::aes(xmin = x - ux, xmax = x + ux), height = 0.004, color = "gray60") +
    ggplot2::geom_point(data = punto, ggplot2::aes(x = x, y = y), shape = 21, fill = "red", color = "black", size = 3) +
    ggplot2::geom_errorbar(data = punto, ggplot2::aes(ymin = y - uy, ymax = y + uy), width = 0.1, color = "red", linewidth = 1) +
    geom_xerr(mapping = ggplot2::aes(xmin = x - ux, xmax = x + ux), data = punto, height = 0.1, color = "red", linewidth = 1) +
    ggplot2::labs(
      title = paste0("hexeInvMC — modelo: ", modelo_sel, " | decisión: ", status),
      subtitle = paste0(
        "x0 = ", round(x0_hat, 6),
        " | u(x0) = ", round(u_x0, 6),
        " | IC95% = [", round(ic95[1], 6), ", ", round(ic95[2], 6), "]",
        " | R² = ", ifelse(is.finite(r2), round(r2, 6), NA)
      ),
      x = "Concentración (x)",
      y = "Respuesta (y)"
    )

  coefs <- stats::coef(modelo_usar)
  if (modelo_sel == "grado_1") {
    g_main <- g_main + ggplot2::stat_function(fun = function(x) coefs[1] + coefs[2] * x, color = "blue", linewidth = 1.2)
  } else {
    g_main <- g_main + ggplot2::stat_function(fun = function(x) coefs[1] + coefs[2] * x + coefs[3] * x^2, color = "blue", linewidth = 1.2)
  }

  if (length(x_valid) > 0) {
    g_x0 <- ggplot2::ggplot(data.frame(x0 = x_valid), ggplot2::aes(x = x0)) +
      ggplot2::geom_histogram(bins = 40, color = "black", fill = "gray80") +
      ggplot2::geom_vline(xintercept = x0_hat, linewidth = 1.1) +
      ggplot2::geom_vline(xintercept = ic95, linetype = 2, linewidth = 0.9) +
      ggplot2::labs(
        title = "Distribución Monte Carlo de x0",
        subtitle = paste0("media = ", round(x0_hat, 6), " | IC95% = [", round(ic95[1], 6), ", ", round(ic95[2], 6), "]"),
        x = "x0 simulado",
        y = "Frecuencia"
      )
  } else {
    g_x0 <- ggplot2::ggplot() + ggplot2::labs(
      title = "Distribución Monte Carlo de x0",
      subtitle = "No hubo simulaciones válidas de x0.",
      x = "x0 simulado",
      y = "Frecuencia"
    )
  }

  df_rf <- data.frame(fitted = fitted, resid = resid)
  g_rf <- ggplot2::ggplot(df_rf, ggplot2::aes(x = fitted, y = resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_smooth(se = FALSE, linewidth = 0.9) +
    ggplot2::labs(
      title = "Diagnóstico — residuos vs valores ajustados",
      subtitle = "Se espera dispersión sin patrón claro alrededor de cero.",
      x = "Valores ajustados",
      y = "Residuos"
    )

  df_q <- data.frame(resid = resid)
  g_qq <- ggplot2::ggplot(df_q, ggplot2::aes(sample = resid)) +
    ggplot2::stat_qq(size = 2, color = "black") +
    ggplot2::stat_qq_line(linewidth = 0.9) +
    ggplot2::labs(
      title = "Diagnóstico — QQ-plot de residuos",
      subtitle = "Desviaciones sistemáticas de la línea sugieren colas pesadas o asimetría.",
      x = "Cuantiles teóricos",
      y = "Cuantiles de residuos"
    )

  df_inf <- data.frame(idx = seq_len(n), leverage = hatv, rstandard = rstd, cooks = cooks)
  lev_thr <- (2 * ncol(stats::model.matrix(modelo_usar))) / n
  g_inf <- ggplot2::ggplot(df_inf, ggplot2::aes(x = leverage, y = rstandard)) +
    ggplot2::geom_hline(yintercept = c(-2, 0, 2), linetype = c(3, 2, 3)) +
    ggplot2::geom_vline(xintercept = lev_thr, linetype = 2) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::labs(
      title = "Diagnóstico — influencia (leverage vs residuo estandarizado)",
      subtitle = paste0("Referencia: leverage ~ 2p/n = ", round(lev_thr, 3), " | Cook > 4/n = ", round(cook_thr, 3)),
      x = "Leverage",
      y = "Residuo estandarizado"
    )

  plots <- list(
    main = g_main,
    x0 = g_x0,
    resid_fitted = g_rf,
    qq_resid = g_qq,
    influence = g_inf
  )

  res <- list(
    estimate = list(x0_hat = x0_hat, u_x0 = u_x0, ic95_x0 = ic95, y0 = y0, u_y0 = uy0),
    mc = list(n_sim = n_sim, n_validas = sum(validas), frac_validas = frac_valid, x0_sim = x0_sim),
    model = list(
      seleccionado = modelo_sel,
      formula = stats::formula(modelo_usar),
      coeficientes = stats::coef(modelo_usar),
      tabla_coeficientes = coef_tab,
      r2 = as.numeric(r2),
      mandel = list(Fcalc = mandel_Fcalc, Fcrit_99 = mandel_Fcrit, p_value = mandel_p, decision = mandel_decision),
      unicidad_cuadratico = list(x_star = x_star, no_univoco = no_univoco, rango = c(x_min, x_max)),
      chi2_informativo = list(chi2_obs = chi2_obs, df = chi2_df, chi2_crit_95 = chi2_crit, p_value = chi2_p, reject = chi2_reject)
    ),
    tests = tests,
    status = status,
    alerts = alerts,
    notes = notes,
    plot = g_main,
    plots = plots,
    call = match.call()
  )
  class(res) <- "hexeInvMC_result"

  print(g_main)

  cat("===========================================\n")
  cat("        [RESUMEN FINAL hexeInvMC]\n")
  cat("===========================================\n")
  cat("Modelo seleccionado     :", modelo_sel, "\n\n")
  cat("Coeficientes del modelo (estimación, incertidumbre típica, valor t, valor p):\n")
  print(coef_tab)

  cat("\nR² del ajuste            :", ifelse(is.finite(r2), round(r2, 6), NA), "\n")
  cat("x0 (media MC)            :", round(x0_hat, 6), "\n")
  cat("u(x0) (sd MC)            :", round(u_x0, 6), "\n")
  cat("IC95% para x0            : [", round(ic95[1], 6), ", ", round(ic95[2], 6), "]\n", sep = "")
  cat("y0 observado             :", y0, "\n")
  cat("u(y0)                    :", uy0, "\n")
  cat("Simulaciones válidas     :", sum(validas), "/", n_sim, " (", round(100 * frac_valid, 2), " %)\n", sep = "")

  cat("\nPruebas / chequeos:\n")
  print(tests[, c("prueba", "estadistico", "gl1", "gl2", "p_value", "umbral",
                  "interpretacion_estadistica", "interpretacion_metrologica")],
        row.names = FALSE)

  cat("\nDecisión                 :", status, "\n")
  if (length(alerts) == 0) {
    cat("Alertas                  : Ninguna.\n")
  } else {
    cat("Alertas:\n")
    for (r in alerts) cat(" - ", r, "\n", sep = "")
  }
  if (length(notes) > 0) {
    cat("\nNotas:\n")
    for (r in notes) cat(" - ", r, "\n", sep = "")
  }
  cat("===========================================\n")

  invisible(res)
}
