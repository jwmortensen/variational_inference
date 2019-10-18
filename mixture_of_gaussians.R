library(shiny)
library(tidyverse)

# This is an implementation of the CAVI algorithm for the
# univariate mixture of gaussians outlined in
# "Variational Inference: A Review for Statisticians" by Blei et. al. 2016

simulate_data = function(K = 3, sigma2 = 5, n = 100) {
  mu = rnorm(K, sd = sqrt(sigma2))
  c_i = factor(sample.int(n = K, size = n, replace = TRUE))
  C = Matrix::sparse.model.matrix(~ -1 + c_i)
  x = rnorm(n, mean = as.vector(C %*% mu), sd = 1)

  res = list(
    K = K,
    sigma2 = sigma2,
    n = n,
    mu = mu,
    c_i = c_i,
    C = C,
    x = x
  )

  return(res)
}

make_elbo = function(x, sigma2, K) {

  n = length(x)
  X = t(x %*% t(rep(1, K)))

  elbo = function(m, s2, psi) {

    testthat::expect_equal(length(m), length(s2))
    testthat::expect_equal(dim(X), dim(psi))
    testthat::expect_equal(dim(X)[2], length(x))
    testthat::expect_equal(dim(X)[1], length(m))

    aa = sum(-(s2 + m^2) / (2 * sigma2))
    # bb = -n * log(K); this is a constant, so we don't need to include it in the sum
    cc = sum(psi * (m * X - 0.5 * (s2 + m^2)))
    dd = sum(psi * log(psi))
    ee = sum(-0.5 * (1 + log(s2)))

    res = aa + cc - dd - ee
    return(res)
  }

  return(elbo)
}

make_update_params = function(x, sigma2, K) {

  X = t(x %*% t(rep(1, K)))

  update_params = function(m, s2, psi) {

    testthat::expect_equal(length(m), length(s2))
    testthat::expect_equal(dim(X), dim(psi))
    testthat::expect_equal(dim(X)[2], length(x))
    testthat::expect_equal(dim(X)[1], length(m))

    # update psi
    psi_prop = exp(m * X - 0.5 * (m^2 + s2))
    psi = t(t(psi_prop) / colSums(psi_prop))

    # update m
    m = as.vector(rowSums(psi * X) / (1 / sigma2 + rowSums(psi)))
    s2 = as.vector(1 / (1 / sigma2 + rowSums(psi)))

    res = list(
      m = m,
      s2 = s2,
      psi = psi
    )
    return(res)
  }

  return(update_params)
}

initialize_params = function(K, n, sigma2) {
  res = list(
    m = rnorm(K, sd = sqrt(sigma2)),
    s2 = rep(1, K),
    psi = matrix(1 / K, nrow = K, ncol = n)
  )
  return(res)
}


CAVI = function(K, sigma2, n) {
  dat = simulate_data(K = K, sigma2 = sigma2, n = n)
  elbo = make_elbo(x = dat$x, sigma2 = dat$sigma2, K = dat$K)
  update = make_update_params(x = dat$x, sigma2 = dat$sigma2, K = dat$K)
  params = initialize_params(K = K, n = n, sigma2 = dat$sigma2)
  res = list(
    data = dat,
    params = params,
    elbo = elbo,
    update = update
  )
}

plot_mixture_fit = function(x, c_i, mu, m, s2, psi, mean_estimates_only = TRUE) {
  K = length(m)
  x_df = data.frame(x = x, Component = c_i)
  mu_df = data.frame(mu = mu, Component = factor(1:K))
  x_est_df = data.frame(x = x, Component = factor(apply(psi, 2, which.max)))
  m_df = data.frame(m = m,
                    y = seq(0, 0.01, length = K),
                    s2 = s2,
                    ub = m + 1.96 * s2,
                    lb = m - 1.96 * s2)
  nn = FNN::get.knnx(mu_df$mu, m_df$m, k = 1)
  m_df$Component = factor(nn$nn.index, levels = c(1:K, "Unclear"))
  dupes = names(table(m_df$Component))[table(m_df$Component) > 1]
  if (length(dupes) > 0) {
    m_df$Component[m_df$Component %in% dupes] = "Unclear"
  }
  gg_true = ggplot(x_df, aes(x = x, col = Component, fill = Component)) +
    geom_density() +
    scale_color_viridis_d() +
    scale_fill_viridis_d(alpha = 0.2) +
    geom_vline(data = mu_df, aes(xintercept = mu, col = Component), lwd = 1) +
    theme_bw()
  if (mean_estimates_only) {
    gg_res = gg_true +
      geom_point(data = m_df, aes(x = m, y = y, col = Component)) +
      geom_errorbarh(data = m_df,
                     aes(y = y,
                         xmax = ub,
                         xmin = lb,
                         col = Component,
                         height = 0.05),
                     inherit.aes = FALSE)

  } else {

    gg_est = ggplot(x_est_df, aes(x = x, col = Component, fill = Component)) +
      geom_density() +
      scale_color_viridis_d() +
      scale_fill_viridis_d(alpha = 0.2) +
      geom_vline(data = m_df, aes(xintercept = m, col = Component), lwd = 2) +
      theme_bw() +
      ggtitle("Gaussian mixture - estimate")

    cowplot::plot_grid(gg_true, gg_est)
  }

  return(gg_res)
}

plot_component_assignment = function(x, mu, m, s2, psi) {
  K = length(m)
  mu_df = data.frame(mu = mu, Component = factor(1:K))
  x_est_df = data.frame(x = x, Component = factor(apply(psi, 2, which.max)))
  m_df = data.frame(m = m,
                    y = seq(0, 0.01, length = K),
                    s2 = s2,
                    ub = m + 1.96 * s2,
                    lb = m - 1.96 * s2,
                    Component = factor(1:K))
  nn = FNN::get.knnx(mu_df$mu, m_df$m, k = 1)
  recode_vals = as.character(nn$nn.index)
  names(recode_vals) = 1:K
  m_df = m_df %>%
    mutate(Component = recode(Component, !!!recode_vals),
           Component = factor(Component, levels = as.character(1:K)))
  x_est_df = x_est_df %>%
    mutate(Component = recode(Component, !!!recode_vals),
           Component = factor(Component, levels = as.character(1:K)))

  gg_est = ggplot(x_est_df, aes(x = x, col = Component, fill = Component)) +
    geom_density() +
    scale_color_viridis_d(drop = FALSE) +
    scale_fill_viridis_d(alpha = 0.2, drop = FALSE) +
    geom_vline(data = mu_df, aes(xintercept = mu, col = Component), lwd = 1) +
    geom_point(data = m_df, aes(x = m, y = y, col = Component)) +
    geom_errorbarh(data = m_df,
                   aes(y = y,
                       xmax = ub,
                       xmin = lb,
                       col = Component,
                       height = 0.05),
                   inherit.aes = FALSE) +
    theme_bw()

  return(gg_est)
}

plot_elbo = function(elbo_df, current_run) {
  if (nrow(elbo_df) < 2) {
    return(NULL)
  }
  gg_res = ggplot(elbo_df, aes(x = iter, y = elbo, col = run_id, group = run_id)) +
    geom_line() +
    theme_bw() +
    xlim(0, max(10, elbo_df$iter)) +
    scale_color_viridis_d("Run")
  return(gg_res)
}

ui = fluidPage(
  titlePanel("Variational Inference for Univariate Gaussian Mixture Model"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        inputId = "K",
        label = "Number of Clusters:",
        value = 3,
        step = 1,
        min = 1,
        max = 10
      ),
      numericInput(
        inputId = "sigma2",
        label = "Prior variance:",
        value = 10,
        min = 0.1,
        max = 100
      ),
      sliderInput(
        inputId = "n",
        label = "Number of observations:",
        value = 30,
        min = 10,
        max = 1000,
        step = 10
      ),
      checkboxInput(
        inputId = "display_component_assignment",
        label = "Display component assignment?"
      ),
      actionButton("simulate_data", label = "Simulate Data", width = "100%"),
      actionButton("update", label = "Update params", width = "100%"),
      actionButton("reset", label = "Reset Fitting Process", width = "100%")
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.simulate_data > 0",
        h3("True Data"),
        plotOutput("plot_data"),
        conditionalPanel(
          condition = "input.display_component_assignment",
          h3("Estimated Densities"),
          plotOutput("plot_assigned_components")
        ),
        conditionalPanel(
          condition = "input.update > 0",
          h3("ELBO"),
          h6(verbatimTextOutput("elbo_diff")),
          plotOutput("plot_elbo")
        )
      )
    ),
  )
)

server = function(input, output, session) {

  v = reactiveValues(data = NULL,
                     params = NULL,
                     elbo = NULL,
                     update = NULL,
                     n_iter = 0,
                     n_run = 1,
                     elbo_out = numeric(),
                     elbo_runs = data.frame(iter = numeric(), elbo = numeric(), run_id = numeric()))

  observeEvent(input$simulate_data, {
    cavi_out = CAVI(input$K, input$sigma2, input$n)
    v$data = cavi_out$data
    v$params = cavi_out$params
    v$elbo = cavi_out$elbo
    v$update = cavi_out$update
    v$n_iter = 0
    v$elbo_out = v$elbo(v$params$m, v$params$s2, v$params$psi)
    v$n_run = 1
    v$elbo_runs = data.frame(iter = numeric(), elbo = numeric(), run_id = numeric())
  })

  observeEvent(input$update, {
    v$params = v$update(v$params$m, v$params$s2, v$params$psi)
    v$n_iter = v$n_iter + 1
    v$elbo_out[v$n_iter + 1] = v$elbo(v$params$m, v$params$s2, v$params$psi)
  })

  observeEvent(input$reset, {
    v$params = initialize_params(input$K, input$n, input$sigma2)
    elbo_df = data.frame(iter = 0:v$n_iter, elbo = v$elbo_out, run_id = v$n_run)
    v$n_run = v$n_run + 1
    v$elbo_runs = dplyr::bind_rows(v$elbo_runs, elbo_df)
    v$n_iter = 0
    v$elbo_out = v$elbo(v$params$m, v$params$s2, v$params$psi)
  })

  output$plot_data = renderPlot({
    if (is.null(v$data)) return()
    plot_mixture_fit(v$data$x,
                     v$data$c_i,
                     v$data$mu,
                     v$params$m,
                     v$params$s2,
                     v$params$psi)
  })

  output$plot_assigned_components = renderPlot({
    if (is.null(v$data)) return()
    plot_component_assignment(x = v$data$x,
                              mu = v$data$mu,
                              m = v$params$m,
                              s2 = v$params$s2,
                              psi = v$params$psi)
  })

  output$plot_elbo = renderPlot({
    if (is.null(v$data)) return()
    elbo_df = v$elbo_runs %>%
      bind_rows(data.frame(iter = 0:v$n_iter, elbo = v$elbo_out, run_id = v$n_run)) %>%
      mutate(run_id = factor(run_id))

    plot_elbo(elbo_df = elbo_df, current_run = v$n_run)
  })

  output$elbo_diff = renderText({
    if (v$n_iter == 0) {
      paste0("No ELBO increase yet.")
    } else {
      paste0("Current iteration: ", v$n_iter, ". ELBO increased by ", round(v$elbo_out[v$n_iter + 1] - v$elbo_out[v$n_iter], digits = 6))
    }
  })

}

shinyApp(ui, server)