# load libs, load config, initialise log

source("R/libs.R")

is_html <- knitr::is_html_output()
# is_pdf <- knitr::is_latex_output()
# is_word <- !is_html & !is_pdf

ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(text = element_text(size = 8))
ggplot2::theme_update(legend.position = "bottom")
# ggplot2::theme_update(legend.title = element_blank())
ggplot2::theme_update(axis.text.x = element_text(size = 8))
ggplot2::theme_update(axis.text.y = element_text(size = 8))

# g_cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Config - store of local OS file system
f_cfg <- file.path("./etc", "cfg.yml")
g_cfg <- config::get(file = f_cfg)
stopifnot("Config is null" = !is.null(g_cfg))


# message(Sys.time(), " Config read from ", f_cfg)

# Logs
f_log <- file.path("./logs", "log.txt")
log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")

# pars/effects of interest

g_fx <- c("b_r", "b_r1d", "b_r2d", "b_f")
g_mod4_pars <- c("a0", paste0("m_", 1:2), paste0("b_", 1:8))
g_mod4_qnt <- c(
  "eta_r_0", "eta_r_1", "eta_d_0", "eta_d_1", "eta_f_0", "eta_f_1"
  )



# alt_font <- flextable::fp_text_default(bold = TRUE, font.size = 12, font.family = "Reenie Beanie")


odds <- function(p) p / (1 - p)
expit <- function(x) plogis(x)