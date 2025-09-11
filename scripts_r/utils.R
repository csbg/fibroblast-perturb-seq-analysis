library(tidyverse)
library(openxlsx)
library(fs)


# Exporting functions -----------------------------------------------------

#' Save a plot with sensible defaults.
#'
#' @param filename Filename, will be saved in subfolder `plots/`. May contain
#'   additional subfolders, which are possibly created. If `NULL`, exit without
#'   creating a plot.
#' @param type Type of image file.
#' @param plot If `NULL`, save the `last_plot()` via `ggsave()`. Otherwise, save
#'   the graphics object in `plot` via the `png(); print(); dev.off()` workflow.
#' @param width Width in mm.
#' @param height Height in mm.
#' @param dip Resolution in dots per inch.
#' @param crop If `TRUE`, remove white margins from the saved plot.
#' @param bg Background color of saved image.
#' @param ... Other parameters passed to `ggsave()` or `png()`.
#'
#' @return The filename, invisibly.
ggsave_default <- function(filename,
                           type = "png",
                           plot = NULL,
                           width = 297,
                           height = 210,
                           dpi = 300,
                           crop = TRUE,
                           bg = "white",
                           ...) {
  if (is.null(filename))
    return()

  filename <- stringr::str_glue("plots/{filename}.{type}")
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()

  if (is.null(plot)) {
    ggplot2::ggsave(filename, dpi = dpi, units = "mm", limitsize = FALSE,
                    width = width, height = height, bg = bg, ...)
  } else {
    plot_args <- list(
      png = list(res = dpi, units = "mm")
    )
    rlang::exec(type, filename, !!!pluck(plot_args, type),
                width = width, height = height,  ...)
    print(plot)
    dev.off()
  }

  if (crop)
    knitr::plot_crop(filename)

  invisible(filename)
}



#' Save a (list of) data frame as well-formatted XLSX file.
#'
#' @param tables Data frames to be saved. A list of dataframes is saved into
#'   separate worksheets, whose names equal the list names. If the list is
#'   unnamed, worksheets will be named "Sheet1" etc.
#' @param filename Filename without extension; will be saved to `tables/`.
#' @param sheet_name Sheet name if a single unnamed table should be saved.
#'
#' @return Nothing.
save_table <- function(tables, filename, sheet_name = "Sheet1") {
  filename <- str_glue("tables/{filename}.xlsx")
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()
  wb <- createWorkbook()

  # ensure that tables is a named list
  if (inherits(tables, "list")) {
    if (is.null(names(tables)))
      tables <- set_names(tables, paste0("Sheet", seq_along(tables)))
  } else if (inherits(tables, "data.frame")) {
    tables <- list(tables) %>% set_names(sheet_name)
  } else {
    stop("'tables' must be a data frame or list of data frames.")
  }

  # populate Excel file with worksheets
  iwalk(
    tables,
    function(table, sheet_name) {
      addWorksheet(wb, sheet_name)
      writeData(
        wb,
        sheet_name,
        table,
        headerStyle = createStyle(textDecoration = "bold")
      )
      freezePane(wb, sheet_name, firstRow = TRUE)
      setColWidths(wb, sheet_name, 1:ncol(table), "auto")
    }
  )

  saveWorkbook(wb, filename, overwrite = TRUE)
}

# Logging -----------------------------------------------------------------

default_logger <- log4r::logger(threshold = "DEBUG")

debug <- function(..., .envir = parent.frame()) {
  log4r::debug(default_logger, glue::glue(..., .envir = .envir))
}

info <- function(..., .envir = parent.frame()) {
  log4r::info(default_logger, glue::glue(..., .envir = .envir))
}

warn <- function(..., .envir = parent.frame()) {
  log4r::warn(default_logger, glue::glue(..., .envir = .envir))
}



# Common definitions ------------------------------------------------------

CELL_TYPE_COLORS = c(
  "Inflammatory (injury response)" = "#b4c8eb",
  "Proliferative" = "#98d876",
  "Deactivated myofibroblasts" = "#37a53e",
  "Quiescent" = "#f68191",
  "Lamp1-fibroblasts" = "#874943",
  "Tissue-repair" = "#c8a9da",
  "Transitory" = "#8f56b1",
  "Stress-fiber fibroblasts" = "#f46e16",
  "Inflammatory & fibrotic" = "#3974bc",
  "Myofibroblasts" = "#fab162",
  "Proliferative myofibroblast" = "#d82326"
)



# Styling -----------------------------------------------------------------

BASE_TEXT_SIZE_MM = 1.76  # mm, corresponds to 5 pt, use in geom_text()
BASE_TEXT_SIZE_PT = 5 # pt, use in theme()
BASE_LINEWIDTH = 0.25 # pt
BASE_BOXPLOT_SIZE = 0.5

#' Common theme for publication-quality figures.
#'
#' This theme bases upon `theme_bw()` and ensures
#' - common line widths of `BASE_LINEWIDTH`
#' - common text sizes of `BASE_TEXT_SIZE_PT`
#' - a uniform plot margin of 1 mm
#' - a medium strip text, an empty strip background, and
#'   1 mm padding between strip text and panel
#'
#' @param rotate_x_labels If `TRUE`, rotate x-axis tick labels by 90°.
#' @param ... Other parameters passed to `theme_bw()`.
#'
#' @return A theme object.
theme_pub <- function(rotate_x_labels = FALSE, ...){
  res <-
    theme_bw(...) +
    theme(
      line = element_line(linewidth = BASE_LINEWIDTH),
      axis.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      axis.title = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.background = element_blank(),
      legend.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.title = element_text(size = BASE_TEXT_SIZE_PT),
      panel.border = element_rect(linewidth = BASE_LINEWIDTH * 2),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      strip.background = element_blank(),
      strip.text = element_text(
        color = "black",
        size = BASE_TEXT_SIZE_PT
      ),
      strip.text.x = element_text(margin = margin(b = 1, unit = "mm")),
      strip.text.y = element_text(margin = margin(l = 1, unit = "mm")),
      plot.title = element_text(
        size = BASE_TEXT_SIZE_PT,
        face = "plain",
        hjust = 0.5
      )
    )

  if (rotate_x_labels)
    res <-
      res +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  res
}
