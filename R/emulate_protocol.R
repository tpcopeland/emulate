#' Generate a target trial protocol table
#'
#' Creates the 7-component target trial protocol table as described by
#' Hernan and Robins (2016, 2020). This table is the cornerstone of any
#' target trial emulation: it forces the researcher to explicitly define
#' every element of the hypothetical randomized trial before touching the
#' data.
#'
#' This is a \strong{standalone function} -- it does not require a \code{emulate}
#' object or any data. You can (and should) write the protocol before
#' starting the analysis pipeline.
#'
#' @section The 7 components:
#' Hernan and Robins define seven components that fully specify a target
#' trial:
#' \enumerate{
#'   \item \strong{Eligibility criteria} -- Who is eligible to be enrolled?
#'     (e.g., "Adults aged 40-80 with newly diagnosed hypertension and no
#'     prior cardiovascular events.")
#'   \item \strong{Treatment strategies} -- What interventions are being
#'     compared? (e.g., "Initiate Drug A within 30 days of diagnosis vs.
#'     no treatment.")
#'   \item \strong{Treatment assignment} -- How is treatment assigned?
#'     (e.g., "At each eligible time point, individuals are assigned based
#'     on their observed treatment decision.")
#'   \item \strong{Start of follow-up} -- When does follow-up begin?
#'     (e.g., "At the time of treatment eligibility (time zero).")
#'   \item \strong{Outcome} -- What is the primary endpoint?
#'     (e.g., "First myocardial infarction or death from cardiovascular
#'     causes, whichever occurs first.")
#'   \item \strong{Causal contrast} -- What is the estimand?
#'     (e.g., "Per-protocol effect: difference in 5-year risk of the
#'     outcome between treatment strategies.")
#'   \item \strong{Analysis plan} -- How will the data be analyzed?
#'     (e.g., "Sequential trials with clone-censor-weight, stabilized
#'     IPTW, pooled logistic regression, and G-formula standardization.")
#' }
#'
#' @param eligibility A single character string describing the eligibility
#'   criteria for the target trial.
#' @param treatment A single character string describing the treatment
#'   strategies (interventions) being compared.
#' @param assignment A single character string describing how treatment is
#'   assigned in the emulation (e.g., at each eligible time point based on
#'   observed treatment decisions).
#' @param followup_start A single character string describing when follow-up
#'   begins (i.e., the definition of "time zero").
#' @param outcome A single character string defining the primary outcome
#'   (endpoint) of the trial.
#' @param causal_contrast A single character string describing the causal
#'   contrast or estimand (e.g., intention-to-treat effect, per-protocol
#'   effect, risk difference at a specific time horizon).
#' @param analysis A single character string describing the statistical
#'   analysis plan (e.g., "clone-censor-weight sequential trials, stabilized
#'   IPTW, pooled logistic regression, G-formula predictions").
#' @param format A character string specifying the output format. One of:
#'   \describe{
#'     \item{\code{"display"}}{(Default) Prints the protocol table to the
#'       console via \code{message()}.}
#'     \item{\code{"csv"}}{Exports the table to a CSV file. Requires the
#'       \code{export} argument.}
#'     \item{\code{"excel"}}{Exports the table to an Excel (.xlsx) file with
#'       formatting. Requires the \pkg{openxlsx} package and the \code{export}
#'       argument.}
#'     \item{\code{"latex"}}{Exports the table as a LaTeX \code{tabular}
#'       environment. Requires the \code{export} argument.}
#'   }
#' @param export A character string giving the file path for CSV, Excel, or
#'   LaTeX export. Required when \code{format} is not \code{"display"};
#'   ignored when \code{format = "display"}. Default is \code{NULL}.
#' @param title An optional character string used as the table title. Defaults
#'   to \code{"Target Trial Protocol"}.
#'
#' @return A \code{data.frame} with two columns (\code{Component} and
#'   \code{Specification}) and seven rows, returned invisibly. Each row
#'   corresponds to one of the seven protocol components.
#'
#' @references
#' Hernan MA, Robins JM (2016). "Using Big Data to Emulate a Target Trial
#' When a Randomized Trial Is Not Available." \emph{American Journal of
#' Epidemiology}, 183(8), 758-764.
#'
#' Hernan MA, Robins JM (2020). \emph{Causal Inference: What If}.
#' Chapman & Hall/CRC.
#'
#' @examples
#' # Define and display a protocol for a statin trial emulation
#' protocol <- emulate_protocol(
#'   eligibility     = "Adults 40-80, newly diagnosed hyperlipidemia, no prior CVD",
#'   treatment       = "Initiate statin within 6 months vs. no statin",
#'   assignment      = "Observational: assigned by physician prescribing decision",
#'   followup_start  = "Date of first eligible visit",
#'   outcome         = "Composite of MI, stroke, or CV death",
#'   causal_contrast = "Per-protocol effect on 5-year risk",
#'   analysis        = "Sequential trials, IPTW, pooled logistic regression"
#' )
#'
#' # The result is a data.frame
#' str(protocol)
#'
#' \donttest{
#' # Export to CSV
#' tmp <- tempfile(fileext = ".csv")
#' emulate_protocol(
#'   eligibility     = "Adults 40-80, newly diagnosed hyperlipidemia",
#'   treatment       = "Statin vs. no statin",
#'   assignment      = "Observational",
#'   followup_start  = "First eligible visit",
#'   outcome         = "MI, stroke, or CV death",
#'   causal_contrast = "PP effect on 5-year risk",
#'   analysis        = "IPTW + pooled logistic regression",
#'   format = "csv", export = tmp
#' )
#' }
#'
#' @export
emulate_protocol <- function(eligibility, treatment, assignment,
                         followup_start, outcome, causal_contrast,
                         analysis, format = "display", export = NULL,
                         title = NULL) {
  # All 7 components required
  components <- c("Eligibility criteria" = eligibility,
                   "Treatment strategies" = treatment,
                   "Treatment assignment" = assignment,
                   "Start of follow-up" = followup_start,
                   "Outcome" = outcome,
                   "Causal contrast" = causal_contrast,
                   "Analysis plan" = analysis)

  protocol <- data.frame(
    Component = names(components),
    Specification = unname(components),
    stringsAsFactors = FALSE
  )

  if (is.null(title)) title <- "Target Trial Protocol"

  if (format == "display") {
    message(title)
    message(strrep("=", 70))
    for (i in seq_len(nrow(protocol))) {
      message(sprintf("%-25s %s", protocol$Component[i],
                      protocol$Specification[i]))
    }
    message(strrep("=", 70))
  } else if (format == "csv") {
    if (is.null(export)) stop("export path required for csv format", call. = FALSE)
    utils::write.csv(protocol, export, row.names = FALSE)
    message("Protocol exported to: ", export)
  } else if (format == "excel") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("openxlsx package required for Excel export", call. = FALSE)
    }
    if (is.null(export)) stop("export path required for excel format", call. = FALSE)

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Protocol")
    openxlsx::writeData(wb, "Protocol", protocol, startRow = 2)
    openxlsx::writeData(wb, "Protocol", title, startRow = 1, startCol = 1)

    # Formatting
    bold_style <- openxlsx::createStyle(textDecoration = "bold")
    wrap_style <- openxlsx::createStyle(wrapText = TRUE)
    openxlsx::addStyle(wb, "Protocol", bold_style, rows = 2, cols = 1:2)
    openxlsx::addStyle(wb, "Protocol", wrap_style,
                        rows = 3:9, cols = 2, gridExpand = TRUE)
    openxlsx::setColWidths(wb, "Protocol", cols = 1, widths = 28)
    openxlsx::setColWidths(wb, "Protocol", cols = 2, widths = 60)

    openxlsx::saveWorkbook(wb, export, overwrite = TRUE)
    message("Protocol exported to: ", export)
  } else if (format == "latex") {
    if (is.null(export)) stop("export path required for latex format", call. = FALSE)

    lines <- c(
      "\\begin{table}[htbp]",
      paste0("\\caption{", title, "}"),
      "\\begin{tabular}{lp{10cm}}",
      "\\hline",
      "\\textbf{Component} & \\textbf{Specification} \\\\",
      "\\hline"
    )
    for (i in seq_len(nrow(protocol))) {
      lines <- c(lines, paste0(protocol$Component[i], " & ",
                                protocol$Specification[i], " \\\\"))
    }
    lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
    writeLines(lines, export)
    message("Protocol exported to: ", export)
  }

  invisible(protocol)
}
