TCGAanalyze_survival_override <- function(data,
                                          path,
                                          clusterCol = NULL,
                                          legend = "Legend",
                                          labels = NULL,
                                          risk.table = TRUE,
                                          xlim = NULL,
                                          main = "Kaplan-Meier Overall Survival Curves",
                                          ylab = "Probability of survival",
                                          xlab = "Time since diagnosis (days)",
                                          genename = "",
                                          color = NULL,
                                          height = 8,
                                          width = 12,
                                          dpi = 300,
                                          pvalue = TRUE,
                                          conf.int = TRUE,
                                          ...) {
  .e <- environment()


  if (!all(
    c(
      "vital_status",
      "days_to_death",
      "days_to_last_follow_up"
    ) %in% colnames(data)
  )) {
    stop(
      "Columns vital_status, days_to_death and  days_to_last_follow_up should be in data frame"
    )
  }

  if (is.null(color)) {
    color <- rainbow(length(unique(data[, clusterCol])))
  }

  group <- NULL
  if (is.null(clusterCol)) {
    stop("Please provide the clusterCol argument")
  } else if (length(unique(data[, clusterCol])) == 1) {
    stop(
      paste0(
        "Sorry, but I'm expecting at least two groups\n",
        "  Only this group found: ",
        unique(data[, clusterCol])
      )
    )
  }
  notDead <- is.na(data$days_to_death)

  if (any(notDead == TRUE)) {
    data[notDead, "days_to_death"] <-
      data[notDead, "days_to_last_follow_up"]
  }
  if (length(data[which((data[, "days_to_death"] < 0) == T), "sample"]) > 0 &
    "sample" %in% colnames(data)) {
    message(
      "Incosistencies in the data were found. The following samples have a negative days_to_death value:"
    )
    message(paste(data[which((data[, "days_to_death"] < 0) == T), "sample"], collapse = ", "))
  }
  if (any(is.na(data[, "days_to_death"])) &
    "sample" %in% colnames(data)) {
    message(
      "Incosistencies in the data were found. The following samples have a NA days_to_death value:"
    )
    message(paste(data[is.na(data[, "days_to_death"]), "sample"], collapse = ", "))
  }

  # create a column to be used with survival package, info need
  # to be TRUE(DEAD)/FALSE (ALIVE)
  data$s <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE)

  # Column with groups
  data$type <- as.factor(data[, clusterCol])
  data <- data[, c("days_to_death", "s", "type")]
  # create the formula for survival analysis
  f.m <-
    formula(survival::Surv(as.numeric(data$days_to_death), event = data$s) ~ data$type)
  fit <- do.call(survival::survfit, list(formula = f.m, data = data))

  label.add.n <- function(x) {
    na.idx <- is.na(data[, "days_to_death"])
    negative.idx <- data[, "days_to_death"] < 0
    idx <- !(na.idx | negative.idx)
    return(paste0(
      x, " (n = ",
      sum(data[idx, "type"] == x), ")"
    ))
  }

  if (is.null(labels)) {
    d <- survminer::surv_summary(fit, data = data)
    order <-
      unname(sapply(levels(d$strata), function(x) {
        unlist(str_split(x, "="))[2]
      }))
    labels <- sapply(order, label.add.n)
  }
  if (length(xlim) == 1) {
    xlim <- c(0, xlim)
  }
  suppressWarnings({
    surv <- survminer::ggsurvplot(
      fit,
      # survfit object with calculated statistics.
      risk.table = risk.table,
      # show risk table.
      pval = pvalue,
      # show p-value of log-rank test.
      conf.int = conf.int,
      # show confidence intervals for point estimaes of survival curves.
      xlim = xlim,
      # present narrower X axis, but not affect survival estimates.
      main = main,
      # Title
      xlab = xlab,
      # customize X axis label.
      legend.title = legend,
      # Legend title
      legend.labs = labels,
      # change legend labels.
      palette =  color,
      # custom color palettes.
      ...
    )
  })


  dir.create(paste0(path, "km/"), showWarnings = FALSE)
  if (surv_pvalue(fit)$pval > 0.05) {
    dir.create(paste0(path, "km/non_sig/"), showWarnings = FALSE)
    filename <- paste0(path, "km/non_sig/", surv_pvalue(fit)$pval, "-", genename, ".pdf")
  } else {
    dir.create(paste0(path, "km/sig/"), showWarnings = FALSE)
    filename <- paste0(path, "km/sig/", surv_pvalue(fit)$pval, "-", genename, ".pdf")
  }

  ggsave(
    surv$plot,
    filename = filename,
    width = width,
    height = height,
    dpi = dpi
  )
  message(paste0("File saved as: ", filename))
  # if (risk.table) {
  #     g1 <- ggplotGrob(surv$plot)
  #     g2 <- ggplotGrob(surv$table)
  #     min_ncol <- min(ncol(g2), ncol(g1))
  #     g <-
  #         gridExtra::gtable_rbind(g1[, 1:min_ncol], g2[, 1:min_ncol], size = "last")
  #     ggsave(
  #         g,
  #         filename = filename,
  #         width = width,
  #         height = height,
  #         dpi = dpi
  #     )
  # }
}


TCGAanalyze_SurvivalKM_override <- function(clinical_patient,
                                            dataGE,
                                            Genelist,
                                            path,
                                            Survresult = FALSE,
                                            ThreshTop = 0.67,
                                            ThreshDown = 0.33,
                                            p.cut = 0.05,
                                            group1,
                                            group2) {

  # check_package("survival")
  # Check which genes we really have in the matrix
  Genelist <- intersect(rownames(dataGE), Genelist)

  # Split gene expression matrix btw the groups
  dataCancer <- dataGE[Genelist, group2, drop = FALSE]
  dataNormal <- dataGE[Genelist, group1, drop = FALSE]

  colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)

  cfu <- clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% substr(colnames(dataCancer), 1, 12), ]

  if ("days_to_last_followup" %in% colnames(cfu)) {
    colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <-
      "days_to_last_follow_up"
  }
  cfu <-
    as.data.frame(subset(
      cfu,
      select = c(
        "bcr_patient_barcode",
        "days_to_death",
        "days_to_last_follow_up",
        "vital_status"
      )
    ))

  # Set alive death to inf
  if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0) {
    cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <-
      "-Inf"
  }

  # Set dead follow up to inf
  if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0) {
    cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <-
      "-Inf"
  }

  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]

  followUpLevel <- FALSE

  # FC_FDR_table_mRNA
  tabSurv_Matrix <-
    matrix(0, nrow(as.matrix(rownames(dataNormal))), 8)
  colnames(tabSurv_Matrix) <- c(
    "mRNA",
    "pvalue",
    "Cancer Deaths",
    "Cancer Deaths with Top",
    "Cancer Deaths with Down",
    "Mean Tumor Top",
    "Mean Tumor Down",
    "Mean Normal"
  )

  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)

  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <-
    as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"] # mod1

  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]

  cfu_complete <- cfu
  ngenes <- nrow(as.matrix(rownames(dataNormal)))


  # Evaluate each gene
  for (i in 1:nrow(as.matrix(rownames(dataNormal)))) {
    cat(paste0((ngenes - i), "."))
    mRNAselected <- as.matrix(rownames(dataNormal))[i]
    mRNAselected_values <-
      dataCancer[rownames(dataCancer) == mRNAselected, ]
    mRNAselected_values_normal <-
      dataNormal[rownames(dataNormal) == mRNAselected, ]
    if (all(mRNAselected_values == 0)) {
      next
    } # All genes are 0
    tabSurv_Matrix[i, "mRNA"] <- mRNAselected


    # Get Thresh values for cancer expression
    mRNAselected_values_ordered <-
      sort(mRNAselected_values, decreasing = TRUE)
    mRNAselected_values_ordered_top <-
      as.numeric(quantile(as.numeric(mRNAselected_values_ordered), ThreshTop)[1])
    mRNAselected_values_ordered_down <-
      as.numeric(quantile(as.numeric(mRNAselected_values_ordered), ThreshDown)[1])

    mRNAselected_values_newvector <- mRNAselected_values


    if (!is.na(mRNAselected_values_ordered_top)) {
      # How many samples do we have
      numberOfSamples <- length(mRNAselected_values_ordered)

      # High group (above ThreshTop)
      lastelementTOP <-
        max(which(
          mRNAselected_values_ordered > mRNAselected_values_ordered_top
        ))

      # Low group (below ThreshDown)
      firstelementDOWN <-
        min(
          which(
            mRNAselected_values_ordered <= mRNAselected_values_ordered_down
          )
        )

      samples_top_mRNA_selected <-
        names(mRNAselected_values_ordered[1:lastelementTOP])
      samples_down_mRNA_selected <-
        names(mRNAselected_values_ordered[firstelementDOWN:numberOfSamples])

      # Which samples are in the intermediate group (above ThreshLow and below ThreshTop)
      samples_UNCHANGED_mRNA_selected <-
        names(mRNAselected_values_newvector[which((mRNAselected_values_newvector) > mRNAselected_values_ordered_down &
          mRNAselected_values_newvector < mRNAselected_values_ordered_top)])

      cfu_onlyTOP <-
        cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_top_mRNA_selected, ]
      cfu_onlyDOWN <-
        cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_down_mRNA_selected, ]
      cfu_onlyUNCHANGED <-
        cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, ]

      cfu_ordered <- NULL
      cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
      cfu <- cfu_ordered

      ttime <- as.numeric(cfu[, "days_to_death"])

      sum(status <- ttime > 0) # morti
      deads_complete <- sum(status <- ttime > 0)

      ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
      deads_top <- sum(ttime_only_top > 0)

      if (dim(cfu_onlyDOWN)[1] >= 1) {
        ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
        deads_down <- sum(ttime_only_down > 0)
      } else {
        deads_down <- 0
      }

      tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
      tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
      tabSurv_Matrix[i, "Cancer Deaths with Down"] <-
        deads_down
      tabSurv_Matrix[i, "Mean Normal"] <-
        mean(as.numeric(mRNAselected_values_normal))
      dataCancer_onlyTop_sample <-
        dataCancer[, samples_top_mRNA_selected, drop = FALSE]
      dataCancer_onlyTop_sample_mRNASelected <-
        dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == mRNAselected, ]
      dataCancer_onlyDown_sample <-
        dataCancer[, samples_down_mRNA_selected, drop = FALSE]
      dataCancer_onlyDown_sample_mRNASelected <-
        dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == mRNAselected, ]
      tabSurv_Matrix[i, "Mean Tumor Top"] <-
        mean(as.numeric(dataCancer_onlyTop_sample_mRNASelected))
      tabSurv_Matrix[i, "Mean Tumor Down"] <-
        mean(as.numeric(dataCancer_onlyDown_sample_mRNASelected))

      ttime[!status] <-
        as.numeric(cfu[!status, "days_to_last_follow_up"])
      ttime[which(ttime == -Inf)] <- 0

      ttime <- survival::Surv(ttime, status)
      rownames(ttime) <- rownames(cfu)
      legendHigh <- paste(mRNAselected, "High")
      legendLow <- paste(mRNAselected, "Low")

      tabSurv_pvalue <- tryCatch(
        {
          tabSurv <-
            survival::survdiff(ttime ~ c(rep(
              "top", nrow(cfu_onlyTOP)
            ), rep(
              "down", nrow(cfu_onlyDOWN)
            )))
          tabSurv_chis <- unlist(tabSurv)$chisq
          tabSurv_pvalue <-
            as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
        },
        error = function(e) {
          return(Inf)
        }
      )
      tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue

      if (Survresult == TRUE) {
        if (tabSurv_pvalue > 0.05) {
          dir.create(paste0(path, "km/non_sig/"), showWarnings = FALSE)
          pdf(paste0(path, "km/non_sig/", format(tabSurv_pvalue, scientific = FALSE), "-", mRNAselected, ".pdf"))
        } else {
          dir.create(paste0(path, "km/sig/"), showWarnings = FALSE)
          pdf(paste0(path, "km/sig/", format(tabSurv_pvalue, scientific = FALSE), "-", mRNAselected, ".pdf"))
        }


        titlePlot <-
          paste(
            "Kaplan-Meier Survival analysis, pvalue=",
            tabSurv_pvalue
          )
        plot(
          survival::survfit(ttime ~ c(
            rep("low", nrow(cfu_onlyTOP)), rep("high", nrow(cfu_onlyDOWN))
          )),
          col = c("green", "red"),
          main = titlePlot,
          xlab = "Days",
          ylab = "Survival"
        )
        legend(
          100,
          1,
          legend = c(paste0(mRNAselected, " Low"), paste0(mRNAselected, " High")),
          col = c("green", "red"),
          text.col = c("green", "red"),
          pch = 15
        )

        dev.off()


        # survminer::ggsurvplot(
        #     survival::survfit(ttime ~ c(rep("low", nrow(cfu_onlyTOP)), rep("high", nrow(cfu_onlyDOWN)))),
        #     data=cfu,
        #     legend.labs=c(legendLow, legendHigh)
        # )
        # pdf(paste0("/cwd/data/km/", mRNAselected,".pdf"), width=14, height=7)

        print(tabSurv)
      }
    } # end if
  } # end for

  tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0

  tabSurvKM <- tabSurv_Matrix

  # Filtering by selected pvalue < 0.01
  tabSurvKM <- tabSurvKM[tabSurvKM$mRNA != 0, ]
  tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < p.cut, ]
  tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
  rownames(tabSurvKM) <- tabSurvKM$mRNA
  tabSurvKM <- tabSurvKM[, -1]
  tabSurvKM <-
    tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE), ]

  colnames(tabSurvKM) <-
    gsub("Cancer", "Group2", colnames(tabSurvKM))
  colnames(tabSurvKM) <-
    gsub("Tumor", "Group2", colnames(tabSurvKM))
  colnames(tabSurvKM) <-
    gsub("Normal", "Group1", colnames(tabSurvKM))


  return(tabSurvKM)
}
