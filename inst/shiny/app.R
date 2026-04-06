library(shiny)
library(dtasamplesize)

ui <- fluidPage(
  titlePanel("dtasamplesize: Sample Size for DTA Studies"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "method_tab",
        tabPanel("Buderer MC",
          numericInput("bud_se", "Sensitivity", 0.85, 0.01, 0.99, 0.01),
          numericInput("bud_d", "Precision (d)", 0.07, 0.01, 0.20, 0.01),
          numericInput("bud_B", "MC reps", 1000, 100, 10000, 100),
          actionButton("btn_buderer", "Calculate")
        ),
        tabPanel("BAM",
          numericInput("bam_a_se", "Prior Se alpha", 17, 1, 100),
          numericInput("bam_b_se", "Prior Se beta", 3, 1, 100),
          numericInput("bam_delta", "Target width", 0.14, 0.01, 0.50, 0.01),
          numericInput("bam_B", "MC reps", 1000, 100, 10000, 100),
          actionButton("btn_bam", "Calculate")
        ),
        tabPanel("Joint",
          numericInput("jnt_se", "Sensitivity", 0.85, 0.01, 0.99, 0.01),
          numericInput("jnt_sp", "Specificity", 0.90, 0.01, 0.99, 0.01),
          numericInput("jnt_auc", "AUC", 0.80, 0.50, 0.99, 0.01),
          numericInput("jnt_prev", "Prevalence", 0.20, 0.01, 0.99, 0.01),
          numericInput("jnt_B", "MC reps", 1000, 100, 10000, 100),
          actionButton("btn_joint", "Calculate")
        ),
        tabPanel("Imperfect Ref",
          numericInput("imp_se_ref", "Ref Se", 0.90, 0.50, 1.00, 0.01),
          numericInput("imp_sp_ref", "Ref Sp", 0.95, 0.50, 1.00, 0.01),
          numericInput("imp_loss", "Loss rate", 0.10, 0.00, 0.50, 0.01),
          actionButton("btn_imperfect", "Calculate")
        )
      )
    ),
    mainPanel(
      h4("Results"),
      verbatimTextOutput("results_output")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveVal(NULL)

  observeEvent(input$btn_buderer, {
    rv(mc_validate_buderer(
      Se = input$bud_se, d = input$bud_d, B = input$bud_B, seed = 2026
    ))
  })

  observeEvent(input$btn_bam, {
    rv(bam_sample_size(
      prior_se = c(input$bam_a_se, input$bam_b_se),
      delta_se = input$bam_delta,
      B = input$bam_B, seed = 2026
    ))
  })

  observeEvent(input$btn_joint, {
    rv(joint_sample_size(
      Se = input$jnt_se, Sp = input$jnt_sp, AUC = input$jnt_auc,
      prev = input$jnt_prev, B = input$jnt_B, seed = 2026,
      N_range = seq(100, 800, by = 10)
    ))
  })

  observeEvent(input$btn_imperfect, {
    rv(ss_imperfect_ref(
      Se_ref = input$imp_se_ref, Sp_ref = input$imp_sp_ref,
      loss_rate = input$imp_loss, B = 0
    ))
  })

  output$results_output <- renderPrint({
    req(rv())
    print(rv())
  })
}

shinyApp(ui, server)
