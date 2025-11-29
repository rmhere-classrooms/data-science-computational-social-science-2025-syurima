# III Rozprzestrzenianie się informacji w sieciach - dane rzeczywiste + ShinyApps

# Wybrana miara nr 5: aktywacja węzłów o średnich outdegree.
# Miara ta wybiera środkowe 5% węzłów pod względem wychodzących krawędzi - badamy, 
# jak informacja rozprzestrzenia się, gdy zaczynamy od mniej wpływowych węzłów 
# (można porównać do sytuacji z metody 1).

library(shiny)
library(bslib)
library("lubridate") 
library(igraph)

# Define UI ----
ui <- page_sidebar(
  title = "Information spread in emails - simulation",

  sidebar = sidebar(
    # Input: Slider for wij multiplier ----
    sliderInput(
      inputId = "wijMultiplier",
      label = "Mnożnik prawdopodobieńsstwa aktywacji:",
      min = 0.1,
      max = 2,
      value = 1
    ),
    # Input: Slider for number of iterations ----
    sliderInput(
      inputId = "iterations",
      label = "Liczba iteracji:",
      min = 1,
      max = 50,
      value = 10
    ),
    # Input: Number of simulation runs (fixed to 100 for averaging) ----
    numericInput(
      inputId = "runs",
      label = "Liczba uruchomień symulacji:",
      value = 10,
      min = 1,
      max = 100
    )
  ),
  # Output: plots illustrating information spread
  # X axis: iteration number
  # Y axis: number of activated nodes in that iteration
  # five different data series - one for each method of selecting initial nodes
  plotOutput(outputId = "spreadPlot")
)

# Define server logic required to draw a plot ----
server <- function(input, output) {
  # ------ Import data and create graph ------

  # Import email data - skipping first two lines, using only first two columns
  data <- read.csv2("https://bergplace.org/share/out.radoslaw_email_email", skip=2, sep= " ")[, 1:2]
  # create directed graph from data frame
  g <- graph_from_data_frame(data, directed = TRUE)

  # Simplify the graph - remove loops and multiple edges
  # After simplification, the graph should have 167 nodes and 5783 edges.
  g <- simplify(g)
  cat("Graph: nodes =", vcount(g), "edges =", ecount(g), "\n")

  # ------ Set edge weights (use original data frame) ------

  add_weights <- function(graph, df) {
    # cntij: number of emails sent from i to j
    cntij <- aggregate(list(cnt = rep(1, nrow(df))), 
      by = list(from = df[,1], to = df[,2]), FUN = sum)

    # cnti: total emails sent by each sender i
    cnti <- aggregate(cnt ~ from, data = cntij, FUN = sum)
    names(cnti) <- c("from", "cnt_total")

    # compute weights - wij = cntij / cnti
    cntij <- merge(cntij, cnti, by = "from")
    cntij$wij <- cntij$cnt / cntij$cnt_total

    # map weights to graph edges
    edge_df <- as.data.frame(ends(graph, E(graph)), stringsAsFactors = FALSE)
    names(edge_df) <- c("from", "to")
    key_edges <- paste(edge_df$from, edge_df$to, sep = ":::")
    key_cnts  <- paste(cntij$from, cntij$to, sep = ":::")
    edge_df$wij <- cntij$wij[match(key_edges, key_cnts)]

    # assign weights
    E(graph)$weight <- edge_df$wij
    
    return(graph)
  }
  g <- add_weights(g, data)
  cat("Edge weights added.\n")
  cat("Graph: nodes =", vcount(g), "edges =", ecount(g), "\n")


  # ------ Activate nodes (various methods) -----

  activate_high_outdegree <- function(graph, percentage) {
    # set all nodes activated to FALSE
    V(graph)$activated <- FALSE
    
    # calculate outdegree
    outdegree_scores <- degree(graph, mode = "out")
    
    # select top percentage of nodes based on outdegree
    num_to_activate <- ceiling(percentage * vcount(graph))
    top_nodes <- order(outdegree_scores, decreasing = TRUE)[1:num_to_activate]
    V(graph)[top_nodes]$activated <- TRUE
    
    return(graph)
  }

  activate_betweenness <- function(graph, percentage) {
    # set all nodes activated to FALSE
    V(graph)$activated <- FALSE
    
    # calculate betweenness centrality
    betweenness_scores <- betweenness(graph, directed = TRUE)
    
    # select top percentage of nodes based on betweenness
    num_to_activate <- ceiling(percentage * vcount(graph))
    top_nodes <- order(betweenness_scores, decreasing = TRUE)[1:num_to_activate]
    V(graph)[top_nodes]$activated <- TRUE
    
    return(graph)
  }

  activate_closeness <- function(graph, percentage) {
    # set all nodes activated to FALSE
    V(graph)$activated <- FALSE
    
    # calculate closeness centrality
    closeness_scores <- closeness(graph, mode = "out")
    
    # select top percentage of nodes based on closeness
    num_to_activate <- ceiling(percentage * vcount(graph))
    top_nodes <- order(closeness_scores, decreasing = TRUE)[1:num_to_activate]
    V(graph)[top_nodes]$activated <- TRUE
    
    return(graph)
  }

  activate_random <- function(graph, percentage) {
    # set all nodes activated to FALSE
    V(graph)$activated <- FALSE
    
    # randomly select given percentage of nodes to activate
    num_to_activate <- ceiling(percentage * vcount(graph))
    initial_activated <- sample(V(graph), num_to_activate)
    V(graph)[initial_activated]$activated <- TRUE
    
    return(graph)
  }

  activate_mid_outdegree <- function(graph, percentage) {
    # set all nodes activated to FALSE
    V(graph)$activated <- FALSE
    
    # calculate outdegree
    outdegree_scores <- degree(graph, mode = "out")
    
    # sort nodes by outdegree
    sorted_nodes <- order(outdegree_scores)
    
    # select middle percentage of nodes based on outdegree
    n <- vcount(graph)
    num_to_activate <- ceiling(percentage * n)
    start_index <- floor((n - num_to_activate) / 2) + 1
    end_index <- start_index + num_to_activate - 1
    
    middle_nodes <- sorted_nodes[start_index:end_index]
    V(graph)[middle_nodes]$activated <- TRUE
    
    return(graph)
  }

  # ------ Simulate info spread ------

  independent_cascade <- function(graph, wij_multiplier = 1, max_iter = 10) {
    total_activated <- which(V(graph)$activated)
    newly_activated <- total_activated # nodes that were activated in current iteration (to avoid more than one activation attempt) - start with initial nodes
    activated_counts <- numeric() # store count after initial activation
    activated_counts[1] <- length(total_activated)

    while (length(newly_activated) > 0 && length(activated_counts) < max_iter) {
      current_activated <- list()
      for (node in newly_activated) {
        neighbors <- neighbors(graph, node, mode = "out")
        for (neighbor in neighbors) {
          if (!V(graph)[neighbor]$activated) {
            wij <- E(graph)[.from(node) & .to(neighbor)]$weight
            if (is.na(wij)) {
              wij <- 0
            }
            # wij * wij_multiplier is the probability of activation
            prob <- wij * wij_multiplier
            if (prob >= 1 || runif(1) < prob) {
              V(graph)[neighbor]$activated <- TRUE
              current_activated <- c(current_activated, neighbor)
            }
          }
        }
      }
      newly_activated <- unique(current_activated)
      total_activated <- unique(c(total_activated, newly_activated))
      activated_counts[length(activated_counts) + 1] <- length(total_activated) # store count after this iteration
    }

    return(activated_counts)
  }

  # function to run the simulation 100 times for a graph and return average number of activated nodes
  run_simulation <- function(graph, activation_func, wij_multiplier, max_iter = 10, runs = 100) {
    # Store activated counts for each run and each iteration
    all_counts <- vector("list", runs)
    g_activated <- activation_func(graph, percentage = 0.05) # activate 5% nodes

    for (i in 1:runs) {
      g_activated_clone <- g_activated # clone graph to reset activations
      activated_counts <- independent_cascade(g_activated_clone, wij_multiplier, max_iter)
      # Pad activated_counts to max_iter length
      if (length(activated_counts) < max_iter) {
        activated_counts <- c(activated_counts, rep(activated_counts[length(activated_counts)], max_iter - length(activated_counts)))
      } else if (length(activated_counts) > max_iter) {
        activated_counts <- activated_counts[1:max_iter]
      }
      all_counts[[i]] <- activated_counts
    }

    # Calculate average activated nodes for each iteration
    avg_activated_per_iteration <- colMeans(do.call(rbind, all_counts))

    cat ("Simulation completed with wijMultiplier =", input$wijMultiplier, "and iterations =", input$iterations, "runs = ", runs,  "\n")
    return(avg_activated_per_iteration)
  }

  output$spreadPlot <- renderPlot({
    num_activated_random <- run_simulation(g, activate_random, input$wijMultiplier, max_iter = input$iterations, runs = input$runs)
    num_activated_high_outdegree <- run_simulation(g, activate_high_outdegree, input$wijMultiplier, max_iter = input$iterations, runs = input$runs)
    num_activated_betweenness <- run_simulation(g, activate_betweenness, input$wijMultiplier, max_iter = input$iterations, runs = input$runs)
    num_activated_closeness <- run_simulation(g, activate_closeness, input$wijMultiplier, max_iter = input$iterations, runs = input$runs)
    num_activated_mid_outdegree <- run_simulation(g, activate_mid_outdegree, input$wijMultiplier, max_iter = input$iterations, runs = input$runs)
    
    plot(1:input$iterations, num_activated_random, type = "o", col = "blue", ylim = c(0, vcount(g)),
         xlab = "Iteration", ylab = "Number of Activated Nodes", 
         main = "Information Spread Simulation")

    lines(1:input$iterations, num_activated_high_outdegree, type = "o", col = "red")
    lines(1:input$iterations, num_activated_betweenness, type = "o", col = "green")
    lines(1:input$iterations, num_activated_closeness, type = "o", col = "purple")
    lines(1:input$iterations, num_activated_mid_outdegree, type = "o", col = "orange")
    legend("topleft", legend = c("Random", "High Outdegree", "Betweenness", "Closeness", "Med Outdegree"),
           col = c("blue", "red", "green", "purple", "orange"), lty = 1, pch = 1)
  })
}

shinyApp(ui = ui, server = server)
