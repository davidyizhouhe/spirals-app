library(shiny)
library(plotrix)

#-------------------------------------------------------------------------------
# GENERAL FUNCTIONS
#-------------------------------------------------------------------------------

# convenience function to convert from degrees to radians
deg_to_rad <- function(degrees) {
  # each degree is pi/180 radians
  return(degrees * pi/180)
}

# can obtain a rational approximation to a given level of accuracy with a
# simple continued fraction approximation of number in list format ([a;b,c...])
rational_approx <- function(number, iterations = 20) {
  # begin with integer portion of original number
  step <- 1
  terms <- c(floor(number))
  # reduce number by integer part
  number <- number - terms[1]
  
  # to obtain next term in CF, simply take integer portion of reciprocal
  while(step < iterations && 
        abs(number - floor(number)) > .Machine$double.eps^0.5) {
    number <- 1 / number
    step <- step + 1
    terms[step] <- floor(number)
    # reduce reciprocal by integer part
    number <- number - terms[step]
  }
  
  return(terms)
}

# Given a simple continued fraction, with its terms in list notation, generates
# a complete table of numerators/denominators for the convergents of the CF
convergents <- function(terms) {
  # simple continued fraction should have at least one term
  if(length(terms) == 0) {
    return()
  }
  
  # single term is just the term over a denominator of 1
  if(length(terms) == 1) {
    convergents <- rbind(terms, terms, c(1))
    rownames(convergents) <- c("Simple CF Terms", "Numerators", "Denominators")
    return(convergents)
  }
  
  # resolve terms left to right (advantaged for computing future convergents),
  # store complete computation matrix of numerators/denominators of convergents
  numerators <- c(1, terms[1])
  denominators <- c(0, 1)
  for(index in 2:length(terms)) {
    # next numerator is the next term multiplied by previous numerator 
    # plus previous numerator before that
    numerators[index + 1] <- 
      terms[index] * numerators[index] + numerators[index - 1]
    # next denominator is the next term multiplied by previous denominator
    # plus previous denominator before that
    denominators[index + 1] <- 
      terms[index] * denominators[index] + denominators[index - 1]
  }
  
  convergents <- rbind(terms, numerators[-1], denominators[-1])
  rownames(convergents) <- c("Simple CF Terms", "Numerators", "Denominators")
  return(convergents)
}

# A given divergence angle can be considered as a ratio of the circumference
# of a circle to the arc length on the circle defined by that divergence angle,
# and that ratio in turn can be considered in terms of the extended Euclidean
# algorithm, with the numerator of the ratio (i.e. the arc length of the 
# divergence angle) being b, and the denominator of the ratio (i.e. the 
# circumference) being a, with a = qb + r, 0 <= r < b. Each step of the extended
# algorithm gives the linear combination for successively smaller remainders.
gen_ext_euclidean <- function(div_angle) {
  # ratio b/a is defined by ratio of circumference to arc length:
  # we can calculate this through the simple CF approx of one rotation/div angle
  # calculating the convergents of the scf gives us linear combinations of EEA
  coeff <- rational_approx(360/div_angle) |>
    convergents()
  # initial quotient (a) and divisor (b) given by final convergent
  len <- ncol(coeff)
  quotients <- as.vector(coeff["Numerators", len])
  divisors <- as.vector(coeff["Denominators", len])
  
  # all remainders can be expressed as linear combinations of a and b
  remainders <- as.vector(abs(coeff["Numerators",] * divisors[1] - 
                                coeff["Denominators",] * quotients[1]))
  
  if(len < 2) {
    full_euclidean <- rbind(coeff, quotients, divisors, remainders)
    return(full_euclidean)
  }
  
  # new quotient given by previous divisor, new divisor by previous remainder
  for(index in 2:len) {
    quotients[index] <- divisors[index-1]
    divisors[index] <- remainders[index-1]
  }
  
  full_euclidean <- rbind(coeff, quotients, divisors, remainders)
  return(full_euclidean)
}

# using pythagorean and angle sum identities, we see dist is given by cosine law
# i.e. given triangle ABC and angle c (opposite to side C):
# C = (A^2 + B^2 - 2AB cos(c))^(1/2)
polar_dist <- function(r_1, theta_1, r_2, theta_2) {
  # note that cosine is an even function (symmetric about y-axis)
  dist <- (r_1^2 + r_2^2 - 2*r_1*r_2 * cos(theta_1 - theta_2))^(1/2)
  
  return(dist)
}

#-------------------------------------------------------------------------------
# APPLICATION SPECIFIC FUNCTIONS
#-------------------------------------------------------------------------------

# Given the matrix encoding the full steps of the extended Euclidean algorithm,
# gives string descriptions of each step
gen_euclidean_strs <- function(euc) {
  rownames(euc) <- c("m", "b_m", "b_n", "q", "d", "r")
  a <- euc["q",1]
  b <- euc["d",1]
  
  # convenience for formatting coefficients of linear combinations
  format_c <- function(coeff, index) {
    if(index %% 2 == 0) {
      return(coeff)
    }else {
      return(paste("(",-1*coeff,")",sep=""))
    }
  }
  
  # each step of the algorithm is euclidean division: q = md + r
  algo_strs <- apply(euc, 2,
                     {\(x) paste(x["q"], "=", x["m"], "∙", 
                                 x["d"], "+", x["r"])})
  # each remainder r_n as a linear combination of initial quotient and divisor
  linear_comb_strs <- sapply(1:ncol(euc), 
                             {\(x) paste(euc["r", x], "=", 
                                         format_c(euc["b_n", x],(x + 1)), "∙", a, 
                                         "+ <b>", format_c(euc["b_m", x],x), "</b>∙", b)})
  strs <- rbind(algo_strs, linear_comb_strs)
  return(strs)
}

# Given a divergence angle and a count of nodes, highlights which remainders
# are visible at this given count
format_strs <- function(euclidean, num_nodes) {
  e_strs <- gen_euclidean_strs(euclidean)
  combined <- sapply(1:ncol(euclidean), 
                     {\(x) 
                       if(euclidean["Numerators", x] < num_nodes) {
                         paste("<h5 style=\"color:green\">",e_strs["algo_strs", x],
                               ", remainder as a linear combination: ",
                               e_strs["linear_comb_strs", x], "</h5>", sep="")
                       }else {
                         paste("<h5>", e_strs["algo_strs", x],
                               ", remainder as a linear combination: ",
                               e_strs["linear_comb_strs", x], "</h5>", sep="")
                       }})
  
  return(combined)
}

# arithmetic/archimedean spiral: nodes grow radially at fixed speed, distances
# between turnings are constant, r = a + b*theta
arith_nodes <- function(num_nodes, radial_dist, log_fact, meristem_radius,
                        interval = 1) {
  # oldest nodes will be furthest out, calculate from # nodes + radial dist
  n <- seq(from = (num_nodes - 1), 
           to = 0,
           by = -interval)
  dists <- meristem_radius + radial_dist * n
  return(dists)
}


# logarithmic spiral: nodes grow radially at exponential rate, distances between
# turnings perserve the spiral's shape, r = a + b*e^(c*theta)
log_nodes <- function(num_nodes, radial_dist, log_fact, meristem_radius,
                      interval = 1) {
  # oldest nodes are furthest out
  n <- seq(from = (num_nodes - 1), 
           to = 0,
           by = -interval)
  dists <- meristem_radius + radial_dist * exp(log_fact*n)
  return(dists)
}

# fermat spiral: areas between turnings is kept constant, r = a + b*theta^(1/2)
fermat_nodes <- function(num_nodes, radial_dist, log_fact, meristem_radius, 
                         interval = 1) {
  # oldest nodes are furthest out
  n <- seq(from = (num_nodes - 1), 
           to = 0,
           by = -interval)
  dists <- meristem_radius + radial_dist * n^(1/2)
  return(dists)
}

# function to generate cartesian + polar co-ordinates of nodes based on number 
# of nodes, divergence angle, exponent co-eff, and radial growth rate.
gen_nodes <- function(gen_func, num_nodes, div_angle, radial_dist, 
                      log_fact, meristem_radius, interval = 1) {
  
  # to describe each node in polar co-ordinates, need angles and lengths:
  
  # each new node's angle differs from the previous one by divergence angle
  n <- seq(from = 0, 
           to = (num_nodes - 1),
           by = interval)
  angles <- n * div_angle
  
  # use given spiral generation function to generate polar lengths of nodes
  dists <- gen_func(num_nodes, radial_dist, log_fact, meristem_radius)
  
  # convert from polar to cartesian: x = r cos theta, y = r sin theta
  x_coords <- dists * cos(deg_to_rad(angles))
  y_coords <- dists * sin(deg_to_rad(angles))
  
  # 4 x num_nodes 2d matrix will encode everything (x, y, len, angle)
  nodes <- rbind(x_coords, y_coords, dists, angles)
  return(nodes)
}

# separate computation of labels and label co-ordinates for faster loading
# when labels are not desired
gen_labels <- function(num_nodes, nodes) {
  # labels themselves are simply the order in which a given node was created
  label_text <- 1:num_nodes
  
  # labels should extend in same direction as the node, but a bit further out
  x_coords <- nodes["x_coords", ] * 1.1
  y_coords <- nodes["y_coords", ] * 1.1
  
  # size [3 x num_nodes] 2d matrix will encode everything
  labels <- rbind(label_text, x_coords, y_coords)
  return(labels)
}

# INITIALIZATION VALUES
init_num <- 16
init_denom <- 45
init_angle <- (360*init_num)/init_denom

#-------------------------------------------------------------------------------
# SHINY CODE
#-------------------------------------------------------------------------------

# define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # app title ----
  titlePanel("Plant Spirals"),
  
  # sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # sidebar panel for inputs ----
    sidebarPanel(
      
      # input: numberbox for number of primordia
      numericInput(inputId = "num_nodes", 
                   label = "Number of nodes:",
                   min = 1,
                   max = 1000,
                   value = 1),  
      
      # NOTE: we allow the divergence angle to either be defined directly,
      # or as a ratio of the angle to a whole rotation
      
      # input: numerator of ratio defining divergence angle
      div(style="display:inline-block", numericInput(inputId = "div_numerator", 
                                                     label = "Numerator", 
                                                     value = init_num)),
      
      # input: denominator of ratio defining divergence angle
      div(style="display:inline-block", numericInput(inputId = "div_denominator", 
                                                     label = "Denominator", 
                                                     min = 1,
                                                     value = init_denom)),
      
      # input: button to update divergence angle based on provided ratio
      actionButton("update_from_ratio", "Update Angle From Ratio"),
      
      # input: input box for number defining the divergence angle
      numericInput(inputId = "divergence", 
                   label = "Divergence angle:",
                   min = 1,
                   max = 360,
                   value = init_angle),  
      
      # input: button to update divergence angle by directly supplying it
      actionButton("update_from_angle", "Update Ratio From Angle"),
      
      # input: input for radial distance, allows animating nodes expanding 
      sliderInput(inputId = "radial", 
                  label = "Radial distance:",
                  min = 0,
                  max = 0.5,
                  value = 0,
                  step = 0.001,
                  round = -3,
                  animate = animationOptions(interval = 100, loop = FALSE)),  
      
      # input: dropdown menu to choose which kind of spiral to generate
      selectInput(inputId = "spiral_type",
                  label = "Type of Spiral:",
                  choices = c("Arithmetic" = "arith",
                              "Logarithmic" = "log",
                              "Fermat" = "fermat"),
                  selected = "arith"),
      
      # input: input for geometric growth factor
      sliderInput(inputId = "log_fact", 
                  label = "Growth factor:",
                  min = 0,
                  max = 0.1,
                  value = 0.01,
                  step = 0.01,
                  round = -2,
                  animate = animationOptions(interval = 100, loop = FALSE)),  
      
      # input: checkbox as to whether to label nodes with their creation order
      checkboxInput(inputId = "labels", label = "Label Nodes", value = TRUE),
      
      # input: checkbox to display continuous spiral modelling node creation
      checkboxInput(inputId = "show_spirals", 
                    label = "Show Spiral", value = FALSE),
      
      # thinner side panel
      width = 3
    ),
    
    # Main panel for displaying output ----
    mainPanel(
      tabsetPanel(
        type = "tabs",
        
        # first panel: plot of meristem and nodes ----
        tabPanel("Plot", plotOutput(outputId = "draw_nodes"),
                 htmlOutput(outputId = "euclidean_algo")),
        
        # second panel: display of distance function and distances between nodes
        tabPanel("Node Distances", 
                 # graph of distance from chosen node as function of node index
                 plotOutput(outputId = "plot_distances"),
                 
                 # input: dropdown menu to choose which node to measure from
                 selectInput(inputId = "select_node",
                             label = "Node to Measure Distance From:",
                             choices = NULL), # reactively updates
                 
                 # general distance function for spiral type, and specific for
                 # these given parameters
                 htmlOutput(outputId = "distance_function"),
                 # table of all the distances from selected node and other nodes
                 tableOutput(outputId = "distance_table"))
        
        # debugging
        #tabPanel("Debug", textOutput(outputId = "debug")),
      )
    )
  )
)

# Define server logic required to plot meristem + nodes and display other info
server <- function(input, output, session) {
  # radius of meristem, hardcoded for now
  m_radius <- 1
  
  # polar + cartesian co-ordinates of nodes (auto-update on input changes)
  node_coords <- reactive({
    # how the nodes are generated depends on the type of growth model chosen
    gen_func <- switch(input$spiral_type,
                       "arith" = arith_nodes,
                       "log" = log_nodes,
                       "fermat" = fermat_nodes)
    gen_nodes(gen_func, n_nodes(), divergence(), 
              input$radial, input$log_fact, m_radius)})
  
  # only want the divergence angle to update on button presses (manual update)
  divergence <- reactiveVal(init_angle)
  
  # ensure that the number of nodes is always valid (1 <= num <= 1000)
  n_nodes <- reactive( {
    n <- input$num_nodes
    # if values are completely invalid, reset to 1
    if(is.na(n) || !is.numeric(n) || n < 1) {
      return(1)
    }
    
    # round it if it isn't an int already
    if(!is.integer(n)) {
      n <- round(n)
    }
    
    # is the number too large?
    if(n > 1000) {
      return(1000)
    }
    
    return(n)
  })
  
  # also update input box to valid number of nodes (see above)
  observe({
    updateNumericInput(inputId = "num_nodes", value = n_nodes())
  })
  
  # update divergence angle based on ratio/direct value when buttons pressed
  # TODO add input validation for this as well
  
  # update from ratio of angle to entire rotation
  observeEvent(input$update_from_ratio, {
    div_angle <- (360 * input$div_numerator)/input$div_denominator
    divergence(div_angle)
    # also update the value of angle input field
    updateNumericInput(inputId = "divergence", 
                       value = div_angle)
    
  })
  
  # update from direct value
  observeEvent(input$update_from_angle, {
    divergence(input$divergence)
    # also update numerator/denominator input fields
    # note: last value is denom of final convergent, 2nd last is numerator
    c <- rational_approx(input$divergence/360) |> convergents()
    updateNumericInput(inputId = "div_numerator",
                       value = c[length(c)-1])
    updateNumericInput(inputId = "div_denominator",
                       value = c[length(c)])
  })
  
  # Reference node selection options updates with number of nodes
  observeEvent(n_nodes(), {
    updateSelectInput(inputId = "select_node", choices = 1:n_nodes())
  })
  
  # FIRST TAB OUTPUT
  
  # plot of meristem and nodes
  output$draw_nodes <- renderPlot({
    # node co-ordinates should be kept up to date automatically
    nodes <- node_coords()
    # scale of plot depends on radial distance of oldest node
    max_dist <- nodes[3]
    axes_scale <- c(-max_dist * 1.1, max_dist*1.1)
    
    # plot on cartesian grid, with meristem as circle centred at origin
    plot(nodes["x_coords", ], nodes["y_coords", ],
         xlim = axes_scale, ylim = axes_scale,
         axes = FALSE, ylab = '', xlab = '',
         asp = 1);
    draw.circle(0,0,m_radius)
    
    # draw node labels if desired (scale label size based on label count)
    if(input$labels) {
      labels <- gen_labels(n_nodes(), nodes)
      text(labels["x_coords",], labels["y_coords",],
           labels["label_text",], 
           cex = 3/log(n_nodes()))
    }
    
    # draw generative spiral if desired
    if(input$show_spirals) {
      # construct interpolated points at interval of 100 between nodes
      gen_func <- switch(input$spiral_type,
                         "arith" = 
                           {\(n,d,f,r) arith_nodes(n,d,f,r,interval=0.01)},
                         "log" = 
                           {\(n,d,f,r) log_nodes(n,d,f,r,interval=0.01)},
                         "fermat" = 
                           {\(n,d,f,r) fermat_nodes(n,d,f,r,interval=0.01)})
      interpolated <- gen_nodes(gen_func, 
                                n_nodes(), divergence(), 
                                input$radial, input$log_fact, m_radius,
                                interval = 0.01)
      lines(interpolated["x_coords",], interpolated["y_coords",])
    }
  }) 
  
  # text output relating extended euclidean to nodes, updates on div + nodes
  output$euclidean_algo <- renderUI({
    strs <- gen_ext_euclidean(divergence()) |>
      format_strs(n_nodes())
    
    HTML(strs)
  })
  
  # SECOND TAB OUTPUT
  
  # plot of distance function from other nodes to reference node
  output$plot_distances <- renderPlot({
    #TODO plot with function instead?
    
    # construct interpolated points at intervals of 100 between nodes
    gen_func <- switch(input$spiral_type,
                       "arith" = 
                         {\(n,d,f,r) arith_nodes(n,d,f,r,interval=0.01)},
                       "log" = 
                         {\(n,d,f,r) log_nodes(n,d,f,r,interval=0.01)},
                       "fermat" = 
                         {\(n,d,f,r) fermat_nodes(n,d,f,r,interval=0.01)})
    interpolated <- gen_nodes(gen_func, 
                              n_nodes(), divergence(), 
                              input$radial, input$log_fact, m_radius,
                              interval = 0.01)
    
    # retrieve chosen reference node
    nodes <- node_coords()
    ref_n <- nodes[,as.integer(input$select_node)]
    
    # calculate distances between chosen node and interpolated points
    i_dists <- apply(interpolated, 2, 
                       {\(x) polar_dist(ref_n[3], ref_n[4] |> deg_to_rad(), 
                                        x[3], x[4] |> deg_to_rad())})
    indices <- seq(from = 1, 
                   to = n_nodes(),
                   by = 0.01)
    
    # distances between nodes and reference node
    n_dists <- apply(nodes, 2, 
    {\(x) polar_dist(ref_n[3], ref_n[4] |> deg_to_rad(), 
                     x[3], x[4] |> deg_to_rad())})
    help(plot)
    plot(x = 1:n_nodes(), y = n_dists,
         ylim = c(0,max(i_dists)),
         xlab = "Node", ylab = "Distance from Reference Node")
    lines(indices, i_dists)
  })
  
  # TODO fix later :(
  # output$distance_function <- renderUI({
  #   # general variable declarations
  #   vars = paste("Note: n<sub>1</sub> = index of ref node, ",
  #                "n<sub>2</sub> = index of compared node, ",
  #                "r = radial growth distance, ",
  #                "&#952; = divergence angle, ",
  #                "m = meristem radius")
  #   
  #   # reference node length and angle
  #   ref_n <- node_coords()[,as.integer(input$select_node)]
  #   ref_l <- ref_n[3]
  #   ref_a <- ref_n[4] %% 360
  #   
  #   # distances depend on type of growth pattern
  #   if(input$spiral_type == "arith") {
  #     #general function
  #     general_f <- paste("(for arithmetic growth) general distance = &radic;",
  #                        "[((n<sub>1</sub>-1)r + m)<sup>2</sup> +", 
  #                        "((n<sub>2</sub>-1)r + m)<sup>2</sup> + ",
  #                        "2((n<sub>1</sub>-1)r + m)((n<sub>2</sub>-1)r + m)",
  #                        "cos((n<sub>1</sub>-n<sub>2</sub>)&#952;)]", sep = "")
  #     
  #     # we have the length and distance of the reference node, can simplify
  #     specific_f <- paste("(for arithmetic growth) general distance = &radic;",
  #                        "[((n<sub>1</sub>-1)r + m)<sup>2</sup> +", 
  #                        "((n<sub>2</sub>-1)r + m)<sup>2</sup> + ",
  #                        "2((n<sub>1</sub>-1)r + m)((n<sub>2</sub>-1)r + m)",
  #                        "cos((n<sub>1</sub>-n<sub>2</sub>)&#952;)]", sep = "")
  #     
  #     combined <- paste(vars, "<br>", general_f, sep = "")
  #     HTML(combined)
  #   }else if (input$spiral_type == "log") {
  #     #general function
  #     general_f <- paste("&radic;[((n<sub>1</sub>-1)r + m)<sup>2</sup> +", 
  #                        "((n<sub>2</sub>-1)r + m)<sup>2</sup> + ",
  #                        "2((n<sub>1</sub>-1)r + m)((n<sub>2</sub>-1)r + m)",
  #                        "cos((n<sub>1</sub>-n<sub>2</sub>)&#952;)]", sep = "")
  #     combined <- paste(vars, "<br>", general_f, sep = "")
  #     HTML(combined)
  #   }else {
  #     #general function
  #     general_f <- paste("&radic;[((n<sub>1</sub>-1)r + m)<sup>2</sup> +", 
  #                        "((n<sub>2</sub>-1)r + m)<sup>2</sup> + ",
  #                        "2((n<sub>1</sub>-1)r + m)((n<sub>2</sub>-1)r + m)",
  #                        "cos((n<sub>1</sub>-n<sub>2</sub>)&#952;)]", sep = "")
  #     combined <- paste(vars, "<br>", general_f, sep = "")
  #     HTML(combined)
  #   }
  # })
  
  # TODO allow sorting options, better scrolling
  output$distance_table <- renderTable({
    # retrieve chosen reference node
    nodes <- node_coords()
    ref_n <- nodes[,as.integer(input$select_node)]
    
    # calculate distances between chosen node and all other nodes
    distances <- apply(nodes, 2, 
                       {\(x) polar_dist(ref_n[3], ref_n[4] |> deg_to_rad(), 
                                        x[3], x[4] |> deg_to_rad())})
    table_data <- cbind(1:n_nodes(), distances)
    colnames(table_data) <- c("Node Index", "Distance")
    return(table_data)
  })
  
  # debug text (should not show up in production!)
  output$debug <- renderText({
    nodes <- gen_nodes(input$num_nodes, divergence(), input$radial, m_radius)
    paste("Co-ordinates for", input$num_nodes, "nodes: ")
    #paste(nodes)
    #labels <- gen_labels(input$num_nodes, nodes)
    #paste(labels)
  })
  
  # TODO remove in production, for my sanity running local rstudio app 
  # gosh-forsaken backend DB connectors!
  session$onSessionEnded(stopApp)
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

#-------------------------------------------------------------------------------
# TESTING CODE
#-------------------------------------------------------------------------------


# anonymous functions and pipes
# test_arr <- 0:10
# test_arr |> sapply(\(x) return(x + 1))
# test_arr |> sapply((\(x) return(x + 1)))
# test_arr |> sapply({\(x) return(x + 1)})
# 
# #intermediate function output formats
# rational_approx(2.8125)
# rational_approx(2.8125) |> convergents()
# gen_ext_euclidean(360/2.8125)
# 
# 
# # fermat
# nums <- 1:10
# (nums-1)^(1/2)
# 
# # node polar distance generation
# seq(from = 0, to = 0, by = 128)
# n <- 20
# a <- 128
# d <- 0
# f <- 1
# m <- 1
# 
# test_nodes <- gen_nodes(arith_nodes, n, a, d, f, m)
# test_nodes[,0][3]
# apply(test_nodes, 2, {\(x) polar_dist(1,0,x[3],x[4])})
# 
# # plot on cartesian grid, with meristem as circle centred at origin
# max_dist <- test_nodes[3]
# axes_scale <- c(-max_dist * 1.1, max_dist*1.1)
# plot(test_nodes["x_coords", ], test_nodes["y_coords", ],
#      xlim = axes_scale, ylim = axes_scale,
#      axes = FALSE, ylab = '', xlab = '',
#      asp=1);
# draw.circle(0,0,1)
# 
# # spiral plotting
# test_spiral <- gen_nodes({\(n,d,f,r) arith_nodes(n,d,f,r,interval=0.1)}, 
#                          n, a, d, f, m, interval = 0.1)
# 
# lines(test_spiral["x_coords", ],test_spiral["y_coords", ])
# 
# 
# # distance between nodes
# 
# polar_dist(1,0,1.1,deg_to_rad(30))
# polar_dist(1,0,2,deg_to_rad(300))
# 
# test_dists <- seq(from=1, to=3, by=0.1)
# test_angles <- seq(from=0, to=600, by=30) |> deg_to_rad()
# test_polar_points <- rbind(test_dists, test_angles)
# test_polar_points
# apply(test_polar_points, 2, {\(x) polar_dist(1,0,x[1],x[2])})
# polar_dist(0,deg_to_rad(30),2,deg_to_rad(90))
# polar_dist(1,deg_to_rad(30),2,deg_to_rad(90))
# 
# # distance data
# n <- 20
# a <- 128
# d <- 0.1
# f <- 1
# m <- 1
# 
# test_nodes <- gen_nodes(arith_nodes, n, a, d, f, m)
# test_nodes
# test_spiral <- gen_nodes({\(n,d,f,r) arith_nodes(n,d,f,r,interval=0.1)}, 
#                          n, a, d, f, m, interval = 0.1)
# test_spiral
# apply(test_nodes, 2, {\(x) polar_dist(2.9,0,x[1],x[2])})
# dists <- apply(test_spiral, 2, {\(x) polar_dist(2.9,0,x[1],x[2])})
# max(dists)
# 
# 
# 
# # find nearest neighbours
# testing_dist_nodes <- gen_nodes(100, 137.5, 0.05, 1)
# 
# dist(testing_dist_nodes)
# 
# count <- 20
# spacing <- 0.1
# start <- 10
# 
# dists <- seq(from = start + count * spacing, 
#              to = start, 
#              by = -spacing)
# dists
# 
# log_increments <- (spacing * 0.9^(seq(0,count-2,1))) |>
#   cumsum()
# log_increments
# 
# log_dists <- (log_increments + start) |> 
#   append(start, 0) |>
#   sort(decreasing = TRUE)
# log_dists
# 
# 
# x <- cbind(c(1:5), c(10,2,6,8,11))
# x
# rev(x)

