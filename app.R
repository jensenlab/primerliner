#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

fromto <- function(from, to) if (to < from) numeric(0) else seq(from, to, by=1)

complement <- function(dnas) chartr("acgtACTG", "tgcaTGAC", dnas)
reverse <- function(dnas) sapply(lapply(strsplit(dnas, NULL), rev), paste, collapse="")
rc <- function(dnas) reverse(complement(dnas))

str2char <- function(x) strsplit(x, NULL)[[1]]
char2str <- function(x) paste(x, collapse="")
str2mat <- function(x) t(str2char(x))

match_length <- function(pattern, x) attr(regexpr(pattern, x), "match.length")
first_notna <- function(x) which(!is.na(x))[1]

align_vec <- function(x1, x2) {
  l1 <- length(x1)
  l2 <- length(x2)
  X <- matrix(NA, nrow=2, ncol=2*l1+l2-1)
  matches <- numeric(l1+l2)
  X[2,(l1+1):(l1+l2)] <- x2
  for (i in 1:(l1+l2)) {
    X[1,i:(l1+i-1)] <- x1
    matches[i] <- sum(apply(X, 2, function(c) !any(is.na(c)) && tolower(c[1]) == tolower(c[2])))
    X[1,i:(l1+i-1)] <- NA
  }
  pos <- which.max(matches)
  return(pos - (l1+1))
}

align_mat <- function(M1, M2) {
  r1 <- dim(M1)[[1]]
  c1 <- dim(M1)[[2]]
  r2 <- dim(M2)[[1]]
  c2 <- dim(M2)[[2]]
  
  offset <- align_vec(M1[r1, ], M2[1, ])
  offset1 <- ifelse(offset < 0, 0, offset)
  offset2 <- ifelse(offset < 0, abs(offset), 0)
  M <- matrix(NA, nrow=r1+r2, max(offset1+c1, offset2+c2))
  M[1:r1,(offset1+1):(offset1+c1)] <- M1
  M[(r1+1):(r1+r2),(offset2+1):(offset2+c2)] <- M2
  return(M)
}

align_chars <- function(chars) {
  Reduce(align_mat, lapply(chars[2:length(chars)], str2mat), str2mat(chars[1]))
}

align_oligos <- function(forward, reverse, remove_spaces=TRUE) {
  if (remove_spaces) {
    forward <- gsub("\\s", "", forward)
    reverse <- gsub("\\s", "", reverse)
  }
  nf <- length(forward)
  nr <- length(reverse)
  strs <- character(0)
  if (nf > 0) {
    strs <- c(strs, forward)
  } 
  if (nr > 0) {
    strs <- c(strs, sapply(reverse, rc))
  }
  
  A <- align_chars(strs)
  ar <- dim(A)[[1]]
  ac <- dim(A)[[2]]
  decs <- matrix(NA, nrow=ar-1, ncol=ac)
  for (i in fromto(1, nf-1)) {
    same <- !is.na(A[i, ]) & !is.na(A[i+1, ]) & tolower(A[i, ]) == tolower(A[i+1, ])
    decs[i,same] <- "v"
  }
  if (nr > 0) {
    same <- !is.na(A[nf, ]) & !is.na(A[nf+1, ]) & tolower(A[nf, ]) == tolower(A[nf+1, ])
    decs[nf,same] <- "|"
  }
  for (i in fromto(1, nr-1)) {
    same <- !is.na(A[nf+i, ]) & !is.na(A[nf+i+1, ]) & tolower(A[nf+i, ]) == tolower(A[nf+i+1, ])
    decs[nf+i,same] <- "^"
  }
  
  for (i in fromto(1, nr)) {
    A[nf+i, ] <- sapply(A[nf+i, ], complement)
  }
  
  final <- matrix(NA, nrow=2*ar-1, ncol=ac)
  for (i in fromto(1, ar-1)) {
    final[2*i-1, ] <- A[i, ]
    final[2*i, ] <- decs[i, ]
  }
  final[2*ar-1, ] <- A[ar, ]
  final[is.na(final)] <- " "
  
  return(apply(final, 1, paste, collapse=""))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Primerliner: the simple oligo aligner"),
   
   textInput("f1", "Forward Oligos"),
   textInput("f2", NULL),
   textInput("f3", NULL),
   textInput("r1", "Reverse Oligos"),
   textInput("r2", NULL),
   textInput("r3", NULL),
   verbatimTextOutput("value")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$value <- renderText({
     forward <- c(input$f1, input$f2, input$f3)
     reverse <- c(input$r1, input$r2, input$r3)
     forward <- forward[nzchar(forward)]
     reverse <- reverse[nzchar(reverse)]
     if (length(forward) + length(reverse) > 1) {
       paste(align_oligos(forward, reverse), collapse="\n")
     } else {
       paste(c(forward, reverse), collapse="\n")
     }
      
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

