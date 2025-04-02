# normal equations in matrix form

library(matlib)

XPX <- matrix(c(
  "n", "\\sum x",
  "\\sum x", "\\sum x^2"), 2, 2)

XPY <- matrix(c("\\sum y", "\\sum x y"), 2,1)

b <- matrix(c("a", "b"), 2, 1)

# Eqn(latexMatrix(XPX, matrix = "bmatrix"),
#     latexMatrix(b),
#     "=",
#     latexMatrix(XPY))
#
# Eqn(
#   underbrace(latexMatrix(XPX, matrix = "bmatrix"), label ="\\mathbf{X}^\\top\\mathbf{X}"),
#   underbrace(latexMatrix(b), label = "\\mathbf{b}") ,
#     "=",
#   underbrace(latexMatrix(XPY), label = "\\mathbf{X}^\\top\\mathbf{y}")
# )


Eqn(latexMatrix(XPX, matrix = "bmatrix"),
    latexMatrix(b),
    "& =",
    latexMatrix(XPY),
    Eqn_newline(),
    "\\mathbf{X}^\\top \\mathbf{X}",
    ' & =',
    "\\mathbf{X}^\\top \\mathbf{y}",
    align = TRUE
)



# Letting $\\mathbf{X} = [ \mathbf{1} \\mathbf{x}]$ and $\mathbf{b} = ( a / b)$

latexMatrix("\\mathbf{1} \\mathbf{x}", 1, 1, matrix = "bmatrix")

latexMatrix(matrix(c("a", "b")), 2,1)

latexMatrix(matrix(c("\\mathbf{1}", "\\mathbf{x}"), 1, 2), matrix = "bmatrix")

