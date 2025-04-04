# quadratic for potential energy

k <- 1/2
x <- seq(-6, 6, 0.5)
Px <- k * x^2

linepts <- function(a, b, x) {
  y <- a + b * x
  data.frame(x,y)
}

# # line through (4, 16), extending from x=2 to x=6
# #   slope = d k X^2 / dx = 2 k X
# #   -> y - 16 = 2*k*4 * (x - 4) -> y = 8 X -16
# pts <- linepts(0, 8, c(2,4,6))


op <- par(mar = c(4,4,2,1)+.5)
plot(x, Px,
     type = "l",
     lwd = 2,
     axes = FALSE,
     ylab = "",
     cex.lab = 2)
axis(1)
axis(2)
abline(v=0, lwd = 1.2)
text(0, 18, "P(x)", cex = 2, xpd = TRUE)

# line through (3, 4.5), extending from x=2 to x=4
#   slope = d k X^2 / dx = X = 3
#   -> y - 4.5 = 3 * (x - 3) -> y = 3 X - 9
pts <- linepts(-4.5, 3, c(1.5, 3, 4.5))

lines(pts, col = "blue", lwd = 3)
points(pts[2,], pch = 16, cex= 2.3, col = "red")

text(-3, 14, expression(P(x) == frac(1,2) * X^2), cex = 1.2)
text(-2.8, 12, expression(F(x) == 2 * X), cex = 1.2)
#text(-2.8, 12, expression(paste(F(x) == -P), "'", sep=""), cex = 1.2)

#plot(1,1,xlab=expression(paste(K,"'",(t),sep="")))

points(0, 0, pch = 16, cex= 1.5, col = "red")
segments(x0 = -1.2, x1 = 1.2,
         y0 = 0, y1 = 0, col = "blue", lwd=2)


text(4 + c(0, .5), 6 + c(-0.7, -3),
     c(expression(F(x) == -dP/dx),
       expression(F(x) == - slope[x])),
     cex = 1.4,
     srt = 60,
     col = "blue",
     xpd = TRUE)
par(op)

