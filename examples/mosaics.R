library(vcdExtra)
# two-way mosaic displays
HairEye <- margin.table(HairEyeColor, c(1,2))
HairEye <- HairEye[,c(1, 3, 4, 2)]
HairEye

vcd::mosaic(HairEye, shade=TRUE)

# How to fill the mosaic with points?

library(ggmosaic)
library(forcats)

ggplot(data = titanic) +
  geom_mosaic(aes(x = product(Class), fill = Survived), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Class), color = Survived))

str(titanic)

haireye <- expand.dft(HairEye)
str(haireye)

ggplot(data = haireye) +
  geom_mosaic(aes(x = product(Hair), fill = Eye), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Hair), color = Eye)) +
  coord_flip()

# print freq in cells -- doesn't work
ggplot(data = haireye) +
  geom_mosaic(aes(x = product(Hair), fill = Eye), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Hair), color = Eye)) +
  geom_mosaic_text(aes(x = product(Eye, Hair),
                       label = after_stat(.wt)),
                   size = 5) +
  coord_flip()

geom_mosaic_text(aes(x = product(do_you_recline, rude_to_recline),
                     label = after_stat(.wt)))
# expected frequencies?
ggplot(data = haireye) +
  geom_mosaic(aes(x = product(Hair), fill = Eye, weight = 1), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Hair), color = Eye)) +
  coord_flip()


# How to keep order the same?

# doesn't work

# doesn't work
ggplot(data = haireye) +
  geom_mosaic(aes(x = product(fct_inorder(Hair)), fill = Eye), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(fct_inorder(Hair)), color = Eye)) +
  coord_flip()

haireye2 <- haireye |>
  mutate(Hair = fct_inorder(as.factor(Hair)),
         Eye = fct_inorder(as.factor(Eye)))


ggplot(data = haireye2) +
  geom_mosaic(aes(x = product(Hair), fill = Eye), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Hair), color = Eye)) +
  coord_flip()


# Warning messages:
#   1: Computation failed in `stat_mosaic()`.
# Caused by error in `[.data.frame`:
#   ! undefined columns selected

# try factors

# no effect
haireye2 <- haireye |>
  mutate(Hair = as.factor(Hair),
         Eye = as.factor(Eye))
ggplot(data = haireye2) +
  geom_mosaic(aes(x = product(Hair), fill = Eye), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Hair), color = Eye)) +
  coord_flip()

haireye3 <- haireye |>
  mutate(Hair = ordered(Hair, levels = c("Black", "Brown", "Red", "Blond")),
         Eye = ordered(Eye, levels = c("Brown", "Hazel", "Green", "Blue"))
  )
ggplot(data = haireye2) +
  geom_mosaic(aes(x = product(Hair), fill = Eye), alpha = 0.3) +
  geom_mosaic_jitter(aes(x = product(Hair), color = Eye)) +
  coord_flip()




