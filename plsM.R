library(mdatools)
data(people)

idx = seq(4, 32, 4)
Xc = people[-idx, -4]
yc = people[-idx, 4, drop = FALSE]
Xt = people[idx, -4]
yt = people[idx, 4, drop = FALSE]

m = pls(Xc, yc, 7, scale = TRUE, info = "Shoesize prediction model")

# Gives an error, needs prper validation
# One way to do validation is to do the "full cross-validation"
m = pls(Xc, yc, 7, scale = TRUE, cv = 1, info = "Shoesize prediction model")
m = selectCompNum(m, 3)

summary(m)
plot(m)

par(mfrow = c(2, 2))
plotRegcoeffs(m)
plotRegcoeffs(m, ncomp = 2)
plot(m$coeffs, ncomp = 3, type = "b", show.labels = TRUE)
plot(m$coeffs, ncomp = 2)

