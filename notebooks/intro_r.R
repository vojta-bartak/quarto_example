# Generate 100 random numbers from normal distribution
x <- rnorm(100)
summary(x)

# Plot histogram
hist(x, main = "Histogram of x", col = "skyblue", border = "white")
