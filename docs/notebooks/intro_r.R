# Create a sequence of values from 0 to 10 with steps of 0.1
x <- seq(0, 10, by = 0.1)

# Calculate the sine of x
y <- sin(x)

# Plot the sine wave
plot(x, y, type="l", main="Sine Wave")
