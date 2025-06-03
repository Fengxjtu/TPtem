rm(list = ls())
library(stats)
library(forecast)
library(pracma)
library(ggplot2)
library(dplyr)
library(zoo)
library(Rssa)


pre <- read.csv("D:/SSAyece/TPtem.csv")
colnames(pre) <- c("Year", "TEM")
data <- pre

ts_data <- ts(data$PRE, start = min(data$Year), end = max(data$Year), frequency = 1)
#  Calculate sliding average (aligned to the right)
ma21 <- rollmean(ts_data, k = 1, align = "right")
ma_years <- (min(data$Year) + 0):(max(data$Year))
ma_ts <- ts(ma21, start = min(ma_years), frequency = 1)

filtered_ma_ts <- window(ma_ts, start = 1352, end = 2010)
train_data <- window(filtered_ma_ts, start = min(time(filtered_ma_ts)), end = max(time(filtered_ma_ts)))


L <- floor(length(train_data)/2)
ssa_train <- ssa(train_data, L = L, kind = "1d-ssa",neig = L)
plot(ssa_train, type = "values", main = "Eigenvalue map")
plot(ssa_train, type = "vectors", main = "Principal Component Diagram")

wcor <- wcor(ssa_train)
plot(wcor, main = "W-Correlation Matrix")

eigenvalues <- ssa_train$sigma^2
total_variance <- sum(eigenvalues)
variance_explained <- eigenvalues / total_variance * 100
cumulative_variance <- cumsum(variance_explained)
components <- 1:length(eigenvalues)
variance_df <- data.frame(
  Component = components,
  Variance = variance_explained,
  Cumulative = cumulative_variance
)

top_10 <- data.frame(
  Component = 1:20,
  Variance_Explained = round(variance_explained[1:20], 2),
  Cumulative_Variance = round(cumulative_variance[1:20], 2)
)
write.csv(top_10,"PCA_variance.csv",row.names = F)


#Set the number and grouping of groups
n_components <- 100
groups <- list(
  trend = 1:4,            
  seasonal = 5:8,        
  noise = 9:n_components 
)


recon <- reconstruct(ssa_train, groups = groups)

trend <- recon$trend
seasonal <- recon$seasonal
noise <- recon$noise
reconstructed <- trend + seasonal

plot(train_data, main = "Comparison between raw and reconstructed", col = "black", lwd = 2)
lines(reconstructed, col = "red", lwd = 2, lty = 2)
lines(seasonal,col = "green", lwd = 2, lty = 2)
lines(trend,col = "purple", lwd = 2, lty = 2)
#legend("topright", legend = c("RAW", "REC","seasonal","trend"),
#col = c("black", "red","green","purple"), lty = c(1, 2), lwd = 2)

train_results <- data.frame(
  Year = time(filtered_ma_ts),
  Trend = as.numeric(trend),
  Seasonal = as.numeric(seasonal),
  Noise = as.numeric(noise)
)

write.csv(train_results, "train_results.csv", row.names = FALSE)

residuals <- train_data - reconstructed
rmse <- sqrt(mean(residuals^2))

residual_sd <- sd(residuals)
confidence_level <- 0.95
z_value <- qnorm(1 - (1 - confidence_level)/2) # 


forecast_years_future <- 2011:2100
forecast_len <- length(forecast_years_future)
future_all <- rforecast(ssa_train,groups = groups,
                        base = c("reconstructed"),
                        len = forecast_len,
                        only.new = TRUE)

future = future_all$trend+future_all$seasonal

future_ts <- window(ma_ts, start = 2011, end = 2100)
futre_data <- window(future_ts, start = min(time(future_ts)), end = max(time(future_ts)))
#residuals2 <- futre_data - future
#rmse2 <- sqrt(mean(residuals2^2))

#residual_sd2 <- sd(residuals2)
#lower_bound <- future - z_value * residual_sd2
#upper_bound <- future + z_value * residual_sd2

lower_bound <- future - z_value * residual_sd
upper_bound <- future + z_value * residual_sd


plot(ma_ts, 
     main = "SSA Forecast with Confidence Intervals",
     ylab = "PRE",
     xlab = "Year",ylim = c(20,22),
     xlim = c(min(ma_years), max(forecast_years_future)),
     col = "black", lwd = 2)
lines(reconstructed, col = "red", lwd = 2, lty = 2)
lines(ts(future, start = min(time(future)), frequency = 1), 
      col = "green", lwd = 2)
lines(ts(lower_bound, start = min(time(future)), frequency = 1), 
      col = "blue", lwd = 1, lty = 2)
lines(ts(upper_bound, start = min(time(future)), frequency = 1), 
      col = "blue", lwd = 1, lty = 2)
legend("topleft", 
       legend = c("Historical", "Reconstructed", "Forecast", "95% Confidence Interval"),
       col = c("black","red",  "green", "blue"),
       lty = c(1, 2, 1, 2), lwd = c(2,2,2,1))


forecast_df <- data.frame(
  Year = forecast_years_future,
  Forecast = as.numeric(future),
  Lower_95 = as.numeric(lower_bound),
  Upper_95 = as.numeric(upper_bound)
)
write.csv(forecast_df, "full_analysis_results.csv", row.names = FALSE)
########################################################################################




