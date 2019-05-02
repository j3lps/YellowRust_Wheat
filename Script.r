# Need an odeSolver, here we use the deSolve package, if it isn't installed use 'install.packages("deSolve") to install it'
library(deSolve)

# Set the working directory to the location of your repository
#setwd("C:/Users/...")

# Read in the model equations
source("./ModelEquations.r")
# Read in the parameters
source("./Parameters.r")
# Read in the weather data - contains ECN hourly weather data
weatherData = read.csv("./WData.csv")

# Need to adjust the weather data from hourly weather into daily weather, and format appropriately:

# Convert date to date class
weatherData$SDATE = as.Date(weatherData$SDATE,format="%Y-%m-%d")
weatherData$Year = as.integer(format(weatherData$SDATE,"%Y"))

# Cut the weather data between the start and harvest
startDay = which(weatherData$SDATE==paste(2001,"-10-01",sep="") & weatherData$SHOUR=="1")
endDay = which(weatherData$SDATE==paste(2002,"-08-31",sep="") & weatherData$SHOUR=="24")
if(length(startDay) == 0 | length(endDay) == 0){
	print(paste("Year",2002,"cannot find a start or end date"))
	next
}
weatherData = weatherData[startDay:endDay,]

# Convert from hourly weather to daily weather
weatherDataDaily = aggregate(SHOUR ~ SDATE, weatherData, FUN = function(x) length(unique(x)))
weatherDataDaily = weatherDataDaily[order(as.Date(weatherDataDaily$SDATE, format="%Y-%m-%d")),]
DailyData = aggregate(weatherData$DRYTMP,by=list(weatherData$SDATE),mean)
SurwetDaily = aggregate(weatherData$SURWET,by=list(weatherData$SDATE),sum)
DailyData$SURWET = SurwetDaily$x/60
DailyData$Radiation = aggregate(weatherData$SOLAR,by=list(weatherData$SDATE),mean)$x * 2.3 / 10
names(DailyData) = c("SDATE","Tmean","SurWet","Light")

# Plot the temperature to check everything works
plot(DailyData$Tmean~DailyData$SDATE,type="l",xlab="Time (days)",ylab="Temperature")

# Now we can run the model using the parameters, model and weather data

# Run the model for a single year
output = runModel(Parameters,DailyData)

# Plot the severity:
plot(output$severity~output$time,type="l",xlab="Time (Julian Days after planting)",ylab="Severity (%)",ylim=c(0,100))
