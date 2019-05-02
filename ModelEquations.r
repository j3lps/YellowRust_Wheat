
# Calculate the AUDPC when passed an object output from runModel
calculateAUDPC<-function(resultDataFrame){
	
	# Check that the names severity and time exist
	if(sum(names(resultDataFrame)=="severity") != 1 | sum(names(resultDataFrame)=="time") != 1){
		stop("Severity or time column not found in results object when calculating AUDPC")
	}
	
	# Check for any na's in severity or time:
	if(sum(is.na(resultDataFrame$severity)) > 0 | sum(is.na(resultDataFrame$time))>0){
		stop("NAs found in severity or time when calculating AUDPC")
	}
	
	AUDPC = 0.0
	
	for(i in 2:nrow(resultDataFrame)){
		AUDPC = AUDPC + 0.5 * (resultDataFrame$severity[i]+resultDataFrame$severity[i-1])/(resultDataFrame$time[i]-resultDataFrame$time[i-1])
	}
	
	return(AUDPC)
	
}

# Calculate the AUDPC from the area of the infectious lesions (rather than from severity)
calculateAUDPC2<-function(resultDataFrame){
	
	AUDPC = 0.0
	
	for(i in 2:nrow(resultDataFrame)){
		AUDPC = AUDPC + 0.5 * (resultDataFrame$I[i]+resultDataFrame$I[i-1])/(resultDataFrame$time[i]-resultDataFrame$time[i-1])
	}
	
	return(AUDPC)
	
}

# Calculate the latent period from the temperature
getLatentPeriod<-function(tempNow,parms){
	
	with(as.list(parms),{
	
		currentLP = a * exp(-b * tempNow)
		
		return(currentLP)
	
	})
	
}

# Calculate the infectious efficiency  
getInfectionEfficiency<-function(tempNow,wetnessNow,lightNow,parms){
	
	with(as.list(parms),{
	
		# On temperature
		IEtemp = p * ((tempNow - T_min) / (T_max - T_min))^n * 
					((T_max - tempNow) / (T_max - T_min))^m
		IEtemp[tempNow >= T_max | tempNow <= T_min] = 0.0
		IEtemp[IEtemp < 0] = 0
		
		# On wetness duration  
		intInfRate = f + (g*tempNow) - (h * (tempNow^2))
		Wmin = j - (k * tempNow) + (r * (tempNow^2))
		IEwetness = 1 - exp(- intInfRate * (wetnessNow - Wmin))
		IEwetness[IEwetness < 0] = 0
		
		# On light quantity
		IElight = 1 - exp(-c * (lightNow + d))
		
		# Relative infectious efficiency
		currentIE = I_max * IEtemp * IEwetness * IElight
		
		return(currentIE)
	})
}

# Calculate the total spore production of a single lesion from temperature
getTotalSporulation<-function(tempNow,parms){
	
	with(as.list(parms),{
		currentSR = a * ((tempNow-T_min)/(T_max-T_min))^n * 
						((T_max-tempNow)/(T_max-T_min))^m

		# If the temperature is above the maximum or below the minimum, then set currentSR to zero
		currentSR[tempNow >= T_max | tempNow <= T_min] = 0
		# If the function has gone negative, set to zero
		currentSR[currentSR < 0] = 0

		return(currentSR)
	})
}

# Calculate the infectious period from temperature
getInfectiousPeriod<-function(tempNow,parms){
	
	with(as.list(parms),{
		# Calculate the infectious period (4 is the number of infectious stages)
		currentIP = 4 / (m * tempNow + c)
		# If the infectious period would be negative, make it very very small (if zero, then end up with 1/0)
		currentIP[currentIP < 0.0001] = 0.0001

		return(currentIP)
	})
}

# This is the overall model equations.
modelEquations<-function(t,y,Parameters,wData){
	
	# Create an environment with y and parms unwrapped
	with(as.list(c(y,Parameters)),{
		
		# Create an empty vector to store the derivatives
		dy = rep(0,length(y))
		
		nLatentStages = Pathogen[["nLS"]]
		nInfecStages = Pathogen[["nIS"]]
		
		# If we have some pathogen, then check for weather variables, otherwise we don't need to
		if(sum(y[(2+1):(2+nLatentStages+nInfecStages)])>0){
			
			# Get today's weather
			todaysWeather = wData[floor(t+1):floor(t+2),]
			difference = t-floor(t)
			temp=todaysWeather$Tmean[1] * difference + todaysWeather$Tmean[2] * (1-difference)
			surwet=todaysWeather$SurWet[1] * difference + todaysWeather$SurWet[2] * (1-difference)
			light=todaysWeather$Light[1] * difference + todaysWeather$Light[2] * (1-difference)
			if(length(temp) > 1 | length(surwet) > 1 | length(light) > 1) print(paste("Warning: Weather parameter error: Site:",site,"Year:",year,"Day:",currentDay))
			
			# Calculate the latent period, infectious period, infection efficiency and total sporulation rate from the weather
			p = getLatentPeriod(temp,parms=Parameters$PathWeath$LatentPeriod)
			TSR = getTotalSporulation(temp,parms=Parameters$PathWeath$SporulationRate)
			IE = getInfectionEfficiency(temp,surwet,light,parms=Parameters$PathWeath$InfectionEfficiency)
			i = getInfectiousPeriod(temp,parms=Parameters$PathWeath$InfectiousPeriod)
			
		} else {
			
			# Since there's no pathogen these values don't matter
			p = Inf
			TSR = 0
			IE = 0
			i = Inf
			
		}
		
		# Convert the sporulation rate and infection efficiency into a infection rate
		rB = Pathogen[["alpha"]] * IE * TSR / i
		
		# Model system
		# dA/dt
		dy[1] = Crop[["growthRate"]] * A * (1 - A / Crop[["maxA"]])
		# dH/dt
		dy[2] = Crop[["growthRate"]] * A * (1 - A / Crop[["maxA"]]) - rB * y[2+nLatentStages+nInfecStages] * nInfecStages * H - H * exp(Crop[["senesRate"]] * (t-Crop[["tSen"]]))
		# dL(j)/dt
		for(j in 1:nLatentStages){
			if(j == 1){
				dy[3] = rB * sum(y[2+nLatentStages+nInfecStages]) * nInfecStages * H
			} else { 
				dy[2+j] = y[2+j-1] * nLatentStages / p 
			}
			dy[2+j] = dy[2+j] - y[2+j] * nLatentStages / p - y[2+j] * exp(Crop[["senesRate"]] * (t-Crop[["tSen"]]))
		}
		# dI(j)/dt - Calculate the derivative of each infectious stage
		for(j in 1:nInfecStages){
			if(j==1){
				# The derivative of the first infectious stage comes from the last latent stage
				dy[2 + nLatentStages+1] = y[2 + nLatentStages] * nLatentStages / p
			} else { 
				dy[2 + nLatentStages + j] = y[2 + nLatentStages + j - 1] * nInfecStages / i 
			}
			dy[2 + nLatentStages + j] = dy[2+nLatentStages+j] - y[2+nLatentStages+j] * nInfecStages / i - y[2+nLatentStages+j] * exp(Crop[["senesRate"]] * (t-Crop[["tSen"]]))
		}
		# dR/dt
		dy[2+nLatentStages+nInfecStages+1] = y[2+nLatentStages+nInfecStages] * nInfecStages / i - exp(Crop[["senesRate"]] * (t-Crop[["tSen"]])) * y[2+nLatentStages+nInfecStages+1]
		
		# Return the derivatives
		names(dy) = names(y)
		return(list(dy))
		
	})
}

runModel<-function(Parameters,weatherData){
	
	# Run the model up to the time of inoculation
	output<-as.data.frame(ode(y=c(A=Parameters$Crop[["LAI0"]],H=Parameters$Crop[["LAI0"]],L=rep(0,Parameters$Pathogen[["nLS"]]),I=rep(0,Parameters$Pathogen[["nIS"]]),R=0.0),
							  parms=Parameters,
							  times=seq(0,Parameters$Pathogen[["tInoc"]]),
							  func=modelEquations,wData=weatherData)
							)
	
	# Inoculate with some latent tissue
	init=c(output[nrow(output),-1]);init=unlist(init);init[3]=init[2]*Parameters$Pathogen[["inocSev"]]; init[2] = init[2]*(1-Parameters$Pathogen[["inocSev"]]);
	output<-rbind(output[-nrow(output),],
				  as.data.frame(ode(y=init,
									parms=Parameters,
									times=seq(Parameters$Pathogen[["tInoc"]],Parameters$Crop[["tHarvest"]]),
									func=modelEquations,wData=weatherData)
				  )
	)
	
	# Get the latent and infectious columns
	output$L = rowSums(output[,grep("L",names(output))])
	output$I = rowSums(output[,grep("I",names(output))])
	
	# Calculate severity
	output$severity = output$I / (output$H + output$L + output$I) * 100
	
	# Return the results as a data.frame
	return(output)
}
