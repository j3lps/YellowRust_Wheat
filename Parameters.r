# This contains all the parameters needed for model simulation

Parameters = list()

# ********* Crop parameters ***********#

Parameters$Crop = c(
	# Growth rate of wheat canopy
	growthRate = 0.06456,
	# Senescence coefficient of wheat
	senesRate = 0.05,
	# Maximum leaf area index
	maxA = 7.565,
	# The time (in Julian days after planting) at which rate of senescence = 1
	tSen = 347,
	# Initial leaf area index at emergence
	LAI0 = 1.395e-5,
	# Time at which harvest occurs
	tHarvest = 333
)

# ********* Pathogen-weather parameters ***********#
Parameters$PathWeath = list()
Parameters$PathWeath$LatentPeriod = c(a=28.24 ,b=0.06)
Parameters$PathWeath$InfectiousPeriod = c(m=0.06,c=-0.02)
Parameters$PathWeath$InfectionEfficiency = c(T_max=19.8,T_min=2.37,p=2.24,n=0.87,m=0.41,c=0.045,d=-0.065,I_max=0.421,f=-0.023,g=0.0246,h=0.00101,j=10.14,k=1.024,r=0.0427)
Parameters$PathWeath$SporulationRate = c(a=8799,T_min=1.7,T_max=20.84,n=0.84,m=0.58)

# ********* Pathogen parameters***********#
Parameters$Pathogen = c(
	# Transmission coefficient
	alpha = 0.003,
	# Number of latent stages
	nLS = 4,
	# Number of infectious stages - determines the shape of the sporulation curve, see paper
	nIS = 4,
	# Time of the inoculuation after planting
	tInoc = 180,
	# Severity at inoculation (between 0 and 1)
	inocSev = 0.1
)
