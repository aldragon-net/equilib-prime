class Species:
    """class of mixture components"""
    def __init__(self, name, molefrac, weight, low_T, threeshold_T,  high_T, a_low, a_high):
        self.name = name
        self.molefrac = molefrac
        self.weight = weight
        self.low_T = low_T
        self.threeshold_T = threeshold_T
        self.high_T = high_T
        self.a_low = a_low
        self.a_high = a_high

		
def read_thermodat(composition, file)
    mixture = ()
    thermdat = open('therm.dat', 'r')
    for i in len(composition):
	    while True:
            line = readline(thermdat)					#1st line
            if line[:18].strip = composition[i]:
               break			
	    name = line[:18].strip
		low_T = line[45:55]
		threeshold_T = line[55:65]
		high_T = line[65:73]
		line = readline(thermdat)						#2nd line
		a_high[0] = line[0:15]
		a_high[1] = line[15:30]
		a_high[2] = line[30:45]
		a_high[3] = line[45:60]
		a_high[4] = line[60:75]
		line = readline(thermdat)						#3rd line
		a_high[5] = line[0:15]
		a_high[6] = line[15:30]
		a_high[7] = line[30:45]
		a_low[1] = line[45:60]
		a_low[2] = line[60:75]
		line = readline(thermdat)						#4th line
		a_low[3] = line[0:15]
		a_low[4] = line[15:30]
		a_low[5] = line[30:45]
		a_low[6] = line[45:60]
		a_low[7] = line[60:75]
		
	return mixture

#x = Species('Ar', 1, 40, 200, 2000, 6000, (0,1,2,3,4,5), (6,7,8,9,10))
#print(x.a_low)
#print(x.a_high[0])
