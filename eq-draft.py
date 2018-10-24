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
        
def readinpfile(inppath):
    """read input parameters from equilib.inp"""
    with open(inppath, 'r') as inpfile:
        while True:
                line = inpfile.readline()
                if line.strip() == '':
                    continue
                if line.strip() == 'MIXTURE':
                    break
                if line.split()[1]=='T1': T1 = float(line.split()[0])
                if line.split()[1]=='P1': P1 = float(line.split()[0])
                if line.split()[1]=='dt': dt = float(line.split()[0])
                if line.split()[1]=='L': L = float(line.split()[0])
        u1 = L/dt
        composition = []
        while True:
                line = inpfile.readline()              
                if line == '':
                    break
                composition.append([line.split()[0], float(line.split()[1])])
        return T1, P1, u1, composition
        
def readlist(thermpath):
    """makes list of available gaseous species from therm.dat"""
    specset = set()
    with open(thermpath, 'r') as thermdat:
        while True:
            line = thermdat.readline()
            if len(line) >= 80:
                if line[79] == '1' and line [44] == 'G':
                    specset.add(line[:18].strip())                  
            if line == '':
                break
    thermdat.closed
    #print(elemset)
    return specset

        
def weightcalc(stringofelements):
    """calculates weight from elements string from therm.dat"""
    elements = {'H': 1.0079, 'D': 2.0141, 'C': 12.011, 'O': 15.9994,
                'HE': 4.0026, 'N': 14.0067, 'NE': 20.179, 'AR': 39.948, 'KR': 83.80, 'XE': 131.30,
                'F': 18.9984, 'CL': 35.453, 'BR': 79.904, 'I': 126.904, 'P': 30.9738, 'S': 32.06,
                'AL': 26.9815, 'FE': 55.847, 'CR': 51.996, 'MO': 95.94, 'B': 10.81,  'GA': 69.72}
    weight = 0.0
    for i in range(4):
        s = stringofelements[(0+5*i):(5+5*i)]
        if not (s[0] == ' '):
            weight = weight + elements[(s[0:2]).strip()]*int((s[4:5]).strip())
    return weight
 
def read_thermodat(composition, file):
    """read thermodynamical properties of mixture"""
    mixset = set()
    for i in range(len(composition)):
        mixset.add(composition[i][0])
    thermdat = open('therm.dat', 'r')                       #инициализировать mixture
    for i in range(len(composition)):
        while True:
            line = thermdat.readline()                   #1st line
            if line[:18].strip() in mixset:                         ####!!!
                break			
        name = line[:18].strip()
        weight = weightcalc(line[24:45])                #weight from elements
        low_T = line[45:55]
        threeshold_T = line[55:65]
        high_T = line[65:73]
        a_high = [0,0,0,0,0,0,0]
        a_low = [0,0,0,0,0,0,0]
        line = thermdat.readline() 						#2nd line
        a_high[0] = float(line[0:15])
        a_high[1] = float(line[15:30])
        a_high[2] = float(line[30:45])
        a_high[3] = float(line[45:60])
        a_high[4] = float(line[60:75])
        line = thermdat.readline() 						#3rd line
        a_high[5] = float(line[0:15])
        a_high[6] = float(line[15:30])
        a_low[0] = float(line[30:45])
        a_low[1] = float(line[45:60])
        a_low[2] = float(line[60:75])
        line = thermdat.readline() 						#4th line
        a_low[3] = float(line[0:15])
        a_low[4] = float(line[15:30])
        a_low[5] = float(line[30:45])
        a_low[6] = float(line[45:60])
        print(a_low[5])
        print(len(composition))
    return mixture

weight = weightcalc('C   2H   2AL  1H   10')
print(weight)

specset = readlist('therm.dat')
print(specset)

(T1, P1, u1, composition) = readinpfile('equilib.inp')
print(composition)
print(T1)
print(P1)
print(u1)

mixture = read_thermodat(composition, 'therm.dat')

print(mixture)
