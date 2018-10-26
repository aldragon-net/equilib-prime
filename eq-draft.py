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

R = 8.31446 #ideal gas constant
        
        
def readinpfile(inppath):
    """read input parameters from equilib.inp"""
    with open(inppath, 'r') as inpfile:
        while True:
                line = inpfile.readline()
                if line.strip() == '':
                    continue
                if line.strip() == 'MIXTURE':
                    break
                if line.split()[1]=='T0': T0 = 273.15+float(line.split()[0])
                if line.split()[1]=='P0': P0 = float(line.split()[0])
                if line.split()[1]=='dt': dt = float(line.split()[0])
                if line.split()[1]=='L': L = float(line.split()[0])
        u_isw = 1e3*L/dt
        composition = []
        while True:
                line = inpfile.readline()              
                if line == '':
                    break
                composition.append([line.split()[0], float(line.split()[1])])
        return T0, P0, u_isw, composition
        
def readlist(thermpath):
    """makes list of available gaseous species from therm.dat"""
    specset = set()
    with open(thermpath, 'r') as thermdat:
        while True:
            line = thermdat.readline()
            if len(line) >= 80:
                if line[79] == '1' and line [44] == 'G':
                    specset.add(line[:18].split()[0])                  
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
    mixture = []
    for i in range(len(composition)):
        thermdat = open('therm.dat', 'r')
        while True:
            line = thermdat.readline()                   #1st line
            if line[:18].split()[0] == composition[i][0]:
                break			
        name = composition[i][0]
        molefrac = 0.01*composition[i][1]
        weight = weightcalc(line[24:45])                #weight from elements
        low_T = float(line[45:55])
        high_T = float(line[55:65])
        threeshold_T = float(line[65:73])
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
        component = Species(name, molefrac, weight, low_T, threeshold_T,  high_T, a_low, a_high)
        mixture.append(component)
    return mixture

def Cp_calc(mixture, T):
    """heat capacity of given mixture at given T"""
    Cp = 0.0
    for i in range(len(mixture)):
        if T <= mixture[i].threeshold_T:
            (a1, a2, a3, a4, a5) = mixture[i].a_low[:5]
        else:
            (a1, a2, a3, a4, a5) = mixture[i].a_high[:5]
        Cp = Cp + mixture[i].molefrac*(R*(a1+a2*T+a3*T**2+a4*T**3+a5*T**4))
    return Cp

def gamma(mixture, T):
    """gamma of given mixture at given T"""
    Cp = Cp_calc(mixture, T)
    gamma = Cp/(Cp-R)
    return gamma
    
def enthalpy(mixture, T):
    """enthalpy of given mixture at given T"""
    H = 0.0
    for i in range(len(mixture)):
        if T <= mixture[i].threeshold_T:
            (a1, a2, a3, a4, a5, a6) = mixture[i].a_low[:6]
        else:
            (a1, a2, a3, a4, a5, a6) = mixture[i].a_high[:6]
        H = H + mixture[i].molefrac*(R*(a1*T+(a2/2)*T**2+(a3/3)*T**3+(a4/4)*T**4+(a5/5)*T**5+a6))
    return H
    
def T_from_enthalpy(mixture, H):
    """T corresponding given enthalpy of given mixture"""
    T = 250
    dT = 200
    sign = 1
    while abs(H-enthalpy(mixture, T)) >= 0.01:
        if H-enthalpy(mixture, T) > 0:
            if sign == -1:
               sign = 1
               dT = dT/2
            T = T + dT               
        else:
            if sign == 1:
               sign = -1
               dT = dT/2
            T = T - dT
    return T
    
def mixweight(mixture):
    weight = 0
    for i in range(len(mixture)):
        weight = weight + mixture[i].molefrac*mixture[i].weight
    return weight

def ISW(mixture, P0, T0, u_isw):
    """parameters behind incident wave"""
    gamma1 = gamma(mixture, T0)
    weight = mixweight(mixture)/1000
    a0 = (gamma1*T0*R/weight)**0.5;
    M_isw = u_isw/a0;
    rratio_isw = (gamma1+1)/(gamma1-1)
    Pa, Pb = 1, 2
    js = 1
    dr = 0.01*rratio_isw;
    print(enthalpy(mixture, T0))
    while abs(Pa-Pb)>=1E-8:
        H2 = enthalpy(mixture, T0)+(0.5*u_isw**2*(1-(1/rratio_isw)**2))*weight
        T_isw = T_from_enthalpy(mixture, H2)
        Pa = rratio_isw*T_isw/T0
        Pb = 1+(u_isw**2*weight/(T0*R))*(1-1/rratio_isw);
        if ((Pa>Pb) and (js<0)) or ((Pa<Pb) and (js>0)):
            dr = dr/4
            js = -js
        if Pa < Pb:
            rratio_isw = rratio_isw + dr
        else:
            rratio_isw = rratio_isw - dr  
                    
    Pratio=0.5*(Pa+Pb)
    P_isw = P0*Pratio
    
    return P_isw,T_isw, u_isw, M_isw, rratio_isw
    
    
def RSW(mixture, T0, P0, u_isw, T_isw):
    """parameters behind reflected wave"""
    gamma1 = gamma(mixture, T0)
    weight = mixweight(mixture)/1000
    a1 = (gamma1*T0*R/weight)**0.5;
    print ('a1=', a1, gamma1, T0, R, weight)
    M_isw = u_isw/a1;
    Ta = T0*(2*(gamma1-1)*M_isw**2+3-gamma1)*(M_isw**2*(3*gamma1-1)-2*gamma1+2)/((gamma1+1)*M_isw)**2
    Tb = T_isw
    dT = 1
    js = 1
    H2 = enthalpy(mixture, T_isw)
    while abs(Ta-Tb)>=1E-5:
        H5 = enthalpy(mixture, Ta)
        s0 = 0.5*(u_isw*(1-1/rratio_isw))**2
        s1 = (H5-H2)/weight
        Tb = (s1-s0)*(T_isw/(s1+s0)+weight/R)
        if ((Ta>Tb) and (js<0)) or ((Ta<Tb) and (js>0)):
            dT = dT/10
            js = -js
        if Ta > Tb:
            Ta = Ta + dT
        else:
            Ta = Ta - dT
    T_rsw = 0.5*(Ta+Tb);
    rratio_rsw = (s1+s0)/(s1-s0);
    P_rsw = P0*rratio_rsw*rratio_isw*T_rsw/T0
    u_rsw = u_isw*(1-1/rratio_isw)/(1-1/rratio_rsw)
    u_rsw_lab = u_rsw - u_isw*(1-1/rratio_isw)
    return P_rsw, T_rsw, rratio_rsw, u_rsw_lab
            
    
specset = readlist('therm.dat')
(T0, P0, u_isw, composition) = readinpfile('equilib.inp')
mixture = read_thermodat(composition, 'therm.dat')
(P_isw,T_isw, u_isw, M_isw, rratio_isw) = ISW(mixture, P0, T0, u_isw)
(P_rsw,T_rsw, rratio_rsw, u_rsw_lab) = RSW(mixture, T0, P0, u_isw, T_isw)

print('ISW:', P_isw, T_isw, u_isw, M_isw, rratio_isw)
print('RSW:', P_rsw, T_rsw, rratio_rsw, u_rsw_lab)


