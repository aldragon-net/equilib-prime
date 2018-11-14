from tkinter import *
from tkinter import ttk
from tkinter import font

class Species:
    """class of mixture components"""
    def __init__(self, name, molefrac, weight, low_T,
                 threeshold_T,  high_T, a_low, a_high):
        self.name = name
        self.molefrac = molefrac
        self.weight = weight
        self.low_T = low_T
        self.threeshold_T = threeshold_T
        self.high_T = high_T
        self.a_low = a_low
        self.a_high = a_high

R = 8.31446 #ideal gas constant
k = 1.3806e-23 #Boltzmann consant

def readlist(thermpath):
    """makes set of available gaseous species from therm.dat"""
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
    return specset

def readinpfile(inppath):
    """read input parameters from equilib.inp"""
    with open(inppath, 'r') as inpfile:
        u = 0
        while True:
                line = inpfile.readline()
                if line.strip() == '':
                    continue
                if line.strip() == 'MIXTURE':
                    break
                if line.split()[1]=='T1': T1 = 273.15+float(line.split()[0])
                if line.split()[1]=='P1': P1 = 0.001*float(line.split()[0])
                if line.split()[1]=='dt': dt = float(line.split()[0])
                if line.split()[1]=='L': L = float(line.split()[0])
                if line.split()[1]=='u': u = float(line.split()[0])
        inpfile.closed
        if u == 0:
            u_isw = 1e3*L/dt
        else:
            u_isw = u
        composition = []
        while True:
                line = inpfile.readline()              
                if line == '':
                    break
                composition.append([line.split()[0], float(line.split()[1])])
        return T1, P1, dt, L, u_isw, composition
        
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
        weight = weightcalc(line[24:45])     #weight from elements
        low_T = float(line[45:55])
        high_T = float(line[55:65])
        threeshold_T = float(line[65:73])
        a_high = [0,0,0,0,0,0,0]
        a_low = [0,0,0,0,0,0,0]
        line = thermdat.readline() 			#2nd line
        a_high[0] = float(line[0:15]);  a_high[1] = float(line[15:30])
        a_high[2] = float(line[30:45]); a_high[3] = float(line[45:60])
        a_high[4] = float(line[60:75])
        line = thermdat.readline() 			#3rd line
        a_high[5] = float(line[0:15]); a_high[6] = float(line[15:30])
        a_low[0] = float(line[30:45]); a_low[1] = float(line[45:60])
        a_low[2] = float(line[60:75])
        line = thermdat.readline() 			#4th line
        a_low[3] = float(line[0:15]); a_low[4] = float(line[15:30])
        a_low[5] = float(line[30:45]);a_low[6] = float(line[45:60])
        component = Species(name, molefrac, weight, low_T,
                            threeshold_T,  high_T, a_low, a_high)
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
    T_limit = 6000
    for i in range(len(mixture)):
        if mixture[i].high_T < T_limit:
            T_limit = mixture[i].high_T
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
        if T>=T_limit:
            return T_limit
    return T
    
def mixweight(mixture):
    """mean molecular weight of mixture in a.m.u."""
    weight = 0
    for i in range(len(mixture)):
        weight = weight + mixture[i].molefrac*mixture[i].weight
    return weight

def ISW(mixture, P0, T0, u_isw):
    """parameters behind incident wave"""
    gamma1 = gamma(mixture, T0)
    weight = mixweight(mixture)/1000    #a.m.u. to kg/mole (SI)
    a0 = (gamma1*T0*R/weight)**0.5  
    M_isw = u_isw/a0;
    if M_isw<=1:
        status_label.set('ERROR: ISW velocity too low, no shock')
        return 1, 1, 1, 1, 1, 1, 1, 1, False
    rratio_isw = (gamma1+1)/(gamma1-1)
    Pa, Pb = 1, 2
    js = 1
    dr = 0.01*rratio_isw;
    T_limit = 6000
    for i in range(len(mixture)):
        if mixture[i].high_T < T_limit:
            T_limit = mixture[i].high_T
    while abs(Pa-Pb)>=5E-6:
        H2 = enthalpy(mixture, T0)+(0.5*u_isw**2*(1-(1/rratio_isw)**2))*weight
        T_isw = T_from_enthalpy(mixture, H2)
        if T_isw >= T_limit:
            status_label.set('ERROR: SW too fast, T₂ exceed limit')
            return 1, 1, 1, 1, 1, 1, 1, 1, False
        Pa = rratio_isw*T_isw/T0
        Pb = 1+(u_isw**2*weight/(T0*R))*(1-1/rratio_isw);
        if ((Pa>Pb) and (js<0)) or ((Pa<Pb) and (js>0)):
            dr = dr/4
            js = -js
        if Pa < Pb:
            rratio_isw = rratio_isw + dr
        elif not (Pa==Pb):
            rratio_isw = rratio_isw - dr               
    Pratio=0.5*(Pa+Pb)
    P_isw = P0*Pratio
    n_isw = (P_isw*1e-1)/(k*T_isw)  #bar to Pa, m3 to cm3 => 0.1 factor
    gamma2 = gamma(mixture, T_isw)
    a_isw = (gamma2*T_isw*R/weight)**0.5
    tratio_isw = rratio_isw
    return T_isw, P_isw, n_isw, rratio_isw, u_isw, a_isw, M_isw, tratio_isw, True
    
    
def RSW(mixture, T1, P1, u_isw, T_isw, rratio_isw):
    """parameters behind reflected wave"""
    ok_flag = True
    gamma1 = gamma(mixture, T1)
    weight = mixweight(mixture)/1000
    a1 = (gamma1*T1*R/weight)**0.5;
    M_isw = u_isw/a1;
    Ta = T1*(2*(gamma1-1)*M_isw**2+3-gamma1)*(M_isw**2*(3*gamma1-1)-2*gamma1+2)/((gamma1+1)*M_isw)**2
    Tb = T_isw
    dT = 1
    js = 1
    H2 = enthalpy(mixture, T_isw)
    T_limit = 6000
    for i in range(len(mixture)):
        if mixture[i].high_T < T_limit:
            T_limit = mixture[i].high_T
    while abs(Ta-Tb)>=1E-5:
        if Ta > T_limit:
            status_label.set('ERROR: wave too fast, T₅ exceed limit')
            return 1, 1, 1, 1, 1, 1, 1, False
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
    P_rsw = P1*rratio_rsw*rratio_isw*T_rsw/T1
    n_rsw = (P_rsw*1e-1)/(k*T_rsw)  #bar to Pa, m3 to cm3 => 0.1 factor
    u_rsw_self = u_isw*(1-1/rratio_isw)/(1-1/rratio_rsw)
    u_rsw_lab = u_rsw_self - u_isw*(1-1/rratio_isw)
    gamma2 = gamma(mixture, T_isw)
    gamma5 = gamma(mixture, T_rsw)
    a_isw = (gamma2*T_isw*R/weight)**0.5
    a_rsw = (gamma5*T_rsw*R/weight)**0.5
    M_rsw = u_rsw_self/a_isw 
    return T_rsw, P_rsw, n_rsw, rratio_rsw, u_rsw_lab, a_rsw, M_rsw, ok_flag

#interface_initialization
root = Tk()
root.title("EQUILIB Prime")
root.iconbitmap('oivticon.ico')

capFont = font.Font(family="Helvetica",size=16)
majFont = font.Font(family="Helvetica",size=14,weight="bold")
comFont = font.Font(family="Helvetica",size=14)

mainframe = ttk.Frame(root, padding="10 10 10 10")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

L_label = StringVar(); dt_label = StringVar(); u_label = StringVar();
P1_label = StringVar();T1_label = StringVar(); N_exp_label = StringVar()
spec = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]
frac = [StringVar(), StringVar(), StringVar(), StringVar(), StringVar()]

ttk.Label(mainframe, text="input", font=capFont).grid(column=1, row=1, sticky=W)
ttk.Label(mainframe, text="P1 = ").grid(column=1, row=2, sticky=E)
ttk.Label(mainframe, text="dt = ").grid(column=1, row=3, sticky=E)
ttk.Label(mainframe, text="L = ").grid(column=1, row=4, sticky=E)
ttk.Label(mainframe, text="u = ").grid(column=1, row=5, sticky=E)
ttk.Label(mainframe, text="T1 = ").grid(column=1, row=6, sticky=E)
ttk.Label(mainframe, text="mbar").grid(column=3, row=2, sticky=W)
ttk.Label(mainframe, text="mks").grid(column=3, row=3, sticky=W)
ttk.Label(mainframe, text="mm").grid(column=3, row=4, sticky=W)
ttk.Label(mainframe, text="m/s").grid(column=3, row=5, sticky=W)
ttk.Label(mainframe, text="°C").grid(column=3, row=6, sticky=W)

P1_entry = ttk.Entry(mainframe, width=7, textvariable=P1_label)
P1_entry.grid(column=2, row=2, pady=1, sticky=(W, E))
dt_entry = ttk.Entry(mainframe, width=7, textvariable=dt_label)
dt_entry.grid(column=2, row=3, pady=1, sticky=(W, E))
L_entry = ttk.Entry(mainframe, width=7, textvariable=L_label)
L_entry.grid(column=2, row=4, pady=1, sticky=(W, E))
u_entry = ttk.Entry(mainframe, width=7, textvariable=u_label)
u_entry.grid(column=2, row=5, pady=5, sticky=(W, E))
u_entry = ttk.Entry(mainframe, width=7, textvariable=T1_label)
u_entry.grid(column=2, row=6, pady=5, sticky=(W, E))

s1 = ttk.Separator(mainframe, orient=HORIZONTAL)
s1.grid(column=1, columnspan=3, pady = 0, row=7, sticky=(W, E))

ttk.Label(mainframe, text="Mixture Species").grid(column=1, row=8, sticky=W, columnspan =2)
ttk.Label(mainframe, text="Fraction").grid(column=3, row=8, sticky=W)
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec[0], justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=9, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac[0])
frac_entry.grid(column=3, row=9, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec[1], justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=10, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac[1])
frac_entry.grid(column=3, row=10, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec[2], justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=11, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac[2])
frac_entry.grid(column=3, row=11, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec[3], justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=12, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac[3])
frac_entry.grid(column=3, row=12, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec[4], justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=13, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac[4])
frac_entry.grid(column=3, row=13, pady=1, sticky=(W, E))

s2 = ttk.Separator(mainframe, orient=VERTICAL)
s2.grid(column=4, row=1, rowspan=12, padx=5, sticky=(N, S))

mixture_label = StringVar()
T_isw_label = StringVar(); P_isw_label = StringVar(); u_isw_label = StringVar()
M_isw_label = StringVar(); n_isw_label = StringVar(); rratio_isw_label = StringVar()
a_isw_label = StringVar(); a_rsw_label = StringVar(); tratio_isw_label = StringVar()
T_rsw_label = StringVar(); P_rsw_label = StringVar(); u_rsw_label = StringVar()
M_rsw_label = StringVar(); n_rsw_label = StringVar(); rratio_rsw_label = StringVar()
status_label = StringVar(); flag_color = StringVar(); N_exp_label = StringVar()

flag_color.set('black')

ttk.Label(mainframe, textvariable=T_isw_label, font=majFont).grid(column=5, row=2, rowspan =2,  sticky=W)
ttk.Label(mainframe, textvariable=P_isw_label, font=majFont).grid(column=5, row=4, sticky=W)
ttk.Label(mainframe, textvariable=n_isw_label, font=comFont).grid(column=5, row=6, sticky=W)
ttk.Label(mainframe, textvariable=rratio_isw_label, font=comFont).grid(column=5, row=8, sticky=W)
ttk.Label(mainframe, textvariable=u_isw_label, font=comFont).grid(column=5, row=9, sticky=W)
ttk.Label(mainframe, textvariable=M_isw_label, font=comFont).grid(column=5, row=10, sticky=W)
ttk.Label(mainframe, textvariable=a_isw_label, font=comFont).grid(column=5, row=11, sticky=W)
ttk.Label(mainframe, textvariable=tratio_isw_label, font=comFont).grid(column=5, row=12, sticky=W)

s2 = ttk.Separator(mainframe, orient=VERTICAL)
s2.grid(column=6, row=2, rowspan=11, padx=5, sticky=(N, S))

ttk.Label(mainframe, textvariable=T_rsw_label, font=majFont).grid(column=7, row=2, rowspan =2,  sticky=W)
ttk.Label(mainframe, textvariable=P_rsw_label, font=majFont).grid(column=7, row=4, sticky=W)
ttk.Label(mainframe, textvariable=n_rsw_label, font=comFont).grid(column=7, row=6, sticky=W)
ttk.Label(mainframe, textvariable=rratio_rsw_label, font=comFont).grid(column=7, row=8, sticky=W)
ttk.Label(mainframe, textvariable=u_rsw_label, font=comFont).grid(column=7, row=9, sticky=W)
ttk.Label(mainframe, textvariable=M_rsw_label, font=comFont).grid(column=7, row=10, sticky=W)
ttk.Label(mainframe, textvariable=a_rsw_label, font=comFont).grid(column=7, row=11, sticky=W)

mainframe.columnconfigure(5, weight=1)
mainframe.columnconfigure(7, weight=1)

ttk.Label(mainframe, textvariable=mixture_label, font=comFont).grid(column=5, columnspan=3, row=1, sticky=W)

statusbar = ttk.Frame(mainframe, borderwidth=0, relief='sunken', padding="2 0 0 0")
statusbar.grid(column=5, columnspan=3, row=13, sticky=(W,N,S,E))
statusbar.columnconfigure(0, weight=1)
statusbar.rowconfigure(0, weight=1)
ttk.Label(statusbar, textvariable=status_label).grid(column=0, row=0)
ttk.Label(statusbar, text='N = ').grid(column=1, row=0)
N_exp_entry = ttk.Entry(statusbar, width=4, textvariable=N_exp_label)
N_exp_entry.grid(column=2, row=0, sticky=(E))

#end of interface initialization

def input_file_write(inppath, T1, P1, dt, L, u_isw, composition):
    """writes input parameters to equilib.inp"""
    with open(inppath, 'w') as inpfile:
        inpfile.write(str(T1-273.15) + '\tT1 [C]'+'\n')
        inpfile.write(str(P1*1000) + '\tP1 [mbar]'+'\n')
        inpfile.write(str(dt) + '\tdt [mks]'+'\n')
        inpfile.write(str(L) + '\tL [mm]'+'\n')
        inpfile.write('\n')
        inpfile.write('MIXTURE'+'\n')
        for i in range(len(composition)):
            if i+1 < len(composition):
                inpfile.write(composition[i][0] + '\t' + str(composition[i][1])+'\n')
            else:
                inpfile.write(composition[i][0]+ '\t' + str(composition[i][1]))
    inpfile.close()
    return

def input_rewrite(T1, P1, dt, L, u_isw, composition):
    """write input labels from file"""
    def nice_mix(composition):
        for i in range(len(composition)):
            if composition[i][0] in ('AR', 'HE', 'NE', 'KR'):
                composition[i] , composition[-1] = composition[-1], composition[i]
            if composition[i][0] in ('O2'):
                composition[i] , composition[-2] = composition[-2], composition[i]
        mixstring='for SW in '
        for i in range(len(composition)):
            if i>0:
                mixstring=mixstring + ' + ' 
            if composition[i][1] >= 0.1:
                mixstring=mixstring+"{:.1f}".format(composition[i][1])+'% '+composition[i][0]
            elif composition[i][1] >= 0.0001:
                mixstring=mixstring+"{:.0f}".format(composition[i][1]*1e4)+'ppm '+composition[i][0]
            else:
                mixstring=mixstring+"{:.2f}".format(composition[i][1]*1e4)+'ppm '+composition[i][0]
        submixstring = mixstring[0]
        for i in range(len(mixstring)-1):
            if (mixstring[i+1].isdigit()) and (submixstring[i] not in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ' ', '.')):
                if mixstring[i+1] == '0':
                    submixstring = submixstring+'₀'
                if mixstring[i+1] == '1':
                    submixstring = submixstring+'₁'
                if mixstring[i+1] == '2':
                    submixstring = submixstring+'₂'
                if mixstring[i+1] == '3':
                    submixstring = submixstring+'₃'
                if mixstring[i+1] == '4':
                    submixstring = submixstring+'₄'
                if mixstring[i+1] == '5':
                    submixstring = submixstring+'₅'
                if mixstring[i+1] == '6':
                    submixstring = submixstring+'₆'
                if mixstring[i+1] == '7':
                    submixstring = submixstring+'₇'
                if mixstring[i+1] == '8':
                    submixstring = submixstring+'₈'
                if mixstring[i+1] == '9':
                    submixstring = submixstring+'₉'
            else:
                submixstring = submixstring+mixstring[i+1]
        mixstring = submixstring
        return mixstring

    P1_label.set(str(P1*1000))
    dt_label.set(str(dt))
    L_label.set(str(L))
    T1_label.set(str(T1-273.15))
    for i in range(5):
        if i<len(composition):
            spec[i].set(composition[i][0])
            frac[i].set(composition[i][1])
        else:
            spec[i].set('')
            frac[i].set('')
    mixstring = nice_mix(composition)
    mixture_label.set(mixstring)
    return

def output_rewrite(T_isw, P_isw, n_isw, rratio_isw, u_isw, a_isw, M_isw, tratio_isw,
              T_rsw, P_rsw, n_rsw, rratio_rsw, u_rsw, a_rsw, M_rsw):
    """rewrite output labels"""
    T_isw_label.set('T₂ = ' + "{:.0f}".format(T_isw) + ' K')
    P_isw_label.set('P₂ = '+ "{:.3f}".format(P_isw) + ' bar')
    n_isw_label.set('n₂ = '+ "{:.3e}".format(n_isw) + ' cm⁻³')
    rratio_isw_label.set('ρ₂/ρ₁ = '+ "{:.3f}".format(rratio_isw))
    u_isw_label.set('u₂ = '+ "{:.0f}".format(u_isw) + ' m/s')
    M_isw_label.set('M₂ = '+ "{:.3f}".format(M_isw))
    a_isw_label.set('a₂ = '+ "{:.0f}".format(a_isw) + ' m/s')
    tratio_isw_label.set('τ = '+ "{:.3f}".format(tratio_isw) + 't')
    T_rsw_label.set('T₅ = ' + "{:.0f}".format(T_rsw) + ' K')
    P_rsw_label.set('P₅ = '+ "{:.3f}".format(P_rsw) + ' bar')
    n_rsw_label.set('n₅ = '+ "{:.3e}".format(n_rsw) + ' cm⁻³')
    rratio_rsw_label.set('ρ₅/ρ₁ = '+ "{:.3f}".format(rratio_rsw*rratio_isw))
    u_rsw_label.set('u₅ = '+ "{:.0f}".format(u_rsw) + ' m/s')
    M_rsw_label.set('M₅ = '+ "{:.3f}".format(M_rsw))
    a_rsw_label.set('a₅ = '+ "{:.0f}".format(a_rsw) + ' m/s')
    status_label.set('status: ok')
    return

def output_fail():
    """rewrite output labels stubs if calculatuion failed"""
    T_isw_label.set('T₂ =  . . . . .  K'); P_isw_label.set('P₂ =  . . . . .  bar')
    n_isw_label.set('n₂ = ..... cm⁻³'); rratio_isw_label.set('ρ₂/ρ₁ = .....')
    u_isw_label.set('u₂ = ..... m/s'); M_isw_label.set('M₂ = .....')
    a_isw_label.set('a₂ = ..... m/s'); tratio_isw_label.set('τ = ..... t')
    T_rsw_label.set('T₅ =  . . . . .  K'); P_rsw_label.set('P₅ = . . . . . bar')
    n_rsw_label.set('n₅ = ..... cm⁻³'); rratio_rsw_label.set('ρ₅/ρ₁ = .....')
    u_rsw_label.set('u₅ = ..... m/s'); M_rsw_label.set('M₅ = .....')
    a_rsw_label.set('a₅ = ..... m/s')
    return

def closeapp(event, *args):
    global root
    root.destroy()
    return
    
def change_flag(event, *args):
    status_label.set('Input changed. Press Ctrl+Enter to recalculate')
    return

def readGUIinp(T1_label, P1_label, dt_label, L_label, u_label, spec, frac, N_exp_label):
    """read input parameters from GUI"""
    ok_flag = True
    try:
        T1 = 273.15+float(T1_label.get())
    except ValueError:
        status_label.set('ERROR: incorrect T1 value'); ok_flag = False
    if (T1<200) or (T1>6000):
        status_label.set('ERROR: extreme T1 value'); ok_flag = False
    try:
        P1 = 0.001*float(P1_label.get())
    except ValueError:
        status_label.set('ERROR: incorrect P1 value'); ok_flag = False
    if (P1<=0):
        status_label.set('ERROR: incorrect P1 value'); ok_flag = False

    if not (dt_label.get()==''):
        try:
            dt = float(dt_label.get())
        except ValueError:
            status_label.set('ERROR: incorrect dt value'); ok_flag = False
        if (dt<=0):
            status_label.set('ERROR: incorrect dt value'); ok_flag = False

    if not (L_label.get()==''):
        try:
            L = float(L_label.get())
        except ValueError:
            status_label.set('ERROR: incorrect L value'); ok_flag = False
        if (L<=0):
            status_label.set('ERROR: incorrect L value'); ok_flag = False


    if not (u_label.get()==''):
        try:
            u_isw = float(u_label.get())
        except ValueError:
            status_label.set('ERROR: incorrect u value'); ok_flag = False
        if (u_isw<=0):
            status_label.set('ERROR: incorrect u value'); ok_flag = False
    else:
        u_isw = 1e3*L/dt

    composition = []
    for i in range(5):
        if not (spec[i].get()==''):
            if spec[i].get() not in specset:
                status_label.set('ERROR: no \"' + spec[i].get()+ '\" in therm.dat')
                ok_flag = False
                continue
            try:
                composition.append([spec[i].get(), float(frac[i].get())])
            except ValueError:
                status_label.set('ERROR: incorrect \"' + spec[i].get()+ '\" fraction')
                ok_flag = False

      
    sum = 0
    for i in range(len(composition)): sum = sum + composition[i][1]
    for i in range(len(composition)): composition[i][1] = composition[i][1]*100/sum

    N_exp = N_exp_label.get()
    
    return T1, P1, dt, L, u_isw, composition, N_exp, ok_flag


def outfile_write(outpath, N_exp, P1,  u_isw, T1, P_isw, T_isw, P_rsw, T_rsw, rratio_isw, rratio_rsw, composition):
    """writes otuput data to equilib.dat"""
    try:
        outfile = open(outpath, 'r')
    except:
        outfile = open(outpath, 'w')
        outfile.write('N_exp '+' P1   '+'u_isw  '+' T2   '+' P2   '+' T5   '+
                      ' P5   '+'ro2ro1 '+'ro5ro1 '+'   n1    '+'   n2     '+'   n5')
    outfile.close
    with open(outpath, 'a') as outfile:
        outfile.write('\n'+(N_exp+10*' ')[:6])
        outfile.write(('{:.3f}'.format(P1)+10*' ')[:7])
        outfile.write(('{:.0f}'.format(u_isw)+10*' ')[:6])
        outfile.write(('{:.0f}'.format(T_isw)+3*' ')[:6])
        outfile.write(('{:.3f}'.format(P_isw)+3*' ')[:6])
        outfile.write(('{:.0f}'.format(T_rsw)+3*' ')[:6])
        outfile.write(('{:.3f}'.format(P_rsw)+3*' ')[:6])
        outfile.write(('{:.3f}'.format(rratio_isw)+4*' ')[:7])
        outfile.write(('{:.3f}'.format(rratio_rsw*rratio_isw)+4*' ')[:7])
        outfile.write(('{:.3e}'.format(P1*1e-1/(k*T1))+4*' ')[:10])
        outfile.write(('{:.3e}'.format(P_isw*1e-1/(k*T_isw))+4*' ')[:10])
        outfile.write(('{:.3e}'.format(P_rsw*1e-1/(k*T_rsw))+4*' ')[:10])
    outfile.close()
    return    

def restart(event, *args):
    """recalculation on Ctrl+Enter"""
    T1, P1, dt, L, u_isw, composition, N_exp, ok_flag = readGUIinp(T1_label, P1_label, dt_label, L_label, u_label, spec, frac, N_exp_label)
    if not ok_flag:
        return
    input_rewrite(T1, P1, dt, L, u_isw, composition)
    mixture = read_thermodat(composition, 'therm.dat')
    (T_isw, P_isw, n_isw, rratio_isw, u_isw, a_isw, M_isw, tratio_isw, ok_flag) = ISW(mixture, P1, T1, u_isw)
    if ok_flag:
        (T_rsw, P_rsw, n_rsw, rratio_rsw, u_rsw, a_rsw, M_rsw, ok_flag) = RSW(mixture, T1, P1, u_isw, T_isw, rratio_isw)
    if not ok_flag:
        output_fail()
    else:
        output_rewrite(T_isw, P_isw, n_isw, rratio_isw, u_isw, a_isw, M_isw, tratio_isw,
              T_rsw, P_rsw, n_rsw, rratio_rsw, u_rsw, a_rsw, M_rsw)
           
    if not (N_exp == ''):
        outfile_write('equilib.dat',  N_exp, P1,  u_isw, T1, P_isw, T_isw, P_rsw, T_rsw, rratio_isw, rratio_rsw, composition)
        N_exp_label.set('')
        status_label.set('status: Ok\t \''+ N_exp+'\' was written to equilib.out')
                
    input_file_write('equilib.inp', T1, P1, dt, L, u_isw, composition)

    return

#mainbody start

specset = readlist('therm.dat')
(T1, P1, dt, L, u_isw, composition) = readinpfile('equilib.inp')

input_rewrite(T1, P1, dt, L, u_isw, composition)

T1, P1, dt, L, u_isw, composition, N_exp, ok_flag = readGUIinp(T1_label, P1_label, dt_label, L_label, u_label, spec, frac, N_exp_label)

mixture = read_thermodat(composition, 'therm.dat')

(T_isw, P_isw, n_isw, rratio_isw, u_isw, a_isw, M_isw, tratio_isw, ok_flag) = ISW(mixture, P1, T1, u_isw)
if ok_flag:
    (T_rsw, P_rsw, n_rsw, rratio_rsw, u_rsw, a_rsw, M_rsw, ok_flag) = RSW(mixture, T1, P1, u_isw, T_isw, rratio_isw)
if not ok_flag:
    output_fail()
else:
    output_rewrite(T_isw, P_isw, n_isw, rratio_isw, u_isw, a_isw, M_isw, tratio_isw,
              T_rsw, P_rsw, n_rsw, rratio_rsw, u_rsw, a_rsw, M_rsw)

P1_label.trace("w", change_flag); T1_label.trace("w", change_flag)
L_label.trace("w", change_flag); dt_label.trace("w", change_flag)
u_label.trace("w", change_flag); N_exp_label.trace("w", change_flag)
frac[0].trace("w", change_flag); spec[0].trace("w", change_flag)
frac[1].trace("w", change_flag); spec[1].trace("w", change_flag)
frac[2].trace("w", change_flag); spec[2].trace("w", change_flag)
frac[3].trace("w", change_flag); spec[3].trace("w", change_flag)
frac[4].trace("w", change_flag); spec[4].trace("w", change_flag)

root.bind("<Control-Return>", restart)
root.bind("<Escape>", closeapp)

root.mainloop()



