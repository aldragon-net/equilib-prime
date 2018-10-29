from tkinter import *
from tkinter import ttk
from tkinter import font


root = Tk()
root.title("EQUILIB Prime")
root.iconbitmap('oivticon.ico')


capFont = font.Font(family="Helvetica",size=18,weight="bold")
majFont = font.Font(family="Helvetica",size=14,weight="bold")
comFont = font.Font(family="Helvetica",size=14)

mainframe = ttk.Frame(root, padding="10 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

L = StringVar()
dt = StringVar()
u = StringVar()
spec1 = StringVar()
spec2 = StringVar()
spec3 = StringVar()
spec4 = StringVar()
spec5 = StringVar()
frac1 = StringVar()
frac2 = StringVar()
frac3 = StringVar()
frac4 = StringVar()
frac5 = StringVar()


ttk.Label(mainframe, text="ISW", font=capFont).grid(column=1, row=1, sticky=W)
ttk.Label(mainframe, text="dt = ").grid(column=1, row=2, sticky=E)
ttk.Label(mainframe, text="L = ").grid(column=1, row=3, sticky=E)
ttk.Label(mainframe, text="u = ").grid(column=1, row=4, sticky=E)
ttk.Label(mainframe, text="mks").grid(column=3, row=2, sticky=W)
ttk.Label(mainframe, text="mm").grid(column=3, row=3, sticky=W)
ttk.Label(mainframe, text="m/s").grid(column=3, row=4, sticky=W)

L_entry = ttk.Entry(mainframe, width=7, textvariable=L)
L_entry.grid(column=2, row=3, pady=1, sticky=(W, E))
dt_entry = ttk.Entry(mainframe, width=7, textvariable=dt)
dt_entry.grid(column=2, row=2, pady=1, sticky=(W, E))
u_entry = ttk.Entry(mainframe, width=7, textvariable=u)
u_entry.grid(column=2, row=4, pady=5, sticky=(W, E))

s1 = ttk.Separator(mainframe, orient=HORIZONTAL)
s1.grid(column=1, columnspan=3, pady = 10, row=5, sticky=(W, E))

ttk.Label(mainframe, text="Mixture").grid(column=1, row=6, sticky=W, columnspan =3)
ttk.Label(mainframe, text="Species").grid(column=2, row=7, sticky=W)
ttk.Label(mainframe, text="Fraction").grid(column=3, row=7, sticky=W)
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec1, justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=8, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec2, justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=9, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec3, justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=10, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec4, justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=11, pady=1, sticky=(W, E))
spec_entry = ttk.Entry(mainframe, width=7, textvariable=spec5, justify=RIGHT)
spec_entry.grid(column=1, columnspan =2,  row=12, pady=1, sticky=(W, E))

frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac1)
frac_entry.grid(column=3, row=8, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac2)
frac_entry.grid(column=3, row=9, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac3)
frac_entry.grid(column=3, row=10, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac4)
frac_entry.grid(column=3, row=11, pady=1, sticky=(W, E))
frac_entry = ttk.Entry(mainframe, width=7, textvariable=frac5)
frac_entry.grid(column=3, row=12, pady=1, sticky=(W, E))


s2 = ttk.Separator(mainframe, orient=VERTICAL)
s2.grid(column=4, row=1, rowspan=12, padx=5, sticky=(N, S))

T_isw_label = StringVar()
P_isw_label = StringVar()
u_isw_label = StringVar()
M_isw_label = StringVar()
n_isw_label = StringVar()
rratio_isw_label = StringVar()


T_rsw_label = StringVar()
P_rsw_label = StringVar()
u_rsw_label = StringVar()
M_rsw_label = StringVar()
n_rsw_label = StringVar()
rratio_rsw_label = StringVar()

T_isw_label.set('T₂ = 1200 K')
T_rsw_label.set('T₅ = 2500 K')
P_isw_label.set('P₂ = 2.567 bar')
P_rsw_label.set('P₅ = 6.789 bar')
u_isw_label.set('u₂ = 1200 m/s')
u_rsw_label.set('u₅ = 456 m/s')
M_isw_label.set('M₂ = 3.123')
M_rsw_label.set('M₅ = 1.456')
n_rsw_label.set('n₅ = 2.65E19')
n_isw_label.set('n₂ = 1.34E19')
rratio_isw_label.set('ρ₂/ρ₁ = 1.23')
rratio_rsw_label.set('ρ₅/ρ₁ = 3.23')

ttk.Label(mainframe, textvariable=T_isw_label, font=majFont).grid(column=5, row=2, rowspan =2,  sticky=W)

ttk.Label(mainframe, textvariable=P_isw_label, font=majFont).grid(column=5, row=4, sticky=W)

ttk.Label(mainframe, textvariable=u_isw_label, font=comFont).grid(column=5, row=6, sticky=W)
ttk.Label(mainframe, textvariable=M_isw_label, font=comFont).grid(column=5, row=7, sticky=W)
ttk.Label(mainframe, textvariable=n_isw_label, font=comFont).grid(column=5, row=8, sticky=W)
ttk.Label(mainframe, textvariable=rratio_isw_label, font=comFont).grid(column=5, row=9, sticky=W)

s2 = ttk.Separator(mainframe, orient=VERTICAL)
s2.grid(column=6, row=2, rowspan=11, padx=5, sticky=(N, S))

ttk.Label(mainframe, textvariable=T_rsw_label, font=majFont).grid(column=7, row=2, rowspan =2,  sticky=W)

ttk.Label(mainframe, textvariable=P_rsw_label, font=majFont).grid(column=7, row=4, sticky=W)

ttk.Label(mainframe, textvariable=u_rsw_label, font=comFont).grid(column=7, row=6, sticky=W)
ttk.Label(mainframe, textvariable=M_rsw_label, font=comFont).grid(column=7, row=7, sticky=W)
ttk.Label(mainframe, textvariable=n_rsw_label, font=comFont).grid(column=7, row=8, sticky=W)
ttk.Label(mainframe, textvariable=rratio_rsw_label, font=comFont).grid(column=7, row=9, sticky=W)

mainframe.columnconfigure(5, weight=1)
mainframe.columnconfigure(7, weight=1)

ttk.Label(mainframe, text="for SW in 2%Fe(CO)5 + 2%C2H2 + 5%O2 + 80%Ar:", font=comFont).grid(column=5, columnspan=3, row=1, sticky=W)

root.update()
root.minsize(root.winfo_width(), root.winfo_height())
#root.maxsize(root.winfo_width(), root.winfo_height())

root.mainloop()
