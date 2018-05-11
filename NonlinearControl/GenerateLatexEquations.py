import numpy as np
import os
TxtFileName = "equations.txt"
with open(TxtFileName,"r") as mytxtfile:
    Data = mytxtfile.readlines()
for j in range(len(Data)):
    NewEntry = Data[j]
    NewEntry = NewEntry.replace("=","&=")
    NewEntry = NewEntry.replace(" lambda X,U:","")
    NewEntry = NewEntry.replace(" lambda X:","")
    NewEntry = NewEntry.replace("np.cos","\\cos")
    NewEntry = NewEntry.replace("np.sin","\\sin")
    NewEntry = NewEntry.replace("np.tan","\\tan")
    NewEntry = NewEntry.replace("np.exp","\\text{exp}")
    NewEntry = NewEntry.replace("np.log","\\text{ln}")
    NewEntry = NewEntry.replace("np.sign","\\text{sgn}")
    NewEntry = NewEntry.replace("dx","\\dot{x}_")
    NewEntry = NewEntry.replace("d2x","\\ddot{x}_")
    NewEntry = NewEntry.replace("X[","x[")
    NewEntry = NewEntry.replace("X","\\vec{x}")
    NewEntry = NewEntry.replace("U","u")
    for i in reversed(range(12)):
        NewEntry = NewEntry.replace("**"+str(i+1),"^{"+str(i+1)+"}")
        NewEntry = NewEntry.replace("["+str(i)+"]","_{"+str(i+1)+"}")
        NewEntry = NewEntry.replace("g"+str(i+1),"g_{"+str(i+1)+"}")
        NewEntry = NewEntry.replace("c"+str(i+1),"c_{"+str(i+1)+"}")
    NewEntry = NewEntry.replace("*","\\cdot ")
    NewEntry = NewEntry.replace("lambdify([θ_EFE],r1.subs([(θ_PS,np.pi/2)]))(x_{1})","r_{1}(x_{1})")
    NewEntry = NewEntry.replace("lambdify([θ_EFE],dr1.subs([(θ_PS,np.pi/2)]))(x_{1})",\
                                "r'_{1}(x_{1})")
    NewEntry = NewEntry.replace("lambdify([θ_EFE],d2r1.subs([(θ_PS,np.pi/2)]))(x_{1})",\
                                "r''_{1}(x_{1})")
    NewEntry = NewEntry.replace("lambdify([θ_EFE],r2.subs([(θ_PS,np.pi/2)]))(x_{1})","r_{2}(x_{1})")
    NewEntry = NewEntry.replace("lambdify([θ_EFE],dr2.subs([(θ_PS,np.pi/2)]))(x_{1})",\
                                "r'_{2}(x_{1})")
    NewEntry = NewEntry.replace("lambdify([θ_EFE],d2r2.subs([(θ_PS,np.pi/2)]))(x_{1})",\
                                "r''_{2}(x_{1})")
    if "&=" not in NewEntry:
        NewEntry = "&\\hspace{1.5em} " + NewEntry
    NewEntry = NewEntry.replace("\\\n","\\\\")
    NewEntry = NewEntry.replace("\n","\\\\")
    Data[j] = NewEntry

TexFileName = "equation.tex"
latexfile = open(TexFileName,"w")
latexfile.write("\\documentclass[a4paper]{article} \n" + \
                    "\\usepackage[english]{babel} \n" + \
                    "\\usepackage[utf8x]{inputenc} \n" + \
                    "\\usepackage[T1]{fontenc} \n" + \
                    "\\usepackage{ragged2e} \n" + \
                    "\\usepackage{amsmath} \n" + \
                    "\\usepackage[a4paper,top=3cm,bottom=2cm,left=3cm,right=3cm,marginparwidth=1.75cm]{geometry} \n" + \
                    "\\begin{document} \n" + \
                    "\\section*{"+ TexFileName + "} \n" + \
                    "\\allowdisplaybreaks \n" + \
                    "\\begin{flalign*} \n")
[latexfile.write(el) for el in Data]
latexfile.write("\\end{flalign*} \n")
latexfile.write("\\end{document}")
latexfile.close()

os.system("pdflatex " + TexFileName + " &> /dev/null")
