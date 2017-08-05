# Moment arm functions taken from Rankin 2012.

# def MA_function(Coefficients):
#
def F1(q1,Coefficients):
    assert len(Coefficients)>=5, "Insufficient number of coefficients for F1 (5)"
    import numpy as np
    import sympy as sp
    output = (sp.Matrix(Coefficients[:5]).T*sp.Matrix([1,q1,q1**2,q1**3,q1**4]))[0,0]
    return(output)
def F2(q1,q2,Coefficients):
    assert len(Coefficients)>=7, "Insufficient number of coefficients for F2 (7)"
    import numpy as np
    import sympy as sp
    output = (sp.Matrix(Coefficients[:7]).T*sp.Matrix([1,q1,q1**2,q1**3,q2,q2**2,q2**3]))[0,0]
    return(output)
def F3(q1,q2,Coefficients):
    assert len(Coefficients)>=11, "Insufficient number of coefficients for F3 (11)"
    import numpy as np
    import sympy as sp
    output = F2(q1,q2,Coefficients[:7]) + (sp.Matrix(Coefficients[7:11]).T*sp.Matrix([q1**4,q1**5,q2**4,q2**5]))[0,0]
    return(output)
def F4(q1,q2,Coefficients):
    assert len(Coefficients)>=12, "Insufficient number of coefficients for F4 (12)"
    import numpy as np
    import sympy as sp
    output = F2(q1,q2,Coefficients[:7]) + (sp.Matrix(Coefficients[7:12]).T*sp.Matrix([q1*q2,(q1**2)*q2,q1*(q2**2),(q1**3)*q2,q1*(q2**3)]))[0,0]
    return(output)
def F5(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=20, "Insufficient number of coefficients for F5 (20)"
    import numpy as np
    import sympy as sp
    output = (sp.Matrix(Coefficients[:20]).T*\
                sp.Matrix([1,q1,q2,q3,\
                            q1**2,q2**2,q3**2,\
                            q1**3,q2**3,q3**3,\
                            q1*q2,q1*q3,q2*q3,\
                            q1*q2*q3,(q1**2)*q2,(q1**2)*q3,\
                            q1*(q2**2),(q2**2)*q3,q1*(q3**2),\
                            q2*(q3**2)]))[0,0]
    return(output)
def F6(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=27, "Insufficient number of coefficients for F6 (27)"
    import numpy as np
    import sympy as sp
    output =F5(q1,q2,q3,Coefficients) + \
            (sp.Matrix(Coefficients[20:27]).T*\
                sp.Matrix([ (q1**3)*q2,(q1**3)*q3,\
                            q1*(q2**3),(q2**3)*q3,q1*(q3**3),\
                            q2*(q3**3),(q1*q2*q3)**2]))[0,0]
    return(output)
def F7(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=42, "Insufficient number of coefficients for F7 (42)"
    import numpy as np
    import sympy as sp
    output =F6(q1,q2,q3,Coefficients) + \
            (sp.Matrix(Coefficients[27:42]).T*\
                sp.Matrix([ (q1*q2)**2,(q1*q3)**2,(q2*q3)**2,
                            (q1**4)*q2,(q1**4)*q3,q1*(q2**4),\
                            (q2**4)*q3,q1*(q3**4),q2*(q3**4),\
                            (q1**3)*(q2**2),(q1**3)*(q3**2),\
                            (q1**2)*(q2**3),(q2**3)*(q3**2),\
                            (q1**2)*(q3**3),(q2**2)*(q3**3)]))[0,0]
    return(output)
def F8(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=72, "Insufficient number of coefficients for F8 (72)"
    import numpy as np
    import sympy as sp
    output =F6(q1,q2,q3,Coefficients) + \
            (sp.Matrix(Coefficients[27:72]).T*\
                sp.Matrix([ (q1**2)*q2*q3,q1*(q2**2)*q3,q1*q2*(q3**2),\
                            (q1**3)*q2*q3,q1*(q2**3)*q3,q1*q2*(q3**3),\
                            (q1**3)*(q2**2),(q1**3)*(q3**2),\
                            (q1**2)*(q2**3),(q2**3)*(q3**2),\
                            (q1**2)*(q3**3),(q2**2)*(q3**3),\
                            (q1**3)*(q2**2)*q3,(q1**3)*q2*(q3**2),\
                            (q1**2)*(q2**3)*q3,q1*(q2**3)*(q3**2),\
                            (q1**2)*q2*(q3**3),q1*(q2**2)*(q3**3),\
                            (q1**4),(q2**4),(q3**4),\
                            (q1**4)*q2,q1*(q2**4),\
                            q1*(q3**4),(q1**4)*q3,\
                            (q2**4)*q3,q2*(q3**4),\
                            (q1**4)*(q2**2),(q1**2)*(q2**4),\
                            (q1**2)*(q3**4),(q1**4)*(q3**2),\
                            (q2**4)*(q3**2),(q2**2)*(q3**4),\
                            (q1**4)*(q2**3),(q1**3)*(q2**4),\
                            (q1**3)*(q3**4),(q1**4)*(q3**3),\
                            (q2**4)*(q3**3),(q2**3)*(q3**4),\
                            (q1**4)*(q2**2)*q3,(q1**2)*(q2**4)*q3,\
                            (q1**2)*q2*(q3**4),(q1**4)*q2*(q3**2),\
                            q1*(q2**4)*(q3**2),q1*(q2**2)*(q3**4),\
                            ]))[0,0]
    return(output)
def F9(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=54, "Insufficient number of coefficients for F9 (54)"
    import numpy as np
    import sympy as sp
    output =F7(q1,q2,q3,Coefficients) + \
            (sp.Matrix(Coefficients[42:54]).T*\
                sp.Matrix([ (q1**5),(q2**5),(q3**5),\
                            (q1**2)*(q2**2)*q3,(q1**2)*q2*(q3**2),\
                            q1*(q2**2)*(q3**2),(q1**3)*q2*(q3**2),\
                            (q1**2)*(q2**3)*q3,q1*(q2**3)*(q3**2),\
                            q1*(q2**2)*(q3**3),(q1**3)*(q2**3),\
                            (q1**3)*(q3**3)]))[0,0]
    return(output)
def F10(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=88, "Insufficient number of coefficients for F10 (88)"
    import numpy as np
    import sympy as sp
    output =F7(q1,q2,q3,Coefficients) + \
            (sp.Matrix(Coefficients[42:88]).T*\
                sp.Matrix([ (q1**5),(q2**5),(q3**5),\
                            (q1**3)*(q2**3),(q1**3)*q2*(q3**2),\
                            (q1**2)*(q2**3)*q3,q1*(q2**3)*(q2**2),\
                            (q1**2)*q2*(q3**3),q1*(q2**2)*q3**3,\
                            (q1**2)*(q2**2)*q3,\
                            (q1**2)*q2*(q3**2),\
                            q1*(q2**2)*(q3**2),\
                            (q1**3)*(q2**3)*q3,\
                            (q1**3)*q2*(q3**3),\
                            q1*(q2**3)*(q3**3),\
                            (q1**5)*q2,(q1**5)*q3,\
                            q1*(q2**5),(q2**5)*q3,\
                            q1*(q3**5),q2*(q3**5),\
                            (q1**6)*q2,q1*(q2**6)\
                            (q1**7)*q2,q1*(q2**7)\
                            (q1**8)*q2,q1*(q2**8)\
                            (q1**6)*q3,(q2**6)*q3,\
                            (q1**7)*q3,(q2**7)*q3,\
                            (q1**8)*q3,(q2**8)*q3,\
                            q1*(q3**6),q2*(q3**6),\
                            (q1**4)*(q2**2)*(q3**2),(q1**4)*(q2**3)*(q3**2),\
                            (q1**4)*(q2**2)*(q3**3),(q1**2)*(q2**4)*(q3**2),\
                            (q1**2)*(q2**4)*(q3**3),(q1**3)*(q2**4)*(q3**2),\
                            (q1**2)*(q2**2)*(q3**4),(q1**2)*(q2**3)*(q3**4),\
                            (q1**3)*(q2**2)*(q3**4),\
                            (q1**3)*(q2**3)*(q3**3),\
                            (q1**4)*(q2**4)*(q3**4)]))[0,0]
    return(output)
def F11(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=75, "Insufficient number of coefficients for F11 (75)"
    import numpy as np
    import sympy as sp
    output =F9(q1,q2,q3,Coefficients) + \
            (sp.Matrix(Coefficients[54:75]).T*\
                sp.Matrix([ (q1**3)*(q2**3)*q3,(q1**3)*q2*(q3**3),\
                            q1*(q2**3)*(q3**3),\
                            (q1**5)*q2,(q1**5)*q3,\
                            q1*(q2**5),(q2**5)*q3,\
                            q1*(q3**5),q2*(q3**5),\
                            (q1**6)*q2,q1*(q2**6),\
                            (q1**7)*q2,q1*(q2**7),\
                            (q1**8)*q2,q1*(q2**8),\
                            (q1**6)*q3,(q2**6)*q3,\
                            (q1**7)*q3,(q2**7)*q3,\
                            (q1**8)*q3,(q2**8)*q3,\
                            ]))[0,0]
    return(output)
def F12(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=22, "Insufficient number of coefficients for F12 (22)"
    import numpy as np
    import sympy as sp
    output = (sp.Matrix(Coefficients[:22]).T*\
                sp.Matrix([1,q2,q3,\
                            q2**2,q3**2,\
                            q2**3,q3**3,\
                            q1*q2,q2*q3,\
                            q1*q2*q3,(q1**2)*q2,\
                            q1*(q2**2),(q2**2)*q3,q1*(q3**2),\
                            q2*(q3**2),(q1**3)*q2,(q1**3)*q3,\
                            q1*(q2**3),(q2**3)*q3,q1*(q3**2),\
                            q2*(q3**2),(q1*q2*q3)**2]))[0,0]
    return(output)
def F13(q1,q2,q3,Coefficients):
    assert len(Coefficients)>=30, "Insufficient number of coefficients for F13 (30)"
    import numpy as np
    import sympy as sp
    output = F12(q1,q2,q3,Coefficients)+\
            (sp.Matrix(Coefficients[22:30]).T*\
                sp.Matrix([(q1**4),(q2**4),(q3**4),\
                            (q1*q2)**2,(q2*q3)**2,\
                            (q1**2)*q2*q3,q1*(q2**2)*q3,\
                            q1*q2*(q3**2)]))[0,0]
    return(output)
def F14(q1,q2,q3,q4,Coefficients):
    assert len(Coefficients)>=18, "Insufficient number of coefficients for F14 (18)"
    import numpy as np
    import sympy as sp
    output = (sp.Matrix(Coefficients[:18]).T*\
                sp.Matrix([1,q1,(q1**2),(q1**3),\
                            q2,(q2**3),(q2**4),\
                            q3,(q3**2),(q3**3),\
                            q4,(q4**2),(q4**3),(q4**4),\
                            q1*q2,(q1**2)*q2,\
                            q1*(q2**2),q2*(q3**2)]))[0,0]
    return(output)
def F15(q1,q2,q3,q4,q5,Coefficients):
    assert len(Coefficients)>=29, "Insufficient number of coefficients for F15 (29)"
    import numpy as np
    import sympy as sp
    output = (sp.Matrix(Coefficients[:29]).T*\
                sp.Matrix([1,q1,q2,q3,q4,q5,\
                            (q1**2),(q2**2),(q3**2),(q4**2),(q5**2),\
                            (q1**3),(q2**3),(q3**3),(q4**3),(q5**3),\
                            q1*q2,q1*q3,q2*q3,q4*q5,\
                            q1*q2*q3,(q1**2)*q2,(q1**2)*q3,\
                            q1*(q2**2),(q2**2)*q3,\
                            q1*(q3**2),q2*(q3**2),\
                            (q4**2)*q5,q4*(q5**2)]))[0,0]
    return(output)


import numpy as np
import sympy as sp
q1,q2,q3,q4,q5 = sp.symbols('q1'), sp.symbols('q2'), sp.symbols('q3'), sp.symbols('q4'), sp.symbols('q5')

Triceps = {}
Triceps['Coefficients'] = [-0.0202998882726984,-0.0215606016088360,0.0356484421721072,-0.0178502164049808,0.0026245033894127]
theta = np.arange(0,3*np.pi/4,0.01)
Triceps['MA'] = F1(q1,Triceps['Coefficients'])
[Triceps['MA'].subs(q1,Q) for Q in theta]
