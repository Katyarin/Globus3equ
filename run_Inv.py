import EquInv
import matplotlib.pyplot as plt

progPath = 'c:/TkachenkoEE/work/equilibrium/Inv_Globus/Globus3/'

Ip0 = 0.8
Bt = 1.5

'''Rmax0 = 1.215
Rmin0 = 0.335
Rx0 = 0.67
Zx0 = -0.88'''
Zc0 = 0

Rmax0, Rmin0, Rx0, Zx0 = EquInv.geomParametrs(R=0.78, a=0.44, k=1.9, d=0.1)
print(Rmax0, Rmin0, Rx0, Zx0)
#Rmax0, Rmin0, Rx0, Zx0 = 1.215, 0.335, 0.938, -0.825

li0 = 1
betaI0 = 0.434
al1=1.45
al2=1.5
configure='div'

Psiak = -1.0481

addOpc = 0
psito = 0.16

fix_curr = {6: -1e5}
addDots = []
EquInv.iter4(Ip=Ip0, Bt=Bt, Rmax=Rmax0, Rmin=Rmin0, Rx=Rx0, Zx=Zx0, Zc=Zc0, addDots=addDots,
 li=li0, betaI=betaI0, al1=al1, al2=al2, configure=configure, Psiak=Psiak, addOpc=addOpc, psito=psito, fix_curr=fix_curr)

result = EquInv.runCalc()

if result == 0:
    nr, nz, n, rgr, zgr, ugr, tpl, dots_b, Rm, Zm, psi_m, Rx, Zx, psi_b, Ip, psi_ext, betap, li, Lp, Rmin, Rmax, RZmax, Zmax, kx, dx, Nopc, Ncur, curr = EquInv.readOut()
    date, coil, alf0, alf1, alf2, betaI = EquInv.readIter4out()

    print(alf0, alf1, alf2, betaI)
    print(kx, dx)
    R = (Rmax+Rmin)/2
    a = (Rmax-Rmin)/2

    textstr1 = '\n'.join((
        r'Parameters:',
        r'$I_p=%.2f$ MA' % (Ip/1e6,),
        r'$B_T=%.2f$ T' % (Bt,),
        r'$R=%.2f$ m' % (R,),
        r'$a=%.2f$ m' % (a,),
        r'$\kappa=%.2f$' % (kx,),
        r'$\delta=%.2f$' % (dx,),
        r'$\beta_I=%.2f$' % (betaI,),
        r'$l_i=%.2f$' % (li,),
    ))
    curr_tuple = tuple([r'$I_{%i}=%e$ A' %(i, curr[i]) for i in range(1, Ncur+1)])
    textstr2 = '\n'.join(('Currents:', ) + curr_tuple)

    EquInv.plot_equ(rgr, zgr, ugr, psi_b, Rmax0, Rmin0, Zc0, Rx0, Zx0, Rmax, Rmin, Rx, Zx, Rm, Zm, textstr1, textstr2)

    EquInv.shotnDataFileWrite(R, a, kx, dx, Ip/1e6, Bt)
    EquInv.fileForPet(curr, betaI, alf0, alf1, alf2, Ip, Bt, R, Rx, Zx, Zmax, RZmax, Rm, Zm)

    EquInv.end()
    plt.show()



