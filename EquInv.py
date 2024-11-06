import json

import matplotlib.pyplot as plt
import subprocess
import os
import shutil

progPath = 'c:/TkachenkoEE/work/equilibrium/Inv_Globus/Globus3/'
resPath = 'c:/TkachenkoEE/work/Data/equlibrium/Inv_db/'

Ncur = 6

with open('%sshotn.txt' %resPath, 'r') as shtFile:
    shotn = float(shtFile.read())

def geomParametrs(R, a, k, d, xdot='down'):
    Rmax = R + a
    Rmin = R - a
    Rx = R - d*a
    if xdot == 'up':
        Zx = a * k
    else:
        Zx = - a * k

    return  Rmax, Rmin, Rx, Zx


def iter4(Ip, Bt, Rmax, Rmin, Rx, Zx, Zc, addDots=[], li=1, betaI=0.5, al1=1.5, al2=1.5, configure='div', Psiak=0, addOpc=0, psito=0.16, fix_curr={}):
    Nopc = 90
    Nt = 2 + len(addDots)
    if configure == 'lim':
        Nopc = 10
        Nt = 4 + len(addDots)
    if addOpc and Psiak:
        Nopc += addOpc

    cur0 = [0]*Ncur
    cur0[0] = 5.5910E+05
    cur0[1] = -8.6690E+04
    cur0[2] = 5.5972E+04
    cur0[3] = 2.9223E+04
    cur0[4] = 6.0795E+06

    if (Rmax + Rmin)/2 - Rx > 0:
        cur0[5] = 0
        fix_curr[6] = 0
    else:
        cur0[5] = -100E+03


    with open('%siter4.dat' %progPath, 'w') as file:
        #lin1
        file.write('%i  ' %Nopc)
        file.write('%i  ' %Nt)
        file.write('!Globus-3 Warm v3, Ip=%.2fMA, Ics(t=0)=40kA \n\n' )

        #line3
        for i in range(Ncur):
            file.write('%e ' %cur0[i])
        file.write('\n\n')

        #line5
        file.write(' %.1f  ' %li)
        file.write('%.3f  ' %betaI)
        file.write('%.4f  ' %Psiak)
        file.write('%.1fe6  ' %Ip)
        file.write('%.1f  ' %al1)
        file.write('%.1f  ' %al2)
        file.write('0')
        file.write('\n\n')

        #lines dot
        file.write(' %.3f    %.3f \n' %(Rmin, Zc))
        file.write(' %.3f    %.3f \n' %(Rmax, Zc))

        for i in range(len(addDots)):
            file.write(' %.3f    %.3f \n' % (addDots[i][0], addDots[i][1]))
        file.write('\n')

        if configure=='div':
            file.write(' %.3f    %.3f \n' % (Rx, Zx))
            file.write('\n')

        #B line
        file.write(' %.2f  ' % psito)
        file.write('%.2f  ' % (Bt*(Rmax+Rmin)/2))
        file.write('%.3f  ' % (abs(Zx)-0.001))
        file.write('\n\n')

        #fix currents line
        curcount = 0
        for key in fix_curr:
            file.write(' %i  %e   ' %(key, fix_curr[key]))
            curcount+=1

        for i in range(Ncur-curcount-1):
            file.write('0  1   ')
        file.write('\n\n')

        #end
        file.write('-1  0   ')
        file.write('\n\n')


def readOut():
    with open('%sout.wr' %progPath, 'r') as file:
        dataRaw = file.read().splitlines()
        data = []
        for line in dataRaw:
            data.append(line.split())
        nr = int(data[0][0])
        nz = int(data[0][1])
        n = int(data[0][2])

        rgr = [float(data[1][i]) for i in range(nr)]
        zgr = [float(data[1][i]) for i in range(nr, nr+nz)]

        ugr = [[float(data[2][j*nr+i]) for i in range(nr)] for j in range(nz)]
        tpl = [[float(data[3][j*nr+i]) for i in range(nr)] for j in range(nz)]

        dots_b = [[float(data[4][i]), float(data[4][i+n])] for i in range(n)]

        Rm = float(data[5][0])
        Zm = float(data[5][1])
        psi_m = float(data[5][2])
        Rx = float(data[5][3])
        Zx = float(data[5][4])
        psi_b = float(data[5][5])

        Ip = float(data[6][0])
        psi_ext = float(data[6][1])
        betap = float(data[6][2])
        li = float(data[6][3])
        Lp = float(data[6][4])
        Rmin = float(data[6][5])
        Rmax = float(data[6][6])
        RZmax = float(data[6][7])
        Zmax = float(data[6][8])
        kx = float(data[6][9])
        dx = float(data[6][10])

        Ncur = int(data[7][0])
        Nopc = float(data[7][1])

        curr = {}
        for i in range(Ncur):
            curr[i+1] = float(data[8][i])

        return nr, nz, n, rgr, zgr, ugr, tpl, dots_b, Rm, Zm, psi_m, Rx, Zx, psi_b, Ip, psi_ext, betap, li, Lp, Rmin, Rmax, RZmax, Zmax, kx, dx, Nopc, Ncur, curr


def readIter4out():
    with open('%siter4_new.out' %progPath, 'r') as file:
        dataRaw = file.read().splitlines()
        kiter_line = 19
        for i, line in enumerate(dataRaw):
            if line.split():
                if line.split()[0] == 'kiter=':
                    kiter_line = i
        date = dataRaw[0]
        coil = dataRaw[1:kiter_line]
        params0 = dataRaw[kiter_line+2].split()
        params1 = dataRaw[kiter_line+6].split()
        alf0 = float(params0[4])
        alf1 = float(params0[5])
        alf2 = float(params0[6])
        betaI = float(params1[5])

    return date, coil, alf0, alf1, alf2, betaI


def plot_equ(rgr, zgr, ugr, psi_b, Rmax0, Rmin0, Zc0, Rx0, Zx0, Rmax, Rmin, Rx, Zx, Rm, Zm, textstr1, textstr2):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 2.6)
    ax.set_ylim(-1.3, 1.3)
    ax.contour(rgr, zgr, ugr, levels=[psi_b], color='g')
    ax.contour(rgr, zgr, ugr, levels=100, alpha=0.5)

    ax.scatter(Rmin0, Zc0, marker='x', zorder=2)
    ax.scatter(Rmax0, Zc0, marker='x', zorder=2)
    ax.scatter(Rx0, Zx0, marker='x', zorder=2)

    ax.scatter(Rmin, Zm, marker='.', zorder=3)
    ax.scatter(Rmax, Zm, marker='.', zorder=3)
    ax.scatter(Rx, Zx, marker='.', zorder=3)

    plt.scatter(Rm, Zm)

    with open(progPath + 'BLANFW.DAT', 'r') as file3:
        blanfw = []
        for line in file3:
            blanfw.append(line.split())

    with open(progPath + 'LIMPNT.dat', 'r') as file4:
        limpnt = []
        for line in file4:
            limpnt.append(line.split())

    # limitr

    nlim = int(limpnt[0][0])
    #print(nlim)
    nvv = int(blanfw[1][0])
    ax.plot([float(blanfw[i + 2][1]) for i in range(nvv)], [float(blanfw[i + 2][2]) for i in range(nvv)], 'black')
    ax.plot([float(blanfw[i + 2][3]) for i in range(nvv)], [float(blanfw[i + 2][4]) for i in range(nvv)], 'black')

    ax.plot([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])],
            [float(limpnt[i + 1][1]) for i in range(nlim)] + [float(limpnt[1][1])], 'm')
    #print("!!!!!!!!!!!!!!!!!")
    #print(min([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])]))
    #print(max([float(limpnt[i + 1][0]) for i in range(nlim)] + [float(limpnt[1][0])]))
    # ax.grid()
    ax.set_ylabel('z, м')
    ax.set_xlabel('r, м')
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    ax.text(0.55, 0.8, textstr1, transform=ax.transAxes, fontsize=14, verticalalignment='top')
    ax.text(0.55, 0.4, textstr2, transform=ax.transAxes, fontsize=14, verticalalignment='top')

    plt.savefig('%s%i.png' %(resPath, shotn))
    #ax.legend()


def runCalc():
    try:
        process = subprocess.Popen(["%siter1_exe.bat" %progPath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(timeout=30)
        print(stderr)
        shutil.copyfile('%siter4.out' % progPath, '%siter4_new.out' % progPath)
        os.remove('%siter4.out' % progPath)
        return 0
    except subprocess.TimeoutExpired:
        print('time over')
        subprocess.check_call("TASKKILL /F /PID {pid} /T".format(pid=process.pid))
        return -1

def end():
    with open('%sshotn.txt' % resPath, 'w') as shtFile:
        shtFile.write('%i' %(shotn+1))

def fileForPet(curr, betaI, alf0, alf1, alf2, Ip, Bt, R, Rx, Zx, Zmax, RZmax, Rm, Zm):
    data = {}
    data['COIL'] = curr
    data['DURS'] = {'betaI': betaI,
                    'alf0': alf0,
                    'alf1': alf1,
                    'alf2':alf2,
                    'Ip': Ip*1e-6}
    data['DINA'] = {'R': R, 'Bt': Bt}
    if Zx > 0:
        Rx_up = Rx
        Zx_up = Zx
        Rx_d = Rx
        Zx_d = -Zx
    else:
        Rx_up = RZmax
        Zx_up = Zmax
        Rx_d = Rx
        Zx_d = Zx
    data['DATA'] = {'Rx_up': Rx_up, 'Zx_up': Zx_up, 'Rx_down': Rx_d, 'Zx_down': Zx_d,
                    'Rm': Rm, 'Zm': Zm}
    with open('%s%i.json' %(resPath, shotn), 'w') as jsFile:
        json.dump(data, jsFile)

def shotnDataFileWrite(R, a, k, d, Ip, Bt):
    with open('%sshotnData.txt' %resPath, 'a') as file:
        '''file.write('%8s' % 'shotn')
        file.write('%8s' % 'R, m')
        file.write('%8s' % 'a, m')
        file.write('%8s' % 'elong')
        file.write('%8s' % 'triang')
        file.write('%8s' % 'Ip, MA')
        file.write('%8s' % 'Bt, T')
        file.write('\n')'''
        file.write('%8i' %shotn)
        file.write('%8.2f' %R)
        file.write('%8.2f' %a)
        file.write('%8.2f' %k)
        file.write('%8.2f' %d)
        file.write('%8.2f' %Ip)
        file.write('%8.2f' %Bt)
        file.write('\n')









