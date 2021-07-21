import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
from matplotlib.animation import FuncAnimation
from celluloid import Camera
import math
from mpmath import mpf

#________________________________________________
#==============ADJUSTABLE PARAMETER==============

element="Xe" #Tersedia: He, Ne, Ar, Kr, Xe
temperatur=1000 #Kelvin 
numpar=50 #Banyak particles #150
Xl=0 #X low boundary for simulated box
Xh=150 #X high boundary for simulated box
Yl=0 #Y low boundary for simulated box
Yh=150 #Y low boundary for simulated box
durasi = 0.001 #0.03
dt = 0.0001 #timestep

#==============ADJUSTABLE PARAMETER==============
#------------------------------------------------



#Initial data of particles and simulation
if (element=="He"):
    m = mpf(6.64*(10**(-27)))
    r = 1.28
    epsilon = (10.22*1.38*(10**(-23)))
    ballcolor ='w'
if (element=="Ne"):
    m = mpf(3.35*(10**(-26)))
    r = 1.375
    epsilon = (35.6*1.38*(10**(-23)))
    ballcolor ='r'
if (element=="Ar"):
    m = mpf(6.64*(10**(-26)))
    r = 1.7
    epsilon = (120*1.38*(10**(-23)))
    ballcolor ='b'
if (element=="Kr"):
    m = mpf(1.39*(10**(-25)))
    r = 1.8
    epsilon = (171*1.38*(10**(-23)))
    ballcolor ='w'
if (element=="Xe"):
    m = mpf(2.18*(10**(-25)))
    r = 2.05
    epsilon = (220*1.38*(10**(-23)))
    ballcolor ='purple'

t = np.arange(0,durasi,dt) #simulated time
e = 1 #elasticity

Rcutoff=1000 #Rcutoff for Lennard-Jones calculation
numvar=11 #number of variables.

#Vawal dihitung bergantung pada temperatur
Kboltzman=mpf(1.38*(10**-23))
Vawal = np.sqrt((3*Kboltzman*temperatur)/m) #maximum speed

#Initialize the array of datas.
#3d array with time, number of particle and tracked variable.
Whole = np.zeros((len(t),numpar,numvar))
Em = np.zeros(len(t)-1)
Epotensial = np.zeros(len(t)-1)
#at numvar data is as follows
#numvar 0 = X position
#numvar 1 = Y position
#numvar 2 = Vx
#numvar 3 = Vy
#numvar 4 = Ek
#numvar 5 = Fx
#numvar 6 = Fy
#numvar 7 = Ax (acceleration)
#numvar 8 = Ay
#numvar 9 = Epotensial
#numvar 10 = Etotal
collide = False #check if collision calculation had been done
Distnow=0
Distbef=0
theta = 0.0
n = 0
closeparticles = 1
searching = True

for n in range(0,numpar):
    #making sure there aren't any clippings of particles at start
    while searching:
        clip=False
        Whole[0][n][0]= random.uniform((Xl+2*r), (Xh-2*r))
        Whole[0][n][1]= random.uniform((Yl+2*r), (Yh-2*r))
        for b in range(0,n):
            Distnow = np.sqrt((Whole[0][b][1]-Whole[0][n][1])**2+(Whole[0][b][0]-Whole[0][n][0])**2) 
            if (b!=n) and (Distnow<=2.2*r):
                clip=True
        if clip:
            searching=True
        elif not clip:
            break
    Whole[0][n][2]= random.uniform(-Vawal, Vawal)
    if (random.random()>0.5):
        tanda=1
    else:
        tanda=-1
    Whole[0][n][3]= tanda*np.sqrt(Vawal**2-Whole[0][n][2]**2)
    Whole[0][n][4]=0.5*m*(Whole[0][n][2]*Whole[0][n][2]+Whole[0][n][3]*Whole[0][n][3])
    Whole[0][n][5]= 0
    Whole[0][n][6]= 0
    Whole[0][n][7]= 0
    Whole[0][n][8]= 0
    print(n)
print ('Particles initialized!')

progres=0
for n in range(0,len(t)-1):
    progres=progres+1
    if(progres>(len(t)/8)):
       progres=0
       print("Perhitungan: ",int((n/(len(t)))*100),"%")
    for b in range (0,numpar):
        if Yh-r > (Whole[n][b][1] + Whole[n+1][b][3]*(t[n+1] - t[n])) > Yl+r:
            collide=False
            for p in range (0,numpar):
                Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)

                # Recently added program
                if (p!=b) and (Distnow < Rcutoff):
                    theta = (math.asin((Whole[n][p][1]-Whole[n][b][1])/Distnow))
                    Whole[n][b][6] = Whole[n][b][6] + ((24*(epsilon/(2*r))*(((-2*((2*r)/Distnow)**13))+((((2*r)/Distnow)**7))))*(math.sin(theta)))
                    Whole[n][b][9] = Whole[n][b][9] + ((4*epsilon)*((((2*r)/Distnow)**2)-(((2*r)/Distnow)**6)))
                # Recently added program
               
                if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                    collide=True

            # Recently added program
            Whole[n][b][8] = (Whole[n][b][6] / m)
            Whole[n+1][b][3] = Whole[n][b][3] + (Whole[n][b][8]*(t[n+1] - t[n]))
            Whole[n+1][b][1] = Whole[n][b][1] + Whole[n+1][b][3]*(t[n+1] - t[n])
            # Recently added program
            if Xh-r > (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])) > Xl+r:
                collide=False
                for p in range (0,numpar):
                    Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                    Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)

                    # Recently added program
                    if (p!=b) and (Distnow < Rcutoff):
                        theta = (math.acos((Whole[n][p][0]-Whole[n][b][0])/Distnow))
                        Whole[n][b][5] = Whole[n][b][5] + ((24*(epsilon/(2*r))*(((-2*((2*r)/Distnow)**13))+((((2*r)/Distnow)**7))))*(math.cos(theta)))
                    # Recently added program
                    
                    if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                        collide=True

                # Recently added program
                Whole[n][b][7] = (Whole[n][b][5] / m)
                Whole[n+1][b][2] = Whole[n][b][2] + (Whole[n][b][7]*(t[n+1] - t[n]))
                Whole[n+1][b][0] = Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])
                # Recently added program

            elif ((Xh-r < (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n]))) and (Xh-Whole[n][b][0]>Xh-Whole[n-1][b][0])) or ((((Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])) < Xl+r)) and (Whole[n][b][0]-Xl>Whole[n-1][b][0]-Xl)):
                collide=False
                for p in range (0,numpar):
                    Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                    Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)
                    if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                        collide=True
                if not collide :
                    Whole[n+1][b][2] = Whole[n][b][2]
                Whole[n+1][b][0] = Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])
            else:
                Whole[n+1][b][2] = Whole[n+1][b][2]-e*Whole[n][b][2]
                Whole[n+1][b][0] = (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n]))
        elif ((Yh-r < (Whole[n][b][1] + Whole[n+1][b][3]*(t[n+1] - t[n]))) and Yh-Whole[n][b][1]>Yh-Whole[n-1][b][1]) or (((Whole[n][b][1] + Whole[n+1][b][3]*(t[n+1] - t[n])) < Yl+r) and (Whole[n][b][1]-Yl>Whole[n-1][b][1]-Yl)):
            collide=False
            for p in range (0,numpar):
                Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)
                if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                    collide=True
            if not collide :
                Whole[n+1][b][3] = Whole[n][b][3]
            Whole[n+1][b][1] = Whole[n][b][1] + Whole[n+1][b][3]*(t[n+1] - t[n])
            if Xh-r > (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])) > Xl+r:
                collide=False
                for p in range (0,numpar):
                    Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                    Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)
                    if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                        collide=True
                if not collide :
                    Whole[n+1][b][2] = Whole[n][b][2]
                Whole[n+1][b][0] = Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])
            elif ((Xh-r < (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n]))) and (Xh-Whole[n][b][0]>Xh-Whole[n-1][b][0])) or ((((Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])) < Xl+r)) and (Whole[n][b][0]-Xl>Whole[n-1][b][0]-Xl)):
                collide=False
                for p in range (0,numpar):
                    Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                    Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)
                    if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                        collide=True
                if not collide :
                    Whole[n+1][b][2] = Whole[n][b][2]
                Whole[n+1][b][0] = Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])
            else:
                Whole[n+1][b][2] = Whole[n+1][b][2]-e*Whole[n][b][2]
                Whole[n+1][b][0] = (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n]))
        else:
            Whole[n+1][b][3] = Whole[n+1][b][3]-e*(Whole[n][b][3] )
            Whole[n+1][b][1]= (Whole[n][b][1] + Whole[n+1][b][3]*(t[n+1] - t[n]))
            if Xh-r > (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])) > Xl+r:
                collide=False
                for p in range (0,numpar):
                    Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                    Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)
                    if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                        collide=True
                if not collide :
                    Whole[n+1][b][2] = Whole[n][b][2]
                Whole[n+1][b][0] = Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])
            elif ((Xh-r < (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n]))) and (Xh-Whole[n][b][0]>Xh-Whole[n-1][b][0])) or ((((Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])) < Xl+r)) and (Whole[n][b][0]-Xl>Whole[n-1][b][0]-Xl)):
                collide=False
                for p in range (0,numpar):
                    Distnow = np.sqrt((Whole[n][b][1]-Whole[n][p][1])**2+(Whole[n][b][0]-Whole[n][p][0])**2)
                    Distbef = np.sqrt((Whole[n-1][b][1]-Whole[n-1][p][1])**2+(Whole[n-1][b][0]-Whole[n-1][p][0])**2)
                    if (p!=b) and ((Distnow < 2*r and Distbef>Distnow)):
                        collide=True
                if not collide :
                    Whole[n+1][b][2] = Whole[n][b][2]
                Whole[n+1][b][0] = Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n])
            else:
                Whole[n+1][b][2] = Whole[n+1][b][2]-e*Whole[n][b][2]
                Whole[n+1][b][0] = (Whole[n][b][0] + Whole[n+1][b][2]*(t[n+1] - t[n]))
        Whole[n][b][4]=0.5*m*(Whole[n][b][2]*Whole[n][b][2]+Whole[n][b][3]*Whole[n][b][3])
        Emaxis= np.arange(0,(durasi-dt),dt)
        Em[n]=Em[n]+Whole[n][b][4]        
##Whole[-1][b][4]=0.5*m*(Whole[-1][b][2]*Whole[-1][b][2]+Whole[-1][b][3]*Whole[-1][b][3])
##Em[-1]=Em[-1]+Whole[-1][b][4]
fig = plt.figure()
plt.xlim(Xl, Xh)
plt.ylim(Yl, Yh)
camera = Camera(fig)
previoussecs = int(5/dt)
snaptimer = 0
progres = 0
time=0
for i in range(0,len(t)):
    progres=progres+1
    if(progres>(len(t)/8)):
       progres=0
       print("Animasi: ",int((i/(len(t))*100)),"%")
    snaptimer = snaptimer+1
    rect=patches.Rectangle((0,0),Xh,Yh, ec='r')
    plt.gcf().gca().add_artist(rect)
    for b in range (0,numpar):
        draw_circle = plt.Circle((Whole[i][b][0],Whole[i][b][1]), radius=r, color=ballcolor)
        plt.gcf().gca().add_artist(draw_circle)    
    texty=plt.text(105,5,"Waktu: "+("%.5f" %time)+"dtk",color='b',backgroundcolor='w')
    xlabely=plt.xlabel("Lebar Kotak (Angstrom)")
    ylabely=plt.ylabel("Panjang Kotak (Angstrom)")
    plt.gcf().gca().add_artist(texty)
    time=time+dt
    camera.snap()
animation = camera.animate()
animation.save(("Animasi "+element+", Temp"+str(temperatur)+"K, "+"partikel"+str(numpar)+"buah, "+"durasi"+str(durasi)+"dtk, "+"step"+str(dt)+"dtk"+'.gif'), writer='PillowWriter', fps=20)

plt.close(fig)
plt.plot (Emaxis,Em)
plt.plot(Emaxis,Epotensial,'r')
plt.ylabel('Energi Total(Joule)')
plt.xlabel('Waktu(dtk)')
plt.savefig(("Energi Total "+element+", Temp"+str(temperatur)+"K, "+"partikel"+str(numpar)+"buah, "+"durasi"+str(durasi)+"dtk, "+"step"+str(dt)+"dtk"+'.png'))

allspeed=np.zeros(numpar)
totalspeed=0
for p in range(0,numpar):
    speed = np.sqrt((Whole[len(t)-1][p][2]**2)+(Whole[len(t)-1][p][3]**2))
    totalspeed+=speed
    allspeed[p]=speed
  
plt.figure()
plt.hist(allspeed,bins=12)
plt.title("Grafik distribusi kecepatan (Kec Rata2: "+"%.2f" %(totalspeed/numpar)+")")
plt.xlabel("Kecepatan(m/s)")
plt.ylabel("Banyak partikel(buah)")
plt.savefig(("Distribusi kecepatan "+element+", Temp"+str(temperatur)+"K, "+"partikel"+str(numpar)+"buah, "+"durasi"+str(durasi)+"dtk, "+"step"+str(dt)+"dtk"+'.jpg'))

plt.show()
