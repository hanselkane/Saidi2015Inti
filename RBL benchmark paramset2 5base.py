import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.constants as const
import math
#-------Cluster decay investigation within a modified Woods-Saxon potential------
#----------------------------------------------------------------------------------------JANGAN LUPA COMMENTTT
print("berikut adalah data - datanya: ")
print("Ra(A=222; Z=86) --> Pb(A=208; Z=82)  +  C(A=14; Z=6) ")

#constant
hbareV = 6.582119569*(10**-16)
heV = 4.135667696*(10**-15)
c = 2.997924588*(10**8)

Ndata = 3000
base = 0
#Bikin array untuk plot data
x = np.empty(Ndata+1)
g=0
x[g] = 4.5
while(g<Ndata):
    x[g+1] = x[g]+0.008
    g = g+1

result = np.empty(23)
arrayerror = np.empty(23)

#-----------------------------Menghitung potensial-----------------------------              
def Ii(A, Z):
    return(((A-Z)-Z)/A)

def a(Ad, Zd):
    return(0.2+(-0.254*Ii(Ad, Zd)))

def Vc(Zc, Zd, r):
    #diubah ke femtometer, ke MeV
    #return((((Zc*Zd*2.56698*9)*(10**1))/r)*6.242*(10^12))
    return(1.44*Zc*Zd/r)

def Vl(l, r, Ad, Ac):
    return(((const.hbar**2)*l*(l+1))/(2*miu(Ad, Ac)*r))

def miu(Ad, Ac):
    return((Ad*Ac)/(Ad+Ac))

def Vo(Ad, Zd, Ac, Zc):
    adad = ((Ad**(1/3))*(Ac**(1/3)))/((Ad**(1/3))+(Ac**(1/3)))
    return(-45.16*(1-(-0.11*(Ii(Ad, Zd)-Ii(Ac, Zc))))*adad)

def Ri(Ai):
    #base99 min 1.25
    return(1.17*(Ai**(1/3)))

def Ro(Ac, Ad):
    return(Ri(Ac)+Ri(Ad)-1.37)

def Vn(Ad, Zd, Ac, Zc, r):
    asdf = (r-Ro(Ac, Ad))/a(Ad, Zd)
    return(Vo(Ad, Zd, Ac, Zc)/(1+np.exp(asdf)))

def Vtotal(Ad, Zd, Ac, Zc, r, l):
    return(Vn(Ad, Zd, Ac, Zc, r)+Vc(Zc, Zd, r)+Vl(l, r, Ad, Ac))

#-----------------------------Menghitung potensial-----------------------------

#------------------------Menghitung Half-life--------------------------
#P nilainya terlalu kecil, digunakan sifat logaritma
#Tahap : hitung ln P, hitung pers:: ln lambda = ln v + ln P
#Didapat ln lambda, hitung pers:: ln lambda + ln T = ln 2
#Didapat ln T. Konversi ln T

def Ev(Ap, Zp, Ac, Zc, Ad, Zd):  
    return(Qtheory*(0.056+(0.039*math.exp((4-Ad)/2.5))))

def v(Ap, Zp, Ac, Zc, Ad, Zd):
    return(2*Ev(Ap, Zp, Ac, Zc, Ad, Zd)*(10**6)/heV)

def integrating(Ad, Zd, Ac, Zc, l, Rin, Rout):
    #Qtheory harusnya Qresult
    return(integrate.quad(lambda r: math.sqrt(Vtotal(Ad, Zd, Ac, Zc, r, l)-Qtheory), Rin, Rout))
           
    #Digunakan sifat logaritma
def lnP(Ad, Zd, Ac, Zc, l, Rin, Rout):
    #unitedconst adalah konstanta gabungan di dalam exp P, ini terdiri dari
    #-2 ; 1/hbareV ; konversi hbareV ke Mev; 1/c(light speed) ;
    #konversi c(m/s) ke (fm/s); akar(931.5 MeV*2*miu); dr(m) ke (fm)
    unitedconst = math.sqrt(2*931.5*miu(Ad,Ac))*(-2/(6.58212*2.998))*(10**-1)
    return(unitedconst*integrating(Ad, Zd, Ac, Zc, l, Rin, Rout)[0])

def lnlambda(Ap, Zp, Ad, Zd, Ac, Zc, l, Rin, Rout):
    return(math.log(v(Ap, Zp, Ac, Zc, Ad, Zd))+lnP(Ad, Zd, Ac, Zc, l, Rin, Rout))

def lnT(Ap, Zp, Ad, Zd, Ac, Zc, l, Rin, Rout):
    return(math.log(0.693)-lnlambda(Ap, Zp, Ad, Zd, Ac, Zc, l , Rin, Rout))
#------------------------Menghitung Half-life--------------------------

def error(base, result):
    if(base==0):
        kunjaw=11.14
    if(base==1):
        kunjaw=16.333
    if(base==2):
        kunjaw=21.821
    if(base==3):
        kunjaw=18.123
    if(base==4):
        kunjaw=21.122
    if(base==5):
        kunjaw=20.095
    if(base==6):
        kunjaw=21.068
    if(base==7):
        kunjaw=19.147
    if(base==8):
        kunjaw=25.08
    if(base==9):
        kunjaw=30.398
    if(base==10):
        kunjaw=24.792
    if(base==11):
        kunjaw=29.854
    if(base==12):
        kunjaw=23.864
    if(base==13):
        kunjaw=29.933
    if(base==14):
        kunjaw=28.034
    if(base==15):
        kunjaw=23.918
    if(base==16):
        kunjaw=23.943
    if(base==17):
        kunjaw=29.79
    if(base==18):
        kunjaw=27.44
    if(base==19):
        kunjaw=18.862
    if(base==20):
        kunjaw=24.541
    if(base==21):
        kunjaw=23.664
    return((result-kunjaw)/kunjaw)*100

def average(values):
    sum = 0
    i=0
    #Tambahinnnn
    while(i<22):
        sum = sum+(abs(values[i]))
        i=i+1
    return(sum/(i+1))    

#--------------------------MAIN----------------------
#Tambahinnn
while(base<22):
    if(base==0):
        Ap = 222; Zp = 88; Ac = 14; Zc = 6;
        Ad = 208; Zd = 82; Qtheory = 33.049;
    if(base==1):
        Ap = 224; Zp = 88; Ac = 14; Zc = 6
        Ad = 210; Zd = 82; Qtheory = 30.535;
    if(base==2):
        Ap = 226; Zp = 88; Ac = 14; Zc = 6
        Ad = 212; Zd = 82; Qtheory = 28.196;
    if(base==3):
        Ap = 226; Zp = 90; Ac = 18; Zc = 8;
        Ad = 208; Zd = 82; Qtheory = 45.726;
    if(base==4):
        Ap = 228; Zp = 90; Ac = 20; Zc = 8
        Ad = 208; Zd = 82; Qtheory = 44.722;
    if(base==5):
        Ap = 230; Zp = 92; Ac = 22; Zc = 10
        Ad = 208; Zd = 82; Qtheory = 61.387;
    if(base==6):
        Ap = 230; Zp = 92; Ac = 24; Zc = 10;
        Ad = 206; Zd = 82; Qtheory = 61.35;
    if(base==7):
        Ap = 232; Zp = 92; Ac = 24; Zc = 10
        Ad = 208; Zd = 82; Qtheory = 62.309;
    if(base==8):
        Ap = 234; Zp = 92; Ac = 24; Zc = 10
        Ad = 210; Zd = 82; Qtheory = 58.825;
    if(base==9):
        Ap = 236; Zp = 92; Ac = 24; Zc = 10
        Ad = 212; Zd = 82; Qtheory = 55.944;
    if(base==10):
        Ap = 234; Zp = 92; Ac = 26; Zc = 10
        Ad = 208; Zd = 82; Qtheory = 59.464;
    if(base==11):
        Ap = 236; Zp = 92; Ac = 26; Zc = 10
        Ad = 210; Zd = 82; Qtheory = 56.744;
    if(base==12):
        Ap = 230; Zp = 90; Ac = 24; Zc = 10
        Ad = 206; Zd = 80; Qtheory = 57.761;
    if(base==13):
        Ap = 232; Zp = 90; Ac = 24; Zc = 10
        Ad = 208; Zd = 80; Qtheory = 54.509;
    if(base==14):
        Ap = 232; Zp = 90; Ac = 26; Zc = 10
        Ad = 206; Zd = 80; Qtheory = 55.964;
    if(base==15):
        Ap = 232; Zp = 92; Ac = 28; Zc = 12
        Ad = 204; Zd = 80; Qtheory = 74.318;
    if(base==16):
        Ap = 234; Zp = 92; Ac = 28; Zc = 12
        Ad = 206; Zd = 80; Qtheory = 74.11;
    if(base==17):
        Ap = 236; Zp = 92; Ac = 28; Zc = 12
        Ad = 208; Zd = 80; Qtheory = 70.564;
    if(base==18):
        Ap = 236; Zp = 92; Ac = 30; Zc = 12
        Ad = 206; Zd = 80; Qtheory = 72.303;
    if(base==19):
        Ap = 236; Zp = 94; Ac = 28; Zc = 12
        Ad = 208; Zd = 82; Qtheory = 79.669;
    if(base==20):
        Ap = 238; Zp = 94; Ac = 28; Zc = 12
        Ad = 210; Zd = 82; Qtheory = 75.911;
    if(base==21):
        Ap = 238; Zp = 94; Ac = 30; Zc = 12
        Ad = 208; Zd = 82; Qtheory = 76.823;
    y = np.empty(Ndata+1)
    i = 0
    Mroot = np.empty(3)  #Matrix untuk titik r di potensial yg = Q value
    iroot=-1
    locked=0
        
    #Mengisi matrix y(Vtotal), dan jika ad yg nilainya dekat dngn Q, maka dicatat di
    #dicatat di Mroot
    while(i<Ndata+1):
        y[i]=Vtotal(Ad, Zd, Ac, Zc, (x[i]), 0)
        if (y[i]-Qtheory<= 0.7)&(y[i]-Qtheory>=0):
            if(locked==0):
                locked=1
                iroot=iroot+1
                Mroot[iroot] = x[i]
        if (y[i]-Qtheory<0)or((x[i]>Mroot[1]+2)&(x[i]<(Mroot[1]+5))):
            locked=0
        i = i+1
        
    Rin=Mroot[1]
    Rout=Mroot[2]
    result[base]=lnT(Ap, Zp, Ad, Zd, Ac, Zc, 0, Rin, Rout)*(math.log(math.e,10))
    arrayerror[base]=error(base, result[base])
    base=base+1

i=0
while(i<22):
    print('Error',i,' = ',arrayerror[i],'; ')
    i=i+1
print("Error average = ", average(arrayerror))
#--------------------------MAIN----------------------  


#Dibawah ini khusus untuk plot saja-------------------------------
base=1
if(base==0):
    Ap = 222; Zp = 88; Ac = 14; Zc = 6;
    Ad = 208; Zd = 82; Qtheory = 33.049;
if(base==1):
    #base 26206
    Ap = 224; Zp = 88; Ac = 14; Zc = 6
    Ad = 210; Zd = 82; Qtheory = 30.535;
if(base==2):
    #base 99
    Ap = 226; Zp = 88; Ac = 14; Zc = 6
    Ad = 212; Zd = 82; Qtheory = 28.196;
yVc = np.empty(Ndata+1)
yVo = np.empty(Ndata+1)
yVn = np.empty(Ndata+1)

#Mengisi matrix Qtheory(qq) untuk membentuk garis lurus
qq = np.empty(Ndata+1)
i = 0
while(i<Ndata+1):
    qq[i]=Qtheory
    i = i+1

i = 0
while(i<Ndata+1):
    y[i]=Vtotal(Ad, Zd, Ac, Zc, (x[i]), 0)
    i=i+1
    
#Mengisi matrix yVc
i = 0
while(i<Ndata+1):
    yVc[i]=Vc(Zc, Zd, x[i])
    i=i+1

#Mengisi matrix yVo
i = 0
while(i<Ndata+1):
    yVo[i]=Vo(Ad, Zd, Ac, Zc)
    i=i+1

#Mengisi matrix yVn
i = 0
while(i<Ndata+1):
    yVn[i]=Vn(Ad, Zd, Ac, Zc, x[i])
    i=i+1

plt.xlabel('r(fm)')
plt.ylabel('V(MeV)')
plt.xlim([0, 35])
plt.ylim([15, 65])
plt.title('Ra(A=222; Z=86) --> Pb(A=208; Z=82)  +  C(A=14; Z=6)')
plt.plot(x,qq,label='Q value')
plt.plot(x,y,label='Interaction Potential')
plt.plot(x,yVo,label='Vo constant')
plt.plot(x,yVc,label='Vc Coulomb Potential')
plt.plot(x,yVn,label='Vn Potential')
plt.legend()
plt.ion()
plt.show()
plt.pause(0.001)

print(Mroot)








