import math
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import numpy as np

# calculation of the shadow cast by a wall
# context :
#   an existing wall is raised : what is the impact on the shadow ?
# algorithm from 
# http://heliodon.net/downloads/Beckers_2010_Helio_006_fr_2.pdf

def solve_Kepler(M,e):
  tol = 1.e-09
  breakflag = 0
  E1 = M
  while breakflag == 0:
    E=M+e*math.sin(E1)
    if abs(E-E1)<tol :
      breakflag=1
    E1=E
  while E > 2*math.pi:
    E=E-2*math.pi
  while E < 0:
    E=E+2*math.pi
  return E

def heliodon(Y,Mois,D,hloc):

  UT = hloc-fuseau
  d=367*Y-math.floor((7*(Y+math.floor((Mois+9)/12)))/4)+math.floor((275*Mois)/9)+D-730530
  w=282.9404 + ((4.70935)*10**(-5))*d
  e= 0.016709 - 1.151*10**(-9)*d
  M= 356.0470 + 0.9856002585 *d - math.floor((356.0470 + 0.9856002585 *d)/360)*360
  L= w+M - math.floor((w+M)/360)*360
  obl = 23.4393 - 3.563 * 10**(-7)*d
  M_rad= M *math.pi / 180
  E_rad = solve_Kepler(M_rad,e)
  E = E_rad * 180 / math.pi
  cov= (math.cos(E_rad)-e)/(1-e*math.cos(E_rad))
  siv = math.sqrt(1-e**2)*math.sin(E_rad)/(1-e*math.cos(E_rad))
  xv =  cov*(1-e*math.cos(E_rad))
  yv =  siv*(1-e*math.cos(E_rad))
  r=math.sqrt(xv*xv+yv*yv)
  v=math.atan2(siv,cov)
  v_deg = v*180/math.pi
  #print(d,L,M,E,v_deg)
  lon = v_deg + w - math.floor((v_deg+w)/360)*360
  x = r*math.cos(lon*math.pi/180)
  y = r*math.sin(lon*math.pi/180)
  z=0
  xeq = x
  yeq = y*math.cos(obl*math.pi/180)-z*math.sin(obl*math.pi/180)
  zeq = y*math.sin(obl*math.pi/180)+z*math.cos(obl*math.pi/180)
  #print(xeq,yeq,zeq)
  phi_heures = (math.atan2(yeq,xeq))*12/math.pi
  #print(24+phi_heures)
  alp_deg = (math.atan2(zeq,math.sqrt(xeq**2+yeq**2)))*180/math.pi
  GMSTO = (L+180) / 15
  SIDTIME = GMSTO + UT + LON/15
  #print(GMSTO,UT,LON/15,SIDTIME)
  ph_HA = (SIDTIME-phi_heures)*15
  x_sol = math.cos(ph_HA*math.pi/180) * math.cos(alp_deg*math.pi/180)
  y_sol = math.sin(ph_HA*math.pi/180) * math.cos(alp_deg*math.pi/180)
  z_sol = math.sin(alp_deg*math.pi/180)
  xL = x_sol * math.sin(lat*math.pi/180) - z_sol*math.cos(lat*math.pi/180)
  yL= y_sol
  zL = x_sol * math.cos(lat*math.pi/180) + z_sol*math.sin(lat*math.pi/180)
  #print(xL,yL,zL)
  azimut= math.atan2(yL,xL)*180/math.pi+180
  hauteur_deg = math.asin(zL)*180/math.pi
  #print(azimut-180,hauteur_deg)
  return (azimut-180,hauteur_deg)

# ------------------------
# BEGIN DATA
# ------------------------

LON=4.3
lat = 50.8
fuseau=1

azimuthWall  =0
heightWallAfter = 4.06
longueurMur = 9.1

heightWallBefore=2

d0 = datetime(2022, 1, 1)
d1 = datetime(2023, 1, 1)

# ------------------------
# END DATA
# ------------------------

d=d0

tab_d = []
tab_s = []
tab_sb = []
tab_sa = []
tab_m = []

tab_h= []
tab_t0 = []
tab_t1 = []

average=0
average_b=0
average_a=0
nbminute=0

while d<d1:
  Y=d.year
  Mois = d.month
  Day=d.day
  hloc=d.hour+d.minute/60
  t=heliodon(Y,Mois,Day,hloc)
  tab_h.append(d.strftime("%Y-%m-%d %H:%M"))
  tab_t0.append(t[0])
  tab_t1.append(t[1])
  if t[0] > azimuthWall and t[1]>0 :
    nbminute +=1
    surfaceShadow = heightWallAfter / math.tan(t[1]*math.pi/180) * longueurMur
    surfaceShadowBefore = heightWallBefore / math.tan(t[1]*math.pi/180) * longueurMur
    average_b += min(surfaceShadowBefore,500)
    average_a += min(surfaceShadow,500)
    average += min(surfaceShadow-surfaceShadowBefore,500)

  if hloc == 0 and nbminute>0:
    tab_d.append(d+timedelta(days=-1))
    tab_s.append(average/nbminute)
    tab_sb.append(average_b/nbminute)
    tab_sa.append(average_a/nbminute)
    tab_m.append(nbminute)
    average=0
    average_b=0
    average_a=0
    nbminute=0
  d = d + timedelta(minutes=1)

df = pd.DataFrame(data={"date": tab_h, "azimuth": tab_t0, "height": tab_t1})
df.to_csv("sun.csv", sep=';',index=False)

fig, (ax1) = plt.subplots(1, 1)

#ax1.plot(tab_d, tab_s)
ax1.set_ylabel('m2 moyen à l\'ombre')

ax2 = ax1.twinx()
line1, = ax1.plot(tab_d, tab_sb,color="blue",label="before")
ax1.plot(tab_d, tab_sa,color="blue",label="after")
line1.set_dashes([2, 2, 2, 2])  # 2pt line, 2pt break, 10pt line, 2pt break

ax2.plot(tab_d, tab_m, color="red",label="minutes")

ax1.set_xlabel('date')
ax1.set_ylabel('m2 moyen à l\'ombre', color='blue')
ax2.set_ylabel('minutes à l\'ombre', color='red')

ax1.grid(True)

#ax2.plot(tab_d, tab_m)
#ax2.set_xlabel('date')
#ax2.set_ylabel('minutes à l\'ombre')

plt.savefig('wall.png')

plt.show()