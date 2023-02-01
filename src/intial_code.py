import matplotlib.pyplot as plt
import numpy as np
import math as mt

#--------------------------------FUNCTIONS------------------------------------

#simply supported beam intermediate load calc with super position
def simp_im_def (x,a1,b1,a2,b2,w1,w2):
    imconst1 = w1*b1/(6*geo_const*axl)
    imconst2 = w2*b2/(6*geo_const*axl)
    ans1=0
    ans2=0
    if x <= a1:
        ans1= -imconst1*b1*x*(pow(axl,2)-pow(x,2)-pow(b1,2))
    if x > a1:
        ans1= -imconst1*b1*((axl/b1*pow((x-a1),3))+(x*(pow(axl,2)-pow(b1,2)))-pow(x,3))
    if x <= a2:
        ans2= -imconst2*b2*x*(pow(axl,2)-pow(x,2)-pow(b2,2))
    if x > a2:
        ans2= -imconst2*b2*((axl/b2*pow((x-a2),3))+(x*(pow(axl,2)-pow(b2,2)))-pow(x,3))
    return ans1 + ans2

#distance formula
def dist_form(x1,x2):
    return (mt.sqrt(pow(x2[0]-x1[0],2)-pow(x2[1]-x1[1],2)))

#strain
def strain(L1,L2):
    return abs(L1-L2)/L1

#------------------------------CONSTANTS--------------------------------------
counter1 = 0
counter2 = 0
calc_reso=1 #resolution of calcution in inches
reso_count=1/calc_reso
#----------------------------axle porperties----------------------------------
axl=48.0 #axle length in inches`
OD=5.0 #outer diameter
thick=0.25 #axle tubing wall thickness
ID=OD-thick #inner diameter
#---------------leaf spring position of intermediate forces-------------------
sproff = 12 #wheel to leaf spring value
ima1 = sproff
imb1 = axl - sproff
ima2= imb1
imb2 = ima1
#----------------sensor dimension----------------------------------------------
senlen = 3.5
senend1 = axl/2 - senlen/2
senend2 = axl/2 + senlen/2
senpltx = [senend1,senend2]
#-------------------------Geometric Constants----------------------------------
inrt = (3.1415*(pow(OD,4)-pow(ID,4)))/64 #axle moment of inertica as pipe
flexmod = 29000000 #modulus of elasticity of steel, 29 10^6 psi
geo_const = flexmod*inrt #E*I
#------------------------------Forces-------------------------------------
maxweight=16000 #maximum weight on axle in lbs F= maxweight/2
skew = 500 #cannot be greater 8000
 #weight "skew"
wL = (maxweight/2) + skew #L spring weight
wR = maxweight - wL #R spring weight
#---------------------------------BODY----------------------------------------

axlpos=np.arange(0,axl,calc_reso) #set xarray, axle positrion
TR=np.arange(0,axl,calc_reso) #set same length
iTR=np.arange(0,axl,calc_reso) #set same length
calc_len = len(axlpos)

#check with skew
for check in range(0,calc_len):
    pos = axlpos[check]
    TR[check]=simp_im_def(pos,ima1,imb1,ima2,imb2,wL,wR)

#for distance calc
ssnrpos1 = simp_im_def(senend1,ima1,imb1,ima2,imb2,wL,wR)
ssnrpos2 = simp_im_def(senend2,ima1,imb1,ima2,imb2,wL,wR)
sensor1 = [senend1,ssnrpos1]
sensor2 = [senend2,ssnrpos2]
ssenplt1 = [ssnrpos1,ssnrpos2]
nsenl = dist_form(sensor1,sensor2)
ustrain = strain(senlen,nsenl)*1000000
print("Skewed Sensor Compressed Length: "+str(round(nsenl,6)))
print ("Skewed Axle Micostrain:"+ str(round(ustrain,2)))

#check without skew
for check in range(0,calc_len):
    pos = axlpos[check]
    iTR[check]=simp_im_def(pos,ima1,imb1,ima2,imb2,maxweight/2,maxweight/2)

issnrpos1 = simp_im_def(senend1,ima1,imb1,ima2,imb2,maxweight/2,maxweight/2)
issnrpos2 = simp_im_def(senend2,ima1,imb1,ima2,imb2,maxweight/2,maxweight/2)
issenplt1 = [issnrpos1,issnrpos2]
isensor1 = [senend1,issnrpos1]
isensor2 = [senend2,issnrpos2]
insenl = dist_form(isensor1,isensor2)
iustrain = strain(senlen,insenl)*1000000
print ("Ideal Axle Micostrain:"+ str(round(iustrain,2)))
print("Ideal Sensor Compressed Length: "+str(round(insenl,6)))
print("Percent Error: " + str(round(abs((iustrain-ustrain)/iustrain*100),2)) +"%")


#print(ustrain)

#get distance for sensor body
#fig, ax = plt.subplots()
plt.plot(axlpos,TR,label="Total Reaction")
plt.plot(axlpos,iTR,label = "even load")
plt.plot(senpltx,ssenplt1,label="skew sensor",linewidth=4)
plt.plot(senpltx,issenplt1,label="idea sensor",linewidth=4)
plt.title("Axle Position VS Deflection")
plt.xlabel('axle position')
plt.ylabel('deflection')
plt.legend()
plt.figure(num=None,figsize=[10,4])
plt.show()