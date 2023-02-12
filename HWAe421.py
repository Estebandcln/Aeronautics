import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
from scipy import signal

# wing geometric data
b=10.48
S=21.55
A=5.1
lambd=0.5
CC=2.13
e=0.94

# steady flight conditions
g=9.81
h=9144
U0=182
XX_cg=0.3
rho=0.458
M=0.6
m=5440
I_xx=30048.7
I_yy=24415
I_zz=55299.8
I_xz=606
C_D0=0.0295
C_L0=0.299
C_T0=0.0295
qq=rho*(U0**2)/2
a=U0/M
theta0=0
gamma=0
a=303.3
epsilon=0.073

# I. Longitudinal Dynamic Stability########################################
# longitudinal stability derivatives
C_Du=0.0035
C_Ddeltae=0
C_Lalpha=5.16
C_Lalphadot=1.74
C_Lu=0.0881
C_Lq=3.9
C_mu=0.043
C_malpha=-1.09
C_malphadot=-4.98
C_mq=-11.7
C_Ldeltae=0.43
C_Mdeltae=-1.13
C_TdeltaT=0
C_LdeltaT=0
C_MdeltaT=0
X_u=-(qq*S/(m*U0))*(2*C_D0+C_Du)
X_w=(qq*S/(m*U0))*(C_L0*(1-2*C_Lalpha/(np.pi*e*A)))
Z_u=-(qq*S/(m*U0))*(2*C_L0+C_Lu)
Z_w=-(qq*S/(m*U0))*(C_D0+C_Lalpha)
Z_wdot=(qq*S*CC/(2*m*U0*U0))*(C_D0+C_Lalphadot)
Z_q=(qq*S*CC/(2*m*U0))*C_Lq
M_u=(qq*S*CC/(I_yy*U0))*C_mu
M_w=(qq*S*CC/(I_yy*U0))*C_malpha
M_wdot=(qq*S*CC*CC/(2*I_yy*U0*U0))*C_malphadot
M_q=(qq*S*CC*CC/(2*I_yy*U0))*C_mq
X_deltae=(qq*S/(m*U0))*C_Ddeltae
Z_deltae=(qq*S/(m*U0))*C_Ldeltae
M_deltae=(qq*S*CC/(I_yy*U0))*C_Mdeltae
X_deltaT=(qq*S/(m*U0))*C_TdeltaT
Z_deltaT=(qq*S/(m*U0))*C_LdeltaT
M_deltaT=(qq*S*CC/(I_yy*U0))*C_MdeltaT

# aircraft longitudinal matrices
Mat_A_long=np.array([[X_u,X_w,0,
                      -g*np.cos(theta0)],
                     [Z_u,Z_w,U0,-g*np.sin(theta0)],
                     [M_u+Z_u*M_wdot,M_w+Z_w*M_wdot,M_q+U0*M_wdot,0],
                     [0,0,1,0]])
Mat_B_long=np.array([[X_deltae,X_deltaT],[Z_deltae,Z_deltaT],[M_deltae+Z_deltae*M_wdot,M_deltaT+Z_deltaT*M_wdot],[0,0]])
#characteristic equation

Mat_A_long_symb=sp.Matrix(Mat_A_long)
eqlong=Mat_A_long_symb.charpoly()
'''
np.poly(Mat_A_long)
'''

#eigenvalues

eigenvalueslong=np.linalg.eigvals(Mat_A_long)
#factored polynomial= (lambda**2 + 2*0.69484*lambda + (0.69484**2 + 1.3433**2))*(lambda**2 + 2*0.010054*lambda + (0.010054**2 + 0.12482**2))
'''
(λ^2+ 2*ksi_sp*omega_sp "λ +" omega_sp^2 )(λ^2+2ksi_ph omega_ph "λ+ " omega_ph^2 )
'''
#curves of longitudinal motion
#short period mode------------------------------------------------
w0=12
omega_sp=math.sqrt(np.poly(Mat_A_long)[2])
ksi_sp=np.poly(Mat_A_long)[1]/(2*omega_sp)

def spm_w(t):
    return (math.exp(-omega_sp*ksi_sp*t)*(w0/10)*math.cos(omega_sp*math.sqrt(1-ksi_sp**2)*t)+(10*omega_sp*ksi_sp)/(omega_sp*math.sqrt(1-ksi_sp**2))*math.sin(omega_sp*math.sqrt(1-ksi_sp**2)*t))
def spm_q(t):
    return (math.exp(-omega_sp*ksi_sp*t)*(theta0/10)*math.cos(omega_sp*math.sqrt(1-ksi_sp**2)*t)+(theta0*omega_sp*ksi_sp)/(omega_sp*math.sqrt(1-ksi_sp**2))*math.sin(omega_sp*math.sqrt(1-ksi_sp**2)*t))

time_sp=[i for i in range (0,100)]
spmlst_w=[spm_w(t) for t in time_sp]
spmlst_q=[spm_q(t) for t in time_sp]
plt.figure(1,dpi=300)
plt.title('Short period')
plt.plot(time_sp, spmlst_w,label='Short period w')
plt.plot(time_sp, spmlst_q,label='Short period q')
plt.legend()
plt.show()
#phugoid mode----------------------------------------------------
'''regarder tuto 1'''
omega_ph=math.sqrt(np.poly(Mat_A_long)[2])
ksi_ph=omega_ph/2*0.010054
time_ph=[i for i in range (0,60)]

def phm_u(t):
    return (math.exp(-omega_ph*ksi_ph*t)*(U0/10)*math.cos(omega_ph*math.sqrt(1-ksi_ph**2)*t)+(10*omega_ph*ksi_ph)/(omega_ph*math.sqrt(1-ksi_ph**2))*math.sin(omega_ph*math.sqrt(1-ksi_ph**2)*t))
def phm_theta(t):
    return (math.exp(-omega_ph*ksi_ph*t)*(theta0/10)*math.cos(omega_ph*math.sqrt(1-ksi_ph**2)*t)+(theta0*omega_ph*ksi_ph)/(omega_ph*math.sqrt(1-ksi_ph**2))*math.sin(omega_ph*math.sqrt(1-ksi_ph**2)*t))

phmlst_u=[phm_u(t) for t in time_ph]
plt.figure(2,dpi=300)
plt.title('Phugoid')
plt.plot(time_ph, phmlst_u,label='Phugoid u')
phmlst_theta=[phm_theta(t) for t in time_ph]
plt.plot(time_ph, phmlst_theta,label='Phugoid theta')
plt.legend()
plt.show()
#approximation-------------------------------------------------------------------------
#Short period mode approximation--------------------------
#Short period Matrix
Mat_A_sp=np.array([[Z_w,U0],[M_w+Z_w*M_wdot,M_q+U0*M_wdot]])
Mat_B_sp=np.array([[Z_deltae,Z_deltaT],[M_deltae+Z_deltae*M_wdot,M_deltaT+Z_deltaT*M_wdot]])

#characteristic equation
Mat_A_sp_symb=sp.Matrix(Mat_A_sp)
eq_sp=Mat_A_sp_symb.charpoly()
'''
np.poly(Mat_A_sp)
'''
#eigenvalues
eigenvaluessp=np.linalg.eigvals(Mat_A_sp)
omega_approx_sp=math.sqrt(np.poly(Mat_A_sp)[2])
ksi_approx_sp=np.poly(Mat_A_sp)[1]/(2*omega_approx_sp)

#period
T_sp=2*math.pi/(omega_approx_sp*math.sqrt(1-ksi_approx_sp**2))

#transfer functions
#Phugoid mode approximation-----------------------------
#Phugoid Matrix
Mat_A_ph=np.array([[X_u,-g*math.cos(theta0)],[-Z_u/U0,g*math.sin(theta0)/U0]])
Mat_B_ph=np.array([[X_deltae,X_deltaT],[-Z_deltae/U0,-Z_deltaT/U0]])

#characteristic equation
Mat_A_ph_symb=sp.Matrix(Mat_A_ph)
eq_ph=Mat_A_ph_symb.charpoly()
'''
np.poly(Mat_A_ph)
'''
#eigenvalues
eigenvaluesph=np.linalg.eigvals(Mat_A_ph)
omega_approx_ph=math.sqrt(np.poly(Mat_A_ph)[2])
ksi_approx_ph=np.poly(Mat_A_ph)[1]/2*omega_approx_ph

#period
T_ph=2*math.pi/(omega_approx_ph*math.sqrt(1-ksi_approx_ph**2))

'''approximation response tuto 2'''


# II. Lateral Dynamic Stability###########################################

# wing geometric data
b=10.36
S=21.37
A=5.03
lambd=0.5
CC=2.13
e=0.94

# steady flight conditions
g=9.81
h=12192
U0=206.3
XX_cg=0.32
rho=0.303
M=0.7
m=5785
I_xx=37979
I_yy=25500
I_zz=63750
I_xz=1763
C_D0=0.256
C_L0=1.64
C_T0=0.256
qq=rho*(U0**2)/2
a=U0/M
theta0=2.7 
gamma=0
a=303.3
epsilon=0.073

# lateral stability derivatives
C_lbeta=-0.11
C_yp=0
C_lp=-0.45
C_lr=0.16
C_ldeltar=-0.019
C_ldeltaa=0.178
C_nbeta=0.127
C_np=-0.008
C_nr=-0.2
C_ndeltaa=-0.02
C_ndeltar=-0.074
C_ybeta=-0.73
C_yr=0.4
C_ydeltaa=0
C_ydeltar=0.14
Y_v=qq*S*C_ybeta/(m*U0)
Y_p=qq*S*b*C_yp/(2*m*U0)
Y_r=qq*S*b*C_yr/(2*m*U0)
L_v=qq*S*b*C_lbeta/(I_xx*U0)
L_p=qq*S*b*b*C_lp/(2*I_xx*U0)
L_r=qq*S*b*b*C_lr/(2*I_xx*U0)
N_v=qq*S*b*C_nbeta/(I_zz*U0)
N_p=qq*S*b*b*C_np/(2*I_zz*U0)
N_r=qq*S*b*b*C_nr/(2*I_zz*U0)
Y_deltar=qq*S*C_ydeltar/(m*U0)
L_deltar=qq*S*b*C_ldeltar/(I_xx*U0)
N_deltar=qq*S*b*C_ndeltar/(I_zz*U0)
Y_deltaa=qq*S*C_ydeltaa/(m*U0)
L_deltaa=qq*S*b*C_ldeltaa/(I_xx*U0)
N_deltaa=qq*S*b*C_ndeltaa/(I_zz*U0)

# aircraft lateral matrices
Mat_A_lat=np.array([[Y_v, Y_p, (Y_r-U0), g*math.cos(theta0)],[L_v, L_p, L_r,0],[N_v, N_p, N_r, 0],[0, 1, 0, 0]])
Mat_B_lat=np.array([[Y_deltar, Y_deltaa],[L_deltar, L_deltaa],[N_deltar, N_deltaa],[0,0]])

#characteristic equation
Mat_A_lat_symb=sp.Matrix(Mat_A_lat)
eqlat=Mat_A_lat_symb.charpoly()
'''
np.poly(Mat_A_lat)
'''

#eigenvalues
eigenvalueslat=np.linalg.eigvals(Mat_A_lat)

#eigenvectors
eigenvectorslat=np.linalg.eig(Mat_A_lat)

#spiral mode
lambda_s=eigenvalueslat[2].real
eigenvectors_s=eigenvectorslat[1]
time_s=[i for i in range (0,30)]

for i in range(4):
    def spiral_eigen(t):
        return (eigenvectors_s[i]*math.exp(lambda_s*t))
    spiralvalues=[spiral_eigen(j) for j in time_s]
    
    """
    plt.figure(i+3,dpi=300)
    plt.plot(time_s, spiralvalues)
    plt.title(label='Spiral mode '+str(i+1))
    plt.grid()
    """
    
    
omega_s=math.sqrt(np.poly(Mat_A_lat)[3])
ksi_s=np.poly(Mat_A_lat)[1]/(2*omega_s)
plt.figure(12,dpi=300)

def spiral(t):
    return(math.exp(lambda_s*t)*(math.cos(omega_s*t)+ksi_s*math.sin(omega_s*t)))

Svalues=[spiral(j) for j in time_s]
plt.plot(time_s, Svalues)
plt.grid()
plt.title('Spiral')


#time constant
T_s=1/lambda_s

#rolling mode
listerolling=[]
lambda_r=eigenvalueslat[1].real
eigenvectors_r=eigenvectorslat[1]
time_r=[i for i in range (0,100)]

omega_r=math.sqrt(np.poly(Mat_A_lat)[2])   # 2 ?
ksi_r=np.poly(Mat_A_lat)[1]/(2*omega_r)
w0=12


for i in range(4):
    def rolling_eigen(t):
        return (eigenvectors_r[i]*math.exp(lambda_r*t))
    rollingvalues=[rolling_eigen(j) for j in time_r]
    """
    plt.figure(i+7,dpi=300)
    plt.plot(time_r, rollingvalues)
    plt.title(label='Rolling mode '+str(i+1))
    plt.grid()
    """
    

omega_r=math.sqrt(np.poly(Mat_A_lat)[2])
ksi_r=np.poly(Mat_A_lat)[1]/(2*omega_r)
plt.figure(13,dpi=300)

def rolling(t):
    return(math.exp(lambda_r*t)*(math.cos(omega_r*t)+ksi_r*math.sin(omega_r*t)))

Rvalues=[rolling(j) for j in time_r]
plt.plot(time_r, Rvalues)
plt.grid()
plt.title('Rolling')

#time constant
T_r=1/lambda_r

#dutch roll mode

lambdareal_dr=eigenvalueslat[0].real

omega_dr=math.sqrt(0.641176)
ksi_dr=0.293342/(2*omega_dr)

plt.figure(14,dpi=300)

def dutchroll(t):
    return(math.exp(lambdareal_dr*t)*(math.cos(omega_dr*t)+ksi_dr*math.sin(omega_dr*t)))

time_dr=[i for i in range (0,100)]
DRvalues=[dutchroll(j) for j in time_dr]
plt.plot(time_dr, DRvalues)
plt.grid()
plt.title('Dutchroll')



#curves of lateral motion

#transfer functions
'''
def lat(X,t):
    eta=1
    return (Mat_A_lat*X+Mat_B_lat*eta)
t=np.linspace(0,10,100)
y=odeint(lat,0,t)
plt.figure(1)
plt.plot(t,y)

I=np.eye(3)
num=[Mat_B_lat*np.linalg.inv(Mat_A_lat).T]
den=[1]

TF=signal.TransferFunction(num,den)

print(TF)
'''