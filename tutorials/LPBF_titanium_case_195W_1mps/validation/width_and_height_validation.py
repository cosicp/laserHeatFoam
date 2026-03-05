import numpy as np
import matplotlib.pyplot as plt

dilip_et_al_width = np.array([
    [0.4997247475071271, 192.73581949232215],
    [0.7506233785106446, 160.75963647759323],
    [1.0013454888937328, 130.28580382325288],
    [1.2011069937213368, 122.05672203071137]
])

dilip_et_al_height = np.array([
    [0.5009434413435996, 175.3414274398405],
    [0.7508762700589977, 109.0449591652782],
    [1.002224260789795, 80.4520499303488],
    [1.200897790956006, 68.64740105432793]
])

laserMeltFoam_width = np.array([
    [0.5001854962451969, 197.533375117518],
    [0.7511738834963906, 152.24099028517708],
    [1.0012976188949723, 128.54653713569414],
    [1.2008796112272244, 119.16295832027578]
])

laserMeltFoam_height = np.array([
    [0.5014182053745723, 144.175930731201],
    [0.7518714485085366, 110.41872661222037],
    [1.0023064314874635, 87.84572942558248],
    [1.2007912733849544, 77.20644068722515]
])

# laserHeatFoam_width = np.array([
#     [0.5, ],
#     [0.75, ],
#     [1, 132],
#     [1.2,120]
# ])

laserHeatFoam_width = np.array([
    [0.5, 195],
    [0.75, 156],
    [1, 132],
    [1.2,120]
])

# laserHeatFoam_height = np.array([
#     [0.5, ],
#     [0.75, ],
#     [1, 85.5],
#     [1.2,78 ]
# ])

laserHeatFoam_height = np.array([
    [0.5, 110],
    [0.75, 94.5],
    [1, 85.5],
    [1.2,78 ]
])

plt.figure(1,figsize=(8, 6))
plt.plot(dilip_et_al_width[:, 0], dilip_et_al_width[:, 1], 'o-',markersize=10, label='Dilip et al. (2020)',linewidth=3)
plt.plot(laserMeltFoam_width[:, 0], laserMeltFoam_width[:, 1], 'o-',markersize=10,label='laserMeltFoam',linewidth=3,color='red')
plt.plot(laserHeatFoam_width[:, 0], laserHeatFoam_width[:, 1], 's-',markersize=10,label='laserHeatFoam',linewidth=3,color='green')
plt.xlabel('Scanning speed (m/s)',fontsize=14)
plt.ylabel('Average melt pool width (micrometers)',fontsize=14)
plt.title('Average melt pool width vs Scanning speed for 195W power',fontsize=14,pad=16)
plt.axis([0.4, 1.3, 40, 220])
plt.legend()
plt.grid(linestyle='--')

plt.figure(2,figsize=(8, 6))
plt.plot(dilip_et_al_height[:, 0], dilip_et_al_height[:, 1], 'o-',markersize=10, label='Dilip et al. (2020)',linewidth=3)
plt.plot(laserMeltFoam_height[:, 0], laserMeltFoam_height[:, 1], 'o-',markersize=10,label='laserMeltFoam',linewidth=3,color='red')
plt.plot(laserHeatFoam_height[:, 0], laserHeatFoam_height[:, 1], 's-',markersize=10,label='laserHeatFoam',linewidth=3,color='green')
plt.xlabel('Scanning speed (m/s)',fontsize=14)
plt.ylabel('Average melt pool height (micrometers)',fontsize=14)
plt.title('Average melt pool height vs Scanning speed for 195W power',fontsize=14,pad=16)
plt.axis([0.4, 1.3, 40, 220])
plt.legend()
plt.grid(linestyle='--')
plt.show()