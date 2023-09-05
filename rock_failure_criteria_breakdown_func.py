def mud_weight_breakdown(Depth, Svg, SHg, Shg, Pog, C, F_ang, v, Azi, Inc, T):
    import numpy as np
    import math

           

    F_ang = F_ang*math.pi/180
    Azi = Azi*math.pi/180
    Inc = Inc*math.pi/180


# Calculations of in-situ stresses and pore pressure at the depth of interest

    Sv = Svg*Depth
    SH = SHg*Depth
    Sh = Shg*Depth
    Po = Pog*Depth

# Determination of Mogi-Coulomb strength parameters

    a = ((2*2**0.5)/3)* C * math.cos(F_ang)
    b = ((2*2**0.5)/3) * math.sin(F_ang)

# Estimation of in situ stresses in the vicinity of the borehole for each borehole trajectory:

    S=[[SH,0,0], [0,Sh,0], [0,0,Sv]]
    rotation_mat=[[math.cos(Azi)*math.cos(Inc), math.sin(Azi)*math.cos(Inc), -math.sin(Inc)], [-math.sin(Azi), math.cos(Azi), 0], [math.cos(Azi)*math.sin(Inc), math.sin(Azi)*math.sin(Inc), math.cos(Inc)]]

    aa = np.array(S)
    bb = np.array(rotation_mat)
    c = bb.transpose()
    St = np.matmul(bb, aa)
    stress = np.matmul(St, c)

    sx = stress[0,0]
    sy = stress[1,1]
    sz = stress[2,2]
    sxy = stress[0,1]
    syz = stress[1,2]
    sxz = stress[0,2]

# Specifying the location of the maximum stress concentration:
# The orientation of the maximum and minimum tangential stresses:

    if sx == sy:
        teta1 = math.pi/4
    else:
        teta1 = 0.5*math.atan(2*sxy/(sx-sy))

    teta2 = teta1 + math.pi/2

# Identifying the angle that is associated with the maximum tangential stress:

    s_teta1 = sx+sy-2*(sx-sy)*math.cos(2*teta1)-4*sxy*math.sin(2*teta1)
    s_teta2 = sx+sy-2*(sx-sy)*math.cos(2*teta2)-4*sxy*math.sin(2*teta2)

# The location of the maximum stress concentration,

    if s_teta1 < s_teta2:
        teta = teta1
    else:
        teta = teta2

# The axial and shear stresses in θ-z plane at θ max:
    
    S_teta = sx+sy-2*(sx-sy)*math.cos(2*teta)-4*sxy*math.sin(2*teta)
   
    Sz = sz-v*(2*(sx-sy)*math.cos(2*teta)+4*sxy*math.sin(2*teta))
  
    S_teta_z = 2*(-sxz*math.sin(teta)+syz*math.cos(teta))

    

    mud_weight_breakdown = S_teta-Po+T
    

    return mud_weight_breakdown



    
        
        


