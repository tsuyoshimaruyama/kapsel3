showRigidCOM=0
type=$constitutive_eq.type
print (type)
if type == "Navier_Stokes" :
	dx=$constitutive_eq.Navier_Stokes.DX
elif type == "Shear_Navier_Stokes" :
	dx=$constitutive_eq.Shear_Navier_Stokes.DX
elif type == "Shear_Navier_Stokes_Lees_Edwards" :
	dx=$constitutive_eq.Shear_Navier_Stokes_Lees_Edwards.DX
elif type == "Electrolyte" :
	dx=$constitutive_eq.Electrolyte.DX
elif type == "Navier_Stokes_FDM" :
	dx=$constitutive_eq.Navier_Stokes_FDM.DX
elif type == "Navier_Stokes_Cahn_Hilliard_FDM" :
	dx=$constitutive_eq.Navier_Stokes_Cahn_Hilliard_FDM.DX
	potential = $constitutive_eq.Navier_Stokes_Cahn_Hilliard_FDM.Potential.type
	if potential == "Landau":
		cval = $constitutive_eq.Navier_Stokes_Cahn_Hilliard_FDM.Potential.Landau.composition_ratio
	elif potential == "Flory_Huggins":
		cval = $constitutive_eq.Navier_Stokes_Cahn_Hilliard_FDM.Potential.Flory_Huggins.composition_ratio
elif type == "Shear_Navier_Stokes_Lees_Edwards_FDM" :
	dx=$constitutive_eq.Shear_Navier_Stokes_Lees_Edwards_FDM.DX
elif type == "Shear_NS_LE_CH_FDM" :
	dx=$constitutive_eq.Shear_NS_LE_CH_FDM.DX
	potential = $constitutive_eq.Shear_NS_LE_CH_FDM.Potential.type
	if potential == "Landau":
		cval = $constitutive_eq.Shear_NS_LE_CH_FDM.Potential.Landau.composition_ratio
	elif potential == "Flory_Huggins":
		cval = $constitutive_eq.Shear_NS_LE_CH_FDM.Potential.Flory_Huggins.composition_ratio
NX=2**$mesh.NPX
NY=2**$mesh.NPY
NZ=2**$mesh.NPZ
LX=dx*NX
LY=dx*NY
LZ=dx*NZ
print (LX,LY,LZ)
objType=$object_type.type
if objType=="spherical_particle":
	Ns=getArray($object_type.spherical_particle.Particle_spec[])
elif objType=="chain":
	Ns=getArray($object_type.chain.Chain_spec[])
elif objType == "rigid":
	Ns=getArray($object_type.rigid.Rigid_spec[])
size_Ns=len(Ns)
RAD=($A*dx)*1.
for i in range(size_Ns):
	print (Ns[i][0])
spat=[
	[1.0, 1.0, 1.0, 1.0, RAD],
	[1.0, 0.0, 0.0, 1.0, RAD],
	[0.0, 1.0, 0.0, 1.0, RAD],
	[0.0, 0.0, 1.0, 1.0, RAD],
	[1.0, 1.0, 0.0, 1.0, RAD],
	[0.0, 1.0, 1.0, 1.0, RAD],
	[1.0, 0.0, 1.0, 1.0, RAD]
	]
n_offset = 0
if objType=="rigid" or objType=="chain":
        if showRigidCOM == 0:
                #species
                for i in range(size_Ns):
                        #chains
                        for m in range(Ns[i][1]):
                                #beads
                                for n in range(Ns[i][0]):
                                        r=$Particles[n_offset+n].R
                                        sphere(r,spat[i%len(spat)])
                                n_offset += Ns[i][0]
        elif showRigidCOM == 1:
                #species
                for i in range(size_Ns):
                        #chains
                        for m in range(Ns[i][1]):
                                r=$RigidParticles[n_offset+m].R
                                sphere(r,spat[i%len(spat)])
                        n_offset+=Ns[i][1]
else:
        #species
        for i in range(size_Ns):
                #particles
                for n in range(Ns[i][0]):
                        r=$Particles[n_offset+n].R
                        sphere(r,spat[i%len(spat)])
                n_offset+=Ns[i][0]

coord_list_list = [[0.0, LX],[0.0, LY],[0.0, LZ]]
div_list = [NX, NY, NZ]
mesh_obj = meshfield("regular", coord_list_list, div_list)
for i in range (NX):
	for j in range (NY):
		for k in range (NZ):
			if(type == "Navier_Stokes_Cahn_Hilliard_FDM" or type == "Shear_NS_LE_CH_FDM"):
				mesh_obj.set([i,j,k],$PSI[i][j][k].psi)
			else:
				mesh_obj.set([i,j,k],0.0)

mesh_obj.isovalue(cval, [0.2, 0.7, 1.0, 0.8])
mesh_obj.draw(frame=[1.0, 1.0, 1.0, 0.5])
