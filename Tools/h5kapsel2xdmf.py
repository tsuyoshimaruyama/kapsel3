#!/usr/bin/python
import sys
import os.path
import xml.etree.ElementTree as et
from xml.dom import minidom

floatType = ['Float', '4']
intType   = ['Int', '4']
#read all lines in conf file and purge all comments and blank lines
def readConf(fname):
    fconf = open(fname, "r")
    lines = fconf.readlines()
    cline = []
    for line in lines:
        if(line.strip() and not line.lstrip().startswith('#')):
            dmy = line.split("#", 1)[0]
            cline.append(dmy.lstrip().rstrip())
    return cline

def printList(l):
    s=str(l[0])
    for i in range(1, len(l)):
        s+=" "+str(l[i])
    return s

#read id tags
def readTag(line):
    idl = line.index("[")
    idr = line.index("]")
    return line[idl+1:idr].lower()

#read argument list for num arguments separated by spaces
def readArg(line, num):
    dmy = line.split(" ")[0:num]
    if(num != len(dmy)):
        print "Format Error!"
        print num
        print line
        exit(-1)
    if num > 1:
        return dmy
    else:
        return dmy[0]

def parseConf(fname):
    h5file =""
    time = [0.0, 0.0, 0.0]
    dims = [0, 0, 0]
    time = [0, 0, 0]
    nump= 0
    scalar_field = []
    vector_field = []
    tensor6_field= []
    tensor_field= []

    scalar_particle =[]
    vector_particle =[]
    tensor6_particle=[]
    tensor_particle=[]
    lines = readConf(fname)
    i = 0
    while i < len(lines):
        tag = readTag(lines[i])
        if(tag == "filename"):
            h5file = readArg(lines[i+1], 1)
        elif(tag == "time series"):
            time = readArg(lines[i+1], 3)
        elif(tag == "grid dimensions"):
            dims  = readArg(lines[i+1], 3)
        elif(tag == "particle number"):
            nump = readArg(lines[i+1], 1)
        elif(tag == "field scalar"):
            scalar_field.append(readArg(lines[i+1], 1))
        elif(tag == "field vector"):
            vector_field.append(readArg(lines[i+1], 4))
        elif(tag == "field tensor6"):
            tensor6_field.append(readArg(lines[i+1], 7))
        elif(tag == "field tensor"):
            tensor_field.append(readArg(lines[i+1], 10))
        elif(tag == "particle grid"):
            particle_grid = readArg(lines[i+1], 1)
        elif(tag == "particle scalar"):
            scalar_particle.append(readArg(lines[i+1], 1))
        elif(tag == "particle vector"):
            vector_particle.append(readArg(lines[i+1], 1))
        elif(tag == "particle tensor6"):
            tensor6_particle.append(readArg(lines[i+1], 1))
        elif(tag == "particle tensor"):
            tensor_particle.append(readArg(lines[i+1], 1))
        else:
            print "Format Error: Unknown Tag"
            print lines[i-1]
            print lines[i]
            exit(-1)
        i+=2

    print "data file  :", h5file
    print "time series:", time
    print "grid dims  :", dims
    print "nump       :", nump

    print "\n***Field Data***"

    print "Scalar Fields:", len(scalar_field)
    for scalar in scalar_field:
        print "\t",scalar

    print "Vector Fields:", len(vector_field)
    for vector in vector_field:
        print "\t",vector[:-1],"@",vector[-1]

    print "Tensor6 Fields:", len(tensor6_field)
    for tensor6 in tensor6_field:
        print "\t",tensor6[:-1], "@", tensor6[-1]

    print "Tensor Fields:", len(tensor_field)
    for tensor in tensor_field:
        print "\t",tensor[:-1], "@", tensor[-1]
    
    print "\n***Particle Data***"

    print "Particle Grid:", particle_grid
    print "Particle Scalar:", len(scalar_particle)
    for scalar in scalar_particle:
        print "\t",scalar
    
    print "Particle Vector:", len(vector_particle)
    for vector in vector_particle:
        print "\t",vector
    
    print "Particle Tensor6:", len(tensor6_particle)
    for tensor6 in tensor6_particle:
        print "\t",tensor6
    
    print "Particle Tensor:", len(tensor_particle)
    for tensor in tensor_particle:
        print "\t",tensor
    return [[h5file, time, dims, nump],\
            [scalar_field, vector_field, tensor6_field, tensor_field],\
            [particle_grid, scalar_particle, vector_particle,\
             tensor6_particle, tensor_particle]]


#add xdmf doctype tag to xml and pretty print
def finalizexdmf(elem):
    rough_string = "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" + et.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def prettyprintxml(elem):
    rough_string = et.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    reparsed = reparsed.toprettyxml(indent="  ")
    print reparsed
    return reparsed

#add data item to xdmf node
def xdmf_addDataItem(root, ddims, dformat, dtype):
    data = et.SubElement(root, 'DataItem')
    data.set('NumberType', dtype[0])
    data.set('Precision', dtype[1])
    data.set('Dimensions', ddims)
    data.set('Format', dformat)
    return data
#add data item function to xdmf node
def xdmf_addDataItemFunc(root, fdims, ffunc):
    func = et.SubElement(root, 'DataItem')
    func.set('ItemType', 'Function')
    func.set('Dimensions', fdims)
    func.set('Function', ffunc)
    return func
#add attribute item to xdmf node
def xdmf_addAttribute(root, atype, aname):
    at = et.SubElement(root, 'Attribute')
    at.set('AttributeType', atype)
    at.set('Name', aname)
    at.set('Center', 'Node')
    return at
#init xdmf xml tree    
def xdmf_init(time, name):
    xdmf = []
    xdmf.append(et.Element('Xdmf'))
    xdmf[-1].set('Version', '2.0')
    xdmf.append(et.SubElement(xdmf[-1], 'Domain'))
    xdmf.append(et.SubElement(xdmf[-1], 'Grid'))
    xdmf[-1].set('CollectionType', 'Temporal')
    xdmf[-1].set('GridType', 'Collection')
    xdmf[-1].set('Name', name)

    xdmf.append(et.SubElement(xdmf[-1], 'Time'))
    xdmf[-1].set('TimeType','HyperSlab')
    xdmf.append(xdmf_addDataItem(xdmf[-1], '3', 'XML', floatType))
    xdmf[-1].text = time
    xdmf.pop()
    xdmf.pop()
    return xdmf

#create xml fluid grid
def xdmf_fluid_grid(coord_loc, field_loc, dims, \
                    scalar_field, vector_field, tensor6_field, tensor_field):
    grid =[] 
    grid.append(et.Element("Grid"))
    grid[-1].set('Name', 'Fluid Frame')
    grid[-1].set('GridType', 'Uniform')

    #grid->topology
    grid.append(et.SubElement(grid[-1], 'Topology'))
    grid[-1].set('TopologyType', '3DSMesh')
    grid[-1].set('NumberOfElements', dims)
    grid.pop()

    #grid->geometry
    grid.append(et.SubElement(grid[-1], 'Geometry'))
    grid[-1].set('GeometryType', 'X_Y_Z')
    #grid->geometry->x axis
    grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
    grid[-1].text = coord_loc+"x"
    grid.pop()
    #grid->geometry->y axis
    grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
    grid[-1].text = coord_loc+"y"
    grid.pop()
    #grid->geometry->z axis
    grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
    grid[-1].text = coord_loc+"z"
    grid.pop()
    grid.pop()

    #attributes
    for scalar in scalar_field:
        grid.append(xdmf_addAttribute(grid[-1], 'Scalar', scalar))
        grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
        grid[-1].text = field_loc+scalar
        grid.pop()
        grid.pop()
    dmy_dims = dims + " 3"
    for vector in vector_field:
        grid.append(xdmf_addAttribute(grid[-1], 'Vector', vector[-1]))
        grid.append(xdmf_addDataItemFunc(grid[-1], dmy_dims, 'JOIN($0, $1, $2)'))
        for i in range(3):
            grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
            grid[-1].text = field_loc+vector[i]
            grid.pop()
        grid.pop()
        grid.pop()
    dmy_dims = dims + " 6"
    for tensor6 in tensor6_field:
        grid.append(xdmf_addAttribute(grid[-1], 'Tensor6', tensor6[-1]))
        grid.append(xdmf_addDataItemFunc(grid[-1], dmy_dims, 'JOIN($0, $1, $2, $3, $4, $5)'))
        for i in range(6):
            grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
            grid[-1].text = field_loc+tensor6[i]
            grid.pop()
        grid.pop()
        grid.pop()
    dmy_dims = dims + " 9"
    for tensor in tensor_field:
        grid.append(xdmf_addAttribute(grid[-1], 'Tensor', tensor[-1]))
        grid.append(xdmf_addDataItemFunc(grid[-1], dmy_dims, 'JOIN($0, $1, $2, $3, $4, $5, $6, $7, $8)'))
        for i in range(9):
            grid.append(xdmf_addDataItem(grid[-1], dims, "HDF", floatType))
            grid[-1].text = field_loc+tensor[i]
            grid.pop()
        grid.pop()
        grid.pop()

    return grid[0]

def xdmf_part_grid(part_loc, coord_loc, nump, scalar_particle, vector_particle, \
                   tensor6_particle, tensor_particle):
    grid=[]
    grid.append(et.Element("Grid"))
    grid[-1].set('Name', 'Particle Frame')
    grid[-1].set('GridType', 'Uniform')
    
    #grid -> topology
    grid.append(et.SubElement(grid[-1], 'Topology'))
    grid[-1].set('TopologyType', 'Polyvertex')
    grid[-1].set('NodesPerElement', nump)
    grid.pop()
    
    #grid -> geometry
    grid.append(et.SubElement(grid[-1], 'Geometry'))
    grid[-1].set('GeometryType', 'XYZ')
    #grid -> geometry -> coordinates
    dmy_dims = nump + " 3"
    grid.append(xdmf_addDataItem(grid[-1], dmy_dims, "HDF", floatType))
    grid[-1].text = part_loc + coord_loc
    grid.pop()
    grid.pop()

    dmy_dims = nump + " 1"
    for scalar in scalar_particle:
        grid.append(xdmf_addAttribute(grid[-1], 'Scalar', scalar))
        grid.append(xdmf_addDataItem(grid[-1], dmy_dims, "HDF", intType))
        grid[-1].text = part_loc + scalar
        grid.pop()
        grid.pop()

    dmy_dims = nump + " 3"
    for vector in vector_particle:
        grid.append(xdmf_addAttribute(grid[-1], 'Vector', vector))
        grid.append(xdmf_addDataItem(grid[-1], dmy_dims, "HDF", floatType))
        grid[-1].text = part_loc + vector
        grid.pop()
        grid.pop()

    dmy_dims = nump + " 6"
    for tensor6 in tensor6_particle:
        grid.append(xdmf_addAttribute(grid[-1], 'Tensor6', tensor6))
        grid.append(xdmf_addDataItem(grid[-1], dmy_dims, "HDF", floatType))
        grid[-1].text = part_loc + tensor6
        grid.pop()
        grid.pop()

    dmy_dims = nump + " 9"
    for tensor in tensor_particle:
        grid.append(xdmf_addAttribute(grid[-1], 'Tensor', tensor))
        grid.append(xdmf_addDataItem(grid[-1], dmy_dims, "HDF", floatType))
        grid[-1].text = part_loc + tensor
        grid.pop()
        grid.pop()
    return grid[0]

##### MAIN ####
if (not len(sys.argv) == 2) or (not os.path.isfile(str(sys.argv[1]))):
    print "Usage: ./kapselxdmf.py conf_file"
    exit(-1)
conf_file = str(sys.argv[1])
[[h5file, time, dims, nump],\
 [scalar_field, vector_field, tensor6_field, tensor_field],\
 [particle_grid, scalar_particle, vector_particle,\
  tensor6_particle, tensor_particle]] = parseConf(conf_file)

fx = xdmf_init(printList(time), "KAPSEL Fluid Trajectory")
px = xdmf_init(printList(time), "KAPSEL Particle Trajectory")
grid_loc=h5file+".h5:/coord_data/"

for i in range(int(time[2])):
    frame = h5file + ".h5:/frame_"+str(i)
    fluid_loc   = frame+"/field_data/"
    particle_loc= frame+"/particle_data/"
    fx[-1].append(xdmf_fluid_grid(grid_loc, fluid_loc, printList(dims), \
                                  scalar_field, vector_field, \
                                  tensor6_field, tensor_field))
    px[-1].append(xdmf_part_grid(particle_loc, particle_grid, nump, \
                                 scalar_particle, vector_particle, \
                                 tensor6_particle, tensor_particle))
    
fxx = open(h5file+'_fluid.xdmf', 'w')
pxx = open(h5file+'_particle.xdmf', 'w')
fxx.write(finalizexdmf(fx[0]))
pxx.write(finalizexdmf(px[0]))
