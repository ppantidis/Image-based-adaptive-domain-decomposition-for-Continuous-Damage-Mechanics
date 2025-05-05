import gmsh
import sys
import os
import matlab.engine

'''# predefined variables for testing. comment out later
filename = "1dom2.txt"
#path = "C:/Users/co978/Desktop/solver/domaindecomp"
path ="/Users/drebbel1z/Desktop/solver/domaindecomp_debug"
modeOfFailure = "1 2 0 5"
# first position: mode of failure(1, 2, or 3)
# will work for only 1 and 3 no mode 2 in 2D
# second position: direction (1 for x and 2 for y)
# third position: minimum value in direction
# fourth position: maximum value in direction
# space is added between values to allow for use of split method
# space is also used to allow float values to be used in third and fourth position
# (the delimiter could be anything)'''
 

p1 =os.path.dirname(os.path.abspath("__file__"))
p2=[]

for i in range (len(p1)):
    
    if p1[i]=="\\":
        p2.append('/')
    else:
        p2.append(p1[i])

for i in range(len(p2)-1,0,-1):
    a=p2.pop(i)
    if a =="/":
        break

dpath=''.join(p2)
 

def read_input_file(filename="1dom2.txt", path=dpath, modeOfFailure="1 2 0 5"):
    """ you can choose to start a new engine or connect to an engine you are already sharing.
    First one is slower
    if you want to use a shared engine, open another matlab instance and share it before you start running the code
    running the code from the same matlab engine as your solver may result in error
    """

    #eng = matlab.engine.start_matlab() #to start a new engine
    eng = matlab.engine.connect_matlab() #to connect with an already shared engine

    #change to directory that has the python files
    os.chdir(path+'/python_files')

    # if 'create_input_file' not in sys.modules:
    #    import create_input_file
    # this is on and off. sometimes when I call from matlab,
    # it doesn't work. the easy fix is to import the function everytime
    # you call this function

    import create_input_file #module for creating input file for the solver

    # open the user defined input file and read its content
    with open(path+'/python_files/'+filename) as f:
        lines = f.readlines()

    count = 0  # we are about to start reading the lines one by one

    # open file to write the model parameters in.
    parameters_files = open(path+'/parameters.txt', 'w')

    while True:  # copy every line into parameters file until you get to Number_of_domains:
        if lines[count] == "\n":
            parameters_files.write(lines[count])
            count += 1
        elif lines[count].split()[0] == "Number_of_domains:":
            break
        else:
            parameters_files.write(lines[count])
            count += 1

    # get number of domains to dictate how many input files to be made
    NumOfDomains = int(lines[count].split()[1])
    # write line to parameters file after you found the numOfDomains
    parameters_files.write(lines[count])
    count += 1

    while True:  # copy every line into parameters file until you get to FinalDisplacement:
        if lines[count] == "\n":
            parameters_files.write(lines[count])
            count += 1
        elif lines[count].split()[0] == "FinalDisplacement:":
            break
        else:
            parameters_files.write(lines[count])
            count += 1

    # get final displacement for whole domain
    FinalDisplacement = float(lines[count].split()[1])
    count += 2

    for Domain in range(NumOfDomains):
        parameters_files.write(lines[count])  # just for consistent styles really
        count += 1

        DomainID = int(lines[count].split()[1]) #get domain ID
        parameters_files.write(lines[count])
        count += 1

        while True:  # copy every line into parameters file until you get to Number_of_points:
            if lines[count] == "\n":
                parameters_files.write(lines[count])
                count += 1
            elif lines[count].split()[0] == "Number_of_points:":
                break
            else:
                parameters_files.write(lines[count])
                count += 1

        # find number of points, copy the points into a python list and use it in gmsh
        NumOfPoints = int(lines[count].split()[1])
        count += 3

        PtsInDomain = [[0 for x in range(4)] for y in range(NumOfPoints)] 
        #4 represents x,y,z,lc

        for pt in range(NumOfPoints):
            line = lines[count].split()
            for cpt in range(4):
                PtsInDomain[pt][cpt] = float(line[cpt])
            count += 1

        count += 1

        # GMSH INPUT FILE

        NameOfDomain = 'd'+str(DomainID)

        # Before using any functions in the Python API, Gmsh must be initialized:
        gmsh.initialize()
        gmsh.model.add(NameOfDomain)

        for i in range(NumOfPoints):
                gmsh.model.geo.addPoint(
                    PtsInDomain[i][0], PtsInDomain[i][1], PtsInDomain[i][2], PtsInDomain[i][3], i+1)

        for i in range(NumOfPoints):
            if i != NumOfPoints-1:
                gmsh.model.geo.addLine(i+1, i+2, i+1)
            else:
                gmsh.model.geo.addLine(i+1, 1, i+1)

        gmsh.model.geo.addCurveLoop([i+1 for i in range(NumOfPoints)], 1)

        #only one plane surface at a time, in the future this can be changed to create multiple plane surfaces

        gmsh.model.geo.addPlaneSurface([1], 1) 

        gmsh.model.geo.synchronize()

        # gmsh.model.mesh.refine()

        gmsh.model.mesh.setRecombine(2, 1)
        # For even better 2D (planar) quadrilateral meshes, you can try the
        # experimental "Frontal-Delaunay for quads" meshing algorithm, which is a
        # triangulation algorithm that enables to create right triangles almost
        # everywhere: J.-F. Remacle, F. Henrotte, T. Carrier-Baudouin, E. Bechet,
        # E. Marchandise, C. Geuzaine and T. Mouton. A frontal Delaunay quad mesh
        # generator using the L^inf norm. International Journal for Numerical Methods
        # in Engineering, 94, pp. 494-512, 2013. Uncomment the following line to try
        # the Frontal-Delaunay algorithms for quads:
        gmsh.option.setNumber("Mesh.Algorithm", 8)

        # The default recombination algorithm might leave some triangles in the mesh, if
        # recombining all the triangles leads to badly shaped quads. In such cases, to
        # generate full-quad meshes, you can either subdivide the resulting hybrid mesh
        # (with `Mesh.SubdivisionAlgorithm' set to 1), or use the full-quad
        # recombination algorithm, which will automatically perform a coarser mesh
        # followed by recombination, smoothing and subdivision.
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)

        gmsh.model.mesh.generate(2)

        gmsh.write(NameOfDomain+".msh")

        #clear file if it's already present
        if os.path.isfile('./' + NameOfDomain+'.txt'):
            os.remove(NameOfDomain+'.txt')

        #change .msh file to .txt file
        os.rename(NameOfDomain+'.msh', NameOfDomain+'.txt')

        #open the renamed file and add points used to make domain
        infile = open(NameOfDomain+'.txt', 'a')

        infile.write("Points_Used_to_Make_Domain: "+str(NumOfPoints))
        infile.write("\n")


        #add points used to make domain to the renamed file
        for y in range(NumOfPoints):
                for x in range(2):
                    infile.write(str(PtsInDomain[y][x])+' ')
                infile.write('\n')
        infile.close()

        # if '-nopopup' not in sys.argv:
        #  gmsh.fltk.run()

        gmsh.finalize()

        DomainFilename = NameOfDomain+'.txt'

        outputFile = path+'/input' + DomainFilename
        create_input_file.create_input_file(
                DomainFilename, FinalDisplacement, outputFile, eng, path + '/python_files', modeOfFailure)

        os.remove(DomainFilename)

    parameters_files.close()


read_input_file(filename, path, modeOfFailure)
