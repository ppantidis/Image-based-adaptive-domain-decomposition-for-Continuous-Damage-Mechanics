import fileinput
import os
import matlab.engine


def create_input_file(mesh_file,FinalDisplacement,output_file, eng, path2use, modeOfFailure):
    
    count = 0
    count2 = 0
    count3=0

    print(f"Running '{__file__}'")

    with fileinput.input(files=(mesh_file)) as f:
        ctr=0
        mnl=open("meshnl.txt", "w")
        mcn=open("meshcn.txt", "w")
        for line in f:
            
            '''When you get to the nodes, start writing meshnl'''
            if line == '$Nodes\n':
                count =count+1
                continue
                
            if count>0:
                to=line.rstrip()
                if(len(to.split(' '))==3):
                    countstr =str(count)+ " "
                    mnl.write(countstr)
                    mnl.write(line)    
                    count = count+1
                    
                

            '''When you get to elements, start writing meshcn'''
            if line =='$Elements\n':
                mnl.close()
                count =0
                count2=1
                continue
            
            if count2 > 0:
                to=line.rstrip()
                if(len(to.split(' '))==5):
                    countstr =str(count2)+ " "
                    mcn.write(countstr)
                    line=line[line.index(' ')+1:len(line)]
                    mcn.write(line)
                    count2=count2+1
                    

            if line.split(' ')[0] =='Points_Used_to_Make_Domain:':
                count2=0
                mcn.close()
                count3=1
                pts=[[0 for x in range(2)] for y in range(int(line.split(' ')[1]))]
                continue

            if count3>0:
                pts[count3-1][0]=float(line.split(' ')[0])
                pts[count3-1][1]=float(line.split(' ')[1])
                count3+=1

    print(f"Creating ", output_file, "in MATLAB'")
    FinalDisplacement=matlab.double([FinalDisplacement])
    modeoffail=matlab.int8([int(modeOfFailure.split()[0])])
    direction=matlab.int8([int(modeOfFailure.split()[1])])
    minInDirection=matlab.double([float(modeOfFailure.split()[2])])
    maxInDirection=matlab.double([float(modeOfFailure.split()[3])])
    
    ptsMat=matlab.double(pts)

    eng.evalc("oldpath=path")
    eng.evalc("path(oldpath,'" + path2use +"')")
    _=eng.create_file_gmsh_contact(FinalDisplacement,ptsMat,output_file,modeoffail,direction,minInDirection,maxInDirection)

    eng.evalc("path(oldpath)")
    
    print(output_file, "has been created")
    os.remove("meshcn.txt")
    os.remove("meshnl.txt")
