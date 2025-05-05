import gmsh
import sys
import os
import matlab.engine
import create_input_file


# predefined variables for testing. comment out later
eng = matlab.engine.connect_matlab()
filename = "2dom.txt"
os.chdir('../')
path = os.getcwd()
interpolate = 0
first_time = 1
modeOfFailure = "1 2 0 51"
# first position: mode of failure(1, 2, or 3); will work for only 1 and 3 no mode 2 in 2D
# second position: direction (1 for x and 2 for y)
# third position: minimum value in direction
# fourth position: maximum value in direction
# space is added between values to allow for use of split method and also to allow float values to be used in third and fourth position
DomainFilename=path+"/python_files/d1.txt"
FinalDisplacement= 0.0014
outputFile=path+'/inputd1.txt'
create_input_file.create_input_file(
                DomainFilename, FinalDisplacement, outputFile, eng, path + '/python_files', modeOfFailure)
