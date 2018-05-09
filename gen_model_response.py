import sys
import random
import numpy as np
import os
import subprocess
import uuid
import shutil
#generates a set of models - one row per model as blob_modelling_parameters
#and corresponding responses - one row per responses

## Michael King

# now add a line of parameters and responses.
parameterFile='./Responses/parameters1.txt'

# open up the files.
f2=open(parameterFile,'a');

for i in range(0,5000):

    #generating unique ID for current blob
    identifier = str(uuid.uuid4())

    #blob parameters
    
    
    background = random.uniform(0.4,0.6)
    #intensity = 0.3
    #strength = 0.9
    #attenuation = 0.9
    #background = random.uniform(0,1)
    intensity = random.uniform(0,1)
    strength = random.uniform(0.01,0.99)
    attenuation = random.uniform(0.02,0.99)
    
    #mesh is 250,000m in x and y direction - 50,000 of which is densley meshed in 16x3125m
    #mesh is 150,000m in z direction - 100,000 of which is meshed at 420*i^1.09088674 for i=0; i<19 (for 20 layers(
    #max zpos is 0.4: max(deltaX,deltaY,deltaZ) = 250km. deltaZ = 150km where last 50km is padding
    #therefore 100km/250km = 0.4
    
    xpos = random.uniform(0.4, 0.6)
    ypos = random.uniform(0.4, 0.6)
    xrad = random.uniform(0.05,0.2)
    yrad = random.uniform(0.05,0.2)
    theta = 0.0
    zpos = (random.uniform(0, 1))**1.2 * 0.4
    zrot = theta
    zrad = random.uniform(0.05,0.2)
    xrot = 0.0
    yrot = 0.0
    '''
    #blob 2
    intensity2 = random.uniform(0,1)
    strength2 = random.uniform(0.01,0.99)
    attenuation2 = random.uniform(0.02,0.99)
    xpos2 = random.uniform(0.4, 0.6)
    ypos2 = random.uniform(0.4, 0.6)
    xrad2 = random.uniform(0.05,0.2)
    yrad2 = random.uniform(0.05,0.2)
    theta2 = 0.0
    zpos2 = (random.uniform(0, 1))**1.2 * 0.4
    zrad2 = random.uniform(0.05,0.2)
    xrot2 = 0.0
    yrot2 = 0.0
    
    #blob 3
    intensity3 = random.uniform(0,1)
    strength3 = random.uniform(0.01,0.99)
    attenuation3 = random.uniform(0.02,0.99)
    xpos3 = random.uniform(0.4, 0.6)
    ypos3 = random.uniform(0.4, 0.6)
    xrad3 = random.uniform(0.05,0.2)
    yrad3 = random.uniform(0.05,0.2)
    theta3 = 0.0
    zpos3 = (random.uniform(0, 1))**1.2 * 0.4
    zrad3 = random.uniform(0.05,0.2)
    xrot3 = 0.0
    yrot3 = 0.0
    
    #blob 4
    intensity4 = random.uniform(0,1)
    strength4 = random.uniform(0.01,0.99)
    attenuation4 = random.uniform(0.02,0.99)
    xpos4 = random.uniform(0.4, 0.6)
    ypos4 = random.uniform(0.4, 0.6)
    xrad4 = random.uniform(0.05,0.2)
    yrad4 = random.uniform(0.05,0.2)
    theta4 = 0.0
    zpos4 = (random.uniform(0, 1))**1.2 * 0.4
    zrad4 = random.uniform(0.05,0.2)
    xrot4 = 0.0
    yrot4 = 0.0
'''
    #b=[0.5, 1.0, 0.9, 0.9, 0.9, 0.5, 0.15, 0.15, 0.0, 0.1,0.15, 0.0, 0]
    
    #combining all 4 blobs into one array
    b=[background, intensity, strength, attenuation, xpos, ypos, xrad, yrad, theta, zpos, zrad, xrot, yrot]
    ''' intensity2, strength2, attenuation2, xpos2, ypos2, xrad2, yrad2, theta2, zpos2, zrad2, xrot2, yrot2, intensity3, strength3, attenuation3, xpos3, ypos3, xrad3, yrad3, theta3, zpos3, zrad3, xrot3, yrot3, intensity4, strength4, attenuation4, xpos4, ypos4, xrad4, yrad4, theta4, zpos4, zrad4, xrot4, yrot4]'''
    
    #b1=[x*0.95 for x in b]
    c=' '.join(map(str,b))
    print(c)

    model_template = './Startup_Files/20x20x20.model'
    target_model = './Current_Model/targ.model'

    #Taking blob parameters and generating model
    callString='./saveModel '+ model_template + ' ' + c + ' > ' + target_model
    os.system(callString);


    if '-nan' not in open(target_model).read():
        print "ok"

        p1 = subprocess.Popen(['./wsinv3dmt'])
        p1.wait()

        #storing response and model files with unique identifier
        response='./Current_Model/blk3d_resp.02';
        shutil.copy2(response, "./Responses/"+identifier)

        #callString = 'python ws2vtk.py ./Models/' + identifier + ' ./Responses/' + identifier
        #os.system(callString)

        #storing a copy VTK of model for viewing (has unique identifier)
        #shutil.copy2('./VTKResistivityGrid.vtr', './Grids/' + identifier+'.vtr')

        # now write the parameters and the response
        f2.write(identifier + ': ' + c+'\n')
        print "done"

f2.close()



