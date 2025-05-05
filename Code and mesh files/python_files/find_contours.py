import sys
import numpy as np
import cv2
from matplotlib import pyplot as plt
import find_where
import find_contours_funcs
import matlab.engine
import os

######################################################################################
# ------------------------------------------------------------------------------------
# Define the kernel size (depends on the image size) 
kernel_size     = 5 # 11 

# ------------------------------------------------------------------------------------
# Define the domain dimensions
min_x_in_domain = 0
min_y_in_domain = 0
max_x_in_domain = 25
max_y_in_domain = 51 

######################################################################################
# ------------------------------------------------------------------------------------
def find_contours(path, coords, num_cont, old_small_cnt, lc, weightx, weighty, first_time, elastic_inc, centroids, coordinates_redundant_border_stored, increment_py):

    repeat  = False
    output  = []  # Coordinate points in the damaged domain
    coords  = np.transpose(coords) # Stylistic preference, could have used it in its original form
    scale   = 2 # User defined scaling factor
    d2b     = 1 # Distance to boundary
    eng     = matlab.engine.connect_matlab()
    
    # ------------------------------------------------------------------------------------
    # Load the image
    img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_0_loading.jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 100])
    
    # ------------------------------------------------------------------------------------
    # Reverse the colormap
    h, w = img.shape[:2]
    img = abs(255 * np.ones([h,w]) - img)
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_1_reversing.jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 100])
    
    # ------------------------------------------------------------------------------------
    # Draw the border 
    pad_size = 100
    # print(img.shape)
    img      = cv2.copyMakeBorder(img, pad_size, pad_size, pad_size, pad_size, cv2.BORDER_CONSTANT, value=[np.median(img),np.median(img),np.median(img)])
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_2_making_border.jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 100]) 
    
    # ------------------------------------------------------------------------------------
    # While we are still in the elastic zone, identify the coordinates of the pixels that comprise the redundant border (MATLAB boundary, notches)
    if elastic_inc == 1:
        coordinates_redundant_border = np.where(img < [np.median(img)]) #  100
        elastic_inc = 0
    else:
        coordinates_redundant_border = coordinates_redundant_border_stored
    
    img[coordinates_redundant_border] = np.median(img)

    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_3_crack_removal.jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 100])
    
    # ------------------------------------------------------------------------------------
    # Apply thresholding
    _, img = cv2.threshold(img, np.median(img)-1, 255, cv2.THRESH_BINARY) # +cv2.THRESH_OTSU  np.median(img)-5
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_4_otsu_threshold.jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 100])

    # ------------------------------------------------------------------------------------
    # Apply morphological opening
    num_iterations  = 1 # 5
    kernel          = cv2.getStructuringElement(cv2.MORPH_RECT,(kernel_size,kernel_size))

    img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel, iterations=num_iterations)

    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_5_morphological_opening.jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 100])

    # ------------------------------------------------------------------------------------
    # Determine contour of all blobs found
    img = img.astype(np.uint8)
    contours, _ = cv2.findContours(img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    h, w = img.shape[:2]
    vis  = np.ones((h, w, 3), np.uint8) * 255
    cv2.drawContours(vis, contours, -1, (0, 0, 0), 3, cv2.LINE_AA)
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_6_contours1.jpg", vis, [cv2.IMWRITE_JPEG_QUALITY, 100])

    # Approximate countours as a closed polynomial
    contours = [cv2.approxPolyDP(cnt, 0.01*cv2.arcLength(cnt, True), True) for cnt in contours]
    
    # Sort the contours according to area
    contours = sorted(contours, key = cv2.contourArea)

    # Neglect the last contour because that is the boundary of the image
    contours = contours[:-1]

    # If you started n contours and suddenly you see more than n contours, n largest contours detected.
    if first_time == 0 and len(contours) > num_cont:
        contours = contours[-num_cont:]
        print("Hey look at if1")

    # If you started n contours and suddenly you see less than n contours, fill the missing contours with the smallest contour you see. 
    # This will be called especially in the double notch problem or any problem where there is more than 1 seed crack
    if first_time == 0 and len(contours) < num_cont:
        print("Hey look at if2")
        for i in range(num_cont-len(contours)):
            contours.insert(0,contours[0])

    # Display all polygon contours
    h, w = img.shape[:2]
    vis  = np.ones((h, w, 3), np.uint8) * 255
    cv2.drawContours(vis, contours, -1, (0, 0, 0), 3, cv2.LINE_AA)
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_6_contours2.jpg", vis, [cv2.IMWRITE_JPEG_QUALITY, 100])
    h, w = img.shape[:2]

    # Use the centroids to make sure that the sorted contours are in order of the previous sorted contours
    contours, new_centroid = find_contours_funcs.sort_contours_by_centroid(contours,first_time,centroids)
    # new_centroid = centroids
    
    # Make contours rectangular
    for i in range(len(contours)): # -1
        min_dir_y = min((contours[i])[:, 0, 1])
        max_dir_y = max((contours[i])[:, 0, 1])
        min_dir_x = min((contours[i])[:, 0, 0])
        max_dir_x = max((contours[i])[:, 0, 0])
        contours[i] = np.asarray([[[min_dir_x, min_dir_y]],
                                  [[max_dir_x, min_dir_y]],
                                  [[max_dir_x, max_dir_y]],

                                  [[min_dir_x, max_dir_y]]])

    # Display all rectangle contours
    vis  = np.ones((h, w, 3), np.uint8) * 255
    cv2.drawContours(vis, contours, -1, (0, 0, 0), 3, cv2.LINE_AA)        
    # cv2.imwrite("Inc_" + str(int(increment_py)) + "_After_7_contours.jpg", vis, [cv2.IMWRITE_JPEG_QUALITY, 100])
    h, w = img.shape[:2]

    x_min_p = pad_size
    y_min_p = pad_size
    x_max_p = w - pad_size - 1
    y_max_p = h - pad_size - 1

    # Obtain the real coordinate location of the contours
    real_pos_of_contour_points = []

    for i in range(len(contours)):
        real_pos_of_contour_points.append(find_contours_funcs.scale_to_real(np.asarray(contours[i]), x_min_p, y_min_p, x_max_p, y_max_p, min_x_in_domain, min_y_in_domain, max_x_in_domain, max_y_in_domain))
    contours_send_back = real_pos_of_contour_points.copy()

    ######################################################################################
    # If no NEW contours are detected in this step: 
    if(len(contours) == num_cont):

        # If no contours have been detected so far, then return 
        if(num_cont == 0):
            old_small_cnt  = old_small_cnt.tolist()
            new_decompbool = False
            
            print("Exit 1")
            # print("weightx: ", weightx.tolist())
            # print("weighty: ", weighty.tolist())            
            return repeat, output, np.asarray(old_small_cnt), new_decompbool, weightx.tolist(), weighty.tolist(), num_cont, first_time, elastic_inc, contours_send_back, new_centroid.tolist(), coordinates_redundant_border

        # If at least 1 contour has been detected so far, then check whether the new damage is close to boundary of the old damage domain
        else:
            for i in range(len(contours)):
                current_damage = np.squeeze(contours[i])
                current_damager = find_contours_funcs.scale_to_real(current_damage, x_min_p, y_min_p, x_max_p, y_max_p, min_x_in_domain, min_y_in_domain,max_x_in_domain, max_y_in_domain)
                dist = eng.min_dist_between_two_polygons(matlab.double(old_small_cnt[i].tolist()), matlab.double(current_damager.tolist()), matlab.double([weightx[i]]), matlab.double([weighty[i]]))
                # print("i: ", i)
                # print("current_damager: ", current_damager)
                # print("old_small_cnt[i].tolist(): ", old_small_cnt[i].tolist())
                # print("dist", dist)

                if(dist > d2b):
                    no_update_needed = True
                else:
                    print(dist)
                    if(dist == 0):
                        repeat = True
                    no_update_needed = False
                    break

            if no_update_needed:
                new_decompbool = False
            
                print("Exit 2")     
                # print("weightx: ", weightx.tolist())
                # print("weighty: ", weighty.tolist())            
                return repeat, output, np.asarray(old_small_cnt), new_decompbool, weightx.tolist(), weighty.tolist(), num_cont, first_time, elastic_inc, contours_send_back, new_centroid.tolist(), coordinates_redundant_border


                # if(dist > d2b + lc):  # 
                #     print("Check first if")
                #     new_decompbool = False
                #     print("Exit 2")

                #     return repeat, output, np.asarray(old_small_cnt), new_decompbool, weightx.tolist(), weighty.tolist(), num_cont, first_time, elastic_inc, contours_send_back, new_centroid.tolist(), coordinates_redundant_border

                # else:
                #     print("len(contours): ", len(contours))
                #     print("You need an update because: ")
                #     print("current_damager", current_damager)
                #     print("old_small_cnt", old_small_cnt[i])
                #     print(i)

                #     print(dist)
                #     print(d2b + lc) #  
                #     if(dist == 0):
                #         repeat = True
                #         print("YOU ARE HERE")
                #     break               

    ######################################################################################
    # If at least a new contour is detected in this step: 
    # OR 
    # If the existing contours are close to the interface:  

    for i in range(len(contours)):
        
        # Determine how the damaged regions are changing in the major direction
        curr_cont = np.squeeze(real_pos_of_contour_points[i])
        dx = max(curr_cont[:, 0])-min(curr_cont[:, 0])
        dy = max(curr_cont[:, 1])-min(curr_cont[:, 1])

        # If this is the first time damaged regions are detected, scale everything in both directions using the same weight 
        if first_time == 1:
            weightx[i] = 1
            weighty[i] = 1
            real_pos_of_contour_points[i] = find_contours_funcs.scale_contour_radial(real_pos_of_contour_points[i], d2b + 1, 3) # d2b + 1 
            if i == len(contours) - 1:
                first_time = 0

        # If this is not the first time damage has been detected, then look at how the damage is changing in the major directions and scale accordingly
        elif first_time == 0:

            # Growing more in the x direction
            if (dx/dy > dy/dx):
                scalex = scale
                scaley = max((dy/dx) * scale, 2)
                weightx[i] = 1
                weighty[i] = scalex/scaley

            # Growing more in the y direction
            elif(dx/dy < dy/dx):
                scalex = max((dx/dy) * scale, 2)
                scaley = scale
                weighty[i] = 1
                weightx[i] = scaley/scalex

            # Growing equally in both directions                
            else:
                scalex = scale
                scaley = scale
                weightx[i] = 1
                weighty[i] = 1

            real_pos_of_contour_points[i] = find_contours_funcs.scale_contour_radial(real_pos_of_contour_points[i], scaley, 2)
            real_pos_of_contour_points[i] = find_contours_funcs.scale_contour_radial(real_pos_of_contour_points[i], scalex, 1)

    # Check if the new scaled contours are intersecting
    for i in range(len(contours)):
        for j in range(i+1,len(contours)):
            current_scaled1 = np.squeeze(real_pos_of_contour_points[i])
            current_scaled2 = np.squeeze(real_pos_of_contour_points[j])
            dist = eng.min_dist_between_two_polygons(matlab.double(current_scaled1.tolist()), matlab.double(current_scaled2.tolist()), matlab.double([weightx[i]]), matlab.double([weighty[i]]))

            # If they are intersecting, then combine the intersecting domains into one domain
            if (dist == 0):
                min_dir_y = min(min((real_pos_of_contour_points[i])[:, 0, 1]),min((real_pos_of_contour_points[j])[:, 0, 1]))
                max_dir_y = max(max((real_pos_of_contour_points[i])[:, 0, 1]),max((real_pos_of_contour_points[j])[:, 0, 1]))
                min_dir_x = min(min((real_pos_of_contour_points[i])[:, 0, 0]), min((real_pos_of_contour_points[j])[:, 0, 0]))
                max_dir_x = max(max((real_pos_of_contour_points[i])[:, 0, 0]),max((real_pos_of_contour_points[j])[:, 0, 0]))
                real_pos_of_contour_points[i] = np.asarray([[[min_dir_x, min_dir_y]],
                                                [[max_dir_x, min_dir_y]],
                                                [[max_dir_x, max_dir_y]],
                                                [[min_dir_x, max_dir_y]]])
                real_pos_of_contour_points[j] = real_pos_of_contour_points[i]

    # Convert damaged region contours from pixel values to real values
    real_pos_of_contour_points2 = []

    for i in range(len(contours)):
        tmp_c = np.squeeze(real_pos_of_contour_points[i])
        real_pos_of_contour_points2.append(tmp_c)

    real_pos_of_contour_points = real_pos_of_contour_points2

    # Check which coordinates in the mesh are in the damaged region
    for j in range(len(coords)):
        for i in range(len(real_pos_of_contour_points)):
            if find_where.is_inside_polygon(real_pos_of_contour_points[i], tuple(coords[j, :])):
                output.append(j+1)
                break

    # print(real_pos_of_contour_points)
    old_small_cnt   = real_pos_of_contour_points
    num_cont        = len(old_small_cnt) # + 1
    old_small_cnt   = np.asarray(old_small_cnt)
    new_decompbool  = True

    print("Exit 3")
    # print("weightx: ", weightx.tolist())
    # print("weighty: ", weighty.tolist())            
    return repeat, output, old_small_cnt, new_decompbool, weightx.tolist(), weighty.tolist(), num_cont, first_time, elastic_inc, contours_send_back, new_centroid.tolist(), coordinates_redundant_border


# ------------------------------------------------------------------------------------
repeat, output, new_decomp, new_decompbool, weightx, weighty, num_of_conts, first_time, elastic_inc, contours_send_back, new_centroid, coordinates_redundant_border = find_contours(
    path, coords, prev_number_of_contours, old_small_cnt, LC, weightx, weighty, first_time, elastic_inc, centroids, coordinates_redundant_border_stored, increment_py)


