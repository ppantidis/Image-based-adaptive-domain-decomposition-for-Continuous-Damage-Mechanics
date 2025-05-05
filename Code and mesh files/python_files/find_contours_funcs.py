import numpy as np
import cv2

######################################################################################
# ------------------------------------------------------------------------------------
def calc_contour_area(contour):
    # Calculate the contour area
    area = cv2.contourArea(contour)
    
    return area

# ------------------------------------------------------------------------------------
def calc_contour_centroid(contour):
    num = contour.shape[0]
    cx = sum(contour[:,0,0])/num
    cy = sum(contour[:,0,1])/num
    
    return cx, cy

# ------------------------------------------------------------------------------------
def scale_contour_radial(cnt, scale, i):
    cx, cy = calc_contour_centroid(cnt)
    cent=np.asarray([cx,cy]).reshape((1,1,2))
    if(i == 1):
        # scale in x direction
        scaled_cnt = cnt-cent
        scaled_cnt[:, 0, 0] = scaled_cnt[:, 0, 0]*scale
        scaled_cnt = scaled_cnt+cent

    elif(i == 2):
        # scale in y direction
        scaled_cnt = cnt-cent
        scaled_cnt[:, 0, 1] = scaled_cnt[:, 0, 1]*scale
        scaled_cnt = scaled_cnt+cent

    elif(i == 3):
        # scale in radial direction
        scaled_cnt = cnt-cent
        scaled_cnt = scaled_cnt*scale
        scaled_cnt = scaled_cnt+cent

    return scaled_cnt

# ------------------------------------------------------------------------------------
def scale_to_real(X, x_min_p, y_min_p, x_max_p, y_max_p, x_min_r, y_min_r, x_max_r, y_max_r):
    if len(X.shape)==2:
        X2Return = [[0 for i in range(X.shape[1])] for j in range(X.shape[0])]
        X2Return = np.asarray(X2Return,dtype=np.float32)
        for i in range(X.shape[0]):
            x_p = X[i][0]
            y_p = X[i][1]
            X2Return[i][0] = ((x_p-x_max_p)*(x_min_r)/(x_min_p-x_max_p)) + \
                ((x_p-x_min_p)*(x_max_r)/(x_max_p-x_min_p))
            X2Return[i][1] = ((y_p-y_max_p)*(y_max_r)/(y_min_p-y_max_p)) + \
                ((y_p-y_min_p)*(y_min_r)/(y_max_p-y_min_p))
        X2Return = np.asarray(X2Return)
        
    if len(X.shape)==3:
        X2Return = [[[0 for i in range(X.shape[2])] for j in range(X.shape[1])]for k in range(X.shape[0])]
        X2Return = np.asarray(X2Return,dtype=np.float32)
        
        for i in range(X.shape[0]):
            x_p = X[i][0][0]
            y_p = X[i][0][1]
            X2Return[i][0][0] = ((x_p-x_max_p)*(x_min_r)/(x_min_p-x_max_p)) + \
                ((x_p-x_min_p)*(x_max_r)/(x_max_p-x_min_p))
            X2Return[i][0][1] = ((y_p-y_max_p)*(y_max_r)/(y_min_p-y_max_p)) + \
                ((y_p-y_min_p)*(y_min_r)/(y_max_p-y_min_p))
        X2Return = np.asarray(X2Return)

    return X2Return

# ------------------------------------------------------------------------------------
def scale_to_fake(X, x_min_p, y_min_p, x_max_p, y_max_p, x_min_r, y_min_r, x_max_r, y_max_r):
    X2Return = [[0 for i in range(X.shape[1])] for j in range(X.shape[0])]
    for i in range(X.shape[0]):
        X2Return = np.asarray(X2Return,dtype=np.float32)
        x_r = X[i][0]
        y_r = X[i][1]
        X2Return[i][0] = ((x_r-x_max_r)*(x_min_p)/(x_min_r-x_max_r)) + \
            ((x_r-x_min_r)*(x_max_p)/(x_max_r-x_min_r))
        X2Return[i][1] = ((y_r-y_max_r)*(y_max_p)/(y_min_r-y_max_r)) + \
            ((y_r-y_min_r)*(y_min_p)/(y_max_r-y_min_r))
    #X2Return = np.asarray(X2Return)
    X2Return = np.round(X2Return, decimals=0)

    return X2Return

# ------------------------------------------------------------------------------------
def sort_contours_by_centroid(contours,first_time,centroids):
        
    new_centroid = np.zeros((len(contours),2))

    for i in range(len(contours)):
        new_centroid[i,:] = calc_contour_centroid(contours[i])

    if not first_time:
        Big_I = np.zeros(len(new_centroid),dtype=np.int32)

        for i in range (len(new_centroid)):
            distances = np.zeros(len(centroids))
            for j in range(len(centroids)):
                distances[j]=np.linalg.norm(new_centroid[i,:]-centroids[j,:])
            Big_I[i] = np.argmin(distances)
            
        sorted_contours = []
        sorted_centroid = []

        for i in range(len(contours)):
            sorted_contours.append(contours[Big_I[i]])
            sorted_centroid.append(new_centroid[Big_I[i]])
        contours     = sorted_contours
        new_centroid = np.asarray(sorted_centroid)

    return contours, new_centroid




