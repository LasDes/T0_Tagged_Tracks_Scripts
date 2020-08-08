# Write code for plotting Y vs. Z and Y vs. Corrected X Plots of the Anode- And Cathode-Piercing Cosmic Ray Muon Events
import ROOT
import numpy as np
# import matplotlib.pyplot as plt			
import root_numpy 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.ticker
import scipy.interpolate
import math
import offset_correction_algo

# Set the numpy print options so that entire arrays will be printed
#np.set_printoptions(threshold=np.inf)

f = ROOT.TFile('50000_event_gallery_output_mcc83_track_trajectory_information.root')
tree = f.Get('SCEtree')

# Define a variable for the number of events in this file (in float form, so that it divides properly)
num_of_events = 50000.0

# Make the branches in the TTree into an array, and don't impose the cut that the track must enter or exit the cathode JUST YET.  It will make looping through that variable later on much easier.
arr = root_numpy.tree2array(tree, branches = ['matched', 'event_runNum', 'event_subrunNum', 'event_eventNum', 'mctrack_t0', 'track_t0', 'time_flash', 'pe_flash', 'track_length', 'track_startX', 'track_startY', 'track_startZ', 'track_endX', 'track_endY', 'track_endZ', 'anode', 'cathode', 'enters_top', 'enters_front', 'enters_back', 'exits_bottom', 'exits_front', 'exits_back'], selection = 'enters_top == 1 && anode == 1')
df = pd.DataFrame(arr)

print "The number of ALL data tracks in this sample = %d." % df.shape[0]
print "\n"

# Define vectors for the three plots that you are making here: 
x_start     = []
x_end       = []

y_start     = []
y_end       = []

z_start     = []
z_end       = []

y_coord_avg = []
z_coord_avg = []

# Define coordinates here for the coordinates of the track and their average (for organization)
track_x_start     = 0.0
track_x_end       = 0.0
track_y_start     = 0.0
track_y_end       = 0.0
track_z_start     = 0.0
track_z_end       = 0.0
track_y_coord_avg = 0.0
track_z_coord_avg = 0.0

# Start a loop over the 'df'
track_iterator = 0

# Define variables for the increments that you'll be using to make the plots below
dx = 12.0
dy = 12.0
dz = 48.0

# Define variables for min_x (the minimum tick on the x-axis), max_x (the maximum tick on the x-axis), min_y (the minimum tick on the y-axis), max_y (the maximum tick on the y-axis), min_z (the minimum tick on the z-axis), and max_z (the maximum tick on the z-axis)
min_x = -25.0
max_x = 275.0

min_y = -150.0
max_y = 150.0

min_z = -100.0
max_z = 1100.0

# Define the offsets that you'll have to use for the anode case and the cathode case (I'll call them 'opposite_offset' because it's the opposite of what was placed in the algorithm)
# I do not need any offset, because none is provided in the producer.
anode_opposite_offset   = 0.0
cathode_opposite_offset = 0.0
drift_velocity          = 0.1114

# Begin the code for the grid of the yz-plot

# Define numpy grids for the two grids in y and z
# y - y axis
# z - x axis 
z_grid_for_yz_plot, y_grid_for_yz_plot = np.mgrid[slice(min_z, max_z + dz, dz), slice(min_y, max_y + dy, dy)]

print "The outer dimension of the z_grid_for_yz_plot is %d and its inner dimension is %d." %(len(z_grid_for_yz_plot), len(z_grid_for_yz_plot[0]))
print "The outer dimension of the y_grid_for_yz_plot is %d and its inner dimension is %d." %(len(y_grid_for_yz_plot), len(y_grid_for_yz_plot[0]))
print "\n"

# Print out these two grids
print "The grid in z that I'm using for the yz plot = ", z_grid_for_yz_plot
print "\n"

print "The grid in y that I'm using for the yz plot = ", y_grid_for_yz_plot
print "\n"

# Define a numpy array for the output points in z that will be divided by the number of events in the end
# This numpy array should have a dimension of (len(y_grid_for_yz_plot) - 1, len(y_grid_for_yz_plot[0] - 1))
number_of_cosmics_passing_through_total_yz_plot = np.zeros([len(y_grid_for_yz_plot) - 1, len(y_grid_for_yz_plot[0]) - 1])

# Define a variable for the number of coordinates out of range in the yz plot
num_of_coords_out_of_range_yz_plot = 0

# Begin the code for the grid of the yx-plot

# Define numpy grids for the two grids in y and x
# y - y axis
# x - x axis

x_grid_for_yx_plot, y_grid_for_yx_plot = np.mgrid[slice(min_x, max_x + dx, dx), slice(min_y, max_y + dy, dy)]

print "The outer dimension of the x_grid_for_yx_plot is %d and its inner dimension is %d." %(len(x_grid_for_yx_plot), len(x_grid_for_yx_plot[0]))
print "The outer dimension of the y_grid_for_yx_plot is %d and its inner dimension is %d." %(len(y_grid_for_yx_plot), len(y_grid_for_yx_plot[0]))
print "\n"

# Print out the two grids
print "The grid in x that I'm using for the yx plot = ", x_grid_for_yx_plot
print "\n"

print "The grid in y that I'm using for the yx plot = ", y_grid_for_yx_plot
print "\n"

# Define a numpy array for the output points in z that will be divided by the number of events in the end                                                                      
# This numpy array should have a dimension of (len(x_grid_for_yx_plot) - 1, len(x_grid_for_yx_plot[0] - 1))                                                                 
number_of_cosmics_passing_through_total_yx_plot = np.zeros([len(y_grid_for_yx_plot) - 1, len(y_grid_for_yx_plot[0]) - 1])

# Define a variable for the number of coordinates out of range in the yx plot                                                       
num_of_coords_out_of_range_yx_plot = 0

# Define number of x, y, and z coordinates out of range (I think all of the coordinates out of range may be y-coordinates)
x_coords_out_of_range = 0
y_coords_out_of_range = 0
z_coords_out_of_range = 0

while track_iterator < df.shape[0]:

    print "The current track number that we're on in this loop is: ", track_iterator, "."
    print "\n"

    # VERY IMPORTANT: Impose that the track either enters or exits the cathode, which are the tracks that we are looking at in this case
    # This cut should have already been applied, but keep it in there so you do not have to change the indentation in the loops.
    if (df['anode'][track_iterator] == 1):

        # Keep this in here so you do not have to change the indentation in the loops.
        # Only look at those tracks with a 2.0 us difference between the reconstructed t0 and the flash time.
        if ( df['track_t0'][track_iterator] - df['time_flash'][track_iterator] ) > -2.0 or ( df['track_t0'][track_iterator] - df['time_flash'][track_iterator] ) < 2.0:

            # Set the track's starting y-coordinate, the starting z-coordinate, and the starting z-coordinate, and the ending z-coordinate equal to the values of the df at this entry
            track_x_start = df['track_startX'][track_iterator]
            track_x_end   = df['track_endX'][track_iterator]
            track_y_start = df['track_startY'][track_iterator]
            track_y_end   = df['track_endY'][track_iterator]
            track_z_start = df['track_startZ'][track_iterator]
            track_z_end   = df['track_endZ'][track_iterator]

            # Define the reconstructed T0 as 'rc_time' so that it can be used in the calculations for adjusted x-coordinates below
            rc_time = df['track_t0'][track_iterator]

            # Correct 'track_x_start' by the value of the 'rc_time' times 0.1114 cm/us (defined above as 'drift_velocity') (call this value 'track_x_start_rctime')                                 
            track_x_start_rctime = track_x_start - (drift_velocity*rc_time)
            track_x_end_rctime = track_x_end - (drift_velocity*rc_time)

            # Calculate the 'track_y_coord_avg' and the 'track_z_coord_avg' that we're looking at 
            track_y_coord_avg = (track_y_start + track_y_end)/2.0
            track_z_coord_avg = (track_z_start + track_z_start)/2.0

            # Append all of these values to their respective lists
            x_start.append(track_x_start_rctime)
            x_end.append(track_x_end_rctime)
            y_start.append(track_y_start)
            y_end.append(track_y_end)
            z_start.append(track_z_start)
            z_end.append(track_z_end)
            y_coord_avg.append(track_y_coord_avg)
            z_coord_avg.append(track_z_coord_avg)

            print "track_x_start = %f." % track_x_start_rctime
            print "track_x_end = %f." % track_x_end_rctime
            print "track_y_start = %f." % track_y_start
            print "track_y_end = %f." % track_y_end
            print "track_z_start = %f." % track_z_start
            print "track_z_end = %f." % track_z_end
            print "\n"

            # Define the linear functions for the 'yz' and 'yx' cases given the starting and ending points below
            yz_function = scipy.interpolate.interp1d([track_y_start, track_y_end], [track_z_start, track_z_end], kind='linear')
            yx_function = scipy.interpolate.interp1d([track_y_start, track_y_end], [track_x_start_rctime, track_x_end_rctime], kind='linear')
            
            # Find out which of the coordinates are greater
            # Find which values of x, y, and z are greater than the other                                                                                                     
            y_lesser  = track_y_end
            y_greater = track_y_start

            x_corresponding_to_y_lesser  = track_x_end_rctime
            x_corresponding_to_y_greater = track_x_start_rctime

            z_corresponding_to_y_lesser  = track_z_end
            z_corresponding_to_y_greater = track_z_start

            # Switch them if this is not the case                                                                                                                        
            if track_y_start < track_y_end:

                y_lesser  = track_y_start
                y_greater = track_y_end

                x_corresponding_to_y_lesser  = track_x_start_rctime
                x_corresponding_to_y_greater = track_x_end_rctime
        
                z_corresponding_to_y_lesser  = track_z_start
                z_corresponding_to_y_greater = track_z_end

            # Define a variable for the difference between the two points in z (this will be a positive number, based on the code immediately above)
            y_diff = y_greater - y_lesser
                
            # Start iterating from the first point on the line, 'y_lesser'
            y_iterator = y_lesser

            # Define a boolean for if one of the y-coordinates are out of range
            y_coord_out_of_range = False

            # NEW LOOP FOR Y COORDINATES THAT ARE OUT OF RANGE:
            # The only coordinates I'm going to throw out are those for which the 'y_lesser' coordinate is greater than the maximimum y coordinate in the plot and the 'y_greater' coordinate is less than the minimum y-coordinate in the plot.  These are the tracks that don't enter into the image that the plot covers at all.
            if y_lesser > max_y or y_greater < min_y:

                # Set 'y_coord_out_of_range' equal to 'True'
                y_coord_out_of_range = True

                num_of_coords_out_of_range_yz_plot += 1
                num_of_coords_out_of_range_yx_plot += 1
                
                y_coords_out_of_range += 1

            # Define a boolean for if one of the z-coordinates are out of range
            z_coord_out_of_range = False

            # Handle the case in which one of the track's z-coordinates are out of range
            if z_corresponding_to_y_greater > max_z or z_corresponding_to_y_greater < min_z or z_corresponding_to_y_lesser > max_z or z_corresponding_to_y_lesser < min_z:

                # Only increment 'num_of_coords_out_of_range_yz_plot' if 'y_coord_out_of_range' == False; otherwise, that variable has already been incremented for this track.
                if y_coord_out_of_range == False:

                    num_of_coords_out_of_range_yz_plot += 1

                z_coord_out_of_range = True

                # Increment the number of z-coordinates out of range                                                                                       
                z_coords_out_of_range += 1

            # Define a boolean for if one of the x-coordinates are out of range                                                                                  
            x_coord_out_of_range = False

            # Handle the case in which one of the track's x-coordinates are out of range
            if x_corresponding_to_y_greater > max_x or x_corresponding_to_y_greater < min_x or x_corresponding_to_y_lesser > max_x or x_corresponding_to_y_lesser < min_x: 

                # Only increment 'num_of_coords_out_of_range_yx_plot' if 'y_coord_out_of_range' == False; otherwise, that variable has already been incremented for this track.
                if y_coord_out_of_range == False:

                    num_of_coords_out_of_range_yx_plot += 1

                x_coord_out_of_range = True

                # Increment the number of x-coordinates out of range                                                                                                            
                x_coords_out_of_range += 1

            # Only enter this outer loop if the y-coordinate is within range.  Otherwise, you can skip this track (it will not hurt our statistics.)
            if y_coord_out_of_range == False:

                # Begin the loop over 'y_iterator', which is the independent variable in each of the functions above
                while y_iterator < y_greater:

                    # Only enter into the loop for the yz function value if the z-coordinate is within range
                    if z_coord_out_of_range == False:

                        # Include a statement for if the track starts outside the image and leaves the image/starts inside the image and leaves the image
                        if y_iterator < min_y or y_iterator > max_y:

                            # If the loop passes this statement, then increment y_iterator and continue until the loop is terminated/the loop starts (if the track has yet to enter the image)
                            y_iterator += dy
                            continue

                        current_z_func_value = yz_function(y_iterator)

                        print "The current y value on the function line = %f." % y_iterator
                        print "The current z value retrieved from the function line = %f." % current_z_func_value
                    
                        # Find out which indices in the array 'number_of_cosmics_passing_through_total_yz_plot' should be filled
                        y_index_heatmap = int((y_iterator - min_y) / dy) 
                        z_index_heatmap = int((current_z_func_value - min_z) / dz)

                        print "The y-index on the heatmap for this y value = %d." % y_index_heatmap
                        print "The z-index on the heatmap for this z value = %d." % z_index_heatmap
                        print "\n"

                        # Increment 'number_of_cosmics_passing_through_total_yz_plot' at this location in the array
                        number_of_cosmics_passing_through_total_yz_plot[z_index_heatmap][y_index_heatmap] += 1

                    # Only enter into the loop for the xy function value if the z-coordinate is within range
                    if x_coord_out_of_range == False:

                        # Include a statement for if the track starts outside the image and leaves the image/starts inside the image and leaves the image                            
                        if y_iterator < min_y or y_iterator > max_y:

                            # If the loop passes this statement, then increment y_iterator and continue until the loop is terminated/the loop starts (if the track has yet to enter the image)                                                                                                                                     
                            y_iterator += dy
                            continue

                        current_x_func_value = yx_function(y_iterator)

                        # I do not even need to apply this offset.
                        # Correct 'current_x_func_value' using 'offset_correction_algo' to account for the offsets present in the Producer
                        # current_x_func_value = offset_correction_algo.offset_correction_algo(anode_opposite_offset, cathode_opposite_offset, current_x_func_value, drift_velocity)

                        print "The current y value on the function line = %f." % y_iterator
                        print "The current x value retrieved from the function line = %f." % current_x_func_value

                        # Find out which indices in the array 'number_of_cosmics_passing_through_total_yx_plot' should be filled
                        y_index_heatmap = int((y_iterator - min_y) / dy)
                        x_index_heatmap = int((current_x_func_value - min_x) / dx)

                        print "The y-index on the heatmap for this y value = %d." % y_index_heatmap
                        print "The x-index on the heatmap for this x value = %d." % x_index_heatmap
                        print "\n"

                        # Increment 'number_of_cosmics_passing_through_total_yx_plot'
                        number_of_cosmics_passing_through_total_yx_plot[x_index_heatmap][y_index_heatmap] += 1

                    # Increment the y_iterator
                    y_iterator += dy

                # Now, consider the cases on the end of the lines after 'y_iterator' has exceeded 'y_greater'

                # Include a final statement at the end of the 'yz' track for the possibility that the endpoints weren't reached.  This happens if the algorithm "overshoots" the endpoint in the final box, meaning that that box is not counted.

                # Don't assign 'y_index_heatmap' a new value in the first loop, because then you cannot use it for comparison to enter the second loop.

                # Include the 'if' statement for 'z_coord_out_of_range' first (this is located within the loop of if the y-coordinate is within range as well)
                if z_coord_out_of_range == False:

                    # Make sure that the ending y-value is within the boundaries of the image, as well
                    if (y_index_heatmap != int((y_greater - min_y) / dy) or z_index_heatmap != int((z_corresponding_to_y_greater - min_z) / dz)) and y_greater < max_y:

                        # Print out the indices for these ending values
                        print "The y-index on the heatmap for this y value = %d." % int((y_greater - min_y) / dy)
                        print "The z-index on the heatmap for this z value = %d." % int((z_corresponding_to_y_greater - min_z) / dz)

                        print "The greater value of y (among the 'start' and 'end' coordinates) = %f." % y_greater
                        print "The value of z corresponding to the greater value of y (among the 'start' and 'end' coordinates) = %f." % z_corresponding_to_y_greater
                        print "\n"

                        # In this case, increment the indices of the last point on the particle's track (they should stay within the graph's limits; those were set very liberally)
                        number_of_cosmics_passing_through_total_yz_plot[int((z_corresponding_to_y_greater - min_z) / dz)][int((y_greater - min_y) / dy)] += 1

                # Include the 'if' statement for 'x_coord_out_of_range' first (this is located within the loop of if the y-coordinate is within range as well)
                if x_coord_out_of_range == False:

                    # Include a final statement at the end of the 'yx' track for the possibility that the endpoints weren't reached.  This happens if the algorithm "overshoots" the endpoint in the final box, meaning that that box is not counted. 
                    # Make sure that the ending y-value is within the boundaries of the image, as well
                    if (y_index_heatmap != int((y_greater - min_y) / dy) or x_index_heatmap != int((x_corresponding_to_y_greater - min_x) / dx)) and y_greater < max_y:

                        # Correct 'x_corresponding_to_y_greater' using 'offset_correction_algo' to account for the offsets present in the Producer                             
                        x_corresponding_to_y_greater = offset_correction_algo.offset_correction_algo(anode_opposite_offset, cathode_opposite_offset, x_corresponding_to_y_greater, drift_velocity)
                        
                        # Print out the indices for those ending values 
                        print "The y-index for the heatmap for this y value = %d." % int((y_greater - min_y) / dy)
                        print "The x-index on the heatmap for this x value = %d." % int((x_corresponding_to_y_greater - min_x)/ dx)

                        print "The greater value of y (among the 'start' and 'end' coordinates) = %f." % y_greater
                        print "The value of x corresponding to the greater value of y (among the 'start' and 'end' coordinates) = %f." % x_corresponding_to_y_greater
                        print "\n"
                                       
                        number_of_cosmics_passing_through_total_yx_plot[int((x_corresponding_to_y_greater - min_x)/ dx)][int((y_greater - min_y) / dy)] += 1

    # Increment 'track_iterator'
    track_iterator += 1

# BEGIN GENERATING THE YZ PLOTS (The Log Color Scale and the Linear Color Scale)

# Define a new numpy array, 'number_of_cosmics_passing_through_per_event', having the same dimensions as the original array
number_of_cosmics_passing_through_per_event_yz_plot = np.zeros([len(z_grid_for_yz_plot) - 1, len(z_grid_for_yz_plot[0]) - 1])

# Fill this new array with a normalized version of the other array
for outer_index in range(len(y_grid_for_yz_plot) - 1):

    for inner_index in range(len(y_grid_for_yz_plot[0]) - 1):

        number_of_cosmics_passing_through_per_event_yz_plot[outer_index][inner_index] = float(number_of_cosmics_passing_through_total_yz_plot[outer_index][inner_index]) / num_of_events

# Print out the dimensions of this list just to ensure that it is not too large
print "The outer dimension of 'number_of_cosmics_passing_through_per_event_yz_plot' = %d and the inner dimension = %d." %(len(number_of_cosmics_passing_through_per_event_yz_plot), len(number_of_cosmics_passing_through_per_event_yz_plot[0]))
print "\n"

# Begin code to plot Y vs. Z on the log scale

# Pick the desired colormap (this will be used for the linear scale plot as well)
cmap = plt.get_cmap('Reds')

# Declare the first figure and its first plot.  This is the log plot of Y (y-axis) vs. Z (x-axis) with a log color scale.                                               
fig1, ax1 = plt.subplots(nrows=1)

# Plot Y vs. Z on the log scale on the figure just declared
im1 = ax1.pcolormesh(z_grid_for_yz_plot, y_grid_for_yz_plot, number_of_cosmics_passing_through_per_event_yz_plot, cmap=cmap, norm=matplotlib.colors.LogNorm(1.0e-4, 2.0e-1, clip=True))
fig1.colorbar(im1, ax=ax1)
ax1.set_title('MCC8.1 Top/Anode Tracks: Track-Hit Density Per Event With Additional Cuts', fontsize=9, y = 1.02)
ax1.set_xlim([-100, 1100])
ax1.set_ylim([-150, 150])

# Change the size of the font for the axes labels                                                                                                                                
ax1.xaxis.label.set_fontsize(15)
ax1.yaxis.label.set_fontsize(15)

# Set the ticks on the y- and z-axes
plt.xticks(np.arange(-100.0, 1300.0, 200.0)) 
plt.yticks(np.arange(-150.0, 175.0, 25.0))

# Label the two axes
plt.xlabel('Track $z$ Coordinate [cm]')
plt.ylabel('Track $y$ Coordinate [cm]')

# Begin code to plot Y vs. Z on the linear scale

# Define the color levels for the linear track density heatmap                                                                                                                
yz_track_density_levels = matplotlib.ticker.MaxNLocator(nbins=100).tick_values(number_of_cosmics_passing_through_per_event_yz_plot.min(), number_of_cosmics_passing_through_per_event_yz_plot.max())
yz_track_density_norm = matplotlib.colors.BoundaryNorm(yz_track_density_levels, ncolors = cmap.N, clip=True)

# Declare the second figure and its second plot.  This is the linear plot of Y (y-axis) vs. Z (x-axis) with a linear color scale.
fig2, ax2 = plt.subplots(nrows=1)

# Plot Y vs. Z on the linear scale on the figure just declared
im2 = ax2.pcolormesh(z_grid_for_yz_plot, y_grid_for_yz_plot, number_of_cosmics_passing_through_per_event_yz_plot, cmap=cmap, norm = yz_track_density_norm)# matplotlib.colors.SymLogNorm(linthresh = 0.002, linscale = 0.28, vmin = 0.0, vmax = 0.043)) 
yz_lin_colorbar = fig2.colorbar(im2, ax=ax2)
ax2.set_title('MCC8.1 Top/Anode Tracks: Track-Hit Density Per Event With Additional Cuts', fontsize=9, y = 1.02)
ax2.set_xlim([-100, 1100])
ax2.set_ylim([-150, 150])

# yz_lin_colorbar.set_ticks([0.001, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036, 0.043])
# yz_lin_colorbar.set_ticklabels(["{0:.3f}".format(0.001), "{0:.3f}".format(0.004), "{0:.3f}".format(0.008), "{0:.3f}".format(0.012), "{0:.3f}".format(0.016), "{0:.3f}".format(0.020), "{0:.3f}".format(0.024), "{0:.3f}".format(0.028), "{0:.3f}".format(0.032), "{0:.3f}".format(0.036), "{0:.3f}".format(0.043)])

# Change the size of the font for the axes labels                                                                                                                               
ax2.xaxis.label.set_fontsize(15)
ax2.yaxis.label.set_fontsize(15)

# Set the ticks on the y- and z-axes                                                                                                                                           
plt.xticks(np.arange(-100.0, 1300.0, 200.0))
plt.yticks(np.arange(-150.0, 175.0, 25.0))

# Label the two axes
plt.xlabel('Track $z$ Coordinate [cm]')
plt.ylabel('Track $y$ Coordinate [cm]')

# Print the number of coordinates out of range in the yz plot
print "The number of coordinates out of range in the yz plot  = %d." % num_of_coords_out_of_range_yz_plot

# BEGIN GENERATING THE YX PLOTS (The Log Color Scale and the Linear Color Scale)

# Define a new numpy array, 'number_of_cosmics_passing_through_per_event', having the same dimensions as the original array                                                   
number_of_cosmics_passing_through_per_event_yx_plot = np.zeros([len(x_grid_for_yx_plot) - 1, len(x_grid_for_yx_plot[0]) - 1])

# Fill this new array with a normalized version of the other array                                                                                                            
for outer_index in range(len(y_grid_for_yx_plot) - 1):

    for inner_index in range(len(y_grid_for_yx_plot[0]) - 1):

        number_of_cosmics_passing_through_per_event_yx_plot[outer_index][inner_index] = float(number_of_cosmics_passing_through_total_yx_plot[outer_index][inner_index]) / num_of_events

# Begin code to plot Y vs. X on the log scale

# Pick the desired colormap (this will be used for the linear scale plot as well, and I redefine it here just for bookkeeping purposes)                                         
cmap = plt.get_cmap('Reds')

# Declare the third figure and its third plot. This is the log plot of Y (y-axis) vs. X (x-axis) with a log color scale.                                                          
fig3, ax3 = plt.subplots(nrows=1)

 # Plot Y vs. X on the log scale on the figure just declared
im3 = ax3.pcolormesh(x_grid_for_yx_plot, y_grid_for_yx_plot, number_of_cosmics_passing_through_per_event_yx_plot, cmap=cmap, norm=matplotlib.colors.LogNorm(1.0e-4, 2.0e-1, clip=True))
fig3.colorbar(im3, ax=ax3)
ax3.set_title('MCC8.1 Top/Anode Tracks: Track-Hit Density Per Event With Additional Cuts', fontsize=9, y = 1.02)
ax3.set_xlim([-50, 325])
ax3.set_ylim([-150, 150])

# Change the size of the font for the axes labels                                                                                                                                
ax3.xaxis.label.set_fontsize(15)
ax3.yaxis.label.set_fontsize(15)

# Set the ticks on the y- and z-axes                                                                                                                                              
plt.xticks(np.arange(-50.0, 350.0, 50.0)) # Here, the last entry will be at x = 300.0 cm, but the plot will go up to 325.0 cm (the upper limit that I considered)
plt.yticks(np.arange(-150.0, 175.0, 25.0))

# Label the two axes
plt.xlabel('Corrected Track $x$ Coordinate [cm]')
plt.ylabel('Track $y$ Coordinate [cm]')

print "The maximum value in the y vs. x heatmap on the cathode-piercing side = %f." % number_of_cosmics_passing_through_per_event_yx_plot.max()
print "\n"

# Begin code to plot Y vs. X on the linear scale

# Define the color levels for the purity heatmap                                                                                                                               
yx_track_density_levels = matplotlib.ticker.MaxNLocator(nbins=100).tick_values(number_of_cosmics_passing_through_per_event_yx_plot.min(), number_of_cosmics_passing_through_per_event_yx_plot.max())
yx_track_density_norm = matplotlib.colors.BoundaryNorm(yx_track_density_levels, ncolors = cmap.N, clip=True)

# Declare the fourth figure and its fourth plot. This is the linear plot of Y (y-axis) vs. X (x-axis) with a linear color scale.                                                
fig4, ax4 = plt.subplots(nrows=1)

# Plot Y vs. X on the linear scale on the figure just declared 
im4 = ax4.pcolormesh(x_grid_for_yx_plot, y_grid_for_yx_plot, number_of_cosmics_passing_through_per_event_yx_plot, cmap=cmap, norm=yx_track_density_norm) # matplotlib.colors.SymLogNorm(linthresh = 0.00335, linscale = 0.040, vmin = 0.0, vmax = 0.067))
yx_lin_colorbar = fig4.colorbar(im4, ax=ax4)
ax4.set_title('MCC8.1 Top/Anode Tracks: Track-Hit Density Per Event With Additional Cuts', fontsize=9, y = 1.02)
ax4.set_xlim([-25, 275])
ax4.set_ylim([-150, 150])

# Set the tick values on the colorbar
# yx_lin_colorbar.set_ticks([0.001, 0.005, 0.010, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.067])
# yx_lin_colorbar.set_ticklabels(["{0:.3f}".format(0.001), "{0:.3f}".format(0.00500), "{0:.3f}".format(0.010), "{0:.3f}".format(0.020), "{0:.3f}".format(0.025), "{0:.3f}".format(0.030), "{0:.3f}".format(0.035), "{0:.3f}".format(0.040), "{0:.3f}".format(0.045), "{0:.3f}".format(0.050), "{0:.3f}".format(0.055), "{0:.3f}".format(0.060), "{0:.3f}".format(0.067)])

# Change the size of the font for the axes labels                                                                                                                              
ax4.xaxis.label.set_fontsize(15)
ax4.yaxis.label.set_fontsize(15)

# Set the ticks on the y- and z-axes                                                                                                                                           
plt.xticks(np.arange(-25.0, 325.0, 50.0)) # Here, the last entry will be at x = 300.0 cm, but the plot will go up to 325.0 cm (the upper limit that I considered)              
plt.yticks(np.arange(-150.0, 175.0, 25.0))

# Label the two axes
plt.xlabel('Corrected Track $x$ Coordinate [cm]')
plt.ylabel('Track $y$ Coordinate [cm]')

# Print the number of coordinates out of range in the yx plot                                                                                                                    
print "The number of coordinates out of range in the yx plot  = %d." % num_of_coords_out_of_range_yx_plot

# Print the number of coordinates of each type out of range
print "The number of x-coordinates out of range = %d." % x_coords_out_of_range 
print "The number of y-coordinates out of range = %d." % y_coords_out_of_range
print "The number of z-coordinates out of range = %d." % z_coords_out_of_range

plt.show()
