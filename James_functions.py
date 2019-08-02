import numpy as np
import matplotlib.pyplot as plt
import cv2
from matplotlib import cm
from scipy import spatial as sp
import os
import pandas as pd
from G_of_r import KNN_search_pbc


def g_6_over_g(filepath,
               distance_from_center,
               frame=0,
               r_cutoff=200,
               num_bins=100):
    """
    Calculate the ratio g_6(r)/g(r) for a dataframe and plot the histogram
    using the sqrt of the order parameter.
    PhysRevX.5.031025
    NOT FINISHED
    """
    dataframe = pd.read_hdf(filepath)
    # find central particles
    middle = find_central_particles(dataframe, distance_from_center)
    # Set dataframe index to particle
    dataframe = dataframe[(dataframe.frame == frame)].set_index('particle')
    # array of the particles that are in the middle for the right frame
    particle_list = middle.loc[middle['frame'] == frame, 'particle'].as_matrix()
    # data arrays to fill
    bins = np.linspace(0, r_cutoff, num_bins)
    g6 = np.zeros(len(bins) - 1)
    g = np.zeros(len(bins) - 1)
    # for each particle in the middle
    for i in particle_list:
        # distance to all other particles
        r = ((dataframe.get_value(i, 'x') - dataframe['x']) ** 2 +
             (dataframe.get_value(i, 'y') - dataframe['y']) ** 2) ** 0.5
        # sqrt of product of order parameters between starting particle and all other particles
        r6 = np.sqrt(dataframe.get_value(i, 'Order_param') *
                     dataframe['Order_param'])
        r = np.asarray(r)
        r6 = np.asarray(r6)
        for j in range(len(bins) - 1):
            # find which r values are in each bin
            in_bin = np.where((r > bins[j]) & (r <= bins[j + 1]))
            if np.shape(in_bin)[1] > 0:
                # add the population of this bin for this particle to g
                g[j] += np.size(r[in_bin])
                # add the sum of the products of orders for this bin for this particle to g6
                g6[j] += np.sum(r6[in_bin])

    # move bins values to middle of bins and divide by particle diameter
    bins = np.cumsum(np.diff(bins))
    bins = bins / (2 * dataframe['size'].mean())
    plt.figure()
    plt.loglog(bins, g6 / g, 'r-')
    plt.xlabel('r/d')
    plt.ylabel('g$_6$(r)/g(r)')


def find_central_particles(dataframe,
                           distance_from_center):
    """
    Returns a dataframe with just the particles that are in a square with
    x, y coordinates distance_from_center from the mean x,y positions
    """
    center = (dataframe['x'].mean(), dataframe['y'].mean())
    xmin = center[0] - distance_from_center
    xmax = center[0] + distance_from_center
    ymin = center[1] - distance_from_center
    ymax = center[1] + distance_from_center
    return dataframe.query('x>=@xmin and x<=@xmax and y>=@ymin and y<=@ymax')


def add_folder(path):
    "Adds the path if it doesn't exist"
    if not os.path.exists(path):
        os.makedirs(path)


def order_histogram(directory,
                    all_files=False,
                    file_number=0,
                    save_hist_data=False,
                    plot=False,
                    number_of_bins=100):
    """
    Calculates and plots the histogram of the order parameters for a dataframe
    using the sqrt of the order parameter

    Inputs:
        directory: directory holding the pandas dataframes
        all_files: If True process all dataframes in directory
        file_number: if all_files==False the dataframe to process - can be a list or int
        save_hist_data: if True save hist data to text file
        plot: if True plot stacked histograms
        number_of_bins: number of bins for the histogram
    """
    files = find_dataframes_in_directory(directory)
    # select the files and the indices for the files
    if not all_files:
        files = files[file_number]
        if type(file_number) is int:
            file_number = [file_number]
    else:
        file_number = np.arange(0, len(files), 1)
    # make sure directory ends with '/'
    if not directory.endswith("/"):
        directory = directory + '/'
    # make directory for data
    save_data_folder = directory + 'order_histogram_data/'
    add_folder(save_data_folder)

    volume_fraction = np.loadtxt(directory + 'volume_fractions.txt')

    i = -1
    for file in files:
        i += 1
        dataframe = pd.read_hdf(directory + file)
        order = dataframe.loc[dataframe['On_Edge'] == 1, 'Order_param'].as_matrix()  # exclude edge particles
        order = np.sqrt(order)
        freq, bins = produce_histogram_data(order, low_limit=0, upper_limit=1, number_of_bins=number_of_bins)
        if save_hist_data:
            vf = int(volume_fraction[file_number[i]] * 1000) / 1000  # force precision
            np.savetxt(save_data_folder + str(vf) + '_bins.txt', bins)
            np.savetxt(save_data_folder + str(vf) + '_freq.txt', freq)
        if plot:
            plt.figure()
            plt.plot(bins, freq, 'o')
            plt.xlabel('Order')
            plt.ylabel('Probability')


def kinetic_histogram(directory,
                      all_files=False,
                      file_number=0,
                      save_hist_data=False,
                      plot=False,
                      number_of_bins=100):
    """
    Calculates the kinetic histograms for the low and high ordered parts of
    the dataframe.

    Inputs:
        directory: directory holding the pandas dataframes
        all_files: If True process all dataframes in directory
        file_number: if all_files==False the dataframe to process - can be a list or int
        save_hist_data: if True save hist data to text file
        plot: if True plot stacked histograms
        number_of_bins: number of bins for the histogram
    """
    files = find_dataframes_in_directory(directory)
    # select the files and the indices for the files
    if not all_files:
        files = files[file_number]
        if type(file_number) is int:
            file_number = [file_number]
    else:
        file_number = np.arange(0, len(files), 1)
    # make sure directory ends with '/'
    if not directory.endswith("/"):
        directory = directory + '/'
    # make directory for data
    save_data_folder = directory + 'kinetic_histogram_data/'
    add_folder(save_data_folder)
    volume_fraction = np.loadtxt(directory + 'volume_fractions.txt')

    i = -1
    for file in files:
        i += 1
        dataframe = pd.read_hdf(directory + file)
        low = dataframe.query('Order_param<=0.5')  # dataframe with order <= 0.5
        high = dataframe.query('Order_param>0.5')  # dataframe with order > 0.5
        low_speed = low.loc[low['On_Edge'] == 1, 'Speed'].as_matrix()
        high_speed = high.loc[high['On_Edge'] == 1, 'Speed'].as_matrix()  # exclude edge particles
        low_energy = np.power(low_speed, 2)
        high_energy = np.power(high_speed, 2)
        freq_low, bins_low = produce_histogram_data(low_energy,
                                                    low_limit=min(low_energy),
                                                    upper_limit=max(low_energy),
                                                    number_of_bins=number_of_bins)
        freq_high, bins_high = produce_histogram_data(high_energy,
                                                      low_limit=min(high_energy),
                                                      upper_limit=max(high_energy),
                                                      number_of_bins=number_of_bins)
        if save_hist_data:
            vf = int(volume_fraction[file_number[i]] * 1000) / 1000
            np.savetxt(save_data_folder + str(vf) + "_low_bin.txt", bins_low)
            np.savetxt(save_data_folder + str(vf) + "_high_bin.txt", bins_high)
            np.savetxt(save_data_folder + str(vf) + "_low_freq.txt", freq_low)
            np.savetxt(save_data_folder + str(vf) + "_high_freq.txt", freq_high)
        if plot:
            plt.figure()
            plt.plot(bins_low, freq_low, 'o')
            plt.plot(bins_high, freq_high, 'o')
            plt.xlabel('Energy')
            plt.ylabel('Probability')
            plt.legend(['low', 'high'])


def density_histogram(directory,
                      all_files=False,
                      file_number=0,
                      save_hist_data=False,
                      plot=False,
                      number_of_bins=100):
    """
    Calculates the density histograms for the low and high ordered parts of
    the dataframe.

    Inputs:
        directory: directory holding the pandas dataframes
        all_files: If True process all dataframes in directory
        file_number: if all_files==False the dataframe to process - can be a list or int
        save_hist_data: if True save hist data to text file
        plot: if True plot stacked histograms
        number_of_bins: number of bins for the histogram
    """
    files = find_dataframes_in_directory(directory)
    # select the files and the indices for the files
    if not all_files:
        files = files[file_number]
        if type(file_number) is int:
            file_number = [file_number]
    else:
        file_number = np.arange(0, len(files), 1)
    # make sure directory ends with '/'
    if not directory.endswith("/"):
        directory = directory + '/'
    # make directory for data
    save_data_folder = directory + 'density_histogram_data/'
    add_folder(save_data_folder)
    volume_fraction = np.loadtxt(directory + 'volume_fractions.txt')

    i = -1
    for file in files:
        i += 1
        dataframe = pd.read_hdf(directory + file)
        low = dataframe.query('Order_param<=0.5')
        high = dataframe.query('Order_param>0.5')
        low_area = low.loc[low['On_Edge'] == 1, 'Voronoi_area'].as_matrix()
        high_area = high.loc[high['On_Edge'] == 1, 'Voronoi_area'].as_matrix()
        low_density = 1 / low_area
        high_density = 1 / high_area
        freq_low, bins_low = produce_histogram_data(low_density,
                                                    low_limit=min(low_density),
                                                    upper_limit=max(low_density),
                                                    number_of_bins=number_of_bins)
        freq_high, bins_high = produce_histogram_data(high_density,
                                                      low_limit=min(high_density),
                                                      upper_limit=max(high_density),
                                                      number_of_bins=number_of_bins)
        if save_hist_data:
            vf = int(volume_fraction[file_number[i]] * 1000) / 1000
            np.savetxt(save_data_folder + str(vf) + "_low_bin.txt", bins_low)
            np.savetxt(save_data_folder + str(vf) + "_high_bin.txt", bins_high)
            np.savetxt(save_data_folder + str(vf) + "_low_freq.txt", freq_low)
            np.savetxt(save_data_folder + str(vf) + "_high_freq.txt", freq_high)
        if plot:
            plt.figure()
            plt.plot(bins_low, freq_low, 'o')
            plt.plot(bins_high, freq_high, 'o')
            plt.xlabel('Density (pixels^-2)')
            plt.ylabel('Probability')
            plt.legend(['low', 'high'])


def produce_histogram_data(data,
                           low_limit=0,
                           upper_limit=440,
                           number_of_bins=100):
    freq, bin_edges = np.histogram(data, number_of_bins, [low_limit, upper_limit])
    bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    freq = freq.astype(float)
    freq = (freq / sum(freq))
    return freq, bins


def cluster_histogram(directory,
                      all_files=False,
                      file_number=0,
                      save_raw_data=False,
                      save_hist_data=False,
                      plot=False,
                      number_of_bins=100):
    """
    Calculates histograms for the number of clusters and the size of clusters.

    Inputs:
        directory: directory holding the pandas dataframes
        all_files: If True process all dataframes in directory
        file_number: if all_files==False the dataframe to process - can be a list or int
        save_hist_data: if True save hist data to text file
        plot: if True plot stacked histograms
        number_of_bins: number of bins for the histogram
    """
    files = find_dataframes_in_directory(directory)

    if not all_files:
        files = files[file_number]
        if type(file_number) is int:
            file_number = [file_number]
    else:
        file_number = np.arange(0, len(files), 1)
    if not directory.endswith("/"):
        directory = directory + '/'
    else:
        file_number = np.arange(0, len(files), 1)
    save_data_folder = directory + 'cluster_histogram_data/'
    add_folder(save_data_folder)
    volume_fraction = np.loadtxt(directory + 'volume_fractions.txt')
    i = -1
    for file in files:
        i += 1
        dataframe = pd.read_hdf(directory + file)
        number_of_frames = dataframe['frame'].max()
        cluster_size = []
        number_of_clusters = []
        for n in range(int(number_of_frames)):
            # get all clustesr for frame
            cluster = dataframe.loc[dataframe['frame'] == n,
                                    'Cluster'].as_matrix()
            # find the different clusters and their frequency
            cluster_id, cluster_frequency = np.unique(cluster,
                                                      return_counts=True)
            # excluding the 0th cluster (not in a cluster) append the size of each cluster
            # to the empty array
            if len(cluster_frequency) > 1:
                for item in cluster_frequency[1:]:
                    cluster_size.append(item)
            # append the number of clusters excluding cluster number 0
            number_of_clusters.append(len(cluster_id) - 1)
        cluster_size = np.ndarray.flatten(np.asarray(cluster_size))
        number_of_clusters = np.asarray(number_of_clusters)

        if save_raw_data:
            vf = int(volume_fraction[file_number[i]] * 1000) / 1000
            np.savetxt(save_data_folder + str(vf) + '_size_raw.txt', cluster_size)
            np.savetxt(save_data_folder + str(vf) + '_number_raw.txt', number_of_clusters)
        numbers_freq, numbers_bins = produce_histogram_data(number_of_clusters,
                                                            low_limit=min(number_of_clusters),
                                                            upper_limit=max(number_of_clusters),
                                                            number_of_bins=number_of_bins)
        size_freq, size_bins = produce_histogram_data(cluster_size,
                                                      low_limit=min(cluster_size),
                                                      upper_limit=max(cluster_size),
                                                      number_of_bins=number_of_bins)
        if save_hist_data:
            vf = int(volume_fraction[file_number[i]] * 1000) / 1000
            np.savetxt(save_data_folder + str(vf) + '_size_bins.txt', size_bins)
            np.savetxt(save_data_folder + str(vf) + '_size_freq.txt', size_freq)
            np.savetxt(save_data_folder + str(vf) + '_number_bins.txt', numbers_bins)
            np.savetxt(save_data_folder + str(vf) + '_number_freq.txt', numbers_freq)
        if plot:
            plt.figure()
            plt.plot(numbers_bins, numbers_freq, 'o')
            plt.xlabel('Number of clusters')
            plt.ylabel('Probability')
            plt.figure()
            plt.plot(size_bins, size_freq, 'o')
            plt.xlabel('Size of clusters')
            plt.ylabel('Probability')


def lever_rule_directory(directory,
                         volume_fraction_filename='volume_fractions.txt',
                         save_lever=False,
                         plot=False):
    """
    Calculates the area fraction of the solid phase to the total area as a
    function of the volume fraction

    Inputs:
        directory: string directory
        volume_fraction_filename: filename of text file containing
                                  volume_fractions
        save_lever: Boolean. If True save susceptibility as text file
        plot: Boolean. If True plot susceptibility against volume_fraction

    """
    files = find_dataframes_in_directory(directory)
    if not directory.endswith("/"):
        directory = directory + '/'
    crystal_fraction = []
    for file in files:
        print(len(crystal_fraction))
        dataframe = pd.read_hdf(directory + file)
        # dataframe with 'Voronoi_area' and 'On_Edge' for 'Order_param' > 0.5
        crystal_area = dataframe.loc[dataframe['Order_param'] > 0.5,
                                     ['Voronoi_area', 'On_Edge']]
        # total area for all the particles in a crystal that aren't on the edge
        total_crystal_area = crystal_area.loc[crystal_area['On_Edge'] == 1,
                                              'Voronoi_area'].sum()
        # total area for all particles that aren't on the edge
        all_area = dataframe.loc[dataframe['On_Edge'] == 1, 'Voronoi_area'].sum()
        crystal_fraction.append(total_crystal_area / all_area)
    if save_lever:
        np.savetxt(directory + 'crystal_fraction.txt', crystal_fraction)
    if plot:
        vol_fraction = np.loadtxt(directory + volume_fraction_filename)
        plt.figure()
        plt.plot(vol_fraction, crystal_fraction)
        plt.xlabel('$\phi$')
        plt.ylabel('$A_s/A$')


def susceptibility_directory(directory,
                             volume_fraction_filename='volume_fractions.txt',
                             save_susceptibility=False,
                             plot=False):
    """
    Calculates the susceptibility for all the dataframes in a directory

    Inputs:
        directory: string directory
        volume_fraction_filename: filename of text file containing
                                  volume_fractions
        save_susceptibility: Boolean. If True save susceptibility as text file
        plot: Boolean. If True plot susceptibility against volume_fraction
    """
    files = find_dataframes_in_directory(directory)
    if not directory.endswith("/"):
        directory = directory + '/'
    susceptibility = []
    for file in files:
        dataframe = pd.read_hdf(directory + file)
        # average order parameter over whole video
        order_space_time_mean = dataframe.loc[dataframe['On_Edge'] == 1,
                                              'Order_param'].mean()
        frames = max(dataframe['frame'])
        sus_array = np.zeros(int(frames))
        for f in range(int(frames)):
            # dataframe with 'Order_param' and 'On_Edge' columns for the frame
            order = dataframe.loc[dataframe['frame'] == f, ['Order_param', 'On_Edge']]
            # mean order for particle in this frame not on the edge
            mean_order = order.loc[order['On_Edge'] == 1, 'Order_param'].mean()
            # difference between the mean in this frame and the video mean
            sus_array[f] = mean_order - order_space_time_mean
        # square values for video
        sus_array = np.power(sus_array, 2)
        # append the mean of all the values for the video to get the susceptibility
        # value for a single volume fraction
        susceptibility.append(np.mean(sus_array))
    if save_susceptibility:
        np.savetxt(directory + "susceptibility.txt", susceptibility)
    if plot:
        vol_fraction = np.loadtxt(directory + volume_fraction_filename)
        plt.figure()
        plt.plot(vol_fraction, susceptibility)
        plt.xlabel('$\phi')
        plt.ylabel('$\chi_6$')


def phase_transition_directory(directory,
                               volume_fraction_filename='volume_fractions.txt',
                               save_order=False,
                               plot=False):
    """
    Finds the mean order parameter for each dataframe in a directory
    and plots it against the volume fraction for each dataframe.

    Inputs:
        directory: directory of dataframes with forward slashes
        volume_fraction_filename: filename of the text file with volume fractions
        save_order:  Boolean. If True save mean order in text file with filename
                     'mean_order.txt' in directory
        plot: Boolean. If True plot mean order against volume fraction
    """
    files = find_dataframes_in_directory(directory)
    if not directory.endswith("/"):
        directory = directory + '/'
    mean_order = []
    for file in files:
        dataframe = pd.read_hdf(directory + file)
        mean_order.append(dataframe['Order_param'].mean())
    if save_order:
        mean_order = np.asarray(mean_order)
        np.savetxt(directory + 'mean_order.txt', mean_order)
    if plot:
        plt.figure()
        vf = np.loadtxt(directory + volume_fraction_filename)
        plt.plot(vf, mean_order, 'o')
        plt.xlabel('$\phi$')
        plt.ylabel('<$\psi_6$>')


def find_dataframes_in_directory(directory):
    """
    Finds all the dataframes in a directory by looking for filename endings
    'hdf5'

    Inputs:
        directory: full directory with forward slashes
    Outputs:
        dataframe_files = list of filenames for the dataframes
    """
    files = os.listdir(directory)
    dataframe_files = []
    for names in files:
        if names.endswith("hdf5"):
            dataframe_files.append(names)
    return dataframe_files


def add_parameter_to_directory(directory, parameter):
    """
    Adds a parameter to all the dataframes in a directory
    Inputs:
        directory: string of the full directory. All backwards slashes should
                   be replaced with forwards slashes
        parameter: String heading of the dataframe column to be added.
    """
    dataframe_files = find_dataframes_in_directory(directory)
    if not directory.endswith("/"):
        directory = directory + '/'
    for file in dataframe_files:
        dataframe = pd.read_hdf(directory + file)
        if parameter == 'Speed':
            dataframe = add_speed(dataframe, wait_time=10)
        elif parameter == 'Order_param':
            dataframe = add_order_param(dataframe, window_size=5)
        elif parameter == 'x_mean' or parameter == 'y_mean':
            dataframe = add_rolling_mean(dataframe, window_size=5)
        elif parameter == 'Voronoi_area':
            dataframe = add_voronoi_cell_area(dataframe)
        elif parameter == 'Cluster':
            dataframe = add_cluster(dataframe, 'Order_param', 0.5)
        Jpt.save_panda(dataframe, directory + file[:-4])


def add_speed(dataframe, wait_time=1):
    """
    Adds the speed column to the dataframe calculated using the trackpy
    relate_frames function.

    Inputs:
        dataframe: pandas dataframe containing 'x', 'y', 'particle' and 'frame'
                   columns.
        wait_time: The number of frames to calculate the speed using

    Outputs:
        dataframe: pandas dataframe with additional 'Speed' column
    """
    dataframe = dataframe.sort_values(by='frame')
    speed_arr = np.asarray([])
    num_frames = dataframe['frame'].max()
    for n in range(int(num_frames) + 1):
        computed_movement_array = tp.relate_frames(dataframe, n, n + wait_time)
        speed_arr = np.append(speed_arr, computed_movement_array['dr'].as_matrix())
    dataframe['Speed'] = speed_arr
    return dataframe


def visualise_clusters(dataframe,
                       core_filename):
    """
    Overlays circles on top of the video with different colours corresponding
    to different cluster numbers

    Inputs:
        dataframe: pandas dataframe containing at least 'x', 'y', 'Cluster'
                   and 'frame' columns
        core_filename: full filepath of the video/dataframe without the
                       file extension
    """
    video_filename = core_filename + '.mp4'
    save_video_filename = core_filename + '_cluster.avi'
    crop_filename = core_filename + 'crop.txt'
    crop = np.loadtxt(crop_filename)
    cap = cv2.VideoCapture(video_filename)
    number_frames = cap.get(7)
    ret, frame = cap.read()
    fps = cap.get(5)
    if np.size(crop) == 4:
        frames_size = (crop[1, 1] - crop[1, 0], crop[0, 1] - crop[0, 0])
    else:
        frames_size = (int(len(frame)), crop[1] - crop[0])
    fourcc = cv2.VideoWriter_fourcc('X', 'V', 'I', 'D')
    write_obj = cv2.VideoWriter(
        save_video_filename,
        fourcc,
        int(fps),
        (int(frames_size[0]), int(frames_size[1])),
        True)
    for n in np.arange(0, number_frames - 1, 1):
        print(n)
        if n != 0:
            ret, frame = cap.read()
        if np.size(crop) == 4:
            frame = frame[crop[1, 1] - crop[1, 0], crop[0, 0]:crop[0, 1], :]
        else:
            frame = frame[:, crop[0]:crop[1], :]
        frame_current_x = dataframe[dataframe['frame'] == n, 'x'].as_matrix()
        frame_current_y = dataframe[dataframe['frame'] == n, 'y'].as_matrix()
        frame_current_parameter = dataframe[dataframe['frame'] == n,
                                            'Cluster'].as_matrix()
        max_param = np.max(frame_current_parameter)
        for m in np.arange(0, np.shape(frame_current_x)[0], 1):
            x = frame_current_x[m]
            y = frame_current_y[m]
            param = frame_current_parameter[m]
            param2_normalised = frame_current_parameter[m]
            param2_normalised = param2_normalised / max_param
            map_val = param2_normalised * 255
            col = cm.jet(np.floor(map_val).astype(int))[0:3]
            col_test = col[2], col[1], col[0]
            font = cv2.FONT_HERSHEY_SIMPLEX
            frame = cv2.circle(frame,
                               (x, y),
                               11,
                               np.floor(np.multiply(col_test, 255)),
                               -1)
            frame = cv2.putText(frame,
                                str(param),
                                (int(x), int(y)),
                                font,
                                1,
                                (255, 255, 255),
                                2)
        write_obj.write(frame)
    write_obj.release()


def first_nonzero(arr, axis, invalid_val=-1):
    mask = arr != 0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)


def find_bigNN(indx, dist, N=10000):

    index = indx[0::4]

    indx = np.ceil((indx + 1) / 4)
    indx_molecules = indx[0::4] - (np.arange(int(N / 4)).reshape(int(N / 4), 1) + 1)

    first_indx = first_nonzero(indx_molecules, axis=1, invalid_val=100)
    first_indx = np.vstack((np.arange(int(N / 4)), first_indx))
    dist_molecules = dist[0::4]
    param = dist_molecules[first_indx[0, :], first_indx[1, :]]

    indx_param = index[first_indx[0, :], first_indx[1, :]]

    indx_molecules[first_indx[0, :], first_indx[1, :]] = 0

    second_indx = first_nonzero(indx_molecules, axis=1, invalid_val=100)
    second_indx = np.vstack((np.arange(int(N / 4)), second_indx))
    dist_molecules_2 = dist[0::4]
    param_2 = dist_molecules_2[second_indx[0, :], second_indx[1, :]]

    indx_param2 = index[second_indx[0, :], second_indx[1, :]]


    parameters = np.vstack((param, param_2)).T
    parameters_indx = np.vstack((indx_param, indx_param2)).T

    return parameters, parameters_indx, indx


def add_param(dataframe, r, r2, N):

    number_frames = max(dataframe['frame'])
    print(number_frames)

    for f in range(int(number_frames + 1)):
        print(f)
        x_data = dataframe.loc[dataframe['frame'] == f, 'x']
        y_data = dataframe.loc[dataframe['frame'] == f, 'y']

        points = np.vstack((x_data, y_data))
        points = np.transpose(points)

    indx, dist = KNN_search_pbc(points[:, 0], points[:, 1], r, r2, N, srch=7)


    param, param_indx, molecule_indx = find_bigNN(indx, dist, N)
    param_indx = np.floor(param_indx/4)
    dataframe['molecule'] = molecule_indx[:, 0]

    array = np.arange(0, N)
    drop = np.append(array[1::4], array[2::4])
    drop = np.append(drop, array[3::4])
    dataframe = dataframe.drop(drop)

    dataframe['param_dist'] = param.tolist()
    dataframe['param_indx'] = param_indx.tolist()
    return dataframe




def add_cluster(dataframe,
                param,
                threshold,
                greater=True):
    """
    Finds clusters of particles defined under/over the threshold of a
    certain parameter by comparing parameter values for delaunay neighbours.

    Inputs:
        dataframe: panda dataframe including 'x', 'y', 'frame', param column
        param: The title of the dataframe column to define the clusters from
        threshold: The threshold value of param above/below which a
                   cluster is defined
        greater: Boolean. True = clusters defined above threshold.
                          False = clusters below

    Outputs:
    dataframe: pandas dataframe with added column 'Cluster'
    """
    number_frames = max(dataframe['frame'])
    dataframe['Cluster'] = 0
    print(number_frames)
    for f in range(int(number_frames)):
        print(f)
        x_data = dataframe.loc[dataframe['frame'] == f, 'x']
        y_data = dataframe.loc[dataframe['frame'] == f, 'y']
        parameter = dataframe.loc[dataframe['frame'] == f, param]

        points = np.vstack((x_data, y_data))
        points = np.transpose(points)

        tess = sp.Delaunay(points)


        neighbour_vertices = tess.vertex_neighbor_vertices
        list_indices = neighbour_vertices[0]
        point_indices = neighbour_vertices[1]
        """
        list_indices gives the initial and final indices of the
        point_indices list which gives the corresponding neighbours
        for each particle.
        """
        # Boolean array of points that are over the threshold
        if greater:
            parameter_boolean = parameter > threshold
        else:
            parameter_boolean = parameter <= threshold
        parameter_boolean = parameter_boolean.as_matrix()

        # Boolean array of the point_indices array for being over the threshold
        point_indices_bool = parameter_boolean[point_indices]

        # indices of points that are over the threshold
        clustered_points = np.argwhere(parameter_boolean)

        # check each of the points over the threshold
        cluster = np.zeros(len(points), int)
        cluster_no = 1
        for i in range(len(clustered_points)):
            cluster_no_store = cluster_no
            # Check to see if starting position is on a cluster
            if cluster[clustered_points[i]] == 0:
                cluster[clustered_points[i]] = cluster_no
            else:
                cluster_no = cluster[clustered_points[i]]

            # Find the neighbours that meet the threshold
            neighbours_all = point_indices[
                             list_indices[clustered_points[i]]:
                             list_indices[clustered_points[i] + 1]]
            neighbours_bool = point_indices_bool[
                              list_indices[clustered_points[i]]:
                              list_indices[clustered_points[i] + 1]]
            neighbours = neighbours_all[neighbours_bool]

            # Check which neighbours are closer than the cut-off distance
            neighbours_points = points[neighbours, :]
            cut_off_dist = 50

            point = points[clustered_points[i], :]
            point2neighbours_vector = neighbours_points - point
            point2neighbours_magnitude = np.sqrt(
                np.power(point2neighbours_vector[:, 0], 2) +
                np.power(point2neighbours_vector[:, 1], 2))
            indices_to_send = np.where(point2neighbours_magnitude <
                                       cut_off_dist)
            neighbours = neighbours[indices_to_send]

            # Check to see if any of the neighbours are already in a cluster
            neighbours_clusters = cluster[neighbours]

            # Find the equivalent clusters
            equivalent_cluster = neighbours_clusters[
                np.nonzero(neighbours_clusters)]
            equivalent_cluster = np.append(equivalent_cluster, cluster_no)

            # Set the neighbours cluster values to the current cluster
            cluster[neighbours] = cluster_no

            # Overwrite equivalent clusters
            lowest_cluster = np.min(equivalent_cluster)
            for j in range(len(equivalent_cluster)):
                cluster_no_temp = equivalent_cluster[j]
                ints = np.where(cluster == cluster_no_temp)
                cluster[ints] = lowest_cluster

            # Increase cluster value
            cluster_no = cluster_no_store + 1
        # Unique cluster values
        cluster_val_unique = np.unique(cluster)
        for i in range(len(cluster_val_unique)):
            unique_indices = np.nonzero(cluster == cluster_val_unique[i])
            if np.size(unique_indices, 1) == 1:
                cluster[unique_indices] = 0
        cluster_val_unique = np.unique(cluster)
        cluster_vals = np.arange(0, len(cluster_val_unique), 1)

        for k in range(len(cluster_vals)):
            temp = np.where(cluster == cluster_val_unique[k])
            cluster[temp] = cluster_vals[k]

        # Collect Panda data for first frame and add cluster number to the end
        dataframe.loc[dataframe['frame'] == f, 'Cluster'] = cluster

    return dataframe


def add_order_to_vid(dataframe,
                     core_filename):
    """

    This function uses the order paramter from the panda dataframe
    and adds circles of different colours depending on their order

    """
    parameter = 'Order_param'
    visualise_data(dataframe, core_filename, parameter)


def calculate_phase_transition_curve(dataframe,
                                     core_filename,
                                     shape='nut'):
    """
    Calculates and plots the average order parameter against area
    fraction for a single video where the barriers move
    """
    crop_filename = core_filename + 'crop.txt'
    number_of_frames = int(dataframe['frame'].max())
    area_fraction = np.zeros(number_of_frames)
    average_order = np.zeros(number_of_frames)
    crop = np.loadtxt(crop_filename)
    mask_lines = np.loadtxt(core_filename + '_mask_lines.txt')
    if shape == 'nut':
        particle_area = 0.5 * np.sqrt(3) * 4 ** 2

    for f in range(number_of_frames):
        print(f, ' of ', number_of_frames)
        order = dataframe.loc[dataframe['frame'] == f,
                              'Order_param'].as_matrix()
        pixel_width = mask_lines[0, f] - mask_lines[1, f]
        if np.size(crop) == 4:
            pixel_height = abs(crop[1, 0] - crop[1, 1])
        else:
            pixel_height = 1080
        height = 136  # mm
        width = (height / pixel_height) * pixel_width
        if f == 0:
            print("width = ", width)
        area_fraction[f] = particle_area * len(order) / (width * height)
        average_order[f] = np.mean(order)

    plt.figure()
    plt.plot(area_fraction, average_order, 'x')
    plt.xlabel('area fraction')
    plt.ylabel('<order parameter>')
    plt.savefig(core_filename + '_phase_transition.png')
    np.savetxt(core_filename + '_phase_transition.txt',
               (area_fraction, average_order))


def visualise_data(dataframe,
                   core_filename,
                   parameter):
    """
    This function takes a video and a dataframe and annotes the video with the
    data pertaining to the given parameter

    Inputs: P_data = Dataframe with x,y,frame and desired parameter columns
            VidFilename = Path to input video_file
           saveVid_filename = Path for saving annotated video
            crop = crop values used to crop video when initially processed or
                   produced using crop_circle of crop_hexagon classes
           parameter = header of parameter column in submitted Dataframe to
                        use to annotate Video

    """
    video_filename = core_filename + '.mp4'
    crop_filename = core_filename + 'crop.txt'
    save_vid_filename = core_filename + '_order.avi'
    cap = cv2.VideoCapture(video_filename)
    crop = np.loadtxt(crop_filename)
    number_frames = cap.get(7)
    ret, frame = cap.read()
    if np.size(crop) == 4:
        frame_size = (int(crop[1, 1] - crop[1, 0]), int(crop[0, 1] - crop[0, 0]))
    else:
        frame_size = (int(len(frame)), int(crop[1] - crop[0]))
    fourcc = cv2.VideoWriter_fourcc('X', 'V', 'I', 'D')
    write_obj = cv2.VideoWriter(save_vid_filename,
                                fourcc,
                                30.0,
                                (frame_size[0], frame_size[1]),
                                True)

    for n in np.arange(0, number_frames, 1):
        print(n)
        if n != 0:
            ret, frame = cap.read()
        if np.size(crop) == 4:
            frame = frame[crop[1, 0]:crop[0, 1], crop[0, 0]:crop[0, 1], :]
        else:
            frame = frame[:, crop[0]:crop[1], :]

        frame_currx = dataframe.loc[dataframe['frame'] == n,
                                    'x'].as_matrix()
        frame_curry = dataframe.loc[dataframe['frame'] == n,
                                    'y'].as_matrix()
        frame_currparam = dataframe.loc[dataframe['frame'] == n,
                                        parameter].as_matrix()

        pix = 14
        pixperm = pix / 8e-3
        invfps = 1 / 60
        for m in np.arange(0, np.shape(frame_currx)[0], 1):
            x = frame_currx[m]
            y = frame_curry[m]
            param = frame_currparam[m]
            if parameter == 'Order_param':
                map_val = param * 255
                col = cm.jet(np.floor(map_val).astype(int))[0:3]
                col_test = col[2], col[1], col[0]
            if parameter == 'Lindemann_param':
                param = param / np.max(dataframe[parameter])
                map_val = param * 255
                col = cm.jet(np.floor(map_val).astype(int))[0:3]
                col_test = col[2], col[1], col[0]
            if parameter == 'r_velocity':
                param = ((param / pixperm) / invfps) * 1e-3
                print(param)
                if param < -0.0003:
                    col_test = [0, 0, 0]
                else:
                    col_test = [0, 0, 1]

            frame = cv2.circle(frame,
                               (x, y),
                               11,
                               np.floor(np.multiply(col_test, 255)),
                               -1)
        if n == 0:
            plt.imshow(frame)
        write_obj.write(frame)

    write_obj.release()


def add_rolling_mean(dataframe,
                     window_size=5):
    """
    Adds the rolling 'x_mean' and 'y_mean' columns to the dataframe based on
    window_size
    """
    particles = dataframe['particle']
    particles = particles.unique()
    togo = np.shape(particles)
    i = 0
    for m in particles:
        print(i, " of ", togo)
        i += 1
        dataframe.loc[dataframe['particle'] == m, 'x_mean'] = \
            dataframe.loc[dataframe['particle'] == m, 'x'].rolling(min_periods=1,
                                                                   window=window_size,
                                                                   center=False).mean()
        dataframe.loc[dataframe['particle'] == m, 'y_mean'] = \
            dataframe.loc[dataframe['particle'] == m, 'y'].rolling(min_periods=1,
                                                                   window=window_size,
                                                                   center=False).mean()
    return dataframe


def add_order_param(dataframe,
                    window_size=5):
    """

    This function calculates the order parameter for each particle in each
    frame and adds a column to the inputed dataframe with the calculated values

    Inputs: p_data = Dataframe with at least x,y and frame columns
            window_size = number of frames to average x and y over to improve
                          noise cause by particles rattling in cages
    Outputs p_data = Dataframe with Order_param column added

    """

    dataframe = dataframe.sort_values(by='frame')
    dataframe_columns = dataframe.columns.tolist()
    if "x_mean" not in dataframe_columns:
        dataframe = add_rolling_mean(dataframe)
    num_frames = np.max(dataframe['frame'])
    data_points = dataframe['x'].shape
    order_param_arr = np.zeros(data_points)

    count = 0
    separations = np.array([])
    for n in np.arange(0, num_frames + 1, 1):
        print(n)
        frame_x = dataframe.loc[dataframe['frame'] == n, 'x_mean']
        frame_y = dataframe.loc[dataframe['frame'] == n, 'y_mean']

        points = np.vstack((frame_x, frame_y))
        points = np.transpose(points)

        tess = sp.Delaunay(points)

        neighbour_vertices = tess.vertex_neighbor_vertices

        list_indices = neighbour_vertices[0]
        # For each point it gives you the initial
        # and final indices of the point indices list which
        # correspond to its nearest neighbours
        point_indices = neighbour_vertices[1]
        # An organised list of point indices
        # which when sliced using list indices correspond to the
        # nearest neighbour of each point sequentially

        bead_num = len(frame_x)

        if n == 0:
            for m in np.arange(0, bead_num, 1):
                point = points[m, :]

                neighbour_indices = point_indices[list_indices[m]:list_indices[m + 1]]
                neighbour_points = points[neighbour_indices, :]

                point2neighbours_vector = neighbour_points - point
                point2neighbours_magnitude = np.sqrt(
                    np.power(point2neighbours_vector[:, 0], 2) +
                    np.power(point2neighbours_vector[:, 1], 2))
                separations = np.append(separations,
                                        point2neighbours_magnitude)

        cut_off_dist = np.mean(separations) * 10
        for m in np.arange(0, bead_num, 1):
            point = points[m, :]

            neighbour_indices = point_indices[list_indices[m]:list_indices[m + 1]]
            neighbour_points = points[neighbour_indices, :]

            point2neighbours_vector = neighbour_points - point
            point2neighbours_magnitude = np.sqrt(
                np.power(point2neighbours_vector[:, 0], 2) +
                np.power(point2neighbours_vector[:, 1], 2))
            indices_to_send = np.where(point2neighbours_magnitude <
                                       cut_off_dist)

            theta_arr = angle_between(point2neighbours_vector[indices_to_send,
                                      :])

            order_param_exp = np.cos(6 * theta_arr)
            order_param = np.sum(order_param_exp) / len(theta_arr)
            order_param_arr[count] = np.power(order_param, 2)

            count = count + 1

    dataframe['Order_param'] = order_param_arr

    return dataframe


def angle_between(p2):
    """
    This function takes a list of a co-ordinates and computes the angles
    between them.

    Inputs: p2 = numpy array of cordinates (x and y co-ordinates listed)
    Outputs: diffs = numpy array of angular separations
    """

    ang2 = np.arctan2(p2[:, :, 0], p2[:, :, 1])
    sorted_arr = np.transpose(np.sort(-ang2, 0))
    sorted_arr = np.append(sorted_arr, (sorted_arr[0] + (2 * np.pi)))
    diffs = np.diff(sorted_arr)
    return diffs


def density_video(dataframe,
                  core_filename,
                  particle_type,
                  window_size=10):
    video_filename = core_filename + '.mp4'
    crop_filename = core_filename + 'crop.txt'
    save_video_filename = core_filename + '_density.avi'
    cap = cv2.VideoCapture(video_filename)
    number_of_frames = cap.get(7)
    crop = np.loadtxt(crop_filename)
    size = np.mean(dataframe.loc[dataframe['frame'] == 0,
                                 'size'].as_matrix())

    # Calculate the area of a particle
    if particle_type == 'nut':
        particle_area = 0.5 * 6 * size ** 2 * np.sin(2 * np.pi / 6)
    else:
        particle_area = np.pi * size ** 2

    # Crop the frame and write the new video file
    ret, frame = cap.read()
    if np.size(crop) == 4:
        frame_cropped = frame[:, crop[0, 0]:crop[0, 1], [2, 1, 0]]
    else:
        frame_cropped = frame[:, crop[0]:crop[1], [2, 1, 0]]

    fourcc = cv2.VideoWriter_fourcc('X', 'V', 'I', 'D')
    write_obj = cv2.VideoWriter(save_video_filename,
                                fourcc,
                                30.0,
                                (np.shape(frame_cropped)[1],
                                 np.shape(frame_cropped)[0]),
                                True)
    w = window_size // 2

    # Calculate density for each frame and add to video
    for n in range(int(number_of_frames) - 1):
        print(n, " of ", number_of_frames)

        # crop the frame if not the first frame
        if n != 0:
            ret, frame = cap.read()
            if np.size(crop) == 4:
                frame_cropped = frame[:, crop[0, 0]:crop[0, 1], [2, 1, 0]]
            else:
                frame_cropped = frame[:, crop[0]:crop[1], [2, 1, 0]]
        frame_cropped = frame_cropped.copy()

        x_pos = dataframe.loc[dataframe['frame'] == n,
                              'x'].as_matrix()
        voronoi = dataframe.loc[dataframe['frame'] == n,
                                'Voronoi_area'].as_matrix()
        xmin = np.min(x_pos)
        xmax = np.max(x_pos)

        # Calculate the average density along x for windows in x
        x_plot_points = np.arange(int(xmin), int(xmax), 1)
        for x in x_plot_points:
            in_window = np.argwhere((x_pos >= x - w) & (x_pos < x + w))
            average_area = np.nanmedian(voronoi[in_window])
            density = ((particle_area / average_area) *
                       np.shape(frame_cropped)[0])
            if density < 1080 and not np.isnan(density):
                frame_cropped = cv2.circle(frame_cropped,
                                           (x, int(density)),
                                           5, (0, 255, 0), -1)
        write_obj.write(frame_cropped)
    write_obj.release()


def add_voronoi_cell_area(dataframe,
                          scale=1):
    """
    This function calculates the voronoi cell area for each particle in
    each frame and adds a column to the inputed dataframe with the calculated
    values

    Inputs: P_data = Dataframe with at least x,y, xmean, ymean and frame
            window_size = number of frames to average x and y over to improve
                          noise cause by particles rattling in cages
    Outputs P_data = Dataframe with Voronoi_area column added
    """

    dataframe = dataframe.sort_values(by='frame')
    dataframe_columns = dataframe.columns.tolist()
    if "x_mean" not in dataframe_columns:
        dataframe = add_rolling_mean(dataframe)
    num_frames = np.max(dataframe['frame'])

    voronoi_area_param_arr = []  # empty array for voronoi area data
    particle_on_edge = []
    for n in np.arange(0, num_frames + 1, 1):
        print(n)

        # actual points from dataframe
        frame_currx = dataframe.loc[dataframe['frame'] == n,
                                    'x_mean'].as_matrix()
        frame_curry = dataframe.loc[dataframe['frame'] == n,
                                    'y_mean'].as_matrix()
        points = np.vstack((frame_currx, frame_curry))
        points = np.transpose(points)

        # vertices for the square of extra boundary points, first repeated
        corners_x = [min(frame_currx) - 20,
                     max(frame_currx) + 20,
                     max(frame_currx) + 20,
                     min(frame_currx) - 20,
                     min(frame_currx) - 20]
        corners_y = [min(frame_curry) - 20,
                     min(frame_curry) - 20,
                     max(frame_curry) + 20,
                     max(frame_curry) + 20,
                     min(frame_curry) - 20]

        # Add boundary points
        x_edge = []
        y_edge = []
        for m in range(4):
            x_edge.append(np.linspace(corners_x[m], corners_x[m + 1], 50))
            y_edge.append(np.linspace(corners_y[m], corners_y[m + 1], 50))
        x_edge = np.ndarray.flatten(np.array(x_edge))
        y_edge = np.ndarray.flatten(np.array(y_edge))

        # Combine actual and boundary points
        all_x = np.concatenate((frame_currx, np.array(x_edge)))
        all_y = np.concatenate((frame_curry, np.array(y_edge)))
        all_points = np.vstack((all_x, all_y))
        all_points = all_points.transpose()

        # spatial algorithms
        vor = sp.Voronoi(all_points)
        tri = sp.Delaunay(all_points)

        # Calculates the voronoi areas of the actual points by taking
        # the area of the polygon given by the vertices for each point
        points_region = vor.point_region
        regions = vor.regions
        region_vertices = vor.vertices
        for p in np.arange(0, len(points), 1):
            region_indice = points_region[p]
            vertices_indices = regions[region_indice]
            points_region_vertices = region_vertices[vertices_indices]
            voronoi_area_param_arr.append(polygon_area(points_region_vertices))

        # Check the actual points to see if any of their neighbours
        # given by the delaunay tesselation are one of the artifically
        # added boundary points.
        # If true then point is on the edge and it gets a value 1.
        # on_edge = np.zeros(len(points))
        delaunay_neighbours = tri.vertex_neighbor_vertices
        list_indices = delaunay_neighbours[0]
        point_indices = delaunay_neighbours[1]
        for m in range(len(points)):
            on_edge = 0
            neighbours = point_indices[list_indices[m]:list_indices[m + 1]]
            edge = neighbours > len(points)
            if np.any(edge):
                on_edge = 1
            particle_on_edge.append(on_edge)
        print("particle on edge = ", np.shape(particle_on_edge))
        # Check results
        if n == 100:
            plt.figure()
            plt.plot(frame_currx,
                     frame_curry,
                     'rx')
            plt.plot(x_edge,
                     y_edge,
                     'bo')
            edge_indices = np.argwhere(on_edge == 1)
            plt.plot(frame_currx[edge_indices],
                     frame_curry[edge_indices],
                     'go')
            sp.voronoi_plot_2d(vor)
            sp.delaunay_plot_2d(tri)

    # Add dataframe columns
    voronoi_area_param_arr = np.asarray(voronoi_area_param_arr)
    print("particle on edge's shape is ", np.shape(particle_on_edge))
    dataframe['Voronoi_area'] = voronoi_area_param_arr * scale
    dataframe['Edge_particle'] = particle_on_edge

    return dataframe


def polygon_area(corners):
    """
    This function calculates the area of a polygon using the shoelace formula

    Inputs: corners = list of a polygons corner co-ordinates MUST BE ORDERED
                      in either CLOCKWISE or COUNTER-CLOCKWISE
    Outputs area = Polygon area
    """

    n = len(corners)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area
