import numpy as np
import numbers
import scipy.signal as s

def cyclic_correlate(a: np.ndarray, v: np.ndarray):
    '''
    description:
        calculate the cyclic auto correlation funtion.
        The length of the output will be the same as the length of the input

    inputs:
        a                -  The first function to be correlated. The lenth of
                            this function is the same as the lenth of the output
                            function
        v                -  The second function to be correlated.

    returns:
        corr_function    -  The cyclically calculated auto correlation function

    '''
    a_size = np.shape(a)[0]
    begin_of_seq = int(a_size/2)
    end_of_seq = begin_of_seq + a_size-1
    begin_of_constructed_array = begin_of_seq-int(len(v)/2)
    end_of_constructed_array = end_of_seq+int(len(v)/2)+len(v)%2
    constructed_array = np.tile(np.fft.fftshift(a), 2)\
        [begin_of_constructed_array:end_of_constructed_array]
    return np.correlate(constructed_array, v, mode = "valid")

def values_above(sequence: np.ndarray, threshold: numbers.Number):
    '''
    descrpition:
        Find the longest consecutive chain of values above a threshold
        in a sequence, containing the middle element

    inputs:
        sequence          - The array of datapoints
        threshold         - The absolute threshold, the sequence is
                            tested against

    returns:
        max_val           - the longest consecutive chain of ones
                            in the sequence
    '''
    #reorder array and mark values above the threhold as 1 and below as 0
    bool_sequence = np.fft.fftshift(sequence > threshold)
    count_poitive = 0
    count_negative = 0
    middle = int(len(sequence)/2)
    #count values to the right of the maximum value
    for idx, s in enumerate(bool_sequence):
        if s <= threshold or idx == middle:
            count_poitive = idx
            break
    #flip the array
    bool_sequence = np.flip(bool_sequence)
    #count values to the left of the maximum value
    for idx, s in enumerate(bool_sequence):
        if s <= threshold or idx == len(sequence)-middle:
            count_negative = idx
            break
    #return the comined counting sequences
    return count_poitive+count_negative

class Rx_power:
    def __init__(self, power_dBFS: np.ndarray, time: np.ndarray):
        '''
        description:
            Constructor of Rx_power class
        inputs:
            power_dBFS       - A sequence of power level
            time             - The time axis that belong to the power sequence
        '''

        self.power_dBFS = power_dBFS
        self.time = time

        self.average = None
        self.min = None
        self.max = None

    def set_avg_min_max(self):
        '''
        description:
            Find average, min and max vaues of the sequence of power levels
        '''
        self.average = 10* np.log10(np.mean(10**(self.power_dBFS/10)))
        self.min = np.min(self.power_dBFS)
        self.max = np.max(self.power_dBFS)

class Impulse_response:
    def __init__(self, h: np.ndarray, time: np.ndarray, delay : np.ndarray):
        '''
        description:
            Constructor of Impulse_response class
        inputs:
            h                - A 2D-sequence with the values of the impulse response
            time             - The time axis of h
            delay            - The delay axis of h
        '''
        self.h = h
        self.time = time
        self.delay = delay

class Transfer_function:
    def __init__(self, T: np.ndarray, time: np.ndarray, frequency: np.ndarray):
        '''
        description:
            Constructor of Transfer_function class
        inputs: 
            T                - A 2D-sequence with the values of the transfer function
            time             - The time axis of T
            frequency        - The frequency axis of T
        '''

        self.T = T
        self.time = time
        self.frequency = frequency

class Coherence_bandwidth:
    def __init__(self, B_coh_classic: np.ndarray, B_coh_cyclic: np.ndarray, time: np.ndarray, threshold: numbers.Number):
        '''
        description:
            Constructor of Coherence_bandwidth class
        inputs:
            B_coh_classic    - The sequence of coherence bandwidths over time with classic correlation method
            B_coh_cyclic     - The sequence of coherence bandwidths over time with cyclic correlation method
            time             - The time axis of the coherence bandwidths
            threshold        - Threshold level the coherence bandwiths have been evaluated at
        '''

        self.B_coh_classic = B_coh_classic
        self.B_coh_cyclic = B_coh_cyclic
        self.time = time
        self.treshold_level = threshold

        self.average_classic = None
        self.min_classic = None
        self.max_classic = None

        self.average_cyclic = None
        self.min_cyclic = None
        self.max_cyclic = None
    
    def set_cyclic(self, B_coh_cyclic: np.ndarray, time: np.ndarray, threshold: numbers.Number):
        '''
        description:
            set the coherence bandwith computed with cyclic correlation method
        inputs:
            B_coh_cyclic     - A sequence with the coherence bandwidth values
            time             - The time axis of the coherence bandwidths
            threshold        - The threshold level of the coherence bandwidths
        '''
        self.B_coh_cyclic = B_coh_cyclic
        self.time = time
        self.treshold_level = threshold

    def set_classic(self, B_coh_classic: np.ndarray, time: np.ndarray, threshold: numbers.Number):
        '''
        description:
            set the coherence bandwith computed with classic correlation method
        inputs:
            B_coh_classic    - A sequence with the coherence bandwidth values
            time             - The time axis of the coherence bandwidths
            threshold        - The threshold level of the coherence bandwidths
        '''
        self.B_coh_classic = B_coh_classic
        self.time = time
        self.treshold_level = threshold

    def set_avg_min_max_cyclic(self):
        '''
        description:
            Find average, min and max vaues of the sequence of coherence bandwidths with the cyclic correlation method
        '''
        self.average_cyclic = np.mean(self.B_coh_cyclic)
        self.min_cyclic = np.min(self.B_coh_cyclic)
        self.max_cyclic = np.max(self.B_coh_cyclic)

    def set_avg_min_max_classic(self):
        '''
        description:
            Find average, min and max vaues of the sequence of coherence bandwidths with the classic correlation method
        '''
        self.average_classic = np.mean(self.B_coh_classic)
        self.min_classic = np.min(self.B_coh_classic)
        self.max_classic = np.max(self.B_coh_classic)

class Coherence_time:
    def __init__(self, T_coh_classic: np.ndarray, T_coh_cyclic: np.ndarray,time: np.ndarray, \
            threshold: numbers.Number, windowsize_in_sec: numbers.Number, \
            points_per_window: numbers.Number):
        '''
        description:
            Constructor of Coherence_time class
        inputs:
            T_coh_classic    - A sequence of coherence times with classic correlation method
            T_coh_cyclic     - A sequence of coherence times with cyclic correlation method
            time             - The time axis of coherence times
            threshold        - The correlation threshold leve of the coherence times
            windowsize_in_sec- determines the size of the window, that is used to
                               calculate the moving average function. This
                               parameter is only in use in the cyclic method
            points_per_window- Determines the number of points of the T_coh
                               function. The number of points is the number
                               of windows that fit inside the function side
                               by side, muliplied with points_per_window
        '''
        self.T_coh_classic = T_coh_classic
        self.T_coh_cyclic = T_coh_cyclic
        self.time = time
        self.thresold_level = threshold
        self.windowsize_in_sec = windowsize_in_sec
        self.points_per_window = points_per_window

        self.average_classic = None
        self.min_classic = None
        self.max_classic = None

        self.average_cyclic = None
        self.min_cyclic = None
        self.max_cyclic = None

    def set_cyclic(self, T_coh_cyclic: np.ndarray,time: np.ndarray, \
                threshold: numbers.Number, windowsize_in_sec: numbers.Number, \
                points_per_window: numbers.Number):
        '''
        description:
            set the coherence time for the cyclic correlation method
        inputs:
            T_coh_cyclic     - A sequence of coherence times with cyclic correlation method
            time             - The time axis of coherence times
            threshold        - The correlation threshold leve of the coherence times
            windowsize_in_sec- Determines the size of the window, that is used to
                               calculate the moving average function. This
                               parameter is only in use in the cyclic method
            points_per_window- Determines the number of points of the T_coh
                               function. The number of points is the number
                               of windows that fit inside the function side
                               by side, muliplied with points_per_window
        '''
        self.T_coh_cyclic = T_coh_cyclic
        self.time = time
        self.thresold_level = threshold
        self.windowsize_in_sec = windowsize_in_sec
        self.points_per_window = points_per_window

    def set_classic(self, T_coh_classic: np.ndarray, time: np.ndarray, \
                threshold: numbers.Number, windowsize_in_sec: numbers.Number, \
                points_per_window: numbers.Number):
        '''
        description:
            set the coherence time for the classic correlation method
        inputs:
            T_coh_classic    - A sequence of coherence times with classic correlation method
            time             - The time axis of coherence times
            threshold        - The correlation threshold leve of the coherence times
            windowsize_in_sec- Determines the size of the window, that is used to
                               calculate the moving average function. This
                               parameter is only in use in the cyclic method
            points_per_window- Determines the number of points of the T_coh
                               function. The number of points is the number
                               of windows that fit inside the function side
                               by side, muliplied with points_per_window
        '''
        self.T_coh_classic = T_coh_classic
        self.time = time
        self.thresold_level = threshold
        self.windowsize_in_sec = windowsize_in_sec
        self.points_per_window = points_per_window

    def set_avg_min_max_cyclic(self):
        '''
        description:
            Find average, min and max vaues of the sequence of coherence times with the cyclic correlation method
        '''
        self.average_cyclic = np.mean(self.T_coh_cyclic)
        self.min_cyclic = np.min(self.T_coh_cyclic)
        self.max_cyclic = np.max(self.T_coh_cyclic)

    def set_avg_min_max_classic(self):
        '''
        description:
            Find average, min and max vaues of the sequence of coherence times with the classic correlation method
        '''
        self.average_classic = np.mean(self.T_coh_classic)
        self.min_classic = np.min(self.T_coh_classic)
        self.max_classic = np.max(self.T_coh_classic)

class Measurement:
    def __init__(self, name: str, capture_date: str, capture_time: str, raw_data: np.ndarray, \
                capture_interval_in_seconds: numbers.Number, batchsize: numbers.Number, \
                fs_in_Hz: numbers.Number, fc_in_Hz: numbers.Number, zadoff_chu_sequence: np.ndarray):
        '''
        description:
            Construcor of Measurement class
        inputs:
            name                         - Name of the measurement
            capture_date                 - Date of the capture
            capture_time                 - Time of the capture
            raw_data                     - A sequence containing the raw measurement data
            capture_interval_in_seconds  - The used capture interval in seconds
            batchsize                    - The size of a single batch in samples
            fs_in_Hz                     - The sampling frequency in Hz
            fc_in_Hz                     - The carrier frequency in Hz
            zadoff_chu_sequence          - A sequence containing the sent zadoff-chu-sequence
        '''

        self.name = name
        self.caputre_date = capture_date
        self.capture_time = capture_time
        self.raw_data = raw_data
        self.capture_interval = capture_interval_in_seconds
        self.batchsize = batchsize
        self.fs = fs_in_Hz
        self.fc = fc_in_Hz
        self.zadoff_chu_sequence = zadoff_chu_sequence
        self.length_zadoff_chu_sequence = np.shape(self.zadoff_chu_sequence)[0]

        if not np.shape(self.raw_data)[0] % self.batchsize == 0:
            raise ValueError("The length of the captured measutement is not a multiple of the batchsize. \
                            It is likely that something went wrong with the measutement!")

        self.number_of_batches = int(round(np.shape(self.raw_data)[0]/self.batchsize))
        self.batches = np.reshape(self.raw_data, (self.number_of_batches, self.batchsize))
        self.dft_resolution = self.fs/self.length_zadoff_chu_sequence

        self.rx_power = Rx_power(power_dBFS=None, time=None)
        self.impulse_response = Impulse_response(h=None, time=None, delay=None)
        self.transfer_function = Transfer_function(T=None, time=None, frequency=None)
        self.coherence_bandwidths = Coherence_bandwidth(B_coh_classic=None, time=None, threshold=None, B_coh_cyclic=None)
        self.coherence_times = Coherence_time(T_coh_classic=None, time=None, threshold=None, T_coh_cyclic=None, windowsize_in_sec=None, points_per_window=None)

        self.correlation_coefficient_B_coh_rx_power_classic = None
        self.correlation_coefficient_B_coh_rx_power_cyclic = None

    @classmethod
    def from_filepath(cls, path_to_file: str, path_to_zc_sequence: str):
        '''
        description:
            Alternative constructor that uses a file path to create the class. The name of the file must contain the meta data of the measurement
            in the following format: { Name }_{ Date }_{ Time }_{ fc }MHz_{ fs }MSps_{ capture_interval }ms.dat
        inputs:
            path_to_file         - the path to the file containing the raw measurement data
            path_to_zc_sequence  - the path to the file containing the sent zadoff-chu-sequence
        '''
        try:
            filename = path_to_file.split("/")[-1]
            str = filename.split('_')
            name_ = str[0]
            date_ = str[1]
            time_ = str[2]
            fc_ = int( str[3][0:str[3].find('MHz')] ) * 1e6
            fs_ = int( str[4][0:str[4].find('MSps')] ) * 1e6
            batchsize_ = int( str[5][0:str[5].find('S')] )
            capture_interval_ = int( str[6][0:str[6].find('ms')] ) * 1e-3
            raw_data = np.fromfile(open(f"{path_to_file}"), dtype=np.complex64)
            zadoff_chu_sequence=np.load(f'{path_to_zc_sequence}')
            return cls(name=name_, capture_date=date_, capture_time=time_, raw_data=raw_data,\
                capture_interval_in_seconds=capture_interval_, batchsize=batchsize_, \
                fs_in_Hz=fs_, fc_in_Hz=fc_, zadoff_chu_sequence=zadoff_chu_sequence)
        except(IndexError):
            raise NameError (f"Illegal file name detected: '{filename}'."+ " Format of file should be: { Name }_{ Date }_{ Time }_{ fc }MHz_{ fs }MSps_{ capture_interval }ms.dat")

    def compute_rx_power(self):
        '''
        description:
            computes the recieved power in dBFS from the raw measurement data
        '''
        #compute recieved power in a linear scale
        recieved_power = np.sum(abs(self.batches**2), axis = 1) / self.batchsize
        #convert to log scale
        power_dBFS = 10 * np.log10(recieved_power)
        #create time axis
        time = np.array([i*self.capture_interval for i in range(self.number_of_batches)])
        self.rx_power = Rx_power(power_dBFS, time)
        self.rx_power.set_avg_min_max()

    def compute_impulse_response(self):
        '''
        description:
            compute the impulse response from the raw measurement data. 
            A sampling time offset correction(STO correction) is performed for more accurate results.
        '''
        ##correlate batches with zc sequence

        corr = np.zeros(np.shape(self.batches), dtype = np.complex64)
        for idx,c in enumerate(self.batches):
            # correlate every batch with zc-sequence
            c = cyclic_correlate(c, self.zadoff_chu_sequence)/self.batchsize
            corr[idx] = c
        # rearrange into 3D-array
        corr = np.reshape(corr, (self.number_of_batches,
            int(round(self.batchsize/self.length_zadoff_chu_sequence )), self.length_zadoff_chu_sequence))


        ##sto correction

        upsample_factor = 10
        corrected_corr = np.zeros(np.shape(corr), dtype =np.complex64)
        for idx_batch, batch in enumerate(corr.copy()): # for every batch:
            for idx_seq, seq in enumerate(batch): #for every sequence in a batch:
                # Perform sinc interpolation
                upsampled_sequence = s.resample(seq, self.length_zadoff_chu_sequence * upsample_factor)
                offset = np.argmax(abs(upsampled_sequence)) % upsample_factor
                #create axis for downsampling
                resampled_axis = np.array([(offset + i*upsample_factor) \
                % len(upsampled_sequence-1) for i in range(self.length_zadoff_chu_sequence)])
                # resample the data
                resampled_sequence = upsampled_sequence[resampled_axis]
                # crop the function at the end, if its too long
                if len(resampled_sequence) >self.length_zadoff_chu_sequence:
                    resampled_sequence = resampled_sequence[:self.length_zadoff_chu_sequence]
                # add 0s to the beginning of the function, if it's too short
                if len(resampled_sequence) < self.length_zadoff_chu_sequence:
                    resampled_sequence = np.pad(resampled_sequence, \
                    (0,self.length_zadoff_chu_sequence - len(resampled_sequence)), 'constant')
                corrected_corr[idx_batch][idx_seq] = resampled_sequence

        alignment_sequence = np.zeros(self.length_zadoff_chu_sequence, dtype = np.complex64)
        # Allign all maximums of the impulse responses
        for idx_batch, batch in enumerate(corrected_corr.copy()):
            for idx_response, response in enumerate(batch):
                if idx_response == 0:
                    alignment_sequence = response
                    corrected_corr[idx_batch, idx_response] = response
                    continue
                cross_corr = abs(cyclic_correlate(alignment_sequence, response))
                shift = np.argmax(cross_corr)
                response = np.roll(response, int(self.length_zadoff_chu_sequence/2)+shift)
                corrected_corr[idx_batch, idx_response] = response    


        ## get impulse response

        # average over the multiple impulse responses 
        h = np.mean(corrected_corr, axis = 1)
        #data allocation
        corrected_h = np.zeros(np.shape(h),dtype = np.complex64)
        alignment_sequence = np.zeros(self.length_zadoff_chu_sequence ,dtype = np.complex64)
        #pick a batch and allign the rest to this batch
        for idx_h, h_ in enumerate(h):
            if idx_h == 0:
                alignment_sequence = h_
                corrected_h[idx_h] = h_
                continue
            cross_corr = abs(cyclic_correlate(alignment_sequence, h_))
            shift = np.argmax(cross_corr)
            h_ = np.roll(h_, int(self.length_zadoff_chu_sequence/2)+shift)
            corrected_h[idx_h] = h_
        
        
        h = corrected_h 
        time = np.array([i*self.capture_interval for i in range(self.number_of_batches)])
        delay = np.array([i/self.fs for i in range(np.shape(h)[1])])
        self.impulse_response = Impulse_response(h=h, time=time, delay=delay)
     
    def compute_transfer_function(self):
        '''
        description:
            Compute the transfer function from the impulse response
        '''

        # timevariant transfer function
        T = np.zeros(np.shape(self.impulse_response.h),dtype = np.complex64)
        for idx_h, h_ in enumerate(self.impulse_response.h):
            T[idx_h] = np.fft.fftshift(np.fft.fft(h_, norm = "ortho")) # fft
            # cut away lowpass effect
        self.transfer_function.T = T[:,3:-3]

        self.transfer_function.time = np.array([i*self.capture_interval for i in range(self.number_of_batches)])
        cutoff_frequency = (1 - ( (1-np.shape(self.transfer_function.T)[1]) / np.shape(self.impulse_response.h)[1] ) *self.fs)/2
        self.transfer_function.frequency = np.linspace( -cutoff_frequency, cutoff_frequency, np.shape(self.transfer_function.T)[1])

    def compute_coherence_bandwidth(self, threshold, correlation_method):
        '''
        description:
            compute the coherence bandwidth over time from the transfer function
        inputs:
            threshold            -  The threshold value as a factor to the maximum
                                    value of the autocorrelation function up
                                    to which the function is still considered coherent
            correlation_method   -  calculate B_coh either with autocorrelation
                                    function ("classic") or the cyclic corelation
                                    function ("cyclic")
        '''
        if not 0<=threshold<=1:
            raise ValueError("When calculating B_coh: threshold value \
                            has to be between 0 and 1")
        if not correlation_method in ["classic", "cyclic"]:
            raise ValueError("Correlation_method is required to either be 'classic' or 'cyclic'!")

        #data allocation
        B_coh = np.zeros(np.shape(self.transfer_function.T)[0])
        for idx_t, T_ in enumerate(self.transfer_function.T):
            freq_corr_function = []
            #compute frequency correlation function for a point in time
            if correlation_method == "cyclic":
                freq_corr_function = cyclic_correlate(T_, T_)
            if correlation_method == "classic":
                freq_corr_function = np.correlate(T_,T_,mode = "full")
            #max value of frequency correlation function
            freq_corr_max = abs( np.max(freq_corr_function) )
            #counts all values above the given threshold level
            bandwidth_count = values_above(abs(freq_corr_function), \
            threshold * freq_corr_max)
            #computes the coherence bandwidth
            B_coh[idx_t] = bandwidth_count/2 * self.dft_resolution
        time = np.array([i*self.capture_interval for i in range(self.number_of_batches)])
        if correlation_method == "classic":
            self.coherence_bandwidths.set_classic(B_coh_classic=B_coh, time=time, threshold=threshold)
            self.coherence_bandwidths.set_avg_min_max_classic()
        else:
            self.coherence_bandwidths.set_cyclic(B_coh_cyclic=B_coh, time=time, threshold=threshold)
            self.coherence_bandwidths.set_avg_min_max_cyclic()

    def compute_coherence_time(self, threshold, correlation_method, windowsize_in_sec, points_per_window):
        '''
        description:
            compute the coherence time over time from the transfer function
        inputs:
            threshold            - The threshold value as a factor to the maximum
                                   value of the autocorrelation function up
                                   to which the function is still considered coherent
            correlation_method   -  calculate T_coh either with autocorrelation
                                    function("classic") or the cyclic corelation
                                    function ("cyclic")
            windowsize_in_sec    -  determines the size of the window, that is used to
                                    calculate the moving average function and that
                                    slices the transfer function into parts
            points_per_window    -  Determines the number of points of the T_coh
                                    function. The number of points is the number
                                    of windows that fit inside the function side
                                    by side, muliplied with points_per_window
        '''
        if not 0<=threshold<=1:
            raise ValueError("When calculating B_coh: threshold value \
                            has to be between 0 and 1")
        if not correlation_method in ["classic", "cyclic"]:
            raise ValueError("Correlation_method is required to either be 'classic' or 'cyclic'!")

        batches_in_window = int(round(windowsize_in_sec/self.capture_interval))
        T = np.swapaxes(self.transfer_function.T,0,1)
        increment = batches_in_window/points_per_window
        evaluation_points = np.array([ int(round(i*increment)) \
        for i in range(int((self.number_of_batches-batches_in_window)/increment)+1)])
        T_coh = np.zeros(len(evaluation_points))
        for idx_t, t in enumerate(evaluation_points):
            T_coh_f = np.zeros(np.shape(T)[0])
            for idx_f, f in enumerate(T):
                time_corr_max = 0
                time_corr_function = []

                if correlation_method == "classic":
                    time_corr_function = np.correlate( f[t:t+batches_in_window], \
                    f[t:t+batches_in_window], mode = "full" )
                    time_corr_max = abs(np.max(time_corr_function))

                if correlation_method == "cyclic":
                    time_corr_function = cyclic_correlate(f[t:t+batches_in_window],\
                    f[t:t+batches_in_window])
                    time_corr_max = abs( np.max(time_corr_function) )

                time_count = values_above(abs(time_corr_function),\
                threshold * time_corr_max)

                T_coh_f[idx_f]=time_count/2*self.capture_interval
            T_coh[idx_t] = np.mean(T_coh_f)
        time = evaluation_points * self.capture_interval + 0.5*windowsize_in_sec*np.ones(np.shape(evaluation_points))
        if correlation_method == "classic":
            self.coherence_times.set_classic(T_coh_classic=T_coh, time=time, threshold=threshold, windowsize_in_sec=windowsize_in_sec, points_per_window=points_per_window )
            self.coherence_times.set_avg_min_max_classic()
        else:
            self.coherence_times.set_cyclic(T_coh_cyclic=T_coh, time=time, threshold=threshold, windowsize_in_sec=windowsize_in_sec, points_per_window=points_per_window )
            self.coherence_times.set_avg_min_max_cyclic()
