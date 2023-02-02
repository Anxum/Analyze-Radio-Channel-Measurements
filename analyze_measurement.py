from Measurement import *
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os
import argparse

class CoherenceValues():
    def __init__(self) -> None:
        self.coherence_bandwidths_classic = []
        self.coherence_bandwidths_cyclic = []
        self.coherence_times_classic = []
        self.coherence_times_cyclic = []
        self.names = []

    def addMeasurement(self, measurement: Measurement):
        if not isinstance(measurement.coherence_bandwidths.B_coh_classic,type(None)):
            self.coherence_bandwidths_classic.append(measurement.coherence_bandwidths.B_coh_classic)
        if not isinstance(measurement.coherence_bandwidths.B_coh_cyclic,type(None)):
            self.coherence_bandwidths_cyclic.append(measurement.coherence_bandwidths.B_coh_cyclic)
        if not isinstance(measurement.coherence_times.T_coh_classic,type(None)):
            self.coherence_times_classic.append(measurement.coherence_times.T_coh_classic)
        if not isinstance(measurement.coherence_times.T_coh_cyclic,type(None)):
            self.coherence_times_cyclic.append(measurement.coherence_times.T_coh_cyclic)
        self.names.append(measurement.name)

    def evaluateCorrelationMethods(self, outpath):
        if len(self.coherence_bandwidths_classic) == 0:
            raise ValueError("Classic coherence bandwidth must not be an empty array!")
        if len(self.coherence_bandwidths_cyclic) == 0:
            raise ValueError("Cyclic coherence bandwidth must not be an empty array!")
        if len(self.coherence_times_classic) == 0:
            raise ValueError("Classic coherence time must not be an empty array!")
        if len(self.coherence_times_cyclic) == 0:
            raise ValueError("Cyclic coherence time must not be an empty array!")

        mins_coherence_bandwidth_classic = []
        mins_coherence_bandwidth_cyclic = []
        mins_coherence_time_classic = []
        mins_coherence_time_cyclic = []

        maxs_coherence_bandwidth_classic = []
        maxs_coherence_bandwidth_cyclic = []
        maxs_coherence_time_classic = []
        maxs_coherence_time_cyclic = []

        avgs_coherence_bandwidth_classic = []
        avgs_coherence_bandwidth_cyclic = []
        avgs_coherence_time_classic = []
        avgs_coherence_time_cyclic = []

        for b_cl,b_cy,t_cl,t_cy in zip(self.coherence_bandwidths_classic, self.coherence_bandwidths_cyclic, self.coherence_times_classic, self.coherence_times_cyclic):
            mins_coherence_bandwidth_classic.append(np.min(b_cl,axis=0))
            mins_coherence_bandwidth_cyclic.append(np.min(b_cy,axis=0))
            mins_coherence_time_classic.append(np.min(t_cl,axis=0))
            mins_coherence_time_cyclic.append(np.min(t_cy,axis=0))

            maxs_coherence_bandwidth_classic.append(np.max(b_cl,axis=0))
            maxs_coherence_bandwidth_cyclic.append(np.max(b_cy,axis=0))
            maxs_coherence_time_classic.append(np.max(t_cl,axis=0))
            maxs_coherence_time_cyclic.append(np.max(t_cy,axis=0))

            avgs_coherence_bandwidth_classic.append(np.mean(b_cl,axis=0))
            avgs_coherence_bandwidth_cyclic.append(np.mean(b_cy,axis=0))
            avgs_coherence_time_classic.append(np.mean(t_cl,axis=0))
            avgs_coherence_time_cyclic.append(np.mean(t_cy,axis=0))

        mean_b_coh_classic_avg = np.mean(np.array(avgs_coherence_bandwidth_classic))
        mean_b_coh_cyclic_avg = np.mean(np.array(avgs_coherence_bandwidth_cyclic))
        mean_t_coh_classic_avg = np.mean(np.array(avgs_coherence_time_classic))
        mean_t_coh_cyclic_avg = np.mean(np.array(avgs_coherence_time_cyclic))

        mean_b_coh_classic_min = np.mean(np.array(mins_coherence_bandwidth_classic))
        mean_b_coh_cyclic_min = np.mean(np.array(mins_coherence_bandwidth_cyclic))
        mean_t_coh_classic_min = np.mean(np.array(mins_coherence_time_classic))
        mean_t_coh_cyclic_min = np.mean(np.array(mins_coherence_time_cyclic))

        mean_b_coh_classic_max = np.mean(np.array(maxs_coherence_bandwidth_classic))
        mean_b_coh_cyclic_max = np.mean(np.array(maxs_coherence_bandwidth_cyclic))
        mean_t_coh_classic_max = np.mean(np.array(maxs_coherence_time_classic))
        mean_t_coh_cyclic_max = np.mean(np.array(maxs_coherence_time_cyclic))

        s = "Average classic coherence bandwidth: {:.2f}MHz \n".format(mean_b_coh_classic_avg*1e-6)
        s = s+"Average cyclic coherence bandwidth: {:.2f}MHz \n".format(mean_b_coh_cyclic_avg*1e-6)
        s = s+"Average classic coherence time: {:.2f}ms \n".format(mean_t_coh_classic_avg*1e3)
        s = s+"Average cyclic coherence time: {:.2f}ms \n".format(mean_t_coh_cyclic_avg*1e3)
        s = s+"\n"

        s = s+"Average minimum classic coherence bandwidth: {:.2f}MHz \n".format(mean_b_coh_classic_min*1e-6)
        s = s+"Average minimum cyclic coherence bandwidth: {:.2f}MHz \n".format(mean_b_coh_cyclic_min*1e-6)
        s = s+"Average minimum classic coherence time: {:.2f}ms \n".format(mean_t_coh_classic_min*1e3)
        s = s+"Average minimum cyclic coherence time: {:.2f}ms \n".format(mean_t_coh_cyclic_min*1e3)
        s = s+"\n"

        s = s+"Average maximum classic coherence bandwidth: {:.2f}MHz \n".format(mean_b_coh_classic_max*1e-6)
        s = s+"Average maximum cyclic coherence bandwidth: {:.2f}MHz \n".format(mean_b_coh_cyclic_max*1e-6)
        s = s+"Average maximum classic coherence time: {:.2f}ms \n".format(mean_t_coh_classic_max*1e3)
        s = s+"Average maximum cyclic coherence time: {:.2f}ms \n".format(mean_t_coh_cyclic_max*1e3)

        textfile = open(f"{outpath}/min_max_avg_values.txt", "w")
        textfile.write(s)
        textfile.close()





    def computeAbsoluteValues(self):
        raise NotImplementedError
        #TODO check for computed correlation methods
        #TODO find absolute min, max, average over all measurements


        


def b_t_p(measurement: Measurement, correlation_method: str, outpath: str):
    x_axis_allign = (measurement.rx_power.time[-1]-measurement.rx_power.time[-1]*1.05,measurement.rx_power.time[-1]*1.05)

    if correlation_method == "classic":
        B_coh = measurement.coherence_bandwidths.B_coh_classic
        T_coh = measurement.coherence_times.T_coh_classic
    if correlation_method == "cyclic":
        B_coh = measurement.coherence_bandwidths.B_coh_cyclic
        T_coh = measurement.coherence_times.T_coh_cyclic

    fig, axs = plt.subplots(3,1,constrained_layout = True)
    axs[0].plot(measurement.rx_power.time, measurement.rx_power.power_dBFS)
    axs[0].set_xlabel("Time [s]")
    axs[0].set_ylabel("Power [dBFS]")
    axs[0].set_ylim((-60,0))
    axs[0].set_xlim(x_axis_allign)
    axs[0].hlines(measurement.rx_power.average, 0, measurement.rx_power.time[-1], linestyles= "dotted", label = "Mean Power")
    axs[0].legend()

    fig.suptitle(f"Power, cohenece bandwith and time of '{measurement.name}'")
    axs[1].plot(measurement.coherence_bandwidths.time, B_coh*1e-6)
    axs[1].set_xlabel("Time [s]")
    axs[1].set_ylabel(r"$B_{coh}$ [MHz]")
    axs[1].set_ylim((-0.3,15))
    axs[1].set_xlim(x_axis_allign)

    axs[2].plot(measurement.coherence_times.time, T_coh*1e3)
    axs[2].set_xlabel("Time [s]")
    axs[2].set_ylabel(r"$T_{coh}$ [ms]")
    axs[2].set_ylim((-3,500))
    axs[2].set_xlim(x_axis_allign)

    plt.savefig(f"{outpath}/{measurement.name}/rx_power_and_coherence_over_time_{correlation_method}.pdf")
    plt.close()

def covariance_plot(measurement: Measurement, correlation_method: str, outpath: str):

    if correlation_method == "classic":
        B_coh = measurement.coherence_bandwidths.B_coh_classic*1e-6
    if correlation_method == "cyclic":
        B_coh = measurement.coherence_bandwidths.B_coh_cyclic*1e-6
    power = measurement.rx_power.power_dBFS
    f = plt.figure(num="Power vs Coherence Bandwidth", figsize = (6,2.5))
    f.subplots_adjust(bottom=0.17)
    x = np.array(np.append([power], [B_coh], axis = 0))
    if correlation_method == "classic":
        measurement.correlation_coefficient_B_coh_rx_power_classic = round(np.corrcoef(x)[1,0],2)
    else:
        measurement.correlation_coefficient_B_coh_rx_power_cyclic = round(np.corrcoef(x)[1,0],2)
    plt.scatter(power, B_coh, marker = ".")
    plt.xlabel(f"RX power [dBFS]")
    plt.ylabel(f"Coherence Bandwidth [MHz]")
    plt.xlim([-60,0])
    if correlation_method == "classic":
        plt.ylim([-0.2,5])
    if correlation_method == "cyclic":
        plt.ylim([0,15])
    plt.title(f"Coherence Bandwidth over RX power of '{measurement.name}' - {correlation_method}")
    x = np.array([np.min(power)*1.1,np.max(power)*0.9])
    b, m = polyfit(power, B_coh, 1)
    plt.plot(x, b + m * x, '-', color = "r")
    plt.savefig(f"{outpath}/{measurement.name}/scatterplot_{correlation_method}.pdf")
    plt.close()

def location_of_percentile(percentile, hist):
    '''
    description:
        Find the location of the percentile in the plot

    inputs:
        percentile       -  The percentile in use
        hist             -  The histogram

    returns:
        location         -  The exact loction of the percentile inside the plot
    '''
    all = hist[0].sum()
    running_count = 0
    percentile_bin = 0
    for frequency in hist[0]:
        if running_count >= all*percentile*0.01:
            break
        running_count = running_count+frequency
        percentile_bin = percentile_bin+1
    return hist[1][percentile_bin]

def plot_percentiles(percentiles, hist, mean, unit, fading_margin = False):
    y_max = np.max(hist[0])
    colors = ["magenta","red","orange","purple"]
    for p,c in zip(percentiles, colors):
        if not 0<=p<=100:
            raise ValueError("percentiles have to be between 0 and 100!")
        x_location = location_of_percentile(p,hist)

        plt.vlines(x_location, 0, y_max, linestyles = "dotted", label = r"$F^{-1}$"+"({:.0%}) = {:.2f}".format(p/100,x_location)+f"{unit}", color = c)
    if fading_margin:
        plt.vlines(mean, 0, y_max, linestyles = "dotted", label = f"mean = {round(mean,1)}{unit}", color = "black")
        plt.ylim([0, y_max*1.2])
    plt.legend()

def histogramm_plot(measurement: Measurement, correlation_method: str, outpath: str):
    if correlation_method == "classic":
        B_coh = measurement.coherence_bandwidths.B_coh_classic*1e-6
        T_coh = measurement.coherence_times.T_coh_classic*1e3
    if correlation_method == "cyclic":
        B_coh = measurement.coherence_bandwidths.B_coh_cyclic*1e-6
        T_coh = measurement.coherence_times.T_coh_cyclic*1e3

    plt.figure()
    plt.title(f'Histogramm of T_coh - {measurement.name}')
    hist_T = plt.hist(T_coh, 125, (-3,500))
    plot_percentiles([1,10,90,99], hist_T,np.mean(T_coh), "ms")
    plt.xlabel(r"$T_{coh}$[ms]")
    plt.xlim([0,500])
    plt.ylabel("Count of T_coh values")
    plt.savefig(f"{outpath}/{measurement.name}/Histogram_of_T_coh{correlation_method}.pdf")
    plt.close()

    plt.figure()
    plt.title(f'Histogramm of B_coh - {measurement.name}')
    hist_B = plt.hist(np.array(B_coh), 125,(-0.3, 15) )
    plot_percentiles([1,10,90,99], hist_B, np.mean(B_coh), "MHz")
    plt.xlabel(r"$B_{coh}$[MHz]")
    plt.xlim([0, 15.36])
    plt.ylabel("Count of B_coh values")
    plt.savefig(f"{outpath}/{measurement.name}/Histogram_of_B_coh_{correlation_method}.pdf")
    plt.close()

    f = plt.figure(figsize= (7,3))
    f.subplots_adjust(bottom= 0.15)
    plt.title(f"Histogramm of RX power '{measurement.name}'")
    hist_P = plt.hist(np.array(measurement.rx_power.power_dBFS), 125, (-60,0))
    plot_percentiles([1,10], hist_P, measurement.rx_power.average, "dBFS", fading_margin = True)
    plt.xlabel("Received power [dBFS]")
    plt.xlim([-60, 0])
    plt.ylabel("Amount of samples")
    plt.savefig(f"{outpath}/{measurement.name}/Histogram_of_RX_power.pdf")
    plt.close()

def data_out(m: Measurement, opt):
    name= f"Name: {m.name}"
    date= f"Date: {m.caputre_date}"
    time= f"Time: {m.capture_time}"
    fs =  f"fs:   {round(m.fs*1e-6,2)}MHz"
    fc=   f"fc:   {round(m.fc*1e-6,2)}MHz"
    meta_data = f"Meta Data:\n\t{name}\n\t{date}\n\t{time}\n\t{fs}\n\t{fc}\n\n"

    power_avg="Mean RX power: {:.2f}dBFS".format(m.rx_power.average,2)
    power_min="Min RX power:  {:.2f}dBFS".format(m.rx_power.min,2)
    power_max="Max RX power:  {:.2f}dBFS".format(m.rx_power.max,2)
    power = f"Power:\n\t{power_avg}\n\t{power_min}\n\t{power_max}\n\n"

    if opt.corr_method in ["both", "classic"]:
        B_coh_classic_avg="Mean classic B_coh({:.0%}): {:.2f}MHz".format(m.coherence_bandwidths.treshold_level, m.coherence_bandwidths.average_classic*1e-6)
        B_coh_classic_min="Min classic B_coh({:.0%}):  {:.2f}MHz".format(m.coherence_bandwidths.treshold_level, m.coherence_bandwidths.min_classic*1e-6)
        B_coh_classic_Max="Max classic B_coh({:.0%}):  {:.2f}MHz".format(m.coherence_bandwidths.treshold_level, m.coherence_bandwidths.max_classic*1e-6)
        B_coh_classic =f"Coherence bandwidth classic:\n\t{B_coh_classic_avg}\n\t{B_coh_classic_min}\n\t{B_coh_classic_Max}\n\n"

    if opt.corr_method in ["both", "cyclic"]:
        B_coh_cyclic_avg=f"Mean cyclic B_coh({int(m.coherence_bandwidths.treshold_level*100)}%): {round(m.coherence_bandwidths.average_cyclic*1e-6,2)}MHz"
        B_coh_cyclic_min=f"Min cyclic B_coh({int(m.coherence_bandwidths.treshold_level*100)}%):  {round(m.coherence_bandwidths.min_cyclic*1e-6,2)}MHz"
        B_coh_cyclic_Max=f"Max cyclic B_coh({int(m.coherence_bandwidths.treshold_level*100)}%):  {round(m.coherence_bandwidths.max_cyclic*1e-6,2)}MHz"
        B_coh_cyclic = f"Coherence bandwidth cyclic:\n\t{B_coh_cyclic_avg}\n\t{B_coh_cyclic_min}\n\t{B_coh_cyclic_Max}\n\n"

    if opt.corr_method in ["both", "classic"]:
        T_coh_classic_avg=f"Mean classic T_coh({int(m.coherence_times.thresold_level*100)}%): {round(m.coherence_times.average_classic*1e3,2)}ms"
        T_coh_classic_min=f"Min classic T_coh({int(m.coherence_times.thresold_level*100)}%):  {round(m.coherence_times.min_classic*1e3,2)}ms"
        T_coh_classic_Max=f"Max classic T_coh({int(m.coherence_times.thresold_level*100)}%):  {round(m.coherence_times.max_classic*1e3,2)}ms"
        T_coh_classic = f"Coherence time classic:\n\t{T_coh_classic_avg}\n\t{T_coh_classic_min}\n\t{T_coh_classic_Max}\n\n"

    if opt.corr_method in ["both", "cyclic"]:
        T_coh_cyclic_avg=f"Mean cyclic T_coh({int(m.coherence_times.thresold_level*100)}%): {round(m.coherence_times.average_cyclic*1e3,2)}ms"
        T_coh_cyclic_min=f"Min cyclic T_coh({int(m.coherence_times.thresold_level*100)}%):  {round(m.coherence_times.min_cyclic*1e3,2)}ms"
        T_coh_cyclic_Max=f"Max cyclic T_coh({int(m.coherence_times.thresold_level*100)}%):  {round(m.coherence_times.max_cyclic*1e3,2)}ms"
        T_coh_cyclic = f"Coherence time cyclic:\n\t{T_coh_cyclic_avg}\n\t{T_coh_cyclic_min}\n\t{T_coh_cyclic_Max}\n\n"

    if opt.corr_method in ["both", "classic"]:
        Corr_coeff_classic = f"Correlation RX-power - B_coh classic: {round(m.correlation_coefficient_B_coh_rx_power_classic,2)}"
    if opt.corr_method in ["both", "cyclic"]:
        Corr_coeff_cyclic =  f"Correlation RX-power - B_coh cyclic:  {round(m.correlation_coefficient_B_coh_rx_power_cyclic,2)}"

    if opt.corr_method in ["both"]:
        corr = f"Correlation coefficients:\n\t{Corr_coeff_classic}\n\t{Corr_coeff_cyclic}"
    if opt.corr_method in ["classic"]:
        corr = f"Correlation coefficients:\n\t{Corr_coeff_classic}"
    if opt.corr_method in ["cyclic"]:
        corr = f"Correlation coefficients:\n\t{Corr_coeff_cyclic}"
    
    if opt.corr_method == "both":
        str_out = meta_data+power+B_coh_classic+B_coh_cyclic+T_coh_classic+T_coh_cyclic+corr
    if opt.corr_method == "classic":
        str_out = meta_data+power+B_coh_classic+T_coh_classic+corr
    if opt.corr_method == "cyclic":
        str_out = meta_data+power+B_coh_cyclic+T_coh_cyclic+corr

    textfile = open(f"{opt.outpath}/{m.name}/data_of_{m.name}.txt", "w")
    textfile.write(str_out)
    textfile.close()


def compute_measurement(measurement: Measurement, options): 
    measurement.compute_impulse_response()
    measurement.compute_transfer_function()
    measurement.compute_rx_power()
    if not os.path.exists(f"{options.outpath}/{measurement.name}"):
        os.mkdir(f"{options.outpath}/{measurement.name}")
    if options.corr_method in ["both", "classic"]:
        measurement.compute_coherence_bandwidth(threshold=options.thresh, correlation_method="classic")
        measurement.compute_coherence_time(threshold=options.thresh, correlation_method="classic",\
                                windowsize_in_sec=options.window_size,points_per_window=options.pts_per_window)
        b_t_p(measurement, "classic", options.outpath)
        covariance_plot(measurement, "classic", options.outpath)
        histogramm_plot(measurement,"classic", options.outpath)

    if options.corr_method in ["both", "cyclic"]:
        measurement.compute_coherence_bandwidth(threshold=options.thresh, correlation_method="cyclic")
        measurement.compute_coherence_time(threshold=options.thresh, correlation_method="cyclic", \
                                windowsize_in_sec=options.window_size,points_per_window=options.pts_per_window)
        b_t_p(measurement, "cyclic", options.outpath)
        covariance_plot(measurement, "cyclic", options.outpath)
        histogramm_plot(measurement,"cyclic", options.outpath)
    data_out(measurement,options)

def locate_measurements(options):
    if not os.path.exists(options.outpath):
        os.mkdir(options.outpath)

    if os.path.isdir(options.inpath):
        coherence_values = CoherenceValues()
        for n, file in enumerate(os.listdir(options.inpath)):
            if file.endswith(".dat"):
                try:
                    measurement=Measurement.from_filepath(f"{options.inpath}/{file}", options.zc_path)
                    compute_measurement(measurement, options)
                    coherence_values.addMeasurement(measurement)
                except NameError as error:
                    print(error) 
            print("{}/{} measurements are computed ({:.2f}%)".format(n+1,len(os.listdir(options.inpath)), (n+1)/len(os.listdir(options.inpath))*100))
        coherence_values.evaluateCorrelationMethods(options.outpath)
    if os.path.isfile(options.inpath):
        compute_measurement(Measurement.from_filepath(f"{options.inpath}", options.zc_path), options)        
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog= 'analyze_measurements',
                        description= 'Analyze one or multiple radio channel measutements. Either specify a measurement or a folder to multiple measurements.'
                        )

    parser.add_argument("--inpath", type=str, default="input", help="Specifies the path to a single measurement or to a folder containing measurements.")
    parser.add_argument("--zc-path", type=str, default="zc_sequence.npy", help="Specifies the filepath to the zadoff-chu-sequence.")
    parser.add_argument("--outpath", type=str,default="output", help="Specifies the path to store the results.")
    parser.add_argument("--thresh",type=float, default=0.9, help= "Specifies the threshold for coherence calculation. Values are accepted from 0 to 1.")
    parser.add_argument("--corr-method", type=str, default="both", choices=["classic", "cyclic", "both"], help="Specifies the correlation method in use. Options are: 'classic', 'cyclic' or 'both'.")
    parser.add_argument("--window-size", type=float, default=1,  help= "Specifies the window size in seconds to calculate the coherence time.")
    parser.add_argument("--pts-per-window",type=int, default=10, help= "Specifies the amount of points used per window to calculate the coherence time.")

    options = parser.parse_args()
    if options.inpath[-1] == "/":
        options.inpath = options.inpath[:-1]

    if options.outpath[-1] == "/":
        options.outpath = options.outpath[:-1]

    print("Analyzing measurements with the following arguments:")
    print(options)

    locate_measurements(options=options)