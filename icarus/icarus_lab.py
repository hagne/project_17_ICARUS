import pandas as pd
import numpy as np
import matplotlib.pylab as _plt
# import os
from atmPy.aerosols.size_distribution import sizedistribution as sd
# import atmPy
from atmPy.aerosols.instruments.POPS import housekeeping
from atmPy.general import timeseries

def read_pops_raw(fname, bin_edges):
    col_names = pd.read_csv(fname, sep=',', nrows=1, header=None,
                            #             index_col=1,
                            #             usecols=np.arange()
                            ).values[0][:-1].astype(str)
    col_names = np.char.strip(col_names)

    data = pd.read_csv(fname, sep=',', skiprows=1, header=None,
                       #             index_col=1,
                       #             usecols=np.arange()
                       )

    data_hk = data.iloc[:, :27]
    data_hk.columns = col_names
    data_hk.index = pd.to_datetime(data_hk['DateTime'], unit='s')
    data_hk.drop('DateTime', axis=1, inplace=True)
    #     hk = atmPy.general.timeseries.TimeSeries(data_hk, sampling_period = 1)
    hk = housekeeping.POPSHouseKeeping(data_hk, sampling_period=1)
    hk.data['Barometric_pressure'] = hk.data['P']
    hk.get_altitude()

    data_hist = data.iloc[:, 27:]
    data_hist.index = data_hk.index

    dist = sd.SizeDist_TS(data_hist, bin_edges, 'numberConcentration')
    dist._data_period = 1
    dist *= 1 / hk.data.Flow_Set.mean()
    return hk, dist  # , col_names

def read_imet(fname, POPS = 14, column_label_line = 17):
    skip = column_label_line - 1
    labels = pd.read_csv(fname, skiprows = skip, nrows = 1).columns.values
    labels = np.append(labels,np.arange(POPS))
    labels = labels.astype(str)
    labels = np.char.strip(labels)

    radio_sonde = pd.read_csv(fname, skiprows=skip + 1, names = labels)
    date_time = radio_sonde['date [y-m-d GMT]'] + radio_sonde['time [h:m:s GMT]']
    radio_sonde.index = pd.to_datetime(date_time)
    radio_sonde = timeseries.TimeSeries(radio_sonde, 1)
    return radio_sonde


def read_cpc(fname, skiprows=17):
    data = pd.read_csv(fname, skiprows=skiprows, usecols=[0, 1], encoding='ISO-8859-1', skipfooter=3, engine='python')

    rein = open(fname, 'r', encoding='ISO-8859-1')
    header = []
    for e in range(skiprows):
        header.append(rein.readline().split(','))
    rein.close()

    for line in header:
        if line[0] == 'Start Date':
            start_date = pd.to_datetime(line[1])
    start_date = '{}-{:02d}-{:02d} '.format(start_date.year, start_date.month, start_date.day)
    data.index = pd.to_datetime(start_date + data.Time)

    data.drop('Time', axis=1, inplace=True)
    data = timeseries.TimeSeries(data, sampling_period=1)
    return data


class tbs_raw(object):
    def __init__(self,
                 fname_pops_dry,
                 fname_pops_wet,
                 fname_imet,
                 fname_cpc,
                 fname_bins,
                 pops_dry_t_offset = None,
                 pops_wet_t_offset = None,
                 remove_artefacts = False,
                 align_and_merge = False,
                 ):
        bin_edges = pd.read_csv(fname_bins).columns.values.astype(float)

        self.hk_wet, self.dist_wet = read_pops_raw(fname_pops_wet, bin_edges)
        self.adjust_time_offset('POPS_wet', pops_wet_t_offset)

        self.hk_dry, self.dist_dry = read_pops_raw(fname_pops_dry, bin_edges)
        self.adjust_time_offset('POPS_dry', pops_dry_t_offset)

        imet = read_imet(fname_imet, column_label_line=19)
        imet.data['altitude (from iMet PTU) [m]'] = imet.data['altitude (from iMet PTU) [km]'] * 1000
        self.imet_raw = imet
        if remove_artefacts:
            _ = imet.remove_artefacts('altitude (from iMet PTU) [m]', inplace=True)

        self.cpc_raw = read_cpc(fname_cpc)

        if align_and_merge:
            self.align_and_merge()


    def plot2check_timeoffset(self, which, offset_t = (0 ,'s'), offset_alt = 0, xylim = (None, None)):
        """
        Arguments
        ----------
        which: str, ['POPS_dry', POPS_wet', 'CPC']
        """
        f ,a = _plt.subplots()
        self.imet.data['altitude (from iMet PTU) [m]'].plot(ax = a, label = 'alt')
        if which == 'POPS_dry':
            hkw = self.hk_dry
        elif which == 'POPS_wet':
            hkw = self.hk_wet


        a.plot(hkw.data.index + np.timedelta64(offset_t[0], offset_t[1]), hkw.data.Altitude + offset_alt, label = which)
        a.set_xlim(xylim[0])
        a.set_ylim(xylim[1])
        a.legend()
        return a

    def adjust_time_offset(self, which, offset_t):
        if not offset_t:
            return

        if which == 'POPS_dry':
            hk = self.hk_dry
            dist = self.dist_dry
        elif which == 'POPS_wet':
            hk = self.hk_wet
            dist = self.dist_wet
        # hk = hk.copy()
        #         dist = dist.copy()
        hk.data.index = hk.data.index + np.timedelta64(offset_t[0], offset_t[1])
        dist.data.index = dist.data.index + np.timedelta64(offset_t[0], offset_t[1])
        return

    def align_and_merge(self):
        data = self.hk_dry.copy()
        data.data.columns = 'POPSdry_' + data.data.columns.values

        pmd_al = self.dist_dry.particle_mean_diameter.align_to(data)
        pmd_al.data.columns = 'POPSdry_' + pmd_al.data.columns.values
        data = data.merge(pmd_al)

        hk_wet_al = self.hk_wet.align_to(data)
        hk_wet_al.data.columns = 'POPSwet_' + hk_wet_al.data.columns.values
        data = data.merge(hk_wet_al)

        pmd_al = self.dist_wet.particle_mean_diameter.align_to(data)
        pmd_al.data.columns = 'POPSwet_' + pmd_al.data.columns.values
        data = data.merge(pmd_al)

        imet_al = self.imet_raw.align_to(data)
        imet_al.data.columns = 'iMet_' + imet_al.data.columns.values
        data = data.merge(imet_al)

        cpc_al = self.cpc_raw.align_to(data)
        cpc_al.data.columns = 'CPC_' + cpc_al.data.columns.values
        data = data.merge(cpc_al)

        data.data['Altitude'] =  data.data['iMet_altitude (from iMet PTU) [m]']
        #         dist.housekeeping = hk
        self.data = data
        return

    def plot_overview_ts(self, alt_lim = (None, None), nc_lim = (None, 200), cpc_lim = (None, None), figh_hight_scale = 1.5):
        f ,a = _plt.subplots(5, sharex=True, gridspec_kw={'hspace': 0})
        f.set_figheight(f.get_figheight() * figh_hight_scale)
        a_alt = a[0]
        a_nc = a[2]
        a_rh = a[3]
        a_t = a[1]
        a_cpc = a[4]


        # cpc
        self.data.data['CPC_Concentration (#/cm³)'].plot(ax = a_cpc)
        a_cpc.set_ylim(cpc_lim)
        a_cpc.set_ylabel('CPC (#/cm^3)')

        # altitude
        self.data.data.Altitude.plot(ax = a_alt, label = 'iMet')
        self.data.data['POPSdry_Altitude'].plot(ax = a_alt, label = 'POPS_dry')
        self.data.data['POPSwet_Altitude'].plot(ax = a_alt, label = 'POPS_wet')
        a_alt.set_ylabel('Alt. (m)')
        a_alt.set_ylim(alt_lim)
        a_alt.legend()

        # particle number
        self.data.data.POPSdry_PartCon.plot(ax = a_nc, label = 'dry')
        self.data.data.POPSwet_PartCon.plot(ax = a_nc, label = 'wet')
        a_nc.set_ylabel('NC (#/cm^3)')
        a_nc.set_ylim(nc_lim)
        a_nc.legend()

        # temperatur
        self.data.data['iMet_iMet air temperature (corrected) [deg C]'].plot(ax = a_t, label = 'iMet')
        self.data.data.POPSdry_Temp.plot(ax = a_t, label = 'POPS_dry')
        self.data.data.POPSwet_Temp.plot(ax = a_t, label = 'POPS_wet')
        a_t.set_ylabel('Temp (°C)')
        a_t.legend()

        # RH
        self.data.data['iMet_iMet humidity [RH %]'].plot(ax = a_rh)
        a_rh.set_ylabel('RH (%)')
        return a

    def plot_overview_vp(self,
                         alt_lim = (None, None),
                         nc_lim = (None, None),
                         temp_lim = (None, None),
                         rh_lim = (None, None),
                         md_lim = (None, None),
                         cpc_lim = (None, None),
                         avg_time = (30 ,'s'),
                         alpha = 0.3,
                         ):

        f ,a = _plt.subplots(1, 5, sharey=True, gridspec_kw={'wspace': 0})
        f.set_figheight(f.get_figheight() * 1.5)
        #     a_alt = a[0]
        a[0].set_ylabel('Altitude (m)')
        a_nc = a[2]
        a_rh = a[1]
        a_t = a[0]
        a_md = a[3]
        a_cpc = a[4]

        #         hk_dry = dist_dry.housekeeping
        #         hk_wet = dist_wet.housekeeping

        #     hk_dry_avg = hk_dry.average_time(avg_time)
        #     hk_wet_avg = hk_wet.average_time(avg_time)

        dist_dry_avg = self.dist_dry.average_time(avg_time)
        dist_wet_avg = self.dist_wet.average_time(avg_time)

        data_avg = self.data.average_time(avg_time)

        # CPC

        a_cpc.plot(self.data.data['CPC_Concentration (#/cm³)'], self.data.data.Altitude)
        a_cpc.plot(data_avg.data['CPC_Concentration (#/cm³)'], data_avg.data.Altitude)
        a_cpc.set_xlim(cpc_lim)
        a_cpc.set_xlabel('CPC (#/cm^3)')


        # particle number

        a_nc.plot(self.data.data.POPSdry_PartCon ,self.data.data.Altitude, alpha = alpha)
        g = a_nc.get_lines()[-1]
        g.set_label(None)
        col1 = g.get_color()
        a_nc.plot(self.data.data.POPSwet_PartCon ,self.data.data.Altitude, alpha = alpha)
        g = a_nc.get_lines()[-1]
        g.set_label(None)
        col2 = g.get_color()
        a_nc.plot(data_avg.data.POPSdry_PartCon ,data_avg.data.Altitude, label = 'dry', color = col1)
        a_nc.plot(data_avg.data.POPSwet_PartCon ,data_avg.data.Altitude, label = 'wet', color = col2)
        a_nc.set_xlabel('NC (#/cm^3)')
        a_nc.set_xlim(nc_lim)
        a_nc.set_ylim(alt_lim)
        a_nc.legend()

        # temperatur
        a_t.plot(self.data.data['iMet_iMet air temperature (corrected) [deg C]'] ,self.data.data.Altitude, label = 'iMet', alpha = alpha)
        g = a_t.get_lines()[-1]
        col1 = g.get_color()
        a_t.plot(data_avg.data['iMet_iMet air temperature (corrected) [deg C]'] ,data_avg.data.Altitude, label = 'iMet', color = col1)
        a_t.set_xlabel('Temp (°C)')
        a_t.legend()
        a_t.set_xlim(temp_lim)
        a_t.set_ylim(alt_lim)

        # RH
        col = 'iMet_iMet humidity [RH %]'
        a_rh.plot(self.data.data[col] ,self.data.data.Altitude, label = 'iMet', alpha = alpha)
        g = a_t.get_lines()[-1]
        col1 = g.get_color()
        a_rh.plot(data_avg.data[col] ,data_avg.data.Altitude, label = 'iMet', color = col1)
        a_rh.set_xlabel('RH (%)')
        a_rh.set_xlim(rh_lim)

        # mean diameter
        a_md.plot(self.data.data['POPSdry_mean d (nm)'] ,self.data.data.Altitude, alpha = alpha)
        g = a_md.get_lines()[-1]
        g.set_label(None)
        col1 = g.get_color()
        a_md.plot(self.data.data['POPSwet_mean d (nm)'] ,self.data.data.Altitude, alpha = alpha)
        g = a_md.get_lines()[-1]
        g.set_label(None)
        col2 = g.get_color()
        a_md.plot(data_avg.data['POPSdry_mean d (nm)'] ,data_avg.data.Altitude, label = 'dry', color = col1)
        a_md.plot(data_avg.data['POPSwet_mean d (nm)'] ,data_avg.data.Altitude, label = 'wet', color = col2)

        a_md.set_xlabel('mean d (nm)')
        a_md.set_xlim(md_lim)
        a_md.legend()
        return a