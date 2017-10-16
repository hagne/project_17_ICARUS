import pandas as _pd
import numpy as _np
import matplotlib.pylab as _plt
# import os
import atmPy.aerosols.size_distribution.sizedistribution as _sd
# import atmPy
import atmPy.aerosols.instruments.POPS.housekeeping as _housekeeping
import atmPy.general.timeseries as _timeseries
from matplotlib.ticker import MaxNLocator as _MaxNLocator
from matplotlib.ticker import LogLocator as _LogLocator
from matplotlib.ticker import LogFormatterSciNotation as _LogFormatterSciNotation
from matplotlib.colors import LogNorm as _LogNorm
from matplotlib.dates import DateFormatter as _DateFormatter
import atmPy.data_archives.arm as arm
import plt_tools as _plt_tools
import warnings as _warnings


def read_pops_raw(fname, bin_edges):
    if type(fname).__name__ == 'str':
        col_names = _pd.read_csv(fname, sep=',', nrows=1, header=None,
                                 #             index_col=1,
                                 #             usecols=np.arange()
                                 ).values[0][:-1].astype(str)
        col_names = _np.char.strip(col_names)

        data = _pd.read_csv(fname, sep=',', skiprows=1, header=None,
                            #             index_col=1,
                            #             usecols=np.arange()
                            )

        data_hk = data.iloc[:, :27]
        data_hk.columns = col_names
        data_hk.index = _pd.to_datetime(data_hk['DateTime'], unit='s')
        data_hk.drop('DateTime', axis=1, inplace=True)
        #     hk = atmPy.general.timeseries.TimeSeries(data_hk, sampling_period = 1)
        hk = _housekeeping.POPSHouseKeeping(data_hk, sampling_period=1)

        data_hist = data.iloc[:, 27:]
        data_hist.index = data_hk.index

        dist = _sd.SizeDist_TS(data_hist, bin_edges, 'numberConcentration')
        dist._data_period = 1
        dist *= 1 / hk.data.Flow_Set.mean()

    elif type(fname).__name__ == 'tuple':
        dist = _sd.read_csv(fname[0])
        dist._data_period = 1
        hk = _housekeeping.read_file(fname[1])

    hk.data['Barometric_pressure'] = hk.data['P']
    hk.get_altitude()
    return hk, dist  # , col_names

def read_pops_sd_qc(fname, bins='standard'):
    """This is data that has been preprocessed by Fan"""
    df = _pd.read_csv(fname)
    df.index = _pd.to_datetime(df.UTC)
    df.drop('UTC',axis=1, inplace=True)

    sd = df.iloc[:,1:]
    if bins == 'standard':
        bineg = _np.array([135., 159.11916648, 187.54747513, 221.0547994,
         260.54855872, 307.09829254, 361.96462472, 426.63340283,
         502.85593669, 592.69642598, 698.58786133, 823.39791267,
         970.50658925, 1143.89777444, 1348.26711416, 1589.14918076,
         1873.06735601, 2207.71049227, 2602.14113606, 3067.04095292,
         3615.])
    else:
        raise ValueError('not implemented yet')
        # bincenters = _np.array([float(i.replace('Dp_','')) for i in sd.columns.values])
        # bineg, colnames = atmPy.aerosols.size_distribution.diameter_binning.bincenters2binsANDnames(bincenters)
        # sd.columns = colnames
    distts = _sd.SizeDist_TS(sd,bineg, distType = 'numberConcentration', ignore_data_gap_error = False)
    distts._data_period = 1
    return distts

def read_tbs_qc(fname):
    """This is data that has been preprocessed by Fan"""
    # fname = '/Users/htelg/data/2017_ICARUS/Fan_qc_controled_email_170905/oli.tbs.2017_iop1_iop2/oli.tbs.20170523.172028.txt'
    data = _pd.read_csv(fname)
    data.index = _pd.to_datetime(data.UTC)
    data = data.drop('UTC', axis=1)
    data = _timeseries.TimeSeries(data,  sampling_period = 1)
    return data

def read_imet(fname, POPS = 14, column_label_line = 17):
    skip = column_label_line - 1
    labels = _pd.read_csv(fname, skiprows = skip, nrows = 1).columns.values
    labels = _np.append(labels, _np.arange(POPS))
    labels = labels.astype(str)
    labels = _np.char.strip(labels)

    radio_sonde = _pd.read_csv(fname, skiprows=skip + 1, names = labels)
    date_time = radio_sonde['date [y-m-d GMT]'] + radio_sonde['time [h:m:s GMT]']
    radio_sonde.index = _pd.to_datetime(date_time)
    radio_sonde = _timeseries.TimeSeries(radio_sonde, 1)
    return radio_sonde


def read_cpc(fname, skiprows=17):
    data = _pd.read_csv(fname, skiprows=skiprows, usecols=[0, 1], encoding='ISO-8859-1', skipfooter=3, engine='python')

    rein = open(fname, 'r', encoding='ISO-8859-1')
    header = []
    for e in range(skiprows):
        header.append(rein.readline().split(','))
    rein.close()

    for line in header:
        if line[0] == 'Start Date':
            start_date = _pd.to_datetime(line[1])
    start_date = '{}-{:02d}-{:02d} '.format(start_date.year, start_date.month, start_date.day)
    data.index = _pd.to_datetime(start_date + data.Time)

    data.drop('Time', axis=1, inplace=True)
    data = _timeseries.TimeSeries(data, sampling_period=1)
    return data

def save_figure(tbs,f, name_base):
    st = tbs.data_ts.get_timespan()[0]
    ts = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
    f.savefig('/Users/htelg/projecte/17_ICARUS_aerosol_cloud_interaction/all_flights_figures/{}_{}.png'.format(name_base,ts))

def plot_on_clouds(tbs, cpc = True, relative_humidity = True, timestamp = True, altitude_column = 'Altitude_iMet', which_pops = 'wet', cloud_vmin = 1e2, linewidth = 3.5):
    bsts = tbs.ceilometer.backscatter.get_timespan()
    tsts = tbs.data_ts.get_timespan()
    if not ((bsts[1] > tsts[0]) and (tsts[1] > bsts[0])):
        raise ValueError('there is no overlap of the tbs and the ceilometer data!')

    f,aa = _plt.subplots(4, sharex=True, gridspec_kw={'hspace' : 0})
    a_rh, a_t, a_nc, a_cpc= aa
    bs = tbs.ceilometer.backscatter

    ts = tbs.data_ts.average_time((10, 's'))
    time = ts.data.index.values
    altitude = ts.data[altitude_column].values
    ###################
    # relative humidity
    if not relative_humidity:
        txt = 'no RH'
        a_rh.text(0.5, 0.5, txt, transform=a_rh.transAxes, fontsize='xx-large', va='center', ha='center')
        lc,cb = None,None
        a = a_rh
    else:
        f,a,pc,cb = bs.plot(ax = a_rh, cb_kwargs=False)
        # f.set_figwidth(f.get_figwidth() * 2)
        pc.set_cmap(_plt_tools.colormap.my_colormaps.clouds())
        pc.set_norm(_LogNorm())
        pc.set_clim(vmin = cloud_vmin)#, vmax= 1e4)
        # f.colorbar(pc)
    #     a.set_ylim(top=alt_max)



        RH = ts.data['iMet_iMet humidity [RH %]'].values
        a,lc,cb = _plt_tools.plot.plot_gradiant_color(time, altitude, RH, ax = a, colorbar=False)
        lc.set_cmap(_plt_tools.colormap.my_colormaps.relative_humidity())
        lc.set_linewidth(linewidth)
        lc.set_clim(0, 100)

        cb, cax = _plt_tools.colorbar.colorbar_axis_split_off(lc, a)
        cb.locator = _MaxNLocator(4, prune = 'both')
        cb.update_ticks()
        cax.set_ylabel('RH (%)', labelpad=0.5)
        # cax.yaxis.set_major_locator()
    set_rh = (a, lc, cb)

    ###################
    # Temperature
    f,a,pc,cb = bs.plot(ax = a_t, cb_kwargs=False)
    pc.set_cmap(_plt_tools.colormap.my_colormaps.clouds())
    pc.set_norm(_LogNorm())
    pc.set_clim(vmin = cloud_vmin)#, vmax= 1e4)
#     a.set_ylim(top=alt_max)

    temp = ts.data['iMet_iMet air temperature (corrected) [deg C]'].values
    a,lc,cb = _plt_tools.plot.plot_gradiant_color(time, altitude, temp, ax = a, colorbar=False)
    lc.set_cmap(_plt_tools.colormap.my_colormaps.relative_humidity(reverse=True))
    lc.set_linewidth(linewidth)
    lc.set_clim(_np.nanmin(temp),_np.nanmax(temp))

    cb, cax = _plt_tools.colorbar.colorbar_axis_split_off(lc, a)
    cb.locator = _MaxNLocator(5, prune = 'both')
    cb.update_ticks()
    cax.set_ylabel('Temp (°C)', labelpad=0.5)
    set_t = (a,lc,cb)

    ##################
    # POPS number concentration
    f,a,pc,cb = bs.plot(ax = a_nc)
    pc.set_cmap(_plt_tools.colormap.my_colormaps.clouds())
    pc.set_norm(_LogNorm())
    pc.set_clim(vmin = cloud_vmin)#, vmax= 1e4)
#     a.set_ylim(top=alt_max)

    if which_pops == 'wet':
        column = 'POPSwet_PartCon_fromsizedist'
    if which_pops == 'dry':
        column = 'POPSdry_PartCon_fromsizedist'
    nc = ts.data[column]
    a,lc,cb = _plt_tools.plot.plot_gradiant_color(time, altitude, nc, ax = a, colorbar=False)
    lc.set_cmap(_plt_tools.colormap.my_colormaps.relative_humidity(reverse=True))
    lc.set_linewidth(linewidth)
    lc.set_linestyle('-')

    cb, cax = _plt_tools.colorbar.colorbar_axis_split_off(lc, a)
    cb.locator = _MaxNLocator(5, prune = 'both')
    cb.update_ticks()
    cax.set_ylabel('NC$_{POPS}$ (#/cm$^3$)', labelpad=0.5)
    set_nc = (a, lc, cb)

    ###################
    # CPC

    f, a, pc, cb = bs.plot(ax=a_cpc, cb_kwargs=False)
    pc.set_cmap(_plt_tools.colormap.my_colormaps.clouds())
    pc.set_norm(_LogNorm())
    pc.set_clim(vmin=cloud_vmin)  # , vmax= 1e4)
    if cpc:
        cpc = ts.data['CPC_Concentration (#/cm³)'].values

        a, lc, cb = _plt_tools.plot.plot_gradiant_color(time[~ _np.isnan(cpc)], altitude[~ _np.isnan(cpc)], cpc[~ _np.isnan(cpc)], ax=a, colorbar=False)
        lc.set_cmap(_plt_tools.colormap.my_colormaps.relative_humidity(reverse=True))
        lc.set_linewidth(linewidth)
        lc.set_clim(_np.nanmin(cpc), _np.nanmax(cpc))
        lc.set_norm(_LogNorm())

        # try:
        cb, cax = _plt_tools.colorbar.colorbar_axis_split_off(lc, a,
                                                          cb_kwargs = {'format': '%.0e'}
                                                          )
        cb.locator = _LogLocator()
        cb.update_ticks()
        # except:
            # pass
            # from matplotlib.colors import NoNorm
            # lc.set_norm(NoNorm())
            # cb, cax = _plt_tools.colorbar.colorbar_axis_split_off(lc, a,
            #                                                       # cb_kwargs={'format': '%.0e'}
            #                                                       )
        # cb.formatter = _LogFormatterSciNotation()

        # cb.locator = _MaxNLocator(5, prune='both')
        # cb.update_ticks()

        cax.set_ylabel('NC$_{CPC}$ (#/cm$^3$)', labelpad=0.5)
        set_cpc = (a, lc, cb)
    else:
        txt = 'no CPC'
        a_cpc.text(0.5, 0.5, txt, transform=a_cpc.transAxes, fontsize='xx-large', va='center', ha='center')
        # lc,cb = None,None
        # a = a_cpc
        set_cpc = (a_cpc,None,None)
    ##############

    _plt_tools.axes.labels.set_shared_label(aa, 'Altitude (m)', axis='y', labelpad=0.01)

    at = aa[-1]
    at.xaxis.set_major_formatter(_plt.DateFormatter("%H:%M:%S"))
    at.set_xlabel('')

    for at in aa:
        at.yaxis.set_major_locator(_MaxNLocator(5, prune='both'))

    f.set_figheight(f.get_figheight() * 1.2)
    f.tight_layout()
    f.patch.set_alpha(0)

    f.save = lambda: save_figure(tbs, f, 'plot_on_clouds')

    if timestamp:
        at = aa[0]
        st = tbs.data_ts.get_timespan()[0]
        txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
        at.text(0.05, 0.8, txt, transform=at.transAxes)

    sets = (set_rh, set_t, set_nc, set_cpc)

    for s in sets:
        a, lc, _ = s
        try:
            pass
            # lc.set_linewidth(2.5)
        except:
            pass

    return f, sets

def plot_on_clouds_flightpath(tbs, cloud_vmin = 100, altitude_column = 'Altitude_iMet', mode = None):
    """

    Parameters
    ----------
    tbs
    altitude_column
    mode: str ([None], 'overview')

    Returns
    -------

    """

    bsts = tbs.ceilometer.backscatter.get_timespan()
    tsts = tbs.data_ts.get_timespan()
    if not ((bsts[1] > tsts[0]) and (tsts[1] > bsts[0])):
        raise ValueError('there is no overlap of the tbs and the ceilometer data!')

    f,a = _plt.subplots()
    if tbs.kazar:
        tbs.kazr.reflectivity.plot(snr_max=10, ax=a)
        tbs.ceilometer.cloudbase.plot(ax=a)
        g = a.get_lines()[-1]
        g.set_markersize(5)
        g.set_markeredgewidth(1.5)
    else:
        tbs.ceilometer.backscatter.plot(ax = a, vmin=cloud_vmin)
    # tbs.data_ts.data[altitude_column].plot(ax = a)




    a.plot(tbs.data_ts.data.index, tbs.data_ts.data[altitude_column])
    a.xaxis.set_major_formatter(_DateFormatter("%H:%M:%S"))
    a.set_xlabel('')

    if mode == 'overview':
        f.set_figheight(2)
        a.set_ylim((tbs.data_ts.data[altitude_column].min() - 10, tbs.data_ts.data[altitude_column].max() + 50))
        st = tbs.data_ts.get_timespan()[0]
        txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
        a.text(0.05, 0.85, txt, transform=a.transAxes)
        f.save = lambda: save_figure(tbs, f, 'plot_on_clouds_flightpath_overview')
        a.yaxis.set_major_locator(_MaxNLocator(6, prune='both'))

    else:
        timestamp = True
        if timestamp:
            at = a
            st = tbs.data_ts.get_timespan()[0]
            txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
            at.text(0.05, 0.9, txt, transform=at.transAxes)
        f.save = lambda: save_figure(tbs, f, 'plot_on_clouds_flightpath')

    f.patch.set_alpha(0)
    f.tight_layout()
    return f,a


class Sections(object):
    def __init__(self,parent):
#         super().__init__(sections)
        self._parent = parent
        self._sections = None

    @property
    def sections_dict(self):
        return self._sections

    @sections_dict.setter
    def sections_dict(self, sections):
        # self.sections = sections
        self._sections = sections
        self._split_into_sections()

    def _split_into_sections(self):
        for sec in self._sections.keys():
            zt = self._sections[sec]
            sec_data = self._parent.data_ts.zoom_time(zt[0], zt[1])

            if self._parent.dist_dry:
                try:
                    dist_dry = self._parent.dist_dry.zoom_time(zt[0], zt[1])
                except IndexError:
                    dist_dry = None
            else:
                dist_dry = None

            if self._parent.dist_wet:
                try:
                    dist_wet = self._parent.dist_wet.zoom_time(zt[0], zt[1])
                except IndexError:
                    dist_wet = None
            else:
                dist_wet = None

            if self._parent.dist_uhsas:
                try:
                    dist_uhsas = self._parent.dist_uhsas.zoom_time(zt[0], zt[1])
                except IndexError:
                    dist_uhsas = None
            else:
                dist_uhsas = None
            setattr(self, sec, TBS_flight(data_ts = sec_data, dist_dry = dist_dry, dist_wet = dist_wet, dist_uhsas = dist_uhsas))

    def plot(self, show_clouds = True, exclude = [], altitude_column = 'Altitude_iMet', timestamp = True):
        f, a = _plt.subplots()

        maxs = []
        mins = []
        if show_clouds:
            self._parent.ceilometer.backscatter.plot(ax = a, zorder = 0)
        for e,att in enumerate(dir(self)):
            if att[0] == '_' or att in ['sections_dict','plot', 'plot_avg_sd']:
                continue
            if att in exclude:
                continue
            sect = getattr(self, att)
            # sect.data_ts.data[altitude_column].plot(ax=a, label=att)
            a.plot(sect.data_ts.data.index.values, sect.data_ts.data[altitude_column].values, label = att, zorder = e+1)
            a.legend(loc = 1)
            maxs.append(sect.data_ts.data[altitude_column].max())
            mins.append(sect.data_ts.data[altitude_column].min())

        g, = a.plot(self._parent.data_ts.data.index, self._parent.data_ts.data[altitude_column], lw=0.8, ls=':', zorder = 1)
        g.set_color('black')
        # g.set_zorder(0)
        if timestamp:
            at = a
            st = self._parent.data_ts.get_timespan()[0]
            txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
            at.text(0.05, 0.9, txt, transform=at.transAxes)

        min = _np.array(mins).min()
        max = _np.array(maxs).max()
        range = max - min
        toleranz = range * 0.05
        a.set_ylim((min - toleranz,  max + toleranz))
        a.xaxis.set_major_formatter(_DateFormatter('%H:%M:%S'))
        a.set_xlabel('')
        f.patch.set_alpha(0)
        f.autofmt_xdate()
        try:
            f.tight_layout()
        except:
            pass
        f.save = lambda x = 'sections': save_figure(self._parent, f, x)
        return f,a

    def plot_avg_sd(self, which_pops = 'wet', exclude = [],  timestamp = True):
        f, a = _plt.subplots()
        maxs = []
        mins = []

        for att in dir(self):
            if att[0] == '_' or att in ['sections_dict','plot', 'plot_avg_sd']:
                continue
            if att in exclude:
                continue
            sect = getattr(self, att)
            # sect.data_ts.data[altitude_column].plot(ax=a, label=att)
            if which_pops == 'wet':
                if sect.dist_wet:
                    avg = sect.dist_wet.average_overAllTime()
                else:
                    raise ValueError('There does not seam to be data from wet POPS')
            elif which_pops == 'dry':
                if sect.dist_dry:
                    avg = sect.dist_dry.average_overAllTime()
                else:
                    raise ValueError('There does not seam to be data from dry POPS')
            avg = avg.convert2dNdlogDp()
            avg.plot(ax=a, label=att)
            a.legend(loc = 1)
            maxs.append(avg.data.iloc[0,:].values.max())
            values = avg.data.iloc[0,:].values
            values = values[values > 0]
            mins.append(values.min())


        if timestamp:
            at = a
            st = self._parent.data_ts.get_timespan()[0]
            txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
            at.text(0.05, 0.9, txt, transform=at.transAxes)

        min = _np.array(mins).min()
        max = _np.array(maxs).max()
        # range = max - min
        # toleranz = range * 0.05
        # print(min, mins)
        a.set_ylim(min*0.9,max * 1.1)
        # a.set_ylim(top =  max + toleranz)
        a.set_yscale('log')
        # a.xaxis.set_major_formatter(_DateFormatter('%H:%M:%S'))
        # a.set_xlabel('')
        f.patch.set_alpha(0)
        # f.autofmt_xdate()
        a.set_xlim(130, 3500)
        try:
            f.tight_layout()
        except:
            pass
        f.save = lambda x = 'parking_positions_sd': save_figure(self._parent, f, x)
        return f,a

class TBS_flight(object):
    """
    - read data files
    - tools to determine time offsets
    - tools to determine altitude offsets
    - align timestamps
    - convert to vertical profiles
    """

    def __init__(self, data_ts = None, dist_dry = None, dist_wet = None, dist_uhsas = None):
        self.data_ts = data_ts
        if data_ts:
            if len([col for col in data_ts.data.columns if 'CPC_' in col]) > 0:
                self.cpc_raw = True
            else:
                self.cpc_raw = None
            if len([col for col in data_ts.data.columns if 'iMet_' in col]) > 0:
                self.imet_raw = True
            else:
                self.imet_raw = None

        self.dist_dry = dist_dry
        self.dist_wet = dist_wet
        self.dist_uhsas = dist_uhsas
        self.data_vp = None
        self.sections = Sections(self)
        self._resolution = None
        self._std = None


    # @property
    # def sections(self):
    #     if not self._sections
    #         self._sections = Sections()


    def adjust_time_offset(self, which, offset_t):
        if not offset_t:
            return

        if which == 'POPS_dry':
            hk = self.hk_dry
            dist = self.dist_dry
        elif which == 'POPS_wet':
            hk = self.hk_wet
            dist = self.dist_wet
        elif which == 'CPC':
            hk = self.cpc_raw
            dist = None
        # hk = hk.copy()
        #         dist = dist.copy()
        hk.data.index = hk.data.index + _np.timedelta64(offset_t[0], offset_t[1])
        if dist:
            dist.data.index = dist.data.index + _np.timedelta64(offset_t[0], offset_t[1])
        return




    def align_and_merge(self,
                        pops_dry_pos_rel2iMet=None,
                        pops_wet_pos_rel2iMet=None,
                        cpc_pos_rel2iMet = None,
                        alt_imet_col='GPS',
                        lanch_and_landing = None):
        """

        Parameters
        ----------
        pops_dry_pos_rel2iMet
        pops_wet_pos_rel2iMet
        cpc_pos_rel2iMet
        alt_imet_col: {'GPS', 'PTU'}
            Which reading to use for iMet altitude

        Returns
        -------

        """

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

        if self.imet_raw:
            imet_al = self.imet_raw.align_to(data)
            imet_al.data.columns = 'iMet_' + imet_al.data.columns.values
            data = data.merge(imet_al)
            if alt_imet_col == 'PTU':
                data.data['Altitude_iMet'] = data.data['iMet_altitude (from iMet PTU) [m]']
            elif alt_imet_col == 'GPS':
                data.data['Altitude_iMet'] = data.data['iMet_GPS altitude [m]']

            data.data['Altitude_POPSdry'] = data.data['Altitude_iMet'].copy()
            if pops_dry_pos_rel2iMet:
                data.data['Altitude_POPSdry'] += pops_dry_pos_rel2iMet

            data.data['Altitude_POPSwet'] = data.data['Altitude_iMet'].copy()
            if pops_dry_pos_rel2iMet:
                data.data['Altitude_POPSwet'] += pops_wet_pos_rel2iMet

        if self.cpc_raw:
            cpc_al = self.cpc_raw.align_to(data)
            cpc_al.data.columns = 'CPC_' + cpc_al.data.columns.values
            data = data.merge(cpc_al)

            if self.imet_raw:
                data.data['Altitude_CPC'] = data.data['Altitude_iMet'].copy()
                if cpc_pos_rel2iMet:
                    data.data['Altitude_CPC'] += cpc_pos_rel2iMet




        if lanch_and_landing:
            data = data.zoom_time(lanch_and_landing[0], lanch_and_landing[1])
            self.dist_dry = self.dist_dry.zoom_time(lanch_and_landing[0], lanch_and_landing[1])
            self.dist_wet = self.dist_wet.zoom_time(lanch_and_landing[0], lanch_and_landing[1])

        try:
            alt_min = data.data.Altitude_iMet.min()
            colls = ['Altitude_POPSdry', 'Altitude_POPSwet', 'Altitude_CPC']
            for col in colls:
                try:
                    data.data.loc[data.data[col] < 0, col] = alt_min
                except KeyError:
                    pass
        except:
            _warnings.warn('The iMet data seams not to exist')
        self.data_ts = data

        if self.imet_raw:
            hk = self.data_ts._del_all_columns_but('Altitude_POPSdry')
            hk.data['Altitude'] = hk.data.Altitude_POPSdry
            self.dist_dry.housekeeping = hk

            hk = self.data_ts._del_all_columns_but('Altitude_POPSwet')
            hk.data['Altitude'] = hk.data.Altitude_POPSwet
            self.dist_wet.housekeeping = hk

        ######
        # ceilometer
        if self.ceilometer:
            ts = self.data_ts.get_timespan()
            self.ceilometer = self.ceilometer.zoom_time(ts[0], ts[1])

        return


    def create_vertical_profile(self, resolution, std = True):
        if (resolution != self._resolution) or (std != self._std):
            self._resolution = resolution
            self._std = std
            alt_cols = [col for col in self.data_ts.data.columns if 'Altitude_' in col]

            min_alt = min([self.data_ts.data[ac].min() for ac in alt_cols])
            max_alt = max([self.data_ts.data[ac].max() for ac in alt_cols])

            instruments = ['POPSdry', 'POPSwet', 'iMet', 'CPC']
            vplist = []
            stdlist = []
            alt_min = 1e10
            alt_max = -1e10

            for inst in instruments:
                colt = [col for col in self.data_ts.data.columns if inst in col]
                datat = self.data_ts._del_all_columns_but(colt)
                alt_coll = [col for col in datat.data.columns if 'Altitude_' in col]
                if len(alt_coll) != 1:
                    raise ValueError('no or multiple collumns with Altitude in it')
                alt_coll = alt_coll[0]
                vpt = datat.convert2verticalprofile(altitude_column=alt_coll, resolution=(resolution, min_alt, max_alt),
                                                    return_std=std)
                if std:
                    vpt, stdt = vpt
                    stdlist.append(stdt)
                vplist.append(vpt)
            # break

            # vp = vertical_profile(pd.DataFrame())
            vp = vplist[0]
            for vpt in vplist[1:]:
                vp = vp.merge(vpt)
            self.data_vp = vp

            std = stdlist[0]
            for stdt in stdlist[1:]:
                std = std.merge(stdt)

            self.data_vp_std = std
            try:
                self.dist_wet_vp = self.dist_wet.convert2verticalprofile(layer_thickness=resolution)
            except AttributeError or ValueError:
                self.dist_wet_vp = None
            try:
                self.dist_dry_vp = self.dist_dry.convert2verticalprofile(layer_thickness=resolution)
            except AttributeError or ValueError:
                self.dist_dry_vp = None


            # self._resolution = resolution
            # self._show_std = std
            # self.data_vp = self.data_ts.convert2verticalprofile(resolution=resolution, return_std=std)
            # if std:
            #     self.data_vp, self.data_vp_std = self.data_vp



    def plot_avg_dist_at_altitude(self, compare = 'altitudes', which='both', layer_list=[], moment='number'):
        """

        Parameters
        ----------
        compare: ['altitudes, instruments']
            Shell each plot show the different altitudes or the different instruments
        which: ['dry', 'wet', 'both']
        layer_list:
            e.g. [(800, 1000),(400, 600), (100, 400)]
        moment

        Returns
        -------

        """

        def plot_sub(which, ax = None):
            if which == 'dry':
                dist_vp = self.dist_dry_vp
            elif which == 'wet':
                dist_vp = self.dist_wet_vp
            if not ax:
                f, a = _plt.subplots()
            else:
                a = ax
            txt_form = '{} - {} m'
            for lay in layer_list:
                distw_layer = dist_vp.zoom_altitude(*lay)
                if moment == 'number':
                    distw_layer = distw_layer.convert2dNdlogDp()
                elif moment == 'surface':
                    distw_layer = distw_layer.convert2dSdlogDp()
                if moment == 'volume':
                    distw_layer = distw_layer.convert2dVdlogDp()
                distw_layer_avg = distw_layer.average_overAllAltitudes()
                #     dist_avg_list.append(distw_layer_avg)
                distw_layer_avg.plot(ax=a, label=txt_form.format(*lay))
            a.legend()
            a.set_yscale('log')
            a.set_xlim((120, 3500))
            return a

        if compare == 'altitudes':
            if which == 'both':
                f,a = _plt.subplots(2, sharex= True, gridspec_kw = {'hspace':0})
                a_dry, a_wet = a
                plot_sub('dry', ax = a_dry)
                plot_sub('wet', ax = a_wet)
                a_dry.set_xlabel('')
                a[1].legend(loc=9)
                lg = a[0].legend()
                lg.remove()
                txt = 'POPS dry'
                at = a[0]
                at.text(0.95, 0.9, txt,
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=at.transAxes)

                txt = 'POPS wet'
                at = a[1]
                at.text(0.95, 0.9, txt,
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=at.transAxes)
            else:
                a = plot_sub(which)
        elif compare == 'instruments':
            f,a = _plt.subplots(len(layer_list), sharex= True, gridspec_kw = {'hspace':0})
            dist_list = []
            for e,lay in enumerate(layer_list):
                distw_layer_dry = self.dist_dry_vp.zoom_altitude(*lay)
                distw_layer_wet = self.dist_wet_vp.zoom_altitude(*lay)
                if moment == 'number':
                    distw_layer_dry = distw_layer_dry.convert2dNdlogDp()
                    distw_layer_wet = distw_layer_wet.convert2dNdlogDp()
                elif moment == 'surface':
                    distw_layer_dry = distw_layer_dry.convert2dSdlogDp()
                    distw_layer_wet = distw_layer_wet.convert2dSdlogDp()
                if moment == 'volume':
                    distw_layer_dry = distw_layer_dry.convert2dVdlogDp()
                    distw_layer_wet = distw_layer_wet.convert2dVdlogDp()
                distw_layer_avg_dry = distw_layer_dry.average_overAllAltitudes()
                distw_layer_avg_wet = distw_layer_wet.average_overAllAltitudes()
                distw_layer_avg_dry.plot(ax = a[e], label = 'dry')
                distw_layer_avg_wet.plot(ax = a[e], label = 'wet')
                txt = 'altitude: {}-{}m'.format(*lay)
                at = a[e]
                at.text(0.95, 0.9, txt,
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=at.transAxes)
                at.legend(loc = 4)
                dist_list.append((distw_layer_dry, distw_layer_wet))


        return a#, dist_list


    def plot2check_timeoffset(self, which, offset_t = (0 ,'s'), offset_alt = 0, xylim = (None, None), ylim_twin = None):
        """
        Arguments
        ----------
        which: str, ['POPS_dry', POPS_wet', 'CPC']
        """
        f ,a = _plt.subplots()
        self.imet_raw.data['altitude (from iMet PTU) [m]'].plot(ax = a, label = 'imet (PTU)')
        self.imet_raw.data['GPS altitude [m]'].plot(ax = a, label = 'imet (GPS)')
        twin = False

        if which == 'POPS_dry':
            hkw = self.hk_dry
            col = 'Altitude'
        elif which == 'POPS_wet':
            hkw = self.hk_wet
            col = 'Altitude'
        elif which == 'CPC':
            twin = True
            col = 'Concentration (#/cm³)'
            hkw = self.cpc_raw

        if twin:
            at = a.twinx()
        else:
            at = a

        at.plot(hkw.data.index + _np.timedelta64(offset_t[0], offset_t[1]), hkw.data[col] + offset_alt, label = which)
        a.set_xlim(xylim[0])
        a.set_ylim(xylim[1])
        a.legend()
        if twin:
            at.set_ylim(ylim_twin)
            at.legend()
            at.set_yscale('log')
            return a, at
        else:
            return a

    def plot4sectioning(self, column = 'Altitude_iMet', show_clouds = False):
        """Plots the altitude time series (of one of the popses right now). If interactive plotting mode is used
        klicking on the plot will print the position and save it in the axis as "plot_click_positions"

        Example
        -------
        %matplotlib nbagg

        a, ts = self.plot4sectioning()

        for click in ts.plot_click_positions:
            print('{}'.format(click[0]))

        """
        f,a = _plt.subplots()
        if show_clouds:
            self.ceilometer.backscatter.plot(ax = a)
        dalt = self.data_ts._del_all_columns_but(column)
        a = dalt.plot(ax = a,picker=5)
        return a, dalt

    plot_on_clouds = plot_on_clouds
    plot_on_clouds_flightpath = plot_on_clouds_flightpath


    def plot_overview_ts(self, alt_lim = (None, None), cpc_lim = (None, None), figh_hight_scale = 1.5):

        f ,a = _plt.subplots(5, sharex=True, gridspec_kw={'hspace': 0})
        f.set_figheight(f.get_figheight() * figh_hight_scale)
        a_alt = a[0]
        a_nc = a[2]
        a_rh = a[3]
        a_t = a[1]
        a_cpc = a[4]

        at = a[0]
        # at.set_title('Overview time series')

        if 1:
            # cpc
            if self.cpc_raw:
                self.data_ts.data['CPC_Concentration (#/cm³)'].plot(ax = a_cpc)
                a_cpc.set_yscale('log')
                a_cpc.set_ylim(cpc_lim)
            a_cpc.set_ylabel('CPC (#/cm^3)')


            # altitude
            if self.imet_raw:
                self.data_ts.data.Altitude_iMet.plot(ax = a_alt, label ='iMet')
            self.data_ts.data['POPSdry_Altitude'].plot(ax = a_alt, label ='POPS_dry')
            self.data_ts.data['POPSwet_Altitude'].plot(ax = a_alt, label ='POPS_wet')
            a_alt.set_ylabel('Alt. (m)')
            a_alt.set_ylim(alt_lim)
            a_alt.legend(loc = 1)

        # particle number
        if self.dist_dry:
           self.data_ts.data.POPSdry_PartCon_fromsizedist.plot(ax=a_nc, label='dry')
        if self.dist_wet:
            self.data_ts.data.POPSwet_PartCon_fromsizedist.plot(ax=a_nc, label='wet')
        # self.data_ts.data.POPSdry_PartCon.plot(ax = a_nc, label ='dry')
        # self.data_ts.data.POPSwet_PartCon.plot(ax = a_nc, label ='wet')
        a_nc.set_ylabel('POPS NC (#/cm^3)')
        a_nc.set_yscale('log')
        a_nc.legend()

        # temperatur
        if self.imet_raw:
            self.data_ts.data['iMet_iMet air temperature (corrected) [deg C]'].plot(ax = a_t, label ='iMet')
        self.data_ts.data.POPSdry_Temp.plot(ax = a_t, label ='POPS_dry')
        self.data_ts.data.POPSwet_Temp.plot(ax = a_t, label ='POPS_wet')
        a_t.set_ylabel('Temp (°C)')
        a_t.legend()

        # RH
        if self.imet_raw:
            self.data_ts.data['iMet_iMet humidity [RH %]'].plot(ax = a_rh)
        a_rh.set_ylabel('RH (%)')

        # other stuff
        for at in a:
            scale = at.get_yscale()
            if scale == 'linear':
                at.yaxis.set_major_locator(_MaxNLocator(prune='both', nbins=5))
        f.patch.set_alpha(0)
        f.tight_layout()
        f.save = lambda: save_figure(self, f, 'plot_overview_ts')

        timestamp = True
        if timestamp:
            at = a[0]
            st = self.data_ts.get_timespan()[0]
            txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
            at.text(0.05, 0.8, txt, transform=at.transAxes)

        return f,a


    def plot_overview_vp(self, resolution=100,
                               show_std = False,
                               # alt_lim=(None, None),
                               # nc_lim=(None, None),
                               # nc_scale = 'log',
                               # temp_lim=(None, None),
                               # rh_lim=(None, None),
                               # md_lim=(None, None),
                               # cpc_lim=(None, None),
                               # cpc_scale = 'log',
                               fighwidthscale = 1.5
                               ):
        """

        Parameters
        ----------
        resolution
        show_std
        fighwidthscale

        Returns
        -------
        a_t,a_rh,a_nc,a_md,a_cpc
        """

        def plot_ibx(alt, mean, std, scale, a, col):
            low = mean - std
            high = mean + std
            if scale == 'log':
                low[low <= 0] = low[low > 0].min()
            a.fill_betweenx(alt, low, high, color = col, alpha = 0.3)

        self.create_vertical_profile(resolution)

        f ,a = _plt.subplots(1, 5, sharey=True, gridspec_kw={'wspace': 0})
        f.set_figwidth(f.get_figwidth() * fighwidthscale)
        #     a_alt = a[0]
        a[0].set_ylabel('Altitude (m)')
        a_nc = a[2]
        a_rh = a[1]
        a_t = a[0]
        a_md = a[3]
        a_cpc = a[4]

        at = a[2]
        # at.set_title('Overview vertical profile')

        alt_data = self.data_vp.data.index

    # CPC
        if 1:
            if self.cpc_raw:
                coll = 'CPC_Concentration (#/cm³)'
                # alt_coll = 'Altitude_CPC'
                a_cpc.plot(self.data_vp.data[coll], alt_data, label = 'dry')
                g = a_cpc.get_lines()[-1]
                col = g.get_color()

                if show_std:
                    plot_ibx(alt_data, self.data_vp.data[coll], self.data_vp_std.data[coll],'log', a_cpc, col)
                    # a_cpc.fill_betweenx(alt_data.values, (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                    #                    (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color = col, alpha = 0.3)

            # a_cpc.set_xlim(cpc_lim)
            a_cpc.set_xlabel('CPC (#/cm^3)')
            a_cpc.set_xscale('log')


    # particle number
    #     coll = 'POPSdry_PartCon'
        if 1:
            if self.dist_dry:
                coll = 'POPSdry_PartCon_fromsizedist'
                # alt_coll = 'Altitude_POPSdry'
                a_nc.plot(self.data_vp.data[coll], alt_data, label = 'dry')
                g = a_nc.get_lines()[-1]
                col = g.get_color()

                if show_std:
                    plot_ibx(alt_data, self.data_vp.data[coll], self.data_vp_std.data[coll], 'log', a_nc, col)
                    # a_nc.fill_betweenx(alt_data.values, (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                    #                    (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color = col, alpha = 0.3)
            if self.dist_wet:
                # coll = 'POPSwet_PartCon'
                coll = 'POPSwet_PartCon_fromsizedist'
                # alt_coll = 'Altitude_POPSwet'
                a_nc.plot(self.data_vp.data[coll], alt_data, label = 'wet')
                g = a_nc.get_lines()[-1]
                col = g.get_color()

                if show_std:
                    plot_ibx(alt_data, self.data_vp.data[coll], self.data_vp_std.data[coll], 'log', a_nc, col)
                    # a_nc.fill_betweenx(alt_data.values, (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                    #                    (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color = col, alpha = 0.3)

            a_nc.set_xlabel('NC (#/cm^3)')
            # a_nc.set_xlim(nc_lim)
            # a_nc.set_ylim(alt_lim)
            a_nc.legend()
            a_nc.set_xscale('log')

    # temperatur
        if 1:
            coll = 'iMet_iMet air temperature (corrected) [deg C]'
            # alt_coll = 'Altitude_iMet'
            a_t.plot(self.data_vp.data[coll], alt_data)
            g = a_t.get_lines()[-1]
            col = g.get_color()

            if show_std:
                a_t.fill_betweenx(alt_data.values, (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                                   (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color = col, alpha = 0.3)

            a_t.set_xlabel('Temp (°C)')
            # a_t.legend()
            # a_t.set_xlim(temp_lim)
            # a_t.set_ylim(alt_lim)

            # RH
            coll = 'iMet_iMet humidity [RH %]'
            # alt_coll = 'Altitude_iMet'
            a_rh.plot(self.data_vp.data[coll], alt_data)
            g = a_rh.get_lines()[-1]
            col = g.get_color()

            if show_std:
                a_rh.fill_betweenx(alt_data.values, (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                                   (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color = col, alpha = 0.3)

            a_rh.set_xlabel('RH (%)')
            # a_rh.set_xlim(rh_lim)

        # mean diameter
        if 1:
            if self.dist_dry:
                coll = 'POPSdry_mean d (nm)'
                # alt_coll = 'Altitude_POPSdry'
                a_md.plot(self.data_vp.data[coll], alt_data, label = 'dry')
                g = a_md.get_lines()[-1]
                col = g.get_color()

                if show_std:
                    a_md.fill_betweenx(alt_data.values, (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                                       (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color = col, alpha = 0.3)

            if self.dist_wet:
                coll = 'POPSwet_mean d (nm)'
                # alt_coll = 'Altitude_POPSwet'
                a_md.plot(self.data_vp.data[coll], alt_data, label = 'wet')
                g = a_md.get_lines()[-1]
                col = g.get_color()

                if show_std:
                    a_md.fill_betweenx(alt_data.values,
                                       (self.data_vp.data[coll] - self.data_vp_std.data[coll]).values,
                                       (self.data_vp.data[coll] + self.data_vp_std.data[coll]).values, color=col, alpha=0.3)

            a_md.set_xlabel('mean d (nm)')
            # a_md.set_xlim(md_lim)
            a_md.legend()

        # other stuff
        for at in a:
            scale = at.get_xscale()
            if scale == 'linear':
                at.xaxis.set_major_locator(_MaxNLocator(prune='both', nbins=5))
            mtl = at.xaxis.get_majorticklabels()
            _plt.setp(mtl, rotation=45, ha='right')

        f.tight_layout()
        f.patch.set_alpha(0)
        f.save = lambda: save_figure(self, f, 'plot_overview_vp')


        timestamp = True
        if timestamp:
            at = a[0]
            st = self.data_ts.get_timespan()[0]
            txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
            at.text(0.05, 0.9, txt, transform=at.transAxes)
        return f, (a_t,a_rh,a_nc,a_md,a_cpc)


    def plot_overview_vp_check(self,
                         # cpc_scale = 'log',
                         avg_time = (30 ,'s'),
                         alpha = 0.3,
                               fighwidthscale=1.5
                         ):
        """This is for plotting when replacing the index of the timeseries with the altitude column. Usefull to see if
        up and down are different. To see a vertical profile that is produced by binning the data into altitude bins
        use plot_overview_vp"""

        f ,a = _plt.subplots(1, 5, sharey=True, gridspec_kw={'wspace': 0})
        f.set_figwidth(f.get_figwidth() * fighwidthscale)
        #     a_alt = a[0]
        a[0].set_ylabel('Altitude (m)')
        a_nc = a[2]
        a_rh = a[1]
        a_t = a[0]
        a_md = a[3]
        a_cpc = a[4]

        at = a[2]
        # at.set_title('Overview vertical profile - quick check')

        data_avg = self.data_ts.average_time(avg_time)

    # CPC
        if self.cpc_raw:
            altitude = self.data_ts.data.Altitude_CPC
            altitude_avg = data_avg.data.Altitude_CPC
            a_cpc.plot(self.data_ts.data['CPC_Concentration (#/cm³)'], altitude, alpha = alpha)
            g = a_cpc.get_lines()[-1]
            col1 = g.get_color()
            a_cpc.plot(data_avg.data['CPC_Concentration (#/cm³)'], altitude_avg, color = col1)
        # a_cpc.set_xlim(cpc_lim)
        a_cpc.set_xlabel('CPC (#/cm^3)')
        a_cpc.set_xscale('log')


    # particle number
        column_dry = 'POPSdry_PartCon_fromsizedist' #'POPSdry_PartCon'
        column_wet = 'POPSwet_PartCon_fromsizedist' #'POPSwet_PartCon'
        altitude = self.data_ts.data.Altitude_POPSdry
        altitude_avg = data_avg.data.Altitude_POPSdry
        altitudeII = self.data_ts.data.Altitude_POPSwet
        altitude_avgII = data_avg.data.Altitude_POPSwet
        if self.dist_dry:
            a_nc.plot(self.data_ts.data[column_dry], altitude, alpha = alpha)
            g = a_nc.get_lines()[-1]
            g.set_label(None)
            col1 = g.get_color()
        if self.dist_wet:
            a_nc.plot(self.data_ts.data[column_wet], altitudeII, alpha = alpha)
            g = a_nc.get_lines()[-1]
            g.set_label(None)
            col2 = g.get_color()
        if self.dist_dry:
            a_nc.plot(data_avg.data[column_dry] ,altitude_avg, label = 'dry', color = col1)
        if self.dist_wet:
            a_nc.plot(data_avg.data[column_wet] ,altitude_avgII, label = 'wet', color = col2)
        a_nc.set_xlabel('NC (#/cm^3)')
        # a_nc.set_xlim(nc_lim)
        # a_nc.set_ylim(alt_lim)
        a_nc.legend()
        a_nc.set_xscale('log')

    # temperatur
        if self.imet_raw:
            altitude = self.data_ts.data.Altitude_iMet
            altitude_avg = data_avg.data.Altitude_iMet
            a_t.plot(self.data_ts.data['iMet_iMet air temperature (corrected) [deg C]'], altitude, label ='iMet', alpha = alpha)
            g = a_t.get_lines()[-1]
            col1 = g.get_color()
            a_t.plot(data_avg.data['iMet_iMet air temperature (corrected) [deg C]'] ,altitude_avg, label = 'iMet', color = col1)
        a_t.set_xlabel('Temp (°C)')
        # a_t.legend()
        # a_t.set_xlim(temp_lim)
        # a_t.set_ylim(alt_lim)

    # RH
        if self.imet_raw:
            altitude = self.data_ts.data.Altitude_iMet
            altitude_avg = data_avg.data.Altitude_iMet
            col = 'iMet_iMet humidity [RH %]'
            a_rh.plot(self.data_ts.data[col], altitude, label ='iMet', alpha = alpha)
            g = a_t.get_lines()[-1]
            col1 = g.get_color()
            a_rh.plot(data_avg.data[col] ,altitude_avg, label = 'iMet', color = col1)
        a_rh.set_xlabel('RH (%)')
        # a_rh.set_xlim(rh_lim)

    # mean diameter
        altitude = self.data_ts.data.Altitude_POPSdry
        altitude_avg = data_avg.data.Altitude_POPSdry
        altitudeII = self.data_ts.data.Altitude_POPSwet
        altitude_avgII = data_avg.data.Altitude_POPSwet
        if self.dist_dry:
            a_md.plot(self.data_ts.data['POPSdry_mean d (nm)'], altitude, alpha = alpha)
            g = a_md.get_lines()[-1]
            g.set_label(None)
            col1 = g.get_color()
        if self.dist_wet:
            a_md.plot(self.data_ts.data['POPSwet_mean d (nm)'], altitudeII, alpha = alpha)
            g = a_md.get_lines()[-1]
            g.set_label(None)
            col2 = g.get_color()
        if self.dist_dry:
            a_md.plot(data_avg.data['POPSdry_mean d (nm)'] ,altitude_avg, label = 'dry', color = col1)
        if self.dist_wet:
            a_md.plot(data_avg.data['POPSwet_mean d (nm)'] ,altitude_avgII, label = 'wet', color = col2)

        a_md.set_xlabel('mean d (nm)')
        # a_md.set_xlim(md_lim)
        a_md.legend()

        for at in a:
            scale = at.get_xscale()
            if scale == 'linear':
                at.xaxis.set_major_locator(_MaxNLocator(prune='both', nbins=5))
            mtl = at.xaxis.get_majorticklabels()
            _plt.setp(mtl, rotation=45, ha='right')
        f.tight_layout()
        f.patch.set_alpha(0)
        f.save = lambda: save_figure(self, f, 'plot_overview_vp_check')

        timestamp = True
        if timestamp:
            at = a[0]
            st = self.data_ts.get_timespan()[0]
            txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
            at.text(0.05, 0.9, txt, transform=at.transAxes)

        return f,a

    def plot_instrument_intercomparison(self, clim = (1e0, 5e2), show_dist_ts = True, show_nc = True, show_dist_avg = True):
        def plot_sdist_ts(uhsas, dist_dry, dist_wet, clim = (1e0, 5e2)):
            f, a = _plt.subplots(1, 3, sharey=True, gridspec_kw={'wspace': 0})
            auh, ad, aw = a
            pcs = []
            a_used = []
            norm = 'linear'

            f, _, pc, cb = uhsas.plot(ax=auh, norm=norm, colorbar=False)
            pcs.append(pc)
            a_used.append(auh)
            if dist_dry:
                f, _, pc, cb = dist_dry.plot(ax=ad, norm=norm, colorbar=False)
                pcs.append(pc)
                a_used.append(ad)
            else:
                txt = 'no \ndry\nPOPS'
                ad.text(0.5, 0.5, txt, transform=ad.transAxes, fontsize='x-large', va='center', ha='center')

            if dist_wet:
                f, _, pc, cb = dist_wet.plot(ax=aw, norm=norm, colorbar=False)
                pcs.append(pc)
                a_used.append(aw)
            else:
                txt = 'no \nwet \nPOPS'
                aw.text(0.5, 0.5, txt, transform=aw.transAxes, fontsize='x-large', va='center', ha='center')

            for pc in pcs:
                pc.set_clim(clim)

            names = ['UHSAS', 'POPS_dry', 'POPS_wet']
            for e, at in enumerate(a):
                at.set_xlabel('')
                # at.xaxis.set_major_formatter(_DateFormatter('%H:%M:%S'))
                # at.xaxis.set_major_locator(_MaxNLocator(5, prune='both'))
                at.text(0.5, 0.95, names[e], transform=at.transAxes, ha='center')
                at.set_ylim(1e2, 4.8e3)

            for e, at in enumerate(a_used):
                # at.set_xlabel('')
                at.xaxis.set_major_formatter(_DateFormatter('%H:%M:%S'))
                at.xaxis.set_major_locator(_MaxNLocator(5, prune='both'))
                # at.text(0.5, 0.95, names[e], transform=at.transAxes, ha='center')
                # at.set_ylim(1e2, 4.8e3)

            f.autofmt_xdate()

            cb, cax, abg = _plt_tools.colorbar.colorbar_inside_plot(a[0], pc, loc=1, extend=('95%', '30%'),
                                                                   extend_cb=('95%', '15%'),
                                                                   colorbar_kw={'orientation': 'horizontal'},
                                                                   color_bg=[1, 1, 1, 0])
            cb.set_label('dN/dlogDp')
            cb.locator = _MaxNLocator(4, prune='both')
            cb.update_ticks()
            for at in a[1:]:
                at.set_ylabel('')

            timestamp = True
            if timestamp:
                at = a[0]
                st = self.data_ts.get_timespan()[0]
                txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
                at.text(0.05, 0.9, txt, transform=at.transAxes)

            # f.tight_layout()
            f.patch.set_alpha(0)
            f.save = lambda x = 'instrument_interc_sdist_ts': save_figure(self, f, x)
            return f, a

        def plot_particle_conc(uhsas, dist_dry, dist_wet):
            binmin = []
            binmax = []
            for dis in [dist_dry, dist_dry, uhsas]:
                try:
                    binmin.append(dis.bins.min())
                    binmax.append(dis.bins.max())
                except:
                    continue

            binminmax = _np.array(binmin).max()
            binmaxmin = _np.array(binmax).min()

            # binminmax = _np.array([dist_dry.bins.min(), dist_dry.bins.min(), uhsas.bins.min()]).max()
            # binmaxmin = _np.array([dist_dry.bins.max(), dist_dry.bins.max(), uhsas.bins.max()]).min()
            if dist_dry:
                dist_dry_minax = dist_dry.zoom_diameter(binminmax, binmaxmin)
            if dist_wet:
                dist_wet_minax = dist_wet.zoom_diameter(binminmax, binmaxmin)
            uhsas_minax = uhsas.zoom_diameter(binminmax, binmaxmin)

            if dist_dry:
                ddmmpnal = dist_dry_minax.particle_number_concentration.align_to(uhsas_minax.particle_number_concentration)
            if dist_wet:
                dwmmpnal = dist_wet_minax.particle_number_concentration.align_to(uhsas_minax.particle_number_concentration)

            if dist_dry:
                ncratio_uh_d = uhsas_minax.particle_number_concentration / ddmmpnal
            if dist_wet:
                ncratio_uh_w = uhsas_minax.particle_number_concentration / dwmmpnal
            if dist_dry and dist_wet:
                ncratio_d_w = ddmmpnal / dwmmpnal

            f, a = _plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
            anc, ar = a
            ccirc = _plt.rcParams['axes.prop_cycle'].by_key()['color']

            # numberconentrations
            uhsas_minax.particle_number_concentration.plot(ax=anc, label='UHSAS', zorder=20)
            if dist_dry:
                dist_dry_minax.particle_number_concentration.plot(ax=anc, label='POPS_dry')
            if dist_wet:
                dist_wet_minax.particle_number_concentration.plot(ax=anc, label='POPS_wet')
            anc.legend()
            # anc.set_ylim(60, 180)
            anc.set_ylabel('NC #/$cm^3$')
            # ratios
            if dist_dry:
                ncratio_uh_d.plot(ax=ar, label='UH/P_dry ($\mu$ = {:0.2f})'.format(
                    ncratio_uh_d.data['Particle number concentration #/$cm^3$'].mean()), color=ccirc[3])
            if dist_wet:
                ncratio_uh_w.plot(ax=ar, label='UH/P_wet ($\mu$ = {:0.2f})'.format(
                    ncratio_uh_w.data['Particle number concentration #/$cm^3$'].mean()), color=ccirc[4])
            if dist_dry and dist_wet:
                ncratio_d_w.plot(ax=ar, label='P_dry/P_wet ($\mu$ = {:0.2f})'.format(
                    ncratio_d_w.data['Particle number concentration #/$cm^3$'].mean()), color=ccirc[5])
            ar.legend()
            ar.set_ylabel('ratios')

            timestamp = True
            if timestamp:
                at = a[0]
                st = self.data_ts.get_timespan()[0]
                txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
                at.text(0.05, 0.9, txt, transform=at.transAxes)

            at = a[-1]
            at.xaxis.set_major_formatter(_DateFormatter('%H:%M:%S'))
            at.set_xlabel('')

            f.tight_layout()
            f.patch.set_alpha(0)
            f.save = lambda x = 'instrument_interc_particles_conc': save_figure(self, f, x)

            return f, a

        def plot_sdist_avg(uhsas, dist_dry, dist_wet):
            avguh = uhsas.average_overAllTime()
            if dist_dry:
                avgd = dist_dry.average_overAllTime()
            if dist_wet:
                avgw = dist_wet.average_overAllTime()

            try:
                avguhlim = avguh.zoom_diameter(avgd.bins.min(), avgd.bins.max())
            except:
                avguhlim = avguh.zoom_diameter(avgw.bins.min(), avgw.bins.max())

            if dist_dry:
                avgdinter = _np.interp(avguhlim.data.columns.values, avgd.data.columns.values, avgd.data.iloc[0, :].values)
            if dist_wet:
                avgwinter = _np.interp(avguhlim.data.columns.values, avgw.data.columns.values, avgw.data.iloc[0, :].values)

            if dist_dry:
                ratio_uh_d = avguhlim.data.iloc[0, :].values / avgdinter
            if dist_wet:
                ratio_uh_w = avguhlim.data.iloc[0, :].values / avgwinter
            if dist_dry and dist_wet:
                ratio_d_w = avgdinter / avgwinter

            f, a = _plt.subplots(3, sharex=True, gridspec_kw={'hspace': 0})
            alin, alog, ar = a
            ccirc = _plt.rcParams['axes.prop_cycle'].by_key()['color']
            # linear
            avguh.plot(ax=alin, label='UHSAS')
            if dist_dry:
                avgd.plot(ax=alin, label='POPS_dry')
            if dist_wet:
                avgw.plot(ax=alin, label='POPS_wet')
            alin.legend(loc=1)

            # log
            avguh.plot(ax=alog, label='UHSAS')
            if dist_dry:
                avgd.plot(ax=alog, label='POPS_dry')
            if dist_wet:
                avgw.plot(ax=alog, label='POPS_wet')

            alog.set_yscale('log')

            # ratios
            x = avguhlim.data.columns.values
            if dist_dry:
                ar.plot(x, ratio_uh_d, label='UH/P_dry ($\mu$ = {:0.2f})'.format(ratio_uh_d.mean()), color=ccirc[3])
            if dist_wet:
                ar.plot(x, ratio_uh_w, label='UH/P_wet ($\mu$ = {:0.2f})'.format(ratio_uh_w.mean()), color=ccirc[4])
            if dist_dry and dist_wet:
                ar.plot(x, ratio_d_w, label='P_dry/P_wet ($\mu$ = {:0.2f})'.format(ratio_d_w.mean()), color=ccirc[5])
            ar.legend(fontsize='x-small')
            for at in a:
                at.set_xlim(1e2, 3.5e3)

            timestamp = True
            if timestamp:
                at = a[0]
                st = self.data_ts.get_timespan()[0]
                txt = "{}{:02d}{:02d}".format(st.year, st.month, st.day)
                at.text(0.05, 0.9, txt, transform=at.transAxes)

            f.tight_layout()
            f.patch.set_alpha(0)
            f.save = lambda x = 'instrument_interc_sdist_avg': save_figure(self, f, x)

            return f, a

        if len(self.sections.sections_dict.keys()) != 1:
            raise KeyError('You need to define exactly one section to make this work. You currently have {} sections with the names {}.'.format(len(self.sections.sections_dict.keys()), self.sections.sections_dict.keys()))


        sect = getattr(self.sections, list(self.sections.sections_dict.keys())[0])

        if sect.dist_uhsas.data.shape[0] == 0:
            raise ValueError(
                'No UHSAS data in this section. Did the measurement run into the next day? If load next day too.')

        fs = []
        axs = []

        uhsas = sect.dist_uhsas.convert2dNdlogDp()
        if sect.dist_dry:
            dist_dry = sect.dist_dry.convert2dNdlogDp()
        else:
            dist_dry = False

        if sect.dist_wet:
            dist_wet = sect.dist_wet.convert2dNdlogDp()
        else:
            dist_wet = None

        if show_dist_ts:
            f, a = plot_sdist_ts(uhsas, dist_dry, dist_wet, clim = clim)
            fs.append(f)
            axs.append(a)
        if show_nc:
            f, a = plot_particle_conc(uhsas, dist_dry, dist_wet)
            fs.append(f)
            axs.append(a)
        if show_dist_avg:
            f, a = plot_sdist_avg(uhsas, dist_dry, dist_wet)
            fs.append(f)
            axs.append(a)
        return fs, axs





    def read_flight_data_qc(self,
                            fname_tbs = None, #'/Users/htelg/data/2017_ICARUS/Fan_qc_controled_email_170905/oli.tbs.2017_iop1_iop2/oli.tbs.20170523.172028.txt',
                            fname_sizedist_wet = None, #'/Users/htelg/data/2017_ICARUS/Fan_qc_controled_email_170905/SD/20170523_SD_SN14.txt',
                            fname_sizedist_dry = None, #'/Users/htelg/data/2017_ICARUS/Fan_qc_controled_email_170905/SD/20170523_SD_SN18.txt',
                            fname_ceilometer = None, #'/Volumes/HTelg_4TB_Backup/arm_data/OLI/ceilometer/oliceilM1.b1.20170523.000009.nc',
                            fname_kazr = None,
                            fname_uhsas = None,
                            pops_dry_pos_rel2iMet=0,
                            pops_wet_pos_rel2iMet=0,
                            cpc_pos_rel2iMet=0,
                            remove_artefacts = None):

        """This is data that has been preprocessed by Fan"""
        self.imet_raw = True
        self.data_ts = read_tbs_qc(fname_tbs)
        self.data_ts.data['Altitude_iMet'] = self.data_ts.data['Alt'] * 1000
        if remove_artefacts:
            _ = self.data_ts.remove_artefacts('Altitude_iMet', inplace=True)
            self.data_ts.data.loc[self.data_ts.data['Altitude_iMet'] > remove_artefacts, ['Altitude_iMet']] = _np.nan

        self.data_ts.data['iMet_iMet humidity [RH %]'] = _np.nan
        self.data_ts.data['iMet_iMet air temperature (corrected) [deg C]'] = self.data_ts.data['ImetT']
        # self.data_ts.data['POPSdry_PartCon_fromsizedist'] = self.data_ts.data['NPOPS1']
        # self.data_ts.data['POPSwet_PartCon_fromsizedist'] = self.data_ts.data['NPOPS2']
        self.cpc_raw = True
        self.data_ts.data['CPC_Concentration (#/cm³)'] = self.data_ts.data['NCPC']
        self.data_ts.data['POPSdry_Altitude'] = self.data_ts.data['Altitude_iMet'] + pops_dry_pos_rel2iMet #originally that would the altitude measured by POPS not derived from the iMET!!!!
        self.data_ts.data['POPSwet_Altitude'] = self.data_ts.data['Altitude_iMet'] + pops_wet_pos_rel2iMet #originally that would the altitude measured by POPS not derived from the iMET!!!!
        self.data_ts.data['Altitude_POPSdry'] = self.data_ts.data['POPSdry_Altitude']
        self.data_ts.data['Altitude_POPSwet'] = self.data_ts.data['POPSwet_Altitude']
        self.data_ts.data['Altitude_CPC'] = self.data_ts.data['Altitude_iMet'] + cpc_pos_rel2iMet
        self.data_ts.data['POPSdry_Temp'] = self.data_ts.data['POPS1T']
        self.data_ts.data['POPSwet_Temp'] = self.data_ts.data['POPS2T']



        if fname_sizedist_dry:
            self.dist_dry = read_pops_sd_qc(fname_sizedist_dry)
            self.dist_dry.housekeeping = self.data_ts
            self.dist_dry.housekeeping.data['Altitude'] = self.dist_dry.housekeeping.data['Altitude_POPSdry']
            self.data_ts.data['POPSdry_mean d (nm)'] = self.dist_dry.particle_mean_diameter.align_to(self.data_ts).data
            self.data_ts.data['POPSdry_PartCon_fromsizedist'] = self.dist_dry.particle_number_concentration.align_to(self.data_ts).data

        if fname_sizedist_wet:
            self.dist_wet = read_pops_sd_qc(fname_sizedist_wet)
            self.dist_wet.housekeeping = self.data_ts
            self.dist_wet.housekeeping.data['Altitude'] = self.dist_wet.housekeeping.data['Altitude_POPSwet']
            self.data_ts.data['POPSwet_mean d (nm)'] = self.dist_wet.particle_mean_diameter.align_to(self.data_ts).data
            self.data_ts.data['POPSwet_PartCon_fromsizedist'] = self.dist_wet.particle_number_concentration.align_to(self.data_ts).data

        if fname_ceilometer:
            self.ceilometer = arm.read_ceilometer_nc(fname_ceilometer)
            ts = self.data_ts.get_timespan()
            self.ceilometer = self.ceilometer.zoom_time(ts[0], ts[1])

        if fname_kazr:
            self.kazr = arm.read_kazr_nc(fname_kazr)
            ts = self.data_ts.get_timespan()
            self.kazr = self.kazr.zoom_time(ts[0], ts[1])

        if fname_uhsas:
            self.dist_uhsas = arm.read_uhsas(fname_uhsas)
            ts = self.data_ts.get_timespan()
            self.dist_uhsas = self.dist_uhsas.zoom_time(ts[0], ts[1])
        else:
            self.dist_uhsas = None

    def read_flight_data_raw(self,
                             fname_pops_dry=None,
                             fname_pops_wet=None,
                             fname_imet=None,
                             fname_cpc=None,
                             fname_bins=None,
                             fname_ceilometer = None,
                             fname_uhsas=None,
                             fname_kazr = None,
                             pops_dry_t_offset=None,
                             pops_wet_t_offset=None,
                             pops_dry_pos_rel2iMet=None,
                             pops_wet_pos_rel2iMet=None,
                             cpc_t_offset = None,
                             cpc_pos_rel2iMet=None,
                             remove_artefacts=False,
                             align_and_merge=False,
                             alt_imet_col = 'GPS',
                             lanch_and_landing = None
                             ):
        """

        Parameters
        ----------
        fname_pops_dry
        fname_pops_wet
        fname_imet
        fname_cpc
        fname_bins
        fname_ceilometer
        pops_dry_t_offset
        pops_wet_t_offset
        pops_dry_pos_rel2iMet
        pops_wet_pos_rel2iMet
        cpc_pos_rel2iMet
        remove_artefacts: bool,int,float
            removes artefacts in the iMet altitudes. If number, all values that are larger than this will be replaced by nan.
        align_and_merge
        alt_imet_col
        lanch_and_landing:
            e.g. lanch_and_landing=('2017-05-23 18:02:44.859767', '2017-05-23 22:14:07')

        Returns
        -------

        """
        # POPS
        ## bins
        bin_edges = _pd.read_csv(fname_bins).columns.values.astype(float)

        ## wet
        self.hk_wet, self.dist_wet = read_pops_raw(fname_pops_wet, bin_edges)
        self.adjust_time_offset('POPS_wet', pops_wet_t_offset)

        ### align dist and hk and get the aligned hk
        self.dist_wet.housekeeping = self.hk_wet
        self.hk_wet = self.dist_wet.housekeeping

        ### add some extra values
        self.hk_wet.data['PartCon_fromsizedist'] = self.dist_wet.particle_number_concentration.data
        # self.hk_wet.data['PartCon_outsiderange'] = self.dist_wet.particle_number_concentration_outside_range.data

        ## dry
        self.hk_dry, self.dist_dry = read_pops_raw(fname_pops_dry, bin_edges)
        self.adjust_time_offset('POPS_dry', pops_dry_t_offset)

        ### align dist and hk and get the aligned hk
        self.dist_dry.housekeeping = self.hk_dry
        self.hk_dry = self.dist_dry.housekeeping

        ### add some extra values
        self.hk_dry.data['PartCon_fromsizedist'] = self.dist_dry.particle_number_concentration.data
        # self.hk_dry.data['PartCon_outsiderange'] = self.dist_dry.particle_number_concentration_outside_range.data

        # iMet
        if fname_imet:
            for cll in [17, 19]:
                try:
                    imet = read_imet(fname_imet, column_label_line=cll)
                except KeyError:
                    imet = None
                    continue
                else:
                    break

            if not imet:
                raise ValueError('The line number where the data is supposed to start seams to be wrong.')

            imet.data['altitude (from iMet PTU) [m]'] = imet.data['altitude (from iMet PTU) [km]'] * 1000
            imet.data['GPS altitude [m]'] = imet.data['GPS altitude [km]'] * 1000

            if remove_artefacts:
                _ = imet.remove_artefacts('altitude (from iMet PTU) [m]', inplace=True)
                imet.data.loc[imet.data['GPS altitude [m]'] > remove_artefacts, ['GPS altitude [m]']] = _np.nan

            self.imet_raw = imet
        else:
            self.imet_raw = None

        # CPC
        if fname_cpc:
            self.cpc_raw = read_cpc(fname_cpc)
            self.adjust_time_offset('CPC', cpc_t_offset)
        else:
            self.cpc_raw = None

        # Ceilometer
        if fname_ceilometer:
            self.ceilometer = arm.read_ceilometer_nc(fname_ceilometer)
        else:
            self.ceilometer = None

        if align_and_merge:
            self.align_and_merge(pops_dry_pos_rel2iMet=pops_dry_pos_rel2iMet,
                                 pops_wet_pos_rel2iMet=pops_wet_pos_rel2iMet,
                                 cpc_pos_rel2iMet=cpc_pos_rel2iMet,
                                 alt_imet_col=alt_imet_col,
                                 lanch_and_landing=lanch_and_landing)


        if fname_uhsas:
            self.dist_uhsas = arm.read_uhsas(fname_uhsas)
            ts = self.data_ts.get_timespan()
            self.dist_uhsas = self.dist_uhsas.zoom_time(ts[0], ts[1])
        else:
            self.dist_uhsas = None

        if fname_kazr:
            self.kazr = arm.read_kazr_nc(fname_kazr)
            ts = self.data_ts.get_timespan()
            self.kazr = self.kazr.zoom_time(ts[0], ts[1])
        else:
            self.kazr = None


