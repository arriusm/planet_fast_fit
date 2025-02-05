#!/Users/arri/miniconda3/bin/python

#%matplotlib

version='Vers. 2025-02-05 Arno Riffeser (arri@usm.lmu.de)\npython planet_fast_fit.py'


import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.ticker as tkr
# from   matplotlib         import gridspec
# from   matplotlib.widgets import Button, RadioButtons
from   matplotlib.widgets import Slider, RangeSlider
from   astropy.io import fits
import argparse
import batman  # see   https://lkreidberg.github.io/batman/docs/html/quickstart.html
import radvel  # see   https://radvel.readthedocs.io/en/latest/tutorials/CustomModel-tutorial.html


# https://docs.astropy.org/en/latest/constants/index.html
G         = 6.6743e-11    # m^3/(s^2*kg)
Msun      = 1.98840987e30 # kg
Mjup      = 1.8981246e27  # kg
Mearth    = 5.97216787e24 # kg
Rsun      = 695700000     # m
Rjup      = 71492000      # m
Rearth    = 6378100       # m
AU        = 1.49597871e11 # m
# rhosun  = Msun / (4./3.*np.pi*Rsun**3)


# model parameters:
#  t0
#  period
#  a/Rstar
#  rp
#  ecc
#  w
#  u1
#  u2
#  inc
#  b
#  rhostar
#  rhoplan
#  Mstar
#  Mplan
#  Rstar
#  Rplan
#  bg
#  F0
#  K
#  offset


###########################################################################################################

def get_TESS_data(lc_filename) :
    data = fits.getdata(lc_filename)
    colnames = data.names
    #print(colnames)
    if 'PDCSAP_FLUX' in colnames :
        fluxname = 'PDCSAP_FLUX'
    elif 'KSPSAP_FLUX' in colnames :
        fluxname = 'KSPSAP_FLUX'
    elif 'DET_FLUX' in colnames :
        fluxname = 'DET_FLUX'
    else :
        return np.array([0.,1.,2.]),np.array([1.,1.,1.]),np.array([1.,1.,1.])
    # print("  fluxname = ",fluxname)
    #print(data['TIME'].ndim, data[fluxname].ndim, data[fluxname+'_ERR'].ndim)
    idx = np.where( (data[fluxname]!=0.) & ( ~np.isnan(data[fluxname]) ) )[0]
    #print(data[fluxname][idx])
    return data['TIME'][idx]+7000., data[fluxname][idx], data[fluxname+'_ERR'][idx]

def get_TESS_data_tab(lc_filename) :
    tab = np.genfromtxt(lc_filename,dtype=[('H1','f8'),('H2','f8'),('H3','f8')],delimiter=",")
    tdata = np.array( [ x[0] for x in tab ] )
    Fdata = np.array( [ x[1] for x in tab ] )
    edata = np.array( [ x[2] for x in tab ] )
    idx = (np.where((Fdata!=0.)&(~np.isnan(Fdata)))[0])
    return tdata[idx], Fdata[idx], edata[idx]

def get_WST_data(lc_filename) :
    tdata, Fdata, edata  = np.genfromtxt(lc_filename,dtype="f8,f8,f8",unpack=True,comments="#")
    idx = np.where((Fdata!=0.)&(~np.isnan(Fdata)))[0]
    return tdata[idx]-2450000., Fdata[idx], edata[idx]


#############

def calc_lc(t,params,params2) :
  m = batman.TransitModel(params, t)
  f = params2.bg + params2.F0 * m.light_curve(params)
  return f

def calc_rv(t,params,params2) :
    tp = radvel.orbit.timetrans_to_timeperi( tc=params.t0, per=params.per,
                                             ecc=params.ecc, omega=params.w/180.*np.pi )
    orbel_synth = np.array([params.per, tp, params.ecc, params.w/180.*np.pi, params2.K])
    rv = radvel.kepler.rv_drive(t, orbel_synth) + params2.offset
    return rv

#def calc_orbit(phi,params) :
#    r=(params.a*(1-params.ecc**2))/(1+params.ecc*np.cos(phi+(params.w)/180.*np.pi))
#    return r

def calc_orbit_x(phi,params) :
    r=(params.a*(1-params.ecc**2))/(1+params.ecc*np.cos(phi+(params.w)/180.*np.pi))
    return r*np.cos(phi)

def calc_orbit_y(phi,params) :
    r=(params.a*(1-params.ecc**2))/(1+params.ecc*np.cos(phi+(params.w)/180.*np.pi))
    return r*np.sin(phi)


#############

def calc_xlim(t) :
    xmin = np.min(t)
    xmax = np.max(t)
    Dx = xmax-xmin
    return [xmin-0.025*Dx,xmax+0.025*Dx]

def calc_ylim(flux) :
    ymin = np.min(flux)
    ymax = np.max(flux)
    Dy = ymax-ymin
    return [ymin-0.25*Dy,ymax+0.25*Dy]

def calc_rv_xlim(t) :
    xmin_rv = np.min(t)
    xmax_rv = np.max(t)
    Dx = xmax_rv-xmin_rv
    return [xmin_rv-0.025*Dx,xmax_rv+0.025*Dx]

def calc_rv_ylim(y) :
    ymin_rv = np.min(y)
    ymax_rv = np.max(y)
    Dy = ymax_rv-ymin_rv
    return [ymin_rv-0.2*Dy,ymax_rv+0.2*Dy]

#############

def calc_b__aR_inc_ecc_w(params) :
    myb = params.a * np.cos(params.inc/180.*np.pi) * \
        (1.-params.ecc**2)/(1.+params.ecc*np.sin(params.w/180.*np.pi))
    #print('b=',myb)
    return myb

def calc_inc__b_aR_ecc_w(params,params2) :
    #print("params2.b = ",params2.b)
    #print("Nenner = ",(params.a * (1.-params.ecc**2)/(1.+params.ecc*np.sin(params.w/180.*np.pi))))
    N = params2.b
    D = params.a * (1.-params.ecc**2) / ( 1.+params.ecc*np.sin(params.w/180.*np.pi) )
    if N<=D :
        myinc = np.arccos( N / D ) / np.pi*180.
    else :
        myinc = 0.
        #print("N=",N)
        #print("D=",D)
    #print('inc=',myinc)
    return myinc

#############

# def calc_depth(params) :
# return params.rp**2 * 1000.
def calc_depth(params,params2) :
    #(x,y) = model2.get_data()
    #ymax = np.max(y)
    #ymin = np.min(y)
    ymax = params2.bg + params2.F0 
    ymin = calc_lc(np.array([params.t0]),params,params2)[0]
    return (ymax-ymin) * 1000.

def calc_dur(params,params2) :
    # tau_num = 0.
    # (x,y) = model2.get_data()
    # mm = np.max(y)
    # bad = np.where((y>mm-1e-6))
    # xtran = np.delete(x,bad)
    # if len(xtran)>=2 :
    #     tau_num = np.max(xtran)-np.min(xtran)
    tau_ana = 0.
    if (params2.b<1+params.rp) :
        delta = np.arcsin(params2.b-params.rp) # TRICK
        Rtilde = params2.Rstar*np.cos(delta)
        vttrans = (2.*np.pi*params.a*params2.Rstar) / params.per * \
            (1.+params.ecc*np.sin(params.w/180.*np.pi)) / np.sqrt((1.-params.ecc**2))
        tau_ana = 2 * Rtilde*(1+params.rp) / vttrans # plus planet radius
    # print(tau_ana, tau_num)
    return tau_ana*24.


############# rhostar a per

def calc_rhostar__aR_per(params) :
    #myrhostar = params.a**3 * (3.*np.pi) / G / (params.per*(24.*3600.))**2 / rhosun
    myrhostar = params.a**3 * (3.*np.pi) / G / (params.per*(24.*3600.))**2 / 1000.
    #print('myrhostar=',myrhostar)
    return myrhostar

def calc_aR__rhostar_per(params,params2) :
    #mya = (params2.rhostar / (3.*np.pi) * G * (params.per*(24.*3600.))**2 * rhosun)**(1./3.)
    mya = (params2.rhostar / (3.*np.pi) * G * (params.per*(24.*3600.))**2 *1000.)**(1./3.)
    #print('mya=',mya)
    return mya

############# Korrektur mit Mplan

# def calc_rhostar__rhostar_Mplan_Rstar(params2) : # Korrektur mit Mplan
#     return ( params2.rhostar*1000. + params2.Mplan*Mearth / ( 4./3.*np.pi*(params2.Rstar*Rsun)**3 ) ) / 1000.
#
# def calc_rhostar__rhostar_Mplan_Rstar(params2) :# Korrektur mit Mplan
#     return ( params2.rhostar*1000. - params2.Mplan*Mearth / ( 4./3.*np.pi*(params2.Rstar*Rsun)**3 ) ) / 1000.

############# rhostar Mstar Rstar

# def calc_rhostar__Mstar_Rstar(params2) :
#     rhostar = params2.Mstar*Msun / (4./3.*np.pi*(params2.Rstar*Rsun)**3) / 1000.
#     #print("rhostar=",rhostar)
#     return rhostar

def calc_rhostar__Mstar_Rstar(params2) :
    #print("rhostar_Mstar_Rstar")
    #print("Rstar=",params2.Rstar)
    #print("Mstar=",params2.Mstar)
    ### rhostar = (params2.Mstar*Msun) / (4./3.*np.pi*(params2.Rstar*Rsun)**3) / 1000.
    # GENAUER: rhostar = (params2.Mstar*Msun + params2.Mplan*Mearth) / (4./3.*np.pi*(params2.Rstar*Rsun)**3) / 1000.
    rhostar = (params2.Mstar*Msun) / (4./3.*np.pi*(params2.Rstar*Rsun)**3) / 1000.
    #print("rhostar=",rhostar)
    return rhostar

def calc_Mstar__rhostar_Rstar_Mplan(params2) :
    #print("calc_Mstar__rhostar_Rstar_Mplan")
    #print("rhostar=",params2.rhostar)
    #print("Rstar=",params2.Rstar)
    Mstar = (params2.rhostar*1000.) * (4./3.*np.pi*(params2.Rstar*Rsun)**3) / Msun
    # GENAUER: Mstar = ( (params2.rhostar*1000.) * (4./3.*np.pi*(params2.Rstar*Rsun)**3) - params2.Mplan*Mearth) / Msun
    #print("Mstar=",Mstar)
    return Mstar

def calc_Rstar__rhostar_Mstar_Mplan(params2) :
    #print("calc_Rstar__rhostar_Mstar_Mplan")
    #print("rhostar=",params2.rhostar)
    #print("Mstar=",params2.Mstar)
    Rstar = (( params2.Mstar * Msun ) / (params2.rhostar*1000.) / (4./3.*np.pi))**(1/3) / Rsun
    # GENAUER: Rstar = ( ( params2.Mstar*Msun + params2.Mplan*Mearth ) / (params2.rhostar*1000.) / (4./3.*np.pi) )**(1/3) / Rsun
    #print("Rstar=",Rstar)
    return Rstar

############# rhoplan Mplan Rplan

def calc_rhoplan__Mplan_Rplan(params2) :
    #print("rhoplan_Mplan_Rplan")
    #print("Rplan=",params2.Rplan)
    #print("Mplan=",params2.Mplan)
    rhoplan = (params2.Mplan*Mearth) / (4./3.*np.pi*(params2.Rplan*Rearth)**3) / 1000.
    #print("rhoplan=",rhoplan)
    return rhoplan

def calc_Mplan__rhoplan_Rplan(params2) :
    #print("calc_Mplan__rhoplan_Rplan")
    #print("rhoplan=",params2.rhoplan)
    #print("Rplan=",params2.Rplan)
    Mplan = (params2.rhoplan*1000.) * (4./3.*np.pi*(params2.Rplan*Rearth)**3) / Mearth
    #print("Mplan=",Mplan)
    return Mplan

def calc_Rplan__rhoplan_Mplan(params2) :
    #print("calc_Rplan__rhoplan_Mplan")
    #print("rhoplan=",params2.rhoplan)
    #print("Mplan=",params2.Mplan)
    Rplan = (( params2.Mplan * Mearth ) / (params2.rhoplan*1000.) / (4./3.*np.pi))**(1/3) / Rearth
    #print("Rplan=",Rplan)
    return Rplan

#############

def calc_rp__Rplan_Rstar(params,params2) :
    return (params2.Rplan*Rearth/Rsun) / params2.Rstar

def calc_Rplan__rp_Rstar(params,params2) :
    return params.rp * params2.Rstar * Rsun/Rearth

# def calc_rp__Rplan_Rstar(params,params2) :
#     return params2.Rplan / params2.Rstar
#
# def calc_Rplan__rp_Rstar(params,params2) :
#     return params.rp * params2.Rstar

def calc_K__Mplan(params,params2) :
    K  = params2.Mplan * Mearth * ( np.sin(params.inc/180.*np.pi) / ( params2.Mstar * Msun )  *  \
                                    (2*np.pi*G*params2.Mstar*Msun/(params.per*86400.))**(1./3.) * \
                                    1/(1.-params.ecc**2) )
    return K

def calc_Mplan__K(params,params2) :
    Mp = params2.K / ( np.sin(params.inc/180.*np.pi) / ( params2.Mstar * Msun )  *  \
                       (2*np.pi*G*params2.Mstar*Msun/(params.per*86400.))**(1./3.) * \
                       1/(1.-params.ecc**2) ) / Mearth
    return Mp

###########################################################################################################

def update_plot1(params,params2) :
    global myalpha
    global lc_filename
    m1 = batman.TransitModel(params, t1)
    f1 = params2.bg + params2.F0 * m1.light_curve(params)
    model1.set_ydata(f1)
    if lc_filename!='' :
        data1.set_alpha(myalpha)

def update_plot2(params,params2) :
    global myalpha
    global lc_filename
    t2 = np.linspace(params.t0-params.per/2/zoom, params.t0+params.per/2/zoom, prec)
    m2 = batman.TransitModel(params, t2)
    f2 = params2.bg + params2.F0 * m2.light_curve(params)
    model2.set_xdata(t2-params.t0)
    model2.set_ydata(f2)
    if lc_filename!='' :
        data2.set_xdata( (tdata - params.t0 + params.per/2. ) % params.per - params.per/2. )
        data2.set_alpha(myalpha)
    # newlalpha = myalpha*zoom*params.per/(tmin+tmax)
    # if newlalpha>1. :
    #     data2.set_alpha(1.)
    # elif newlalpha<0.01 :
    #     data2.set_alpha(0.01)
    # else :
    #     data2.set_alpha(newlalpha)
    ax2.set_xlim([-params.per/2 / zoom, params.per/2 / zoom])
    #ax2.set_xlim([-3,3])

def update_plot_rv_1(params,params2) :
    global myalpha_rv
    rv1 = calc_rv(t1_rv,params,params2)
    model_rv_1.set_ydata(rv1)
    if rv_filename!='' :
        data_rv_1.set_alpha(myalpha_rv)

def update_plot_rv_2(params,params2) :
    global myalpha_rv
    global rv_filename
    t2_rv = np.linspace(params.t0-params.per/2, params.t0+params.per/2, prec)
    rv2 = calc_rv(t2_rv,params,params2)
    model_rv_2.set_xdata(t2_rv-params.t0)
    model_rv_2.set_ydata(rv2)
    if rv_filename!='' :
        data_rv_2.set_xdata( (x_rv - params.t0 + params.per/2. ) % params.per - params.per/2. )
        data_rv_2.set_alpha(myalpha_rv)
    ax4.set_xlim([-params.per/2 / zoom_rv, params.per/2 / zoom_rv])

def update_transit(params,params2) :
    circ_star.set_radius(params2.Rstar)
    circ_planet.set_radius(params.rp*params2.Rstar) #### NEW
    circ_planet.set_center((0.,params2.b*params2.Rstar)) #### NEW
    #circ_planet.set_radius(params.rp)
    #circ_planet.set_center((0.,params2.b)) 
    size=1.15*(params2.Rstar)
    ax5.axis([-size,size,-size,size])

def update_orbit(params,params2) :
    #phi = orbit_p.get_xdata()
    #orbit_p.set_ydata(calc_orbit(phi,params))
    x = calc_orbit_x(phi,params)*params2.Rstar*Rsun/AU
    y = calc_orbit_y(phi,params)*params2.Rstar*Rsun/AU
    orbit_p.set_xdata(x)
    orbit_p.set_ydata(y)
    orbit_l.set_ydata([-2.*params.a*params2.Rstar*Rsun/AU,-params.rp*params2.Rstar*Rsun/AU])
    circ_star_orbit.set_radius(params2.Rstar*Rsun/AU)
    circ_planet_orbit.set_radius(params.rp*params2.Rstar*Rsun/AU)
    circ_planet_orbit.set_center((0.,calc_orbit_y(-np.pi/2.,params)*params2.Rstar*Rsun/AU))
    size=1.15*(params.a*params2.Rstar*Rsun/AU*(1-params.ecc**2))/(1-params.ecc)
    ax6.axis([-size,size,-size,size])

###########################################################################################################

# def update_sliders(params,params2) :
#     global     slider_aR,slider_rhostar,slider_Rstar,slider_Mstar,slider_rhoplan,slider_Rplan,slider_Mplan,slider_rp,slider_b,slider_inc,slider_u1,slider_u2,slider_ecc,slider_w,slider_F0,slider_bg,slider_rv_K,slider_rv_offset
#     #slider_t0.set_val(params.t0)
#     #slider_t0fine.set_val(params.t0)
#     #slider_per.set_val(params.per)
#     #slider_perfine.set_val(params.per)
#     slider_aR.set_val(params.a)
#     slider_rhostar.set_val(params2.rhostar)
#     slider_Rstar.set_val(params2.Rstar)
#     slider_Mstar.set_val(params2.Mstar)
#     slider_rhoplan.set_val(params2.rhoplan)
#     slider_Rplan.set_val(params2.Rplan)
#     slider_Mplan.set_val(params2.Mplan)
#     slider_rp.set_val(params.rp)
#     slider_b.set_val(params2.b)
#     slider_inc.set_val(params.inc)
#     slider_u1.set_val(params.u[0])
#     slider_u2.set_val(params.u[1])
#     slider_ecc.set_val(params.ecc)
#     slider_w.set_val(params.w)
#     slider_F0.set_val(params2.F0)
#     slider_bg.set_val(params2.bg)
#     slider_rv_K.set_val(params2.K)
#     slider_rv_offset.set_val(params2.offset)
#     params2.depth = calc_depth(params,params2)
#     params2.dur   = calc_dur(params,params2)
#     slider_dur.set_val(params2.dur)
#     slider_depth.set_val(params2.depth)
#     return
#
# def slider_update_xstart(val):
#     global params,params2,ax,xstart,xend,t1,m1
#     if val<xend :
#         xstart = val
#         t1 = np.linspace(xstart,xend, prec)
#         m1 = batman.TransitModel(params, t1)
#         f1 = params2.bg + params2.F0 * m1.light_curve(params)
#         model1.set_xdata(t1)
#         model1.set_ydata(f1)
#         ax1.set_xlim([xstart,xend])
#
# def slider_update_xend(val):
#     global params,params2,ax,xstart,xend,t1,m1
#     if val>xstart :
#         xend = val
#         t1 = np.linspace(xstart,xend, prec)
#         m1 = batman.TransitModel(params, t1)
#         f1 = params2.bg + params2.F0 * m1.light_curve(params)
#         model1.set_xdata(t1)
#         model1.set_ydata(f1)
#         ax1.set_xlim([xstart,xend])
#
# def slider_update_ystart(val):
#     global ax,ystart,yend,b
#     if val<yend :
#         ystart = val
#         ax1.set_ylim([ystart,yend])
#         ax2.set_ylim([ystart,yend])
#
# def slider_update_yend(val):
#     global ax,ystart,yend,b
#     if val>ystart :
#         yend = val
#         ax1.set_ylim([ystart,yend])
#         ax2.set_ylim([ystart,yend])

def slider_update_x(val) :
    global ax1,xstart,xend
    if (xstart,xend) != val :
        xstart,xend = val
        ax1.set_xlim([xstart,xend])

def slider_update_y(val) :
    global ax1,ax2,ystart,yend
    if (ystart,yend) != val :
        ystart,yend = val
        ax1.set_ylim([ystart,yend])
        ax2.set_ylim([ystart,yend])

def slider_update_zoom(val) :
    global zoom
    if zoom != val :
        zoom = val
        update_plot2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_alpha(val) :
    global myalpha
    if myalpha != val :
        myalpha= val
        update_plot1(params,params2)
        update_plot2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_x_rv(val) :
    global ax3,xstart_rv,xend_rv
    if (xstart_rv,xend_rv) != val :
        xstart_rv,xend_rv = val
        ax3.set_xlim([xstart_rv,xend_rv])

def slider_update_y_rv(val) :
    global ax3,ax4,ystart_rv,yend_rv
    if (ystart_rv,yend_rv) != val :
        ystart_rv,yend_rv = val
        ax3.set_ylim([ystart_rv,yend_rv])
        ax4.set_ylim([ystart_rv,yend_rv])

def slider_update_zoom_rv(val) :
    global zoom_rv
    if zoom_rv != val :
        zoom_rv = val
        update_plot_rv_2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_alpha_rv(val) :
    global myalpha_rv
    if myalpha_rv != val :
        myalpha_rv = val
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_t0(val) :
    global params,params2
    if params.t0 != val :
        params.t0 = val
        slider_t0.set_val(params.t0)
        slider_t0fine.set_val(params.t0)
        slider_t0fine.valmin = params.t0-0.1
        slider_t0fine.valmax = params.t0+0.1
        ax_t0fine.set_xlim( params.t0-0.1, params.t0+0.1)
        update_plot1(params,params2)
        update_plot2(params,params2)
        line1.set_xdata([params.t0,params.t0])
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        line_rv_1.set_xdata([params.t0,params.t0])
        plt.gcf().canvas.draw_idle()

def slider_update_t0fine(val) :
    global params,params2
    if params.t0 != val :
        params.t0 = val
        slider_t0.set_val(params.t0)
        slider_t0fine.set_val(params.t0)
        update_plot1(params,params2)
        update_plot2(params,params2)
        line1.set_xdata([params.t0,params.t0])
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        line_rv_1.set_xdata([params.t0,params.t0])
        plt.gcf().canvas.draw_idle()

def slider_update_per(val) :  # per -> a  fix: rhostar |  a -> inc  fix: b ecc w
    global params,params2
    if params.per != val :
        params.per = val
        slider_per.set_val(params.per)
        slider_perfine.set_val(params.per)
        slider_perfine.valmin = params.per-0.1
        slider_perfine.valmax = params.per+0.1
        ax_perfine.set_xlim( params.per-0.1, params.per+0.1)
        # rhostar fixed!!
        params.a    = calc_aR__rhostar_per(params,params2)
        params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_aR.set_val(params.a)
        slider_aAU.set_val(params2.aAU)
        slider_inc.set_val(params.inc)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_perfine(val) : # per -> a  fix: rhostar |  a -> inc  fix: b ecc w
    global params,params2
    if params.per != val :
        params.per = val
        slider_per.set_val(params.per)
        slider_perfine.set_val(params.per)
        # rhostar fixed!!
        params.a    = calc_aR__rhostar_per(params,params2)
        params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_aR.set_val(params.a)
        slider_aAU.set_val(params2.aAU)
        slider_inc.set_val(params.inc)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_rp(val) : # rp -> Rplan  fix: Rstar  |  Rplan -> rhoplan  fix: Mplan
    global params,params2
    if params.rp != val :
        rp_buf = params.rp
        params.rp  = val
        if ( params.ecc > 1 - 1/(params.a/(1+params.rp)) ) : 
            params.rp  = rp_buf
            slider_rp.set_val(params.rp)
        else :
          params.rp  = val
          params2.Rplan = calc_Rplan__rp_Rstar(params,params2)  # fix Rstar
          params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)  # fix Mplan
          slider_Rplan.set_val(params2.Rplan)
          slider_rhoplan.set_val(params2.rhoplan)
          params2.depth = calc_depth(params,params2)
          params2.dur   = calc_dur(params,params2)
          slider_dur.set_val(params2.dur)
          slider_depth.set_val(params2.depth)
          # update_sliders(params,params2)
          update_plot1(params,params2)
          update_plot2(params,params2)
          update_transit(params,params2)
          update_orbit(params,params2)
          plt.gcf().canvas.draw_idle()

def slider_update_Rstar(val) :  # 1) Rstar -> rp  fix: Rplan  |  Rstar -> Mstar    fix: rhostar |  Mstar -> Mplan    fix: K     |  Mplan -> rhoplan  fix: Rplan
                                # 2) Rstar -> rp  fix: Rplan  |  Rstar -> rhostar  fix: Mstar   |  rhostar -> a  fix: P         |  a -> inc  fix: b ecc w
                                # 3) Rstar -> Rplan  fix: rp  |  Rplan -> rhoplan  fix: Mplan   |  Rstar -> rhostar  fix: Mstar   |  rhostar -> a  fix: P         |  a -> inc  fix: b ecc w
    global params,params2
    # if params2.Rstar != val :
    #     params2.Rstar = val
    #     params.rp     = calc_rp__Rplan_Rstar(params,params2)
    #     # params2.calc_rhostar = calc_rhostar__Mstar_Rstar(params2)   # fix Mstar
    #     params2.Mstar = calc_Mstar__rhostar_Rstar_Mplan(params2)     # fix rhostar
    #     # if (params2.Mstar*Msun < params2.Mplan*Mearth) :
    #     #    params2.Mstar = params2.Mplan*Mearth/Msun
    #     if ( params2.Mstar < 0.08 ) :
    #         params2.Mstar = 0.08
    #         params2.Rstar = calc_Rstar__rhostar_Mstar_Mplan(params2)
    #         params.rp     = calc_rp__Rplan_Rstar(params,params2)
    #         slider_Rstar.set_val(params2.Rstar)
    #     params2.Mplan = calc_Mplan__K(params,params2)
    #     params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
    #     slider_Mstar.set_val(params2.Mstar)
    #     slider_rp.set_val(params.rp)
    #     slider_Mplan.set_val(params2.Mplan)
    #     slider_rhoplan.set_val(params2.rhoplan)
    #     # params2.depth = calc_depth(params,params2)
    #     # params2.dur   = calc_dur(params,params2)
    #     # slider_dur.set_val(params2.dur)
    #     # slider_depth.set_val(params2.depth)
    #     # update_sliders(params,params2)
    #     update_plot1(params,params2)
    #     update_plot2(params,params2)
    #     circ_planet.set_radius(params.rp)
    #     circ_planet.set_center((0.,params2.b))
    #     plt.gcf().canvas.draw_idle()
    if params2.Rstar != val : # 3)
        Rstar_buf = params2.Rstar
        rhostar_buf = params2.rhostar
        a_buf = params.a
        params2.Rstar = val
        params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
        params.a = calc_aR__rhostar_per(params,params2)
        if ( params.ecc > 1 - 1/(params.a/(1+params.rp)) ) :
            params2.Rstar = Rstar_buf
            params2.rhostar = rhostar_buf
            params.a = a_buf 
            slider_Rstar.set_val(params2.Rstar)
            # pos 2 : # (adjust other parameters)
            # pos 2 : params.ecc = 1 - 1/(params.a/(1+params.rp))
            # pos 2 : slider_ecc.set_val(params.ecc)
            # pos 2 : params.inc = calc_inc__b_aR_ecc_w(params,params2)
            # pos 2 : slider_inc.set_val(params.inc)
            # pos 2 : params2.K = calc_K__Mplan(params,params2)
            # pos 2 : slider_rv_K.set_val(params2.K)
        else :
            params2.Rplan   = calc_Rplan__rp_Rstar(params,params2)
            params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
            ##### params2.rhostar = calc_rhostar__aR_per(params)
            params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
            params.inc = calc_inc__b_aR_ecc_w(params,params2)
            params2.depth = calc_depth(params,params2)
            params2.dur   = calc_dur(params,params2)
            slider_dur.set_val(params2.dur)
            slider_depth.set_val(params2.depth)
            slider_Rplan.set_val(params2.Rplan)
            slider_rhoplan.set_val(params2.rhoplan)
            slider_rhostar.set_val(params2.rhostar)
            slider_aR.set_val(params.a)
            slider_aAU.set_val(params2.aAU)
            slider_inc.set_val(params.inc)
            # update_sliders(params,params2)
            update_plot1(params,params2)
            update_plot2(params,params2)
            update_plot_rv_1(params,params2)
            update_plot_rv_2(params,params2)
            update_transit(params,params2)
            update_orbit(params,params2)
            plt.gcf().canvas.draw_idle()
            
def slider_update_Rplan(val) :   #  Rplan -> rp, rhoplan  fix: Rstar,Mplan
    global params,params2
    if params2.Rplan != val :
        Rplanbuf = params2.Rplan
        rp_buf = params.rp    
        params2.Rplan   = val
        params.rp       = calc_rp__Rplan_Rstar(params,params2)
        if ( params.ecc > 1 - 1/(params.a/(1+params.rp)) ) :
            params2.Rplan = Rplanbuf
            params.rp = rp_buf
            slider_Rplan.set_val(params2.Rplan)
        else :
            params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
            slider_rp.set_val(params.rp)
            slider_rhoplan.set_val(params2.rhoplan)
            params2.depth = calc_depth(params,params2)
            params2.dur   = calc_dur(params,params2)
            slider_dur.set_val(params2.dur)
            slider_depth.set_val(params2.depth)
            # update_sliders(params,params2)
            update_plot1(params,params2)
            update_plot2(params,params2)
            update_transit(params,params2)
            update_orbit(params,params2)
            plt.gcf().canvas.draw_idle()

def slider_update_aR(val) :  # 1) a -> inc  fix: b ecc w | a -> rhostar  fix: P |  rhostar -> Mstar  fix: Rstar  | Mstar -> Mplan fix: K | Mplan -> rhoplan fix: Rplan
                            # FALSCH 2) a -> inc  fix: b ecc w | a -> rhostar  fix: P |  rhostar -> Rstar  fix: Mstar  | Rstar -> Rplan  fix: rp |  Rplan -> rhoplan  fix: Mplan
    global params,params2
    if params.a != val :
        params.a = val
        ## rhobuf = calc_rhostar__aR_per(params)
        ## if (rhobuf>10.) :
        ##     params2.rhostar = 10.
        ##     params.a = calc_aR__rhostar_per(params,params2)
        ## if ( params.ecc > 1 - 1/(params.a-params.rp) ) :
        ##     params.ecc = 1 - 1/(params.a-params.rp)
        # params2.b = calc_b__aR_inc_ecc_w(params)
        ## params.inc = calc_inc__b_aR_ecc_w(params,params2)
        ## slider_inc.set_val(params.inc)
        params2.rhostar = calc_rhostar__aR_per(params)
        params2.Mstar = calc_Mstar__rhostar_Rstar_Mplan(params2)
        # params2.rhostar  = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
        #if (params2.Mstar*Msun < params2.Mplan*Mearth) :
        #    params2.Mstar = params2.Mplan*Mearth/Msun
        if ( params2.Mstar < 0.08 ) :
            params2.Mstar = 0.08
            params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
            params.a = calc_aR__rhostar_per(params,params2)
            # params2.rhostar = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
            slider_aR.set_val(params.a)
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.Mplan = calc_Mplan__K(params,params2)
        params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
        params2.aAU = params.a*params2.Rstar*Rsun/AU
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_rhostar.set_val(params2.rhostar)
        slider_inc.set_val(params.inc)
        slider_Mstar.set_val(params2.Mstar)
        slider_Mplan.set_val(params2.Mplan)
        slider_rhoplan.set_val(params2.rhoplan)
        slider_aAU.set_val(params2.aAU)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()
    # if params.a != val :
    #     params.a = val
    #     params.inc = calc_inc__b_aR_ecc_w(params,params2)
    #     params2.rhostar = calc_rhostar__aR_per(params)
    #     params2.Rstar = calc_Rstar__rhostar_Mstar_Mplan(params2)
    #     params2.Rplan = calc_Rplan__rp_Rstar(params,params2)
    #     params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
    #     slider_inc.set_val(params.inc)
    #     slider_rhostar.set_val(params2.rhostar)
    #     slider_Rstar.set_val(params2.Rstar)
    #     slider_Rplan.set_val(params2.Rplan)
    #     slider_rhoplan.set_val(params2.rhoplan)
    #     # update_sliders(params,params2)
    #     update_plot1(params,params2)
    #     update_plot2(params,params2)
    #     circ_planet.set_radius(params.rp)
    #     circ_planet.set_center((0.,params2.b))
    #     plt.gcf().canvas.draw_idle()


def slider_update_aAU(val) :  # 1) a -> inc  fix: b ecc w | a -> rhostar  fix: P |  rhostar -> Mstar  fix: Rstar  | Mstar -> Mplan fix: K | Mplan -> rhoplan fix: Rplan
                            # 2) a -> inc  fix: b ecc w | a -> rhostar  fix: P |  rhostar -> Rstar  fix: Mstar  | Rstar -> Rplan  fix: rp |  Rplan -> rhoplan  fix: Mplan
    global params,params2
    if params2.aAU != val :
        params2.aAU = val
        params.a = params2.aAU * AU / (params2.Rstar*Rsun)
        params2.rhostar = calc_rhostar__aR_per(params)
        params2.Mstar = calc_Mstar__rhostar_Rstar_Mplan(params2)
        # params2.rhostar  = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
        #if (params2.Mstar*Msun < params2.Mplan*Mearth) :
        #    params2.Mstar = params2.Mplan*Mearth/Msun
        if ( params2.Mstar < 0.08 ) :
            params2.Mstar = 0.08
            params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
            params.a = calc_aR__rhostar_per(params,params2)
            # params2.rhostar = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
            params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
            slider_aAU.set_val(params2.aAU)
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.Mplan = calc_Mplan__K(params,params2)
        params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_rhostar.set_val(params2.rhostar)
        slider_inc.set_val(params.inc)
        slider_Mstar.set_val(params2.Mstar)
        slider_Mplan.set_val(params2.Mplan)
        slider_rhoplan.set_val(params2.rhoplan)
        slider_aR.set_val(params.a)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()


def slider_update_rhostar(val) :  # rhostar -> Mstar  fix: Rstar  |  rhostar -> a  fix: P  |  a -> inc  fix: b ecc w  |  Mstar -> Mplan fix: K | Mplan -> rhoplan fix: Rplan
                                  # rhostar -> Rstar  fix: Mstar  |  rhostar -> a  fix: P  |  a -> inc  fix: b ecc w  |  Rstar -> Rplan  fix: rp |  Rplan -> rhoplan  fix: Mplan
    global params,params2
    # if params2.rhostar != val :
    #     params2.rhostar = val
    #     # params2.rhostar = calc_rhostar__rhostar_Mplan_Rstar(params2)  # Korrektur mit Mplan
    #     params2.Mstar = calc_Mstar__rhostar_Rstar_Mplan(params2)
    #     if ( params2.Mstar < 0.08 ) :
    #         params2.Mstar = 0.08
    #         params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
    #         # params2.rhostar = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
    #         # slider_rhostar.set_val(params2.rhostar)
    #     params.a = calc_aR__rhostar_per(params,params2)
    #     params2.aAU = params.a*params2.Rstar*Rsun/AU
    #     params.inc = calc_inc__b_aR_ecc_w(params,params2)
    #     params2.Mplan = calc_Mplan__K(params,params2)
    #     params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
    #     slider_aR.set_val(params.a)
    #     slider_aAU.set_val(params2.aAU)
    #     slider_Mstar.set_val(params2.Mstar)
    #     slider_inc.set_val(params.inc)
    #     slider_rv_K.set_val(params2.K)
    #     slider_Mplan.set_val(params2.Mplan)
    #     slider_rhoplan.set_val(params2.rhoplan)
    #     # update_sliders(params,params2)
    #     update_plot1(params,params2)
    #     update_plot2(params,params2)
    #     update_plot_rv_1(params,params2)
    #     update_plot_rv_2(params,params2)
    #     circ_planet.set_radius(params.rp)
    #     circ_planet.set_center((0.,params2.b))
    #     plt.gcf().canvas.draw_idle()
    if params2.rhostar != val :
        rhostar_buf = params2.rhostar
        Rstar_buf = params2.Rstar
        Rplan_buf = params2.Rplan
        a_buf = params.a
        params2.rhostar = val
        params2.Rstar = calc_Rstar__rhostar_Mstar_Mplan(params2)
        params2.Rplan = calc_Rplan__rp_Rstar(params,params2)
        params.a = calc_aR__rhostar_per(params,params2)
        if ( params.ecc > 1 - 1/(params.a/(1+params.rp)) ) :
            params2.rhostar = rhostar_buf
            params2.Rstar = Rstar_buf
            params2.Rplan = Rplan_buf
            params.a = a_buf
            slider_rhostar.set_val(params2.rhostar)
        else :
            params2.aAU = params.a*params2.Rstar*Rsun/AU
            params.inc = calc_inc__b_aR_ecc_w(params,params2)
            params2.Rplan = calc_Rplan__rp_Rstar(params,params2)
            params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
            params2.depth = calc_depth(params,params2)
            params2.dur   = calc_dur(params,params2)
            slider_dur.set_val(params2.dur)
            slider_depth.set_val(params2.depth)
            slider_Rstar.set_val(params2.Rstar)
            slider_aR.set_val(params.a)
            slider_aAU.set_val(params2.aAU)
            slider_inc.set_val(params.inc)
            slider_Rplan.set_val(params2.Rplan)
            slider_rhoplan.set_val(params2.rhoplan)
            # update_sliders(params,params2)
            update_plot1(params,params2)
            update_plot2(params,params2)
            update_plot_rv_1(params,params2)
            update_plot_rv_2(params,params2)
            update_transit(params,params2)
            update_orbit(params,params2)
            plt.gcf().canvas.draw_idle()

# def slider_rv_update_rv_Ms(val):
#     global params,params2
#     if params2.Mstar != val :
#         params2.Mstar = val
#         params2.Mplan = calc_Mplan__K(params,params2)
#         slider_rv_Mp.set_val(params2.Mplan)
#         update_plot_rv_1(params,params2)
#         update_plot_rv_2(params,params2)
#         plt.gcf().canvas.draw_idle()

def slider_update_Mstar(val) : # 1)  Mstar -> rhostar   fix: Rstar  |  rhostar -> a  fix: P  |  a -> inc  fix: b ecc w  |  Mstar -> K     fix: Mplan
                               # 2)  Mstar -> rhostar   fix: Rstar  |  rhostar -> a  fix: P  |  a -> inc  fix: b ecc w  |  Mstar -> Mplan fix: K      |  Mplan -> rhoplan fix: Rplan
                               # 3)  Mstar -> Rstar   fix: rhostar  |  Mstar -> Mplan fix: K |  Rstar -> Rplan  fix: rp |  Rplan -> rhoplan  fix: Mplan
    global params,params2
    # if params2.Mstar != val :
    #     params2.Mstar = val
    #     params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
    #     # params2.rhostar  = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
    #     params.a = calc_aR__rhostar_per(params,params2)
    #     params2.aAU = params.a*params2.Rstar*Rsun/AU
    #     params.inc = calc_inc__b_aR_ecc_w(params,params2)
    #     params2.K = calc_K__Mplan(params,params2)
    #     # params2.Mplan = calc_Mplan__K(params,params2)
    #     # params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
    #     slider_rhostar.set_val(params2.rhostar)
    #     slider_aR.set_val(params.a)
    #     slider_aAU.set_val(params2.aAU)
    #     slider_inc.set_val(params.inc)
    #     slider_rv_K.set_val(params2.K)
    #     # slider_Mplan.set_val(params2.Mplan)
    #     # slider_rhoplan.set_val(params2.rhoplan)
    #     # update_sliders(params,params2)
    #     update_plot1(params,params2)
    #     update_plot2(params,params2)
    #     update_plot_rv_1(params,params2)
    #     update_plot_rv_2(params,params2)
    #     circ_planet.set_radius(params.rp)
    #     circ_planet.set_center((0.,params2.b))
    #     plt.gcf().canvas.draw_idle()
    if params2.Mstar != val :
        params2.Mstar = val
        params2.Rstar = calc_Rstar__rhostar_Mstar_Mplan(params2)
        params2.aAU = params.a*params2.Rstar*Rsun/AU
        params2.Mplan = calc_Mplan__K(params,params2)
        params2.Rplan = calc_Rplan__rp_Rstar(params,params2)
        params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
        slider_Rstar.set_val(params2.Rstar)
        slider_Mplan.set_val(params2.Mplan)
        slider_Rplan.set_val(params2.Rplan)
        slider_rhoplan.set_val(params2.rhoplan)
        slider_aAU.set_val(params2.aAU)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_rhoplan(val) :
    global params,params2
    if params2.rhoplan != val :
        params2.rhoplan = val
        # params2.Rplan = calc_Rplan__rhoplan_Mplan(params2)
        # params.rp       = calc_rp__Rplan_Rstar(params,params2)
        params2.Mplan = calc_Mplan__rhoplan_Rplan(params2)
        params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
        # params2.rhostar = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan
        params.a = calc_aR__rhostar_per(params,params2)
        params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.K = calc_K__Mplan(params,params2)
        slider_Mplan.set_val(params2.Mplan)
        slider_rhostar.set_val(params2.rhostar)
        slider_aR.set_val(params.a)
        slider_aAU.set_val(params2.aAU)
        slider_inc.set_val(params.inc)
        slider_rv_K.set_val(params2.K)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_Mplan(val) :
    global params,params2
    if params2.Mplan != val :
        params2.Mplan = val
        params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
        params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
        #if (params2.rhostar>10.) :
        #    params2.rhostar = 10.
        params.a = calc_aR__rhostar_per(params,params2)
        params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        params2.K = calc_K__Mplan(params,params2)
        slider_rhoplan.set_val(params2.rhoplan)
        slider_rhostar.set_val(params2.rhostar)
        slider_aR.set_val(params.a)
        slider_aAU.set_val(params2.aAU)
        slider_rv_K.set_val(params2.K)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_inc(val) :
    global params,params2
    if params.inc != val :
        params.inc = val
        params2.b = calc_b__aR_inc_ecc_w(params)
        params2.K = calc_K__Mplan(params,params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_b.set_val(params2.b)
        slider_rv_K.set_val(params2.K)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        #update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_ecc(val) :
    global params,params2
    if params.ecc != val :
        params.ecc = val
        # 1-ecc^2 >! 1/a*(2-1/a) -> planet is still outside the star rt > Rstar (see THEORY )
        # a*(1-ecc^2) >! (2-1/a)
        # 2-a*(1-ecc^2) <! 1/a
        # (1-ecc^2) * a^2 -2 a +1 >! 0
        # a = (2 +- sqrt(4-4*(1-ecc^2))) / (2(1-ecc^2))
        # a > (1 + ecc ) / (1-ecc^2) = 1/(1 - ecc)
        # a < (1 - ecc ) / (1-ecc^2) # not possible
        # a = params.a/(1+params.rp) -> Radius um Planetenradius vergroessert!
        if ( params.ecc > 1 - 1/(params.a/(1+params.rp)) ) :
            params.ecc = 1 - 1/(params.a/(1+params.rp))
            slider_ecc.set_val(params.ecc)
        #params2.b = calc_b__aR_inc_ecc_w(params)
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        params2.K = calc_K__Mplan(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_inc.set_val(params.inc)
        slider_rv_K.set_val(params2.K)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_w(val) :
    global params,params2
    if params.w != val :
        params.w = val
        #params2.b = calc_b__aR_inc_ecc_w(params)
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_inc.set_val(params.inc)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_b(val) :
    global params,params2
    if params2.b != val :
        params2.b = val
        params.inc = calc_inc__b_aR_ecc_w(params,params2)
        params2.K = calc_K__Mplan(params,params2)
        params2.depth = calc_depth(params,params2)
        params2.dur   = calc_dur(params,params2)
        slider_dur.set_val(params2.dur)
        slider_depth.set_val(params2.depth)
        slider_inc.set_val(params.inc)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        update_transit(params,params2)
        # update_orbit(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_u1(val) :
    global params,params2
    if params.u[0] != val :
        params.u[0] = val
        if (params.u[0] + params.u[1] > 1) :
            params.u[0] = 1 - params.u[1]
            # pos. 2:  params.u[1] = 1 - params.u[0]
        if (params.u[0] + 2.*params.u[1] < 0.) :
            params.u[0] =0 - 2.*params.u[1]
            # pos. 2: params.u[1] = -params.u[0]/2.
        # pos. 2: slider_u2.set_val(params.u[1])
        params2.depth = calc_depth(params,params2)
        slider_u1.set_val(params.u[0])
        slider_depth.set_val(params2.depth)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_u2(val) :
    global params,params2
    if params.u[1] != val :
        params.u[1] = val
        if (params.u[0] + params.u[1] > 1) :
            # pos. 2: params.u[0] = 1 - params.u[1]
            params.u[1] = 1 - params.u[0]
        if (params.u[0] + 2.*params.u[1] < 0.) :
            # pos. 2: params.u[0] = -2.*params.u[1]
            params.u[1]=-params.u[0]/2.
        # pos. 2: slider_u1.set_val(params.u[0])
        params2.depth = calc_depth(params,params2)
        slider_u2.set_val(params.u[1])
        slider_depth.set_val(params2.depth)
        # update_sliders(params,params2)
        update_plot1(params,params2)
        update_plot2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_F0(val) :
    global params,params2,ax
    if params2.F0 != val :
        params2.F0 = val
        params2.depth = calc_depth(params,params2)
        slider_depth.set_val(params2.depth)
        update_plot1(params,params2)
        update_plot2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_bg(val) :
    global params,params2
    if params2.bg != val :
        params2.bg = val
        # params2.depth = calc_depth(params,params2)
        # slider_depth.set_val(params2.depth)
        update_plot1(params,params2)
        update_plot2(params,params2)
        plt.gcf().canvas.draw_idle()

# def slider_update_rv_xstart(val) :
#     global xstart, t1_rv
#     xstart = val
#     t1_rv = np.linspace(xstart,xend,prec)
#     update_plot_rv_1(params2)
#     model_rv_1.set_xdata(t1_rv)
#     ax3.set_xlim([xstart,xend])
#
# def slider_update_rv_xend(val) :
#     global xend, t1_rv
#     xend = val
#     t1_rv = np.linspace(xstart,xend,prec)
#     update_plot_rv_1(params2)
#     model_rv_1.set_xdata(t1_rv)
#     ax3.set_xlim([xstart,xend])
#
# def slider_update_rv_y(val) :
#     global ystart_rv,yend_rv
#     if val!=(ystart_rv,yend_rv) :
#         ystart_rv,yend_rv = val
#         ax3.set_ylim([ystart_rv,yend_rv])
#         ax4.set_ylim([ystart_rv,yend_rv])
#
# def slider_update_rv_ystart(val) :
#     global ystart
#     ystart = val
#     ax3.set_ylim([ystart,yend])
#
# def slider_update_rv_yend(val) :
#     global yend
#     yend = val
#     ax3.set_ylim([ystart,yend])
#
# def slider_update_rv_t0(val) :
#     global params2
#     if params.t0 != val :
#         params2.t0 = val
#         slider_rv_t0.set_val(params2.t0)
#         update_plot_rv_1(params,params2)
#         update_plot_rv_2(params,params2)
#         line_rv_1.set_xdata([params.t0,params.t0])
#         plt.gcf().canvas.draw_idle()
#
# def slider_update_rv_P(val) :
#     global params2
#     if params.per != val :
#         params.per = val
#         params2.Mplan = calc_Mplan__K(params2)
#         slider_rv_P.set_val(params.per)
#         slider_rv_Mp.set_val(params2.Mplan)
#         update_plot_rv_1(params2)
#         update_plot_rv_2(params2)
#         plt.gcf().canvas.draw_idle()
#
# def slider_update_rv_ecc(val) :
#     global params2
#     if params2.ecc != val :
#         params2.ecc = val
#         params2.Mplan = calc_Mplan__K(params,params2)
#         slider_rv_Mp.set_val(params2.Mplan)
#         update_plot_rv_1(params,params2)
#         update_plot_rv_2(params,params2)
#         plt.gcf().canvas.draw_idle()
#
# def slider_update_rv_w(val) :
#     global params2
#     if params2.w != val :
#         params2.w = val
#         update_plot_rv_1(params,params2)
#         update_plot_rv_2(params,params2)
#
# def slider_update_rv_inc(val) :
#     global params2
#     if params2.inc != val :
#         params2.inc = val
#         params2.Mplan = calc_Mplan__K(params,params2)
#         slider_rv_Mp.set_val(params2.Mplan)
#         update_plot_rv_1(params,params2)
#         update_plot_rv_2(params,params2)
#         plt.gcf().canvas.draw_idle()

def slider_update_rv_K(val) :
    global params2
    if params2.K != val :
        params2.K = val
        params2.Mplan = calc_Mplan__K(params,params2)
        params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
        slider_Mplan.set_val(params2.Mplan)
        slider_rhoplan.set_val(params2.rhoplan)
        # update_sliders(params,params2)
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)
        plt.gcf().canvas.draw_idle()

def slider_update_rv_offset(val) :
    global params2
    if params2.offset != val :
        params2.offset = val
        update_plot_rv_1(params,params2)
        update_plot_rv_2(params,params2)

#  def slider_update_rv_Mp(val) :
#      global params2
#      if params2.Mplan != val :
#          params2.Mplan = val
#          params2.K = calc_K__Mplan(params,params2)
#          slider_rv_K.set_val(params2.K)
#          update_plot_rv_1(params,params2)
#          update_plot_rv_2(params,params2)
#          plt.gcf().canvas.draw_idle()
#
#  def slider_update_rv_Ms(val) :
#      global params2
#      if params2.Mstar != val :
#          params2.Mstar = val
#          params2.Mplan = calc_Mplan__K(params,params2)
#          slider_rv_Mp.set_val(params2.Mplan)
#          update_plot_rv_1(params,params2)
#          update_plot_rv_2(params,params2)
#          plt.gcf().canvas.draw_idle()


############################################################### MAIN ###############################################################

marker_style = dict(linestyle='',  marker='o', fillstyle='full', markeredgewidth=0)
line_style   = dict(linestyle='-', marker='')
onlyplot = False
#prec = 3001
prec = 501

parser = argparse.ArgumentParser(description=version)
parser.add_argument('files',    nargs='*',          help='files')
parser.add_argument('-lc',      dest='lc_filename', type=str,   default='',     help='[%(default)s] lc_filename')
parser.add_argument('-rv',      dest='rv_filename', type=str,   default='',     help='[%(default)s] rv_filename')
parser.add_argument('-save',    dest='save',        type=str,   default='',     help='[%(default)s] save with filename (pdf/png)')
#parser.add_argument('-save',   dest='save',        action='store_false',       help='[%(default)s] save')
parser.add_argument('-tex',     dest='tex',         action='store_true',        help='[%(default)s] tex Titel')
parser.add_argument('-tit',     dest='tit',         type=str,   default='',     help='[%(default)s] Titel')
parser.add_argument('-Pmax',    dest='Pmax',        type=float, default='0',    help='[%(default)s] Pmax')
parser.add_argument('-amax',    dest='amax',        type=float, default='0',    help='[%(default)s] amax')
parser.add_argument('-t0',      dest='t0',          type=float, default='999',  help='[%(default)s] t0')
parser.add_argument('-P',       dest='P',           type=float, default='2.',   help='[%(default)s] Periode P')
parser.add_argument('-Mstar',   dest='Mstar',       type=float, default='0',    help='[%(default)s] Masse Stern  / Mstar')
parser.add_argument('-Mplan',   dest='Mplan',       type=float, default='0',    help='[%(default)s] Masse Planet / Mplan')
parser.add_argument('-rp',      dest='rp',          type=float, default='0.1',  help='[%(default)s] Radiusverhaeltnis rp')
parser.add_argument('-Rstar',   dest='Rstar',       type=float, default='0',    help='[%(default)s] Radius Stern Rstar')
parser.add_argument('-Rplan',   dest='Rplan',       type=float, default='0',    help='[%(default)s] Radius Planet Rplan')
parser.add_argument('-REB',     dest='REB',         type=float, default='0',    help='[%(default)s] Radius EB 2nd star')
parser.add_argument('-a',       dest='a',           type=float, default='0.',   help='[%(default)s] Grosse Halbachse a')
parser.add_argument('-rhostar', dest='rhostar',     type=float, default='1.0',  help='[%(default)s] mittl. Sterndichte rhostar')
parser.add_argument('-rhoplan', dest='rhoplan',     type=float, default='0.0',  help='[%(default)s] mittl. Planetendichte rhoplan')
parser.add_argument('-b',       dest='b',           type=float, default='0.',   help='[%(default)s] Impaktparameter b')
parser.add_argument('-i',       dest='i',           type=float, default='0.',   help='[%(default)s] Inklination i')
parser.add_argument('-u1',      dest='u1',          type=float, default='0.',   help='[%(default)s] Limb-Darkening u1')
parser.add_argument('-u2',      dest='u2',          type=float, default='0.',   help='[%(default)s] Limb-Darkening u2')
parser.add_argument('-e',       dest='e',           type=float, default='0.',   help='[%(default)s] Excentrizitaet e')
parser.add_argument('-w',       dest='w',           type=float, default='90.',  help='[%(default)s] Argument Periastron w')
parser.add_argument('-norm',    dest='norm',        type=float, default='0.',   help='[%(default)s] norm')
parser.add_argument('-x0',      dest='x0',          type=float, default=0,      help='[%(default)s] x0')
parser.add_argument('-x1',      dest='x1',          type=float, default=0,      help='[%(default)s] x1')
parser.add_argument('-y0',      dest='y0',          type=float, default=0,      help='[%(default)s] y0')
parser.add_argument('-y1',      dest='y1',          type=float, default=0,      help='[%(default)s] y1')
parser.add_argument('-zoom',    dest='zoom',        type=float, default=1.,     help='[%(default)s] zoom')
parser.add_argument('-alpha',   dest='alpha',       type=float, default=0.05,   help='[%(default)s] alpha')
parser.add_argument('-F0',      dest='F0',          type=float, default=1.,     help='[%(default)s] F0')
parser.add_argument('-bg',      dest='bg',          type=float, default=0.,     help='[%(default)s] bg')
parser.add_argument('-K',       dest='K',           type=float, default=100.,   help='[%(default)s] K')
parser.add_argument('-off',     dest='offset',      type=float, default=0.,     help='[%(default)s] offset')
args = parser.parse_args()


# reading data

#print(' args.files = ', args.files)
#if len(args.files)==0 and args.lc_filename=='' :
#    print("no light curve specified")
#    exit(-1)

lc_filename=''
if len(args.files)!=0 :
    lc_filename=args.files[0]
if args.lc_filename!='' :
    lc_filename=args.lc_filename

rv_filename = args.rv_filename

if lc_filename=='' and rv_filename=='' :
    print("no data specified")
    exit(-1)


##########################  LC data

print('loading LC data:\n',lc_filename)

norm = 1.
myalpha    = args.alpha
myalpha_rv = args.alpha*3

if lc_filename!='' :
    #if lc_filename.endswith('.fits') and not lc_filename.endswith('tp.fits') and not lc_filename.endswith('dvt.fits') :
    if lc_filename.endswith('.fits') :
        tdata, Fdata, edata  = get_TESS_data(lc_filename)
    elif lc_filename.endswith('.csv') :
        # awk '(NR>1){print $5,$15,$25}' TIC233497719-01_202*.tbl > TOI1877_43cm.tab
        tdata, Fdata, edata  = get_TESS_data_tab(lc_filename)
    elif lc_filename.endswith('.tab') or lc_filename.endswith('.tbl') :
        # awk '(NR>1){print $5,$15,$25}' TIC233497719-01_202*.tbl > TOI1877_43cm.tab
        tdata, Fdata, edata  = get_WST_data(lc_filename)
    else :
        print('  no flux data')
        #continue
        exit(-1)

    #print(tdata.ndim)
    #print(Fdata.ndim)
    #print(edata.ndim)
    #print(tdata[0:10])
    #print(Fdata[0:10])
    #print(edata[0:10])

    if (args.norm==0.) :
        norm = np.median(Fdata)
        #norm = np.percentile(Fdata,95.)
    else :
        norm = args.norm
    Fdata /= norm
    edata /= norm

    [tmin,tmax] = [np.min(tdata),np.max(tdata)]
    [ymin,ymax] = [np.min(Fdata),np.max(Fdata)]



# if args.REB != 0. :
#     if args.Rplan == 0. :
#         args.Rplan = args.REB*Rsun/Rearth
#     else :
#         print("R_EB and R_planet ???")
#         exit(-1)

###################### RV data

if rv_filename!='' :
    print('loading RV data:\n',rv_filename)
    x_rv, y_rv, err_rv  = np.genfromtxt(rv_filename, comments='#',dtype="f8,f8,f8", unpack=True)
    x_rv -= 2450000.
    if lc_filename=='' :
        t0start = (np.min(x_rv)+np.max(x_rv))/2.
        [tmin,tmax]  = [t0start-2.5*args.P,t0start+2.5*args.P]
        [ymin,ymax] = [0.9,1.1]

if args.t0==999 :
    t0start = (tmin+tmax)/2.
else :
    t0start = args.t0



# # simulated data
# tlen = 83
# x_rv = np.sort(np.random.uniform(low=tstart,high=tend, size=tlen))
# yerr_rv = 5.0 #assuming we know the Gaussian uncertainty
# synth_params2 = radvel.Parameters(1,basis='per tc e w k')
# synth_params2['per1'] = radvel.Parameter(value = params2.per)
# synth_params2['tc1']  = radvel.Parameter(value = params.t0)
# synth_params2['e1']   = radvel.Parameter(value = params2.ecc)
# synth_params2['w1']   = radvel.Parameter(value = params2.w)
# synth_params2['k1']   = radvel.Parameter(value = params2.K)
# synth_params2['dvdt'] = radvel.Parameter(value = 0)
# synth_params2['curv'] = radvel.Parameter(value = 0)
# synth_model = radvel.RVModel(params2=synth_params2)
# y_rv = synth_model(x_rv) + yerr_rv * np.random.randn(len(x_rv)) + params2.offset
#
# params2 = RVParams()
# params.t0     = 2456301.6 - 2400000.  # time of inferior conjunction
# params.per    = 200.31                # orbital period
# params2.ecc    = 0.                    # eccentricity
# params2.w      = 0.                    # longitude of periastron (in radians)
# params2.K      = 39.1
# params2.offset = 0.
#
# tstart = 2456200. -  2400000.
# tend   = 2457140. -  2400000.

#rc_filename = '51Peg_RV_WST_2-1m.tab'     # 51 Peg Wendelstein

if lc_filename!='' :
    [tmin_rv,tmax_rv] = [tmin,tmax]
    [ymin_rv,ymax_rv] = [args.offset-2*args.K,args.offset+2*args.K]
if rv_filename!='' :
    [tmin_rv,tmax_rv]  = [np.min(x_rv),np.max(x_rv)]
    [ymin_rv,ymax_rv]  = [np.min(y_rv),np.max(y_rv)]


#print("tmin_rv,tmax_rv = ",tmin_rv,tmax_rv)
#print("ymin_rv,ymax_rv = ",ymin_rv,ymax_rv)


################################ model parameters


class OtherParams :
    def __init__(self) :
        self.b        = 0.
        self.rhostar  = 0.
        self.rhoplan  = 0.
        self.Mstar    = 0.
        self.Mplan    = 0.
        self.Rstar    = 0.
        self.Rplan    = 0.
        self.dur      = 0.
        self.depth    = 0.
        self.bg       = 0.
        self.F0       = 0.
        #self.t0     = 0.   # RV
        #self.P    = 0.     # RV
        #self.ecc    = 0.   # RV
        #self.w      = 0.   # RV
        self.K      = 0.   # RV
        self.offset = 0.   # RV
        self.q      = [0.,0.]
        #self.inc     = 90. # RV
        #self.Mp     = 0.   # RV
        #self.Ms     = 0.   # RV

params  = batman.TransitParams()
params2 = OtherParams()


def print_all(params,params2) :
    q1 = (params.u[0]+params.u[1])**2                    # Kipping 2013 Eq. 17
    if (params.u[0]+params.u[1])==0. :
        q2 = 0.25
    else :
        q2 = params.u[0]/(2.*(params.u[0]+params.u[1]))  # Kipping 2013 Eq. 18
    print('  t0 [d]           =',params.t0   )
    print('  period [d]       =',params.per  )
    print('  a/Rstar          =',params.a    )
    print('  a [AU]           =',params2.aAU  )
    print('  rp = Rplan/Rstar =',params.rp   )
    print('  ecc              =',params.ecc  )
    print('  w                =',params.w    )
    print('  u1               =',params.u[0] )
    print('  u2               =',params.u[1] )
    print('  q1 (juliet)      =', q1 )
    print('  q2 (juliet)      =', q2 )
    print('  inc [deg]        =',params.inc )
    print('  b                =',params2.b       )
    print('  rhostar [g/cm3]  =',params2.rhostar )
    print('  rhoplan [g/cm3]  =',params2.rhoplan )
    print('  Mstar [Msun]     =',params2.Mstar   )
    print('  Mplan [Mjup]     =',params2.Mplan*Mearth/Mjup   )
    print('  Mplan [Mearth]   =',params2.Mplan   )
    print('  Rstar [Rsun]     =',params2.Rstar   )
    print('  Rplan [Rjup]     =',params2.Rplan*Rearth/Rjup   )
    print('  Rplan [Rearth]   =',params2.Rplan   )
    print('  bg               =',params2.bg      )
    print('  F0               =',params2.F0      )
    print('  K [m/s]          =',params2.K       )
    print('  offset [m/s]     =',params2.offset  )
    print('  dur [h]          =',params2.dur     )
    print('  dur [min]        =',params2.dur*60. )
    print('  depth            =',params2.depth   )

params.t0  = t0start                  # time of inferior conjunction
params.per = args.P                   # orbital period
params.rp  = args.rp                  # planet radius (in units of stellar radii)
params.ecc = args.e                   # eccentricity
params.w   = args.w                   # longitude of periastron (in degrees)
params.u   = [args.u1,args.u2]        # limb darkening coefficients [u1, u2]
if (params.u[0] + params.u[1] > 1) :
    params.u[1] = 1 - params.u[0]
if (params.u[0] + 2.*params.u[1] < 0.) :
    params.u[1] = -params.u[0]/2.
params.limb_dark = "quadratic"        # limb darkening model
params2.b  = args.b
params2.bg = args.bg
params2.F0 = args.F0


if args.a > 0. :
    params.a = args.a
    if args.rhostar == 0. and params.per!=0 :
        params2.rhostar = calc_rhostar__aR_per(params)
        if args.Mstar == 0. and args.Rstar != 0. :
            params2.Rstar   = args.Rstar
            params2.Mstar   = calc_Mstar__rhostar_Rstar_Mplan(params2)
        elif args.Mstar != 0. and args.Rstar == 0. :
            params2.Mstar   = args.Mstar
            params2.Rstar   = calc_Rstar__rhostar_Mstar_Mplan(params2)
        elif  args.Mstar == 0. and args.Rstar == 0. :
            params2.Rstar   = 1.
            params2.Mstar   = calc_Mstar__rhostar_Rstar_Mplan(params2)
        else :
            print("too many given parameters:")
            print("args.a       = ",args.a)
            print("args.Mstar   = ",args.Mstar)
            print("args.Rstar   = ",args.Rstar)
    else :
        print("inconsistent given parameters:")
        print("args.a        = ",args.a)
        print("args.rhostar = ",args.rhostar)
        print("args.P        = ",args.P)
else :
    if args.rhostar == 0. :
        if args.Mstar != 0. and args.Rstar != 0. :
            params2.Mstar   = args.Mstar
            params2.Rstar   = args.Rstar
            params2.rhostar = calc_rhostar__Mstar_Rstar(params2)
            params.a        = calc_aR__rhostar_per(params,params2)   # semi-major axis (in units of stellar radii)
            params2.aAU     =  params.a * (params2.Rstar*Rsun) / AU
        else :
            print("inconsistent given parameters:")
            print("args.a       = ",args.a)
            print("args.rhostar = ",args.rhostar)
            print("args.P       = ",args.P)
            exit(-1)
    else :
        params2.rhostar = args.rhostar
        if   args.Mstar == 0. and args.Rstar != 0. :
            params2.Rstar   = args.Rstar
            params2.Mstar   = calc_Mstar__rhostar_Rstar_Mplan(params2)
            params.a        = calc_aR__rhostar_per(params,params2)   # semi-major axis (in units of stellar radii)
            params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        elif args.Mstar != 0. and args.Rstar == 0. :
            params2.Mstar   = args.Mstar
            params2.Rstar   = calc_Rstar__rhostar_Mstar_Mplan(params2)
            params.a        = calc_aR__rhostar_per(params,params2)   # semi-major axis (in units of stellar radii)
            params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        elif  args.Mstar == 0. and args.Rstar == 0. :
            params2.Rstar   = 1.
            params2.Mstar   = calc_Mstar__rhostar_Rstar_Mplan(params2)
            params.a        = calc_aR__rhostar_per(params,params2)   # semi-major axis (in units of stellar radii)
            params2.aAU =  params.a * (params2.Rstar*Rsun) / AU
        else :
            print("too many given parameters:")
            print("args.rhostar = ",args.rhostar)
            print("args.Mstar   = ",args.Mstar)
            print("args.Rstar   = ",args.Rstar)
            exit(-1)

if args.i != 0. :
    params.inc = args.i
    params2.b = calc_b__aR_inc_ecc_w(params)
params.inc = calc_inc__b_aR_ecc_w(params,params2) # orbital inclination (in degrees)

if args.rp == 0. :
    if args.rhoplan == 0. :
        if args.K>0 :
            params2.K      = args.K
            if args.Mplan == 0. and args.Rplan != 0. :
                params2.Rplan   = args.Rplan
                params2.Mplan   = calc_Mplan__K(params,params2)
                params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
                params.rp = calc_rp__Rplan_Rstar(params,params2)
            else :
                print("inconsistent given parameters:")
                print("args.K       = ",args.K)
                print("args.Mplan   = ",args.Mplan)
                print("args.Rplan   = ",args.Rplan)
                print("args.a       = ",args.a)
                print("args.rhoplan = ",args.rhoplan)
                print("args.P       = ",args.P)
                exit(-1)
        else :
            if args.Mplan != 0. and args.Rplan != 0. :
                params2.Mplan   = args.Mplan
                params2.Rplan   = args.Rplan
                params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
                params.rp = calc_rp__Rplan_Rstar(params,params2)
                params2.K   = calc_K__Mplan(params,params2)
            else :
                print("inconsistent given parameters:")
                print("args.K       = ",args.K)
                print("args.Mplan   = ",args.Mplan)
                print("args.Rplan   = ",args.Rplan)
                print("args.a       = ",args.a)
                print("args.rhoplan = ",args.rhoplan)
                print("args.P       = ",args.P)
                exit(-1)
    else :
        params2.rhoplan = args.rhoplan
        if args.K==0 :
            if args.Mplan == 0. and args.Rplan != 0. :
                params2.Rplan   = args.Rplan
                params2.Mplan   = calc_Mplan__rhoplan_Rplan(params2)
            elif args.Mplan != 0. and args.Rplan == 0. :
                params2.Mplan   = args.Mplan
                params2.Rplan   = calc_Rplan__rhoplan_Mplan(params2)
            elif  args.Mplan == 0. and args.Rplan == 0. :
                params2.Rplan   = 1.
                params2.Mplan   = calc_Mplan__rhoplan_Rplan(params2)
            else :
                print("too many given parameters:")
                print("args.rhoplan = ",args.rhoplan)
                print("args.Mplan   = ",args.Mplan)
                print("args.Rplan   = ",args.Rplan)
                exit(-1)
        else :
            params2.K      = args.K
            if  args.Mplan == 0. and args.Rplan == 0. :
                params2.Mplan   = calc_Mplan__K(params,params2)
                params2.Rplan   = calc_Rplan__rhoplan_Mplan(params2)
            else :
                print("too many given parameters:")
                print("args.rhoplan = ",args.rhoplan)
                print("args.Mplan   = ",args.Mplan)
                print("args.Rplan   = ",args.Rplan)
                exit(-1)
else :
    params.rp = args.rp
    # params2.Rplan = calc_Rplan__rp_Rstar(params,params2)
    if args.rhoplan == 0. :
        if args.K > 0 :
            params2.K = args.K
            if args.Mplan == 0. and args.Rplan == 0. :
                params2.Mplan   = calc_Mplan__K(params,params2)
                params2.Rplan   = calc_Rplan__rp_Rstar(params,params2)
                params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
            else :
                print("too many given parameters:")
                print("args.rp      = ",args.rp)
                print("args.Rplan   = ",args.Rplan)
                print("args.K       = ",args.K)
                print("args.Mplan   = ",args.Mplan)
                print("args.Mstar   = ",args.Mstar)
                print("args.rhoplan = ",args.rhoplan)
                exit(-1)
        else :
            if   args.rhoplan == 0. and args.Mplan != 0. :
                params2.Mplan   = args.Mplan
                params2.rhoplan = calc_rhoplan__Mplan_Rplan(params2)
            elif args.rhoplan != 0. and args.Mplan == 0. :
                params2.rhoplan = args.rhoplan
                params2.Mplan   = calc_Mplan__rhoplan_Rplan(params2)
            elif args.rhoplan == 0. and args.Mplan == 0. :
                params2.rhoplan = args.rhoplan
                params2.Mplan   = args.Mplan
            else :
                print("inconsistent given parameters:")
                print("args.a       = ",args.a)
                print("args.rhoplan = ",args.rhoplan)
                print("args.P       = ",args.P)
                exit(-1)

# params.rp = calc_rp__Rplan_Rstar(params,params2)
# params2.rhostar  = calc_rhostar__rhostar_Mplan_Rstar(params2) # Korrektur mit Mplan

params2.aAU   = params.a*params2.Rstar*Rsun/AU
params2.dur   = calc_dur(params,params2)
params2.depth = params.rp**2 * 1000. # first guess
#params2.t0     = t0start   # time of inferior conjunction
#params2.P      = args.P    # orbital period
#params2.ecc    = args.e    # eccentricity
#params2.w      = args.w
params2.offset = args.offset
#params2.inc    = params.inc
#params2.Ms     = 1.
#params2.MP     = calc_Mplan__K(params,params2)


############ plotting areas and slider ranges

Pmax = 20.
amax = 20.
if args.Pmax==0 and args.P!=2.:
    Pmax = 3*args.P
    amax = 3*params.a

if args.x0!=0. and args.x1!=0. :
    [xstart,xend] = [args.x0,args.x1]
else :
    if lc_filename!='' :
        [xstart,xend] = calc_xlim(tdata)
    elif rv_filename!='' :
        # [xstart,xend] = calc_xlim(x_rv)
        [xstart,xend] = [tmin,tmax]

if args.y0!=0. and args.y1!=0. :
    [ystart,yend] = [args.y0,args.y1]
else :
    if lc_filename!='' :
        [ystart,yend] = calc_ylim(Fdata)
    else :
        # [ystart,yend] = calc_ylim(f1)
        [ystart,yend] = [0.9,1.1]

[xstart_rv,xend_rv] = [xstart,xend]
[ystart_rv,yend_rv] = [ymin_rv,ymax_rv]
if rv_filename!='' :
    [xstart_rv,xend_rv] = calc_rv_xlim(x_rv)
    [ystart_rv,yend_rv] = calc_rv_ylim(y_rv)

zoom    = args.zoom
zoom_rv = 1.



#########################################

if onlyplot :
    print('noch falsch...')
    fig = plt.figure(figsize=(8,4),facecolor='w',dpi=100)
    ax1 = plt.axes([0.10, 0.12, 0.88, 0.78], facecolor='w')
    ax1.plot( tdata, Fdata, c='b' ,zorder=0, alpha=myalpha, ms=7, **marker_style)
    ax1.set_xlim([xstart,xend])
    ax1.set_ylim([ystart,yend])
    ax1.set_title(lc_filename)
    ax1.set_xlabel("Julian Date - 2450000 [d]")
    ax1.set_ylabel("normalized flux")
    #ax1.xaxis.set_minor_locator(AutoMinorLocator(n=4))
    #ax1.yaxis.set_minor_locator(AutoMinorLocator(n=4))
    plt.savefig(lc_filename+'.pdf')
    plt.show()
    exit(-1)


# fig = plt.figure(figsize=(11,7),facecolor='w',dpi=100)
# fig = plt.figure(figsize=(13,7),facecolor='w',dpi=100)
fig = plt.figure(figsize=(14,8),facecolor='w',dpi=100)

xx1 = 0.07
#xx2 = 0.475
xx2 = 0.48
lx  = 0.37
yy1 = 0.45
yy2 = 0.75
ly  = 0.2

ax0 = plt.axes([0.00, 0.95,  1.00,  0.05], facecolor='w')        # title
ax1 = plt.axes([xx1,         yy2, lx, ly],        facecolor='w') # lc 1
ax2 = plt.axes([xx2,         yy2, lx, ly],        facecolor='w') # lc 2
ax3 = plt.axes([xx1,         yy1, lx, ly],        facecolor='w') # rv 1
ax4 = plt.axes([xx2,         yy1, lx, ly],        facecolor='w') # rv 2
# ax5 = plt.axes([xx2+lx-0.03, yy2, ly, ly], facecolor='w', xlabel=None, ylabel=None) # star-planet
ax5 = plt.axes([xx2+lx-0.015, yy2, ly, ly], facecolor='w', ylabel=None) # star-planet
# ax6 = plt.axes([xx2+lx-0.03, yy1, ly, ly], facecolor='w', xlabel=None, ylabel=None, polar=True) # orbit
# ax6 = plt.axes([xx2+lx-0.03, yy1, ly, ly], facecolor='w', xlabel=None, ylabel=None) # orbit
ax6 = plt.axes([xx2+lx-0.015, yy1, ly, ly], facecolor='w', ylabel=None) # orbit

ax0.axis('off')
#print(">>>>>>>>>>>>>>.",lc_filename)
ax0.text(0.5,0.5,args.tit+'  '+lc_filename+'  '+rv_filename,size = 12, va='center_baseline', ha='center', color='k')

#ax1.xaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
#ax1.yaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
ax1.tick_params(axis='both', bottom=True, top=False, left=True, right=False)
#ax1.tick_params(which='major', labelsize=12, length=10, axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax1.tick_params(which='minor',               length=5 , axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
ax1.set_xlim([xstart,xend])
ax1.set_ylim([ystart,yend])
ax1.set_xlabel("julian date - 2450000 [d]",size=8)
#ax1.set(xlabel=None)
ax1.set(ylabel=None)

#ax2.xaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
#ax2.yaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
ax2.tick_params(axis='both', bottom=True, top=False, left=True,  right=False)
#ax2.tick_params(which='major', labelsize=12, length=10, axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax2.tick_params(which='minor',               length=5 , axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
ax2.set_xlim([-params.per/2 / zoom, params.per/2 / zoom])
#ax2.set_xlim([-3, 3])
ax2.set_ylim([ystart,yend])
#ax2.set(xlabel=None)
ax2.set_xlabel("(julian date) mod P [d]",size=8)
ax2.set(ylabel=None)
#ax2.set_yticklabels([])

## ax3 = plt.axes([0.07, 0.65, 0.90, 0.21], facecolor='w')
ax3.set_xlim([xstart_rv,xend_rv])
ax3.set_ylim([ystart_rv,yend_rv])
#ax3.set(xlabel=None)
#ax3.set(ylabel=None)
#ax3.set_ylabel('Radial Velocity')
#ax3.xaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
#ax3.yaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
ax3.tick_params(axis='both', bottom=True, top=False, left=True, right=False)
#ax3.tick_params(which='major', labelsize=12, length=10, axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax3.tick_params(which='minor',               length=5 , axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax3.set_xlim([xstart,xend])
#ax3.set_ylim([ystart,yend])
ax3.set(xlabel=None)
ax3.set_xlabel("julian date - 2450000 [d]",size=8)
ax3.set(ylabel=None)

## ax4 = plt.axes([0.07, 0.36, 0.90, 0.21], facecolor='w')
#ax4.xaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
#ax4.yaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
ax4.tick_params(axis='both', bottom=True, top=False, left=True, right=False)
#ax4.tick_params(which='major', labelsize=12, length=10, axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax4.tick_params(which='minor',               length=5 , axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax4.set_xlim([-params.per/2 / zoom, params.per/2 / zoom])
ax4.set_xlim([-params.per/2 / zoom_rv, params.per/2 / zoom_rv])
ax4.set_ylim([ystart_rv,yend_rv])
#ax4.set(xlabel=None)
#ax4.set(ylabel=None)
#ax3.set_xlabel("Julian Date - 2450000 [d]")
#ax4.set_ylabel('Radial Velocity')
#plt.title('RV Model and Data')
#ax4.xaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
#ax4.yaxis.set_minor_locator(tkr.AutoMinorLocator(n=5))
ax4.tick_params(axis='both', bottom=True, top=False, left=True,  right=False)
#ax4.tick_params(which='major', labelsize=12, length=10, axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax4.tick_params(which='minor',               length=5 , axis='both', direction='in',
#                   bottom=True, top=True, left=True, right=True)
#ax4.set_xlim([-params.per/2 / zoom_rv, params.per/2 / zoom_rv])
#ax4.set_xlim([-3, 3])
#ax4.set_ylim([ystart,yend])
#ax4.set(xlabel=None)
ax4.set_xlabel("(julian date) mod P [d]",size=8)
ax4.set(ylabel=None)
#ax4.set_yticklabels([])

ax5.set_aspect(1)
#ax5.set(xlabel=None)
#ax5.set(ylabel=None)
#ax5.set_facecolor('black')
#ax5.set_facecolor('white')
#ax5.set_xticks([])
#ax5.set_xticks([], minor=True)
ax5.set_yticks([])
ax5.set_yticks([], minor=True)
ax5.tick_params(axis='both', bottom=True, top=False, left=False, right=False)
### ax5.axis("off")
ax5.axis([-1.5,1.5,-1.5,1.5])
ax5.text(0.5,1.08,'transit',size=10,transform=ax5.transAxes,va='top',ha='center',color='k')
ax5.set_ylabel("y [Rsun]",size=8)
ax5.set_xlabel("x [Rsun]",size=8)

ax6.set_aspect(1)
#ax6.grid(False)
#ax6.set_xticklabels([])
#ax6.set_yticklabels([])
#ax6.set(xlabel=None)
#ax6.set(ylabel=None)
ax6.set_yticks([])
ax6.set_yticks([], minor=True)
ax6.tick_params(axis='both', bottom=True, top=False, left=False, right=False)
### ax6.axis("off")
ax6.text(0.5,1.08,'orbit',size=10,transform=ax6.transAxes,va='top',ha='center',color='k')
ax6.set_ylabel("z [AU]",size=8)
ax6.set_xlabel("x [AU]",size=8)

plt.text(0.003,yy2+ly,     "LC", horizontalalignment='left', verticalalignment='top', fontsize=16, transform=plt.gcf().transFigure)
plt.text(0.003,yy2+ly-0.03, "[phot]", horizontalalignment='left', verticalalignment='top', fontsize=8, transform=plt.gcf().transFigure)
plt.text(0.003,yy1+ly,     "RV", horizontalalignment='left', verticalalignment='top', fontsize=16, transform=plt.gcf().transFigure)
plt.text(0.003,yy1+ly-0.03, "[m/s]", horizontalalignment='left', verticalalignment='top', fontsize=8, transform=plt.gcf().transFigure)

###################### sliders

axcolor = 'lightgrey'
start1 = 0.07
start2 = 0.56
starty = 0.04
lenx = 0.38
Dfin = 0.025

# ax_xstart  = plt.axes([start1,      0.31, lenx,      0.02], facecolor=axcolor)
# ax_xend    = plt.axes([start2,      0.31, lenx,      0.02], facecolor=axcolor)
# ax_ystart  = plt.axes([start1,      0.29, lenx,      0.02], facecolor=axcolor)
# ax_yend    = plt.axes([start2,      0.29, lenx,      0.02], facecolor=axcolor)
ax_x         = plt.axes([xx1+0.01,    yy2-0.07, lx-0.09,   0.02], facecolor=axcolor)
ax_y         = plt.axes([xx1+0.01,    yy2-0.09, lx-0.09,   0.02], facecolor=axcolor)
ax_zoom      = plt.axes([xx2+0.03,    yy2-0.07, lx-0.06,   0.02], facecolor=axcolor)
ax_alpha     = plt.axes([xx2+0.03,    yy2-0.09, lx-0.06,   0.02], facecolor=axcolor)

ax_x_rv      = plt.axes([xx1+0.01,    yy1-0.07, lx-0.09,   0.02], facecolor=axcolor)
ax_y_rv      = plt.axes([xx1+0.01,    yy1-0.09, lx-0.09,   0.02], facecolor=axcolor)
ax_zoom_rv   = plt.axes([xx2+0.03,    yy1-0.07, lx-0.06,   0.02], facecolor=axcolor)
ax_alpha_rv  = plt.axes([xx2+0.03,    yy1-0.09, lx-0.06,   0.02], facecolor=axcolor)

ax_per       = plt.axes([start1,      starty+0.27, lenx,      0.02], facecolor=axcolor)
ax_perfine   = plt.axes([start1+Dfin, starty+0.25, lenx-Dfin, 0.02], facecolor=axcolor)
ax_aR        = plt.axes([start1,      starty+0.22, lenx,      0.02], facecolor=axcolor)
ax_aAU       = plt.axes([start1,      starty+0.20, lenx,      0.02], facecolor=axcolor)
ax_Rstar     = plt.axes([start1,      starty+0.17, lenx,      0.02], facecolor=axcolor)
ax_Mstar     = plt.axes([start1,      starty+0.15, lenx,      0.02], facecolor=axcolor)
ax_rhostar   = plt.axes([start1,      starty+0.13, lenx,      0.02], facecolor=axcolor)
ax_inc       = plt.axes([start1,      starty+0.10, lenx,      0.02], facecolor=axcolor)
ax_ecc       = plt.axes([start1,      starty+0.08, lenx,      0.02], facecolor=axcolor)
ax_w         = plt.axes([start1,      starty+0.06, lenx,      0.02], facecolor=axcolor)
ax_F0        = plt.axes([start1,      starty+0.03, lenx,      0.02], facecolor=axcolor)
ax_rv_K      = plt.axes([start1,      starty+0.00, lenx,      0.02], facecolor=axcolor)
ax_dur       = plt.axes([start1,      starty-0.03, lenx,      0.02], facecolor=axcolor)

ax_t0        = plt.axes([start2,      starty+0.27, lenx,      0.02], facecolor=axcolor)
ax_t0fine    = plt.axes([start2+Dfin, starty+0.25, lenx-Dfin, 0.02], facecolor=axcolor)
ax_rp        = plt.axes([start2,      starty+0.20, lenx,      0.02], facecolor=axcolor)
ax_Rplan     = plt.axes([start2,      starty+0.17, lenx,      0.02], facecolor=axcolor)
ax_Mplan     = plt.axes([start2,      starty+0.15, lenx,      0.02], facecolor=axcolor)
ax_rhoplan   = plt.axes([start2,      starty+0.13, lenx,      0.02], facecolor=axcolor)
ax_b         = plt.axes([start2,      starty+0.10, lenx,      0.02], facecolor=axcolor)
ax_u1        = plt.axes([start2,      starty+0.08, lenx,      0.02], facecolor=axcolor)
ax_u2        = plt.axes([start2,      starty+0.06, lenx,      0.02], facecolor=axcolor)
ax_bg        = plt.axes([start2,      starty+0.03, lenx,      0.02], facecolor=axcolor)
ax_rv_offset = plt.axes([start2,      starty+0.00, lenx,      0.02], facecolor=axcolor)
ax_depth     = plt.axes([start2,      starty-0.03, lenx,      0.02], facecolor=axcolor)

#############################################################

title_x0       = r'x0'
title_y0       = r'y0'
title_x1       = r'x1'
title_y1       = r'y1'
title_x        = r'x'
title_y        = r'y'
title_zoom     = r'zoom'
title_alpha    = r'alpha'
title_x_rv     = r'x'
title_y_rv     = r'y'
title_zoom_rv  = r'zoom'
title_alpha_rv = r'alpha'
title_t0       = r't0'
title_t0fine   = r'fine:'
title_per      = r'P [d]'
title_perfine  = r'fine:'
title_aR       = r'a / Rs'
title_aAU      = r'a [AU]'
title_rhostar  = r'rhos [g/cm3]'
title_Rstar    = r'Rs [Rsun]'
title_Mstar    = r'Ms [Msun]'
title_rhoplan  = r'rhop [g/cm3]'
title_Rplan    = r'Rp [Rearth]'
title_Mplan    = r'Mp [Mearth]'
title_rp       = r'Rp / Rs'
title_b        = r'b'
title_inc      = r'i [deg]'
title_u1       = r'u1'
title_u2       = r'u2'
title_ecc      = r'e'
title_w        = r'omega'
title_F0       = r'F0'
title_bg       = r's'
title_K        = r'K [m/s]'
title_offset   = r'offset [m/s]'
title_dur      = r'dur [h]'
title_depth    = r'depth [ppt]'

if args.tex :
    plt.rcParams['text.usetex'] = True
    title_xstart  = r'$x0$'
    title_xend    = r'$x1$'
    title_ystart  = r'$y0$'
    title_yend    = r'$y1$'
    title_x       = r'$x$'
    title_y       = r'$y$'
    title_zoom    = r'zoom'
    title_alpha   = r'alpha'
    title_t0      = r'$t_0$'
    title_t0fine  = r'fine:'
    title_per     = r'$P\,[\mathrm{d}]$'
    title_perfine = r'fine:'
    title_aR      = r'$a\,/\,R_\mathrm{s}$'
    title_aAU     = r'$a$ [AU]'
    title_rhostar = r'$\rho_{\mathrm{s}}\,[\mathrm{g}/\mathrm{cm}^3]$'
    title_Rstar   = r'$R_{\mathrm{s}}\,[\,R_\odot]$'
    title_Mstar   = r'$M_{\mathrm{s}}\,[\,M_\odot]$'
    title_rhoplan = r'$\rho_{\mathrm{p}}\,[\mathrm{g}/\mathrm{cm}^3]$'
    title_Rplan   = r'$R_{\mathrm{p}}\,[\,R_\oplus]$'
    title_Mplan   = r'$M_{\mathrm{p}}\,[\,M_\oplus]$'
    title_rp      = r'$R_{\mathrm{p}}\,/\,R_\mathrm{s}$'
    title_b       = r'$b$'
    title_inc     = r'$i\,[^\circ]$'
    title_u1      = r'$u_1$'
    title_u2      = r'$u_2$'
    title_ecc     = r'$e$'
    title_w       = r'$\omega$'
    title_F0      = r'$F_0$'
    title_bg      = r'$s$'
    title_K       = r'$K$ [m/s]'
    title_offset  = r'$offset$ [m/s]'
    title_dur     = r'$\mathrm{dur}\,[\mathrm{h}]$'
    title_depth   = r'$\mathrm{depth}\,[\mathrm{ppt}]$'


# https://matplotlib.org/stable/api/widgets_api.html#matplotlib.widgets.Slider
# class matplotlib.widgets.Slider     (ax, label, valmin, valmax, *, valinit=0.5,  valfmt=None, closedmin=True, closedmax=True, slidermin=None, slidermax=None, dragging=True, valstep=None, orientation='horizontal', initcolor='r', track_color='lightgrey', handle_style=None, **kwargs)[source]
# class matplotlib.widgets.RangeSlider(ax, label, valmin, valmax, *, valinit=None, valfmt=None, closedmin=True, closedmax=True,                                 dragging=True, valstep=None, orientation='horizontal',                track_color='lightgrey', handle_style=None, **kwargs)[source]

# slider_xstart  = Slider(ax_xstart,    title_x0,      valmin=tmin,           valmax=tmax,           valfmt="%1.3f", valinit=xstart,          color='b')
# slider_xend    = Slider(ax_xend,      title_y0,      valmin=tmin,           valmax=tmax,           valfmt="%1.3f", valinit=xend,            color='b')
# slider_ystart  = Slider(ax_ystart,    title_x1,      valmin=0.,             valmax=2.,             valfmt="%1.3f", valinit=ystart,          color='b')
# slider_yend    = Slider(ax_yend,      title_y1,      valmin=0.,             valmax=2.,             valfmt="%1.3f", valinit=yend,            color='b')
slider_x    = RangeSlider(ax_x,         title_x,       valmin=xstart,         valmax=xend,           valfmt="%1.1f", valinit=(xstart,xend),   color='lightblue')
slider_y    = RangeSlider(ax_y,         title_y,       valmin=ystart,         valmax=yend,            valfmt="%1.3f", valinit=(ystart,yend),   color='lightblue')
slider_zoom      = Slider(ax_zoom,      title_zoom,    valmin=0.01,           valmax=20.,            valfmt="%1.3f", valinit=zoom,            color='lightblue')
slider_alpha     = Slider(ax_alpha,     title_alpha,   valmin=0.01,           valmax=1.,             valfmt="%1.3f", valinit=myalpha,         color='lightblue')
slider_x_rv = RangeSlider(ax_x_rv,      title_x,       valmin=xstart_rv,      valmax=xend_rv,        valfmt="%1.1f", valinit=(xstart_rv,xend_rv), color='lightblue')
slider_y_rv = RangeSlider(ax_y_rv,      title_y,       valmin=ystart_rv,      valmax=yend_rv,        valfmt="%1.0f", valinit=(ystart_rv,yend_rv), color='lightblue')
slider_zoom_rv   = Slider(ax_zoom_rv,   title_zoom,    valmin=0.01,           valmax=20.,            valfmt="%1.3f", valinit=zoom_rv,         color='lightblue')
slider_alpha_rv  = Slider(ax_alpha_rv,  title_alpha,   valmin=0.01,           valmax=1.,             valfmt="%1.3f", valinit=myalpha_rv,      color='lightblue')
slider_t0        = Slider(ax_t0,        title_t0,      valmin=t0start-10.,    valmax=t0start+10.,    valfmt="%1.3f", valinit=params.t0,       color='g')
slider_t0fine    = Slider(ax_t0fine,    title_t0fine,  valmin=t0start-0.1,    valmax=t0start+0.1,    valfmt="%1.3f", valinit=params.t0,       color='g')
slider_per       = Slider(ax_per,       title_per,     valmin=0.1,            valmax=Pmax,           valfmt="%1.3f", valinit=params.per,      color='g')
slider_perfine   = Slider(ax_perfine,   title_perfine, valmin=params.per-0.5, valmax=params.per+0.5, valfmt="%1.3f", valinit=params.per,      color='g')
slider_aR        = Slider(ax_aR,        title_aR,      valmin=2.,             valmax=amax,           valfmt="%1.3f", valinit=params.a,        color='g')
slider_aAU       = Slider(ax_aAU,       title_aAU, valmin=2.*params2.Rstar*Rsun/AU, valmax=amax*params2.Rstar*Rsun/AU, valfmt="%1.3f", valinit=params2.aAU,        color='b')
slider_rhostar   = Slider(ax_rhostar,   title_rhostar, valmin=0.01,           valmax=10.,            valfmt="%1.3f", valinit=params2.rhostar, color='b')
slider_Rstar     = Slider(ax_Rstar,     title_Rstar,   valmin=0.1,            valmax=5.,             valfmt="%1.3f", valinit=params2.Rstar,   color='b')
slider_Mstar     = Slider(ax_Mstar,     title_Mstar,   valmin=0.08,           valmax=10.,            valfmt="%1.3f", valinit=params2.Mstar,   color='b')
slider_rhoplan   = Slider(ax_rhoplan,   title_rhoplan, valmin=0.,             valmax=10.,            valfmt="%1.3f", valinit=params2.rhoplan, color='b')
slider_Rplan     = Slider(ax_Rplan,     title_Rplan,   valmin=1,              valmax=100.,           valfmt="%1.2f", valinit=params2.Rplan,   color='b')
slider_Mplan     = Slider(ax_Mplan,     title_Mplan,   valmin=0,              valmax=1000.,          valfmt="%1.3f", valinit=params2.Mplan,   color='b')
slider_rp        = Slider(ax_rp,        title_rp,      valmin=0.01,           valmax=1.,             valfmt="%1.3f", valinit=params.rp,       color='g')
slider_b         = Slider(ax_b,         title_b,       valmin=0.,             valmax=2.,             valfmt="%1.3f", valinit=params2.b,       color='b')
slider_inc       = Slider(ax_inc,       title_inc,     valmin=70.,            valmax=90.,            valfmt="%1.3f", valinit=params.inc,      color='g')
slider_u1        = Slider(ax_u1 ,       title_u1,      valmin=0.,             valmax=2.,             valfmt="%1.3f", valinit=params.u[0],     color='g')
slider_u2        = Slider(ax_u2 ,       title_u2,      valmin=-1.,            valmax=1.,             valfmt="%1.3f", valinit=params.u[1],     color='g')
slider_ecc       = Slider(ax_ecc,       title_ecc,     valmin=0.,             valmax=1.,             valfmt="%1.3f", valinit=params.ecc,      color='g')
slider_w         = Slider(ax_w,         title_w,       valmin=0.,             valmax=360.,           valfmt="%1.3f", valinit=params.w,        color='g')
slider_F0        = Slider(ax_F0,        title_F0,      valmin=0.95,           valmax=1.05,           valfmt="%1.3f", valinit=params2.F0,      color='b')
slider_bg        = Slider(ax_bg,        title_bg,      valmin=-0.01,          valmax=0.01,           valfmt="%1.3f", valinit=params2.bg,      color='b')
slider_rv_K      = Slider(ax_rv_K,      title_K,       valmin=1.,             valmax=ymax_rv-ymin_rv,valfmt="%1.1f", valinit=params2.K,       color='g')
slider_rv_offset = Slider(ax_rv_offset, title_offset,  valmin=ymin_rv,        valmax=ymax_rv,        valfmt="%1.1f", valinit=params2.offset,  color='g')
slider_dur       = Slider(ax_dur,       title_dur,     valmin=0.,             valmax=24,             valfmt="%1.3f", valinit=params2.dur,     color='grey', dragging=False)
slider_depth     = Slider(ax_depth,     title_depth,   valmin=0.,             valmax=100,            valfmt="%1.3f", valinit=params2.depth,   color='grey', dragging=False)

#slider_xstart.label.set_size(8)
#slider_xend.label.set_size(8)
#slider_ystart.label.set_size(8)
#slider_yend.label.set_size(8)
slider_x.label.set_size(8)
slider_y.label.set_size(8)
slider_zoom.label.set_size(8)
slider_alpha.label.set_size(8)
slider_x_rv.label.set_size(8)
slider_y_rv.label.set_size(8)
slider_zoom_rv.label.set_size(8)
slider_alpha_rv.label.set_size(8)
slider_t0.label.set_size(8)
slider_t0fine.label.set_size(8)
slider_per.label.set_size(8)
slider_perfine.label.set_size(8)
slider_aR.label.set_size(8)
slider_aAU.label.set_size(8)
slider_rhostar.label.set_size(8)
slider_Rstar.label.set_size(8)
slider_Mstar.label.set_size(8)
slider_rhoplan.label.set_size(8)
slider_Rplan.label.set_size(8)
slider_Mplan.label.set_size(8)
slider_rp.label.set_size(8)
slider_b.label.set_size(8)
slider_inc.label.set_size(8)
slider_u1.label.set_size(8)
slider_u2.label.set_size(8)
slider_ecc.label.set_size(8)
slider_w.label.set_size(8)
slider_F0.label.set_size(8)
slider_bg.label.set_size(8)
slider_rv_K.label.set_size(8)
slider_rv_offset.label.set_size(8)
slider_dur.label.set_size(8)
slider_depth.label.set_size(8)

#slider_xstart.on_changed(slider_update_xstart)
#slider_xend.on_changed(slider_update_xend)
#slider_ystart.on_changed(slider_update_ystart)
#slider_yend.on_changed(slider_update_yend)
slider_x.on_changed(slider_update_x)
slider_y.on_changed(slider_update_y)
slider_zoom.on_changed(slider_update_zoom)
slider_alpha.on_changed(slider_update_alpha)
slider_x_rv.on_changed(slider_update_x_rv)
slider_y_rv.on_changed(slider_update_y_rv)
slider_zoom_rv.on_changed(slider_update_zoom_rv)
slider_alpha_rv.on_changed(slider_update_alpha_rv)
slider_t0.on_changed(slider_update_t0)
slider_t0fine.on_changed(slider_update_t0fine)
slider_per.on_changed(slider_update_per)
slider_perfine.on_changed(slider_update_perfine)
slider_aR.on_changed(slider_update_aR)
slider_aAU.on_changed(slider_update_aAU)
slider_rhostar.on_changed(slider_update_rhostar)
slider_Rstar.on_changed(slider_update_Rstar)
slider_Mstar.on_changed(slider_update_Mstar)
slider_rhoplan.on_changed(slider_update_rhoplan)
slider_Rplan.on_changed(slider_update_Rplan)
slider_Mplan.on_changed(slider_update_Mplan)
slider_rp.on_changed(slider_update_rp)
slider_b.on_changed(slider_update_b)
slider_inc.on_changed(slider_update_inc)
slider_u1.on_changed(slider_update_u1)
slider_u2.on_changed(slider_update_u2)
slider_ecc.on_changed(slider_update_ecc)
slider_w.on_changed(slider_update_w)
slider_F0.on_changed(slider_update_F0)
slider_bg.on_changed(slider_update_bg)
slider_rv_K.on_changed(slider_update_rv_K)
slider_rv_offset.on_changed(slider_update_rv_offset)


############################### plotting data

if lc_filename!='' :
    data1, = ax1.plot( tdata,
                       Fdata, c='b' ,zorder=0, alpha=myalpha, ms=5, **marker_style)
    data2, = ax2.plot( (tdata - params.t0 + params.per/2. ) % params.per - params.per/2.,
                       Fdata, c='b' ,zorder=0, alpha=myalpha, ms=5, **marker_style)
else :
    data1 = None
    data2 = None


if rv_filename!='' :
    # data_rv_1 = ax3.errorbar( x_rv,
    #                           y_rv, yerr=5*err_rv, c='blue', alpha=0.5,  **marker_style)
    # data_rv_2 = ax4.errorbar( (x_rv - params.t0 + params.per/2. ) % params.per - params.per/2.,
    #                           y_rv, yerr=5*err_rv, c='b' ,zorder=0, alpha=my_alpha, **marker_style)
    data_rv_1, = ax3.plot( x_rv,
                         y_rv, c='b', zorder=0, alpha=myalpha_rv, ms=5, **marker_style)
    data_rv_2, = ax4.plot( (x_rv - params.t0 + params.per/2. ) % params.per - params.per/2.,
                         y_rv, c='b', zorder=0, alpha=myalpha_rv, ms=5, **marker_style)
else :
    data_rv_1 = None
    data_rv_2 = None


line1,     = ax1.plot([params.t0,params.t0], [0,2],        c='g' ,zorder=1, alpha=0.5, **line_style)
line2,     = ax2.plot([0.,0.],               [0,2],        c='g' ,zorder=1, alpha=0.5, **line_style)

line_rv_1, = ax3.plot([params.t0,params.t0], [-1000,1000], c='g' ,zorder=1, alpha=0.5, **line_style)
line_rv_2, = ax4.plot([0.,0.],               [-1000,1000], c='g' ,zorder=1, alpha=0.5, **line_style)


############################### model curves

t1 = np.linspace(tmin, tmax, prec)
f1 = calc_lc(t1,params,params2)

t2 = np.linspace(params.t0-params.per/2/zoom, params.t0+params.per/2/zoom, prec)
f2 = calc_lc(t2,params,params2)

model1, = ax1.plot(t1,           f1, color='r', **line_style, zorder=3, alpha=0.8)
model2, = ax2.plot(t2-params.t0, f2, color='r', **line_style, zorder=3, alpha=0.8)

params2.depth = calc_depth(params,params2)
slider_depth.set_val(params2.depth)

t1_rv = np.linspace(xstart_rv, xend_rv, prec)
rv1 = calc_rv(t1_rv,params,params2)

t2_rv = np.linspace(params.t0-params.per/2, params.t0+params.per/2, prec)
rv2 = calc_rv(t2_rv,params,params2)

model_rv_1, = ax3.plot(t1_rv,             rv1, color='r', alpha=0.3)
model_rv_2, = ax4.plot(t2_rv - params.t0, rv2, color='r', alpha=0.8)

if lc_filename!='' :
  model1.figure.canvas.manager.set_window_title(lc_filename)
else :
  model1.figure.canvas.manager.set_window_title(rv_filename)
# model1.figure.canvas.manager.set_window_title(args.files[0])
# model_rv_1.figure.canvas.manager.set_window_title(filename)

########## transit sketch
circle_star   = plt.Circle((0.,0.),params2.Rstar,color='orange',fill=True,zorder=0)
circle_planet = plt.Circle((0.,params2.b*params2.Rstar),params.rp*params2.Rstar,color='b',fill=True,alpha=1,zorder=1)
circ_star     = ax5.add_patch(circle_star)
circ_planet   = ax5.add_patch(circle_planet)
size          = 1.15*(params2.Rstar)
ax5.axis([-size,size,-size,size])

########## orbit sketch
phi=np.linspace(0,2*np.pi,100)
x = calc_orbit_x(phi,params)*params2.Rstar*Rsun/AU
y = calc_orbit_y(phi,params)*params2.Rstar*Rsun/AU
orbit_p, = ax6.plot(x,y,color='grey',lw=0.5,zorder=1)
orbit_l, = ax6.plot([0,0],[-2.*params.a*params2.Rstar*Rsun/AU,-params.rp*params2.Rstar*Rsun/AU],color='g',lw=0.5,ls='--',zorder=0)
circle_star_orbit   = plt.Circle((0.,0.),params2.Rstar*Rsun/AU,color='orange',fill=True)
circle_planet_orbit = plt.Circle((0.,calc_orbit_y(-np.pi/2,params)*params2.Rstar*Rsun/AU),params.rp*params2.Rstar*Rsun/AU,color='b',fill=True,alpha=1,zorder=2)
circ_star_orbit     = ax6.add_patch(circle_star_orbit)
circ_planet_orbit   = ax6.add_patch(circle_planet_orbit)
size                = 1.15*(params.a*params2.Rstar*Rsun/AU*(1-params.ecc**2))/(1-params.ecc)
ax6.axis([-size,size,-size,size])

if args.save!='' :
    plt.savefig(args.save)
    #plt.savefig(lc_filename+'.pdf')

print('input parameters:')
print_all(params,params2)

plt.show()

print('output parameters:')
print_all(params,params2)

# print( '  norm=',norm)
# print( '  REB [Rsun] = ',params2.Rplan*Rearth/Rsun)


