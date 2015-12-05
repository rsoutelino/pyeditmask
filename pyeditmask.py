#!/usr/bin/env python
######################################################
## Edits ROMS masks using a GUI
## Nov 2014
## rsoutelino@gmail.com
######################################################
import os
import wx

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Navbar
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import scipy.io as sp
import netCDF4 as nc

# TO-DO LIST: ====================================================
#   - improve point selection based in find_lower_left_node
#   - create better icons for mask/unmask area
#   - resolve untoggle/toggle between mask/unmask functions
#   - add support to other models (POM, SWAN, WW3)
#   - move matplotlib toolbar to the lower part
#   - add a wx.TaskBarIcon to show up on the unity launcher
#   - display local depth of the pixel we are looking at
#   - transform mask/unmask and mask_area and unmask_area in the same
#        function, and figure out how to decide wether to mask or unmask
# ================================================================

# NICE TIP TO DEBUG THIS PROGRAM: ================================
#   - comment out app.MainLoop at the last line of this script
#   - ipython --gui=wx
#   - run pyeditmask.py
#   - trigger the events and check out the objects in the shell
# ================================================================


global currentDirectory
currentDirectory = os.getcwd()

PROJECT_DIR = os.path.abspath(os.path.dirname(__file__))
DEFAULT_VMIN = 0 
DEFAULT_VMAX=1.5 
DEFAULT_CMAP = plt.cm.BrBG
DEFAULT_DEPTH_FOR_LAND = -50


class WW3Grid(object):
    """ 
    Stores and manipulates WW3 gridgen output 
    """
    def __init__(self,filename):
        self.filename = filename
        self.read_meta()
        self.read_grid()
        # self.ncfile = nc.Dataset(filename, mode='r+')
        # self.lonr  = self.ncfile.variables['lon_rho'][:]
        # self.latr  = self.ncfile.variables['lat_rho'][:]
        # self.lonu  = self.ncfile.variables['lon_u'][:]
        # self.latu  = self.ncfile.variables['lat_u'][:]
        # self.lonv  = self.ncfile.variables['lon_v'][:]
        # self.latv  = self.ncfile.variables['lat_v'][:]
        # self.h     = self.ncfile.variables['h'][:]
        # self.maskr = self.ncfile.variables['mask_rho'][:]
        # self.masku = self.ncfile.variables['mask_u'][:]
        # self.maskv = self.ncfile.variables['mask_v'][:]     

    def read_meta(self):
        """
        Reads meta data
        """
        tmp=['dummy']
        f=open(self.filename)
        while tmp[0] != "'RECT'":
            tmp=f.readline().split()
        nx,ny=np.array(f.readline().split()).astype(float)
        res=np.array(f.readline().split()).astype(float)
        x0 = np.array(f.readline().split()).astype(float)
        dlon,dlat=res[0]/res[2],res[1]/res[2]
        lon0,lat0=x0[0]/x0[2],x0[1]/x0[2]
        self.lonr = np.linspace(lon0, lon0+(dlon)*ny, ny, endpoint=True)
        self.latr = np.linspace(lat0, lat0+(dlat)*nx, nx, endpoint=True)

    def read_grid(self,extention='mask'):
        """
        Reads filename and returns data in numpy array
        filename -- filename
        """
        dir,meta = os.path.split(self.filename)
        filename = os.path.join(dir,meta.split('.')[0]+'.'+extention)
        if extention == 'obst':
            temp = np.loadtxt(filename)
            shape = temp.shape
            dum = shape[0]/2
            self.maskr = np.array([temp[:dum,],temp[dum:,]])
        else:
            self.maskr = np.loadtxt(filename).transpose()


# ROMS related objects ---------------------------------------------
class RomsGrid(object):
    """ 
    Stores and manipulates netcdf ROMS grid file information
    """
    def __init__(self,filename):
        self.filename = filename    
        self.ncfile = nc.Dataset(filename, mode='r+')
        self.lonr  = self.ncfile.variables['lon_rho'][:]
        self.latr  = self.ncfile.variables['lat_rho'][:]
        self.lonu  = self.ncfile.variables['lon_u'][:]
        self.latu  = self.ncfile.variables['lat_u'][:]
        self.lonv  = self.ncfile.variables['lon_v'][:]
        self.latv  = self.ncfile.variables['lat_v'][:]
        self.h     = self.ncfile.variables['h'][:]
        self.maskr = self.ncfile.variables['mask_rho'][:]
        self.masku = self.ncfile.variables['mask_u'][:]
        self.maskv = self.ncfile.variables['mask_v'][:]     


def uvp_mask(rfield):
    Mp, Lp = rfield.shape
    M      = Mp - 1
    L      = Lp - 1

    vfield = rfield[0:M,:] * rfield[1:Mp,:]
    ufield = rfield[:,0:L] * rfield[:,1:Lp]
    pfield = ufield[0:M,:] * ufield[1:Mp,:]

    return ufield, vfield, pfield
# -------------------------------------------------------------------


class App(wx.App):
    def OnInit(self):
        self.frame = Interface("PyEditMask 0.1.0", size=(1024,800))
        self.frame.Show()
        return True


class Interface(wx.Frame):
    def __init__(self, title=wx.EmptyString, pos=wx.DefaultPosition, 
                       size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE,
                       *args, **kwargs):
        wx.Frame.__init__(self, None, -1, "PyEditMask 0.1.0", pos=pos, 
                          size=size, style=style, *args, **kwargs)
        
        # Initializing toolbar
        self.toolbar = MainToolBar(self)


        # BASIC LAYOUT OF THE NESTED SIZERS ======================
        panel1 = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel.SetBackgroundColour("WHITE")

        # BOX 1 is the main sizer
        box1 = wx.BoxSizer(wx.HORIZONTAL)
        box1.Add(panel1, 1, wx.EXPAND)
        box1.Add(mplpanel, 15, wx.EXPAND)

        # BOX 2 is the inner sizer of the left big control panel
        box2 = wx.BoxSizer(wx.VERTICAL)

        # BOX 3 is the sizer of the right big parent panel(panel1), the one that will
        #    serve as base for two child panels which will hold
        #    the two matplotlib canvas's
        box3 = wx.BoxSizer(wx.VERTICAL)

        # panel 1 content ========================================
        main_label = wx.StaticText(panel1, label=" ")
        box2.Add(main_label, proportion=0, flag=wx.CENTER)

        # set_land = wx.Button(panel1, label="Set Land", style=wx.ID_CANCEL)
        # box2.Add(set_land, proportion=0, flag=wx.CENTER)
        # set_land.Bind(wx.EVT_BUTTON, self.onSetLand)

        # set_water = wx.Button(panel1, label="Set Water", style=wx.ID_CANCEL)
        # box2.Add(set_water, proportion=0, flag=wx.CENTER)
        # set_water.Bind(wx.EVT_BUTTON, self.onSetWater)   

        # mplpanel content ========================================
        self.mplpanel = SimpleMPLCanvas(mplpanel)
        box3.Add(self.mplpanel.canvas, 1, flag=wx.CENTER) 

        # FINAL LAYOUT CONFIGURATIONS ============================
        self.SetAutoLayout(True)
        panel1.SetSizer(box2)
        # panel2.SetSizer(box4)
        mplpanel.SetSizer(box3)

        self.SetSizer(box1)

        self.InitMenu()
        self.Layout()
        self.Centre()
        # self.ShowModal()


    def InitMenu(self):
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_OPEN, u'&Open ROMS grid file')
        fileMenu.Append(wx.ID_OPEN, u'&Open coastline file')
        fileMenu.Append(wx.ID_SAVE, '&Save grid')
        fileMenu.AppendSeparator()

        qmi = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+W')
        opf = wx.MenuItem(fileMenu, wx.ID_OPEN, '&Open\tCtrl+O')
        opc = wx.MenuItem(fileMenu, wx.ID_OPEN, '&Open\tCtrl+O+C')
        svf = wx.MenuItem(fileMenu, wx.ID_SAVE, '&Save\tCtrl+S')
        fileMenu.AppendItem(qmi)
        # fileMenu.AppendItem(svf)

        self.Bind(wx.EVT_MENU, self.OnQuit, qmi)
        self.Bind(wx.EVT_MENU, self.toolbar.OnLoadGrid, opf)
        self.Bind(wx.EVT_MENU, self.toolbar.OnLoadCoastline, opc)
        self.Bind(wx.EVT_MENU, self.toolbar.OnSaveGrid, svf)

        menubar.Append(fileMenu, u'&PyEditMask')
        self.SetMenuBar(menubar)


    def OnQuit(self, e):
        """Fecha o programa"""
        self.Close()
        self.Destroy()

    def OnCloseWindow(self, e):
        self.Destroy()


class SimpleMPLCanvas(object):
    """docstring for SimpleMPLCanvas"""
    def __init__(self, parent):
        super(SimpleMPLCanvas, self).__init__()
        self.parent = parent
        self.plot_properties()
        self.make_navbar()
        
    def make_navbar(self):
        self.navbar = Navbar(self.canvas)   
        self.navbar.SetPosition(wx.Point(0,0)) # this is not working !!


    def plot_properties(self):
        # Create matplotlib figure
        self.fig = Figure(facecolor='w', figsize=(12,8))
        self.canvas = FigureCanvas(self.parent, -1, self.fig)
        
        self.ax   = self.fig.add_subplot(111)
        # tit = self.ax1.set_title("ROMS mask_rho", fontsize=12, fontweight='bold')
        # tit.set_position([0.9, 1.05])


class MainToolBar(object):
    def __init__(self, parent):
        self.currentDirectory = os.getcwd()
        self.parent = parent
        self.toolbar = parent.CreateToolBar(style=1, winid=1,
                                            name="Toolbar")
        self.tools_params ={ 
            'load_grid': (load_bitmap('grid.png'), u"Load grid",
                        "Load ocean_grd.nc ROMS grid netcdf file"),
            'load_coastline': (load_bitmap('coast.png'), u"Load coastline",
                        "Load *.mat coastline file [lon / lat poligons]"),
            'save_grid': (load_bitmap('save.png'), u"Apply and save",
                        "Save changes to ocean_grd.nc ROMS grid netcdf file"),
            'set_land': (load_bitmap('land.png'), u"Set land",
                        "Set grid point to land"),
            'set_land_area': (load_bitmap('land_area.png'), u"Set land area",
                        "Set poligon area to land"),
            'set_water': (load_bitmap('water.png'), u"Set water",
                        "Set grid point to water"),
            'set_water_area': (load_bitmap('water_area.png'), u"Set water area",
                        "Set poligon area to water"),
            'settings': (load_bitmap('settings.png'), u"PyEditMask settings",
                        "PyEditMask configurations"),
            'quit': (load_bitmap('exit.png'), u"Quit",
                        "Quit PyEditMask"),
        }
        
        self.createTool(self.toolbar, self.tools_params['load_grid'], 
                        self.OnLoadGrid)
        self.createTool(self.toolbar, self.tools_params['load_coastline'], 
                        self.OnLoadCoastline)
        self.createTool(self.toolbar, self.tools_params['save_grid'], 
                        self.OnSaveGrid)
        
        self.toolbar.AddSeparator()

        self.mask_tool = self.createTool(self.toolbar, self.tools_params['set_land'], 
                                         self.OnSetLand, isToggle=True)
        self.mask_area_tool = self.createTool(self.toolbar, 
                                              self.tools_params['set_land_area'], 
                                              self.OnSetLandArea, isToggle=True)
        self.unmask_tool = self.createTool(self.toolbar, self.tools_params['set_water'], 
                                           self.OnSetWater, isToggle=True)
        self.unmask_area_tool = self.createTool(self.toolbar, 
                                                self.tools_params['set_water_area'], 
                                                self.OnSetWaterArea, isToggle=True)

        self.toolbar.AddSeparator()

        self.createTool(self.toolbar, self.tools_params['settings'], 
                        self.OnSettings)
        self.createTool(self.toolbar, self.tools_params['quit'], 
                        self.parent.OnQuit)

        self.toolbar.Realize()


    def createTool(self, parent, params, evt, isToggle=False):
        tool = parent.AddTool(wx.NewId(), bitmap=params[0], shortHelpString=params[1],
                    longHelpString=params[2], isToggle=isToggle)
        self.parent.Bind(wx.EVT_TOOL, evt, id=tool.GetId())
        return tool


    def OnLoadGrid(self, evt):
        openFileDialog = wx.FileDialog(self.parent, "Open gridgen meta file [*.meta]",
                                       "/source/ww3/tools/gridgen/output/", " ",
                                       "meta files (*.meta)|*.meta",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...

        filename = openFileDialog.GetPath()
        grd = WW3Grid(filename)

        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        self.pcolor = ax.pcolormesh(grd.lonr, grd.latr, grd.maskr, 
                                   vmin=DEFAULT_VMIN, vmax=DEFAULT_VMAX, 
                                   cmap=DEFAULT_CMAP)
        ax.plot(grd.lonr, grd.latr, 'k', alpha=0.2)
        ax.plot(grd.lonr.transpose(), grd.latr.transpose(), 'k', alpha=0.2)
        ax.set_xlim([grd.lonr.min(), grd.lonr.max()])
        ax.set_ylim([grd.latr.min(), grd.latr.max()])
        ax.set_aspect('equal')

        mplpanel.canvas.draw()
        self.grd = grd
        self.grd.hmin = grd.ncfile.variables['h'][:].min()


    def OnLoadCoastline(self, evt):
        openFileDialog = wx.FileDialog(self.parent, "Open coastline file - NETCDF",
                                       "/metocean/roms/data", " ",
                                       "NETCDF files (*.nc)|*.nc",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...

        filename = openFileDialog.GetPath()
        coast = nc.Dataset(filename)
        lon, lat = coast.variables['lon'][:], coast.variables['lat'][:]

        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        ax.plot(lon, lat, 'k')

        try:
            ax.set_xlim([self.grd.lonr.min(), self.grd.lonr.max()])
            ax.set_ylim([self.grd.latr.min(), self.grd.latr.max()])
        except AttributeError: # just in case a grid was not loaded before
            ax.set_xlim([np.nanmin(lon), np.nanmax(lon)])
            ax.set_ylim([np.nanmin(lat), np.nanmax(lat)])
        
        ax.set_aspect('equal')
        mplpanel.canvas.draw()


    def OnSaveGrid(self, evt):
        maskr = self.grd.maskr
        [masku, maskv, maskp] = uvp_mask(maskr)
        self.grd.ncfile.variables['mask_rho'][:] = maskr
        self.grd.ncfile.variables['mask_u'][:]   = masku
        self.grd.ncfile.variables['mask_v'][:]   = maskv
        self.grd.ncfile.variables['mask_psi'][:] = maskp
        self.grd.ncfile.variables['h'][:] = self.grd.h
        self.grd.ncfile.sync()


    def OnSetLand(self, evt):
        mplpanel = app.frame.mplpanel

        if self.mask_tool.IsToggled():
            self.cid = mplpanel.canvas.mpl_connect('button_press_event', self.mask)
        else:
            mplpanel.canvas.mpl_disconnect(self.cid)


    def OnSetLandArea(self, evt):
        mplpanel = app.frame.mplpanel

        if self.mask_area_tool.IsToggled():
            self.cid = mplpanel.canvas.mpl_connect('button_press_event', 
                                                    self.mask_area)
        else:
            mplpanel.canvas.mpl_disconnect(self.cid)


    def OnSetWater(self, evt):
        mplpanel = app.frame.mplpanel

        if self.unmask_tool.IsToggled():
            self.cid = mplpanel.canvas.mpl_connect('button_press_event', self.unmask)
        else:
            mplpanel.canvas.mpl_disconnect(self.cid)


    def OnSetWaterArea(self, evt):
        mplpanel = app.frame.mplpanel

        if self.unmask_area_tool.IsToggled():
            self.cid = mplpanel.canvas.mpl_connect('button_press_event', 
                                                    self.unmask_area)
        else:
            mplpanel.canvas.mpl_disconnect(self.cid)


    def OnSettings(self, evt):
        pass


    def mask(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
        line, col = find_lower_left_node(self.grd.lonr, self.grd.latr, x, y)
        self.grd.maskr[line, col] = 0 # assigning new value
        self.grd.h[line, col] = self.grd.hmin
        # refilling with new value
        ax.pcolormesh(self.grd.lonr[line:line+2, col:col+2], 
                      self.grd.latr[line:line+2, col:col+2], 
                      self.grd.maskr[line:line+2, col:col+2], 
                      vmin=DEFAULT_VMIN, vmax=DEFAULT_VMAX, cmap=DEFAULT_CMAP)
        mplpanel.canvas.draw()


    def mask_area(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
        button = evt.button
        if button == 1:
            p = ax.plot(x, y, 'ro')
            try:
                self.points.append(p)
                self.area.append( (x, y) )
            except AttributeError:
                self.points = [p]
                self.area = [ (x, y) ]

            mplpanel.canvas.draw()

        elif button == 3:
            grd = self.grd
            path = Path( self.area )
            a, b = grd.lonr.shape
            for i in range(a):
                for j in range(b):
                    if path.contains_point( [grd.lonr[i, j], 
                                             grd.latr[i, j] ] ) == 1:
                        grd.maskr[i,j] = 0
                        grd.h[i,j] = grd.hmin

            ax.clear()
            self.pcolor = ax.pcolormesh(grd.lonr, grd.latr, grd.maskr, 
                                   vmin=DEFAULT_VMIN, vmax=DEFAULT_VMAX, 
                                   cmap=DEFAULT_CMAP)
            ax.plot(grd.lonr, grd.latr, 'k', alpha=0.2)
            ax.plot(grd.lonr.transpose(), grd.latr.transpose(), 'k', alpha=0.2)
            ax.set_xlim([grd.lonr.min(), grd.lonr.max()])
            ax.set_ylim([grd.latr.min(), grd.latr.max()])
            ax.set_aspect('equal')
            mplpanel.canvas.draw()
            del self.points, self.area



    def unmask(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
        line, col = find_lower_left_node(self.grd.lonr, self.grd.latr, x, y)
        self.grd.maskr[line, col] = 1 # assigning new value
        self.grd.h[line, col] = self.grd.hmin
        # refilling with new value
        ax.pcolormesh(self.grd.lonr[line:line+2, col:col+2], 
                      self.grd.latr[line:line+2, col:col+2], 
                      self.grd.maskr[line:line+2, col:col+2], 
                      vmin=DEFAULT_VMIN, vmax=DEFAULT_VMAX, cmap=DEFAULT_CMAP)
        mplpanel.canvas.draw()


    def unmask_area(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax: return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
        button = evt.button
        if button == 1:
            p = ax.plot(x, y, 'ro')
            try:
                self.points.append(p)
                self.area.append( (x, y) )
            except AttributeError:
                self.points = [p]
                self.area = [ (x, y) ]

            mplpanel.canvas.draw()

        elif button == 3:
            grd = self.grd
            path = Path( self.area )
            a, b = grd.lonr.shape
            for i in range(a):
                for j in range(b):
                    if path.contains_point( [grd.lonr[i, j], 
                                             grd.latr[i, j] ] ) == 1:
                        grd.maskr[i,j] = 1
                        grd.h[i,j] = grd.hmin
                        
            ax.clear()
            self.pcolor = ax.pcolormesh(grd.lonr, grd.latr, grd.maskr, 
                                   vmin=DEFAULT_VMIN, vmax=DEFAULT_VMAX, 
                                   cmap=DEFAULT_CMAP)
            ax.plot(grd.lonr, grd.latr, 'k', alpha=0.2)
            ax.plot(grd.lonr.transpose(), grd.latr.transpose(), 'k', alpha=0.2)
            ax.set_xlim([grd.lonr.min(), grd.lonr.max()])
            ax.set_ylim([grd.latr.min(), grd.latr.max()])
            ax.set_aspect('equal')
            mplpanel.canvas.draw()
            del self.points, self.area



def find_lower_left_node(x, y, x0, y0, n=4):
    # need to improve this, not very accurate yet
    dx = np.abs(x - x0); dx = dx / dx.max()
    dy = np.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    
    line, col, lola = [], [], []

    for k in range(n):
        fn = np.where(dn == dn.min())
        f1, f2 = int(fn[0]), int(fn[1])
        line.append(f1)
        col.append(f2)
        lola.append(x[f1, f2] + y[f1, f2])
        dn[f1, f2] = 1e20

    lola = np.array(lola)
    f = np.where(lola == lola.min())[0][0]

    line = line[f]
    col = col[f]

    return line, col


def load_bitmap(filename, direc=None):
    """
    Load a bitmap file from the ./icons subdirectory. 
    The filename parameter should not
    contain any path information as this is determined automatically.

    Returns a wx.Bitmap object
    copied from matoplotlib resources
    """

    if not direc:
        basedir = os.path.join(PROJECT_DIR,'icons')
    else:
        basedir = os.path.join(PROJECT_DIR, direc)

    bmpFilename = os.path.normpath(os.path.join(basedir, filename))
    if not os.path.exists(bmpFilename):
        raise IOError('Could not find bitmap file "%s"; dying'%bmpFilename)

    bmp = wx.Bitmap(bmpFilename)
    return bmp


if __name__ == "__main__":
    app = App(False)
    app.MainLoop()























