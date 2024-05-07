"""This code reads in an input files containing the wells, bimolecular products,
transition states andbarrierless reactions and creates a PES plot
"""
import os
import sys
import math

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pylab as plt  # translate into pyplot.
import matplotlib.image as mpimg
import numpy as np
import numpy.linalg as la
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw.cairoCanvas import Canvas
from openbabel import pybel
from PIL import Image
import networkx as nx
from pyvis import network as net

from pesviewer.gen_resonant_structs import gen_reso_structs

pybel.ob.obErrorLog.SetOutputLevel(0)
# contains all the options for this PES
options = {}

# global parameters for the plot
xlow = 0.0  # lowest x value on x axis
xhigh = 0.0  # highest x value on x axis
xmargin = 0.0  # margin on x axis
ylow = 0.0  # lowest y value on x axis
yhigh = 0.0  # highest y value on x axis
ymargin = 0.0  # margin on y axis
xlen = 1.0  # length of horizontal lines per st pt

wells = []  # list of wells
bimolecs = []  # list of bimolecular products
tss = []  # list of transition states
barrierlesss = []  # list of barrierless reactions

# text dictionary: key is the chemical structure, value is the text
textd = {}
# lines dictionary: key is the chemical structure,
# value is a list of lines for that structure
linesd = {}
# figures dictionary: key is the chemical structure,
# value is the figure of that structure
imgsd = {}
# extents of the images
extsd = {}


class dragimage(object):
    """
    Class to drag an image
    """
    def __init__(self, figure=None):
        if figure is None:
            figure = plt.gcf()
        # simple attibute to store the dragged text object
        self.struct = None
        self.img = None
        # Connect events and callbacks
        figure.canvas.mpl_connect("button_press_event", self.pick_image)
        figure.canvas.mpl_connect("button_release_event", self.release_image)
        figure.canvas.mpl_connect("motion_notify_event", self.move_image)
    # end def

    def pick_image(self, event):
        for key in imgsd.keys():
            self.struct = None
            self.img = None
            if (imgsd[key].get_extent()[0] < event.xdata and
                    imgsd[key].get_extent()[1] > event.xdata and
                    imgsd[key].get_extent()[2] < event.ydata and
                    imgsd[key].get_extent()[3] > event.ydata):
                self.struct = key
                self.img = imgsd[key]
                self.current_pos = (event.xdata, event.ydata)
                break
            # end if
        # end for
    # end def

    def move_image(self, event):
        if self.img is not None:
            old_extent = self.img.get_extent()
            xchange = event.xdata-self.current_pos[0]
            ychange = event.ydata-self.current_pos[1]
            extent_change = (xchange, xchange, ychange, ychange)
            extent = [old_extent[i] + extent_change[i] for i in range(0, 4)]
            self.img.set_extent(extent=extent)
            self.current_pos = (event.xdata, event.ydata)
            plt.draw()
    # end def

    def release_image(self, event):
        if self.img is not None:
            self.struct = None
            self.img = None
            save_im_extent()
            # end if
        # end if
    # end def
# end class


class selecthandler(object):
    """
    Class to select and move stationary points, which highlight its reactions
    and in the future (TODO) renders the 3D images
    """
    def __init__(self, figure=None):
        if figure is None:
            figure = plt.gcf()
        self.struct = None  # stationary point that is selected
        figure.canvas.mpl_connect("button_press_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event",
                                  self.on_release_event)
        figure.canvas.mpl_connect("motion_notify_event",
                                  self.motion_notify_event)
    # end def

    def on_pick_event(self, event):
        self.struct = None
        # TODO: find more efficient way to iterate
        # all stationary points in one loop?
        # create a new list?
        for w in wells:
            if self.is_close(w, event):
                self.struct = w
                highlight_structure(self.struct)
            # end if
        # end for
        for b in bimolecs:
            if self.is_close(b, event):
                self.struct = b
                highlight_structure(self.struct)
            # end if
        # end for
        for t in tss:
            if self.is_close(t, event):
                self.struct = t
                highlight_structure(self.struct)
            # end if
        # end for
        if self.struct is None:
            highlight_structure()
        # end if
    # end def

    def motion_notify_event(self, event):
        if self.struct is not None:
            # a stationary point got selected
            # current position of the stationary point
            old_pos = (self.struct.x, self.struct.y)
            self.struct.x = event.xdata  # set the new position
            # move all the elements(image, text and lines)
            updateplot(self.struct, (event.xdata-old_pos[0]))
            self.current_pos = old_pos
    # end def

    def on_release_event(self, event):
        if self.struct is not None:
            self.struct = None
            # save the x-values of the startionary points to a file
            save_x_values()
            # save the image extents (x and y coordinates) to a file
            save_im_extent()
        # end if
        return True
    # end def

    def is_close(self, struct, event):
        """
        An event is close if it comes within 2% of the stationary point
        both in the x and in the y direction
        """
        xc = math.fabs(event.xdata - struct.x) < (xhigh-xlow)*0.02
        yc = math.fabs(event.ydata - struct.y) < (yhigh-ylow)*0.02
        return xc and yc
# end class


class line:
    """
    A line contains information about the line on the graph
    it is either a line between a reactant and ts, between a ts and product
    or between a reactant and product (for barrierless reactions)
    the chemstructs are the reactant and ts (f orward),
    ts and product (reverse), or reactant and product (barrierless)
    """
    def __init__(self, x1, y1, x2, y2, chemstruct=None, col='black'):
        if x1 <= x2:
            self.xmin = x1
            self.y1 = y1  # y1 corresponds to xmin
            self.xmax = x2
            self.y2 = y2  # y2 corresponds to xmax
        else:
            self.xmin = x2
            self.y1 = y2  # y1 corresponds to xmin
            self.xmax = x1
            self.y2 = y1  # y2 corresponds to xmax
        # end if
        
        
        if x1 == x2 or y1 == y2:
            self.straight_line = True
            self.coeffs = []
        else:
            self.straight_line = False
            self.coeffs = get_polynomial(self.xmin, self.y1,
                                         self.xmax, self.y2)
        # end if
        if chemstruct is None:
            self.chemstruct = []
        else:
            self.chemstruct = chemstruct
        self.color = col
    # end def
# end class


class well:
    # Well class, contains the name, smiles and energy of a well
    def __init__(self, name, energy, smi=None, energy2=None):
        self.name = name
        self.energy = convert_units(energy)
        self.smi = smi
        self.x = 0.
        self.y = 0.
        self.xyz_files = []
        self.energy2 = energy2
        fn = 'xyz/{name}.xyz'.format(name=name)
        if os.path.exists(fn):
            self.xyz_files.append(fn)
    # end def
# end class


class bimolec:
    # Bimolec class, contains the name,
    # both smiles and energy of a bimolecular product
    def __init__(self, name, energy, smi=None, energy2=None):
        self.name = name
        self.energy = convert_units(energy)
        if smi is None:
            self.smi = []
        else:
            self.smi = smi
        self.x = 0.
        self.y = 0.
        self.xyz_files = []
        # this bimolecular product is placed on the right side of the graph
        self.right = False
        self.energy2 = energy2
        i = 1
        fn = 'xyz/{name}{i}.xyz'.format(name=name, i=i)
        while os.path.exists(fn):
            self.xyz_files.append(fn)
            i += 1
            fn = 'xyz/{name}{i}.xyz'.format(name=name, i=i)
        # end for
    # end def
# end class


class ts:
    """
    TS class, contains the name, the names of the
    reactant and product and the energy of the ts
    """
    def __init__(self, name, names, energy, col='black'):
        self.name = name
        self.energy = convert_units(energy)
        self.color = col
        self.xyz_files = []
        fn = 'xyz/{name}.xyz'.format(name=name)
        if os.path.exists(fn):
            self.xyz_files.append(fn)
        self.reactant = next((w for w in wells if w.name == names[0]), None)
        if self.reactant is None:
            list = (b for b in bimolecs if b.name == names[0])
            self.reactant = next(list, None)
        if self.reactant is None:
            e = exceptions.not_recognized('reactant', names[0], name)
            raise Exception(e)
        self.product = next((w for w in wells if w.name == names[1]), None)
        if self.product is None:
            list = (b for b in bimolecs if b.name == names[1])
            self.product = next(list, None)
        if self.product is None:
            e = exceptions.not_recognized('product', names[1], name)
            raise Exception(e)
        self.lines = []
        self.x = 0.
        self.y = 0.
    # end def
# end class


class barrierless:
    """
    Barrierless class, contains the name and the
    names of the reactant and product
    """
    def __init__(self, name, names, col='black'):
        self.name = name
        self.xyz_files = []
        self.color = col
        fn = f'xyz/{name}.xyz'
        if os.path.exists(fn):
            self.xyz_files.append(fn)
        self.reactant = next((w for w in wells if w.name == names[0]), None)
        if self.reactant is None:
            list = (b for b in bimolecs if b.name == names[0])
            self.reactant = next(list, None)
        if self.reactant is None:
            e = exceptions.not_recognized('reactant', names[0], name)
            raise Exception(e)
        self.product = next((w for w in wells if w.name == names[1]), None)
        if self.product is None:
            list = (b for b in bimolecs if b.name == names[1])
            self.product = next(list, None)
        if self.product is None:
            e = exceptions.not_recognized('product', names[1], name)
            raise Exception(e)
        self.energy = convert_units(self.product.energy)
        self.line = None
    # end def
# end class


class exceptions(object):
    """
    Class that stores the exception messages
    """
    def not_recognized(role, ts_name, species_name):
        s = 'Did not recognize {role}'.format(role=role)
        s += ' {prod} '.format(prod=species_name)
        s += 'for the transition state {ts}'.format(ts=ts_name)
        return s


def get_polynomial(x1, y1, x2, y2):
    """
    Method fits a third order polynomial through two points as such
    that the derivative in both points is zero
    This method should only be used if x1 is not equal to x2
    """
    if x1 == x2:
        print('Error, cannot fit a polynomial if x1 equals x2')
        sys.exit()
    else:
        y = np.matrix([[y1], [y2], [0], [0]])
        x = np.matrix([[x1**3, x1**2, x1, 1],
                      [x2**3, x2**2, x2, 1],
                      [3*x1**2, 2*x1, 1, 0],
                      [3*x2**2, 2*x2, 1, 0]])
        xinv = la.inv(x)
        a = np.dot(xinv, y)
        return np.transpose(a).tolist()[0]
# end def


def read_input(fname):
    """
    Method to read the input file
    """
    if not os.path.exists(fname):  # check if file exists
        raise FileNotFoundError(fname + ' does not exist')
    # end if
    with open(fname, 'r') as f:
        input_str = f.read()
    
    input_str = input_str.replace('\r\n', '\n').replace('\r', '\n').replace('\n\n', '\n')
    input_dict = get_sd_prop(input_str)
    options['id'] = input_dict['id'][0]
    # by default, print the graph title
    options['title'] = 1
    # default units
    options['units'] = 'kJ/mol'
    # units to display the results, defaults to whatever is set in units
    options['display_units'] = None
    # rounding of values on plot
    options['rounding'] = 1
    # shifting energies (in original units)
    options['energy_shift'] = 0.
    # use xyz by default, put 0  to switch off
    options['use_xyz'] = 1
    # no rescale as default, put the well or
    # bimolecular name here to rescale to that value
    options['rescale'] = 0
    # default figure height
    options['fh'] = 9.
    # default figure width
    options['fw'] = 18.
    # scale factor for figures
    options['fs'] = 1.
    # change the linewidth of the traditional Pot. vs Reac. Coord. plot.
    options['lw'] = 1.5
    # default margin on the x and y axis
    options['margin'] = 0.2
    # default dpi of the molecule figures
    options['dpi'] = 120
    # does the plot need to be saved (1) or displayed (0)
    options['save'] = 0
    # Whether to plot the V vs RC plot.
    options['plot'] = 1
    # Whether to plot the PES graph
    options['graph_plot'] = 1
    # booleans tell if the ts energy values should be written
    options['write_ts_values'] = 1
    # booleans tell if the well and bimolecular energy values should be written
    options['write_well_values'] = 1
    # color of the energy values for the bimolecular products
    options['bimol_color'] = 'red'
    # color of the energy values of the wells
    options['well_color'] = 'blue'
    # color or the energy of the ts, put to 'none' to use same color as line
    options['ts_color'] = 'none'
    # boolean tells whether the molecule images should be shown on the graph
    options['show_images'] = 1
    # boolean that specifies which code was used for the 2D depiction
    options['rdkit4depict'] = 1
    # font size of the axes
    options['axes_size'] = 10
    # font size of the energy values
    options['text_size'] = 10
    # use linear lines instead of a polynomial
    options['linear_lines'] = 0
    # image interpolation
    options['interpolation'] = 'hanning'
    # graphs edge color, if set to 'energy', will be colored by that
    options['graph_edge_color'] = 'black'
    # graphs edge color, if set to 'energy', will be colored by that
    options['graph_bimolec_color'] = 'blue'
    # enable/disable generation of 2D depictions for resonant structures.
    options['reso_2d'] = 0
    # print report on paths connecting two species. Replace 0 with the two species names if to be activated.
    options['path_report'] = []
    # depth of search
    options['search_cutoff'] = 10
    # Scale graph nodes according to their stability.
    options['node_size_diff'] = 0

    if 'options' in input_dict:
        for line in input_dict['options']:
            if line.startswith('title'):
                options['title'] = int(line.split()[1])
            elif line.startswith('units'):
                options['units'] = line.split()[1]
            elif line.startswith('display_units'):
                options['display_units'] = line.split()[1]
            elif line.startswith('rounding'):
                options['rounding'] = int(line.split()[1])
            elif line.startswith('energy_shift'):
                options['energy_shift'] = float(line.split()[1])
            elif line.startswith('use_xyz'):
                options['use_xyz'] = int(line.split()[1])
            elif line.startswith('rescale'):
                options['rescale'] = line.split()[1]
            elif line.startswith('fh'):
                options['fh'] = float(line.split()[1])
            elif line.startswith('fw'):
                options['fw'] = float(line.split()[1])
            elif line.startswith('fs'):
                options['fs'] = float(line.split()[1])
            elif line.startswith('lw'):
                options['lw'] = float(line.split()[1])
            elif line.startswith('margin'):
                options['margin'] = float(line.split()[1])
            elif line.startswith('dpi'):
                options['dpi'] = int(line.split()[1])
            elif line.startswith('save'):
                if not options['save_from_command_line']:
                    options['save'] = int(line.split()[1])
            elif line.startswith('plot'):
                options['plot'] = int(line.split()[1])
            elif line.startswith('graph_plot'):
                options['graph_plot'] = bool(int(line.split()[1]))
            elif line.startswith('write_ts_values'):
                options['write_ts_values'] = int(line.split()[1])
            elif line.startswith('write_well_values'):
                options['write_well_values'] = int(line.split()[1])
            elif line.startswith('bimol_color'):
                options['bimol_color'] = line.split()[1]
            elif line.startswith('well_color'):
                options['well_color'] = line.split()[1]
            elif line.startswith('ts_color'):
                options['ts_color'] = line.split()[1]
            elif line.startswith('show_images'):
                options['show_images'] = int(line.split()[1])
            elif line.startswith('rdkit4depict'):
                options['rdkit4depict'] = int(line.split()[1])
            elif line.startswith('axes_size'):
                options['axes_size'] = float(line.split()[1])
            elif line.startswith('text_size'):
                options['text_size'] = float(line.split()[1])
            elif line.startswith('linear_lines'):
                options['linear_lines'] = int(line.split()[1])
            elif line.startswith('graph_edge_color'):
                options['graph_edge_color'] = str(line.split()[1])
            elif line.startswith('graph_bimolec_color'):
                options['graph_bimolec_color'] = line.split()[1]
            elif line.startswith('reso_2d'):
                options['reso_2d'] = int(line.split()[1])
            elif line.startswith('path_report'):
                options['path_report'].append(tuple(str(i) 
                                                    for i in line.split()[1:]))
            elif line.startswith('search_cutoff'):
                options['search_cutoff'] = int(line.split()[1])
            elif line.startswith('node_size_diff'):
                options['node_size_diff'] = float(line.split()[1])
            elif line.startswith('#'):
                # comment line, don't do anything
                continue
            else:
                if len(line) > 0:
                    print('Cannot recognize input line:')
                    print(line)
                # end if
            # end if
        # end for
    else:
        print('Warning, the input file arcitecture has changed,' +
              'use an "options" input tag to put all the options')
    # end if

    if options['display_units'] is None:
        options['display_units'] = options['units']

    for w in input_dict['wells']:
        w = w.split()
        name = w[0]
        energy = float(w[1])
        energy2 = None
        smi = None
        if len(w) > 2 and w[2] != '#':
            try:
                energy2 = float(w[2])
            except ValueError:
                smi = w[2]
        # end if
        w = well(name, energy, smi, energy2)
        wells.append(w)
    # end for
    for b in input_dict['bimolec']:
        b = b.split()
        name = b[0]
        energy = float(b[1])
        energy2 = None
        smi = []
        if len(b) > 2 and b[2] != '#':
            try:
                energy2 = float(b[2])
            except ValueError:
                smi = b[2:]
        b = bimolec(name, energy, smi, energy2)
        bimolecs.append(b)
    # end for

    # it is important that the wells and bimolecular products
    # are read prior to the transition states and barrierless
    # reactions because they need to be added as reactants
    # and products
    for t in input_dict['ts']:
        t = t.split()
        name = t[0]
        energy = eval(t[1])
        names = [t[2], t[3]]
        col = 'black'
        if options['graph_edge_color'] != 'energy':
            col = options['graph_edge_color']
        if len(t) > 4:
            col = t[4]
        t = ts(name, names, energy, col=col)
        tss.append(t)
    # end for
    for b in input_dict['barrierless']:
        b = b.split()
        name = b[0]
        names = [b[1], b[2]]
        col = 'gray'
        if len(b) > 3:
            col = b[3]
        b = barrierless(name, names, col=col)
        barrierlesss.append(b)
    # end for
    return input_str
# end def


def get_sd_prop(all_lines):
    """
    The get_sd_prop method interprets the input file
    which has a structure comparable to sdf molecular files
    """
    # split all the lines according to the keywords
    inputs = all_lines.split('> <')
    # do not consider the first keyword, this contains the comments
    inputs = inputs[1:]
    ret = {}
    for inp in inputs:
        inp = inp .split('>', 1)
        kw = inp[0].lower()  # keyword in lower case
        val = inp[1].strip().split('\n')  # values of the keywords
        val = [vi.strip() for vi in val if not vi.startswith('#')]
        val = [vi for vi in val if vi]
        ret[kw] = val
    # end for
    return ret
# end def


def position():
    """
    This method find initial position for all the wells, products and
    transition states. Initially, the wells are put in the middle and
    the products are divided on both sides. The transition states are
    positioned inbetween the reactants and products of the reaction
    """
    # do the rescaling, i.e. find the y values
    # y0 is the energy of the species to which rescaling is done
    y0 = next((w.energy for w in wells if w.name == options['rescale']), 0.)
    if y0 == 0.:
        list = (b.energy for b in bimolecs if b.name == options['rescale'])
        y0 = next(list, 0.)
    for w in wells:
        if np.isnan(w.energy):
            w.y = w.energy2 - y0
        else:
            w.y = w.energy - y0
    # end for
    for b in bimolecs:
        if np.isnan(b.energy):
            b.y = b.energy2 - y0
        else:
            b.y = b.energy - y0
    # end for
    for t in tss:
        t.y = t.energy - y0
    # end for

    # find the x values
    file_name = '{id}_xval.txt'.format(id=options['id'])
    if os.path.exists(file_name):
        # if file exists, read the x values
        fi = open(file_name, 'r')
        a = fi.read()
        fi.close()
        a = a.split('\n')
        for entry in a:
            if len(entry) > 0:
                name = entry.split(' ')[0]
                xval = eval(entry.split(' ')[1])
                for w in wells:
                    if w.name == name:
                        w.x = xval
                # end for
                for b in bimolecs:
                    if b.name == name:
                        b.x = xval
                # end for
                for t in tss:
                    if t.name == name:
                        t.x = xval
                # end for
            # end if
        # end for
    else:
        n = len(bimolecs)  # n is the number of bimolecular products
        n2 = n // 2

        for i, w in enumerate(wells):
            # wells are put in the middle
            w.x = n2 + 1 + i + 0.

        # end for
        for i, b in enumerate(bimolecs):
            # bimolecular products are put on both sides
            if i < n2:
                b.x = 1 + i + 0.
            else:
                b.x = len(wells) + 1 + i + 0.
                b.right = True
            # end if

        # end for
        for i, t in enumerate(tss):
            # transition states are put inbetween the reactant and proudct
            x1 = t.reactant.x
            x2 = t.product.x
            x3 = (x1+x2)/2.
            t.x = x3
        # end for
        save_x_values()  # write the x values to a file
    # end if
# end def


def generate_lines():
    """
    The method loops over the transition states and barrierless reactions
    and creates lines accordingly depending on the x and y coordinates
    """
    for t in tss:
        line1 = line(t.x,
                     t.y,
                     t.reactant.x,
                     t.reactant.y,
                     [t, t.reactant],
                     col=t.color)
        line2 = line(t.x,
                     t.y,
                     t.product.x,
                     t.product.y,
                     [t, t.product],
                     col=t.color)
        t.lines.append(line1)
        t.lines.append(line2)
    # end for
    for b in barrierlesss:
        b.line = line(b.reactant.x,
                      b.reactant.y,
                      b.product.x,
                      b.product.y,
                      [b, b.reactant, b.product],
                      col=b.color)
    # end for
# end def


def get_sizes():
    """
    Get the axis lengths and the sizes of the images
    """
    global xlow, xhigh, xmargin, ylow, yhigh, ymargin
    # TODO: what if wells or bimoleculs is empty,
    # min and max functions will give error
    if len(bimolecs) > 0:
        # x coords of tss are always inbetween stable species
        xlow = min(min([w.x for w in wells]), min([b.x for b in bimolecs]))
        xhigh = max(max([w.x for w in wells]), max([b.x for b in bimolecs]))
        # all tss lie above the lowest well
        ylow = min(min([w.y for w in wells]), min([b.y for b in bimolecs]))
    else:
        # x coords of tss are always inbetween stable species
        xlow = min([w.x for w in wells])
        xhigh = max([w.x for w in wells])
        # all tss lie above the lowest well
        ylow = min([w.y for w in wells])
    xmargin = options['margin']*(xhigh-xlow)
    try:
        yhigh = max([t.y for t in tss])
    except ValueError:
        yhigh = max([b.y for b in bimolecs])
    yhigh = max(yhigh, max([w.y for w in wells]))
    if len(bimolecs) > 0:
        yhigh = max(yhigh, max([b.y for b in bimolecs]))
    ymargin = options['margin']*(yhigh-ylow)
# end def


def plot():
    """Plotter method takes all the lines and plots them in one graph"""
    global xlow, xhigh, xmargin, ylow, yhigh, ymargin, xlen

    def showimage(s):
        """
        Get the extent and show the image on the plot
        s is the structure of which the image is to be
        plotted
        """
        global ymargin
        fn = '{id}_2d/{name}_2d.png'.format(id=options['id'], name=s.name)
        if os.path.exists(fn):
            img = mpimg.imread(fn)
            extent = None
            if s.name in extsd:
                extent = extsd[s.name]
            else:
                if not options['rdkit4depict']:
                    options['dpi'] = 120
                # end if

                imy = len(img) * options['fs']
                imx = len(img[0]) * options['fs']
                imw = (xhigh-xlow+0.)/(options['fw']+0.)*imx/options['dpi']
                imh = (yhigh-ylow+0.)/(options['fh']+0.)*imy/options['dpi']

                if isinstance(s, bimolec):
                    if s.x > (xhigh-xlow)/2.:
                        extent = (s.x + xmargin/5.,
                                  s.x + xmargin/5. + imw,
                                  s.y-imh/2.,
                                  s.y+imh/2.)
                    else:
                        extent = (s.x - xmargin/5. - imw,
                                  s.x - xmargin/5.,
                                  s.y-imh/2.,
                                  s.y+imh/2.)
                    # end if
                else:
                    extent = (s.x - imw/2.,
                              s.x + imw/2.,
                              s.y-ymargin/5. - imh,
                              s.y-ymargin/5.)
                # end if
            # end if
            im = ax.imshow(img, aspect='auto', extent=extent, zorder=-1, interpolation=options['interpolation'])
            # add to dictionary with the well it belongs to as key
            imgsd[s] = im
        # end if
    
    lines = []
    for t in tss:
        lines.append(t.lines[0])
        lines.append(t.lines[1])
    # end for
    for b in barrierlesss:
        lines.append(b.line)
    # end for
    get_sizes()
    plt.rcParams["figure.figsize"] = [options['fw'], options['fh']]
    plt.rcParams["font.size"] = options['axes_size']
    plt.rcParams['figure.dpi'] = options['dpi']

    matplotlib.rc("figure", facecolor="white")
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.set_xticklabels([])

    if options['show_images']:
        for w in wells:
            showimage(w)
        # end for
        for b in bimolecs:
            showimage(b)
        # end for
    save_im_extent()  # save the positions of the images to a file
    
    # draw the lines
    # in case of linear lines, calculate the distance of the horizontal pieces
    if options['linear_lines']:
        xlen = (len(wells) + len(bimolecs)) / (4 * (xhigh - xlow))
    
    for line in lines:
        lw = options['lw']
        alpha = 1.0
        ls = 'solid'
        if line.color == 'dotted':
            ls = 'dotted'
            line.color = 'gray'
#        elif line.color == 'blue' or line.color == 'b':
#            ls = 'dashed'
        if line.straight_line:
            if line.xmin == line.xmax:  # plot a vertical line
                ymin = min(line.y1, line.y2)
                ymax = max(line.y1, line.y2)
                a = ax.vlines(x=line.xmin,
                              ymin=ymin,
                              ymax=ymax,
                              color=line.color,
                              ls=ls,
                              linewidth=lw,
                              alpha=alpha)
                for struct in line.chemstruct:
                    # add to the lines dictionary
                    linesd[struct].append(a)
                # end for
            else:  # plot a horizontal line
                a = ax.hlines(y=line.y1,
                              xmin=line.xmin,
                              xmax=line.xmax,
                              color=line.color,
                              linewidth=lw,
                              alpha=alpha)
                for struct in line.chemstruct:
                    # add to the lines dictionary
                    linesd[struct].append(a)
                # end for
            # end if
        else:
            if options['linear_lines']:
                xlist = [line.xmin, line.xmin + xlen / 2,
                         line.xmax - xlen / 2, line.xmax]
                y = [line.y1, line.y1, line.y2, line.y2]
            else:
                xlist = np.arange(line.xmin,
                                  line.xmax,
                                  (line.xmax-line.xmin) / 1000)
                a = line.coeffs
                y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
            pl = ax.plot(xlist,
                         y,
                         color=line.color,
                         ls=ls,
                         linewidth=lw,
                         alpha=alpha)
            for struct in line.chemstruct:
                # add to the lines dictionary
                linesd[struct].append(pl)
            # end for
        # end if
    # end for
    ax.set_xlim([xlow-xmargin, xhigh+xmargin])
    ax.set_ylim([ylow-ymargin, yhigh+ymargin])

    # write the name and energies to the plot
    for w in wells:
        if options['write_well_values']:
            t = ax.text(w.x,
                        w.y-ymargin/10,
                        '{:.{}f}'.format(w.y, options['rounding']),
                        fontdict={'size': options['text_size']},
                        ha='center', va='top',
                        color=options['well_color'],
                        picker=True)
            # add to dictionary with the well it belongs to as key
            textd[w] = t
        # end if
    # end for
    for b in bimolecs:
        if options['write_well_values']:
            # write the text values below the line:
            t = ax.text(b.x,
                        b.y-ymargin/10,
                        '{:.{}f}'.format(b.y, options['rounding']),
                        fontdict={'size': options['text_size']},
                        ha='center',
                        va='top',
                        color=options['bimol_color'],
                        picker=True)
            # add to dictionary with the well it belongs to as key
            textd[b] = t
        # end if
    # end for
    for t in tss:
        if options['write_ts_values']:
            color = t.color
            if options['ts_color'] != 'none':
                color = options['ts_color']
            te = ax.text(t.x,
                         t.y+ymargin/30,
                         '{:.{}f}'.format(t.y, options['rounding']),
                         fontdict={'size': options['text_size']},
                         ha='center',
                         va='bottom',
                         color=color,
                         picker=True)
            # add to dictionary with the well it belongs to as key
            textd[t] = te
        # end if
    # end for

    sel = selecthandler()
    dr = dragimage()

    if options['title']:
        plt.title('Potential energy surface of {id}'.format(id=options['id']))
    # defualts to units if not specified
    plt.ylabel('Energy ({display_units})'.format(display_units=options['display_units']))
    if options['save']:
        plt.savefig(f'{options["id"]}_pes_plot.png', bbox_inches='tight')
    else:
        plt.show()
# end def


def generate_2d_depiction():
    """
    2D depiction is generated (if not yet available) and
    stored in the directory join(input_id, '_2d')
    This is only done for the wells and bimolecular products,
    2D of tss need to be supplied by the user
    """
    def get_smis(m, smis, files):
        # name and path of png file
        if len(smis) > 0:
            return smis
        elif files and all([os.path.exists(f) for f in files]):
            smis = []
            for f in files:
                try:
                    obmol = next(pybel.readfile('xyz', f))
                    smi = obmol.write("smi").split()[0]
                    smi = smi.replace('=[CH]=C', '=C[C]')
                    smi = smi.replace('N(=O)=O', 'N(=O)[O]')
                    smi = smi.replace('[NH]([CH2])(O)[O]',
                                      '[NH+]([CH2])(O)[O-]')
                    smis.append(smi)
                except NameError:
                    print('Could not generate smiles for {n}'.format(n=m.name))
                # end try
            # end for
        else:
            raise FileNotFoundError(f'No xyz file found for well {m.name} ('
                                    f'{m.smi})')
        return smis
    # end def

    def reaction_smi():
        """Given a list of transition states and barrierless channels
        write a reaction_smi.out file that contains the reactions as, e.g.:
        OOC[CH2] = [OH] + O1CC1i  0.0  10.72  -16.11
        The three numbers are the energy of the reactant, transition state and product.
        The barrier height for barrierless reactions is given as max(Ereact,Eprod)
        """

        with open('reaction_smi.out', 'w') as f:
            for t in tss:
                if not t.reactant.smi:
                    t.reactant.smi = get_smis(t.reactant, [], t.reactant.xyz_files)
                elif isinstance(t.reactant.smi, str):
                    t.reactant.smi = [t.reactant.smi]
                f.write(t.reactant.smi[0])
                for r in t.reactant.smi[1:]:
                    f.write(' + ')
                    f.write(r)
                f.write(' = ')
                if not t.product.smi:
                    t.product.smi = get_smis(t.product, [], t.product.xyz_files)
                elif isinstance(t.product.smi, str):
                    t.product.smi = [t.product.smi]
                f.write(t.product.smi[0])
                for p in t.product.smi[1:]:
                    f.write(' + ')
                    f.write(p)
                f.write('  {:.2f}  {:.2f}   {:.2f}'.format(t.reactant.energy, t.energy, t.product.energy))
                f.write('\n')

            for b in barrierlesss:
                if not b.reactant.smi:
                    b.reactant.smi = get_smis(b.reactant, [], b.reactant.xyz_files)
                elif isinstance(b.reactant.smi, str):
                    b.reactant.smi = [b.reactant.smi]
                f.write(b.reactant.smi[0])
                for r in b.reactant.smi:
                    f.write(' + ')
                    f.write(r)
                f.write(' = ')
                if not b.product.smi:
                    b.product.smi = get_smis(b.product, [], b.product.xyz_files)
                elif isinstance(b.product.smi, str):
                    b.product.smi = [b.product.smi]
                f.write(b.product.smi[0])
                for p in b.product.smi[1:]:
                    f.write(' + ')
                    f.write(p)
                f.write('  {:.2f}  {:.2f}  {:.2f}'.format(b.reactant.energy, max(b.reactant.energy, b.product.energy), b.product.energy))
                f.write('\n')

    def generate_2d(m, smis):
        
        png_filename = '{id}_2d/{name}_2d{confid}.png'
        if len(smis) == 0:
            print('Could not generate 2d for {name}'.format(name=m.name))
            return
        smi = '.'.join(smis)
        if os.path.isfile(png_filename.format(id=options['id'], name=m.name,
                                              confid='')):
            return
        try:
            if options['reso_2d']:
                try:
                    reson_mols = gen_reso_structs(smi, min_rads=True)
                except RuntimeError:
                    print(f'Warning: Unable to generate resonant structure for '
                          f'{smi}.')
                    options['reso_2d'] = 0
            else:
                mol = Chem.MolFromSmiles(smi, sanitize=False)
                mol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                                 | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                                 | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                                 | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                                 | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                                 | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                                 catchErrors=True)
                reson_mols = [mol]

            resol = 5
            size_x = 100 * resol

            # NEW METHOD
            # opts = Draw.rdMolDraw2D.MolDrawOptions()  # New way
            # opts.minFontSize = 30 * resol
            # opts.maxFontSize = 9 * resol
            # opts.bondLineWidth = 1 * resol
            # opts.padding = 0.15
            # opts.noAtomLabels = True

            # OLD METHOD (cannot use Draw.MolToImage or Draw.MolToFile):
            opts = Draw.DrawingOptions()
            opts.dotsPerAngstrom = 15 * resol
            opts.atomLabelFontSize = 9 * resol
            opts.bondLineWidth = 1 * resol
            opts.radicalSymbol = '•'
            for i, mol in enumerate(reson_mols):
                AllChem.Compute2DCoords(mol)
                cc = mol.GetConformer()
                xx = []
                yy = []
                for j in range(cc.GetNumAtoms()):
                    pos = cc.GetAtomPosition(j)
                    xx.append(pos.x)
                    yy.append(pos.y)
                sc = 50
                dx = round(max(xx) - min(xx)) * sc
                dy = round(max(yy) - min(yy)) * sc
                size_x = round(size_x * (1 + (max(dx, dy) - 200) / 500))
                size = (size_x,) * 2
                # NEW METHOD
                # img = Draw.MolToImage(mol, kekulize=False, wedgeBonds=False,
                #                       options=opts, size=size)

                # OLD METHOD (cannot use Draw.MolToImage or Draw.MolToFile):
                img = Image.new("RGBA", tuple(size))
                canvas = Canvas(img)
                drawer = Draw.MolDrawing(canvas=canvas, drawingOptions=opts)
                drawer.AddMol(mol)
                canvas.flush()

                # Convert each white pixel to transparent.
                pixels = img.getdata()
                new_pixels = []
                for pix in pixels:
                    if pix == (255, 255, 255, 255):
                        new_pixels.append((255, 255, 255, 0))
                    else:
                        new_pixels.append(pix)
                img.putdata(new_pixels)

                if i == 0:
                    img.save(png_filename.format(id=options['id'], name=m.name,
                                                 confid=''))
                else:
                    img.save(png_filename.format(id=options['id'], name=m.name,
                                                 confid=f'_{i}'))
        except (NameError, RuntimeError):
            try:
                options['rdkit4depict'] = 0
                obmol = pybel.readstring("smi", smi)
                obmol.draw(show=False, filename=png_filename.format(id=options['id'],
                                                                    name=m.name,
                                                                    confid=''))
                img = Image.open(png_filename.format(id=options['id'],
                                                     name=m.name,
                                                     confid=''))
                new_size = (280, 280)
                im_new = Image.new("RGB", new_size, 'white')
                im_new.paste(img, ((new_size[0] - img.size[0]) // 2,
                                   (new_size[1] - img.size[1]) // 2))
                im_new.save(png_filename.format(id=options['id'], name=m.name,
                                                confid=''))
            except NameError:
                print('Could not generate 2d for {n}'.format(n=m.name))
                return

    # end def
    # make the directory with the 2d depictions, if not yet available
    dir = options['id'] + '_2d'
    try:
        os.stat(dir)
    except FileNotFoundError:
        os.mkdir(dir)
    for w in wells:
        smis = []
        if w.smi is not None:
            smis.append(w.smi)
        else:
            smis = get_smis(w, smis, w.xyz_files)
        generate_2d(w, smis)
    # end for
    for b in bimolecs:
        generate_2d(b, get_smis(b, b.smi, b.xyz_files))
    # end for

    reaction_smi()  # write file with reaction string

# end def


def updateplot(struct, x_change):
    """
    Update the plot after the drag event of a stationary point (struct),
    move all related objects by x_change in the x direction,
    and regenerate the corresponding lines
    """
    global xlow, xhigh, xmargin, ylow, yhigh, ymargin, xlen
    xlow_old = xlow
    xhigh_old = xhigh
    # set the new sizes of the figure
    get_sizes()
    plt.gca().set_xlim([xlow-xmargin, xhigh+xmargin])
    plt.gca().set_ylim([ylow-ymargin, yhigh+ymargin])
    ratio = (xhigh - xlow) / (xhigh_old - xlow_old)
    # generate new coordinates for the images
    if struct in imgsd:
        old_extent = imgsd[struct].get_extent()
        extent_change = (x_change, x_change, 0, 0)
        extent = [old_extent[i] + extent_change[i] for i in range(0, 4)]
        imgsd[struct].set_extent(extent=extent)
    # need to scale all other images as well, but only doing it if 
    # the overall width has changed
    if ratio != 1:
        for key in imgsd:
            oe = imgsd[key].get_extent()  # old extent
            extent = [oe[0], oe[0] + (oe[1] - oe[0]) * ratio, oe[2], oe[3]]
            imgsd[key].set_extent(extent=extent)
    # generate new coordinates for the text
    if struct in textd:
        old_pos = textd[struct].get_position()
        new_pos = (old_pos[0]+x_change, old_pos[1])
        textd[struct].set_position(new_pos)
    # generate new coordinates for the lines
    for t in tss:
        if (struct == t or struct == t.reactant or struct == t.product):
            t.lines[0] = line(t.x,
                              t.y,
                              t.reactant.x,
                              t.reactant.y,
                              [t, t.reactant],
                              col=t.color)
            t.lines[1] = line(t.x,
                              t.y,
                              t.product.x,
                              t.product.y,
                              [t, t.product],
                              col=t.color)
            for i in range(0, 2):
                li = t.lines[i]
                if li.straight_line:
                    print('straight line')
                else:
                    if options['linear_lines']:
                        xlist = [li.xmin, li.xmin + xlen / 2,
                                 li.xmax - xlen / 2, li.xmax]
                        y = [li.y1, li.y1, li.y2, li.y2]
                    else:
                        xlist = np.arange(li.xmin,
                                          li.xmax,
                                          (li.xmax-li.xmin) / 1000)
                        a = li.coeffs
                        y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
                    linesd[t][i][0].set_xdata(xlist)
                    linesd[t][i][0].set_ydata(y)
                # end if
            # end for
        # end if
    # end for
    for b in barrierlesss:
        if (struct == b.reactant or struct == b.product):
            b.line = line(b.reactant.x,
                          b.reactant.y,
                          b.product.x,
                          b.product.y,
                          [b.reactant, b.product],
                          col=b.color)
            li = b.line
            if li.straight_line:
                print('straight line')
            else:
                if options['linear_lines']:
                    xlist = [li.xmin, li.xmin + xlen / 2,
                             li.xmax - xlen / 2, li.xmax]
                    y = [li.y1, li.y1, li.y2, li.y2]
                else:
                    xlist = np.arange(li.xmin, li.xmax, (li.xmax-li.xmin) / 1000)
                    a = li.coeffs
                    y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
                linesd[b][0][0].set_xdata(xlist)
                linesd[b][0][0].set_ydata(y)
            # end if
        # end if
    # end for
    plt.draw()
# end def


def highlight_structure(struct=None):
    """
    for all the lines, text and structures that are not
    directly connected to struct, set alpha to 0.15
    """
    if struct is None:
        alpha = 1.
    else:
        alpha = 0.15
    # end if
    # get all the tss and barrierlesss with this struct
    highlight = []
    lines = []
    for t in tss:
        if t == struct or t.reactant == struct or t.product == struct:
            highlight.append(t)
            if t.reactant not in highlight:
                highlight.append(t.reactant)
            if t.product not in highlight:
                highlight.append(t.product)
            lines = lines + linesd[t]
        # end if
    # end for
    for b in barrierlesss:
        if b.reactant == struct or b.product == struct:
            highlight.append(b)
            if b.reactant not in highlight:
                highlight.append(b.reactant)
            if b.product not in highlight:
                highlight.append(b.product)
            lines = lines + linesd[b]
        # end if
    # end for
    for struct in linesd:
        for li in linesd[struct]:
            if li in lines:
                li[0].set_alpha(1.)
            else:
                li[0].set_alpha(alpha)
            # end if
        # end for
    # end for
    for struct in textd:
        if struct in highlight:
            textd[struct].set_alpha(1.)
        else:
            textd[struct].set_alpha(alpha)
        # end if
    # end for
    for struct in imgsd:
        if struct in highlight:
            imgsd[struct].set_alpha(1.)
        else:
            imgsd[struct].set_alpha(alpha)
        # end if
    # end for
    plt.draw()
# end def


def save_x_values():
    """
    save the x values of the stationary points to an external file
    """
    fi = open('{id}_xval.txt'.format(id=options['id']), 'w')
    if len(wells) > 0:
        lines = ['{n} {v:.2f}'.format(n=w.name, v=w.x) for w in wells]
        fi.write('\n'.join(lines)+'\n')
    if len(bimolecs) > 0:
        lines = ['{n} {v:.2f}'.format(n=b.name, v=b.x) for b in bimolecs]
        fi.write('\n'.join(lines)+'\n')
    if len(tss) > 0:
        lines = ['{n} {v:.2f}'.format(n=t.name, v=t.x) for t in tss]
        fi.write('\n'.join(lines))
    fi.close()
# end def


def save_im_extent():
    """Save the x values of the stationary points to an external file"""
    fi = open(f'{options["id"]}_im_extent.txt', 'w')
    for key in imgsd:
        e = imgsd[key].get_extent()
        vals = '{:.2f} {:.2f} {:.2f} {:.2f}'.format(e[0], e[1], e[2], e[3])
        fi.write('{name} {vals}\n'.format(name=key.name, vals=vals))
    fi.close()
# end def


def read_im_extent():
    """Read the extents of the images if they are present in a file_name."""
    fname = f'{options["id"]}_im_extent.txt'
    if os.path.exists(fname):
        fi = open(fname, 'r')
        a = fi.read()
        fi.close()
        a = a.split('\n')
        for entry in a:
            pieces = entry.split(' ')
            if len(pieces) == 5:
                extsd[pieces[0]] = [eval(pieces[i]) for i in range(1, 5)]
            # end if
        # end for
    # end if
# end def

def convert_units(energy):
    """Converts energy from 'units' to 'display_units' and apply 'energy_shift'"""
    # apply shift in original units
    energy += options['energy_shift']
    # convert to kJ/mol first
    if options['units'] == 'kcal/mol':
        energy = energy * 4.184
    elif options['units'] == 'eV':
        energy = energy * 96.4852912
    elif options['units'] == 'Hartree':
        energy = energy * 2625.498413
    # convert to display_units
    if options['display_units'] == 'kcal/mol':
        energy = energy / 4.184
    elif options['display_units'] == 'eV':
        energy = energy / 96.4852912
    elif options['display_units'] == 'Hartree':
        energy = energy / 2625.498413
    return energy


def create_interactive_graph(meps):
    """Create an interactive graph with pyvis."""
    
    g = net.Network(height='1000px', width='90%', heading='', directed=True)

    base_energy = next((species.energy for species in wells + bimolecs 
                        if species.name == options['rescale']), 0)
    for well in wells:
        rel_energy = round(well.energy - base_energy, options['rounding'])
        if np.isnan(rel_energy):
            label = '❌'
        else:
            label = str(rel_energy)
        if not well.energy2 is None:
            rel_energy2 = round(well.energy2 - base_energy, options["rounding"])
            if np.isnan(rel_energy2):
                label += ' ❌'
            else:
                label += f' {rel_energy2}'
        g.add_node(well.name, label=label, borderWidth=3, title=f'{well.name}', 
                   shape='circularImage', image=f'{options["id"]}_2d/{well.name}_2d.png', 
                   size=80, font='30',
                   color={'background': '#FFFFFF', 'border': 'black',
                          'highlight': {'border': '#FF00FF', 'background': '#FFFFFF'}})
    for bim in bimolecs:
        border_color = options['graph_bimolec_color']
        rel_energy = round(bim.energy - base_energy, options['rounding'])
        if np.isnan(rel_energy):
            label = '❌'
        else:
            label = str(rel_energy)
        if not bim.energy2 is None:
            rel_energy2 = round(bim.energy2 - base_energy, options["rounding"])
            if np.isnan(rel_energy2):
                label += ' ❌'
            else:
                label += f' {rel_energy2}'
        g.add_node(bim.name, label=label, borderWidth=3, title=f'{bim.name}', 
                   shape='circularImage', image=f'{options["id"]}_2d/{bim.name}_2d.png', 
                   size=80, font='30', 
                   color={'background': '#FFFFFF', 'border': border_color,
                          'highlight': {'border': '#FF00FF', 
                          'background': '#FFFFFF'}})

    min_ts_energy = min([ts.energy for ts in tss])
    max_ts_energy = max([ts.energy for ts in tss])
    ts_energy_range = max_ts_energy - min_ts_energy
    cmap = plt.get_cmap('viridis')

    for ts in tss:
        norm_energy = (ts.energy - min_ts_energy) / ts_energy_range
        if ts in [mep['bottle_neck'] for mep in meps]:
            color = '#E9C46A'
        elif ts in [rxn for mep in meps for rxn in mep['rxns']]:
            color = '#2A9D8F'
        elif options['graph_edge_color'] == 'energy':
            red, green, blue = np.array(cmap.colors[int(norm_energy * 255)]) * 255 
            color = f'rgb({red},{green},{blue})'
        else:  
            color = ts.color
        rel_energy = round(ts.energy - base_energy, options["rounding"])
        g.add_edge(ts.reactant.name, ts.product.name, 
                   title=f'{rel_energy} {options["display_units"]}',
                   color={'highlight': '#FF00FF', 'color': color}, 
                   width=(1 - norm_energy) * 20 + 1, arrows='')

    for bless in barrierlesss:
        norm_energy = (bless.product.energy - min_ts_energy) / ts_energy_range
        if options['graph_edge_color'] == 'energy':
            red, green, blue = np.array(cmap.colors[int(norm_energy * 255)]) * 255 
            color = f'rgb({red},{green},{blue})'
        else:  
            color = bless.color
        rel_energy = round(bless.product.energy - base_energy, options["rounding"])
        g.add_edge(bless.reactant.name, bless.product.name, 
                   title=f'{rel_energy} {options["display_units"]}', 
                   color={"highlight": "#FF00FF", 'color': color},
                   width=(1 - norm_energy) * 20 + 1, arrows='')

    g.set_edge_smooth('dynamic')
    g.show_buttons(filter_=['physics'])
    g.save_graph(f'{options["id"]}.html')

    return 0


def is_path_valid(path):
    return all(['_' not in name for name in path[1:-1]])


def write_section(f, input_lines, stopsign, start, path):
    for ll, line in enumerate(input_lines[start:]):
        if not line.startswith(stopsign):
            if stopsign == '> <wells>':
                if line.startswith('> <id>'):
                    f.write(f'> <id> aux_{path[0]}_{path[-1]}\n')
                elif line.startswith('plot'):
                    f.write(f'plot 1\n')
                elif line.startswith('path_report') or line.startswith('search_cutoff'):
                    continue
                else:
                    f.write(f'{line}\n')
            elif stopsign == '> <bimolec>' or stopsign == '> <ts>':
                if line.split()[0] in path:
                    f.write(f'{line}\n')
                elif line.startswith('>'):
                    f.write(f'{line}\n')
            elif stopsign == '> <barrierless>':
                if len(line.split()) < 4:
                    f.write(f'{line}\n')
                elif line.split()[2] in path and line.split()[3] in path:
                    if abs(path.index(line.split()[2]) - path.index(line.split()[3])) != 1:
                        continue
                    else:
                        f.write(f'{line}\n')
            elif stopsign == '> <help>':
                if len(line.split()) < 3: 
                    f.write(f'{line}\n')
                elif line.split()[1] in path and line.split()[2] in path:
                    if abs(path.index(line.split()[1]) - path.index(line.split()[2])) != 1:
                        continue
                    else:
                        f.write(f'{line}\n')
            else:
                f.write(f'{line}\n')
        else:
            return ll+start


def gen_graph():
    """Generate a networkx graph object to work with

    Returns:
        networkx.Graph: A Graph object representation of the PES.
    """
    graph = nx.Graph()
    base_energy = next((species.energy for species in wells + bimolecs 
                        if species.name == options['rescale']), 0)
    
    for reac in tss + barrierlesss:
        rname = reac.reactant.name
        pname = reac.product.name
        graph.add_edge(rname, pname)
        graph[rname][pname]['energy'] = round(reac.energy - base_energy, 1)

    return graph


def find_mep(graph, user_input):
    """Find the minimum energy path between two species in a PES.

    Args:
        graph (networkx.Graph): A graph object containing the information of 
            the PES
        user_input (str): The input file passed to pesviewer.

    Returns:
        NoneType: None
    """
    meps = []
    for species_pair in options['path_report']:
        meps.append({})
        current_mep = meps[-1]
        spec_1 = species_pair[0]
        spec_2 = species_pair[1]
        max_length = options['search_cutoff']
        paths = nx.all_simple_paths(graph, spec_1, spec_2, cutoff=max_length)
        max_barr = np.inf
        for valid_path in [path for path in paths if is_path_valid(path)]:
            path_energies = [graph[valid_path[i]][valid_path[i+1]]['energy'] 
                             for i in range(len(valid_path)-1)]    
            if (max_barr == max(path_energies) and len(valid_path) < len(mep_species)) \
                    or max_barr > max(path_energies):
                max_barr = max(path_energies)  # the bottle neck
                mep_species = valid_path
                bottle_neck_idx = path_energies.index(max_barr)
        current_mep['species'] = mep_species
        current_mep['energies'] = path_energies
        current_mep['rxns'] = []
        for i, sp in enumerate(mep_species[:-1]):
            for rxn in tss + barrierlesss:
                rname = rxn.reactant.name
                pname = rxn.product.name
                if (rname == sp and pname == mep_species[i+1]) \
                        or (pname == sp and rname == mep_species[i+1]):
                    current_mep['rxns'].append(rxn)
        current_mep['bottle_neck'] = current_mep['rxns'][bottle_neck_idx]

        print(f'The bottle neck barrier between {spec_1} and {spec_2} with a ' \
              f'{options["search_cutoff"]} depth search is {max_barr} ' \
              f'{options["display_units"]} high.')

        # write new pesviewer input for just this path
        input_lines = user_input.split('\n')
        with open(f'{mep_species[0]}_{mep_species[-1]}.inp', 'w') as f:
            stop = write_section(f, input_lines, '> <wells>', 0, mep_species)
            stop = write_section(f, input_lines, '> <bimolec>', stop, mep_species)
            stop = write_section(f, input_lines, '> <ts>', stop, mep_species)
            stop = write_section(f, input_lines, '> <barrierless>', stop, mep_species)
            stop = write_section(f, input_lines, '> <help>', stop, mep_species)

    return meps


def main(fname=None):
    """Main method to run the PESViewer"""
    if fname is None and len(sys.argv) > 1:
        fname = sys.argv[1]
        options['save'] = 0
        options['save_from_command_line'] = 0
        if len(sys.argv) > 2 and sys.argv[2] == 'save':  # Save the plot.
            options['save'] = 1
            options['save_from_command_line'] = 1
    elif fname is None and len(sys.argv) == 1:
        print('To use the pesviewer, supply an input file as argument.')
        sys.exit(-1)
    user_input = read_input(fname)  # read the input file
    # initialize the dictionaries
    for w in wells:
        linesd[w] = []
    
    for b in bimolecs:
        linesd[b] = []
    
    for t in tss:
        linesd[t] = []
    
    for b in barrierlesss:
        linesd[b] = []
    
    read_im_extent()  # read the position of the images, if known
    position()  # find initial positions for all the species on the graph
    generate_lines()  # generate all the line
    # generate 2d depiction from the smiles or 3D structure,
    # store them in join(input_id, '_2d')
    generate_2d_depiction()
    graph = gen_graph()
    meps = find_mep(graph, user_input)
    if options['plot']:
        plot()
        if meps:
            print('To draw 2D plots for the individual MEPs type:')
            for mep in meps:
                print(f'\tpesviewer {mep["species"][0]}_{mep["species"][-1]}.inp')
    
    if options['graph_plot']:
        create_interactive_graph(meps)


def pesviewer(fname=None):
    options['save'] = 0
    options['save_from_command_line'] = 0
    main(fname)


if __name__ == "__main__":
    main()
