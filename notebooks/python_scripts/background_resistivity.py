from SimPEG import Mesh, Utils
import numpy as np
from simpegEM1D import diffusion_distance
from pymatsolver import Pardiso
from SimPEG import EM
from scipy.constants import mu_0
from scipy.interpolate import interp1d
from simpegEM1D.Waveforms import piecewise_pulse_fast
from pyMKL import mkl_set_num_threads
from simpegskytem import ProblemSkyTEM
from mesh_utils import refineTree, meshBuilder

num_threads = 4
mkl_set_num_threads(num_threads)

tmin, tmax = 1e-6, 1e-2
sigma_for_padding = 1./100.
# print (diffusion_distance(1e-2, sigma_for_padding) * 2)
# print (diffusion_distance(1e-6, sigma_for_padding) / 2.)
# print (diffusion_distance(1e-6, sigma_for_padding) * 2)
padding_distance = np.round(diffusion_distance(1e-2, sigma_for_padding) * 2)

sigma_halfspace = 1./20.

data_dir = "../"
waveform_hm = np.loadtxt(data_dir+"HM_butte_312.txt")
time_gates_hm = np.loadtxt(data_dir+"HM_butte_312_gates")[7:,:] * 1e-6
waveform_lm = np.loadtxt(data_dir+"LM_butte_312.txt")
time_gates_lm = np.loadtxt(data_dir+"LM_butte_312_gates")[8:,:] * 1e-6

time_input_currents_HM = waveform_hm[:,0]
input_currents_HM = waveform_hm[:,1]
time_input_currents_LM = waveform_lm[:,0]
input_currents_LM = waveform_lm[:,1]

time_LM = time_gates_lm[:,3] - waveform_lm[:,0].max()
time_HM = time_gates_hm[:,3] - waveform_hm[:,0].max()
source_area = 536.36

base_frequency_HM = 30.
base_frequency_LM = 210.


x = np.linspace(-250, 250)
y = np.linspace(-250, 250)
z = np.array([0.])
dem = Utils.ndgrid(x,y,z)


result_dir = "../background_resistivity/"

maxLevel = 11
h = [0.0625, 0.0625, 0.0625]
octreeLevel = [0, 0, 0, 0, 0, 1, 1, 4, 4, 10]
n_tower = 2
ground_resistance = 20.
layer_thickness = 4
resistivity_near = 20.
resistivity_backgrounds = [5, 20, 80, 320, 1280]
for i_res, resistivity_background in enumerate(resistivity_backgrounds):

    padDist = np.ones((3, 2)) * padding_distance
    mesh = meshBuilder(dem, h, padDist,
                                        meshType='TREE',
                                        verticalAlignment='center')
    # Refine the mesh around topographyttb
    mesh = refineTree(mesh, dem, dtype='surface',
                                       octreeLevels=octreeLevel, finalize=False)
    n_segment = int(n_tower)
    l_copper = 80.
    ys = np.arange(n_segment) * l_copper
    shift = -ys.max()/2.
    ys += shift

    x1 = 0.
    z1, z2 = -3., 10.
    ys_corr = []
    for y_temp in ys:
        ys_corr.append(mesh.vectorNy[np.argmin(abs(mesh.vectorNy-y_temp))])

    ys = np.hstack(ys_corr)
    x1 = mesh.vectorNx[np.argmin(abs(mesh.vectorNx-x1))]
    z1 = mesh.vectorNz[np.argmin(abs(mesh.vectorNz-z1))]
    z2 = mesh.vectorNz[np.argmin(abs(mesh.vectorNz-z2))]

    pts_top = []
    pts_bottom = []
    for ii in range(n_segment-1):
        ind_y_top = np.logical_and(mesh.vectorNy>=ys[ii], mesh.vectorNy<=ys[ii+1])
        ex = np.ones(ind_y_top.sum()) * x1
        ez = np.ones(ind_y_top.sum()) * z2
        pts_top.append(np.c_[ex, mesh.vectorNy[ind_y_top], ez])
        ez = np.ones(ind_y_top.sum()) * z1
        pts_bottom.append(np.c_[ex, mesh.vectorNy[ind_y_top], ez])

    pts_tower = []
    for ii in range(n_segment):
        ind_z_side = np.logical_and(mesh.vectorNz>=z1, mesh.vectorNz<=z2)
        ex = np.ones(ind_z_side.sum())*x1
        ey = np.ones(ind_z_side.sum())*ys[ii]
        pts_tower.append(np.c_[ex, ey, mesh.vectorNz[ind_z_side]])

    # pts = np.vstack((np.vstack(pts_top), np.vstack(pts_bottom), np.vstack(pts_tower)))
    pts = np.vstack((np.vstack(pts_top), np.vstack(pts_tower)))

    mesh = refineTree(mesh, pts, dtype='point',
                                       octreeLevels=[1, 0, 0], finalize=False)
    survey_length = 400.
    dx = 4
    n_src = survey_length / dx
    x = np.arange(n_src) * dx
    x -= x.max()/2.

    y = np.array([abs(mesh.vectorCCy).min()])
    z_src = 40.
    z = np.array([z_src])

    xyz = Utils.ndgrid(x, y, z)
    mesh = refineTree(mesh, xyz, dtype='point',
                                       octreeLevels=[1, 0, 0], finalize=True, maxLevel=maxLevel)
    # sigma_ground = 1e3
    sigma = np.ones(mesh.nC) * 1./resistivity_background
    inds_air = mesh.gridCC[:, 2] > 0.
    sigma[mesh.gridCC[:, 2] > 0.] = 1e-8

    indArr, levels = mesh.__getstate__()
    inds = levels == levels.max()

    temp = np.unique(levels)
    cell_size = 2**(temp.max()-temp)
    print (n_tower)

    radius_copper = 3.264 * 1e-3/2.
    area_copper = np.pi * radius_copper **2
    radius_rod = 15.87 * 1e-3 / 2.
    area_rod = np.pi * radius_copper **2
    sigma_copper = 6e7
    sigma_rod = 1e8
    area = (mesh.hx.min() * 4)**2
    sigma[np.logical_and(inds, inds_air)] = sigma_copper * area_copper / area
    inds_layer_near = np.logical_and(mesh.gridCC[:,2]<0., mesh.gridCC[:,2]>-layer_thickness)
    sigma[inds_layer_near] = 1./resistivity_near
    inds_layer = np.logical_and(mesh.gridCC[:,2]<0., mesh.gridCC[:,2]>-3)
    sigma[np.logical_and(inds, ~inds_air) & (inds_layer)] = sigma_rod * area_rod / area

    from pymatsolver import Pardiso
    from SimPEG import EM
    from scipy.constants import mu_0

    def compute_response(sigma):
        srcList = []
        z_src = 40.
        z_offset = 0.
        x_offset = 0.
        radius = 13.25
        for x_src in x:
            z_src = z_src
            rxloc = np.array([x_src+x_offset, 0., z_src+z_offset])
            srcloc = np.array([x_src, 0., z_src])
            rx = EM.TDEM.Rx.Point_dbdt(rxloc, np.logspace(np.log10(1e-5), np.log10(1e-2), 31), 'z')
            src = EM.TDEM.Src.CircularLoop([rx], waveform=EM.TDEM.Src.StepOffWaveform(), loc=srcloc, radius=radius)
            srcList.append(src)
        survey = EM.TDEM.Survey(srcList)
        prb = ProblemSkyTEM(mesh, verbose=False, sigma=sigma)
        prb.timeSteps = [(3e-7, 6),(1e-6, 5),(2e-6, 5),(5e-6, 5),(1e-5, 5),(2e-5, 5),(5e-5, 5),(1e-4, 5),(2e-4, 5),(5e-4, 5),(1e-3, 8),(2e-3, 8),(5e-3, 5),(1e-2, 4)]
        prb.Solver = Pardiso
        prb.pair(survey)
        data = prb.simulate(
            [],
            time_HM,
            time_LM,
            time_input_currents_HM,
            input_currents_HM,
            time_input_currents_LM,
            input_currents_LM,
        )
        return xyz, data

    xyz, data = compute_response(sigma)
    print(i_res, resistivity_background)
    np.save(result_dir+'xyz', xyz)
    np.save(result_dir+'data' + str(i_res), data)
