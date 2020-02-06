from SimPEG.Utils.matutils import mkvc, ndgrid
import numpy as np
import discretize as Mesh
from scipy.spatial import cKDTree, Delaunay
from scipy.interpolate import LinearNDInterpolator

def meshBuilder(xyz, h, padDist, meshGlobal=None,
                expFact=1.3,
                meshType='TENSOR',
                verticalAlignment='top'):
    """
        Function to quickly generate a Tensor mesh
        given a cloud of xyz points, finest core cell size
        and padding distance.
        If a meshGlobal is provided, the core cells will be centered
        on the underlaying mesh to reduce interpolation errors.

        :param numpy.ndarray xyz: n x 3 array of locations [x, y, z]
        :param numpy.ndarray h: 1 x 3 cell size for the core mesh
        :param numpy.ndarray padDist: 2 x 3 padding distances [W,E,S,N,Down,Up]
        [OPTIONAL]
        :param numpy.ndarray padCore: Number of core cells around the xyz locs
        :object SimPEG.Mesh: Base mesh used to shift the new mesh for overlap
        :param float expFact: Expension factor for padding cells [1.3]
        :param string meshType: Specify output mesh type: "TensorMesh"

        RETURNS:
        :object SimPEG.Mesh: Mesh object

    """

    assert meshType in ['TENSOR', 'TREE'], ('Revise meshType. Only ' +
                                            ' TENSOR | TREE mesh ' +
                                            'are implemented')

    # Get extent of points
    limx = np.r_[xyz[:, 0].max(), xyz[:, 0].min()]
    limy = np.r_[xyz[:, 1].max(), xyz[:, 1].min()]
    limz = np.r_[xyz[:, 2].max(), xyz[:, 2].min()]

    # Get center of the mesh
    midX = np.mean(limx)
    midY = np.mean(limy)

    if verticalAlignment == 'center':
        midZ = np.mean(limz)
    else:
        midZ = limz[0]

    nCx = int(limx[0]-limx[1]) / h[0]
    nCy = int(limy[0]-limy[1]) / h[1]
    nCz = int(np.max([
            limz[0]-limz[1],
            int(np.min(np.r_[nCx*h[0], nCy*h[1]])/3)
            ]) / h[2])

    if meshType == 'TENSOR':
        # Make sure the core has odd number of cells for centereing
        # on global mesh
        if meshGlobal is not None:
            nCx += 1 - int(nCx % 2)
            nCy += 1 - int(nCy % 2)
            nCz += 1 - int(nCz % 2)

        # Figure out paddings
        def expand(dx, pad):
            L = 0
            nC = 0
            while L < pad:
                nC += 1
                L = np.sum(dx * expFact**(np.asarray(range(nC))+1))

            return nC

        # Figure number of padding cells required to fill the space
        npadEast = expand(h[0], padDist[0, 0])
        npadWest = expand(h[0], padDist[0, 1])
        npadSouth = expand(h[1], padDist[1, 0])
        npadNorth = expand(h[1], padDist[1, 1])
        npadDown = expand(h[2], padDist[2, 0])
        npadUp = expand(h[2], padDist[2, 1])

        # Create discretization
        hx = [(h[0], npadWest, -expFact),
              (h[0], nCx),
              (h[0], npadEast, expFact)]
        hy = [(h[1], npadSouth, -expFact),
              (h[1], nCy), (h[1],
              npadNorth, expFact)]
        hz = [(h[2], npadDown, -expFact),
              (h[2], nCz),
              (h[2], npadUp, expFact)]

        # Create mesh
        mesh = Mesh.TensorMesh([hx, hy, hz], 'CC0')

        # Re-set the mesh at the center of input locations
        # Set origin
        if verticalAlignment == 'center':
            mesh.x0 = [midX-np.sum(mesh.hx)/2., midY-np.sum(mesh.hy)/2., midZ - np.sum(mesh.hz)/2.]
        elif verticalAlignment == 'top':
            mesh.x0 = [midX-np.sum(mesh.hx)/2., midY-np.sum(mesh.hy)/2., midZ - np.sum(mesh.hz)]
        else:
            assert NotImplementedError("verticalAlignment must be 'center' | 'top'")

    elif meshType == 'TREE':

        # Figure out full extent required from input
        extent = np.max(np.r_[nCx * h[0] + padDist[0, :].sum(),
                              nCy * h[1] + padDist[1, :].sum(),
                              nCz * h[2] + padDist[2, :].sum()])

        maxLevel = int(np.log2(extent/h[0]))+1

        # Number of cells at the small octree level
        # equal in 3D
        nCx, nCy, nCz = 2**(maxLevel), 2**(maxLevel), 2**(maxLevel)

        # nCy = 2**(int(np.log2(extent/h[1]))+1)
        # nCz = 2**(int(np.log2(extent/h[2]))+1)

        # Define the mesh and origin
        mesh = Mesh.TreeMesh([np.ones(nCx)*h[0],
                              np.ones(nCx)*h[1],
                              np.ones(nCx)*h[2]])

        # Shift mesh if global mesh is used
        center = np.r_[midX, midY, midZ]
        if meshGlobal is not None:

            tree = cKDTree(meshGlobal.gridCC)
            _, ind = tree.query(center, k=1)
            center = meshGlobal.gridCC[ind, :]

        # Set origin
        if verticalAlignment == 'center':
            mesh.x0 = np.r_[center[0] - (nCx)*h[0]/2., center[1] - (nCy)*h[1]/2., center[2] - (nCz)*h[2]/2.]
        elif verticalAlignment == 'top':
            mesh.x0 = np.r_[center[0] - (nCx)*h[0]/2., center[1] - (nCy)*h[1]/2., center[2] - (nCz)*h[2]]
        else:
            assert NotImplementedError("verticalAlignment must be 'center' | 'top'")

    return mesh


def refineTree(
            mesh, xyz,
            finalize=False, dtype="radial",
            octreeLevels=[1, 1, 1],
            octreeLevels_XY=None,
            maxDist=np.inf,
            maxLevel=None
):

    if octreeLevels_XY is not None:

        assert len(octreeLevels_XY) == len(octreeLevels), "Arguments 'octreeLevels' and 'octreeLevels_XY' must be the same length"

    else:

        octreeLevels_XY = np.zeros_like(octreeLevels)

    if maxLevel is None:
        maxLevel = int(np.log2(mesh.hx.shape[0]))

    tree = cKDTree(xyz)

    if dtype == "point":

        mesh.insert_cells(xyz, np.ones(xyz.shape[0])*maxLevel, finalize=False)

        stencil = np.r_[
                np.ones(octreeLevels[0]),
                np.ones(octreeLevels[1])*2,
                np.ones(octreeLevels[2])*3
        ]

        # Reflect in the opposite direction
        vec = np.r_[stencil[::-1], 1, stencil]
        vecX, vecY, vecZ = np.meshgrid(vec, vec, vec)
        gridLevel = np.maximum(np.maximum(np.abs(vecX),
                               np.abs(vecY)), np.abs(vecZ))
        gridLevel = np.kron(np.ones(xyz.shape[0]), mkvc(gridLevel))

        # Grid the coordinates
        vec = np.r_[-np.cumsum(stencil)[::-1], 0, np.cumsum(stencil)]
        vecX, vecY, vecZ = np.meshgrid(vec, vec, vec)
        offset = np.c_[
            mkvc(np.sign(vecX)*np.abs(vecX) * mesh.hx.min()),
            mkvc(np.sign(vecY)*np.abs(vecY) * mesh.hy.min()),
            mkvc(np.sign(vecZ)*np.abs(vecZ) * mesh.hz.min())
        ]

        # Replicate the point locations in each offseted grid points
        newLoc = (
            np.kron(xyz, np.ones((offset.shape[0], 1))) +
            np.kron(np.ones((xyz.shape[0], 1)), offset)
        )

        # Apply max distance
        r, ind = tree.query(newLoc)

        mesh.insert_cells(
            newLoc[r < maxDist], maxLevel-mkvc(gridLevel)[r < maxDist]+1, finalize=False
        )

        if finalize:
            mesh.finalize()

    elif dtype == "radial":

        # Compute the outer limits of each octree level
        rMax = np.cumsum(
            mesh.hx.min() *
            np.asarray(octreeLevels) *
            2**np.arange(len(octreeLevels))
        )

        def inBall(cell):
            xyz = cell.center
            r, ind = tree.query(xyz)

            for ii, nC in enumerate(octreeLevels):

                if r < rMax[ii]:

                    return maxLevel-ii

            return 0

        mesh.refine(inBall, finalize=finalize)

    elif dtype == 'surface':

        # Compute centroid and
        centroid = np.mean(xyz, axis=0)

        # Largest outer point distance
        rOut = np.linalg.norm(np.r_[
            np.abs(centroid[0]-xyz[:, 0]).max(),
            np.abs(centroid[1]-xyz[:, 1]).max()
            ]
        )

        # Compute maximum depth of refinement
        zmax = np.cumsum(
            mesh.hz.min() *
            np.asarray(octreeLevels) *
            2**np.arange(len(octreeLevels))
        )
        padWidth = np.cumsum(
            mesh.hx.min() *
            np.asarray(octreeLevels_XY) *
            2**np.arange(len(octreeLevels_XY))
        )

        depth = zmax[-1]

        # Increment the vertical offset
        zOffset = 0
        xyPad = -1
        # Cycle through the Tree levels backward
        for ii in range(len(octreeLevels)-1, -1, -1):

            dx = mesh.hx.min() * 2**ii
            dy = mesh.hy.min() * 2**ii
            dz = mesh.hz.min() * 2**ii

            # Increase the horizontal extent of the surface
            if xyPad != padWidth[ii]:
                xyPad = padWidth[ii]

                # Calculate expansion for padding XY cells
                expFactor = (rOut + xyPad) / rOut
                xLoc = (xyz - centroid)*expFactor + centroid

                # Create a new triangulated surface
                tri2D = Delaunay(xLoc[:, :2])
                F = LinearNDInterpolator(tri2D, xLoc[:, 2])

            limx = np.r_[xLoc[:, 0].max(), xLoc[:, 0].min()]
            limy = np.r_[xLoc[:, 1].max(), xLoc[:, 1].min()]

            nCx = int(np.ceil((limx[0]-limx[1]) / dx))
            nCy = int(np.ceil((limy[0]-limy[1]) / dy))

            # Create a grid at the octree level in xy
            CCx, CCy = np.meshgrid(
                np.linspace(
                    limx[1], limx[0], nCx
                    ),
                np.linspace(
                    limy[1], limy[0], nCy
                    )
            )

            xy = np.c_[CCx.reshape(-1), CCy.reshape(-1)]

            # Only keep points within triangulation
            indexTri = tri2D.find_simplex(xy)

            z = F(xy[indexTri != -1])

            newLoc = np.c_[xy[indexTri != -1], z]

            # Only keep points within maxDist
            # Apply max distance
            r, ind = tree.query(newLoc)

            # Apply vertical padding for current octree level
            zOffset = 0
            while zOffset < depth:
                indIn = r < (maxDist + padWidth[ii])
                nnz = int(np.sum(indIn))
                if nnz > 0:
                    mesh.insert_cells(
                        np.c_[
                            newLoc[indIn, :2],
                            newLoc[indIn, 2]-zOffset],
                        np.ones(nnz)*maxLevel-ii,
                        finalize=False
                    )

                zOffset += dz

            depth -= dz * octreeLevels[ii]

        if finalize:
            mesh.finalize()

    elif dtype == 'box':

        # Define the data extend [bottom SW, top NE]
        bsw = np.min(xyz, axis=0)
        tne = np.max(xyz, axis=0)

        hx = mesh.hx.min()

        if mesh.dim == 2:
            hz = mesh.hy.min()
        else:
            hz = mesh.hz.min()

        # Pre-calculate max depth of each level
        zmax = np.cumsum(
            hz * np.asarray(octreeLevels) *
            2**np.arange(len(octreeLevels))
        )

        if mesh.dim == 2:
            # Pre-calculate outer extent of each level
            padWidth = np.cumsum(
                    mesh.hx.min() *
                    np.asarray(octreeLevels_XY) *
                    2**np.arange(len(octreeLevels_XY))
                )

            # Make a list of outer limits
            BSW = [
                bsw - np.r_[padWidth[ii], zmax[ii]]
                for ii, (octZ, octXY) in enumerate(
                        zip(octreeLevels, octreeLevels_XY)
                )
            ]

            TNE = [
                tne + np.r_[padWidth[ii], zmax[ii]]
                for ii, (octZ, octXY) in enumerate(
                    zip(octreeLevels, octreeLevels_XY)
                )
            ]

        else:
            hy = mesh.hy.min()

            # Pre-calculate outer X extent of each level
            padWidth_x = np.cumsum(
                    hx * np.asarray(octreeLevels_XY) *
                    2**np.arange(len(octreeLevels_XY))
                )

            # Pre-calculate outer Y extent of each level
            padWidth_y = np.cumsum(
                    hy * np.asarray(octreeLevels_XY) *
                    2**np.arange(len(octreeLevels_XY))
                )

            # Make a list of outer limits
            BSW = [
                bsw - np.r_[padWidth_x[ii], padWidth_y[ii], zmax[ii]]
                for ii, (octZ, octXY) in enumerate(
                        zip(octreeLevels, octreeLevels_XY)
                )
            ]

            TNE = [
                tne + np.r_[padWidth_x[ii], padWidth_y[ii], zmax[ii]]
                for ii, (octZ, octXY) in enumerate(
                    zip(octreeLevels, octreeLevels_XY)
                )
            ]

        def inBox(cell):

            xyz = cell.center

            for ii, (nC, bsw, tne) in enumerate(zip(octreeLevels, BSW, TNE)):

                if np.all([xyz > bsw, xyz < tne]):
                    return maxLevel-ii

            return cell._level

        mesh.refine(inBox, finalize=finalize)

    else:
        NotImplementedError(
            "Only dtype= 'surface' | 'points' "
            " | 'radial' | 'box' has been implemented"
        )

    return mesh
