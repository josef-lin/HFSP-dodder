import cv2
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import heapq
from collections import deque


def load_file(g, t):
    xz_path = (f"../../data/synthetic/synthetic_curves/gen_{g:02d}/plant_gen{g:02d}_t{t:02d}_xz.png")
    yz_path = (f"../../data/synthetic/synthetic_curves/gen_{g:02d}/plant_gen{g:02d}_t{t:02d}_yz.png")
    xy_path = (f"../../data/synthetic/synthetic_curves/gen_{g:02d}/plant_gen{g:02d}_t{t:02d}_xy.png")


    # load images
    img_xz = cv2.imread(xz_path, cv2.IMREAD_GRAYSCALE)
    img_yz = cv2.imread(yz_path, cv2.IMREAD_GRAYSCALE)
    img_xy = cv2.imread(xy_path, cv2.IMREAD_GRAYSCALE)

    if img_xz is None:
        raise FileNotFoundError(xz_path)

    if img_yz is None:
        raise FileNotFoundError(yz_path)

    if img_xy is None:
        raise FileNotFoundError(xy_path)

    # convert to [0,1]
    # rho_xz and rho_yz are the probability distributions along the xz and yz planes
    # (x,z) and (y,z) respectively. The .T is to transpose the arrays so that the z-axis is the first axis.
    rho_xz = (img_xz > 0).astype(np.float32).T
    rho_yz = (img_yz > 0).astype(np.float32).T
    rho_xy = (img_xy > 0).astype(np.float32).T

    # optional denoising / uncertainty model
    # rho_xz = cv2.GaussianBlur(rho_xz, (0, 0), sigmaX=2)
    # rho_yz = cv2.GaussianBlur(rho_yz, (0, 0), sigmaX=2)

    # construct R0
    R0 = rho_xz[:, None, :] * rho_yz[None, :, :] * rho_xy[:, :, None]

    return R0

def create_mask(R, threshold):
    return R > threshold

def distance_transform(mask):
    labels, n = ndimage.label(mask)

    sizes = ndimage.sum(mask, labels, range(1, n+1))
    largest = np.argmax(sizes) + 1

    mask = labels == largest

    return ndimage.distance_transform_edt(mask)

# 26-neighbour connectivity
NEIGHBOURS = np.array([
    [dx,dy,dz]
    for dx in (-1,0,1)
    for dy in (-1,0,1)
    for dz in (-1,0,1)
    if not (dx == dy == dz == 0)
], dtype=int)


def bfs_farthest(mask, start):
    """
    BFS restricted to foreground voxels.
    
    Returns:
        farthest_voxel
        distance_array
        parent_dictionary
    """
    
    dist = {}
    parent = {}

    q = deque([tuple(start)])
    dist[tuple(start)] = 0

    farthest = tuple(start)

    nx, ny, nz = mask.shape

    while q:

        p = q.popleft()

        if dist[p] > dist[farthest]:
            farthest = p

        x, y, z = p

        for dx, dy, dz in NEIGHBOURS:

            xn = x + dx
            yn = y + dy
            zn = z + dz

            if (
                0 <= xn < nx and
                0 <= yn < ny and
                0 <= zn < nz and
                mask[xn, yn, zn]
            ):

                qn = (xn, yn, zn)

                if qn not in dist:
                    dist[qn] = dist[p] + 1
                    parent[qn] = p
                    q.append(qn)

    return farthest, dist, parent

def centered_geodesic(mask, D, start, goal):

    pq = []

    heapq.heappush(pq, (0.0, tuple(start)))

    dist = {tuple(start): 0.0}
    parent = {}

    nx, ny, nz = mask.shape

    while pq:

        curr_cost, p = heapq.heappop(pq)

        if p == tuple(goal):
            break

        if curr_cost > dist[p]:
            continue

        x,y,z = p

        for dx,dy,dz in NEIGHBOURS:

            xn = x + dx
            yn = y + dy
            zn = z + dz

            if (
                0 <= xn < nx and
                0 <= yn < ny and
                0 <= zn < nz and
                mask[xn,yn,zn]
            ):

                q = (xn,yn,zn)

                step = np.sqrt(dx*dx + dy*dy + dz*dz)

                edge_cost = step/(D[q]**2 + 1e-3)

                new_cost = curr_cost + edge_cost

                if q not in dist or new_cost < dist[q]:

                    dist[q] = new_cost
                    parent[q] = p

                    heapq.heappush(
                        pq,
                        (new_cost, q)
                    )

    path = []

    p = tuple(goal)

    while p != tuple(start):

        path.append(p)
        p = parent[p]

    path.append(tuple(start))

    path.reverse()

    return np.array(path)


def find_endpoints(mask):
    # any foreground voxel
    start = np.argwhere(mask)[0]

    # first BFS
    p1, _, _ = bfs_farthest(mask, start)

    # second BFS
    p2, dist, parent = bfs_farthest(mask, p1)
    return p1, p2 

def voxel_to_physical(g, R, path):
    camera_params_path = f"../../data/synthetic/synthetic_curves/gen_{g:02d}/camera_params.csv"

    camera_params = pd.read_csv(camera_params_path)

    scale = camera_params['scale'].values[0]
    x_offset = camera_params['x_offset'].values[0]
    y_offset = camera_params['y_offset'].values[0]
    z_offset = camera_params['z_offset'].values[0]
    imgH = camera_params['imgH'].values[0]
    imgW = camera_params['imgW'].values[0]

    Nx, Ny, Nz = R.shape

    x_phys = (np.arange(Nx) - x_offset) / scale
    y_phys = (np.arange(Ny) - y_offset) / scale
    z_phys = ((Nz - 1) - np.arange(Nz) - z_offset) / scale

    path_phys = np.column_stack([
        x_phys[path[:,0].astype(int)],
        y_phys[path[:,1].astype(int)],
        z_phys[path[:,2].astype(int)]
    ])

    # Then recenter the curve to be at origin 
    x0, y0, z0 = path_phys[0]
    path_phys -= np.array([x0, y0, z0])

    return np.asarray(path_phys, dtype=float)

def smooth_centerline(path):
    return np.column_stack([
        savgol_filter(path[:,0], 21, 3),
        savgol_filter(path[:,1], 21, 3),
        savgol_filter(path[:,2], 21, 3)
    ])

def arclength_resample(centerline, N=100):

    ds = np.linalg.norm(
        np.diff(centerline, axis=0),
        axis=1
    )

    keep = np.concatenate([[True], ds > 1e-10])
    centerline = centerline[keep]

    ds = np.linalg.norm(
        np.diff(centerline, axis=0),
        axis=1
    )

    s = np.concatenate([[0], np.cumsum(ds)])

    s_uniform = np.linspace(0, s[-1], N)

    fx = interp1d(s, centerline[:,0], kind='cubic')
    fy = interp1d(s, centerline[:,1], kind='cubic')
    fz = interp1d(s, centerline[:,2], kind='cubic')

    return np.column_stack([
        fx(s_uniform),
        fy(s_uniform),
        fz(s_uniform)
    ])

def extract_centerline(g, t, N):
    R = load_file(g, t)
    mask = create_mask(R, 0)
    D = distance_transform(mask)
    p1, p2 = find_endpoints(mask)
    centerline = centered_geodesic(mask, D, p1, p2)
    centerline = voxel_to_physical(g, R, centerline)
    centerline = smooth_centerline(centerline)
    return arclength_resample(centerline, N)