import vtk
import pyvista as pv
import numpy as np
import vtkmodules.numpy_interface.dataset_adapter as dsa
import matplotlib.pyplot as plt
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser( prog='Post-treat ccpm',
            description='project phase contact from vertices to cell in VTK format')
    parser.add_argument('filename',nargs='+')
    parser.add_argument('-s','--stat',action='store_true')
    args = parser.parse_args()

    if args.stat:
        # 510 : tag for S+L+G, 425 for L+G, 340 for S+L and 255 for S+G
        fig,axs = plt.subplots(len(args.filename)+1,2)
        vsum = np.ndarray((0,6))
        for i,file in enumerate(args.filename):

            tag = None 
            # t = np.loadtxt(file,skiprows=1,delimiter=',')
            v = np.loadtxt(file,skiprows=1,delimiter=',')
            # u = np.loadtxt(file,skiprows=1,delimiter=',')
            if tag is not None:
                v510 = v[v[:,5]==tag,:]
            else:
                v510 = v
            v510 = v510[v510[:,4]!=0,:]
            v510 = v510[v510[:,3]!=0,:]
            vsum = np.concatenate((vsum, v510), axis=0)
            # t510 = t#[t[:,5]==tag,:]
            # t510 = t510[t510[:,4]!=0,:]
            # t510 = t510[t510[:,3]!=0,:]
            # u510 = u#[u[:,5]==tag,:]
            # u510 = u510[u510[:,4]!=0,:]
            # u510 = u510[u510[:,3]!=0,:]

            xh,yh,h = axs[i,0].hist(v510[:,3],bins=100,label=file)
            axs[i,0].legend()
            axs[i,0].set_xlim([-0.5,1.])
            # xf,yf,f = axs[1,0].hist(t510[:,3],bins=yh,label=args.filename[0])
            # axs[1,0].legend()
            # xg,yg,g = axs[2,0].hist(u510[:,3],bins=yh,label=args.filename[2])
            # axs[2,0].legend()

            xh,yh,h = axs[i,1].hist(v510[:,4],bins=100)
            # xf,yf,f = axs[1,1].hist(t510[:,4],bins=yh)
            # xg,yg,g = axs[2,1].hist(u510[:,4],bins=yh)
            axs[i,1].set_xlim([-0.5,1.])
        xh,yh,h = axs[-1,0].hist(vsum[:,3],bins=100,label='sum',color='r')
        axs[-1,0].legend()
        axs[-1,0].set_xlim([-0.5,1.])
        xh,yh,h = axs[-1,1].hist(vsum[:,4],bins=100,label='sum',color='r')
        axs[-1,1].legend()
        axs[-1,1].set_xlim([-0.5,1.])

        plt.show()

    else:
        r = pv.read(args.filename[0])
        t = np.loadtxt(args.filename[0].split('.')[0]+'.csv',delimiter=',')
        #sorting csv and stl points to match
        ixc = t[:,0].argsort()
        pts = dsa.numpy_support.vtk_to_numpy(r.GetPoints().GetData())
        ixm = pts[:,0].argsort()
        u = np.ndarray(t.shape)
        for i,ix in enumerate(ixm):
            u[ix,:] = t[ixc[i],:]

        grid = r.cast_to_unstructured_grid()
        grid.GetPointData().AddArray( dsa.numpyTovtkDataArray(u[:,-1]) )
        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputData(grid)
        p2c.Update()
        w = vtk.vtkXMLUnstructuredGridWriter()
        w.SetInputData(p2c.GetOutput())
        w.SetFileName(args.filename[0].split('.')[0]+'.vtu')
        w.Update()
        w.Write()
        #w = vtk.vtkXMLPolyDataWriter()
