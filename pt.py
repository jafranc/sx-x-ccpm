from ctypes.wintypes import tagMSG

import vtk
import pyvista as pv
import numpy as np
import vtkmodules.numpy_interface.dataset_adapter as dsa
import matplotlib.pyplot as plt
import argparse

def formula(f : float, a: float = 2.*np.pi, off: float = 0.):

    return np.arccos(((4*np.pi - f - off )/a-1)) * 180/np.pi

if __name__ == "__main__":
    parser = argparse.ArgumentParser( prog='Post-treat ccpm',
            description='project phase contact from vertices to cell in VTK format')
    parser.add_argument('filename',nargs='+')
    parser.add_argument('-s','--stat',action='store_true')
    parser.add_argument('-n','--ndisc',nargs=1, type=int)
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
        n = args.ndisc[0]
        r = pv.read(args.filename[0])
        t = np.loadtxt(args.filename[0].split('.')[0]+'.csv',delimiter=',',skiprows=1)
        #sorting csv and stl points to match
        ixc = t[:,0].argsort()
        pts = dsa.numpy_support.vtk_to_numpy(r.GetPoints().GetData())
        ixm = pts[:,0].argsort()
        u = np.ndarray(t.shape)
        for i,ix in enumerate(ixm):
            u[ix,:] = t[ixc[i],:]

        grid = r.cast_to_unstructured_grid()
        tag = dsa.numpyTovtkDataArray(u[:, -2])
        tag.SetName("boundary")
        grid.GetPointData().AddArray(tag)
        #
        meanC = dsa.numpyTovtkDataArray( .5*(u[:, -5] * u[:, -6]) )
        meanC.SetName("meanCurvature")
        grid.GetPointData().AddArray(meanC)
        gaussianC = dsa.numpyTovtkDataArray( u[:, -5] * u[:, -6] )
        gaussianC.SetName("gaussianCurvature")
        grid.GetPointData().AddArray(gaussianC)
        mC = dsa.numpyTovtkDataArray( u[:, -4] )
        mC.SetName("mc")
        grid.GetPointData().AddArray(mC)
        gC = dsa.numpyTovtkDataArray( u[:, -3] )
        gC.SetName("gc")
        grid.GetPointData().AddArray(gC)
        #
        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputData(grid)
        p2c.Update()
        w = vtk.vtkXMLUnstructuredGridWriter()
        w.SetInputData(p2c.GetOutput())
        w.SetFileName(args.filename[0].split('.')[0]+'.vtu')
        w.Update()
        w.Write()
        #w = vtk.vtkXMLPolyDataWriter()

        qual = vtk.vtkMeshQuality()
        qual.SetInputData(p2c.GetOutput())
        qual.SetTriangleQualityMeasureToArea()
        qual.Update()
        tag = 4.8
        area = dsa.numpy_support.vtk_to_numpy(qual.GetOutput().GetCellData().GetArray(5))
        bound = dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('boundary'))
        print(f"global angle : {4*np.pi} : { np.sum(area*dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gc'))) } : "
              f"{ np.sum(area*( dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gaussianCurvature')) + 0.*dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('meanCurvature'))))}")
        print(f"ring value : { (f1:=np.sum(area[bound>=tag]*dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gc'))[bound>=tag])) } : "
              f"{ (f2:=np.sum(area[bound>=tag]*( dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gaussianCurvature'))[bound>=tag] + dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('meanCurvature'))[bound>=tag]))) }")
        print(f"formula : {formula(f1)} : {formula(f2)}")

        print(f"formula from geometry")
        cc = (vtk.vtkCellCenters())
        cc.SetInputData(p2c.GetOutput())
        cc.Update()
        pts = dsa.numpy_support.vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
        eps = 2/128*n
        # eps = 4
        ymin = np.max(pts[:,1])
        idx = np.intersect1d(np.where( pts[:,1]>ymin-eps)[0],np.where(pts[:,1]<ymin+eps )[0])
        print(f"ring value : { (f1:=np.sum(area[idx]*dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gc'))[idx])) } : "
              f"{ (f2:=np.sum(area[idx]*( dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gaussianCurvature'))[idx] + dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('meanCurvature'))[idx]))) }")
        print(f"formula : {formula(f1)} : {formula(f2)}")
        # p = .25
        # lastIndex = int(area[bound >= 4.1].shape[0] * p)
        # f1p = np.sum(area[bound>=4.1][:lastIndex] * dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gc'))[bound >= 4.1][:lastIndex])
        # f2p = np.sum(area[bound>=4.1][:lastIndex] * (dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gaussianCurvature'))[bound >= 4.1][:lastIndex]
        #                                         + dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('meanCurvature'))[bound>=4.1][:lastIndex]))
        # print(f"ring value (p={p}) : {f1p} : {f2p} ")
        # print(f"formula (p={p}) ): {formula(f1p, np.pi, +3/2*np.pi)} : {formula(f2p,np.pi,0)}")

        #np.arccos(((4*np.pi - 0.88 )/2/np.pi-1)) * 180/np.pi