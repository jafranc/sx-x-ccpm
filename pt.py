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

        qual = vtk.vtkMeshQuality()
        qual.SetInputData(p2c.GetOutput())
        qual.SetTriangleQualityMeasureToArea()
        qual.Update()

        tag = 5.1
        #selection work
        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL)
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
        ids = vtk.vtkIdTypeArray()
        for icell in range(qual.GetOutput().GetNumberOfCells()):
            if(qual.GetOutput().GetCellData().GetArray('boundary').GetValue(icell)>=tag):
                ids.InsertNextValue(icell)
        selectionNode.SetSelectionList(ids)
        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)
        extract = vtk.vtkExtractSelection()
        extract.SetInputConnection(0,qual.GetOutputPort())
        extract.SetInputData(1, selection)
        extract.Update()
        print(f"there is {extract.GetOutput().GetNumberOfCells()} cells in the selection")

        conn = vtk.vtkConnectivityFilter()
        conn.SetInputData(extract.GetOutput())
        conn.SetExtractionModeToAllRegions()
        conn.ColorRegionsOn()
        conn.Update()
        print(f"there is {conn.GetNumberOfExtractedRegions()} contacts detected")

        garea = dsa.numpy_support.vtk_to_numpy(qual.GetOutput().GetCellData().GetArray(3))
        region = dsa.numpy_support.vtk_to_numpy(conn.GetOutput().GetCellData().GetArray('RegionId'))
        area = dsa.numpy_support.vtk_to_numpy(conn.GetOutput().GetCellData().GetArray(3))

        print(f"global angle : {4*np.pi} : { np.sum(garea*dsa.numpy_support.vtk_to_numpy(p2c.GetOutput().GetCellData().GetArray('gc'))) }")

        number_of_valid: int = 0
        s,e = conn.GetOutput().GetCellData().GetArray('RegionId').GetRange()
        for i in np.arange(s,e+1):

            #eliminate minor contacts
            if(region[region==i].shape[0]>15):
                number_of_valid += 1

                print(f"{i} : number of elements {np.sum(region==i)} for a global area of {np.sum(area[region==i])}")
                print(f"{i} : ring value : { (f1:=np.sum(area[region == i]*dsa.numpy_support.vtk_to_numpy(conn.GetOutput().GetCellData().GetArray('gc'))[region == i])) }")
                print(f"{i} : formula : {formula(f1)}")

        print(f"A : ring value : { (f1:=np.sum(area*dsa.numpy_support.vtk_to_numpy(conn.GetOutput().GetCellData().GetArray('gc')))/number_of_valid) }")
        print(f"A : formula : {formula(f1)}")