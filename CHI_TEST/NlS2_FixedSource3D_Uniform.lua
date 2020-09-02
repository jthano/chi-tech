--############################################### Setup mesh
chiMeshHandlerCreate()

nodes={}
N=3
ds=1.0/N
for i=0,N do
    nodes[i+1] = -0.5 + i*ds
end
--surf_mesh,region1 = chiMeshCreate1DSlabMesh(nodes)
surf_mesh,region1 = chiMeshCreate3DOrthoMesh(nodes,nodes,nodes)

--chiSurfaceMesherSetProperty(PARTITION_X,2)
--chiSurfaceMesherSetProperty(PARTITION_Y,2)
--chiSurfaceMesherSetProperty(CUT_X,0.0)
--chiSurfaceMesherSetProperty(CUT_Y,0.0)

chiVolumeMesherExecute();

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add material
material0 = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)
num_groups = 1
chiPhysicsMaterialSetProperty(material0,
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,
                              num_groups,     --Num grps
                              1.0,   --Sigma_t
                              0.2)   --Scattering ratio

--############################################### Setup Physics
phys1 = chiLBSCreateSolver_nlS2()
chiSolverAddRegion(phys1,region1)

for k=1,num_groups do
    chiLBSCreateGroup(phys1)
end

pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2,2)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
chiLBSGroupsetAddGroups(phys1,gs0,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,gs0,pquad)

--========== Boundary conditions
bsrc = {}
bsrc[1] = 0.5
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
                  YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)

--############################################### Setup Output
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

-- cline = chiFFInterpolationCreate(LINE)
-- chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,-1.0,1.0)
-- chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0, 1.0,1.0)
-- chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)
--
-- chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])
-- chiFFInterpolationInitialize(cline)
--
-- chiFFInterpolationExecute(cline)
-- chiFFInterpolationExportPython(cline)

chiExportFieldFunctionToVTK(fflist[1],"NlS2_FixedSource3D","Phi")

-- if (chi_location_id == 0) then
--     local handle = io.popen("python ZLFFI00.py")
-- end
