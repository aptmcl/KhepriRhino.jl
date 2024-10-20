export rhino

### PLugin
#=
Whenever the plugin is updated, run this function and commit the plugin files.
upgrade_plugin()
=#

julia_khepri = dirname(dirname(abspath(@__FILE__)))

dlls = ["KhepriBase.dll", "KhepriRhinoceros.rhp"]
local_plugin = joinpath(dirname(dirname(abspath(@__FILE__))), "Plugin")

upgrade_plugin() =
  let # 1. The dlls are updated in VisualStudio after compilation of the plugin, and they are stored in the folder.
      # contained inside the Plugins folder, which has a specific location regarding this file itself
      plugin_folder = joinpath(dirname(dirname(dirname(dirname(abspath(@__FILE__))))), "Plugins", "KhepriRhinoceros", "KhepriRhinoceros", "bin")
      # 2. The bundle needs to be copied to the current folder
      local_folder = joinpath(julia_khepri, "Plugin")
      # 3. Now we copy the dlls to the local folder
      for dll in dlls
          src = joinpath(plugin_folder, dll)
          dst = joinpath(local_folder, dll)
          rm(dst, force=true)
          cp(src, dst)
      end
  end

##############################################
# We need some additional Encoders
# RH is a subtype of CS
parse_signature(::Val{:RH}, sig::T) where {T} = parse_signature(Val(:CS), sig)
encode(::Val{:RH}, t::Val{T}, c::IO, v) where {T} = encode(Val(:CS), t, c, v)
decode(::Val{:RH}, t::Val{T}, c::IO) where {T} = decode(Val(:CS), t, c)

# We need some additional Encoders
@encode_decode_as(:RH, Val{:Options}, Val{:Dict})
@encode_decode_as(:RH, Val{:Brep}, Val{:Guid})
@encode_decode_as(:RH, Val{:RhinoObject}, Val{:Guid})

encode(::Val{:RH}, t::Union{Val{:Point3d},Val{:Vector3d}}, c::IO, p) =
  encode(Val(:CS), Val(:double3), c, raw_point(p))
decode(::Val{:RH}, t::Val{:Point3d}, c::IO) =
  xyz(decode(Val(:CS), Val(:double3), c)..., world_cs)
decode(::Val{:RH}, t::Val{:Vector3d}, c::IO) =
  vxyz(decode(Val(:CS), Val(:double3), c)..., world_cs)

encode(ns::Val{:RH}, t::Val{:Plane}, c::IO, p) =
  let o = in_world(p),
      x = in_world(vx(1, p.cs)),
      y = in_world(vy(1, p.cs))
    encode(ns, Val(:Point3d), c, o)
    encode(ns, Val(:Vector3d), c, x)
    encode(ns, Val(:Vector3d), c, y)
  end
decode(ns::Val{:RH}, t::Val{:Plane}, c::IO) =
  loc_from_o_vx_vy(
    decode(ns, Val(:Point3d), c),
    decode(ns, Val(:Vector3d), c),
    decode(ns, Val(:Vector3d), c))

rhino_api = @remote_functions :RH """
public void RunScript(string script)
public void SetView(Point3d position, Point3d target, double lens, bool perspective, string style)
public void View(Point3d position, Point3d target, double lens)
public void ViewTop()
public Point3d ViewCamera()
public Point3d ViewTarget()
public double ViewLens()
public string DefaultMaterialsFolder()
public MatId LoadRenderMaterialFromPath(string path)
public MatId CreateMaterial(string name, Color diffuse, Color specular, Color ambient, Color emission, Color reflection, double reflectivity, double reflectionGlossiness, double refractionIndex, double transparency)
public MatId CreateColorMaterial(Color color)
public Guid Point(Point3d p)
public Point3d PointPosition(Guid ent)
public Guid PolyLine(Point3d[] pts)
public Point3d[] LineVertices(RhinoObject id)
public Guid Spline(Point3d[] pts)
public Guid SplineTangents(Point3d[] pts, Vector3d start, Vector3d end)
public Guid ClosedPolyLine(Point3d[] pts)
public Guid ClosedSpline(Point3d[] pts)
public Guid Circle(Point3d c, Vector3d n, double r)
public Point3d CircleCenter(RhinoObject obj)
public Vector3d CircleNormal(RhinoObject obj)
public double CircleRadius(RhinoObject obj)
public Guid Ellipse(Plane p, double radiusX, double radiusY)
public Guid Arc(Plane p, double radius, double startAngle, double endAngle)
public Guid Text(string str, Plane p, double height, string fontName, bool bold, bool italic)
public Guid JoinCurves(Guid[] objs)
public Guid Mesh(Point3d[] pts, int[][] faces, MatId mat)
public Guid NGon(Point3d[] pts, Point3d pivot, bool smooth, MatId mat)
public Guid SurfaceFrom(Guid[] ids, MatId mat)
public Guid SurfaceFromCurvesPoints(Point3d[][] ptss, bool[] smooths, MatId mat)
public Guid SurfaceCircle(Point3d c, Vector3d n, double r, MatId mat)
public Guid SurfaceEllipse(Plane p, double radiusX, double radiusY, MatId mat)
public Guid SurfaceArc(Plane p, double radius, double startAngle, double endAngle, MatId mat)
public Guid SurfaceClosedPolyLine(Point3d[] pts, MatId mat)
public Guid Sphere(Point3d c, double r, MatId mat)
public Guid Torus(Point3d c, Vector3d vz, double majorRadius, double minorRadius, MatId mat)
public Guid Cylinder(Point3d bottom, double radius, Point3d top, MatId mat)
public Guid Cone(Point3d bottom, double radius, Point3d top, MatId mat)
public Guid ConeFrustum(Point3d bottom, double bottom_radius, Point3d top, double top_radius, MatId mat)
public Guid Box(Point3d corner, Vector3d vx, Vector3d vy, double dx, double dy, double dz, MatId mat)
public Guid XYCenteredBox(Point3d corner, Vector3d vx, Vector3d vy, double dx, double dy, double dz, MatId mat)
public Guid IrregularPyramid(Point3d[] pts, Point3d apex, MatId mat)
public Guid IrregularPyramidFrustum(Point3d[] bpts, Point3d[] tpts, MatId mat)
public Guid PrismWithHoles(Point3d[][] ptss, bool[] smooths, Vector3d dir, MatId mat)
public Guid SurfaceFromGrid(int nU, int nV, Point3d[] pts, bool closedU, bool closedV, int degreeU, int degreeV, MatId mat)
public Brep[] Thicken(RhinoObject obj, double thickness)
public double[] CurveDomain(RhinoObject obj)
public double CurveLength(RhinoObject obj)
public Plane CurveFrameAt(RhinoObject obj, double t)
public Plane CurveFrameAtLength(RhinoObject obj, double l)
public double[] SurfaceDomain(RhinoObject obj)
public Plane SurfaceFrameAt(RhinoObject obj, double u, double v)
public Point3d[] BoundingBox(ObjectId[] ids)
public Brep Extrusion(RhinoObject obj, Vector3d dir)
public Brep[] SweepPathProfile(RhinoObject path, RhinoObject profile, double rotation, double scale)
public Guid[] Unite(RhinoObject[]objs)
public Guid[] Intersect(RhinoObject obj0, RhinoObject obj1)
public Guid[] Subtract(RhinoObject obj0, RhinoObject obj1)
public Guid[] Slice(RhinoObject obj, Point3d p, Vector3d n)
public Guid Move(Guid id, Vector3d v)
public Guid Scale(Guid id, Point3d p, double s)
public Guid Rotate(Guid id, Point3d p, Vector3d n, double a)
public Guid Mirror(Guid id, Point3d p, Vector3d n, bool copy)
public bool IsPoint(RhinoObject e)
public bool IsCircle(RhinoObject e)
public bool IsPolyLine(RhinoObject e)
public bool IsSpline(RhinoObject e)
public bool IsClosedPolyLine(RhinoObject e)
public bool IsClosedSpline(RhinoObject e)
public bool IsEllipse(RhinoObject e)
public bool IsArc(RhinoObject e)
public bool IsText(RhinoObject e)
public byte ShapeCode(RhinoObject id)
public Guid ClosedPathCurveArray(Point3d[] pts, double[] angles)
public Brep[] PathWall(RhinoObject obj, double lThickness, double rThickness, double height)
public Brep RectangularTable(Point3d c, double angle, double length, double width, double height, double top_thickness, double leg_thickness)
public Brep[] BaseRectangularTable(double length, double width, double height, double top_thickness, double leg_thickness)
public int CreateRectangularTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)
public Guid Table(Point3d c, double angle, int family)
public Brep[] BaseChair(double length, double width, double height, double seat_height, double thickness)
public int CreateChairFamily(double length, double width, double height, double seat_height, double thickness)
public Guid Chair(Point3d c, double angle, int family)
public GeometryBase[] BaseRectangularTableAndChairs(int tableFamily, int chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)
public int CreateRectangularTableAndChairsFamily(int tableFamily, int chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)
public Guid TableAndChairs(Point3d c, double angle, int family)
public Guid CreateLayer(String name, bool active, Color color)
public Guid CurrentLayer()
public void SetCurrentLayer(Guid id)
public void SetLayerActive(Guid id, bool active) {
public Guid ShapeLayer(RhinoObject objId)
public void SetShapeLayer(RhinoObject objId, Guid layerId)
public Point3d[] CurvePointsAt(RhinoObject obj, double[] ts)
public Vector3d[] CurveTangentsAt(RhinoObject obj, double[] ts)
public Vector3d[] CurveNormalsAt(RhinoObject obj, double[] ts)
public int DeleteAll()
public int DeleteAllInLayer(Guid id)
public void Delete(Guid id)
public void DeleteMany(Guid[] ids)
public void Select(Guid id)
public void SelectMany(Guid[] ids)
public void Highlight(RhinoObject o)
public void Unhighlight(RhinoObject o)
public void UnhighlightAll()
public Guid Clone(RhinoObject e)
public Point3d[] GetPosition(string prompt)
public Guid[] GetPoint(string prompt)
public Guid[] GetPoints(string prompt)
public Guid[] GetCurve(string prompt)
public Guid[] GetCurves(string prompt)
public Guid[] GetSurface(string prompt)
public Guid[] GetSurfaces(string prompt)
public Guid[] GetSolid(string prompt)
public Guid[] GetSolids(string prompt)
public Guid[] GetShape(string prompt)
public Guid[] GetShapes(string prompt)
public Guid[] GetAllShapes()
public Guid[] GetAllShapesInLayer(Guid id)
public Guid CreateLeaderDimension(String text, Point3d p0, Point3d p1, double scale, String mark, Options props)
public Guid CreateDiametricDimension(String text, Plane c, Point3d p, Point3d pdim, double scale, String mark, Options props)
public Guid CreateAngularDimension(String text, Plane p, Point3d p0, Point3d p1, Point3d ptdim, double scale, String mark, Options props) {
public void SunLight(DateTime dt, double latitude, double longitude, int meridian, float turbidity, bool hasSun)
public Guid PointLight(Point3d p, Color c, double power)
public void EnableGrasshopperSolver()
public void DisableGrasshopperSolver()
public void RunGrasshopperSolver()
public void Render(int width, int height, int quality, string path)
public void RenderLoadKhepriEnvironment(string name, string path)
public void RenderUseHDRiEnvironment(string name, double rotation)
public void ClayRenderBlack(string env, double rotation, int width, int height, int quality, string path)
public void ClayRenderWhite(string env, double rotation, int width, int height, int quality, string path)
public void SetBackgroundHDRi(string path, double angle)
public void SaveView(int width, int height, string path)
"""
#public Guid CreateAngularDimension(string text, Plane p, double radius, double startAngle, double endAngle, double scale, double offset, Options props)

#=
public byte Sync()
public byte Disconnect()
public Entity InterpSpline(Point3d[] pts, Vector3d tan0, Vector3d tan1)
public Entity InterpClosedSpline(Point3d[] pts)
public Guid SurfaceFromCurve(Entity curve)
public Guid IrregularPyramidMesh(Point3d[] pts, Point3d apex)
public Entity SolidFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level, double thickness)
public Guid NurbSurfaceFrom(ObjectId id)
public double[] SurfaceDomain(Entity ent)
public Frame3d SurfaceFrameAt(Entity ent, double u, double v)
public Guid Sweep(ObjectId pathId, ObjectId profileId, double rotation, double scale)
public Guid Loft(ObjectId[] profilesIds, ObjectId[] guidesIds, bool ruled, bool closed)
public Guid Revolve(ObjectId profileId, Point3d p, Vector3d n, double startAngle, double amplitude)
public Point3d[] GetPoint(string prompt)
public void ZoomExtents()
public void SetSystemVariableInt(string name, int value)
public int Command(string cmd)
public void DisableUpdate()
public void EnableUpdate()
public BIMLevel FindOrCreateLevelAtElevation(double elevation)
public BIMLevel UpperLevel(BIMLevel currentLevel, double addedElevation)
public double GetLevelElevation(BIMLevel level)
public Entity MeshFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN)

public FloorFamily FloorFamilyInstance(double totalThickness, double coatingThickness)
public Entity LightweightPolyLine(Point2d[] pts, double[] angles, double elevation)
public Entity SurfaceLightweightPolyLine(Point2d[] pts, double[] angles, double elevation)
public ObjectId CreatePathFloor(Point2d[] pts, double[] angles, BIMLevel level, FloorFamily family)
=#


abstract type RHKey end
const RHId = Guid
const RHRef = GenericRef{RHKey, RHId}
const RHRefs = Vector{RHRef}
const RHNativeRef = NativeRef{RHKey, RHId}
const RHNativeRefs = NativeRefs{RHKey, RHId}
const RH = SocketBackend{RHKey, RHId}

KhepriBase.void_ref(::RH) = 0 % UInt128

KhepriBase.failed_connecting(b::RH) =
  begin
	@info("Please, ensure Khepri's Rhinoceros plugin is installed.")
  	@info("Go to Rhinoceros: Tools->Options->Plug-ins->Install")
  	@info("Location: $(joinpath(local_plugin, dlls[2]))")
  end

KhepriBase.after_connecting(b::RH) =
  let mat_folder = @remote(b, DefaultMaterialsFolder()),
  	  m = match(r".+?Rhinoceros\\(.+?)\\Localization", mat_folder)
  	if m === nothing
	  error("Unrecognized Rhino materials folder: $(mat_folder)")
  	else
	  let v = VersionNumber(m[1])
	    if v >= v"6.0"
			# Apparently, materials are not yet working in Rhino 6!!!!
	      set_material(rhino, material_metal, rhino_default_material(raw"Metal\Matte\Matte Silver"))
	      set_material(rhino, material_glass, rhino_default_material(raw"Glass\Clear Glass"))
	      set_material(rhino, material_wood, rhino_default_material(raw"Wood\Pear polished"))
	      set_material(rhino, material_concrete, rhino_default_material(raw"Ceramics\Stoneware"))
	      set_material(rhino, material_plaster, rhino_default_material(raw"White Matte"))
	      set_material(rhino, material_grass, rhino_default_material(raw"Textures\Grass"))
	    elseif v >= v"5.0"
          set_material(rhino, material_metal, rhino_default_material(raw"Metal\Silver"))
          set_material(rhino, material_glass, rhino_default_material(raw"Transparent\Glass"))
          set_material(rhino, material_wood, rhino_default_material(raw"Wood\Scots Pine"))
          set_material(rhino, material_concrete, rhino_default_material(raw"Stone\Granite"))
          set_material(rhino, material_plaster, rhino_default_material(raw"White Matte"))
          set_material(rhino, material_grass, rhino_default_material(raw"Special\Grass"))
	    else
	      error("Unrecognized Rhino version: ")
	    end
	  end
    end
  end

const rhino = RH("Rhino", rhino_port, rhino_api)

KhepriBase.has_boolean_ops(::Type{RH}) = HasBooleanOps{false}()
KhepriBase.backend(::RHRef) = rhino

# Primitives
KhepriBase.b_point(b::RH, p, mat) =
  @remote(b, Point(p))

KhepriBase.b_line(b::RH, ps, mat) =
  @remote(b, PolyLine(ps))

KhepriBase.b_polygon(b::RH, ps, mat) =
  @remote(b, ClosedPolyLine(ps))

KhepriBase.b_spline(b::RH, ps, v0, v1, mat) =
  if (v0 == false) && (v1 == false)
    @remote(b, Spline(ps))
  elseif (v0 != false) && (v1 != false)
    @remote(b, SplineTangents(ps, v0, v1))
  else
    @remote(b, SplineTangents(ps,
                 v0 == false ? ps[2]-ps[1] : v0,
                 v1 == false ? ps[end-1]-ps[end] : v1))
  end

KhepriBase.b_closed_spline(b::RH, ps, mat) =
  @remote(b, ClosedSpline(ps))

KhepriBase.b_circle(b::RH, c, r, mat) =
  @remote(b, Circle(c, vz(1, c.cs), r))

KhepriBase.b_arc(b::RH, c, r, α, Δα, mat) =
  if r == 0
    @remote(b, Point(c))
  elseif Δα == 0
    @remote(b, Point(c + vpol(r, α, c.cs)))
  elseif abs(Δα) >= 2*pi
    @remote(b, Circle(c, vz(1, c.cs), r))
  else
	let β = α + Δα
  	  if β > α
  	  	@remote(b, Arc(c, r, α, β))
  	  else
  	  	@remote(b, Arc(c, r, β, α))
  	  end
    end
  end

KhepriBase.b_ellipse(b::RH, c, rx, ry, mat) =
  @remote(b, Ellipse(c, rx, ry))

#=
Rhino supports meshes and nurbs surfaces. Meshes are more efficient... but 
many operations, such as boolean operations, sweep, extrusion, etc., need 
surfaces (or polysurfaces).
This means that although we might implement trigs and quads using meshes, 
that is not what Rhino would give us if we were creating 3D objects by hand.
Therefore, we will use meshes as a fallback but most of the surface and solids
stuff needs Rhino's surfaces and polysurfaces.
=#

KhepriBase.b_trig(b::RH, p1, p2, p3, mat) =
  @remote(b, Mesh([p1, p2, p3], [[0, 1, 2, 2]], mat))

KhepriBase.b_quad(b::RH, p1, p2, p3, p4, mat) =
 	@remote(b, Mesh([p1, p2, p3, p4], [[0, 1, 2, 2], [0, 2, 3, 3]], mat))

KhepriBase.b_ngon(b::RH, ps, pivot, smooth, mat) =
 	@remote(b, NGon(ps, pivot, smooth, mat))

############################################################
# Second tier: surfaces

KhepriBase.b_surface_rectangle(b::RH, c, dx, dy, mat) =
  @remote(b, SurfaceClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy), c], mat))

KhepriBase.b_surface_regular_polygon(b::RH, edges, c, r, angle, inscribed, mat) =
  b_surface_polygon(b, regular_polygon_vertices(edges, c, r, angle, inscribed), mat)

KhepriBase.b_surface_polygon(b::RH, ps, mat) =
  @remote(b, SurfaceClosedPolyLine([ps...,ps[1]], mat))

KhepriBase.b_surface_polygon_with_holes(b::RH, ps, qss, mat) =
  @remote(b, SurfaceFromCurvesPoints([ps, qss...], falses(1 + length(qss)), mat))

KhepriBase.b_surface_circle(b::RH, c, r, mat) =
	@remote(b, SurfaceCircle(c, vz(1, c.cs), r, mat))

KhepriBase.b_surface_arc(b::RH, c, r, α, Δα, mat) =
  @remote(b, SurfaceArc(c, r, α, α + Δα, mat))

KhepriBase.b_surface_ellipse(b::RH, c, rx, ry, mat) =
   @remote(b, SurfaceEllipse(c, rx, ry, mat))

KhepriBase.b_surface(b::RH, frontier::Shapes, mat) =
  and_mark_deleted(b, @remote(b, SurfaceFrom(collect_ref(b, frontier), mat)), frontier)


KhepriBase.b_surface_grid(b::RH, ptss, closed_u, closed_v, smooth_u, smooth_v, mat) =
  let (nu, nv) = size(ptss),
      order(n) = min(2*ceil(Int,n/16) + 1, 5)
    @remote(b, SurfaceFromGrid(nu, nv,
                               reshape(permutedims(ptss), :),
                               closed_u, closed_v,
                               smooth_u ? order(nu) : 1,
                               smooth_v ? order(nv) : 1,
							   mat))
  end

KhepriBase.b_surface_mesh(b::RH, vertices, faces, mat) =
  let ensure_4(v) = length(v) == 4 ? v : [v..., v[3]]
    @remote(b, Mesh(vertices, map(face->ensure_4(face.-1), faces), mat))
  end

############################################################
# Third tier: solids

KhepriBase.b_solidify(b::RH, refs) = refs
  

KhepriBase.b_generic_prism(b::RH, bs, smooth, v, bmat, tmat, smat) =
	@remote(b, PrismWithHoles([bs], [smooth], v, tmat))

KhepriBase.b_generic_prism_with_holes(b::RH, bs, smooth, bss, smooths, v, bmat, tmat, smat) =
  @remote(b, PrismWithHoles([bs, bss...], [smooth, smooths...], v, tmat))

KhepriBase.b_pyramid_frustum(b::RH, bs, ts, bmat, tmat, smat) =
  @remote(b, IrregularPyramidFrustum(bs, ts, tmat))

KhepriBase.b_pyramid(b::RH, bs, t, bmat, smat) =
  @remote(b, IrregularPyramid(bs, t, smat))

KhepriBase.b_cylinder(b::RH, c, r, h, bmat, tmat, smat) =
  @remote(b, Cylinder(c, r, c + vz(h, c.cs), smat))

#=
See comment above regarding the use of Meshes vs Surfaces or Polysurfaces!

KhepriBase.b_cuboid(b::RH, pb0, pb1, pb2, pb3, pt0, pt1, pt2, pt3, mat) =
  @remote(b, Mesh([pb0, pb1, pb2, pb3, pt0, pt1, pt2, pt3], 
                  [[0, 1, 2, 3],
                   [0, 1, 5, 4],
                   [1, 2, 6, 5],
                   [2, 3, 7, 6],
                   [3, 0, 4, 7],
                   [4, 5, 6, 7]],
                  mat))
=# 

KhepriBase.b_cuboid(b::RH, pb0, pb1, pb2, pb3, pt0, pt1, pt2, pt3, mat) =
  b_pyramid_frustum(b, [pb0, pb1, pb2, pb3], [pt0, pt1, pt2, pt3], mat, mat, mat)

KhepriBase.b_box(b::RH, c, dx, dy, dz, mat) =
  @remote(b, Box(c, vx(1, c.cs), vy(1, c.cs), dx, dy, dz, mat))

KhepriBase.b_sphere(b::RH, c, r, mat) =
	@remote(b, Sphere(c, r, mat))

KhepriBase.b_cone(b::RH, cb, r, h, bmat, smat) =
  @remote(b, Cone(add_z(cb, h), r, cb, smat))

KhepriBase.b_cone_frustum(b::RH, cb, rb, h, rt, bmat, tmat, smat) =
  @remote(b, ConeFrustum(cb, rb, cb + vz(h, cb.cs), rt, smat))

KhepriBase.b_text(b::RH, str, c, size, mat) =
  @remote(b, Text(str, c, size, "ISOCP", false, false))

## Materials
# To encode materials, we use a Guid but convert it to an int,
# taking into consideration that Guid(int) = int + 1
encode(ns::Val{:RH}, t::Val{:MatId}, c::IO, v) =
  encode(ns, Val(:int), c, v%Int32 - 1)
decode(ns::Val{:RH}, t::Val{:MatId}, c::IO) =
  (decode(ns, Val(:int), c) + 1)%Guid

KhepriBase.b_get_material(b::RH, path::String) =
  @remote(b, LoadRenderMaterialFromPath(path))

struct RhinoDefaultMaterial
  name::String
end

export rhino_default_material
rhino_default_material = RhinoDefaultMaterial

KhepriBase.b_get_material(b::RH, mat::RhinoDefaultMaterial) =
  b_get_material(b, joinpath(@remote(b, DefaultMaterialsFolder()), mat.name*".rmtl"))

KhepriBase.b_new_material(b::RH, name, color, specularity, roughness, transmissivity, transmitted_specular) =
  Guid(@remote(b, new_material(path, convert(RGBA, color), specularity, roughness)) + 1)

#=
KhepriBase.b_new_material(b::RH, name,
						  base_color,
						  metallic, specular, roughness,
	                 	  clearcoat, clearcoat_roughness,
						  ior,
						  transmission, transmission_roughness,
	                 	  emission_color,
						  emission_strength) =
  @remote(b, new_material(name,
  						  convert(RGBA, base_color),
						  metallic, specular, roughness,
  						  clearcoat, clearcoat_roughness,
  				  		  ior,
  				  		  transmission, transmission_roughness,
						  convert(RGBA, emission_color), emission_strength))

KhepriBase.b_plastic_material(b::RH, name, color, roughness) =
  @remote(b, new_material(name, convert(RGBA, color), 0.0, 1.0, roughness, 0.0, 0.0, 1.4, 0.0, 0.0, RGBA(0.0, 0.0, 0.0, 1.0), 0.0))

KhepriBase.b_metal_material(b::RH, name, color, roughness, ior) =
  @remote(b, new_metal_material(name, convert(RGBA, color), roughness, ior))

KhepriBase.b_glass_material(b::RH, name, color, roughness, ior) =
  @remote(b, new_glass_material(name, convert(RGBA, color), roughness, ior))

KhepriBase.b_mirror_material(b::RH, name, color) =
  @remote(b, new_mirror_material(name, convert(RGBA, color)))
=#

backend_map_division(b::RH, f::Function, s::Shape1D, n::Int) =
  let conn = connection(b),
      r = ref(b, s).value,
      (t1, t2) = @remote(b, CurveDomain(r)),
      ti = division(t1, t2, n),
      ps = @remote(b, CurvePointsAt(r, ti)),
      ts = @remote(b, CurveTangentsAt(r, ti)),
      #ns = ACADCurveNormalsAt(conn, r, ti),
      frames = rotation_minimizing_frames(@remote(b, CurveFrameAt(r, t1)), ps, ts)
    map(f, frames)
  end

#=
realize(b::RH, s::Text) =
  @remote(b, Text(s.str, s.c, vx(1, s.c.cs), vy(1, s.c.cs), s.h))
=#

backend_surface_domain(b::RH, s::Shape2D) =
    tuple(@remote(b, SurfaceDomain(ref(b, s).value)...))

backend_map_division(b::RH, f::Function, s::Shape2D, nu::Int, nv::Int) =
    let conn = connection(b)
        r = ref(b, s).value
        (u1, u2, v1, v2) = @remote(b, SurfaceDomain(r))
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(@remote(b, SurfaceFrameAt(r, u, v)))
            end
        end
    end

KhepriBase.b_torus(b::RH, c, ra, rb, mat) =
  @remote(b, Torus(c, vz(1, c.cs), ra, rb, mat))

KhepriBase.b_extruded_curve(b::RH, s::Shape1D, v, cb, mat) =
  and_mark_deleted(b,
      map_ref(b, s) do r
          @remote(b, Extrusion(r, v))
      end,
      s)

KhepriBase.b_extruded_surface(b::RH, s::Shape2D, v, cb, mat) =
   and_mark_deleted(b,
      map_ref(b, s) do r
          @remote(b, Extrusion(r, v))
      end,
      s)

#
#KhepriBase.b_sweep_curve(b::RH, path::Path, profile::Path, rotation, scale, mat) =

rhino_sweep(b, path, profile, rotation, scale, mat) =
  and_mark_deleted(b,
    map_ref(b, profile) do profile_r
        map_ref(b, path) do path_r
          @remote(b, SweepPathProfile(path_r, profile_r, rotation, scale))
        end
      end,
    [path, profile])

KhepriBase.b_swept_curve(b::RH, path::Shape1D, profile::Shape1D, rotation, scale, mat) =
  rhino_sweep(b, path, profile, rotation, scale, mat)
KhepriBase.b_swept_surface(b::RH, path::Shape1D, profile::Shape2D, rotation, scale, mat) =
  rhino_sweep(b, path, profile, rotation, scale, mat)

KhepriBase.b_subtract_ref(b::RH, s, r) =
  @remote(b, Subtract(s, r))

KhepriBase.b_intersect_ref(b::RH, s, r) =
  @remote(b, Intersect(s, r))

KhepriBase.b_unite_refs(b::RH, rs) =
  @remote(b, Unite(rs))

KhepriBase.b_slice_ref(b::RH, r, p, v) =
  @remote(b, Slice(r, p, v))

#=
realize(b::RH, s::Revolve) =
  and_delete_shape(
    map_ref(b, s.profile) do r
      @remote(b, Revolve(r, s.p, s.n, s.start_angle, s.amplitude))
    end,
    s.profile)

backend_loft_points(b::Backend, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  let f = (ruled ? (closed ? polygon : line) : (closed ? closed_spline : spline))
    and_delete_shapes(ref(b, f(map(point_position, profiles), backend=b)),
                      vcat(profiles, rails))
  end

backend_loft_curves(b::RH, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  and_delete_shapes(RHLoft(connection(b),
                             collect_ref(b, profiles),
                             collect_ref(b, rails),
                             ruled, closed),
                    vcat(profiles, rails))

backend_loft_surfaces(b::RH, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  backend_loft_curves(b, profiles, rails, ruled, closed)

backend_loft_curve_point(b::RH, profile::Shape, point::Shape) =
  and_delete_shapes(RHLoft(connection(b),
                             vcat(collect_ref(b, profile), collect_ref(b, point)),
                             [],
                             true, false),
                    [profile, point])

backend_loft_surface_point(b::RH, profile::Shape, point::Shape) =
  backend_loft_curve_point(b, profile, point)
=#

#=
slice_ref(b::RH, r::RHNativeRef, p::Loc, v::Vec) =
  (@remote(b, Slice(r.value, p, v); r))

=#
realize(b::RH, s::Move) =
  and_mark_deleted(b,
	  map_ref(b, s.shape) do r
      @remote(b, Move(r, s.v))
      r
    end,
    s.shape)

realize(b::RH, s::Scale) =
  and_mark_deleted(b,
    map_ref(b, s.shape) do r
      @remote(b, Scale(r, s.p, s.s))
      r
    end,
	  s.shape)

realize(b::RH, s::Rotate) =
  and_mark_deleted(b,
  	map_ref(b, s.shape) do r
      @remote(b, Rotate(r, s.p, s.v, s.angle))
      r
	  end,
	  s.shape)

realize(b::RH, s::Mirror) =
  and_mark_deleted(b, 
    map_ref(b, s.shape) do r
      @remote(b, Mirror(r, s.p, s.n, false))
    end,
    s.shape)
  
realize(b::RH, s::Thicken) =
  and_mark_deleted(b,
    map_ref(b, s.shape) do r
      @remote(b, Thicken(r, s.thickness))
    end,
    s.shape)
#=
backend_frame_at(b::RH, s::Shape2D, u::Real, v::Real) =
  @remote(b, SurfaceFrameAt(ref(b, s).value, u, v))

=#

# BIM

# Families

realize(b::RH, f::TableFamily) =
    @remote(b, CreateRectangularTableFamily(f.length, f.width, f.height, f.top_thickness, f.leg_thickness))
realize(b::RH, f::ChairFamily) =
    @remote(b, CreateChairFamily(f.length, f.width, f.height, f.seat_height, f.thickness))
realize(b::RH, f::TableChairFamily) =
    @remote(b, CreateRectangularTableAndChairsFamily(
        realize(b, f.table_family), realize(b, f.chair_family),
        f.table_family.length, f.table_family.width,
        f.chairs_top, f.chairs_bottom, f.chairs_right, f.chairs_left,
        f.spacing))

KhepriBase.b_table(b::RH, c, angle, family) =
    @remote(b, Table(c, angle, realize(b, family)))

KhepriBase.b_chair(b::RH, c, angle, family) =
    @remote(b, Chair(c, angle, realize(b, family)))

KhepriBase.b_table_and_chairs(b::RH, c, angle, family) =
    @remote(b, TableAndChairs(c, angle, realize(b, family)))

# Lights
KhepriBase.b_pointlight(b::RH, loc::Loc, color::RGB, intensity::Real, range::Real) =
  @remote(b, PointLight(loc, color, intensity))



# Dimensioning

no_props = Dict()
base_props = Dict("TextHeight"=>Float64(0.2))
base_props = Dict("DimensionScale"=>Float64(0.2))
label_props = merge(base_props)
angular_props = merge(base_props)
diametric_props = merge(base_props)

KhepriBase.b_labels(b::RH, p, strs, mats, mat) =
  [with(current_layer, mat.layer) do
    @remote(b, CreateLeaderDimension(str, p, p+vpol(annotation_scale(), ϕ), annotation_scale(), mark, label_props))
   end
   for (str, ϕ, mark, mat)
     in zip(strs,
            division(-π/4, 7π/4, length(strs), false),
            Iterators.flatten((Iterators.repeated("_DOT", 1), Iterators.repeated("_NONE"))),
            mats)]

#
KhepriBase.b_radii_illustration(b::RH, c, rs, rs_txts, mats, mat) =
  [with(current_layer, mat.layer) do
    @remote(b, CreateDiametricDimension(r_txt, c, c+vpol(r, ϕ), c+vpol(0.1, ϕ+pi/2), annotation_scale(), "", diametric_props))
   end
   for (r, r_txt, ϕ, mat) in zip(rs, rs_txts, division(π/6, 2π+π/6, length(rs), false), mats)]

# Maybe merge the texts when the radii are the same.
KhepriBase.b_vectors_illustration(b::RH, p, a, rs, rs_txts, mats, mat) =
  [with(current_layer, mat.layer) do
    @remote(b, CreateDiametricDimension(r_txt, p, p+vpol(r, a), p+vpol(0.1, a+pi/2), annotation_scale(), "", diametric_props))
   end
   for (r, r_txt, mat) in zip(rs, rs_txts, mats)]

KhepriBase.b_angles_illustration(b::RH, c, rs, ss, as, r_txts, s_txts, a_txts, mats, mat) =
  let refs = new_refs(b),
      maxr = maximum(rs),
      n = length(rs),
      ars = division(0.2maxr, 0.7maxr, n, false),
      idxs = sortperm(as),
      (rs, ss, as, r_txts, s_txts, a_txts) = (rs[idxs], ss[idxs], as[idxs], r_txts[idxs], s_txts[idxs], a_txts[idxs])
    for (r, ar, s, a, r_txt, s_txt, a_txt, mat) in zip(rs, ars, ss, as, r_txts, s_txts, a_txts, mats)
      with(current_layer, mat.layer) do
        if !(r ≈ 0.0)
          if !(s ≈ 0.0)
          #  push!(refs, @remote(b, CreateAngularDimension(s_txt, c, ar, 0, s, annotation_scale(), 5, no_props)))
            push!(refs, @remote(b, CreateAngularDimension(s_txt, c, c+vpol(0.8*ar, 0), c+vpol(0.8*ar, s), c+vpol(ar, s/2), annotation_scale(), "", angular_props)))
          end
          if !(a ≈ 0.0)
          #  push!(refs, @remote(b, CreateAngularDimension(a_txt, c, ar, s, s + a, annotation_scale(), 5, no_props)))
          push!(refs, @remote(b, CreateAngularDimension(a_txt, c, c+vpol(0.5*ar, s), c+vpol(0.5*ar, s + a), c+vpol(ar, s+a/2), annotation_scale(), "", angular_props)))
          else
            #push!(refs, @remote(b, CreateDiametricDimension(r_txt, c, c+vpol(maxr, s + a), 0.0, annotation_scale(), "")))
          end
        end
        #push!(refs, @remote(b, CreateDiametricDimension(r_txt, c, c+vpol(maxr, s + a), 0.0, annotation_scale(), "", diametric_props)))
      end
    end
    refs
  end

KhepriBase.b_arcs_illustration(b::RH, c, rs, ss, as, r_txts, s_txts, a_txts, mats, mat) =
  let refs = new_refs(b),
      maxr = maximum(rs),
      n = length(rs),
      ars = division(0.2maxr, 0.7maxr, n, false),
      idxs = sortperm(ss),
      (rs, ss, as, r_txts, s_txts, a_txts) = (rs[idxs], ss[idxs], as[idxs], r_txts[idxs], s_txts[idxs], a_txts[idxs])
    for (i, r, ar, s, a, r_txt, s_txt, a_txt, mat) in zip(1:n, rs, ars, ss, as, r_txts, s_txts, a_txts, mats)
      with(current_layer, mat.layer) do
        if !(r ≈ 0.0)
          if !(s ≈ 0.0) && ((i == 1) || !(s ≈ ss[i-1] + as[i-1]))
            push!(refs, @remote(b, CreateAngularDimension(s_txt, c, c+vpol(0.8*ar, 0), c+vpol(0.8*ar, s), c+vpol(ar, s/2), annotation_scale(), "", angular_props)))
          end
          if !(a ≈ 0.0)
            push!(refs, @remote(b, CreateAngularDimension(a_txt, c, c+vpol(0.8*ar, s), c+vpol(0.8*ar, s + a), c+vpol(ar, s+a/2), annotation_scale(), "", angular_props)))
          end
          #push!(refs, @remote(b, CreateDiametricDimension(r_txt, c, c+vpol(maxr, s + a), 0.0, annotation_scale(), "", diametric_props)))
        end
      end
    end
    refs
  end

#=
KhepriBase.b_spotlight(b::ACAD, loc::Loc, dir::Vec, hotspot::Real, falloff::Real) =
    @remote(b, SpotLight(loc, hotspot, falloff, loc + dir))

KhepriBase.b_ieslight(b::ACAD, file::String, loc::Loc, dir::Vec, alpha::Real, beta::Real, gamma::Real) =
    @remote(b, IESLight(file, loc, loc + dir, vxyz(alpha, beta, gamma)))
=#

############################################

# KhepriBase.b_bounding_box(b::RH, shapes::Shapes) =
#   @remote(b, BoundingBox(collect_ref(b, shapes)))

KhepriBase.b_set_view(b::RH, camera, target, lens, aperture) =
  @remote(b, View(camera, target, lens))

KhepriBase.b_get_view(b::RH) =
  @remote(b, ViewCamera()), @remote(b, ViewTarget()), @remote(b, ViewLens())

KhepriBase.b_zoom_extents(b::RH) =
  @remote(b, ZoomExtents())

KhepriBase.b_set_view_top(b::RH) =
  @remote(b, ViewTop())

KhepriBase.b_delete_refs(b::RH, refs::Vector{RHId}) =
  @remote(b, DeleteMany(refs))

KhepriBase.b_all_refs(b::RH) =
  @remote(b, GetAllShapes())

KhepriBase.b_delete_all_refs(b::RH) =
  @remote(b, DeleteAll())

KhepriBase.b_highlight_ref(b::RH, r::RHId) =
  @remote(b, Highlight(r))

KhepriBase.b_unhighlight_ref(b::RH, r::RHId) =
  @remote(b, Unhighlight(r))

KhepriBase.b_unhighlight_all_refs(b::RH) =
  @remote(b, UnhighlightAll())

# Layers
KhepriBase.b_current_layer(b::RH) =
  @remote(b, CurrentLayer())

KhepriBase.b_current_layer(b::RH, layer) =
  @remote(b, SetCurrentLayer(layer))

KhepriBase.b_layer(b::RH, name, active, color) =
  @remote(b, CreateLayer(name, true, color))

KhepriBase.b_delete_all_shapes_in_layer(b::RH, layer) =
  @remote(b, DeleteAllInLayer(layer))

KhepriBase.b_shape_from_ref(b::RH, r) =
  let code = @remote(b, ShapeCode(r)),
      ref = DynRefs(b=>RHNativeRef(r))
    if code == 1
      point(@remote(b, PointPosition(r)), ref=ref)
    elseif code == 2
      circle(loc_from_o_vz(@remote(b, CircleCenter(r)), @remote(b, CircleNormal(r))),
             @remote(b, CircleRadius(r)),
             ref=ref)
    elseif 3 <= code <= 6
      line(@remote(b, LineVertices(r)), ref=ref)
    elseif code == 7
      spline([xy(0,0)], false, false, #HACK obtain interpolation points
             ref=ref)
    elseif 103 <= code <= 106
      polygon(@remote(b, LineVertices(r)), ref=ref)
    elseif code == 40
      # FIXME: frontier is missing
      surface(frontier=[], ref=ref)
    elseif code == 41
      surface(frontier=[], ref=ref)
    elseif code == 81
      sphere(ref=ref)
    elseif code == 82
      cylinder(ref=ref)
    elseif code == 83
      cone(ref=ref)
    elseif code == 84
      torus(ref=ref)
    else
      error("Unknown shape")
    end
  end

#

KhepriBase.b_select_position(b::RH, prompt::String) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = @remote(b, GetPosition(prompt))
      length(ans) > 0 ? ans[1] : nothing
    end
  end

KhepriBase.b_select_positions(b::RH, prompt::String) =
  let ps = Loc[]
      p = nothing
    @info "$(prompt) on the $(b) backend."
    while ((p = @remote(b, GetPosition(prompt)) |> length)) > 0
        push!(ps, p...)
    end
    ps
  end

# HACK: The next operations should receive a set of shapes to avoid re-creating already existing shapes

KhepriBase.b_select_point(b::RH, prompt::String) =
  select_one_with_prompt(prompt, b, @get_remote b GetPoint)
KhepriBase.b_select_points(b::RH, prompt::String) =
  select_many_with_prompt(prompt, b, @get_remote b GetPoints)

KhepriBase.b_select_curve(b::RH, prompt::String) =
  select_one_with_prompt(prompt, b, @get_remote b GetCurve)
KhepriBase.b_select_curves(b::RH, prompt::String) =
  select_many_with_prompt(prompt, b, @get_remote b GetCurves)

KhepriBase.b_select_surface(b::RH, prompt::String) =
  select_one_with_prompt(prompt, b, @get_remote b GetSurface)
KhepriBase.b_select_surfaces(b::RH, prompt::String) =
  select_many_with_prompt(prompt, b, @get_remote b GetSurfaces)

KhepriBase.b_select_solid(b::RH, prompt::String) =
  select_one_with_prompt(prompt, b, @get_remote b GetSolid)
KhepriBase.b_select_solids(b::RH, prompt::String) =
  select_many_with_prompt(prompt, b, @get_remote b GetSolids)

KhepriBase.b_select_shape(b::RH, prompt::String) =
  select_one_with_prompt(prompt, b, @get_remote b GetShape)
KhepriBase.b_select_shapes(b::RH, prompt::String) =
  select_many_with_prompt(prompt, b, @get_remote b GetShapes)

backend_captured_shape(b::RH, handle) =
  backend_shape_from_ref(b, handle)
#
backend_captured_shapes(b::RH, handles) =
  map(handles) do handle
      backend_shape_from_ref(b, handle)
  end

backend_generate_captured_shape(b::RH, s::Shape) =
    println("captured_shape(rhino, $(ref(b, s).value))")
backend_generate_captured_shapes(b::RH, ss::Shapes) =
  begin
    print("captured_shapes(rhino, [")
    for s in ss
      print(ref(b, s).value)
      print(", ")
    end
    println("])")
  end


backend_all_shapes_in_layer(b::RH, layer) =
  Shape[backend_shape_from_ref(b, r) for r in @remote(b, GetAllShapesInLayer(layer))]

# Khepri render quality ranges from -1 to 1
# Rhino render quality ranges from 0 to 3
KhepriBase.b_render_and_save_view(b::RH, path::String) =
  let kind = render_kind()
    if kind == :realistic
      let quality = round(Int, (render_quality() + 1)*3/2)
        @remote(b, Render(render_width(), render_height(), quality, path))
      end
    else
      @remote(b, RenderLoadEnvironment("KhepriStudio", joinpath(@__DIR__, "KhepriStudio.renv")))
      let (camera, target) = (@remote(b, ViewCamera()), @remote(b, ViewTarget())),
          rot = 11π/6 + π/2 - (camera - target).ϕ
        if render_kind() == :white
          @remote(b, ClayRenderWhite("KhepriStudio", rot, render_width(), render_height(), 3, path))
        elseif render_kind() == :black
          @remote(b, ClayRenderBlack("KhepriStudio", rot, render_width(), render_height(), 3, path))
        else
          error("Unknown render kind $kind")
        end
      end
    end
  end

#=
BIM families for Rhino
=#

abstract type RhinoFamily <: Family end

struct RhinoLayerFamily <: RhinoFamily
  name::String
  color::RGB
  ref::Parameter{Any}
end

rhino_layer_family(name, color::RGB=rgb(1,1,1)) =
  RhinoLayerFamily(name, color, Parameter{Any}(nothing))

backend_get_family_ref(b::RH, f::Family, af::RhinoLayerFamily) =
  backend_create_layer(b, af.name, true, af.color)

KhepriBase.b_realistic_sky(b::RH, date, latitude, longitude, elevation, meridian, turbidity, sun) =
  @remote(b, SunLight(date, latitude, longitude, meridian, turbidity, true))

export run_grasshopper_solver
run_grasshopper_solver() =
  @remote(rhino, RunGrasshopperSolver())
