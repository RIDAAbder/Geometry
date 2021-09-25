using Autodesk.Revit.DB;
using Autodesk.Revit.UI;
using RevitUtils.Models.Geometry;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RevitUtils.Logic
{
    public static class Geometry
    {
        const double _offset = 0.1;

        public static List<XYZ> GetPoints(this Face face)
        {
            var list = new List<XYZ>();
            foreach (EdgeArray array in face.EdgeLoops)
            {
                foreach (Edge edge in array)
                {
                    list.Add(edge.AsCurve().GetEndPoint(0));
                    list.Add(edge.AsCurve().GetEndPoint(1));
                }
            }
            return list.Distinct().ToList();
        }
        public static LineXYZ ToSlice(this Curve curve)
        {
            return new LineXYZ(curve.GetEndPoint(0).ToPointXYZ(), curve.GetEndPoint(1).ToPointXYZ());
        }
        public static LineXY ToLineXY(this Curve curve)
        {
            return new LineXY(curve.GetEndPoint(0).ToPointXY(), curve.GetEndPoint(1).ToPointXY());
        }
        public static Curve ToCurve2D(this Curve curve)
        {
            return Line.CreateBound(curve.GetEndPoint(0).ToPoint2D(), curve.GetEndPoint(1).ToPoint2D());
        }
        public static BoundingBoxXYZ GetBounds(this Curve curve)
        {
            return curve.Tessellate().ToList().GetBounds();
        }
        public static BoundingBoxXYZ Offset(this BoundingBoxXYZ bounds, double offset)
        {
            return new BoundingBoxXYZ()
            {
                Min = new XYZ(bounds.Min.X - offset, bounds.Min.Y - offset, bounds.Min.Z - offset),
                Max = new XYZ(bounds.Max.X + offset, bounds.Max.Y + offset, bounds.Max.Z + offset),
            };
        }
        public static bool Intersects(this BoundingBoxXYZ bounds1, BoundingBoxXYZ bounds2)
        {
            if ((bounds1.Max.X < bounds2.Min.X || bounds2.Max.X < bounds1.Min.X)
                && (bounds1.Max.Y < bounds2.Min.Y || bounds2.Max.Y < bounds1.Min.Y)
                && (bounds1.Max.Z < bounds2.Min.Z || bounds2.Max.Z < bounds1.Min.Z))
                return false;
            return true;
        }
        public static double GetLength2D(this Curve curve)
        {
            return curve.GetEndPoint(0).ToPoint2D().DistanceTo(curve.GetEndPoint(1).ToPoint2D());
        }
        public static List<List<XYZ>> GetFloorHoles(Floor floor)
        {
            List<List<XYZ>> polygons = new List<List<XYZ>>();

            var geo = floor.get_Geometry(new Options());
            foreach (GeometryObject obj in geo)
            {
                Solid solid = obj as Solid;
                if (solid != null)
                {
                    GetBoundary(polygons, solid);
                }
            }

            //remove the outer polygon
            var diagonals = polygons.Select(o => o.GetDiagonalLength()).ToList();
            var max = diagonals.Max();
            var index = diagonals.FindIndex(o => o == max);
            polygons.RemoveAt(index);

            return polygons;
        }
        public static bool GetBoundary(List<List<XYZ>> polygons, Solid solid)
        {
            PlanarFace lowest = null;
            FaceArray faces = solid.Faces;
            foreach (Face f in faces)
            {
                PlanarFace pf = f as PlanarFace;
                if (null != pf && pf.IsHorizontal())
                {
                    if ((null == lowest) || (pf.Origin.Z < lowest.Origin.Z))
                    {
                        lowest = pf;
                    }
                }
            }
            if (null != lowest)
            {
                XYZ p, q = XYZ.Zero;
                bool first;
                int i, n;
                EdgeArrayArray loops = lowest.EdgeLoops;
                foreach (EdgeArray loop in loops)
                {
                    List<XYZ> vertices = new List<XYZ>();
                    first = true;
                    foreach (Edge e in loop)
                    {
                        var points = e.Tessellate();
                        p = points.ElementAt(0);
                        if (!first)
                        {
                            Debug.Assert(p.IsAlmostEqualTo(q), "expected subsequent start point" + " to equal previous end point");
                        }
                        n = points.Count;
                        q = points.ElementAt(n - 1);
                        for (i = 0; i < n - 1; ++i)
                        {
                            XYZ v = points.ElementAt(i);
                            v = new XYZ(v.X, v.Y, v.Z - _offset);
                            vertices.Add(v);
                        }
                    }
                    q = new XYZ(q.X, q.Y, q.Z - _offset);
                    Debug.Assert(q.IsAlmostEqualTo(vertices[0]), "expected last end point to equal" + " first start point");
                    polygons.Add(vertices);
                }
            }
            return null != lowest;
        }
        public static bool Contains(this Curve line, XYZ pt)
        {
            var totalLength = line.Length;
            var len1 = pt.DistanceTo(line.GetEndPoint(0));
            var len2 = pt.DistanceTo(line.GetEndPoint(1));
            if (Math.Abs(totalLength - (len1 + len2)) > 0.001)
                return false;
            return true;
        }
        public static bool IsAlmostEqual(this double n1, double n2)
        {
            return n1.IsAlmostEqual(n2, 0.0001);
        }
        public static bool IsAlmostEqual(this double n1, double n2,double tolerance)
        {
            return Math.Abs(n1 - n2) < tolerance;
        }
        public static bool IsVerticalFace(this Face f)
        {
            var normal = f.ComputeNormal(new UV()).Normalize();
            if (normal.Z.IsAlmostEqual(0))
                return true;
            return false;
        }
        public static bool IsHorizontalFace(this Face f)
        {
            var normal = f.ComputeNormal(new UV()).Normalize();
            if (normal.X.IsAlmostEqual(0) && normal.Y.IsAlmostEqual(0))
                return true;
            return false;
        }
        public static bool IsHorizontal(this PlanarFace face)
        {
            double firstZ = 0;
            for (int i = 0; i < face.EdgeLoops.Size; i++)
            {
                var loop = face.EdgeLoops.get_Item(i);
                for (int j = 0; j < loop.Size; j++)
                {
                    if (i == 0 && j == 0)
                        firstZ = loop.get_Item(j).AsCurve().GetEndPoint(0).Z;

                    var curve = loop.get_Item(j).AsCurve();
                    if (curve.GetEndPoint(0).Z != firstZ || curve.GetEndPoint(1).Z != firstZ)
                        return false;
                }
            }
            return true;
        }
        public static bool IsHorizontal3D(this Curve curve)
        {
            return curve.GetEndPoint(0).Z == curve.GetEndPoint(1).Z;
        }
        public static bool IsVertical3D(this Curve curve)
        {
            return curve.GetEndPoint(0).ToPoint2D().IsAlmostEqualTo(curve.GetEndPoint(1).ToPoint2D());
        }
        public static bool IsHorizontal2D(this Curve curve)
        {
            return curve.GetEndPoint(0).Y == curve.GetEndPoint(1).Y;
        }
        public static bool IsVertical2D(this Curve curve)
        {
            return curve.GetEndPoint(0).X == curve.GetEndPoint(1).X;
        }
        public static void MovePoint(this Curve line, XYZ pt1, XYZ vec)
        {
            var newPosition = pt1 + vec;
            if (line.GetEndPoint(0).IsAlmostEqualTo(pt1))
            {
                var par1 = pt1.DistanceTo(newPosition);
                var par2 = line.Length;
                line.MakeBound(par1, par2);
            }
            else
            {
                var par2 = line.GetEndPoint(0).DistanceTo(newPosition);
                line.MakeBound(0, par2);
            }
        }
        public static XYZ GetCentroid(this BoundingBoxXYZ bBox)
        {
            var max = bBox.Max;
            var min = bBox.Min;
            var vec = max - min;
            return min + (0.5 * vec);
        }
        public static double GetDiagonalLength(this List<XYZ> points)
        {
            return PolyLine.Create(points).GetOutline().GetDiagonalLength();
        }
        public static BoundingBoxXYZ GetBounds(this List<XYZ> points)
        {
            var minX = points.Min(o => o.X);
            var minY = points.Min(o => o.Y);
            var minZ = points.Min(o => o.Z);

            var maxX = points.Max(o => o.X);
            var maxY = points.Max(o => o.Y);
            var maxZ = points.Max(o => o.Z);

            return new BoundingBoxXYZ() { Min = new XYZ(minX, minY, minZ), Max = new XYZ(maxX, maxY, maxZ) };
        }
        public static List<Curve> GetBounds2D(this List<XYZ> points)
        {
            var minX = points.Min(o => o.X);
            var minY = points.Min(o => o.Y);
            var maxX = points.Max(o => o.X);
            var maxY = points.Max(o => o.Y);
            var l1 = Line.CreateBound(new XYZ(minX, minY, 0), new XYZ(minX, maxY, 0));
            var l2 = Line.CreateBound(new XYZ(minX, maxY, 0), new XYZ(maxX, maxY, 0));
            var l3 = Line.CreateBound(new XYZ(maxX, maxY, 0), new XYZ(maxX, minY, 0));
            var l4 = Line.CreateBound(new XYZ(maxX, minY, 0), new XYZ(minX, minY, 0));

            return new List<Curve>() { l1, l2, l3, l4 };
        }
        public static List<Curve> GetBounds(this BoundingBoxXYZ box)
        {
            var minPt = box.Min;
            var maxPt = box.Max;
            var minX = minPt.X;
            var minY = minPt.Y;
            var maxX = maxPt.X;
            var maxY = maxPt.Y;
            var l1 = Line.CreateBound(new XYZ(minX, minY, 0), new XYZ(minX, maxY, 0));
            var l2 = Line.CreateBound(new XYZ(minX, maxY, 0), new XYZ(maxX, maxY, 0));
            var l3 = Line.CreateBound(new XYZ(maxX, maxY, 0), new XYZ(maxX, minY, 0));
            var l4 = Line.CreateBound(new XYZ(maxX, minY, 0), new XYZ(minX, minY, 0));

            return new List<Curve>() { l1, l2, l3, l4 };
        }
        public static bool RectangleContain2D(this List<Curve> rectangle, XYZ pt)
        {
            var allPoints = rectangle.SelectMany(o => new List<XYZ>() { o.GetEndPoint(0), o.GetEndPoint(1) }).ToList();
            if (allPoints.Max(o => o.X) > pt.X && allPoints.Min(o => o.X) < pt.X
                    && allPoints.Max(o => o.Y) > pt.Y && allPoints.Min(o => o.Y) < pt.Y)
                return true;
            return false;
        }
        public static bool Intersects2D(this Curve l1, Curve l2)
        {
            var l1Temp = Line.CreateBound(new XYZ(l1.GetEndPoint(0).X, l1.GetEndPoint(0).Y, 0), new XYZ(l1.GetEndPoint(1).X, l1.GetEndPoint(1).Y, 0));
            var l2Temp = Line.CreateBound(new XYZ(l2.GetEndPoint(0).X, l2.GetEndPoint(0).Y, 0), new XYZ(l2.GetEndPoint(1).X, l2.GetEndPoint(1).Y, 0));
            IntersectionResultArray intersection = null;
            l1Temp.Intersect(l2Temp, out intersection);
            if (intersection != null)
                return true;
            return false;
        }
        public static XYZ GetIntersection2D(this Curve l1, Curve l2)
        {
            var l1Temp = Line.CreateBound(new XYZ(l1.GetEndPoint(0).X, l1.GetEndPoint(0).Y, 0), new XYZ(l1.GetEndPoint(1).X, l1.GetEndPoint(1).Y, 0));
            var l2Temp = Line.CreateBound(new XYZ(l2.GetEndPoint(0).X, l2.GetEndPoint(0).Y, 0), new XYZ(l2.GetEndPoint(1).X, l2.GetEndPoint(1).Y, 0));
            IntersectionResultArray intersection = null;
            l1Temp.Intersect(l2Temp, out intersection);
            if (intersection != null)
                return intersection.get_Item(0).XYZPoint;
            return null;
        }

        public static XYZ ToPoint2D(this XYZ pt)
        {
            return new XYZ(pt.X, pt.Y, 0);
        }
        public static XYZ GetNormalVector2D(this XYZ pt)
        {
            return new XYZ(-pt.Y, pt.X, 0);
        }

        public static PointXYZ ToPointXYZ(this XYZ pt)
        {
            return new PointXYZ(pt.X, pt.Y, pt.Z);
        }
        public static PointXY ToPointXY(this XYZ pt)
        {
            return new PointXY(pt.X, pt.Y);
        }
        public static XYZ GetMidPoint(this Curve curve)
        {
            double param1 = curve.GetEndParameter(0);
            double param2 = curve.GetEndParameter(1);

            double paramCalc = param1 + ((param2 - param1) * 0.5); //* requiredDist / curveLength);

            XYZ evaluatedPoint = null;

            if (curve.IsInside(paramCalc))
            {
                double normParam = curve
                  .ComputeNormalizedParameter(paramCalc);

                evaluatedPoint = curve.Evaluate(normParam, true);
            }
            return evaluatedPoint;
        }
        public static XYZ GetProjection(this XYZ p, XYZ start, XYZ end)
        {
            var v1 = new XYZ(p.X - start.X, p.Y - start.Y, 0);
            var v2 = new XYZ(end.X - start.X, end.Y - start.Y, 0);
            var l2 = v2.X * v2.X + v2.Y * v2.Y;
            var dotproduct = v1.X * v2.X + v1.Y * v2.Y;
            var ndistance = dotproduct / l2;
            return new XYZ(start.X + v2.X * ndistance, start.Y + v2.Y * ndistance, 0);
        }

        public static List<XYZ> GetCurvePoints(this Curve curve)
        {
            return new List<XYZ>() { curve.GetEndPoint(0), curve.GetEndPoint(1), curve.GetMidPoint() };
        }
        public static Curve GetCurve(this Element element)
        {
            return (element.Location as LocationCurve).Curve;
        }
        public static XYZ GetVector(this Curve l)
        {
            return l.GetEndPoint(1) - l.GetEndPoint(0);
        }
        public static bool IsParallel(this Curve l1, Curve l2)
        {
            var vec1 = l1.GetVector().Normalize();
            var vec2 = l2.GetVector().Normalize();
            if (vec1.IsAlmostEqualTo(vec2) || vec1.IsAlmostEqualTo(-vec2))
                return true;
            return false;
        }
        public static List<Curve> GetEdges(this Face f)
        {
            var list = new List<Curve>();
            foreach (EdgeArray array in f.EdgeLoops)
                foreach (Edge edge in array)
                    list.Add(edge.AsCurve());
            return list;
        }
        public static RectangleXY GetBounds(this List<LineXY> lines)
        {
            var minX = lines.Min(o => Math.Min(o.Start.X, o.End.X));
            var maxX = lines.Max(o => Math.Max(o.Start.X, o.End.X));
            var width = Math.Abs(maxX - minX);
            var minY = lines.Min(o => Math.Min(o.Start.Y, o.End.Y));
            var maxY = lines.Max(o => Math.Max(o.Start.Y, o.End.Y));
            var height = Math.Abs(maxY - minY);
            return new RectangleXY(minX, minY, width, height);
        }
        public static Face GetMinFace(this List<Face> faces)
        {
            var hozFaces = faces.Where(o => o.ComputeNormal(new UV()).Normalize().IsAlmostEqualTo(XYZ.BasisZ) || o.ComputeNormal(new UV()).Normalize().IsAlmostEqualTo(-XYZ.BasisZ)).ToList();
            if (hozFaces.Count == 0) return null;
            return hozFaces.OrderBy(o => o.EdgeLoops.get_Item(0).get_Item(0).AsCurve().GetEndPoint(0).Z).ElementAt(0);
        }
        public static bool IsApproxEqual(this double num1, double num2)
        {
            return Math.Abs(num1 - num2) < 0.0001;
        }
        public static void PaintElemnets(this List<ElementId> elements, Document doc, Color edgesColor, Color fillColor)
        {
            OverrideGraphicSettings ogs = new OverrideGraphicSettings();
            ogs.SetProjectionLineColor(edgesColor);
            Element solidFill = new FilteredElementCollector(doc).OfClass(typeof(FillPatternElement)).Where(q => q.Name.Contains("Solid")).FirstOrDefault();

            if (solidFill != null)
            {
                if (doc.Application.VersionNumber != "2020")
                {
                    SetFill(fillColor, ogs, solidFill);
                }
                foreach (var item in elements)
                {
                    try
                    {
                        doc.ActiveView.SetElementOverrides(item, ogs);
                    }
                    catch (Exception)
                    {
                    }
                }
            }
        }

        private static void SetFill(Color fillColor, OverrideGraphicSettings ogs, Element solidFill)
        {
            ogs.SetProjectionFillPatternId(solidFill.Id);
            ogs.SetProjectionFillColor(fillColor);
        }

        public static void ChangeTransparency(this List<ElementId> elements, Document doc, int transparent)
        {
            OverrideGraphicSettings ogs = new OverrideGraphicSettings();
            ogs.SetSurfaceTransparency(transparent);
            foreach (var item in elements)
            {
                try
                {
                    doc.ActiveView.SetElementOverrides(item, ogs);
                }
                catch (Exception)
                {
                }
            }
        }
        public static Curve MovePointTo(this Curve curve, XYZ pt, XYZ newPosition)
        {
            if (curve.GetEndPoint(0).IsAlmostEqualTo(pt))
                return Line.CreateBound(newPosition, curve.GetEndPoint(1));
            else
                return Line.CreateBound(curve.GetEndPoint(0), newPosition);
        }
        public static double GetHeight(this BoundingBoxXYZ bBox)
        {
            return bBox.Max.Z - bBox.Min.Z;
        }
        public static BoundingBox ToCustomBBox(this BoundingBoxXYZ bBox)
        {
            return new BoundingBox(bBox);
        }
        public static List<Curve> GetShownCurves(this Element elem, Options op)
        {
            var geom = elem.get_Geometry(op);
            var edges = new List<Curve>();
            if (geom == null)
                return edges;
            foreach (GeometryObject obj in geom)
            {
                if (obj is GeometryInstance)
                {
                    var geomInst = obj as GeometryInstance;
                    var geom2 = geomInst.GetInstanceGeometry();
                    foreach (GeometryObject obj2 in geom2)
                    {
                        if (obj2 is Curve)
                            edges.Add(obj2 as Curve);
                    }
                }
                else if (obj is Curve)
                    edges.Add(obj as Curve);
            }

            return edges;
        }

        /// <summary>
        /// Check the intersection between a plane and line
        /// https://rosettacode.org/wiki/Find_the_intersection_of_a_line_with_a_plane#C.23
        /// </summary>
        /// <param name="rayVector"></param>
        /// <param name="rayPoint"></param>
        /// <param name="planeNormal"></param>
        /// <param name="planePoint"></param>
        /// <returns></returns>
        public static XYZ IntersectPoint(XYZ rayVector, XYZ rayPoint, XYZ planeNormal, XYZ planePoint)
        {
            var diff = rayPoint - planePoint;
            var prod1 = diff.DotProduct(planeNormal);
            var prod2 = rayVector.DotProduct(planeNormal);
            var prod3 = prod1 / prod2;
            return rayPoint - rayVector * prod3;
        }
        public static bool IsIntersect(this Face f1, Face f2)
        {
            Curve intersection = null;
            f1.Intersect(f2, out intersection);
            if (intersection == null)
                return false;
            return true;
        }

        public static XYZ GetUnitVector(this Curve line)
        {
            return (line.GetEndPoint(1) - line.GetEndPoint(0)).Normalize();
        }
        public static BoundingBoxXYZ GetBoundingBoxXYZ(this List<XYZ> points)
        {
            var minX = points.Min(o => o.X);
            var minY = points.Min(o => o.Y);
            var minZ = points.Min(o => o.Z);
            var maxX = points.Max(o => o.X);
            var maxY = points.Max(o => o.Y);
            var maxZ = points.Max(o => o.Z);
            return new BoundingBoxXYZ() { Min = new XYZ(minX, minY, minZ), Max = new XYZ(maxX, maxY, maxZ) };
        }
        /// <summary>
        /// check that the point M is inside the triangle ABC
        /// </summary>
        /// <param name="m"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public bool CheckPointMInsideTraingle(Point m, Point a, Point b, Point c)
        {
            /
            double condVal = -tolerance * 400;
            double mX = m.cordonnes.X;
            double mY = m.cordonnes.Y;
            double aX = a.cordonnes.X;
            double aY = a.cordonnes.Y;
            double bX = b.cordonnes.X;
            double bY = b.cordonnes.Y;
            double cX = c.cordonnes.X;
            double cY = c.cordonnes.Y;
            //Test1 det((AB) , (AM)) x det((AM) , (AC)) > 0
            //Test2 det((BA) ,(BM)) x det((BM) ,(BC)) > 0 [pemutation A'<-B /  B'<-A ]
            //Test3 det((CA) , (CM)) x det((CM) , (CB)) > 0 [pemutation B'<-C / C<-B' ]
            double t1 = (((aX - bX) * (aY - mY)) - ((aY - bY) * (aX - mX))) * (((aX - mX) * (aY - cY)) - ((aY - mY) * (aX - cX)));
            double t2 = (((bX - aX) * (bY - mY)) - ((bY - aY) * (bX - mX))) * (((bX - mX) * (bY - cY)) - ((bY - mY) * (bX - cX)));
            double t3 = (((cX - aX) * (cY - mY)) - ((cY - aY) * (cX - mX))) * (((cX - mX) * (cY - bY)) - ((cY - mY) * (cX - bX)));

            return (t1 >= condVal && t2 >= condVal && t3 >= condVal) || (t1 <= condVal && t2 <= condVal && t3 <= condVal);
        }
        /// <summary>
        /// Get the list of points of a room
        /// </summary>
        /// <param name="room"></param>
        /// <returns></returns>
        public List<XYZ> GetRoomPoints(Room room)
        {
            var points = new List<XYZ>();
            var xyzcomparator = new XyzComparator();
            var opt = new SpatialElementBoundaryOptions();
            opt.SpatialElementBoundaryLocation = SpatialElementBoundaryLocation.Center;
            var loops = room.GetBoundarySegments(opt);
            foreach (var seg in loops[0])
            {
                var curve = seg.GetCurve();
                XYZ p;
                XYZ q = null;
                p = curve.GetEndPoint(0);
                q = curve.GetEndPoint(1);
                if (!points.Contains(p, xyzcomparator)) points.Add(p);
                if (!points.Contains(q, xyzcomparator)) points.Add(q);
            }
            return points.Distinct(xyzcomparator).ToList();
        }

    }
}
