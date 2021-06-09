"""
code taken from:
https://gist.github.com/lvngd/54a26a748c073d35e269f90419c0f629

Computes the Convex Hull with the Graham Scan algorithm

Use:
	h = ConvexHull(points)
	print(h.hull)
"""
class ConvexHull:
    def __init__(self, points):
        self.points = points
        self.hull, self.hull_ids = self.compute_convex_hull()
        #self.hull_ids = []

    def get_cross_product(self,p1, p2, p3):
        return ((p2[0] - p1[0])*(p3[1] - p1[1])) - ((p2[1] - p1[1])*(p3[0] - p1[0]))

    def get_slope(self,p1, p2):
        if p1[0] == p2[0]:
            return float('inf')
        else:
            return 1.0*(p1[1]-p2[1])/(p1[0]-p2[0])

    def compute_convex_hull(self):
        hull = []
        hull_ids = []

        self.points.sort(key=lambda x:[x[0],x[1]])
        start = self.points.pop(0)

        hull.append(start)
        hull_ids.append(start[2])

        self.points.sort(key=lambda p: (self.get_slope(p,start), -p[1],p[0]))

        for pt in self.points:
            hull.append(pt)
            hull_ids.append(pt[2])

            while len(hull) > 2 and self.get_cross_product(hull[-3],hull[-2],hull[-1]) < 0:
                hull.pop(-2)
                hull_ids.pop(-2)

        return hull, hull_ids