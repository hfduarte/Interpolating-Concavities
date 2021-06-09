import numpy as np
import math
import sys
import time

import shapely.wkt

from shapely.geometry import Point, Polygon, LineString, MultiPolygon
from shapely.geometry import GeometryCollection, mapping

from shapely.wkt import dumps, loads

from MSeg import MSeg
from MSeg2 import MSeg2

from convex_hull import ConvexHull
from skgeom import *

OPT = 0
TESTING = False
NTESTS = 100
Option = 1
type_transf = 1
chull_compute_time = 0
epsilon = 0.000001
PLOT_RES = False

"""
import matplotlib.pyplot as plt

t1 = float(5)
n1 = 1000
n2 = 2000
n3 = 3000
n4 = 10000

t2 = (n2 / n1)**2 * t1
t3 = (n3 / n1)**2 * t1
t4 = (n4 / n1)**2 * t1

plt.plot([n1, n2, n3, n4], [t1, t2, t3, t4], 'ro')
plt.plot([n1, n2, n3, n4], [t1, t2, t3, t4], 'r--')

plt.show()
print(t1, t2, t3, t4)
sys.exit()
"""

if PLOT_RES:
    import matplotlib.pyplot as plt

    t1 = float(0.001)
    n1 = 4
    n2 = 16
    n3 = 256

    t2 = (n2 / n1)**2 * t1
    t3 = (n3 / n1)**2 * t1

    plt.plot([ 4, 16, 256], [0.001, 0.005, 0.102], 'ko')
    plt.plot([ 4, 16, 256], [0.001, 0.005, 0.102], 'k--')

    plt.plot([ 4, 16, 256], [t1, t2, t3], 'ro')
    plt.plot([ 4, 16, 256], [t1, t2, t3], 'r--')

    plt.ylabel('execution time (s)')
    plt.xlabel('number of points (n)')
    plt.show()
    
    print(t1, t2, t3)
    sys.exit()

"""
#import numpy as np

# 100 linearly spaced numbers
x = np.linspace(-5,5,100)

# the function, which is y = x^2 here
y = x**2

# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# plot the function
plt.plot(x,y, 'r')

# show the plot
plt.show()
"""

"""

    Seg >> Concavity OV : N Tests: 1000
    mET (sec);MET (sec);AVGET (sec);NV;NC;NS;CHT;
    0.001;0.002;0.001;4;1;1;9.298324584960938e-06;
    0.004;0.007;0.005;16;1;1;2.86102294921875e-05;
    0.094;0.230;0.102;256;1;1;0.00044465065002441406;

    5, 20.4

    python lcip_opt.py 'LineString (-0.01196172248803817 -0.36602870813397126, -1.05023923444976064 -0.24880382775119614)' 'LineString (-1.01196172248803817 -0.07416267942583721, -0.96650717703349276 0.03349282296650724, -0.88516746411483238 0.11004784688995217, -0.90430622009569372 0.15311004784688997, -0.97129186602870798 0.11004784688995217, -0.96650717703349276 0.07655502392344504, -1.00717703349282273 0.05502392344497609, -1.03588516746411474 0.09569377990430628, -1.08133971291866016 0.09090909090909094, -1.11961722488038262 0.16028708133971292, -1.02631578947368407 0.17942583732057416, -0.92344497607655485 0.22009569377990434, -0.83492822966507174 0.21531100478468901, -0.78468899521531088 0.16507177033492826, -0.6866028708133971 0.18181818181818188, -0.76315789473684204 0.22248803827751196, -0.78229665071770316 0.30622009569377995, -0.8875598086124401 0.32296650717703357, -0.97607655502392343 0.33014354066985652, -1.06220095693779903 0.27033492822966509, -1.16746411483253576 0.22488038277511968, -1.25119617224880386 0.23684210526315796, -1.24880382775119614 0.16985645933014359, -1.25358851674641136 0.10287081339712922, -1.22248803827751185 0.01435406698564601, -1.31100478468899517 0.03827751196172258, -1.3157894736842104 0.10287081339712922, -1.34688995215310992 0.16507177033492826, -1.36363636363636354 0.21052631578947367, -1.32535885167464107 0.30143540669856461, -1.20574162679425823 0.30861244019138756, -1.10765550239234445 0.31818181818181823, -1.09330143540669855 0.38277511961722488, -1.22009569377990434 0.40669856459330145, -1.16507177033492804 0.46650717703349287, -1.09090909090909083 0.45215311004784686, -1.11004784688995217 0.51913875598086134, -1.04066985645932997 0.50717703349282295, -1.02392344497607635 0.44736842105263164, -0.94976076555023914 0.43540669856459335, -0.94497607655502391 0.53349282296650724, -0.98564593301435388 0.59569377990430628, -1.07177033492822948 0.61722488038277512, -1.16267942583732031 0.59330143540669855, -1.24880382775119614 0.59808612440191389, -1.31339712918660267 0.63397129186602874, -1.33253588516746402 0.71531100478468901, -1.24401913875598069 0.78229665071770338, -1.11961722488038262 0.8133971291866029, -0.97607655502392343 0.83014354066985652, -0.89952153110047828 0.77511961722488043, -0.93540669856459302 0.70334928229665072, -0.86363636363636354 0.65789473684210531, -0.85406698564593286 0.49282296650717705, -0.80622009569377973 0.42583732057416268, -0.83014354066985629 0.53110047846889952, -0.74880382775119614 0.52631578947368429, -0.70574162679425823 0.58851674641148333, -0.66985645933014348 0.66267942583732053, -0.6004784688995215 0.69617224880382778, -0.64832535885167464 0.74641148325358853, -0.71770334928229662 0.75598086124401909, -0.72727272727272729 0.82775119617224879, -0.64832535885167464 0.87320574162679421, -0.55023923444976064 0.87799043062200954, -0.48325358851674638 0.80622009569377995, -0.46172248803827753 0.74401913875598091, -0.5023923444976075 0.67703349282296654, -0.50717703349282295 0.60765550239234445, -0.57655502392344493 0.5861244019138756, -0.61244019138755967 0.53110047846889952, -0.53349282296650702 0.48086124401913866, -0.45933014354066981 0.43301435406698563, -0.38995215311004783 0.4856459330143541, -0.37559808612440193 0.57655502392344493, -0.30382775119617222 0.60287081339712922, -0.24641148325358841 0.55263157894736836, -0.17942583732057416 0.47368421052631582, -0.13397129186602874 0.38038277511961727, -0.12918660287081329 0.299043062200957, -0.11004784688995217 0.19377990430622016, -0.21052631578947367 0.16746411483253587, -0.23444976076555024 0.22248803827751196, -0.3014354066985645 0.24162679425837319, -0.36124401913875581 0.24641148325358853, -0.34449760765550241 0.3133971291866029, -0.27990430622009566 0.30622009569377995, -0.22727272727272729 0.29665071770334928, -0.24880382775119614 0.35406698564593309, -0.24641148325358841 0.4138755980861244, -0.29186602870813383 0.48803827751196172, -0.33971291866028697 0.40909090909090906, -0.42583732057416257 0.32057416267942584, -0.5023923444976075 0.28708133971291872, -0.48325358851674638 0.16507177033492826, -0.42583732057416257 0.21052631578947367, -0.36124401913875581 0.18421052631578949, -0.2870813397129186 0.09330143540669855, -0.1866028708133971 0.09330143540669855, -0.17703349282296643 0.03110047846889952, -0.22248803827751185 -0.02870813397129179, -0.06220095693779903 0.04545454545454553, 0.02392344497607679 0.10765550239234456, -0.00478468899521522 0.18181818181818188, 0.08133971291866038 0.2129186602870814, 0.03110047846889952 0.27751196172248804, 0.0837320574162681 0.35167464114832536, -0.02392344497607657 0.37081339712918659, -0.07177033492822971 0.4282296650717704, -0.09569377990430605 0.48086124401913877, -0.13397129186602874 0.55023923444976086, -0.16507177033492826 0.62200956937799046, -0.22009569377990434 0.68181818181818188, -0.21052631578947367 0.74162679425837319, -0.13157894736842102 0.76794258373205748, -0.07894736842105265 0.73923444976076558, -0.09090909090909083 0.69138755980861244, -0.07416267942583721 0.62918660287081341, -0.05023923444976064 0.56937799043062198, -0.02392344497607657 0.50717703349282295, 0.0143540669856459 0.45454545454545459, 0.07177033492822993 0.43779904306220097, 0.09330143540669855 0.48803827751196172, 0.07177033492822993 0.55023923444976086, 0.05741626794258359 0.64832535885167464, 0.10765550239234445 0.72248803827751196, 0.18899521531100483 0.75119617224880386, 0.30622009569378017 0.77272727272727271, 0.33971291866028697 0.69377990430622005, 0.3133971291866029 0.65550239234449759, 0.24401913875598069 0.70574162679425845, 0.17703349282296665 0.67703349282296654, 0.13157894736842124 0.60526315789473684, 0.15311004784688986 0.52392344497607657, 0.16507177033492848 0.45933014354066992, 0.15789473684210531 0.38038277511961727, 0.23205741626794252 0.49043062200956944, 0.18899521531100483 0.55023923444976086, 0.20334928229665072 0.62440191387559807, 0.26794258373205748 0.62200956937799046, 0.32296650717703335 0.4856459330143541, 0.299043062200957 0.36363636363636365, 0.39952153110047872 0.28229665071770338, 0.42822966507177052 0.17464114832535893, 0.31578947368421062 0.0861244019138756, 0.22488038277511979 0.05023923444976075, 0.11961722488038262 0.04545454545454553, 0.02631578947368407 -0.01674641148325362, -0.06459330143540654 -0.04784688995215314, -0.1004784688995215 -0.12440191387559807)' '0' '100' '0' '1' 
    python lcip_opt.py 'LineString (-0.01196172248803817 -0.36602870813397126, -1.05023923444976064 -0.24880382775119614)' 'LineString (-1.01196172248803817 -0.07416267942583721, -0.6866028708133971 0.18181818181818188, -0.36124401913875581 0.18421052631578949, -0.1004784688995215 -0.12440191387559807)' '0' '100' '0' '1' 


    An optimized version

    Changes:
        convex hull is now computed using the Graham scan algorithm:
        
        convex_hull_2 >> ConvexHull(points)
        visibility >> O(n)
        time is computed in a different way

    Examples:
    
    python lcip_opt.py 'LINESTRING(-140.50924455177292316 -54.26986651320897437, -140.06526900166338123 -54.25620572705175704)' 'LINESTRING(-140.38800507462761402 -54.16228782222089677, -140.33848472480769942 -54.09569148970446406, -140.39312786943656874 -54.0393407468059479, -140.34360751961665414 -54.02738755891838451, -140.30774795595397109 -54.06495472085072862, -140.20016926496589349 -54.05812432777211995, -140.16089450476388834 -54.12130546374925189, -140.10795895840468006 -54.08715349835620856, -140.07039179647230753 -54.10252188278307273, -140.13698812898874735 -54.15033463433333338, -140.07039179647233595 -54.18960939453533143)' 0 100 0 1 
    
    python lcip_opt.py 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)' 'LineString (-1.1775510204081634 -0.00408163265306127, -0.81428571428571428 0.19183673469387752, -0.61020408163265305 -0.00816326530612255, -0.42653061224489797 0.15102040816326523, -0.6020408163265305 0.32244897959183672, -0.5122448979591836 0.42040816326530606, -0.19387755102040805 0.16326530612244894, -0.083673469387755 -0.05714285714285716)' '0' '100' '0' '1' 
    
    Tests:
    
    python lcip_opt.py '10'

    Normal Transformation:
    
    python lcip_opt.py 'LineString (-1.06326530612244907 -0.46530612244897962, 0.92448979591836755 -0.42448979591836755)' 'LineString (-1.08775510204081627 -0.07755102040816331, -0.76938775510204072 -0.13469387755102047, -0.71224489795918355 0.00408163265306116, -0.61020408163265305 0.13061224489795908, -0.42653061224489797 0.22448979591836726, -0.23061224489795906 0.19999999999999996, -0.0428571428571427 0.08571428571428563, 0.0346938775510206 -0.03673469387755102, -0.01836734693877551 -0.08571428571428585, -0.07551020408163245 -0.00816326530612255, -0.15714285714285703 0.04081632653061218, -0.21428571428571419 0, -0.23877551020408161 0.02040816326530603, -0.17755102040816317 0.08979591836734691, -0.34897959183673466 0.13061224489795908, -0.45510204081632644 0.08979591836734691, -0.57755102040816331 -0.02040816326530615, -0.45510204081632644 -0.03673469387755102, -0.393877551020408 -0.08979591836734713, -0.40612244897959182 -0.16734693877551021, -0.44693877551020411 -0.13877551020408174, -0.44693877551020411 -0.106122448979592, -0.53673469387755102 -0.07346938775510203, -0.63877551020408152 -0.07346938775510203, -0.63877551020408152 -0.14693877551020429, -0.2795918367346939 -0.23673469387755119, 0.14489795918367365 -0.1959183673469389, 0.2551020408163267 0.00816326530612232, 0.1530612244897962 0.14693877551020396, -0.23469387755102034 0.39183673469387748, -0.47551020408163258 0.38775510204081631, -0.72448979591836737 0.27346938775510199, -0.80204081632653068 0.05306122448979589, -0.86326530612244889 0.05306122448979589, -0.83469387755102042 0.19999999999999996, -0.71632653061224483 0.39183673469387748, -0.49591836734693873 0.48571428571428565, -0.23469387755102034 0.49387755102040809, -0.15714285714285703 0.61632653061224485, -0.02244897959183678 0.63265306122448983, 0.11632653061224518 0.55918367346938769, 0.13265306122448983 0.4408163265306122, 0.06326530612244907 0.4408163265306122, 0.0591836734693878 0.51020408163265296, -0.03469387755102016 0.55102040816326525, -0.07142857142857117 0.4408163265306122, 0.07959183673469417 0.32244897959183672, 0.30000000000000027 0.18775510204081625, 0.47142857142857153 0.06530612244897949, 0.4918367346938779 -0.03265306122448997, 0.40612244897959204 -0.05714285714285716, 0.40204081632653077 0.04897959183673462, 0.35306122448979593 0.04897959183673462, 0.34489795918367339 -0.10204081632653073, 0.47142857142857153 -0.22448979591836737, 0.81836734693877577 -0.08163265306122458)' 0 100 0  1

    python lcip_opt.py 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)' 'LineString (-0.73497709287796686 1.13821740941274596, 0.28299875052061729 1.63827571845064712, 0.5062390670553949 1.29002082465639445, 0.24728029987505318 1.25430237401082989, 0.24728029987505318 1.42396501457726066, -0.0027488546438974 1.29002082465639445, 0.29192836318200843 1.04892128279883501, 0.72054977092878092 1.07571012078300843, 0.63125364431486997 1.45075385256143408, 0.20263223656809748 1.83472719700125109, 0.90807163681799397 1.71864223240316694, 0.88128279883382055 1.34359850062474084, 1.25632653061224664 1.1024989587671814, 2.05999167013744477 1.3793169512703054)' '0' '100' '0' '2' 

    python lcip_opt.py 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)' 'LineString (-0.73497709287796686 1.13821740941274596, 0.28299875052061729 1.63827571845064712, 0.5062390670553949 1.29002082465639445, 0.24728029987505318 1.25430237401082989, 0.24728029987505318 1.42396501457726066, -0.0027488546438974 1.29002082465639445, 0.29192836318200843 1.04892128279883501, 0.72054977092878092 1.07571012078300843, 0.63125364431486997 1.45075385256143408, 0.20263223656809748 1.83472719700125109, 0.90807163681799397 1.71864223240316694, 0.88128279883382055 1.34359850062474084, 1.25632653061224664 1.1024989587671814, 2.05999167013744477 1.3793169512703054)' '0' '100' '0' 

"""

# --------------------------------------------------------------------------------------------------------------------------
# Other Tests
# --------------------------------------------------------------------------------------------------------------------------

def get_input():
    global TESTING
    global NTESTS

    if len(sys.argv) == 2:
        NTESTS = int(str(sys.argv[1]))
        TESTING = True
        return None, None, 0, 0, False, 1

    p_wkt = str(sys.argv[1])
    q_wkt = str(sys.argv[2])
    option = int(str(sys.argv[3]))
    n_obs = int(str(sys.argv[4]))
    debug = int(str(sys.argv[5]))
    type_transf = int(str(sys.argv[6]))

    p = loads(p_wkt)
    q = loads(q_wkt)

    return p, q, option, n_obs, debug, type_transf

# --------------------------------------------------------------------------------------------------------------------------
# Auxiliary Functions
# --------------------------------------------------------------------------------------------------------------------------

def get_transform4(i, j, keys, points, nt, line_coords):
    global epsilon
    viz = []
    isegs = []
    #_isegs = []
    msegs = []
    
    ML = nt
    
    # concavity is a 'triangle'
    if j - i == 2:
        viz = [i, i + 1, j]        
    # concavity is a 'polygon'
    else:
        conc = []
        _v = j

        # 1 point visibility test

        p0 = points[i]
        p1 = points[i+1]
        p2 = points[j]

        conc += [Segment2(p0, p2)]

        _p0 = None
        _p1 = None

        while _v > i:
            _p0 = points[_v]
            _v -= 1
            _p1 = points[_v]

            conc += [Segment2(_p0, _p1)]

        arr = arrangement.Arrangement()

        for s in conc:
            arr.insert(s)

        he = arr.halfedges
        edge = None

        e = next(he)
        edge = e.prev()

        if edge.source().point() != p1 or edge.target().point() != p0:
            edge = e.next()
            edge = edge.twin()

            if edge.source().point() != p1 or edge.target().point() != p0:
                edge = None

        if edge == None:
            print(i, j, 'NULL_EDGE')
            sys.exit()

        ts = TriangularExpansionVisibility(arr)
        tx = ts.compute_visibility(p0, edge)

        # 2 point visibility test

        vpoints = []

        for v in tx.vertices:
            vpoints += [Point2(v.point().x(), v.point().y())]

        conc = []

        _i = 0
        n = len(vpoints) - 1

        while _i < n :
            conc += [Segment2(vpoints[_i], vpoints[_i+1])]
            _i += 1

        conc += [Segment2(vpoints[n], vpoints[0])]

        arr = arrangement.Arrangement()

        for s in conc:
            arr.insert(s)

        p0 = vpoints[0]
        p1 = vpoints[1]

        he = arr.halfedges
        edge = None
        
        # always the first or the second halfedge
        for v in he:
            if v.source().point() == p0 and v.target().point() == p1:
                edge = v
                break

        ts = TriangularExpansionVisibility(arr)
        tx = ts.compute_visibility(p1, edge)        

        # get only visible points that are vertices on the concavity
    
        for v in tx.vertices:
            _key = float('%.9f' % v.point().x())
    
            k = keys.get(_key)
            if k is not None:
                if len(k) == 1:
                    """
                    _id = k[0][1]
                    if _id != i and _id != j:
                        viz = [_id] + viz
                    """
                    viz = [k[0][1]] + viz
                else:
                    for key in k:
                        if v.point().y() - epsilon < key[0] and v.point().y() + epsilon > key[0]:
                            """
                            _id = key[1]
                            if _id != i and _id != j:
                                viz = [_id] + viz
                            """
                            viz = [key[1]] + viz
                            break

    # get the transformations for the sub-concavities if any

    #i = iseg[0][0]
    #j = iseg[m - 1][1]
    
    nsegs = len(viz) - 1

    SA = Point(line_coords[i][0], line_coords[i][1])
    SB = Point(line_coords[j][0], line_coords[j][1])
        
    f = float(1) / (nsegs)
    
    #segs = []
    
    dx = (SB.x - SA.x)
    dy = (SB.y - SA.y)
    
    xstep = dx * f
    ystep = dy * f
    
    #h = 0
    #tup = (SA.x, SA.y)
    
    """
    while h < nsegs:
        if h == m - 1:
            segs += [[tup, (B.x, B.y)]]
        else:
            x = A.x + xstep * (h+1)
            y = A.y + ystep * (h+1)
            
            segs += [[tup, (x, y)]]
            tup = (x, y)
        
        h += 1
    """

    #print(viz, 'viz')

    nnt = nt + 1
    
    # visible points form a trriangle
    if nsegs == 2:
            A = (SA.x, SA.y)
            B = (SA.x + xstep, SA.y + ystep)
            
            D = (line_coords[viz[1]][0], line_coords[viz[1]][1])
                        
            msegs += [MSeg2(B, B, D, D, None, None, nt)]
            
            if viz[0] + 1 == viz[1]:
                msegs += [MSeg2(D, D, D, D, None, 1, nt)]

                #print([viz[0], viz[1], nt])
                if viz[1] + 1 != viz[2]:
                    #print([viz[1], viz[2], nt])
                    _msegs, ml = get_transform4(viz[1], viz[2], keys, points, nnt, line_coords)
                    msegs += _msegs

                    if ML < ml:
                        ML = ml
            else:
                msegs += [MSeg2(A, A, A, A, None, 1, nt)]
                #print([viz[0], viz[1], nt])
                _msegs, ml = get_transform4(viz[0], viz[1], keys, points, nnt, line_coords)
                msegs += _msegs

                if ML < ml:
                    ML = ml
                
                msegs += [MSeg2(D, D, D, D, None, 1, nt)]
                                
                if viz[1] + 1 != viz[2]:
                    #print([viz[1], viz[2], nt])
                    _msegs, ml = get_transform4(viz[1], viz[2], keys, points, nnt, line_coords)
                    msegs += _msegs

                    if ML < ml:
                        ML = ml

            D = (line_coords[viz[2]][0], line_coords[viz[2]][1])
            msegs += [MSeg2(D, D, D, D, None, 1, nt)]    
    else:
        g = 0
        M = nsegs - 1
        while g < nsegs:
                _next = g + 1
                #print([viz[g], viz[_next], nt])
                C = (line_coords[viz[g]][0], line_coords[viz[g]][1])
                D = (line_coords[viz[_next]][0], line_coords[viz[_next]][1])

                if g == 0:
                    A = (SA.x, SA.y)
                    B = (SA.x + xstep, SA.y + ystep)

                    msegs += [MSeg2(B, B, C, D, None, None, nt)]

                    if viz[0] + 1 != viz[1]:
                        msegs += [MSeg2(A, A, A, A, None, 1, nt)]

                        _msegs, ml = get_transform4(viz[0], viz[1], keys, points, nnt, line_coords)
                        msegs += _msegs

                        if ML < ml:
                            ML = ml

                    msegs += [MSeg2(D, D, D, D, None, 1, nt)]
                elif g == M:
                    A = (SA.x + xstep * (M), SA.y + ystep * (M))
                    B = (SB.x, SB.y)

                    # final corner
                    msegs += [MSeg2(A, A, C, D, None, None, nt)]
                    # marl initial point
                    msegs += [MSeg2(C, C, C, C, None, 1, nt)]

                    if viz[g] + 1 != viz[_next]:
                        _msegs, ml = get_transform4(viz[g], viz[_next], keys, points, nnt, line_coords)
                        msegs += _msegs

                        if ML < ml:
                            ML = ml

                    # marl final point
                    msegs += [MSeg2(D, D, D, D, None, 1, nt)]
                else:
                    A = (SA.x + xstep * (g), SA.y + ystep * (g))
                    B = (SA.x + xstep * (_next), SA.y + ystep * (_next))

                    msegs += [MSeg2(A, A, C, D, None, None, nt)]
                    msegs += [MSeg2(A, B, D, D, None, None, nt)]
                    
                    msegs += [MSeg2(C, C, C, C, None, 1, nt)]
                    
                    if viz[g] + 1 != viz[_next]:
                        _msegs, ml = get_transform4(viz[g], viz[_next], keys, points, nnt, line_coords)
                        msegs += _msegs

                        if ML < ml:
                            ML = ml

                    msegs += [MSeg2(D, D, D, D, None, 1, nt)]
                
                g += 1

    """
    k = 0
    nvp = len(viz) - 1
    
    while k < nvp:
        if viz[k+1] - viz[k] > 1:

            isegs += [(viz[k], viz[k+1], nt)]
            print([[viz[k], viz[k+1], nt]])
                        
            __isegs, ml = get_transform4(viz[k], viz[k+1], keys, points, nnt, line_coords, e_t)
            #_isegs += __isegs
            isegs += __isegs
            
            if ML < ml:
                ML = ml
        else:
            isegs += [(viz[k], viz[k+1], nt)]
            print([[viz[k], viz[k+1], nt]])

        k += 1
    """
    
    """
    _i = i
    for e in viz:
        #isegs += [[_i, e]]
        if e - _i > 1:
            _nnt = nt + 1
            isegs += [[_i, e, nt]]
            print([[_i, e, nt]])
                        
            __isegs, __ml = get_transform4(_i, e, keys, points, _nnt)
            _isegs += __isegs
            
            if ML < __ml:
                ML = __ml
        else:
            isegs += [[_i, e, nt]]
            print([[_i, e, nt]])

        _i = e
    """

    """
    if j - _i > 1:
        _nnt = nt + 1
        isegs += [[_i, j, nt]]
        print( [[_i, j, nt]])
        
        __isegs, __ml = get_transform4(_i, j, keys, points, _nnt)
        _isegs += __isegs
        
        if ML < __ml:
            ML = __ml
    else:
        isegs += [[_i, j, nt]]
        print([[_i, j, nt]])
    """
    
    return msegs, ML

def get_isegs4(i, j, keys, points, nt, line_coords):
    return get_transform4(i, j, keys, points, nt, line_coords)

# --------------------------------------------------------------------------------------------------------------------------
# Computational Geometry
# --------------------------------------------------------------------------------------------------------------------------

"""
https://www.geeksforgeeks.org/convex-hull-set-1-jarviss-algorithm-or-wrapping/
"""
def orientation(p, q, r): 
    ''' 
    To find orientation of ordered triplet (p, q, r).  
    The function returns the following values:
        0 --> p, q and r are colinear  
        1 --> Clockwise  
        2 --> Counterclockwise  
    '''
    val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y) 
    
    #A, B, C
    #return (B.x - A.x) * (C.y - A.y) > (B.y - A.y) * (C.x - A.x)
  
    if val == 0: 
        return 0
    elif val > 0: 
        return 1
    else: 
        return 2

"""
    O(n2)
"""
def convex_hull_2(points, n): 
      
    # There must be at least 3 points  
    if n < 3: 
        return
  
    # Find the leftmost point 
    #l = 12
    
    # O(n), n the number of points.
    l = Left_index(points) 
    
    #print(l)
    
    hull = [] 
      
    ''' 
    Start from leftmost point, keep moving counterclockwise  
    until reach the start point again. This loop runs O(h)  
    times where h is number of points in result or output.  
    '''
    p = l 
    q = 0
    
    while(True): 
          
        # Add current point to result  
        hull.append(p) 
  
        ''' 
        Search for a point 'q' such that orientation(p, x,  
        q) is counterclockwise for all points 'x'. The idea  
        is to keep track of last visited most counterclock-  
        wise point in q. If any point 'i' is more counterclock-  
        wise than q, then update q.  
        '''
        q = (p + 1) % n 
  
        for i in range(n): 
              
            # If i is more counterclockwise than current q, then update q  
            if(orientation(points[p], points[i], points[q]) == 2): 
                q = i 
  
        ''' 
        Now q is the most counterclockwise with respect to p  
        Set p as q for next iteration, so that q is added to  
        result 'hull'  
        '''
        p = q 
  
        # While we don't come to first point 
        if(p == l): 
            break
    
    return hull

''' 
    Finding the left most point 
'''
def Left_index(points): 
    minn = 0
    
    for i in range(1, len(points)): 
        if points[i].x < points[minn].x: 
            minn = i 
        elif points[i].x == points[minn].x: 
            if points[i].y > points[minn].y: 
                minn = i
    
    return minn 

''' 
    Finding the right most point 
'''
def right_index(points): 
    maxn = 0
    
    for i in range(1, len(points)): 
        if points[i].x > points[maxn].x: 
            maxn = i 
        elif points[i].x == points[maxn].x: 
            if points[i].y < points[maxn].y: 
                maxn = i
    
    return maxn 

# --------------------------------------------------------------------------------------------------------------------------
# Other Tests
# --------------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------------------------------------------------------

def tests(seg_wkt, line_wkt, num_tests, op, face_wkt = None):
    k = 0
    
    min_exec_time = sys.float_info.max
    max_exec_time = sys.float_info.min
    
    precision3 = '.3f'
    precision = '.2f'
    
    t_avg = 0

    p = loads(seg_wkt)
    q = loads(line_wkt)

    if op == 2:
        f = loads(face_wkt)
    elif op == 3:
        objs = face_wkt

    NC = 0
    while k < num_tests:
        begin_exec_time = time.time()
        ch_ex_t = 0
        n_segs = 0
    
        if op == 1:
            msegs, n_segs, ch_ex_t, _n_concavities, start_time, step = get_seg_to_concavity_msegs_opt(p, q)
            NC = _n_concavities
            #n_msegs, n_segs, ch_ex_t = seg_to_concavity_testing(p, q)
        elif op == 2:
            #msegs = seg_to_concavity_with_face(p, q, f)
            
            # pre procc
        
            l_coords = q.coords
            p_coords = f.exterior.coords[:-1]

            line_coords = []
    
            l_id = 15
            o_id = 7
    
            dxx = 0.00001
    
            r = 0
            while r < len(l_coords):
                if r < l_id:
                    line_coords += [(l_coords[r][0], l_coords[r][1])]
                elif r == l_id:
                    line_coords += [(l_coords[r][0] - dxx, l_coords[r][1])]
            
                    a = o_id
                    while a >= 0:
                        if a == o_id:
                            line_coords += [(p_coords[a][0] - dxx, p_coords[a][1])]
                        else:
                            line_coords += [(p_coords[a][0], p_coords[a][1])]
                
                        a -= 1
            
                    a = len(p_coords) - 1

                    while a >= o_id:
                        if a == o_id:
                            line_coords += [(p_coords[a][0] + dxx, p_coords[a][1])]
                        else:
                            line_coords += [(p_coords[a][0], p_coords[a][1])]
                
                        a -= 1
            
                    line_coords += [(l_coords[r][0] + dxx, l_coords[r][1])]
                else:
                    line_coords += [(l_coords[r][0], l_coords[r][1])]
    
                r += 1

            nline = LineString(line_coords)
    
            # end pre procc
        
            msegs, nsegs, ch_et, _n_concavities, start_time, step = get_seg_to_concavity_msegs_opt(p, nline)            
        elif op == 3:
            #reg_to_reg_with_concavity(p, q, objs)
            msegs, nsegs, ch_et, _n_concavities, start_time, step = get_seg_to_concavity_msegs_opt(p, q)
        
            for obj in objs:
                a = obj[0]
                b = obj[1]
        
                A = (a.coords[0][0], a.coords[0][1])
                B = (a.coords[1][0], a.coords[1][1])
        
                C = (b.coords[0][0], b.coords[0][1])
                D = (b.coords[1][0], b.coords[1][1])
        
                msegs += [MSeg2(A, A, C, D, 0, 1, 0, True)]
                msegs += [MSeg2(A, B, D, D, 0, 1, 0, True)]
            
        end_exec_time = time.time()
        exec_time = end_exec_time - begin_exec_time

        min_exec_time = min(min_exec_time, exec_time)
        max_exec_time = max(max_exec_time, exec_time)
        
        t_avg += exec_time
        k += 1
    
    NV = len(loads(line_wkt).coords)
    t_avg = t_avg / num_tests
    
    return (format(min_exec_time, precision3), format(max_exec_time, precision3), format(t_avg, precision3), str(NV), str(NC), str(n_segs), str(ch_ex_t))

# --------------------------------------------------------------------------------------------------------------------------
# Transformation Algorithm.
# --------------------------------------------------------------------------------------------------------------------------

def get_seg_to_concavity_msegs_opt(seg, line):
    seg_coords = seg.coords
    d_seg = seg.length
    
    line_coords = line.coords
    
    line_min_dist = 0.2
    msegs = []
        
    d_max = 0
    
    ix = sys.float_info.max
    fx = sys.float_info.min

    id_x = -1
    fd_x = -1
    
    n = len(line_coords)
    k = 0
    
    s_t = 0
    e_t = 0.15
        
    # --------------------------------------------------------------------------------------------------------------------------
    # Line convex-hull points ids.
    # --------------------------------------------------------------------------------------------------------------------------
    
    keys = {}
    cpoints = []
    points = []
    J = 0
    for coord in line_coords:
        x = coord[0]
        y = coord[1]
    
        _key = float('%.9f' % x)

        k = keys.get(_key)
        if k is None:
            keys[_key] = [(y, J)]
        else:
            k += [(y, J)]
    
        #points.append((coord[0], coord[1], J)) 
        points.append((x, y, J)) 
        #if J > 0:
        #    cpoints = [Point2(x, y)] + cpoints
        cpoints.append(Point2(x, y))
        J += 1
    
    # Graham scan algorithm: O(nlogn)
    t_0 = time.time()
    H = ConvexHull(points)
    t_1 = time.time()

    chull_ids = H.hull_ids
    
    k = len(chull_ids)
    
    """
    print(chull_ids)
    print(k)
    sys.exit()
    """
        
    # --------------------------------------------------------------------------------------------------------------------------
    # Test visibility between seg endpointa and line endpoints.
    # --------------------------------------------------------------------------------------------------------------------------
    
    _ik = 0
    _fk = n - 1
    
    # s_0
    A = Point(seg_coords[0][0], seg_coords[0][1])
    B = Point(line_coords[_ik][0], line_coords[_ik][1])
        
    _ll = LineString([(A.x, A.y), (B.x, B.y)])
                
    if _ll.crosses(line):
        print_error('The initial endpoint of the target is not visible to the source!"')
        sys.exit()
    
    # s_1
    A = Point(seg_coords[1][0], seg_coords[1][1])
    B = Point(line_coords[_fk][0], line_coords[_fk][1])
        
    _ll = LineString([(A.x, A.y), (B.x, B.y)])
                
    if _ll.crosses(line):
        print_error('The final endpoint of the target is not visible to the source!"')
        sys.exit()
    
    ### End Visibility Test.

    # --------------------------------------------------------------------------------------------------------------------------
    # Choose the ids in the CHull between the initial and final endpoints of the line.
    # --------------------------------------------------------------------------------------------------------------------------

    new_chull_ids = []
    f = 0
    
    for _id in chull_ids:
        if _id == _ik:
            new_chull_ids += [_id]
            f = 1
        elif f == 1 and _id == _fk:
            new_chull_ids += [_id]
            f = 2
            break
        elif f == 1:
            new_chull_ids += [_id]
        
    _n_segs = 0
    
    if f != 2:
        chull_ids = [_ik, _fk]
        
        i_vis_idx = _ik
        f_vis_idx = _fk
        
        _n_segs = 1
    else:
        chull_ids = new_chull_ids
    
        i_vis_idx = _ik
        f_vis_idx = _fk
        
        _n_segs = len(chull_ids) - 1
    
    # --------------------------------------------------------------------------------------------------------------------------    
    
    k = len(chull_ids)
    
    """
    print(chull_ids)
    print(k)
    sys.exit()
    """
    
    i_last_idx = 0
    f_last_idx = len(chull_ids) - 1
    
    # --------------------------------------------------------------------------------------------------------------------------    
    # Initial seg > line interpolation (optional).
    # --------------------------------------------------------------------------------------------------------------------------        
    """
    factor = 0.25
    
    while factor > 0:
        id_x = i_vis_idx
        fd_x = f_vis_idx
        
        i_seg_x = line_coords[id_x][0]
        i_seg_y = line_coords[id_x][1]
            
        f_seg_x = line_coords[fd_x][0]
        f_seg_y = line_coords[fd_x][1]
            
        i_x = seg_coords[0][0] + factor * (i_seg_x - seg_coords[0][0]);
        i_y = seg_coords[0][1] + factor * (i_seg_y - seg_coords[0][1]);

        f_x = seg_coords[1][0] + factor * (f_seg_x - seg_coords[1][0]);
        f_y = seg_coords[1][1] + factor * (f_seg_y - seg_coords[1][1]);

        iseg = LineString([(i_x, i_y), (f_x, f_y)]);
            
        if iseg.distance(line) > line_min_dist:
            A = (seg_coords[0][0], seg_coords[0][1])
            B = (f_x, f_y)
                
            msegs += [MSeg2(A, A, (i_x, i_y), (f_x, f_y), s_t, e_t, 0)]
            msegs += [MSeg2((seg_coords[0][0], seg_coords[0][1]), (seg_coords[1][0], seg_coords[1][1]), B, B, s_t, e_t)]
                
            seg_coords = iseg.coords
            s_t = e_t
            break
            
        factor -= 0.1
    """
    # End principal line interpolation (optional).

    # --------------------------------------------------------------------------------------------------------------------------        
    # Divide the current segment in the transformation in n segments. O(k)    
    # -------------------------------------------------------------------------------------------------------------------------- 
    
    segs = []
    
    if _n_segs == 1:
        segs += [[(seg_coords[0][0],seg_coords[0][1]), (seg_coords[1][0], seg_coords[1][1])]]
    else:
        f = float(1) / _n_segs
        
        dx = (seg_coords[1][0] - seg_coords[0][0])
        dy = (seg_coords[1][1] - seg_coords[0][1])
        
        xstep = dx * f
        ystep = dy * f
        
        h = 0
        v = 1
        tup = (seg_coords[0][0], seg_coords[0][1])
        
        while h < k - 1:
        
            if h < i_last_idx or h >= f_last_idx:
                segs += [[]]
            else:
                x = seg_coords[0][0] + xstep * v
                y = seg_coords[0][1] + ystep * v
                    
                if h == k - 2 or h == f_last_idx - 1:
                    segs += [[(x, y), (seg_coords[1][0], seg_coords[1][1])]]
                else:
                    segs += [[tup, (x, y)]]
                    tup = (x, y)
                v += 1
            
            h += 1
    
    # --------------------------------------------------------------------------------------------------------------------------        
    # Get the transformations and the info to compute the time interval associated with each set of transformations.
    # --------------------------------------------------------------------------------------------------------------------------        
    
    e_t = 0.5

    _n_concavities = 0

    ML = 0
    i = 0
    while i < k - 1:
        j = i + 1
        
        # Seg > Concavity Transformation.        
        if chull_ids[i] + 1 != chull_ids[j]:
            A = (segs[i][0][0], segs[i][0][1])
            B = (segs[i][1][0], segs[i][1][1])
                
            C = (line_coords[chull_ids[i]][0], line_coords[chull_ids[i]][1])
            D = (line_coords[chull_ids[j]][0], line_coords[chull_ids[j]][1])

            msegs += [MSeg2(A, A, C, D, s_t, e_t, 0)]
            msegs += [MSeg2(A, B, D, D, s_t, e_t, 0)]

            #A = (line_coords[chull_ids[i]][0], line_coords[chull_ids[i]][1])
            #msegs += [MSeg2(A, A, A, A, e_t, 1, 0)]
            msegs += [MSeg2(C, C, C, C, e_t, 1, 0)]

             # Seg-concavity transformation.
             
            _msegs, ml = get_isegs4(chull_ids[i], chull_ids[j], keys, cpoints, 0, line_coords)
            msegs += _msegs

            if ML < ml:
                ML = ml

            if i == k - 2:
                msegs += [MSeg2(D, D, D, D, e_t, 1, 0)] 
            
            _n_concavities += 1
        # Direct Seg > Seg Transformation.
        else:
            A = (segs[i][0][0], segs[i][0][1])
            B = (segs[i][1][0], segs[i][1][1])
                
            _i = chull_ids[i]
            _j = chull_ids[j]
                
            C = (line_coords[_i][0], line_coords[_i][1])
            D = (line_coords[_j][0], line_coords[_j][1])
                
            msegs += [MSeg2(A, A, C, D, s_t, e_t, 0, True)]
            msegs += [MSeg2(A, B, D, D, s_t, e_t, 0, True)]

        i += 1
    
    start_time = e_t
    step = (1 - e_t) / (ML + 1)
    
    return msegs, _n_segs, (t_1 - t_0), _n_concavities, start_time, step

# --------------------------------------------------------------------------------------------------------------------------
# Print Functions
# --------------------------------------------------------------------------------------------------------------------------

def get_in_between_observations2(seg, line, msegs, num_samples, start_time, step, is_region = False, face = None):
   
    i = 0
    t = 0
    n = num_samples - 1
    num_invalid_geoms = 0
    num_complex_geoms = 0
    geoms = []
    
    while i < num_samples:
        t = float(i) / n

        if t == 0:
            geoms += [seg]
        elif t == 1:
            if face != None:
                s = GeometryCollection([line, face])
                geoms += [s]
            else:
                geoms += [line]
        else:
            coords = []
            M = 0

            for mseg in msegs:
                if i == 1:
                    xi, yi, xf, yf = mseg.at2(t, start_time, step)
                else:
                    xi, yi, xf, yf = mseg.at(t)
                
                if xi == None:
                    continue
                
                # All Points (No Filter).
                _n = len(coords)

                if _n > 1:
                    _xi = coords[_n - 2][0]
                    _yi = coords[_n - 2][1]
                        
                    _xf = coords[_n - 1][0]
                    _yf = coords[_n - 1][1]
                        
                    if _xi == xi and _yi == yi and _xf == xf and _yf == yf:
                        continue

                coords += [[xi, yi]]
                coords += [[xf, yf]]               
            
            g = LineString(coords)
            g = g.simplify(0.000000001)
            
            #
            _DX = 0.000000001
            _C = g.coords
            _N = len(_C)
            _I = 1
            _Coords = [(_C[0][0], _C[0][1])]
            while _I < _N:
                _X0 = _C[_I-1][0]
                _Y0 = _C[_I-1][1]
                
                _X1 = _C[_I][0]
                _Y1 = _C[_I][1]
                
                if _X1 == _X0 and _Y1 == _Y0:
                    pass
                #elif _X1 == _X0 and _Y0 - _DX <= _Y1 and _Y1 <= _Y0 + _DX:
                #    pass
                elif _X0 - _DX <= _X1 and _X1 <= _X0 + _DX and _Y0 - _DX <= _Y1 and _Y1 <= _Y0 + _DX:
                    pass
                else:
                    _Coords += [(_C[_I][0], _C[_I][1])]
                
                _I += 1
            
            #

            if is_region:
                _Coords += [(_Coords[0][0], _Coords[0][1])]
                g = shapely.geometry.Polygon(_Coords)
            else:
                g = LineString(_Coords)
            
            geoms += [g]
            
            if not g.is_valid:
                num_invalid_geoms += 1
                        
            if not g.is_simple:
                num_complex_geoms += 1
  
        i += 1

    return geoms, num_invalid_geoms, num_complex_geoms

def print_error(message):
    print(0)
    print(0)
    print(0)
    print(message)

# --------------------------------------------------------------------------------------------------------------------------

"""
h = loads("LineString (-1.01196172248803817 -0.07416267942583721, -0.96650717703349276 0.03349282296650724, -0.88516746411483238 0.11004784688995217, -0.90430622009569372 0.15311004784688997, -0.97129186602870798 0.11004784688995217, -0.96650717703349276 0.07655502392344504, -1.00717703349282273 0.05502392344497609, -1.03588516746411474 0.09569377990430628, -1.08133971291866016 0.09090909090909094, -1.11961722488038262 0.16028708133971292, -1.02631578947368407 0.17942583732057416, -0.92344497607655485 0.22009569377990434, -0.83492822966507174 0.21531100478468901, -0.78468899521531088 0.16507177033492826, -0.6866028708133971 0.18181818181818188, -0.76315789473684204 0.22248803827751196, -0.78229665071770316 0.30622009569377995, -0.8875598086124401 0.32296650717703357, -0.97607655502392343 0.33014354066985652, -1.06220095693779903 0.27033492822966509, -1.16746411483253576 0.22488038277511968, -1.25119617224880386 0.23684210526315796, -1.24880382775119614 0.16985645933014359, -1.25358851674641136 0.10287081339712922, -1.22248803827751185 0.01435406698564601, -1.31100478468899517 0.03827751196172258, -1.3157894736842104 0.10287081339712922, -1.34688995215310992 0.16507177033492826, -1.36363636363636354 0.21052631578947367, -1.32535885167464107 0.30143540669856461, -1.20574162679425823 0.30861244019138756, -1.10765550239234445 0.31818181818181823, -1.09330143540669855 0.38277511961722488, -1.22009569377990434 0.40669856459330145, -1.16507177033492804 0.46650717703349287, -1.09090909090909083 0.45215311004784686, -1.11004784688995217 0.51913875598086134, -1.04066985645932997 0.50717703349282295, -1.02392344497607635 0.44736842105263164, -0.94976076555023914 0.43540669856459335, -0.94497607655502391 0.53349282296650724, -0.98564593301435388 0.59569377990430628, -1.07177033492822948 0.61722488038277512, -1.16267942583732031 0.59330143540669855, -1.24880382775119614 0.59808612440191389, -1.31339712918660267 0.63397129186602874, -1.39952153110047806 0.61722488038277534, -1.44736842105263119 0.66746411483253609, -1.45933014354066959 0.74401913875598102, -1.46889952153110004 0.80861244019138778, -1.46411483253588481 0.87320574162679443, -1.48803827751196138 0.93062200956937824, -1.49043062200956888 1.0119617224880384, -1.46889952153110004 1.06459330143540676, -1.42583732057416235 1.09090909090909105, -1.37081339712918604 1.09330143540669877, -1.34210526315789425 1.05741626794258403, -1.36363636363636331 1.03349282296650724, -1.39712918660287055 0.98803827751196194, -1.40430622009569328 0.92344497607655518, -1.38755980861243966 0.84210526315789491, -1.38995215311004738 0.76076555023923464, -1.39234449760765511 0.69377990430622027, -1.33253588516746402 0.71531100478468901, -1.24401913875598069 0.78229665071770338, -1.24641148325358797 0.85645933014354092, -1.27990430622009521 0.90669856459330167, -1.27990430622009521 0.95933014354067003, -1.26315789473684159 1.01435406698564612, -1.20334928229665028 1.04066985645933041, -1.13636363636363602 1.01674641148325384, -1.11244019138755945 0.97607655502392365, -1.11483253588516718 0.93779904306220119, -1.15789473684210487 0.9114832535885169, -1.1339712918660283 0.87081339712918682, -1.11961722488038262 0.8133971291866029, -0.97607655502392343 0.83014354066985652, -0.89952153110047828 0.77511961722488043, -0.93540669856459302 0.70334928229665072, -0.86363636363636354 0.65789473684210531, -0.85406698564593286 0.49282296650717705, -0.80622009569377973 0.42583732057416268, -0.83014354066985629 0.53110047846889952, -0.74880382775119614 0.52631578947368429, -0.70574162679425823 0.58851674641148333, -0.66985645933014348 0.66267942583732053, -0.6004784688995215 0.69617224880382778, -0.64832535885167464 0.74641148325358853, -0.71770334928229662 0.75598086124401909, -0.72727272727272729 0.82775119617224879, -0.64832535885167464 0.87320574162679421, -0.55023923444976064 0.87799043062200954, -0.42105263157894735 0.91148325358851678, -0.33492822966507174 0.91148325358851678, -0.27272727272727271 0.90430622009569372, -0.19138755980861233 0.8995215311004785, -0.21531100478468868 0.96650717703349298, -0.25358851674641114 1.01913875598086134, -0.3205741626794254 1.01435406698564612, -0.36602870813397104 1.01913875598086134, -0.406698564593301 1.05741626794258403, -0.48564593301435366 1.07894736842105288, -0.54306220095693747 1.05263157894736858, -0.63636363636363602 1.06459330143540676, -0.67703349282296621 1.098086124401914, -0.68660287081339666 1.16985645933014371, -0.68899521531100438 1.22966507177033524, -0.67703349282296621 1.27033492822966521, -0.6196172248803824 1.2822966507177036, -0.56937799043062154 1.23923444976076569, -0.5334928229665068 1.17703349282296665, -0.45215311004784642 1.16028708133971303, -0.37559808612440149 1.11722488038277534, -0.32296650717703312 1.13636363636363669, -0.26555023923444931 1.13157894736842124, -0.22248803827751162 1.09569377990430628, -0.18421052631578916 1.04066985645933041, -0.16507177033492781 0.97607655502392365, -0.1052631578947365 0.96650717703349298, -0.12200956937799035 0.9138755980861244, -0.11004784688995217 0.83971291866028708, -0.23684210526315774 0.82296650717703357, -0.29904306220095678 0.84210526315789469, -0.36363636363636354 0.83014354066985652, -0.48325358851674638 0.80622009569377995, -0.46172248803827753 0.74401913875598091, -0.5023923444976075 0.67703349282296654, -0.50717703349282295 0.60765550239234445, -0.57655502392344493 0.5861244019138756, -0.61244019138755967 0.53110047846889952, -0.53349282296650702 0.48086124401913866, -0.45933014354066981 0.43301435406698563, -0.38995215311004783 0.4856459330143541, -0.37559808612440193 0.57655502392344493, -0.30382775119617222 0.60287081339712922, -0.24641148325358841 0.55263157894736836, -0.17942583732057416 0.47368421052631582, -0.13397129186602874 0.38038277511961727, -0.12918660287081329 0.299043062200957, -0.11004784688995217 0.19377990430622016, -0.21052631578947367 0.16746411483253587, -0.23444976076555024 0.22248803827751196, -0.3014354066985645 0.24162679425837319, -0.36124401913875581 0.24641148325358853, -0.34449760765550241 0.3133971291866029, -0.27990430622009566 0.30622009569377995, -0.22727272727272729 0.29665071770334928, -0.24880382775119614 0.35406698564593309, -0.24641148325358841 0.4138755980861244, -0.29186602870813383 0.48803827751196172, -0.33971291866028697 0.40909090909090906, -0.42583732057416257 0.32057416267942584, -0.5023923444976075 0.28708133971291872, -0.48325358851674638 0.16507177033492826, -0.42583732057416257 0.21052631578947367, -0.36124401913875581 0.18421052631578949, -0.2870813397129186 0.09330143540669855, -0.1866028708133971 0.09330143540669855, -0.17703349282296643 0.03110047846889952, -0.22248803827751185 -0.02870813397129179, -0.06220095693779903 0.04545454545454553, 0.02392344497607679 0.10765550239234456, -0.00478468899521522 0.18181818181818188, 0.08133971291866038 0.2129186602870814, 0.03110047846889952 0.27751196172248804, 0.0837320574162681 0.35167464114832536, -0.02392344497607657 0.37081339712918659, -0.07177033492822971 0.4282296650717704, -0.09569377990430605 0.48086124401913877, -0.13397129186602874 0.55023923444976086, -0.16507177033492826 0.62200956937799046, -0.22009569377990434 0.68181818181818188, -0.21052631578947367 0.74162679425837319, -0.13157894736842102 0.76794258373205748, -0.07894736842105265 0.73923444976076558, -0.09090909090909083 0.69138755980861244, -0.07416267942583721 0.62918660287081341, -0.05023923444976064 0.56937799043062198, -0.02392344497607657 0.50717703349282295, 0.0143540669856459 0.45454545454545459, 0.07177033492822993 0.43779904306220097, 0.09330143540669855 0.48803827751196172, 0.07177033492822993 0.55023923444976086, 0.05741626794258359 0.64832535885167464, 0.10765550239234445 0.72248803827751196, 0.18899521531100483 0.75119617224880386, 0.21770334928229684 0.81100478468899539, 0.1937799043062205 0.86842105263157909, 0.14832535885167508 0.83732057416267958, 0.10765550239234467 0.80861244019138778, 0.05741626794258425 0.83253588516746435, 0.09090909090909105 0.87081339712918682, 0.07416267942583787 0.92344497607655518, 0.04306220095693836 0.95693779904306242, 0.04066985645933063 1.04066985645933041, -0.02392344497607612 1.06459330143540676, 0.01913875598086157 1.098086124401914, 0.02870813397129202 1.16507177033492848, 0.08851674641148377 1.14354066985645941, 0.16985645933014371 1.12918660287081352, 0.18660287081339733 1.09090909090909105, 0.10765550239234467 1.08133971291866038, 0.14354066985645964 1.04306220095693791, 0.09090909090909105 1.01435406698564612, 0.12200956937799057 0.94497607655502414, 0.16267942583732098 0.9114832535885169, 0.18899521531100505 0.94976076555023947, 0.23684210526315819 0.96411483253588537, 0.27511961722488087 0.94497607655502414, 0.27511961722488087 0.88755980861244044, 0.29904306220095722 0.83014354066985663, 0.30622009569378017 0.77272727272727271, 0.33971291866028697 0.69377990430622005, 0.3133971291866029 0.65550239234449759, 0.24401913875598069 0.70574162679425845, 0.17703349282296665 0.67703349282296654, 0.13157894736842124 0.60526315789473684, 0.15311004784688986 0.52392344497607657, 0.16507177033492848 0.45933014354066992, 0.15789473684210531 0.38038277511961727, 0.23205741626794252 0.49043062200956944, 0.18899521531100483 0.55023923444976086, 0.20334928229665072 0.62440191387559807, 0.26794258373205748 0.62200956937799046, 0.32296650717703335 0.4856459330143541, 0.299043062200957 0.36363636363636365, 0.397129186602871 0.40909090909090906, 0.43779904306220097 0.45933014354066992, 0.42822966507177052 0.50717703349282295, 0.44497607655502414 0.53588516746411485, 0.45454545454545459 0.61722488038277512, 0.54066985645932997 0.64593301435406703, 0.598086124401914 0.64354066985645941, 0.65550239234449759 0.59808612440191389, 0.71770334928229662 0.53827751196172247, 0.72248803827751207 0.50239234449760772, 0.6866028708133971 0.45933014354066992, 0.65789473684210531 0.51196172248803828, 0.602870813397129 0.56459330143540676, 0.53827751196172269 0.58851674641148333, 0.5143540669856459 0.5430622009569378, 0.56698564593301448 0.52153110047846885, 0.59330143540669855 0.48086124401913877, 0.55263157894736858 0.43062200956937802, 0.47846889952153093 0.41148325358851678, 0.44258373205741641 0.35167464114832536, 0.39952153110047872 0.28229665071770338, 0.42105263157894757 0.23684210526315819, 0.42822966507177052 0.17464114832535893, 0.38277511961722532 0.12200956937799057, 0.31578947368421062 0.0861244019138756, 0.22488038277511979 0.05023923444976075, 0.11961722488038262 0.04545454545454553, 0.02631578947368407 -0.01674641148325362, -0.06459330143540654 -0.04784688995215314, -0.1004784688995215 -0.12440191387559807)")
print(len(h.coords))
sys.exit()
"""

num_invalid_geoms = 0
num_complex_geoms = 0
geoms = []
    
precision = '.4f'

p, q, op, n_obs, debug, type_transf = get_input()

# Transformation.

s_exec_time = 0
e_exec_time = 0

if TESTING:
    tests_r = []
    N = NTESTS
    
    # 4, 16, 256
    # test 00
    seg_wkt = 'LineString (-1.05023923444976064 -0.24880382775119614, -0.01196172248803817 -0.36602870813397126)'
    line_wkt = 'LineString (-1.01196172248803817 -0.07416267942583721, -0.6866028708133971 0.18181818181818188, -0.36124401913875581 0.18421052631578949, -0.1004784688995215 -0.12440191387559807)'
				
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 01
    seg_wkt = 'LineString (-1.05023923444976064 -0.24880382775119614, -0.01196172248803817 -0.36602870813397126)'
    line_wkt = 'LineString (-1.01196172248803817 -0.07416267942583721, -0.88516746411483238 0.11004784688995217, -0.6866028708133971 0.18181818181818188, -0.76315789473684204 0.22248803827751196, -0.78229665071770316 0.30622009569377995, -0.67464114832535871 0.32296650717703357, -0.5885167464114831 0.38277511961722488, -0.5023923444976075 0.28708133971291872, -0.48325358851674638 0.16507177033492826, -0.42583732057416257 0.21052631578947367, -0.36124401913875581 0.18421052631578949, -0.2870813397129186 0.09330143540669855, -0.1866028708133971 0.09330143540669855, -0.17703349282296643 0.03110047846889952, -0.27511961722488021 -0.01913875598086112, -0.1004784688995215 -0.12440191387559807)'
				
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 02
    seg_wkt = 'LineString (-1.05023923444976064 -0.24880382775119614, -0.01196172248803817 -0.36602870813397126)'
    line_wkt = 'LineString (-1.01196172248803817 -0.07416267942583721, -0.96650717703349276 0.03349282296650724, -0.88516746411483238 0.11004784688995217, -0.90430622009569372 0.15311004784688997, -0.97129186602870798 0.11004784688995217, -0.96650717703349276 0.07655502392344504, -1.00717703349282273 0.05502392344497609, -1.03588516746411474 0.09569377990430628, -1.08133971291866016 0.09090909090909094, -1.11961722488038262 0.16028708133971292, -1.02631578947368407 0.17942583732057416, -0.92344497607655485 0.22009569377990434, -0.83492822966507174 0.21531100478468901, -0.78468899521531088 0.16507177033492826, -0.6866028708133971 0.18181818181818188, -0.76315789473684204 0.22248803827751196, -0.78229665071770316 0.30622009569377995, -0.8875598086124401 0.32296650717703357, -0.97607655502392343 0.33014354066985652, -1.06220095693779903 0.27033492822966509, -1.16746411483253576 0.22488038277511968, -1.25119617224880386 0.23684210526315796, -1.24880382775119614 0.16985645933014359, -1.25358851674641136 0.10287081339712922, -1.22248803827751185 0.01435406698564601, -1.31100478468899517 0.03827751196172258, -1.3157894736842104 0.10287081339712922, -1.34688995215310992 0.16507177033492826, -1.36363636363636354 0.21052631578947367, -1.32535885167464107 0.30143540669856461, -1.20574162679425823 0.30861244019138756, -1.10765550239234445 0.31818181818181823, -1.09330143540669855 0.38277511961722488, -1.22009569377990434 0.40669856459330145, -1.16507177033492804 0.46650717703349287, -1.09090909090909083 0.45215311004784686, -1.11004784688995217 0.51913875598086134, -1.04066985645932997 0.50717703349282295, -1.02392344497607635 0.44736842105263164, -0.94976076555023914 0.43540669856459335, -0.94497607655502391 0.53349282296650724, -0.98564593301435388 0.59569377990430628, -1.07177033492822948 0.61722488038277512, -1.16267942583732031 0.59330143540669855, -1.24880382775119614 0.59808612440191389, -1.31339712918660267 0.63397129186602874, -1.39952153110047806 0.61722488038277534, -1.44736842105263119 0.66746411483253609, -1.45933014354066959 0.74401913875598102, -1.46889952153110004 0.80861244019138778, -1.46411483253588481 0.87320574162679443, -1.48803827751196138 0.93062200956937824, -1.49043062200956888 1.0119617224880384, -1.46889952153110004 1.06459330143540676, -1.42583732057416235 1.09090909090909105, -1.37081339712918604 1.09330143540669877, -1.34210526315789425 1.05741626794258403, -1.36363636363636331 1.03349282296650724, -1.39712918660287055 0.98803827751196194, -1.40430622009569328 0.92344497607655518, -1.38755980861243966 0.84210526315789491, -1.38995215311004738 0.76076555023923464, -1.39234449760765511 0.69377990430622027, -1.33253588516746402 0.71531100478468901, -1.24401913875598069 0.78229665071770338, -1.24641148325358797 0.85645933014354092, -1.27990430622009521 0.90669856459330167, -1.27990430622009521 0.95933014354067003, -1.26315789473684159 1.01435406698564612, -1.20334928229665028 1.04066985645933041, -1.13636363636363602 1.01674641148325384, -1.11244019138755945 0.97607655502392365, -1.11483253588516718 0.93779904306220119, -1.15789473684210487 0.9114832535885169, -1.1339712918660283 0.87081339712918682, -1.11961722488038262 0.8133971291866029, -0.97607655502392343 0.83014354066985652, -0.89952153110047828 0.77511961722488043, -0.93540669856459302 0.70334928229665072, -0.86363636363636354 0.65789473684210531, -0.85406698564593286 0.49282296650717705, -0.80622009569377973 0.42583732057416268, -0.83014354066985629 0.53110047846889952, -0.74880382775119614 0.52631578947368429, -0.70574162679425823 0.58851674641148333, -0.66985645933014348 0.66267942583732053, -0.6004784688995215 0.69617224880382778, -0.64832535885167464 0.74641148325358853, -0.71770334928229662 0.75598086124401909, -0.72727272727272729 0.82775119617224879, -0.64832535885167464 0.87320574162679421, -0.55023923444976064 0.87799043062200954, -0.42105263157894735 0.91148325358851678, -0.33492822966507174 0.91148325358851678, -0.27272727272727271 0.90430622009569372, -0.19138755980861233 0.8995215311004785, -0.21531100478468868 0.96650717703349298, -0.25358851674641114 1.01913875598086134, -0.3205741626794254 1.01435406698564612, -0.36602870813397104 1.01913875598086134, -0.406698564593301 1.05741626794258403, -0.48564593301435366 1.07894736842105288, -0.54306220095693747 1.05263157894736858, -0.63636363636363602 1.06459330143540676, -0.67703349282296621 1.098086124401914, -0.68660287081339666 1.16985645933014371, -0.68899521531100438 1.22966507177033524, -0.67703349282296621 1.27033492822966521, -0.6196172248803824 1.2822966507177036, -0.56937799043062154 1.23923444976076569, -0.5334928229665068 1.17703349282296665, -0.45215311004784642 1.16028708133971303, -0.37559808612440149 1.11722488038277534, -0.32296650717703312 1.13636363636363669, -0.26555023923444931 1.13157894736842124, -0.22248803827751162 1.09569377990430628, -0.18421052631578916 1.04066985645933041, -0.16507177033492781 0.97607655502392365, -0.1052631578947365 0.96650717703349298, -0.12200956937799035 0.9138755980861244, -0.11004784688995217 0.83971291866028708, -0.23684210526315774 0.82296650717703357, -0.29904306220095678 0.84210526315789469, -0.36363636363636354 0.83014354066985652, -0.48325358851674638 0.80622009569377995, -0.46172248803827753 0.74401913875598091, -0.5023923444976075 0.67703349282296654, -0.50717703349282295 0.60765550239234445, -0.57655502392344493 0.5861244019138756, -0.61244019138755967 0.53110047846889952, -0.53349282296650702 0.48086124401913866, -0.45933014354066981 0.43301435406698563, -0.38995215311004783 0.4856459330143541, -0.37559808612440193 0.57655502392344493, -0.30382775119617222 0.60287081339712922, -0.24641148325358841 0.55263157894736836, -0.17942583732057416 0.47368421052631582, -0.13397129186602874 0.38038277511961727, -0.12918660287081329 0.299043062200957, -0.11004784688995217 0.19377990430622016, -0.21052631578947367 0.16746411483253587, -0.23444976076555024 0.22248803827751196, -0.3014354066985645 0.24162679425837319, -0.36124401913875581 0.24641148325358853, -0.34449760765550241 0.3133971291866029, -0.27990430622009566 0.30622009569377995, -0.22727272727272729 0.29665071770334928, -0.24880382775119614 0.35406698564593309, -0.24641148325358841 0.4138755980861244, -0.29186602870813383 0.48803827751196172, -0.33971291866028697 0.40909090909090906, -0.42583732057416257 0.32057416267942584, -0.5023923444976075 0.28708133971291872, -0.48325358851674638 0.16507177033492826, -0.42583732057416257 0.21052631578947367, -0.36124401913875581 0.18421052631578949, -0.2870813397129186 0.09330143540669855, -0.1866028708133971 0.09330143540669855, -0.17703349282296643 0.03110047846889952, -0.22248803827751185 -0.02870813397129179, -0.06220095693779903 0.04545454545454553, 0.02392344497607679 0.10765550239234456, -0.00478468899521522 0.18181818181818188, 0.08133971291866038 0.2129186602870814, 0.03110047846889952 0.27751196172248804, 0.0837320574162681 0.35167464114832536, -0.02392344497607657 0.37081339712918659, -0.07177033492822971 0.4282296650717704, -0.09569377990430605 0.48086124401913877, -0.13397129186602874 0.55023923444976086, -0.16507177033492826 0.62200956937799046, -0.22009569377990434 0.68181818181818188, -0.21052631578947367 0.74162679425837319, -0.13157894736842102 0.76794258373205748, -0.07894736842105265 0.73923444976076558, -0.09090909090909083 0.69138755980861244, -0.07416267942583721 0.62918660287081341, -0.05023923444976064 0.56937799043062198, -0.02392344497607657 0.50717703349282295, 0.0143540669856459 0.45454545454545459, 0.07177033492822993 0.43779904306220097, 0.09330143540669855 0.48803827751196172, 0.07177033492822993 0.55023923444976086, 0.05741626794258359 0.64832535885167464, 0.10765550239234445 0.72248803827751196, 0.18899521531100483 0.75119617224880386, 0.21770334928229684 0.81100478468899539, 0.1937799043062205 0.86842105263157909, 0.14832535885167508 0.83732057416267958, 0.10765550239234467 0.80861244019138778, 0.05741626794258425 0.83253588516746435, 0.09090909090909105 0.87081339712918682, 0.07416267942583787 0.92344497607655518, 0.04306220095693836 0.95693779904306242, 0.04066985645933063 1.04066985645933041, -0.02392344497607612 1.06459330143540676, 0.01913875598086157 1.098086124401914, 0.02870813397129202 1.16507177033492848, 0.08851674641148377 1.14354066985645941, 0.16985645933014371 1.12918660287081352, 0.18660287081339733 1.09090909090909105, 0.10765550239234467 1.08133971291866038, 0.14354066985645964 1.04306220095693791, 0.09090909090909105 1.01435406698564612, 0.12200956937799057 0.94497607655502414, 0.16267942583732098 0.9114832535885169, 0.18899521531100505 0.94976076555023947, 0.23684210526315819 0.96411483253588537, 0.27511961722488087 0.94497607655502414, 0.27511961722488087 0.88755980861244044, 0.29904306220095722 0.83014354066985663, 0.30622009569378017 0.77272727272727271, 0.33971291866028697 0.69377990430622005, 0.3133971291866029 0.65550239234449759, 0.24401913875598069 0.70574162679425845, 0.17703349282296665 0.67703349282296654, 0.13157894736842124 0.60526315789473684, 0.15311004784688986 0.52392344497607657, 0.16507177033492848 0.45933014354066992, 0.15789473684210531 0.38038277511961727, 0.23205741626794252 0.49043062200956944, 0.18899521531100483 0.55023923444976086, 0.20334928229665072 0.62440191387559807, 0.26794258373205748 0.62200956937799046, 0.32296650717703335 0.4856459330143541, 0.299043062200957 0.36363636363636365, 0.397129186602871 0.40909090909090906, 0.43779904306220097 0.45933014354066992, 0.42822966507177052 0.50717703349282295, 0.44497607655502414 0.53588516746411485, 0.45454545454545459 0.61722488038277512, 0.54066985645932997 0.64593301435406703, 0.598086124401914 0.64354066985645941, 0.65550239234449759 0.59808612440191389, 0.71770334928229662 0.53827751196172247, 0.72248803827751207 0.50239234449760772, 0.6866028708133971 0.45933014354066992, 0.65789473684210531 0.51196172248803828, 0.602870813397129 0.56459330143540676, 0.53827751196172269 0.58851674641148333, 0.5143540669856459 0.5430622009569378, 0.56698564593301448 0.52153110047846885, 0.59330143540669855 0.48086124401913877, 0.55263157894736858 0.43062200956937802, 0.47846889952153093 0.41148325358851678, 0.44258373205741641 0.35167464114832536, 0.39952153110047872 0.28229665071770338, 0.42105263157894757 0.23684210526315819, 0.42822966507177052 0.17464114832535893, 0.38277511961722532 0.12200956937799057, 0.31578947368421062 0.0861244019138756, 0.22488038277511979 0.05023923444976075, 0.11961722488038262 0.04545454545454553, 0.02631578947368407 -0.01674641148325362, -0.06459330143540654 -0.04784688995215314, -0.1004784688995215 -0.12440191387559807)'
				
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]

    # >>>>>
    
    # test 0
    seg_wkt = 'LINESTRING(-140.50924455177292316 -54.26986651320897437, -140.06526900166338123 -54.25620572705175704)'
    line_wkt = 'LINESTRING(-140.38800507462761402 -54.16228782222089677, -140.33848472480769942 -54.09569148970446406, -140.39312786943656874 -54.0393407468059479, -140.34360751961665414 -54.02738755891838451, -140.30774795595397109 -54.06495472085072862, -140.20016926496589349 -54.05812432777211995, -140.16089450476388834 -54.12130546374925189, -140.10795895840468006 -54.08715349835620856, -140.07039179647230753 -54.10252188278307273, -140.13698812898874735 -54.15033463433333338, -140.07039179647233595 -54.18960939453533143)'
				
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 1
    seg_wkt = 'LineString (-1.06326530612244907 -0.46530612244897962, 0.92448979591836755 -0.42448979591836755)'
    line_wkt = 'LineString (-1.08775510204081627 -0.07755102040816331, -0.76938775510204072 -0.13469387755102047, -0.71224489795918355 0.00408163265306116, -0.61020408163265305 0.13061224489795908, -0.42653061224489797 0.22448979591836726, -0.23061224489795906 0.19999999999999996, -0.0428571428571427 0.08571428571428563, 0.0346938775510206 -0.03673469387755102, -0.01836734693877551 -0.08571428571428585, -0.07551020408163245 -0.00816326530612255, -0.15714285714285703 0.04081632653061218, -0.21428571428571419 0, -0.23877551020408161 0.02040816326530603, -0.17755102040816317 0.08979591836734691, -0.34897959183673466 0.13061224489795908, -0.45510204081632644 0.08979591836734691, -0.57755102040816331 -0.02040816326530615, -0.45510204081632644 -0.03673469387755102, -0.393877551020408 -0.08979591836734713, -0.40612244897959182 -0.16734693877551021, -0.44693877551020411 -0.13877551020408174, -0.44693877551020411 -0.106122448979592, -0.53673469387755102 -0.07346938775510203, -0.63877551020408152 -0.07346938775510203, -0.63877551020408152 -0.14693877551020429, -0.2795918367346939 -0.23673469387755119, 0.14489795918367365 -0.1959183673469389, 0.2551020408163267 0.00816326530612232, 0.1530612244897962 0.14693877551020396, -0.23469387755102034 0.39183673469387748, -0.47551020408163258 0.38775510204081631, -0.72448979591836737 0.27346938775510199, -0.80204081632653068 0.05306122448979589, -0.86326530612244889 0.05306122448979589, -0.83469387755102042 0.19999999999999996, -0.71632653061224483 0.39183673469387748, -0.49591836734693873 0.48571428571428565, -0.23469387755102034 0.49387755102040809, -0.15714285714285703 0.61632653061224485, -0.02244897959183678 0.63265306122448983, 0.11632653061224518 0.55918367346938769, 0.13265306122448983 0.4408163265306122, 0.06326530612244907 0.4408163265306122, 0.0591836734693878 0.51020408163265296, -0.03469387755102016 0.55102040816326525, -0.07142857142857117 0.4408163265306122, 0.07959183673469417 0.32244897959183672, 0.30000000000000027 0.18775510204081625, 0.47142857142857153 0.06530612244897949, 0.4918367346938779 -0.03265306122448997, 0.40612244897959204 -0.05714285714285716, 0.40204081632653077 0.04897959183673462, 0.35306122448979593 0.04897959183673462, 0.34489795918367339 -0.10204081632653073, 0.47142857142857153 -0.22448979591836737, 0.81836734693877577 -0.08163265306122458)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]

    # test 2
    seg_wkt = 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)'
    line_wkt = 'LineString (-1.1775510204081634 -0.00408163265306127, -0.81428571428571428 0.19183673469387752, -0.61020408163265305 -0.00816326530612255, -0.42653061224489797 0.15102040816326523, -0.6020408163265305 0.32244897959183672, -0.5122448979591836 0.42040816326530606, -0.19387755102040805 0.16326530612244894, -0.083673469387755 -0.05714285714285716)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]

    # test 3
    seg_wkt = 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)'
    line_wkt = 'LineString (-0.73497709287796686 1.13821740941274596, 0.28299875052061729 1.63827571845064712, 0.5062390670553949 1.29002082465639445, 0.24728029987505318 1.25430237401082989, 0.24728029987505318 1.42396501457726066, -0.0027488546438974 1.29002082465639445, 0.29192836318200843 1.04892128279883501, 0.72054977092878092 1.07571012078300843, 0.63125364431486997 1.45075385256143408, 0.20263223656809748 1.83472719700125109, 0.90807163681799397 1.71864223240316694, 0.88128279883382055 1.34359850062474084, 1.25632653061224664 1.1024989587671814, 2.05999167013744477 1.3793169512703054)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]

    # test 4
    seg_wkt = 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)'
    line_wkt = 'LineString (-0.73497709287796686 1.13821740941274596, -0.2815076711234259 1.38834907224030912, -0.23017151017008153 1.57943367134442392, -0.17883534921673727 1.45109326896106317, 0.23451459850912548 1.63257170056694223, 0.5062390670553949 1.29002082465639445, 0.24728029987505318 1.25430237401082989, 0.24728029987505318 1.42396501457726066, -0.0027488546438974 1.29002082465639445, 0.29192836318200843 1.04892128279883501, 0.72054977092878092 1.07571012078300843, 0.63125364431486997 1.45075385256143408, 0.20263223656809748 1.83472719700125109, 0.90807163681799397 1.71864223240316694, 0.88128279883382055 1.34359850062474084, 1.25632653061224664 1.1024989587671814, 2.05999167013744477 1.3793169512703054)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]

    # test 5
    seg_wkt = 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)'
    line_wkt = 'LineString (-0.73497709287796686 1.13821740941274596, -0.27865566218157345 1.43112920636809582, -0.34710387678603249 1.49672541203070253, -0.28435968006527834 1.57658166240257147, -0.19309539392599961 1.69351402901852222, -0.2045034296934094 1.4596492957866205, 0.0350653214221972 1.47390934049588296, 0.23451459850912548 1.63257170056694223, 0.5062390670553949 1.29002082465639445, 0.24728029987505318 1.25430237401082989, 0.24728029987505318 1.42396501457726066, -0.0027488546438974 1.29002082465639445, 0.29192836318200843 1.04892128279883501, 0.72054977092878092 1.07571012078300843, 0.63125364431486997 1.45075385256143408, 0.20263223656809748 1.83472719700125109, 0.90807163681799397 1.71864223240316694, 0.88128279883382055 1.34359850062474084, 1.25632653061224664 1.1024989587671814, 2.05999167013744477 1.3793169512703054)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 6
    seg_wkt = 'LineString (-0.75283631820074914 0.70959600166597325, 2.35466888796335061 0.64708871303623572)'
    line_wkt = 'LineString (-0.73497709287796686 1.13821740941274596, -0.27865566218157345 1.43112920636809582, -0.34710387678603249 1.49672541203070253, -0.37384146061589896 1.44467624884189605, -0.42802963051109572 1.46178830249301073, -0.43088163945294816 1.51312446344635498, -0.40236155003442353 1.55020057969043701, -0.36243342484848912 1.59298071381822393, -0.29113320130217757 1.5851376892281297, -0.24906606940985382 1.65073389489073641, -0.19309539392599961 1.69351402901852222, -0.18133085704085783 1.58870270040544526, -0.23694503140698076 1.53166252156839611, -0.2045034296934094 1.4596492957866205, 0.0350653214221972 1.47390934049588296, 0.10387003714438825 1.5559045975741419, 0.06251590748752767 1.64574287924249441, 0.15948421151051129 1.60724075852748616, 0.23451459850912548 1.63257170056694223, 0.5062390670553949 1.29002082465639445, 0.24728029987505318 1.25430237401082989, 0.24728029987505318 1.42396501457726066, -0.0027488546438974 1.29002082465639445, 0.2022643456382982 1.18656943960424877, 0.29352863177757693 1.17516140383683876, 0.33060474802165896 1.20368149325536344, 0.37053287320759343 1.19227345748795366, 0.39334894474241322 1.15234533230201919, 0.31919671225424917 1.10386118029052738, 0.29192836318200843 1.04892128279883501, 0.72054977092878092 1.07571012078300843, 0.63125364431486997 1.45075385256143408, 0.20263223656809748 1.83472719700125109, 0.90807163681799397 1.71864223240316694, 0.88128279883382055 1.34359850062474084, 1.25632653061224664 1.1024989587671814, 2.05999167013744477 1.3793169512703054)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]

    # test 7
    seg_wkt = 'LineString (-0.96530612244897962 -0.416326530612245, 0.82653061224489832 -0.33061224489795937)'
    line_wkt = 'LineString (-0.9408163265306122 -0.18367346938775531, -0.63469387755102047 0.25306122448979584, -0.25918367346938775 -0.21632653061224505, -0.00204081632653041 0.12653061224489792, 0.21836734693877569 -0.15918367346938789, 0.38979591836734695 0.00408163265306116, 0.54489795918367356 -0.08163265306122458, 0.63877551020408196 -0.01632653061224509, 0.68367346938775508 -0.05714285714285716, 0.75714285714285712 0.01632653061224476)'
    
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 8
    seg_wkt = 'LineString (-1.14489795918367343 -0.60408163265306136, 1.41836734693877586 -0.5346938775510206)'
    line_wkt = 'LineString (-1.25102040816326543 0.05306122448979589, -0.67959183673469381 -0.17142857142857149, -0.82653061224489788 0.08163265306122436, -0.62653061224489792 0.33877551020408159, -0.3816326530612244 0.44897959183673464, -0.13673469387755088 0.45306122448979591, -0.10816326530612241 0.59591836734693882, -0.21020408163265292 0.62040816326530601, -0.19387755102040805 0.74285714285714288, -0.05102040816326525 0.67755102040816317, 0.03877551020408188 0.54693877551020398, 0.01836734693877551 0.41632653061224489, 0.14489795918367365 0.32653061224489788, 0.35306122448979593 0.22448979591836726, 0.51632653061224509 0.26530612244897955, 0.54081632653061229 0.43673469387755093, 0.32040816326530619 0.43673469387755093, 0.33265306122449001 0.53469387755102038, 0.63061224489795942 0.52244897959183667, 0.59795918367346967 0.2857142857142857, 0.75306122448979629 0.34285714285714275, 0.73265306122448992 0.47346938775510194, 0.83061224489795915 0.36734693877551017, 0.85918367346938807 0.19591836734693868, 0.74897959183673501 0.24489795918367341, 0.48775510204081662 0.11020408163265294, 0.66734693877551043 -0.06122448979591844, 0.74081632653061247 -0.23673469387755119, 1.32448979591836746 0.05306122448979589)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 9
    seg_wkt = 'LineString (-1.14489795918367343 -0.60408163265306136, 1.41836734693877586 -0.5346938775510206)'
    line_wkt = 'LineString (-1.21836734693877546 -0.26122448979591839, -0.80612244897959173 -0.34693877551020424, -0.74081632653061225 -0.21224489795918378, -0.62653061224489792 -0.20408163265306123, -0.52040816326530615 -0.08571428571428585, -0.3040816326530611 -0.106122448979592, -0.14489795918367343 -0.08979591836734713, 0.01428571428571423 -0.13469387755102047, 0.25918367346938798 -0.08979591836734713, 0.41020408163265332 -0.14285714285714302, 0.43469387755102051 -0.01632653061224509, 0.21428571428571441 0.04897959183673462, -0.21020408163265292 0.03673469387755091, -0.3693877551020408 0.1020408163265305, -0.01428571428571423 0.10612244897959178, -0.08775510204081627 0.28979591836734686, -0.01020408163265296 0.28979591836734686, 0.00612244897959213 0.17551020408163254, 0.11224489795918391 0.17551020408163254, 0.09183673469387754 0.59999999999999998, 0.14897959183673493 0.59183673469387754, 0.16938775510204085 0.12653061224489792, 0.2551020408163267 0.12244897959183665, 0.24693877551020416 0.80408163265306121, -0.71224489795918355 0.80000000000000004, -0.54081632653061229 0.17959183673469381, -0.87551020408163271 -0.05306122448979611, -0.87142857142857144 0.03265306122448974, -0.62653061224489792 0.18367346938775508, -0.78571428571428559 0.8571428571428571, 0.32857142857142874 0.8571428571428571, 0.32448979591836746 0.14693877551020396, 0.54489795918367356 0.00408163265306116, 0.53265306122448974 -0.106122448979592, 0.57755102040816331 -0.20000000000000018, 0.67142857142857171 -0.30204081632653068, 1.41020408163265332 -0.21224489795918378)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 10
    seg_wkt = 'LineString (-1.25102040816326543 -0.73877551020408183, 0.93673469387755137 -0.59183673469387776)'
    line_wkt = 'LineString (-1.3040816326530611 -0.17551020408163276, -0.83469387755102042 -0.26938775510204094, -0.77346938775510199 -0.06938775510204098, -0.86734693877551017 0.02448979591836731, -0.94489795918367347 -0.01224489795918382, -0.97346938775510194 0.04081632653061218, -0.86734693877551017 0.1020408163265305, -0.78163265306122454 0.02857142857142847, -0.71632653061224483 0.00816326530612232, -0.67959183673469381 0.10612244897959178, -0.76530612244897966 0.19999999999999996, -1.00612244897959191 0.17551020408163254, -1.1489795918367347 0.06122448979591832, -1.15714285714285703 0.09387755102040807, -1.04285714285714293 0.19183673469387752, -0.89183673469387759 0.2326530612244897, -1.01428571428571423 0.33877551020408159, -0.9285714285714286 0.35510204081632646, -0.83061224489795915 0.23673469387755097, -0.66326530612244894 0.2204081632653061, -0.6020408163265305 0.14285714285714279, -0.51632653061224487 0.20408163265306112, -0.41428571428571415 0.30204081632653057, -0.25102040816326521 0.26938775510204072, -0.17755102040816317 0.33061224489795915, -0.0428571428571427 0.39591836734693875, 0.03061224489795933 0.25714285714285712, -0.02653061224489806 0.25306122448979584, -0.03877551020408143 0.31836734693877544, -0.15714285714285703 0.28979591836734686, -0.18571428571428572 0.19999999999999996, -0.09591836734693882 0.06122448979591832, 0.00612244897959213 -0.05714285714285716, 0.16530612244897958 0.03265306122448974, 0.21020408163265314 0.1673469387755101, 0.16938775510204085 0.30204081632653057, 0.10000000000000009 0.39999999999999991, -0.11224489795918369 0.47755102040816322, -0.29999999999999982 0.42040816326530606, -0.29591836734693877 0.47346938775510194, -0.13673469387755088 0.52244897959183667, 0.12040816326530601 0.45714285714285707, 0.27551020408163263 0.25714285714285712, 0.24285714285714288 0.01632653061224476, 0.02653061224489806 -0.1183673469387756, -0.07551020408163245 -0.21632653061224505, 0.14897959183673493 -0.33061224489795937, 0.83061224489795915 -0.05306122448979611)'
    
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 11
    seg_wkt = 'LineString (-1.25102040816326543 -0.73877551020408183, 0.93673469387755137 -0.59183673469387776)'
    line_wkt = 'LineString (-1.3040816326530611 -0.17551020408163276, -0.83469387755102042 -0.26938775510204094, -0.77346938775510199 -0.06938775510204098, -0.86734693877551017 0.02448979591836731, -0.94489795918367347 -0.01224489795918382, -0.97346938775510194 0.04081632653061218, -0.86734693877551017 0.1020408163265305, -0.78163265306122454 0.02857142857142847, -0.71632653061224483 0.00816326530612232, -0.67959183673469381 0.10612244897959178, -0.76530612244897966 0.19999999999999996, -0.66326530612244894 0.2204081632653061, -0.6020408163265305 0.14285714285714279, -0.51632653061224487 0.20408163265306112, -0.41428571428571415 0.30204081632653057, -0.25102040816326521 0.26938775510204072, -0.17755102040816317 0.33061224489795915, -0.15714285714285703 0.28979591836734686, -0.18571428571428572 0.19999999999999996, -0.09591836734693882 0.06122448979591832, 0.00612244897959213 -0.05714285714285716, 0.16530612244897958 0.03265306122448974, 0.24285714285714288 0.01632653061224476, 0.02653061224489806 -0.1183673469387756, -0.07551020408163245 -0.21632653061224505, 0.14897959183673493 -0.33061224489795937, 0.83061224489795915 -0.05306122448979611)'

    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test 12
    seg_wkt = 'LineString (-0.54489795918367334 -0.44081632653061242, -0.32040816326530597 -0.44081632653061242)'
    line_wkt = 'LineString (-0.51632653061224487 -0.02448979591836742, -0.71632653061224483 -0.18367346938775531, -0.92040816326530606 -0.15918367346938789, -0.81836734693877555 0.08163265306122436, -1.28775510204081645 0.11020408163265294, -1.30816326530612237 0.21632653061224483, -1.52040816326530615 0.30612244897959173, -1.04693877551020398 0.86122448979591837, -0.95714285714285707 0.75102040816326532, -1.14489795918367343 0.62040816326530601, -1.0714285714285714 0.48571428571428565, -1.13265306122448983 0.45306122448979591, -1.18571428571428572 0.54693877551020398, -1.35306122448979593 0.3755102040816326, -1.18571428571428572 0.19183673469387752, -1.0591836734693878 0.36326530612244889, -0.97755102040816322 0.33469387755102031, -1.10000000000000009 0.19999999999999996, -0.81836734693877555 0.1551020408163265, -0.69999999999999996 0.38775510204081631, -0.83469387755102042 0.46530612244897951, -0.76122448979591839 0.55918367346938769, -0.71224489795918355 0.52244897959183667, -0.74489795918367352 0.47755102040816322, -0.65510204081632639 0.43673469387755093, -0.49591836734693873 0.70204081632653059, -0.41428571428571415 0.59183673469387754, -0.59387755102040818 0.28163265306122442, -0.50408163265306127 0.25306122448979584, -0.2918367346938775 0.67755102040816317, -0.65918367346938767 0.82448979591836735, -0.80612244897959173 0.62040816326530601, -0.87959183673469377 0.62857142857142856, -0.7204081632653061 0.93061224489795924, -0.35714285714285698 0.78367346938775506, -0.31224489795918364 0.90612244897959182, -0.19795918367346932 0.87755102040816324, -0.25510204081632648 0.77551020408163263, -0.13673469387755088 0.72653061224489801, -0.35714285714285698 0.24081632653061213, -0.27142857142857135 0.1551020408163265, -0.12448979591836729 0.13061224489795908, 0.05102040816326525 0.32244897959183672, -0.10816326530612241 0.49387755102040809, 0.14897959183673493 0.72653061224489801, 0.58979591836734713 0.43673469387755093, 0.55306122448979611 0.35510204081632646, 0.16938775510204085 0.61632653061224485, 0.00612244897959213 0.50204081632653064, 0.10816326530612264 0.34285714285714275, 0.18571428571428594 0.33877551020408159, 0.13265306122448983 0.27346938775510199, 0.04285714285714315 0.23673469387755097, 0.01428571428571423 0.08571428571428563, 0.28775510204081645 -0.28979591836734708, 0.04693877551020398 -0.35918367346938784, -0.03469387755102016 -0.23265306122448992, 0.11224489795918391 -0.21224489795918378, 0.10816326530612264 -0.26530612244897966, 0.06734693877551035 -0.26530612244897966, 0.06734693877551035 -0.30612244897959195, 0.18571428571428594 -0.25714285714285734, -0.00612244897959169 0.02448979591836731, -0.29591836734693877 0.04081632653061218, -0.12448979591836729 -0.10204081632653073, -0.23469387755102034 -0.17551020408163276, -0.36122448979591826 -0.07346938775510203)'
    
    tests_r += [tests(seg_wkt, line_wkt, N, 1)]
    
    # test concavity with a face

    seg_wkt = 'LineString (-1.25102040816326543 -0.73877551020408183, 0.93673469387755137 -0.59183673469387776)'
    line_wkt = 'LineString (-1.3040816326530611 -0.17551020408163276, -0.83469387755102042 -0.26938775510204094, -0.77346938775510199 -0.06938775510204098, -0.86734693877551017 0.02448979591836731, -0.94489795918367347 -0.01224489795918382, -0.97346938775510194 0.04081632653061218, -0.86734693877551017 0.1020408163265305, -0.78163265306122454 0.02857142857142847, -0.72123661799388505 0.03752937718888849, -0.67959183673469381 0.10612244897959178, -0.74519746809134857 0.19218577327251635, -0.67549317689872757 0.24446399166698207, -0.6020408163265305 0.14285714285714279, -0.51632653061224487 0.20408163265306112, -0.41628034402616809 0.25971180536536792, -0.2790500207406954 0.25317702806605974, -0.18571428571428572 0.19999999999999996, -0.09591836734693882 0.06122448979591832, -0.06775888806306307 -0.0539575050014266, -0.05904585166398535 -0.12366179619404771, -0.07551020408163245 -0.21632653061224505, 0.14897959183673493 -0.33061224489795937, 0.83061224489795915 -0.05306122448979611)'
    face_wkt = 'Polygon ((-0.65918367346938767 -0.19183673469387763, -0.65510204081632639 -0.00408163265306127, -0.59795918367346923 0.07346938775510192, -0.56530612244897949 -0.02448979591836742, -0.5244897959183672 -0.00816326530612255, -0.55306122448979589 0.1020408163265305, -0.42244897959183669 0.18775510204081625, -0.26734693877551008 0.18775510204081625, -0.15306122448979576 0.05306122448979589, -0.25102040816326521 -0.01224489795918382, -0.27142857142857135 0.06938775510204076, -0.32448979591836724 0.06938775510204076, -0.29999999999999982 -0.04897959183673484, -0.25510204081632648 -0.11020408163265305, -0.15306122448979576 -0.09387755102040818, -0.16530612244897958 -0.24081632653061247, -0.34897959183673466 -0.25306122448979607, -0.41428571428571415 -0.17142857142857149, -0.50408163265306127 -0.13877551020408174, -0.58571428571428563 -0.22448979591836737, -0.65918367346938767 -0.19183673469387763))'

    tests_r += [tests(seg_wkt, line_wkt, N, 2, face_wkt)] 
    
    # test region with concavity

    #seg_wkt = 'LINESTRING(-140.38800507462761402 -54.16228782222089677, -140.07039179647233595 -54.18960939453533143)'
    line_wkt = 'LINESTRING(-140.38800507462761402 -54.16228782222089677, -140.33848472480769942 -54.09569148970446406, -140.39312786943656874 -54.0393407468059479, -140.34360751961665414 -54.02738755891838451, -140.30774795595397109 -54.06495472085072862, -140.20016926496589349 -54.05812432777211995, -140.16089450476388834 -54.12130546374925189, -140.10795895840468006 -54.08715349835620856, -140.07039179647230753 -54.10252188278307273, -140.13698812898874735 -54.15033463433333338, -140.07039179647233595 -54.18960939453533143)'
    seg_wkt = 'LINESTRING(-140.50924455177292316 -54.26986651320897437, -140.06526900166338123 -54.25620572705175704)'
    
    a_wkt = 'LINESTRING(-140.06526900166338123 -54.25620572705175704, -139.86377240584442916 -54.10935227586167429)'
    b_wkt = 'LINESTRING(-139.86377240584442916 -54.10935227586167429, -139.78522288544041885 -53.90444048350342854)'
    c_wkt = 'LINESTRING(-139.78522288544041885 -53.90444048350342854, -140.05160821550614969 -53.74051104961682768)'
    d_wkt = 'LINESTRING(-140.05160821550614969 -53.74051104961682768, -140.4409406209868223 -53.76100222885264657)'
    e_wkt = 'LINESTRING(-140.4409406209868223 -53.76100222885264657, -140.61511564449133971 -53.89761009042481987)'
    f_wkt = 'LINESTRING(-140.61511564449133971 -53.89761009042481987, -140.70732595105255314 -54.11618266894028295)'
    g_wkt = 'LINESTRING(-140.70732595105255314 -54.11618266894028295, -140.50924455177292316 -54.26986651320897437)'

    sa = loads(a_wkt)
    sb = loads(b_wkt)
    sc = loads(c_wkt)
    sd = loads(d_wkt)
    se = loads(e_wkt)
    sf = loads(f_wkt)
    sg = loads(g_wkt)

    a_wkt = 'LINESTRING(-140.07039179647233595 -54.18960939453533143, -139.92695354182154688 -54.07690790873829201)'
    b_wkt = 'LINESTRING(-139.92695354182154688 -54.07690790873829201, -139.91670795220363743 -53.91468607312134509)'
    c_wkt = 'LINESTRING(-139.91670795220363743 -53.91468607312134509, -140.10795895840468006 -53.815645373481523)'
    d_wkt = 'LINESTRING(-140.10795895840468006 -53.815645373481523, -140.4289874330992518 -53.80539978386360644)'
    e_wkt = 'LINESTRING(-140.4289874330992518 -53.80539978386360644, -140.56730289294108616 -53.9317620558178632)'
    f_wkt = 'LINESTRING(-140.56730289294108616 -53.9317620558178632, -140.59633206352518187 -54.07520031046864517)'
    g_wkt = 'LINESTRING(-140.59633206352518187 -54.07520031046864517, -140.38800507462761402 -54.16228782222089677)'

    ta = loads(a_wkt)
    tb = loads(b_wkt)
    tc = loads(c_wkt)
    td = loads(d_wkt)
    te = loads(e_wkt)
    tf = loads(f_wkt)
    tg = loads(g_wkt)

    objs = [[sa, ta]]
    objs += [[sb, tb]]
    objs += [[sc, tc]]
    objs += [[sd, td]]
    objs += [[se, te]]
    objs += [[sf, tf]]
    objs += [[sg, tg]]
    
    tests_r += [tests(seg_wkt, line_wkt, N, 3, objs)] 
    
    # print (output) results
    
    print('Seg >> Concavity OV : N Tests: ' + str(NTESTS))
    print('mET (sec);MET (sec);AVGET (sec);NV;NC;NS;CHT;')

    for el in tests_r:
        str = ''
        for e in el:
            str += e + ';'

        print(str)

    sys.exit()
else:
    # seg to concavity.
    if type_transf == 1:
        s_exec_time = time.time()
        msegs, nsegs, ch_et, _n_concavities, start_time, step = get_seg_to_concavity_msegs_opt(p, q)
        e_exec_time = time.time()
        
        geoms, num_invalid_geoms, num_complex_geoms = get_in_between_observations2(p, q, msegs, n_obs, start_time, step)
    # region with a concavity.
    elif type_transf == 2:
        line_wkt = 'LINESTRING(-140.38800507462761402 -54.16228782222089677, -140.33848472480769942 -54.09569148970446406, -140.39312786943656874 -54.0393407468059479, -140.34360751961665414 -54.02738755891838451, -140.30774795595397109 -54.06495472085072862, -140.20016926496589349 -54.05812432777211995, -140.16089450476388834 -54.12130546374925189, -140.10795895840468006 -54.08715349835620856, -140.07039179647230753 -54.10252188278307273, -140.13698812898874735 -54.15033463433333338, -140.07039179647233595 -54.18960939453533143)'
        seg_wkt = 'LINESTRING(-140.50924455177292316 -54.26986651320897437, -140.06526900166338123 -54.25620572705175704)'
    
        a_wkt = 'LINESTRING(-140.06526900166338123 -54.25620572705175704, -139.86377240584442916 -54.10935227586167429)'
        b_wkt = 'LINESTRING(-139.86377240584442916 -54.10935227586167429, -139.78522288544041885 -53.90444048350342854)'
        c_wkt = 'LINESTRING(-139.78522288544041885 -53.90444048350342854, -140.05160821550614969 -53.74051104961682768)'
        d_wkt = 'LINESTRING(-140.05160821550614969 -53.74051104961682768, -140.4409406209868223 -53.76100222885264657)'
        e_wkt = 'LINESTRING(-140.4409406209868223 -53.76100222885264657, -140.61511564449133971 -53.89761009042481987)'
        f_wkt = 'LINESTRING(-140.61511564449133971 -53.89761009042481987, -140.70732595105255314 -54.11618266894028295)'
        g_wkt = 'LINESTRING(-140.70732595105255314 -54.11618266894028295, -140.50924455177292316 -54.26986651320897437)'

        sa = loads(a_wkt)
        sb = loads(b_wkt)
        sc = loads(c_wkt)
        sd = loads(d_wkt)
        se = loads(e_wkt)
        sf = loads(f_wkt)
        sg = loads(g_wkt)

        a_wkt = 'LINESTRING(-140.07039179647233595 -54.18960939453533143, -139.92695354182154688 -54.07690790873829201)'
        b_wkt = 'LINESTRING(-139.92695354182154688 -54.07690790873829201, -139.91670795220363743 -53.91468607312134509)'
        c_wkt = 'LINESTRING(-139.91670795220363743 -53.91468607312134509, -140.10795895840468006 -53.815645373481523)'
        d_wkt = 'LINESTRING(-140.10795895840468006 -53.815645373481523, -140.4289874330992518 -53.80539978386360644)'
        e_wkt = 'LINESTRING(-140.4289874330992518 -53.80539978386360644, -140.56730289294108616 -53.9317620558178632)'
        f_wkt = 'LINESTRING(-140.56730289294108616 -53.9317620558178632, -140.59633206352518187 -54.07520031046864517)'
        g_wkt = 'LINESTRING(-140.59633206352518187 -54.07520031046864517, -140.38800507462761402 -54.16228782222089677)'

        ta = loads(a_wkt)
        tb = loads(b_wkt)
        tc = loads(c_wkt)
        td = loads(d_wkt)
        te = loads(e_wkt)
        tf = loads(f_wkt)
        tg = loads(g_wkt)

        objs = [[sa, ta]]
        objs += [[sb, tb]]
        objs += [[sc, tc]]
        objs += [[sd, td]]
        objs += [[se, te]]
        objs += [[sf, tf]]
        objs += [[sg, tg]]
        
        p = loads(seg_wkt)
        q = loads(line_wkt)
    
        s_exec_time = time.time()
        #msegs = reg_to_reg_with_concavity(p, q, objs)
        
        msegs, nsegs, ch_et, _n_concavities, start_time, step = get_seg_to_concavity_msegs_opt(p, q)
        
        for obj in objs:
            a = obj[0]
            b = obj[1]
        
            A = (a.coords[0][0], a.coords[0][1])
            B = (a.coords[1][0], a.coords[1][1])
        
            C = (b.coords[0][0], b.coords[0][1])
            D = (b.coords[1][0], b.coords[1][1])
        
            msegs += [MSeg2(A, A, C, D, 0, 1, 0, True)]
            msegs += [MSeg2(A, B, D, D, 0, 1, 0, True)]
        
        e_exec_time = time.time()
        
        p = loads('Polygon ((-140.61511564449133971 -53.89761009042481987, -140.70732595105255314 -54.11618266894028295, -140.50924455177292316 -54.26986651320897437, -140.06526900166338123 -54.25620572705175704, -139.86377240584442916 -54.10935227586167429, -139.78522288544041885 -53.90444048350342854, -140.05160821550614969 -53.74051104961682768, -140.4409406209868223 -53.76100222885264657, -140.61511564449133971 -53.89761009042481987))')
        q = loads('Polygon ((-140.56730289294108616 -53.9317620558178632, -140.59633206352518187 -54.07520031046864517, -140.38800507462761402 -54.16228782222089677, -140.33848472480769942 -54.09569148970446406, -140.39312786943656874 -54.0393407468059479, -140.34360751961665414 -54.02738755891838451, -140.30774795595397109 -54.06495472085072862, -140.20016926496589349 -54.05812432777211995, -140.16089450476388834 -54.12130546374925189, -140.10795895840468006 -54.08715349835620856, -140.07039179647230753 -54.10252188278307273, -140.13698812898874735 -54.15033463433333338, -140.07039179647233595 -54.18960939453533143, -139.92695354182154688 -54.07690790873829201, -139.91670795220363743 -53.91468607312134509, -140.10795895840468006 -53.815645373481523, -140.4289874330992518 -53.80539978386360644, -140.56730289294108616 -53.9317620558178632))')

        #geoms, num_invalid_geoms, num_complex_geoms = get_in_between_observations_reg_to_reg(p, q, msegs, n_obs)
        geoms, num_invalid_geoms, num_complex_geoms = get_in_between_observations2(p, q, msegs, n_obs, start_time, step, True)
    # concavity with a face
    elif type_transf == 3:
        seg = loads('LineString (-1.25102040816326543 -0.73877551020408183, 0.93673469387755137 -0.59183673469387776)')
        line = loads('LineString (-1.3040816326530611 -0.17551020408163276, -0.83469387755102042 -0.26938775510204094, -0.77346938775510199 -0.06938775510204098, -0.86734693877551017 0.02448979591836731, -0.94489795918367347 -0.01224489795918382, -0.97346938775510194 0.04081632653061218, -0.86734693877551017 0.1020408163265305, -0.78163265306122454 0.02857142857142847, -0.72123661799388505 0.03752937718888849, -0.67959183673469381 0.10612244897959178, -0.74519746809134857 0.19218577327251635, -0.67549317689872757 0.24446399166698207, -0.6020408163265305 0.14285714285714279, -0.51632653061224487 0.20408163265306112, -0.41628034402616809 0.25971180536536792, -0.2790500207406954 0.25317702806605974, -0.18571428571428572 0.19999999999999996, -0.09591836734693882 0.06122448979591832, -0.06775888806306307 -0.0539575050014266, -0.05904585166398535 -0.12366179619404771, -0.07551020408163245 -0.21632653061224505, 0.14897959183673493 -0.33061224489795937, 0.83061224489795915 -0.05306122448979611)')
        face = loads('Polygon ((-0.65918367346938767 -0.19183673469387763, -0.65510204081632639 -0.00408163265306127, -0.59795918367346923 0.07346938775510192, -0.56530612244897949 -0.02448979591836742, -0.5244897959183672 -0.00816326530612255, -0.55306122448979589 0.1020408163265305, -0.42244897959183669 0.18775510204081625, -0.26734693877551008 0.18775510204081625, -0.15306122448979576 0.05306122448979589, -0.25102040816326521 -0.01224489795918382, -0.27142857142857135 0.06938775510204076, -0.32448979591836724 0.06938775510204076, -0.29999999999999982 -0.04897959183673484, -0.25510204081632648 -0.11020408163265305, -0.15306122448979576 -0.09387755102040818, -0.16530612244897958 -0.24081632653061247, -0.34897959183673466 -0.25306122448979607, -0.41428571428571415 -0.17142857142857149, -0.50408163265306127 -0.13877551020408174, -0.58571428571428563 -0.22448979591836737, -0.65918367346938767 -0.19183673469387763))')
    
        s_exec_time = time.time()
    
        # pre procc
        
        l_coords = line.coords
        p_coords = face.exterior.coords[:-1]

        line_coords = []
    
        l_id = 15
        o_id = 7
    
        dxx = 0.00001
    
        r = 0
        while r < len(l_coords):
            if r < l_id:
                line_coords += [(l_coords[r][0], l_coords[r][1])]
            elif r == l_id:
                line_coords += [(l_coords[r][0] - dxx, l_coords[r][1])]
            
                a = o_id
                while a >= 0:
                    if a == o_id:
                        line_coords += [(p_coords[a][0] - dxx, p_coords[a][1])]
                    else:
                        line_coords += [(p_coords[a][0], p_coords[a][1])]
                
                    a -= 1
            
                a = len(p_coords) - 1

                while a >= o_id:
                    if a == o_id:
                        line_coords += [(p_coords[a][0] + dxx, p_coords[a][1])]
                    else:
                        line_coords += [(p_coords[a][0], p_coords[a][1])]
                
                    a -= 1
            
                line_coords += [(l_coords[r][0] + dxx, l_coords[r][1])]
            else:
                line_coords += [(l_coords[r][0], l_coords[r][1])]
    
            r += 1

        nline = LineString(line_coords)
    
        # end pre procc
        
        msegs, nsegs, ch_et, _n_concavities, start_time, step = get_seg_to_concavity_msegs_opt(seg, nline)
        
        #msegs = seg_to_concavity_with_face(seg, line, face)
        e_exec_time = time.time()
        
        geoms, num_invalid_geoms, num_complex_geoms = get_in_between_observations2(seg, line, msegs, n_obs, start_time, step, False, face)
        #geoms, num_invalid_geoms, num_complex_geoms = get_in_between_observations_conc_with_face(seg, line, face, msegs, n_obs)

exec_time = e_exec_time - s_exec_time
exec_time = format(exec_time, precision)

# Output.

print(n_obs)
print(str(num_invalid_geoms) + ',' + str(exec_time))
print(num_complex_geoms)

for g in geoms:
    print(g.wkt)