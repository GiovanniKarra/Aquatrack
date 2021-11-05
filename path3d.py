import numpy as np
import scipy.interpolate as spip


def path_points(points, steps=None):
    """
    Calcule un chemin courbe à partir de points de passage.
    Les courbes sont des splines cubiques.
    
    Paramètres:
        points: array[3,N]
            N points de passage en 3 dimensions
        steps: int
            nombre de points à générer, par défaut 10 * N
    Retourne:
        XPath: array[3,steps]
            coordonnées des points sur le chemin en 3 dimensions
    """

    # paramètre: racine carrée de la corde accumulée:
    deltaTau = np.sqrt(np.sum(np.diff(points) ** 2, axis=0))
    tauPoints = np.hstack((0, np.cumsum(deltaTau)))
    tauEnd = tauPoints[-1]

    # Interpoler avec des spline cubiques: 
    spline = spip.splprep(points, u=tauPoints, s=0)[0]

    # Echantillonner à intervalles réguliers:
    if not steps: steps = 10 * points.shape[1]
    tau = np.linspace(0, tauEnd, steps)
    XPath = np.array(spip.splev(tau, spline))

    return XPath


def path_vectors(XPath):
    """
    Calcule la distance curvilinéaire, le vecteur tangent et
    le vecteur de courbure le long du chemin.
    
    Paramètres:
        XPath: array[3,N]
            coordonnées de N points en 3 dimensions
    Retourne: (sPath, TPath, CPath)
        sPath: array[N]
            distance curvilinéaire des points
        TPath: array[3,N]
            vecteur tangent aux points (|T| = 1)
        CPath: array[3,N]
            vecteur de courbure aux points (|C| = courbure = 1/rayon)
    """
    ds = np.sqrt(np.sum(np.diff(XPath) ** 2, axis=0))
    s = np.hstack((0, np.cumsum(ds)))

    dX_ds = np.gradient(XPath, s, axis=1, edge_order=2)
    d2X_ds = np.gradient(dX_ds, s, axis=1, edge_order=2)

    return s, dX_ds, d2X_ds


def path(points, steps=None):
    """
    Calcule les éléments d'un chemin courbe à partir de points de passage.
    
    Paramètres:
        points: array[3,N]
            N points de passage en 3 dimensions
        steps: int
            nombre de points à générer, par défaut 10 * N
    Retourne: (sPath, TPath, CPath)
        XPath: array[3,steps]
            coordonnées des points sur le chemin en 3 dimensions
        sPath: array[N]
            distance curvilinéaire des points
        TPath: array[3,N]
            vecteur tangent aux points (|T| = 1)
        CPath: array[3,N]
            vecteur de courbure aux points (|C| = courbure = 1/rayon)
    """
    XPath = path_points(points, steps)
    sPath, TPath, CPath = path_vectors(XPath)
    return sPath, XPath, TPath, CPath


def ainterp(x, xp, ap, **kwargs):
    """
    Interpolation en x sur plusieurs fonctions xp -> ap[k].
    
    Paramètres:
        x: float ou array[M]
            abscisses x_m où les fonctions sont évaluées
        xp: array[N]
            abscisses x_n des points des fonctions
        ap: array[K, N]
            ordonnées f_k(x_n) des points des fonctions
            
    Résultat:
        a: array[K] ou array[K, M]
            ordonnées évaluées f_k(x_m)
    """
    a = np.array([np.interp(x, xp, fp, **kwargs) for fp in ap])
    return a


def path_at(s, path, **kwargs):
    """
    Retourne les éléments en un point donné d'un chemin,
    par interpolation des données du chemin.
    
    Paramètres:
        s: float
            distance curviligne du point
        path: (sPath, XPath, TPath, CPath)
            éléments du chemin, tels que retournés par path().
            
    Retourne: X, T, C
        X: array[3]
            coordonnées du point
        T: array[3]
            vecteur tangent au point
        C: array[3]
            vecteur de courbure au point
    """
    sPath, XPath, TPath, CPath = path
    X = ainterp(s, sPath, XPath, **kwargs)
    T = ainterp(s, sPath, TPath, **kwargs)
    C = ainterp(s, sPath, CPath, **kwargs)
    return X, T, C


def tilt_vectors(T, beta=0.):
    """
    Retourne les vecteurs d'inclinaison le long d'un chemin.
    
    Paramètres:
        T: array[3]
            vecteur tangent au chemin
        beta: float, défaut=0.
            angle d'inclinaison par rapport à l'horizontale
            (> 0 = penche à droite dans le sens de T, en radians)
            
    Retourne: B, N
        B: array[3]
            vecteur unitaire normal parallèle à l'inclinaison
            (à gauche dans le sens de T)
        N: array[3]
            vecteur unitaire normal perpendiculaire à l'inclinaison
            (vers z > 0)
    """
    T /= np.sqrt(np.sum(T ** 2, axis=0))  # unit
    Z = np.array((0, 0, 1))

    D = np.cross(Z, T, axis=0)
    D /= np.sqrt(np.sum(D ** 2, axis=0))  # unit normale horizontale

    U = np.cross(T, D, axis=0)  # unit normale verticale

    B = D * np.cos(beta) + U * np.sin(beta)  # unit normale parallèle
    N = np.cross(T, B, axis=0)  # unit normale perpendiculaire
    return B, N


def looping_points(steps=21):
    """
    Calcule les points de passage pour un chemin simple avec
    un looping.  Le chemin a pour équations:
        
        x = t*(t-1)**2 + 0.5*t
        y = 0.5*t*(t-1)**2
        z = 2*t**2*(t-0.9)**2
        
    pour -1 ≤ t ≤ 1
    
    Paramètres:
        steps: int (défaut 21)
            le nombre de points
    Retourne: array[3, steps]
        coordonnées des point de passage
    """
    tPoints = np.linspace(-1, 1, 21)
    Xpoints = np.vstack((
        0.5 * tPoints + tPoints * (tPoints - 1) * (tPoints + 1),
        0.5 * tPoints * (tPoints - 1) * (tPoints + 1),
        2 * tPoints * tPoints * (tPoints - 0.9) * (tPoints + 0.9)
    ))
    return Xpoints


#### Code de test ####

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # points de passage
    xyzPoints = looping_points()

    # chemin et vecteurs
    xyzPath = path_points(xyzPoints)
    sPath, TPath, CPath = path_vectors(xyzPath)

    # points jalons à afficher sur le graphique
    length = sPath[-1]
    sMarks = np.linspace(0, length, 30)
    xyzMarks = np.array([np.interp(sMarks, sPath, uPath) for uPath in xyzPath])
    TMarks = np.array([np.interp(sMarks, sPath, uPath) for uPath in TPath])
    CMarks = np.array([np.interp(sMarks, sPath, uPath) for uPath in CPath])

    # graphique 3D : chemin et vecteurs
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect(np.ptp(xyzPath, axis=1))
    ax.plot(xyzPoints[0], xyzPoints[1], xyzPoints[2], 'bo', label='points')
    ax.plot(xyzPath[0], xyzPath[1], xyzPath[2], 'k.', ms=0.5)
    ax.plot(xyzMarks[0], xyzMarks[1], xyzMarks[2], 'r.', ms=2)
    scale = 0.1
    ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
              scale * TMarks[0], scale * TMarks[1], scale * TMarks[2],
              color='r', linewidth=0.5, label='T')
    ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
              scale * CMarks[0], scale * CMarks[1], scale * CMarks[2],
              color='g', linewidth=0.5, label='C')
    ax.legend()
    plt.show()

    # Vecteurs normaux
    # tilt = 2*xyzMarks[1]
    tilt = 0.
    BMarks, NMarks = tilt_vectors(TMarks, tilt)

    # graphique 3D : vecteurs normaux
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect(np.ptp(xyzPath, axis=1))
    ax.plot(xyzPath[0], xyzPath[1], xyzPath[2], 'k.', ms=0.5)
    # ax.plot(xyzMarks[0],xyzMarks[1],xyzMarks[2],'r.', ms=2)
    scale = 0.1
    ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
              scale * TMarks[0], scale * TMarks[1], scale * TMarks[2],
              color='r', linewidth=0.5, label='T')
    ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
              scale * BMarks[0], scale * BMarks[1], scale * BMarks[2],
              color='b', linewidth=0.5, label='B')
    ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
              scale * NMarks[0], scale * NMarks[1], scale * NMarks[2],
              color='g', linewidth=0.5, label='N')
    plt.legend()
    plt.show()
