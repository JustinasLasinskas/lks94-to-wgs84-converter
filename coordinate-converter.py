import math

def wgs2lks(latitude, longitude):
    lat_rad = math.radians(latitude)
    lon_rad = math.radians(longitude)

    a = 6378137
    f = 1 / 298.257223563
    b = a * (1 - f)
    e2 = (a * a - b * b) / (a * a)
    n = (a - b) / (a + b)

    G = a * (1 - n) * (1 - n * n) * (1 + (9 / 4) * n * n + (255 / 64) * math.pow(n, 4)) * (math.pi / 180)

    w = math.radians(longitude - 24)
    t = math.tan(lat_rad)
    rho = a * (1 - e2) / math.pow(1 - (e2 * math.sin(lat_rad) * math.sin(lat_rad)), (3 / 2))
    nu = a / math.sqrt(1 - (e2 * math.sin(lat_rad) * math.sin(lat_rad)))

    psi = nu / rho
    cos_lat = math.cos(lat_rad)

    A0 = 1 - (e2 / 4) - (3 * e2 * e2 / 64) - (5 * math.pow(e2, 3) / 256)
    A2 = (3 / 8) * (e2 + (e2 * e2 / 4) + (15 * math.pow(e2, 3) / 128))
    A4 = (15 / 256) * (e2 * e2 + (3 * math.pow(e2, 3) / 4))
    A6 = 35 * math.pow(e2, 3) / 3072
    m = a * ((A0 * lat_rad) - (A2 * math.sin(2 * lat_rad)) + (A4 * math.sin(4 * lat_rad)) - (A6 * math.sin(6 * lat_rad)))

    eterm1 = (w * w / 6) * cos_lat * cos_lat * (psi - t * t)
    eterm2 = (math.pow(w, 4) / 120) * math.pow(cos_lat, 4) * (4 * math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (1 + 8 * t * t) - psi * 2 * t * t + math.pow(t, 4))
    eterm3 = (math.pow(w, 6) / 5040) * math.pow(cos_lat, 6) * (61 - 479 * t * t + 179 * math.pow(t, 4) - math.pow(t, 6))
    dE = 0.9998 * nu * w * cos_lat * (1 + eterm1 + eterm2 + eterm3)
    long = round(500000 + (dE / 1), 3)

    nterm1 = (w * w / 2) * nu * math.sin(lat_rad) * cos_lat
    nterm2 = (math.pow(w, 4) / 24) * nu * math.sin(lat_rad) * math.pow(cos_lat, 3) * (4 * psi * psi + psi - t * t)
    nterm3 = (math.pow(w, 6) / 720) * nu * math.sin(lat_rad) * math.pow(cos_lat, 5) * (8 * math.pow(psi, 4) * (11 - 24 * t * t) - 28 * math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (1 - 32 * t * t) - psi * 2 * t * t + math.pow(t, 4))
    nterm4 = (math.pow(w, 8) / 40320) * nu * math.sin(lat_rad) * math.pow(cos_lat, 7) * (1385 - 3111 * t * t + 543 * math.pow(t, 4) - math.pow(t, 6))
    dN = 0.9998 * (m + nterm1 + nterm2 + nterm3 + nterm4)
    lat = round(0 + (dN / 1), 3)

    return lat, long 


def lks2wgs(y, x):
    a = 6378137.0
    f = 1 / 298.257222101
    b = a * (1 - f)
    e = math.sqrt((a**2 - b**2) / a**2)
    n = (a - b) / (a + b)
    A1 = a / (1 + n) * (1 + 1/4*n**2 + 1/64*n**4)
    h1 = 1/2*n - 2/3*n**2 + 37/96*n**3 - 1/360*n**4
    h2 = 1/48*n**2 + 1/15*n**3 - 437/1440*n**4
    h3 = 17/480*n**3 - 37/840*n**4
    h4 = 439/161280*n**4
    lon0 = 24 / 180 * math.pi
    k0 = 0.9998
    x0 = 500000
    dy = y
    dx = x - x0
    M = dy / k0
    mu = M / A1
    fp = mu + h1*math.sin(2*mu) + h2*math.sin(4*mu) + h3*math.sin(6*mu) + h4*math.sin(8*mu)
    e1 = e / math.sqrt(1 - e**2)
    sinfp = math.sin(fp)
    tanfp = math.tan(fp)
    cosfp = math.cos(fp)
    t = tanfp * (1 - e**2 * sinfp**2)
    N = a / math.sqrt(1 - e**2 * sinfp**2)
    R = N * (1 - e**2) / (1 - e**2 * sinfp**2)
    D = dx / (N * k0)
    Q1 = N * tanfp / R
    Q2 = D**2 / 2
    Q3 = (5 + 3*t + 10*t**2 - 4*t**3 - 9*e**2) * D**4 / 24
    Q4 = (61 + 90*t + 298*t**2 + 45*t**3 - 252*e**2 - 3*e**2**2) * D**6 / 720
    lat = fp - Q1 * (Q2 - Q3 + Q4)
    Q5 = D
    Q6 = (1 + 2*t + t**2) * D**3 / 6
    Q7 = (5 - 2*t + 28*t**2 - 3*t**3 + 8*e**2 + 24*t**4) * D**5 / 120
    lon = lon0 + (Q5 - Q6 + Q7) / cosfp
    lat = lat * 180 / math.pi
    lon = lon * 180 / math.pi

    return lat, lon