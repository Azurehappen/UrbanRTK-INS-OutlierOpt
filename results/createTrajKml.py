import numpy as np
from scipy.io import loadmat
from simplekml import Kml, Color, IconStyle, AltitudeMode, LookAt
from pyproj import Proj, Transformer

def rgb_to_abgr(r, g, b, alpha=1.0):
    """Converts RGB color values (scaled 0-1) to KML (ABGR) hex format."""
    a = int(alpha * 255)
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    return Color.rgb(b, g, r, a)  # Using simplekml's Color.rgb to create ABGR formatted color

def ecef_to_latlon(ecef):
    # Define the projection for ECEF and geographic coordinates
    transformer = Transformer.from_proj(
        proj_from=Proj(proj='geocent', ellps='WGS84', datum='WGS84'),
        proj_to=Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    )

    # Perform the transformation
    lon, lat, alt = transformer.transform(ecef[0], ecef[1], ecef[2])
    return lat, lon, alt

def create_kml_from_mat(mat_file, kml_file):
    data = loadmat(mat_file)
    pos_ecef = data['output']['pos_ecef'][0,0]
    hor_err = data['output']['hor_err'][0,0].flatten()

    kml = Kml()

    first_lat, first_lon, first_alt = ecef_to_latlon(pos_ecef[:, 0])
    # Set the LookAt to a specific point
    lookat = LookAt(
        longitude=first_lon,  # Longitude for Austin, TX (not San Francisco)
        latitude=first_lat,    # Latitude for Austin, TX
        altitude=1000,         # Altitude in meters
        range=500,             # Distance in meters from the point
        tilt=45,               # Angle from the vertical
        heading=0              # Direction of the camera
    )
    kml.document.lookat = lookat

    # Define MATLAB-like colors, adjusted if necessary for intensity
    colors = [
        rgb_to_abgr(0, 1, 0),  # Green
        Color.cyan,  # Cyan
        rgb_to_abgr(0.5, 0, 0.5),  # Purple
        rgb_to_abgr(1, 0, 0)  # Red
    ]
    thresholds = [1, 3, 20]

    for i in range(pos_ecef.shape[1]):
        lat, lon, alt = ecef_to_latlon(pos_ecef[:, i])
        # Determine the color based on the horizontal error
        error = hor_err[i]
        if error < thresholds[0]:
            color = colors[0] # < 1 m, green
        elif error < thresholds[1]:
            color = colors[1] # 1 - 3 m, blue
        elif error < thresholds[2]:
            color = colors[2] # 3 - 20 m, purple
        else:
            color = colors[3] # > 20 m, red

        pnt = kml.newpoint(coords=[(lon, lat, alt)])
        pnt.style.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/shapes/road_shield3.png'  # Use a small circle icon
        pnt.style.iconstyle.scale = 0.5  # Scale down the icon to make it less prominent
        pnt.style.iconstyle.color = color

    kml.save(kml_file)

# Usage example
# You may have to define the correct path of the mat file.
create_kml_from_mat('RAPS-INS-RTK.mat', 'raps_trajectory.kml')