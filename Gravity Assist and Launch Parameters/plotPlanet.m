function [globe] = plotPlanet(planetID, position, handle, scaleFactor)
%
% Solar System Sun and Planets Plot
%
%DESCRIPTION:
%This code plot every planet of the Solar System and also the Sun
%For the planetID the list of possible choices is presentes herebelow:
%    1:   Mercury
%    2:   Venus
%    3:   Earth
%    4:   Mars
%    5:   Jupiter
%    6:   Saturn
%    7:   Uranus
%    8:   Neptune
%    9:   Pluto
%    10:  Sun
%--------------------------------------------------------------------------
% INPUTS:
%   planetID     [1x1]       Planet ID                 [-]
%   position     [3x1]       Body Position             [km]
%   handle       [-]         Figure Handle             [-]
%   scaleFactor  [1x1]       Body Scaling Factor       [-]
%--------------------------------------------------------------------------
% OUTPUTS:
%   globe        [chart]     Body 3D Representation    [-]
%--------------------------------------------------------------------------
%
%NOTES:
% - The "handle" input can be fullfilled by giving "gca".
% - The "scaleFactor" acts multiplying the sphere radius by its quantity.
%   Hence if it appears too small in the orbital rapresentation the user
%   can set it. (e.g., for a Earth-Mars Transfer a value of 20 for the
%   Sun and 15 for Earth and Mars gives a good result)
% - In order to have the texture working it is necessary to have the folder
%   named "textures" in the main directory. The code automatically extracts
%   the wanted texture from it.
%
%CALLED FUNCTIONS:
% astroConstants
%
%UPDATES:
% 25/12/2020 - Updated the code description and aspect
%
%REFERENCES:
% (none)
%
%AUTHOR(s):
%Luigi De Maria, 2020
%
%% Pre-Processing
%Figure Handle Setting
if nargin<3
    HAXIS = gca;
elseif ishandle(handle)==0
    msg = ['The figure handle is not valid'];
    error(msg)
else
    try
        HAXIS=gca(handle);
    catch
        HAXIS=handle;  
    end
    hold on
end
%Scale Factor Setting
if nargin<4
    if planetID == 10
        scaleFactor = 1;
    else
        scaleFactor = astroConstants(20+planetID)/astroConstants(3);
    end
end
%Planet Names
planetNames = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', ...
               'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Sun'};
    
%Planet  Radius (w.r.t. Sun, [km])
Rplanet = astroConstants(3)*scaleFactor;
%% Plotting Routine
%Pre-Settings
npanels = 360;      %Number of Globe Panels Around the Equator [deg/panel] = [360/npanels]
alpha = 0.9;        %Alpha (i.e. Transparency Level) of the Globe
erad=Rplanet;       %Equatorial Radius [km]
prad=Rplanet;       %Polar Radius [km]
hold on;
%Axis Settings
axis equal;
axis vis3d;
%3D Mesh Creation
% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x,y,z] = ellipsoid(position(1), position(2), position(3), erad, erad, prad, npanels);
globe = surf(HAXIS, x,y,z,'FaceColor','none','EdgeColor',0.5*[1 1 1], 'HandleVisibility','off');
% RMK.: HandleVisibility=off removes unuseful legends for the plotted globe
%Texturing
cdata=imread(sprintf('./textures/%s.jpg',planetNames{planetID})); % Load Body image for texture map
% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata.
globe.FaceColor = 'texturemap';
globe.CData = flip(cdata); % W/o flip() the earth texture looks upside down
globe.EdgeColor = 'none';
%Environment Settings
globe.FaceLighting = 'gouraud';
globe.AmbientStrength = 0.5;
if planetID ~= 10
    globe.FaceAlpha = 0.9;
else
    globe.FaceAlpha = 1.0;
    globe.FaceLighting = 'none';
end
%Position and Orientation Application
rotate(globe, [0 0 1], 180, position);
end