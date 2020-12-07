
# Trinity

<a href="https://portsmouth.github.io/Trinity/">Trinity</a> is a programmable 3D GPU (WebGL) fluid simulator. The fluid simulation logic is fully customizable via GLSL programs coded directly into the browser (via [CodeMirror](https://codemirror.net)). 

<a href='https://portsmouth.github.io/Trinity/?preset="Basic plume"'><img src="./thumbs/Basic-plume.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Plume + sphere collider I"'><img src="./thumbs/Plume-sphere.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Plume + walls"'><img src="./thumbs/Plume-walls.png" width="33%"/></a>
<a href='https://portsmouth.github.io/Trinity/?preset="Nuke III"'><img src="./thumbs/nuke.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Nuke"'><img src="./thumbs/nuke-II.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Moving fireball III"'><img src="./thumbs/fireball.png" width="33%"/></a>
<a href='https://portsmouth.github.io/Trinity/?preset="Dust devil"'><img src="./thumbs/dust-devil.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Dye collision"'><img src="./thumbs/Dye-collision.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Vortex street"'><img src="./thumbs/vortex-street.png" width="33%"/></a>

## Overview

### Simulation and Rendering

Trinity is a WebGL application which solves the [Navier–Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations) of fluid/gas dynamics on GPU, and volume renders the resulting fields. 

Only the core simulation logic is hard-coded, while most of the dynamics is determined by user-written GLSL programs (written in the browser) which specify the injection of fluid velocity and temperature, application of external forces, and the presence of solid obstacles which the fluid collides with. Hot fluid is simulated by injection of a scalar field representing temperature, which is then passively advected and made to affect the dynamics according to buoyancy forces. In general, up to four scalar fields (collectively referred to as "the temperature") may be passively advected and used to drive the dynamics.

For rendering, two color fields representing the extinction (i.e. density) and albedo of an absorbing/scattering/emitting medium, such as dust or ink, are injected and passively advected. These are volume rendered via raymarching, illuminated by a single distant light (the "sun"). The map from the temperature field to emission radiance, to simulate blackbody radiation for example, is provided by the user.

The following 6 user-written GLSL programs (as well as various in-built parameters of the [solver](#solver-parameters) and [renderer](#renderer-parameters)) specify the dynamics and rendering:

  - <a href="#common">Common</a>: declare uniform float and vec3 quantities, and bind them to UI sliders and color pickers.
  - <a href="#initial">Initial</a>: specify initial conditions for velocity, temperature and medium density and albedo.
  - <a href="#inject">Inject</a>: inject velocity, heat or media into the simulation.
  - <a href="#influence">Influence</a>: apply external forces (due to e.g. buoyancy, wind).
  - <a href="#collide">Collide</a>: specify collision geometry via an SDF.
  - <a href="#render">Render</a>: specify how temperature maps to emission, and the phase-function.

### Grid geometry

The simulation is done on a fixed size Eulerian grid.
In all programs, the variable `vec3 wsP` refers to the world space position in coordinates which range from the origin to `vec3 L`, where `L` is in units of voxels.
For example a grid of resolution `(128, 512, 128)` has its lower left corner at `(0, 0, 0)` and its upper right corner at `L=(128.0, 512.0, 128.0)`.
The center of the grid is at `L/2`.

### Technical details

 - As WebGL does not currently support writing to 3D textures from within fragment shaders, the 3D grid has to be represented via 2D textures.
   This is done similarly to the ["flat 3D textures"](https://dl.acm.org/doi/10.5555/844174.844189) of Harris et al (2003). See the commentary [here](https://github.com/portsmouth/Trinity/blob/master/js/solver.js#L144) for a detailed description of the scheme used.
 - Pressure projection is currently rather simplistic and done via Jacobi iteration. The scheme for handling pressure projection in the presence of solid boundaries is taken from ["Real-Time Simulation and Rendering of 3D Fluids"](http://developer.download.nvidia.com/books/gpu_gems_3/samples/gems3_ch30.pdf) by Crane et al (2008).
 - Semi-Lagrangian advection is done via a 4th order Runge-Kutta method.
 - Colliders are currently assumed to be static (i.e. if the SDF is time-dependent, the velocity of the walls will not be transferred to the fluid).
 - Diffusion of the advected terms, as well as fluid viscosity, is ignored (as is common in CFD for graphics).
 - "Neumann" boundary conditions are applied at the edges of the grid (i.e. derivatives are zero at the boundaries, so material flows out freely).
 - Trinity is named after the code name of the [first nuclear weapon test](https://en.wikipedia.org/wiki/Trinity_(nuclear_test)).


## Programs

The simulation is specified in detail by the following six GLSL programs (here with example implementations):

### Common

To allow for real-time interaction with the simulation, there is a simple system for declaring "uniform" float and vec3 variables in the user-code,
and binding them to UI sliders and color-pickers which control them:

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Bind UI parameters to uniforms used in the various programs
// The metadata after the // on each line containing a uniform declaration is a JSON object which is used to
// generate a uniform variable for the shader, which is "bound" to (i.e. driven by) a UI slider or color picker
// (depending on whether "default" is a number or array).
//////////////////////////////////////////////////////////////////////////////////////////////////////

// "Physics"
uniform float gravity;          // {"name":"gravity",          "min":0.0, "max":1.0,   "step":0.01, "default":0.05}
uniform float buoyancy;         // {"name":"buoyancy",         "min":0.0, "max":0.5,   "step":0.01, "default":0.5}
uniform float radiationLoss;    // {"name":"radiationLoss",    "min":0.9, "max":1.0,   "step":0.01,  "default":0.99}

// Blast geometry
uniform float blast_height;     // {"name":"blast_height",     "min":0.1, "max":0.9,   "step":0.01, "default":0.25}
uniform float blast_radius;     // {"name":"blast_radius",     "min":0.0, "max":0.1,   "step":0.01, "default":0.1}
uniform float blast_velocity;   // {"name":"blast_velocity",   "min":0.0, "max":100.0, "step":0.1,  "default":50.0}
uniform float blast_heat_flux;  // {"name":"blast_heat_flux",  "min":0.0, "max":100.0, "step":1.0,  "default":100.0}

// Dust
uniform float dust_inflow_rate; // {"name":"dust_inflow_rate", "min":0.0, "max":10.0,   "step":0.01, "default":1.0}
uniform vec3  dust_absorption;  // {"name":"dust_absorption",  "default":[0.5,0.5,0.5], "scale":1.0}
uniform vec3  dust_scattering;  // {"name":"dust_scattering",  "default":[0.5,0.5,0.5], "scale":1.0}

// Rendering
uniform float TtoKelvin;        // {"name":"TtoKelvin",        "min":0.0, "max":300.0,  "step":0.01, "default":10.0}

float Tambient;

/******************************************************/
/*                 mandatory function                 */
/******************************************************/
void init()
{
    // Any global constants defined here are available in all functions
    Tambient = 1.0;
}
```

### Initial

Specify initial conditions for velocity, temperature and medium density and albedo.

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Specify the initial conditions for the simulation,
// i.e. populate all the relevant fields (velocity, temperature, debris density/albedo) at time 0.0
//////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************************************/
/*                 mandatory function                 */
/******************************************************/
void initial_conditions(in vec3 wsP,               // world space center of current voxel
                        in vec3 L, in float dL,    // world-space extents of grid, and voxel-size
                        inout vec3 v,              // initial velocity
                        inout vec4 T,              // initial temperature
                        inout vec3 medium,         // initial per-channel medium density (extinction)
                        inout vec3 mediumAlbedo)   // initial per-channel medium albedo
{
    v = vec3(0.0);
    T = vec4(Tambient);
    medium = vec3(0.0);
    mediumAlbedo = vec3(0.0);
}
```

### Inject

Inject velocity, heat or media into the simulation.

The temperature field `vec4 T` is in general 4 arbitrary scalar fields which are advected with the flow, which can be used to influence the dynamics in whatever fashion. (For example, one of these fields could be used to represent density of fuel, if simulating combustion). The medium density and albedo fields are "special" in the sense that they determine the absorption and scattering of the rendered media.

The velocity and temperature fields may be injected either by specifying the inflow rate (i.e. supplying a "source" term), or by directly setting the value of the field (i.e like a Dirichlet boundary condition). While only the inflow version is available for the medium field (as using Dirichlet boundary conditions for medium density would be odd).

Note that the simulation time starts at `0.0`, incrementing by one timestep each frame
(looping on reaching max_timesteps).

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Update the velocity, temperature via either:
//  - specification of volumetric inflow/outflow rate due to sources/sinks (vInflow, Tinflow)
//  - modification in-place, i.e. Dirichlet boundary conditions (v, T)
// Also specify the injected medium density inflow rate, and its scattering albedo.
//////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************************************/
/*                 mandatory function                 */
/******************************************************/

void inject(in vec3 wsP,                 // world space center of current voxel
            in float time,               // time (i.e. frame count times timestep value)
            in vec3 L, in float dL,      // world-space extents of grid, and voxel-size
            inout vec3 v,                // modify velocity in-place (defaults to no change)
            inout vec3 vInflow,          // velocity inflow rate (defaults to zero)
            inout vec4 T,                // modify temperature in-place (defaults to no change)
            inout vec4 Tinflow,          // temperature inflow rate (defaults to zero)
            inout vec3 mediumInflow,     // medium density inflow rate (defaults to zero)
            inout vec3 mediumAlbedo)     // medium albedo
{
    vec3 blast_center = vec3(0.5*L.x, blast_height*L.y, 0.5*L.z);
    vec3 dir = wsP - blast_center;
    float r = length(dir);
    dir /= r;
    float rt = r/(blast_radius*L.y);
    if (rt <= 1.0 && time<400.0)
    {
        // Within blast radius: inject velocity and temperature
        float radial_falloff = max(0.0, 1.0 - rt*rt*(3.0 - 2.0*rt));
        vInflow = dir * blast_velocity * radial_falloff;
        Tinflow.r = blast_heat_flux * radial_falloff;

        // Also inject absorbing/scattering "dust"
        vec3 dust_extinction = dust_absorption + dust_scattering;
        mediumInflow = dust_extinction * dust_inflow_rate * radial_falloff;
        mediumAlbedo = dust_scattering / dust_extinction;
    }
    else
    {
        // Apply thermal relaxation due to "radiation loss"
        T.r *= radiationLoss;
    }
}

```

### Influence

Apply external forces (due to e.g. buoyancy, wind), i.e. return the sum of forces to be applied at the given voxel.

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Apply any external forces to the fluid
//////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************************************/
/*                 mandatory function                 */
/******************************************************/

vec3 externalForces(in vec3 wsP,                       // world space center of current voxel
                    in float time,                     // time
                    in vec3 L, in float dL,            // world-space extents of grid, and voxel-size
                    in vec3 v, in float P, in vec4 T,  // velocity, pressure, temperature at current voxel
                    in vec3 medium)                    // medium density at current voxel
{
    // Boussinesq approximation (a la Fedkiw & Stam)
    float densityAvg = (medium.r + medium.g + medium.b)/3.0;
    float buoyancy_force = -densityAvg*gravity + buoyancy*(T.r - Tambient);
    return vec3(0.0, buoyancy_force, 0.0);
}
```

### Collide

Specify collision geometry via an SDF (that is, the volume in which the returned value is negative is taken to be in the interior of a solid obstacle).

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Specify regions which contain impenetrable, static collider material
//////////////////////////////////////////////////////////////////////////////////////////////////////

float sdSphere(vec3 X, in vec3 C, float r) { return length(X-C) - r; }

/******************************************************/
/*                 mandatory function                 */
/******************************************************/

float collisionSDF(in vec3 wsP,            // world space center of current voxel
                   in float time,          // time
                   in vec3 L, in float dL) // world-space extents of grid, and voxel-size
{
    // Return SDF of the collider surface.
    // (where the interior with SDF < 0.0 is a solid obstacle)
    return sdSphere(wsP, vec3(L.x/2.0, L.y/3.0, L.z/2.0), 0.5*L.x*collider_radius);
}
```

### Render

Specify how temperature maps to emission, and phase-function.

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Specify the fluid emission field and phase function
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Approximate map from temperature in Kelvin to blackbody emission color.
// Valid from 1000 to 40000 K (and additionally 0 for pure full white)
vec3 colorTemperatureToRGB(const in float temperature)
{
    // (implementation omitted here)
}

/******************************************************/
/*                 mandatory functions                */
/******************************************************/

// Specify how the temperature is mapped to the local emission radiance
vec3 temperatureToEmission(in vec4 T)
{
    vec3 emission = colorTemperatureToRGB(T.r * TtoKelvin) * pow(T.r/100.0, 4.0);
    return emission;
}

// Optionally remap the medium density (extinction) and albedo
void mediumRemap(inout vec3 medium,
                 inout vec3 mediumAlbedo)
{}

// Specify phase function of medium
float phaseFunction(float mu,         // cosine of angle between incident and scattered ray
                    float anisotropy) // anisotropy coefficient
{
    const float pi = 3.141592653589793;
    float g = anisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*pi)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}
```

## Parameters

For reference, these are the parameters of the solver and volume renderer:

### Solver parameters

The fluid solver has the following parameters:

- *Nx*, *Ny*, *Nz*: the voxel resolution on each axis.
- *NprojSteps*: the number of Jacobi iterations for the pressure projection step
- *max_timesteps*: the maximum timestep count, after which the simulation loops
- *vorticity_scale*: controls the amount of "vorticity confinement" applied, which enhances detail at the expense of simulation stability and correctness (see [Fedkiw & Stam](https://dl.acm.org/doi/10.1145/383259.383260)).
- *expansion*: simulates local fluid expansion due to heating.
- *timestep*: timestep value, normally fixed at 1.0

### Renderer parameters

Volume rendering parameters:

- *extinctionScale*: scale the medium density (extinction) value
- *emissionScale*: scale the emission radiance
- *anisotropy*: set the anisotropy of the phase-function
- *exposure*: tonemapping exposure
- *gamma*: tonemapping gamma
- *saturation*: tonemapping color saturation
- *sunLatitude*: latitude angle of directional light
- *sunLongitude*: longitude angle of directional light
- *sunPower*: power of directional light
- *sunColor*: color of directional light
- *skyColor*: color of (uniform) sky dome light
- *Nraymarch*: number of raymarch steps
- *max_spp*: maximum number of samples to average the jittered raymarch over
- *spp_per_frame*: number of spp to run for each frame
- *show_bounds*: show wireframe grid bounds
- *colliderDiffuse*: collision SDF diffuse reflection color
- *colliderSpec*: collision SDF specular reflection color
- *colliderRoughness*: collision SDF specular roughness

### Load/Save scene

 - *save scene*: save scene out to a JSON file (filename auto-generated with a timestamp)
 - *load scene*: load a previously saved JSON scene file

### Keys

  - Space bar: toggle the simulation play/pause state
  - ESC key: restart the simulation from time 0.0
  - F11 key: go fullscreen
  - O key: dump JSON application state to console, used to generate new presets.
  - F key: move camera to the standard orientation
  - H key: toggle the menus and HUD
  - P key: export current framebuffer to PNG image


