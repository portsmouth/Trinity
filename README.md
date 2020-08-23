
# Trinity

<a href="https://portsmouth.github.io/Trinity/">Trinity</a> is a programmable WebGL 3d fluid simulator.

<a href='https://portsmouth.github.io/Trinity/?preset="Basic plume"'><img src="./thumbs/Basic-plume.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Plume + sphere collider I"'><img src="./thumbs/Plume-sphere.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Plume + walls"'><img src="./thumbs/Plume-walls.png" width="33%"/></a>
<a href='https://portsmouth.github.io/Trinity/?preset="Nuke III"'><img src="./thumbs/nuke.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Nuke"'><img src="./thumbs/nuke-II.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Moving fireball III"'><img src="./thumbs/fireball.png" width="33%"/></a>
<a href='https://portsmouth.github.io/Trinity/?preset="Dust devil"'><img src="./thumbs/dust-devil.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Dye collision"'><img src="./thumbs/Dye-collision.png" width="33%"/></a><a href='https://portsmouth.github.io/Trinity/?preset="Vortex street"'><img src="./thumbs/vortex-street.png" width="33%"/></a>

## Simulation

Trinity solves the Navier-Stokes equations of fluid/gas dynamics (in the zero viscosity limit).
The simulation is Eulerian on a fixed size grid.


## Solver parameters

- *Nx*, *Ny*, *Nz*: the voxel resolution on each axis.
- *tubeWidth*: the radius of the rendered solution tubes, relative to the box maximum extents



The simulation is specified by the following six programs.

### Common

In this program, arbitrary numeric and color quantities which are used in other programs may be declared in the form:
```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Bind UI parameters to uniforms used in the various programs
//////////////////////////////////////////////////////////////////////////////////////////////////////
uniform float blast_height;               // {"name":"blast_height",     "min":0.1, "max":0.9,    "step":0.001, "default":0.25}
uniform float blast_radius;               // {"name":"blast_radius",     "min":0.0, "max":0.1,    "step":0.001, "default":0.1}
uniform float dust_inflow_rate;           // {"name":"dust_inflow_rate", "min":0.0, "max":10.0,   "step":0.01, "default":1.0}
uniform vec3  dust_absorption;            // {"name":"dust_absorption",  "default":[0.5,0.5,0.5], "scale":1.0}
uniform vec3  dust_scattering;            // {"name":"dust_scattering",  "default":[0.5,0.5,0.5], "scale":1.0}
```
The metadata after the `//` is a JSON object which is used to generate a uniform variable for the shader, driven by a UI slider or color picker (depending on whether "default" is a number or array).
```glsl
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
            in float time,               // time
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

```glsl
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Specify the fluid emission field and phase function
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Approximate map from temperature in Kelvin to blackbody emission color.
// Valid from 1000 to 40000 K (and additionally 0 for pure full white)
vec3 colorTemperatureToRGB(const in float temperature)
{
  mat3 m = (temperature <= 6500.0) ? mat3(vec3(0.0, -2902.1955373783176, -8257.7997278925690),
	                                      vec3(0.0, 1669.5803561666639, 2575.2827530017594),
	                                      vec3(1.0, 1.3302673723350029, 1.8993753891711275)) :
	 								 mat3(vec3(1745.0425298314172, 1216.6168361476490, -8257.7997278925690),
   	                                      vec3(-2666.3474220535695, -2173.1012343082230, 2575.2827530017594),
	                                      vec3(0.55995389139931482, 0.70381203140554553, 1.8993753891711275));
  return mix(clamp(vec3(m[0] / (vec3(clamp(temperature, 1000.0, 40000.0)) + m[1]) + m[2]), vec3(0.0), vec3(1.0)),
             vec3(1.0),
             smoothstep(1000.0, 0.0, temperature));
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


