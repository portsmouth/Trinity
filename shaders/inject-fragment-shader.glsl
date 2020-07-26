precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform int W;
uniform int H;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;

// Physics
uniform float time;

/////// input buffers ///////
uniform sampler2D Vair_sampler;      // 0, vec3 air velocity field
uniform sampler2D Tair_sampler;      // 2, float air temperature field
uniform sampler2D debris_sampler;    // 3, float debris density field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Tair_output;
layout(location = 2) out vec4 debris_output;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////

float Tambient;                     // {"label":"Tambient",                   "min":0.0, "max":1000.0, "step":0.01, "default":300.0}
float blast_height;                 // {"label":"blast_height",               "min":0.0, "max":1.0,    "step":0.01, "default":0.25}
float blast_radius;                 // {"label":"blast_radius",               "min":0.0, "max":1.0,    "step":0.01, "default":0.1}
float blast_velocity;               // {"label":"blast_velocity",             "min":0.0, "max":100.0,  "step":0.01, "default":50.0}
float blast_temperature_contrast;   // {"label":"blast_temperature_contrast", "min":0.0, "max":1000.0, "step":0.01, "default":100.0}
float debris_inflow_rate;           // {"label":"debris_inflow_rate",         "min":0.0, "max":10.0,   "step":0.01, "default":1.0}
vec3 debris_albedo;                 // {"label":"debris_albedo",              "default":[0.5, 0.5, 0.5], "scale":1.0}

void init(in vec3 wsP,   // world space point of current voxel
          in float time, // time in units of frames
          in vec3 L)     // world-space extents of grid
{
    Tambient = 300.0;
    blast_height = 0.2;
    blast_radius = 0.05;
    blast_velocity = 5.0;
    blast_temperature_contrast = 100.0;
    debris_inflow_rate = 1.0;
    debris_albedo = vec3(0.5, 0.5, 0.5);
}

void inject(in vec3 wsP,        // world space point of current voxel
            in float time,      // time in units of frames
            in vec3 L,          // world-space extents of grid
            inout vec3 dv,      // injected velocity in voxels/frame (defaults to 0)
            inout float T,      // temperature updated in-place for heating/cooling (no update by default)
            inout vec3 drho,    // injected per-channel debris extinction (defaults to 0)
            inout vec3 albedo)  // albedo of the injected debris (defaults to grey)
{
    vec3 blast_center = vec3(0.5*L.x, blast_height*L.y, 0.5*L.z);
    vec3 dir = wsP - blast_center;
    float r = length(dir);
    dir /= r;
    float rt = r/(blast_radius*L.y);
    if (rt <= 1.0)
    {
        // Within blast radius: add velocity and set temperature
        float radial_falloff = max(0.0, 1.0 - rt*rt*(3.0 - 2.0*rt));
        dv = dir * blast_velocity * radial_falloff;
        T = Tambient * (1.0 + blast_temperature_contrast*radial_falloff);
        drho = vec3(1.0) * debris_inflow_rate * radial_falloff;
        albedo = debris_albedo;
    }
}

///////////////////////////////////////////////////////////////////////////////////

bool isSolidCell(in ivec3 vsPi)
{
    vec3 vsP = vec3(float(vsPi.x)+0.5, float(vsPi.y)+0.5,float(vsPi.z)+0.5);
    vec3 wsP = vsP*dL;
    return isSolid(wsP, L);
}

vec3 mapFragToVs(in ivec2 frag)
{
    // map fragment coord in [W, H] to continuous position of corresponding voxel center in voxel space
    int iu = frag.x;
    int iv = frag.y;
    int row = int(floor(float(iv)/float(Nz)));
    int col = int(floor(float(iu)/float(Nx)));
    int i = iu - col*Nx;
    int j = col + row*Ncol;
    int k = iv - row*Nz;
    return vec3(ivec3(i, j, k)) + vec3(0.5);
}

void main()
{
    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsP = mapFragToVs(frag);
    ivec3 vsPi = ivec3(floor(vsP));
    vec3 wsP = vsP*dL;

    if (isSolidCell(vsPi))
    {
        Vair_output   = texelFetch(Vair_sampler, frag, 0);
        Tair_output   = texelFetch(Tair_sampler, frag, 0);
        debris_output = texelFetch(debris_sampler, frag, 0);
    }
    else
    {
        // Get current velocity, temperature, and debris fields:
        vec3  v = texelFetch(Vair_sampler, frag, 0).rgb;
        float T = texelFetch(Tair_sampler, frag, 0).x;
        float E = texelFetch(debris_sampler, frag, 0).x;

        init(wsP, time, L); // @todo: will remove once uniform binding done

        // Inject mass and modify temperature:
        vec3 dv;
        vec3 drho, albedo;
        inject(wsP, time, L,
               dv, T, drho, albedo);

        v += dv;
        E += drho.r; // @todo: mix injected debris extinction/albedo into current medium of voxel

        Vair_output   = vec4(v, 0.0);
        Tair_output   = vec4(T);
        debris_output = vec4(E);
    }
}