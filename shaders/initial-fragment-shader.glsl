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

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Pair_output;
layout(location = 2) out vec4 Tair_output;
layout(location = 3) out vec4 debris_output;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////

vec3 blast_center;
float blast_radius;
float blast_velocity;
float blast_temperature_contrast;
float T_ambient;
vec3 debris_density;
vec3 debris_albedo;

// Specify velocity, temperature, and debris density/albedo at time=0
void initial_conditions(in vec3 wsP,                // world space point of current voxel
                        in vec3 L,                  // world-space extents of grid
                        inout vec3 velocity,        // initial velocity
                        inout float temperature,    // initial temperature
                        inout vec3 density,         // initial per-channel debris extinction
                        inout vec3 albedo)          // initial per-channel debris albedo
{
    velocity = vec3(0.0);

    float T_ambient = 300.0;
    temperature = T_ambient;

    density = vec3(0.0);
    albedo = vec3(0.0);
}


///////////////////////////////////////////////////////////////////////////////////


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

    vec3 v;
    float T;
    vec3 E; // debris extinction
    vec3 A; // debris albedo
    initial_conditions(wsP, L,
                       v, T, E, A);

    Vair_output   = vec4(v, 0.0);
    Pair_output   = vec4(0.0);
    Tair_output   = vec4(T);
    debris_output = vec4(E, 0.0);
}