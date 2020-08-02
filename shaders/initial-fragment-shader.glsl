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
layout(location = 3) out vec4 absorption_output;
layout(location = 4) out vec4 scattering_output;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////

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
    init();

    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsP = mapFragToVs(frag);
    ivec3 vsPi = ivec3(floor(vsP));
    vec3 wsP = vsP*dL;

    vec3 v;
    float T;
    vec3 absorption;
    vec3 scattering;
    initial_conditions(wsP, L,
                       v, T, absorption, scattering);

    Vair_output   = vec4(v, 0.0);
    Pair_output   = vec4(0.0);
    Tair_output   = vec4(T);
    absorption_output = vec4(absorption, 0.0);
    scattering_output = vec4(scattering, 0.0);
}