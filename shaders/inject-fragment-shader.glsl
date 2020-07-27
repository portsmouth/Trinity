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
    init();

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