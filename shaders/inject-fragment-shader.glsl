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
uniform float timestep;

/////// input buffers ///////
uniform sampler2D Vair_sampler;       // 0, vec3 air velocity field
uniform sampler2D Tair_sampler;       // 1, float air temperature field
uniform sampler2D absorption_sampler; // 2, vec3 absorption field
uniform sampler2D scattering_sampler; // 3, vec3 scattering field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Tair_output;
layout(location = 2) out vec4 absorption_output;
layout(location = 3) out vec4 scattering_output;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////

bool isSolidCell(in ivec3 vsPi)
{
    vec3 vsP = vec3(float(vsPi.x)+0.5, float(vsPi.y)+0.5,float(vsPi.z)+0.5);
    vec3 wsP = vsP*dL;
    return collisionSDF(wsP, time, L, dL) < 0.0;
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
        absorption_output = texelFetch(absorption_sampler, frag, 0);
        scattering_output = texelFetch(scattering_sampler, frag, 0);
    }
    else
    {
        // Get current velocity, temperature, and debris fields:
        vec3 v = texelFetch(Vair_sampler, frag, 0).rgb;
        vec4 T = texelFetch(Tair_sampler, frag, 0);
        vec3 absorption = texelFetch(absorption_sampler, frag, 0).rgb;
        vec3 scattering = texelFetch(scattering_sampler, frag, 0).rgb;

        // Inject mass and modify temperature:
        vec3 vInflow      = vec3(0.0);
        vec4 Tinflow      = vec4(0.0);
        vec3 mediumInflow = vec3(0.0);
        vec3 mediumAlbedo = vec3(0.5);

        inject(wsP, time, L, dL,
               v, vInflow,
               T, Tinflow,
               mediumInflow, mediumAlbedo);

        v += vInflow * timestep;
        T += Tinflow * timestep;

        vec3 scatteringInflow = mediumAlbedo * mediumInflow;
        vec3 absorptionInflow = mediumInflow - scatteringInflow;
        absorption += absorptionInflow * timestep;
        scattering += scatteringInflow * timestep;

        Vair_output       = vec4(v, 0.0);
        Tair_output       = T;
        absorption_output = vec4(absorption, 0.0);
        scattering_output = vec4(scattering, 0.0);
    }
}