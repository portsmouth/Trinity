precision highp float;

uniform int N;
uniform float dL;
uniform float timestep;

uniform float buoyancy;
uniform float gravity;
uniform float radiationLoss;
uniform float T0;       // reference temperature for buoyancy
uniform float Tambient; // temperature which generates buoyancy which balances gravity

in vec2 v_texcoord;

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field
uniform sampler2D Pair_sampler; // 1, float pressure field
uniform sampler2D Tair_sampler; // 2, float temperature field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Pair_output;
layout(location = 2) out vec4 Tair_output;


vec3 mapFragToVs(in ivec2 frag)
{
    // map fragment coord in [N*N, N] to continuous position of corresponding voxel center in voxel space
    int iu = frag.x;
    int iv = frag.y;
    int k = iv;
    int j = int(floor(float(iu)/float(N)));
    int i = iu - N*j;
    return vec3(ivec3(i, j, k)) + vec3(0.5);
}

vec4 interp(in sampler2D S, in vec3 wsP)
{
    vec3 vsP = wsP / dL;
    float pY = vsP.y - 0.5;
    int jlo = clamp(int(floor(pY)), 0, N-1); // lower j-slice
    int jhi = clamp(         jlo+1, 0, N-1); // upper j-slice
    float flo = float(jhi) - pY;             // lower j fraction
    float fhi = 1.0 - flo;                   // upper j fraction
    vec2 resolution = vec2(float(N*N), float(N)); // @todo: precompute
    float v = vsP.z / resolution.y;
    float ulo = (vsP.x + float(jlo)*float(N)) / resolution.x;
    float uhi = (vsP.x + float(jhi)*float(N)) / resolution.x;
    vec4 Slo = texture(S, vec2(ulo, v));
    vec4 Shi = texture(S, vec2(uhi, v));
    return flo*Slo + fhi*Shi;
}

void main()
{
    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsX = mapFragToVs(frag);
    vec3 wsX = vsX*dL;

    // Apply semi-Lagrangian advection
    vec3 v0  = texelFetch(Vair_sampler, frag, 0).xyz;
    vec3 wsXp = wsX - v0*timestep; // @todo: implement RK4
    vec3  v = interp(Vair_sampler, wsXp).xyz;
    float P = interp(Pair_sampler, wsXp).x;
    float T = interp(Tair_sampler, wsXp).x;

    // Apply thermal relaxation to ambient temperature due to "radiation loss"
    T *= exp(-radiationLoss);

    // Apply gravity + buoyancy force
    float buoyancy_expansion = buoyancy*(T - T0);
    v.y += timestep * gravity * (1.0 - buoyancy_expansion);

    // Apply solid boundary condition at y=0
    ivec3 vsI = ivec3(vsX);
    if (vsI.y==0) v.y = abs(v.y);

    Vair_output = vec4(v, 0.0);
    Pair_output = vec4(P);
    Tair_output = vec4(T);
}





