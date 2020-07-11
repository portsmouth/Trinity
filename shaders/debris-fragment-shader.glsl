precision highp float;

uniform int N;
uniform float dL;
uniform float timestep;

in vec2 v_texcoord;

/////// input buffers ///////
uniform sampler2D debris_sampler; // 0, float debris density field
uniform sampler2D Vair_sampler;   // 1, vec3 velocity field

/////// output buffers ///////
layout(location = 0) out vec4 debris_output;

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

    // Current voxel air velocity
    vec3 v0 = texelFetch(Vair_sampler, frag, 0).rgb;

    // Apply semi-Lagrangian advection
    float h = timestep;
    vec3 wsXp = wsX - v0*h; // @todo: implement RK4
    vec3 debris_p = interp(debris_sampler, wsXp).rgb;
    debris_output = vec4(debris_p, 1.0);
}





