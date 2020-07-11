precision highp float;

uniform int N;
uniform float dL;

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field

/////// output buffers ///////
layout(location = 0) out vec4 divVair_output;


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

ivec2 mapVsToFrag(in ivec3 vsP)
{
    // map integer voxel space coords to the corresponding fragment coords
    int i = vsP.x;
    int j = vsP.y;
    int k = vsP.z;
    int ui = N*j + i;
    int vi = k;
    return ivec2(ui, vi);
}

void main()
{
    ivec2 frag = ivec2(gl_FragCoord.xy);
    ivec3 vsX = ivec3(mapFragToVs(frag));
    int ix = vsX.x;
    int iy = vsX.y;
    int iz = vsX.z;

    // Apply Neumann boundary conditions
    ivec2 X_ip = mapVsToFrag(ivec3(min(ix+1, N-1), iy, iz));
    ivec2 X_in = mapVsToFrag(ivec3(max(ix-1, 0),   iy, iz));
    ivec2 X_jp = mapVsToFrag(ivec3(ix, min(iy+1, N-1), iz));
    ivec2 X_jn = mapVsToFrag(ivec3(ix, max(iy-1, 0),   iz));
    ivec2 X_kp = mapVsToFrag(ivec3(ix, iy, min(iz+1, N-1)));
    ivec2 X_kn = mapVsToFrag(ivec3(ix, iy, max(iz-1, 0)  ));

    vec4 V_xp = texelFetch(Vair_sampler, X_ip, 0);
    vec4 V_xn = texelFetch(Vair_sampler, X_in, 0);
    vec4 V_yp = texelFetch(Vair_sampler, X_jp, 0);
    vec4 V_yn = texelFetch(Vair_sampler, X_jn, 0);
    vec4 V_zp = texelFetch(Vair_sampler, X_kp, 0);
    vec4 V_zn = texelFetch(Vair_sampler, X_kn, 0);

    float divVair = 0.5 * (V_xp.x - V_xn.x + V_yp.y - V_yn.y + V_zp.z - V_zn.z) / dL;
    divVair_output = vec4(divVair);
}





