precision highp float;

uniform int N;
uniform float dL;
uniform float timestep;
uniform float expansion;

/////// input buffers ///////
uniform sampler2D Pair_sampler;    // 0, float pressure field
uniform sampler2D Tair_sampler;    // 1, float temperature field
uniform sampler2D divVair_sampler; // 2, float divergence field

/////// output buffers ///////
layout(location = 0) out vec4 Pair_output;


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

    // Get air pressure, temperature and velocity divergence at voxel
    float P    = texelFetch(Pair_sampler, frag, 0).x;
    float T    = texelFetch(Tair_sampler, frag, 0).x;
    float divv = texelFetch(divVair_sampler, frag, 0).x;

    // Introduce local expansion due to heated fluid
    float phi = timestep/float(N) * expansion * T;
    divv -= phi;

    // Get pressure values at local stencil
    float P_xp = texelFetch(Pair_sampler, X_ip, 0).r;
    float P_xn = texelFetch(Pair_sampler, X_in, 0).r;
    float P_yp = texelFetch(Pair_sampler, X_jp, 0).r;
    float P_yn = texelFetch(Pair_sampler, X_jn, 0).r;
    float P_zp = texelFetch(Pair_sampler, X_kp, 0).r;
    float P_zn = texelFetch(Pair_sampler, X_kn, 0).r;

    // Thus compute Jacobi-relaxation solution of Poisson equation:
    //   Laplacian[p] = Divergence[v] / timestep
    float rhs = divv; // / timestep;
    float avgp = (P_xp + P_xn + P_yp + P_yn + P_zp + P_zn) / 6.0;
    float P_relax = avgp - dL*dL*rhs/6.0; // (See Numerical Recipes 19.5.5, extended to 3d)

    Pair_output = vec4(P_relax);
}





