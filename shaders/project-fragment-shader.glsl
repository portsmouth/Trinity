precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;

// Physics
uniform float timestep;
uniform float expansion;

/////// input buffers ///////
uniform sampler2D Pair_sampler;    // 0, float pressure field
uniform sampler2D Tair_sampler;    // 1, float temperature field
uniform sampler2D divVair_sampler; // 2, float divergence field

/////// output buffers ///////
layout(location = 0) out vec4 Pair_output;

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

ivec2 mapVsToFrag(in ivec3 vsP)
{
    // map integer voxel space coords to the corresponding fragment coords in [W, H]
    int i = vsP.x;
    int j = vsP.y;
    int k = vsP.z;
    int row = int(floor(float(j)/float(Ncol)));
    int col = j - row*Ncol;
    int iu = col*Nx + i;
    int iv = row*Nz + k;
    return ivec2(iu, iv);
}

bool isSolidCell(in ivec3 vsPi)
{
    vec3 vsP = vec3(float(vsPi.x)+0.5, float(vsPi.y)+0.5,float(vsPi.z)+0.5);
    vec3 wsP = vsP*dL;
    return isSolid(wsP, L);
}

void main()
{
    init();

    ivec2 frag = ivec2(gl_FragCoord.xy);
    ivec3 vsX = ivec3(floor(mapFragToVs(frag)));
    int ix = vsX.x;
    int iy = vsX.y;
    int iz = vsX.z;

    // Apply Neumann boundary conditions at grid boundaries
    ivec3 X_ip = ivec3(min(ix+1, Nx-1), iy, iz);
    ivec3 X_in = ivec3(max(ix-1,    0), iy, iz);
    ivec3 X_jp = ivec3(ix, min(iy+1, Ny-1), iz);
    ivec3 X_jn = ivec3(ix, max(iy-1,    0), iz);
    ivec3 X_kp = ivec3(ix, iy, min(iz+1, Nz-1));
    ivec3 X_kn = ivec3(ix, iy, max(iz-1,    0));

    // Get air pressure, temperature and velocity divergence at voxel
    float P     = texelFetch(Pair_sampler, frag, 0).x;
    float T     = texelFetch(Tair_sampler, frag, 0).x;
    float div_v = texelFetch(divVair_sampler, frag, 0).x;

    // Introduce local expansion due to heated fluid
    div_v -= timestep * expansion * T;

    // Get pressure values at local stencil
    float P_xp = texelFetch(Pair_sampler, mapVsToFrag(X_ip), 0).r;
    float P_xn = texelFetch(Pair_sampler, mapVsToFrag(X_in), 0).r;
    float P_yp = texelFetch(Pair_sampler, mapVsToFrag(X_jp), 0).r;
    float P_yn = texelFetch(Pair_sampler, mapVsToFrag(X_jn), 0).r;
    float P_zp = texelFetch(Pair_sampler, mapVsToFrag(X_kp), 0).r;
    float P_zn = texelFetch(Pair_sampler, mapVsToFrag(X_kn), 0).r;

    // Apply pressure Neumann boundary-condition for solid cells:
    if (isSolidCell(X_ip)) P_xp = P;
    if (isSolidCell(X_in)) P_xn = P;
    if (isSolidCell(X_jp)) P_yp = P;
    if (isSolidCell(X_jn)) P_yn = P;
    if (isSolidCell(X_kp)) P_zp = P;
    if (isSolidCell(X_kn)) P_zn = P;

    // Thus compute Jacobi-relaxation solution of Poisson equation:
    //   Laplacian[p] = Divergence[v]
    float avgp = (P_xp + P_xn + P_yp + P_yn + P_zp + P_zn) / 6.0;

    // For the pressure projection formula, see:
    //  - "Fast Fluid Dynamics Simulation on the GPU", Harris
    //  - "Real-Time Simulation and Rendering of 3D Fluids", Crane et. al
    //  - Numerical Recipes 19.5.5, extended to 3d
    float P_relax = avgp - dL*dL*div_v/6.0;

    Pair_output = vec4(P_relax);
}





