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

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field
uniform sampler2D Pair_sampler; // 1, float pressure field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;

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

_USER_CODE_

bool isSolidCell(in ivec3 vsPi)
{
    vec3 vsP = vec3(float(vsPi.x)+0.5, float(vsPi.y)+0.5,float(vsPi.z)+0.5);
    vec3 wsP = vsP*dL;
    return isSolid(wsP, L);
}

void main()
{
    // Setup local stencil:
    ivec2 frag = ivec2(gl_FragCoord.xy);
    ivec3 vsXi = ivec3(floor(mapFragToVs(frag)));
    int ix = vsXi.x;
    int iy = vsXi.y;
    int iz = vsXi.z;

    // Apply Neumann boundary conditions at grid boundaries
    ivec3 X_ip = ivec3(min(ix+1, Nx-1), iy, iz);
    ivec3 X_in = ivec3(max(ix-1,    0), iy, iz);
    ivec3 X_jp = ivec3(ix, min(iy+1, Ny-1), iz);
    ivec3 X_jn = ivec3(ix, max(iy-1,    0), iz);
    ivec3 X_kp = ivec3(ix, iy, min(iz+1, Nz-1));
    ivec3 X_kn = ivec3(ix, iy, max(iz-1,    0));

    // air velocity at voxel
    vec3 v = texelFetch(Vair_sampler, frag, 0).rgb;

    // Get pressure values at local stencil
    float  P   = texelFetch(Pair_sampler, frag, 0).r;
    float P_xp = texelFetch(Pair_sampler, mapVsToFrag(X_ip), 0).r;
    float P_xn = texelFetch(Pair_sampler, mapVsToFrag(X_in), 0).r;
    float P_yp = texelFetch(Pair_sampler, mapVsToFrag(X_jp), 0).r;
    float P_yn = texelFetch(Pair_sampler, mapVsToFrag(X_jn), 0).r;
    float P_zp = texelFetch(Pair_sampler, mapVsToFrag(X_kp), 0).r;
    float P_zn = texelFetch(Pair_sampler, mapVsToFrag(X_kn), 0).r;

    // Apply pressure Neumann boundary-condition for solid cells:
    vec3 vMask = vec3(1.0, 1.0, 1.0);
    if (isSolidCell(X_ip)) { P_xp = P; vMask.x = 0.0; }
    if (isSolidCell(X_in)) { P_xn = P; vMask.x = 0.0; }
    if (isSolidCell(X_jp)) { P_yp = P; vMask.y = 0.0; }
    if (isSolidCell(X_jn)) { P_yn = P; vMask.y = 0.0; }
    if (isSolidCell(X_kp)) { P_zp = P; vMask.z = 0.0; }
    if (isSolidCell(X_kn)) { P_zn = P; vMask.z = 0.0; }

    // Update air velocity accordingly (see "Real-Time Simulation and Rendering of 3D Fluids", Crane et. al)
    vec3 gradp = 0.5*vec3(P_xp - P_xn, P_yp - P_yn, P_zp - P_zn)/dL;
    vec3 vnew = v - gradp;
    vnew = vMask*vnew;

    Vair_output = vec4(vnew, 0.0);
}





