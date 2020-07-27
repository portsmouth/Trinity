precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field

/////// output buffers ///////
layout(location = 0) out vec4 divVair_output;

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

    float V_xp = texelFetch(Vair_sampler, mapVsToFrag(X_ip), 0).x;
    float V_xn = texelFetch(Vair_sampler, mapVsToFrag(X_in), 0).x;
    float V_yp = texelFetch(Vair_sampler, mapVsToFrag(X_jp), 0).y;
    float V_yn = texelFetch(Vair_sampler, mapVsToFrag(X_jn), 0).y;
    float V_zp = texelFetch(Vair_sampler, mapVsToFrag(X_kp), 0).z;
    float V_zn = texelFetch(Vair_sampler, mapVsToFrag(X_kn), 0).z;

    // Apply solid no-slip boundary conditions
    if (isSolidCell(X_ip)) V_xp = 0.0;
    if (isSolidCell(X_in)) V_xn = 0.0;
    if (isSolidCell(X_jp)) V_yp = 0.0;
    if (isSolidCell(X_jn)) V_yn = 0.0;
    if (isSolidCell(X_kp)) V_zp = 0.0;
    if (isSolidCell(X_kn)) V_zn = 0.0;

    float divVair = 0.5 * ((V_xp - V_xn) + (V_yp - V_yn) + (V_zp - V_zn)) / dL;
    divVair_output = vec4(divVair);
}





