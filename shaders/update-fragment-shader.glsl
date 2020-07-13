precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
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

void main()
{
    // Setup local stencil:
    ivec2 frag = ivec2(gl_FragCoord.xy);
    ivec3 vsX = ivec3(mapFragToVs(frag));
    int ix = vsX.x;
    int iy = vsX.y;
    int iz = vsX.z;

    // Apply Neumann boundary conditions
    ivec2 X_ip = mapVsToFrag(ivec3(min(ix+1, Nx-1), iy, iz));
    ivec2 X_in = mapVsToFrag(ivec3(max(ix-1,    0), iy, iz));
    ivec2 X_jp = mapVsToFrag(ivec3(ix, min(iy+1, Ny-1), iz));
    ivec2 X_jn = mapVsToFrag(ivec3(ix, max(iy-1,    0), iz));
    ivec2 X_kp = mapVsToFrag(ivec3(ix, iy, min(iz+1, Nz-1)));
    ivec2 X_kn = mapVsToFrag(ivec3(ix, iy, max(iz-1,    0)));

    // air velocity at voxel
    vec3 v = texelFetch(Vair_sampler, frag, 0).rgb;

    // Compute local gradient of pressure field
    float P_xp = texelFetch(Pair_sampler, X_ip, 0).r;
    float P_xn = texelFetch(Pair_sampler, X_in, 0).r;
    float P_yp = texelFetch(Pair_sampler, X_jp, 0).r;
    float P_yn = texelFetch(Pair_sampler, X_jn, 0).r;
    float P_zp = texelFetch(Pair_sampler, X_kp, 0).r;
    float P_zn = texelFetch(Pair_sampler, X_kn, 0).r;
    float dpdx = 0.5*(P_xp - P_xn)/dL;
    float dpdy = 0.5*(P_yp - P_yn)/dL;
    float dpdz = 0.5*(P_zp - P_zn)/dL;
    vec3 gradp = vec3(dpdx, dpdy, dpdz);

    // Update air velocity accordingly
    Vair_output = vec4(v - timestep*gradp, 0.0);
}





