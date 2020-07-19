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
layout(location = 0) out vec4 vorticity_output;

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

vec3 curl_v(in ivec3 vsXi)
{
    int ix = vsXi.x;
    int iy = vsXi.y;
    int iz = vsXi.z;
    ivec3 X_ip = ivec3(min(ix+1, Nx-1), iy, iz);
    ivec3 X_in = ivec3(max(ix-1,    0), iy, iz);
    ivec3 X_jp = ivec3(ix, min(iy+1, Ny-1), iz);
    ivec3 X_jn = ivec3(ix, max(iy-1,    0), iz);
    ivec3 X_kp = ivec3(ix, iy, min(iz+1, Nz-1));
    ivec3 X_kn = ivec3(ix, iy, max(iz-1,    0));

    // Get velocity values on stencil
    vec3 v_xp = texelFetch(Vair_sampler, mapVsToFrag(X_ip), 0).rgb;
    vec3 v_xn = texelFetch(Vair_sampler, mapVsToFrag(X_in), 0).rgb;
    vec3 v_yp = texelFetch(Vair_sampler, mapVsToFrag(X_jp), 0).rgb;
    vec3 v_yn = texelFetch(Vair_sampler, mapVsToFrag(X_jn), 0).rgb;
    vec3 v_zp = texelFetch(Vair_sampler, mapVsToFrag(X_kp), 0).rgb;
    vec3 v_zn = texelFetch(Vair_sampler, mapVsToFrag(X_kn), 0).rgb;

    // Construct components of curl(v)
    float inv2dL = 0.5/dL;
    float dvy_dx = (v_xp.y - v_xn.y) * inv2dL;
    float dvz_dx = (v_xp.z - v_xn.z) * inv2dL;
    float dvx_dy = (v_yp.x - v_yn.x) * inv2dL;
    float dvz_dy = (v_yp.z - v_yn.z) * inv2dL;
    float dvx_dz = (v_zp.x - v_zn.x) * inv2dL;
    float dvy_dz = (v_zp.y - v_zn.y) * inv2dL;
    vec3 curl = vec3(dvz_dy - dvy_dz, dvx_dz - dvz_dx, dvy_dx - dvx_dy);
    return curl;
}

void main()
{
    // Compute vorticity
    ivec2 frag = ivec2(gl_FragCoord.xy);
    ivec3 vsX = ivec3(floor(mapFragToVs(frag)));
    vec3 vorticity = curl_v(vsX);
    vorticity_output = vec4(vorticity, 0.0);
}





