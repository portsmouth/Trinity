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
uniform float vorticity_scale;

in vec2 v_texcoord;

/////// input buffers ///////
uniform sampler2D Vair_sampler;       // 0, vec3 air velocity field
uniform sampler2D Pair_sampler;       // 1, float air pressure field
uniform sampler2D Tair_sampler;       // 2, float air temperature field
uniform sampler2D absorption_sampler; // 3, vec3 absorption field
uniform sampler2D scattering_sampler; // 4, vec3 scattering field
uniform sampler2D vorticity_sampler;  // 5, vec3 air vorticity field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Pair_output;
layout(location = 2) out vec4 Tair_output;
layout(location = 3) out vec4 absorption_output;
layout(location = 4) out vec4 scattering_output;

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

vec2 slicetoUV(int j, vec3 vsP)
{
    // Given y-slice index j, and continuous voxel space xz-location,
    // return corresponding continuous frag UV for interpolation within this slice
    int row = int(floor(float(j)/float(Ncol)));
    int col = j - row*Ncol;
    vec2 uv_ll = vec2(float(col*Nx)/float(W), float(row*Nz)/float(H));
    float du = vsP.x/float(W);
    float dv = vsP.z/float(H);
    return uv_ll + vec2(du, dv);
}

vec4 interp(in sampler2D S, in vec3 wsP)
{
    vec3 vsP = wsP / dL;
    float pY = vsP.y - 0.5;
    int jlo = clamp(int(floor(pY)), 0, Ny-1); // lower j-slice
    int jhi = clamp(         jlo+1, 0, Ny-1); // upper j-slice
    float flo = float(jhi) - pY;              // lower j fraction
    float fhi = 1.0 - flo;                    // upper j fraction
    vec2 uv_lo = slicetoUV(jlo, vsP);
    vec2 uv_hi = slicetoUV(jhi, vsP);
    vec4 Slo = texture(S, uv_lo);
    vec4 Shi = texture(S, uv_hi);
    return flo*Slo + fhi*Shi;
}

vec3 clampToBounds(in vec3 wsP)
{
    vec3 voxel = vec3(dL);
    return clamp(wsP, voxel, L-voxel);
}

vec3 back_advect(in vec3 wsP, in vec3 vX, float h)
{
    // RK4 integration for position wsP advected backwards through time h:
    vec3 k1 = vX;
    vec3 k2 = interp(Vair_sampler, clampToBounds(wsP - 0.5*h*k1)).xyz;
    vec3 k3 = interp(Vair_sampler, clampToBounds(wsP - 0.5*h*k2)).xyz;
    vec3 k4 = interp(Vair_sampler, clampToBounds(wsP -     h*k3)).xyz;
    return clampToBounds(wsP - h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0);
}

vec3 vorticityConfinementForce(ivec2 frag, in ivec3 vsPi)
{
    int ix = vsPi.x;
    int iy = vsPi.y;
    int iz = vsPi.z;
    ivec3 X_ip = ivec3(min(ix+1, Nx-1), iy, iz);
    ivec3 X_in = ivec3(max(ix-1,    0), iy, iz);
    ivec3 X_jp = ivec3(ix, min(iy+1, Ny-1), iz);
    ivec3 X_jn = ivec3(ix, max(iy-1,    0), iz);
    ivec3 X_kp = ivec3(ix, iy, min(iz+1, Nz-1));
    ivec3 X_kn = ivec3(ix, iy, max(iz-1,    0));

    // Take gradient of vorticity magnitude field
    float omega_xp = length(texelFetch(vorticity_sampler, mapVsToFrag(X_ip), 0).rgb);
    float omega_xn = length(texelFetch(vorticity_sampler, mapVsToFrag(X_in), 0).rgb);
    float omega_yp = length(texelFetch(vorticity_sampler, mapVsToFrag(X_jp), 0).rgb);
    float omega_yn = length(texelFetch(vorticity_sampler, mapVsToFrag(X_jn), 0).rgb);
    float omega_zp = length(texelFetch(vorticity_sampler, mapVsToFrag(X_kp), 0).rgb);
    float omega_zn = length(texelFetch(vorticity_sampler, mapVsToFrag(X_kn), 0).rgb);

    vec3 eta = vec3(omega_xp - omega_xn, omega_yp - omega_yn, omega_zp - omega_zn);
    vec3 N = eta / (length(eta) + 1.0e-6);
    vec3 omega = texelFetch(vorticity_sampler, frag, 0).rgb;
    return vorticity_scale * cross(N, omega);
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

    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsP = mapFragToVs(frag);
    ivec3 vsPi = ivec3(floor(vsP));

    if (isSolidCell(vsPi))
    {
        Vair_output       = texelFetch(Vair_sampler, frag, 0);
        Pair_output       = texelFetch(Pair_sampler, frag, 0);
        Tair_output       = texelFetch(Tair_sampler, frag, 0);
        absorption_output = texelFetch(absorption_sampler, frag, 0);
        scattering_output = texelFetch(scattering_sampler, frag, 0);
    }
    else
    {
        // Apply semi-Lagrangian advection
        vec3 v0 = texelFetch(Vair_sampler, frag, 0).xyz;
        vec3 wsP = vsP*dL;
        vec3 wsPp = back_advect(wsP, v0, timestep);
        vec3  v         = interp(Vair_sampler, wsPp).rgb;
        float P         = interp(Pair_sampler, wsPp).x;
        float T         = interp(Tair_sampler, wsPp).x;
        vec3 absorption = interp(absorption_sampler, wsPp).rgb;
        vec3 scattering = interp(scattering_sampler, wsPp).rgb;

        // Apply external forces
        v += timestep * externalForces(wsP, time, v, P, T, L);

        // Apply vorticity confinement:
        if (vorticity_scale > 0.0)
            v += timestep * vorticityConfinementForce(frag, vsPi);

        Vair_output       = vec4(v, 0.0);
        Pair_output       = vec4(vec3(P), 0.0);
        Tair_output       = vec4(vec3(T), 0.0);
        absorption_output = vec4(absorption, 0.0);
        scattering_output = vec4(scattering, 0.0);
    }
}





