var Shaders = {

'advect-fragment-shader': `#version 300 es
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
    return collisionSDF(wsP, time, L, dL) < 0.0;
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
        vec3 wsP_back = back_advect(wsP, v0, timestep);
        vec3 v          = interp(Vair_sampler, wsP_back).rgb;
        float P         = interp(Pair_sampler, wsP_back).r;
        vec4 T          = interp(Tair_sampler, wsP_back);
        vec3 absorption = interp(absorption_sampler, wsP_back).rgb;
        vec3 scattering = interp(scattering_sampler, wsP_back).rgb;
        vec3 medium = absorption + scattering;

        // Apply external forces
        v += timestep * externalForces(wsP, time, L, dL, v, P, T, medium);

        // Apply vorticity confinement:
        if (vorticity_scale > 0.0)
            v += timestep * vorticityConfinementForce(frag, vsPi);

        Vair_output       = vec4(v, 0.0);
        Pair_output       = vec4(vec3(P), 0.0);
        Tair_output       = T;
        absorption_output = vec4(absorption, 0.0);
        scattering_output = vec4(scattering, 0.0);
    }
}
`,

'advect-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 v_texcoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    v_texcoord = TexCoord;
}
`,

'copy-fragment-shader': `#version 300 es
precision highp float;

/////// input buffers ///////
uniform sampler2D Qin;

/////// output buffers ///////
layout(location = 0) out vec4 Qcopy;


void main()
{
    ivec2 frag = ivec2(gl_FragCoord.xy);
    Qcopy = texelFetch(Qin, frag, 0);
}
`,

'copy-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'div-fragment-shader': `#version 300 es
precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;
uniform float time;

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
    return collisionSDF(wsP, time, L, dL) < 0.0;
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
`,

'div-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'initial-fragment-shader': `#version 300 es
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

void main()
{
    init();

    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsP = mapFragToVs(frag);
    ivec3 vsPi = ivec3(floor(vsP));
    vec3 wsP = vsP*dL;

    vec3 v = vec3(0.0);
    vec4 T = vec4(0.0);
    vec3 mediumDensity = vec3(0.0);
    vec3 mediumAlbedo = vec3(0.0);
    initial_conditions(wsP, L, dL,
                       v, T,
                       mediumDensity, mediumAlbedo);

    vec3 scattering = mediumDensity * mediumAlbedo;
    vec3 absorption = mediumDensity - scattering;

    Vair_output   = vec4(v, 0.0);
    Pair_output   = vec4(0.0);
    Tair_output   = T;
    absorption_output = vec4(absorption, 0.0);
    scattering_output = vec4(scattering, 0.0);
}
`,

'initial-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'inject-fragment-shader': `#version 300 es
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

/////// input buffers ///////
uniform sampler2D Vair_sampler;       // 0, vec3 air velocity field
uniform sampler2D Tair_sampler;       // 1, float air temperature field
uniform sampler2D absorption_sampler; // 2, vec3 absorption field
uniform sampler2D scattering_sampler; // 3, vec3 scattering field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Tair_output;
layout(location = 2) out vec4 absorption_output;
layout(location = 3) out vec4 scattering_output;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////

bool isSolidCell(in ivec3 vsPi)
{
    vec3 vsP = vec3(float(vsPi.x)+0.5, float(vsPi.y)+0.5,float(vsPi.z)+0.5);
    vec3 wsP = vsP*dL;
    return collisionSDF(wsP, time, L, dL) < 0.0;
}

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

void main()
{
    init();

    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsP = mapFragToVs(frag);
    ivec3 vsPi = ivec3(floor(vsP));
    vec3 wsP = vsP*dL;

    if (isSolidCell(vsPi))
    {
        Vair_output   = texelFetch(Vair_sampler, frag, 0);
        Tair_output   = texelFetch(Tair_sampler, frag, 0);
        absorption_output = texelFetch(absorption_sampler, frag, 0);
        scattering_output = texelFetch(scattering_sampler, frag, 0);
    }
    else
    {
        // Get current velocity, temperature, and debris fields:
        vec3 v = texelFetch(Vair_sampler, frag, 0).rgb;
        vec4 T = texelFetch(Tair_sampler, frag, 0);
        vec3 absorption = texelFetch(absorption_sampler, frag, 0).rgb;
        vec3 scattering = texelFetch(scattering_sampler, frag, 0).rgb;

        // Inject mass and modify temperature:
        vec3 vInflow      = vec3(0.0);
        vec4 Tinflow      = vec4(0.0);
        vec3 mediumInflow = vec3(0.0);
        vec3 mediumAlbedo = vec3(0.5);

        inject(wsP, time, L, dL,
               v, vInflow,
               T, Tinflow,
               mediumInflow, mediumAlbedo);

        v += vInflow * timestep;
        T += Tinflow * timestep;

        vec3 scatteringInflow = mediumAlbedo * mediumInflow;
        vec3 absorptionInflow = mediumInflow - scatteringInflow;
        absorption += absorptionInflow * timestep;
        scattering += scatteringInflow * timestep;

        Vair_output       = vec4(v, 0.0);
        Tair_output       = T;
        absorption_output = vec4(absorption, 0.0);
        scattering_output = vec4(scattering, 0.0);
    }
}
`,

'inject-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'line-fragment-shader': `#version 300 es
precision highp float;

out vec4 outputColor;
uniform vec4 color;

void main() 
{
	outputColor = color;
}
`,

'line-vertex-shader': `#version 300 es
precision highp float;

uniform mat4 u_projectionMatrix;
uniform mat4 u_modelViewMatrix;

in vec3 Position;

void main()
{
	gl_Position = u_projectionMatrix * u_modelViewMatrix * vec4(Position, 1.0);
}
`,

'project-fragment-shader': `#version 300 es
precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;

// Physics
uniform float time;
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
    return collisionSDF(wsP, time, L, dL) < 0.0;
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
`,

'project-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'tonemap-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Radiance;
in vec2 vTexCoord;

uniform float exposure;
uniform float invGamma;
uniform float saturation;

out vec4 g_outputColor;

float toneMap(float L)
{
  return L / (1.0 + L);
}

void main()
{
    vec3 RGB = texture(Radiance, vTexCoord).rgb;

    // apply gamma correction to convert linear RGB to sRGB
    RGB = pow(RGB, vec3(invGamma));

    // deal with out-of-gamut RGB.
    float delta = -min(0.0, min(min(RGB.r, RGB.g), RGB.b));
    RGB.r += delta;
    RGB.g += delta;
    RGB.b += delta;

    // apply tonemapping
    RGB *= pow(2.0, exposure);
    float R = RGB.r;
    float G = RGB.g;
    float B = RGB.b;
    R = toneMap(R);
    G = toneMap(G);
    B = toneMap(B);

    // apply saturation
    float mean = (R + G + B)/3.0;
    float dR = R - mean;
    float dG = G - mean;
    float dB = B - mean;
    R = mean + sign(dR)*pow(abs(dR), 1.0/saturation);
    G = mean + sign(dG)*pow(abs(dG), 1.0/saturation);
    B = mean + sign(dB)*pow(abs(dB), 1.0/saturation);

    g_outputColor = vec4(vec3(R,G,B), 1.0);
}
`,

'tonemap-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;
out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'update-fragment-shader': `#version 300 es
precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;

// Physics
uniform float time;
uniform float timestep;

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field
uniform sampler2D Pair_sampler; // 1, float pressure field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;

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
    return collisionSDF(wsP, time, L, dL) < 0.0;
}

void main()
{
    init();

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
`,

'update-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'volume-fragment-shader': `#version 300 es
precision highp float;

in vec2 vTexCoord;

// Camera
uniform vec2 resolution;
uniform vec3 camPos;
uniform vec3 camDir;
uniform vec3 camX;
uniform vec3 camY;
uniform float camFovy; // degrees
uniform float camAspect;
uniform vec3 volMin;
uniform vec3 volMax;
uniform vec3 volCenter;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform int W;
uniform int H;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;
uniform int Nraymarch;

// Physics
uniform float time;
uniform float extinctionScale;
uniform float emissionScale;
uniform float anisotropy;

// Lighting
uniform vec3 skyColor;
uniform vec3 sunColor;
uniform float sunPower;
uniform vec3 sunDir;
uniform vec3 colliderDiffuse;
uniform vec3 colliderSpec;
uniform float colliderRoughness;

// Progressive rendering
uniform int spp;

// Tonemapping
uniform float exposure;
uniform float invGamma;

/////// input buffers ///////
uniform sampler2D Radiance;           // 0 (last frame radiance)
uniform sampler2D RngData;            // 1 (last frame rng seed)
uniform sampler2D absorption_sampler; // 2, vec3 absorption field
uniform sampler2D scattering_sampler; // 3, vec3 scattering field
uniform sampler2D Tair_sampler;       // 4, float temperature field

/////// output buffers ///////
layout(location = 0) out vec4 gbuf_rad;
layout(location = 1) out vec4 gbuf_rng;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////

#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }
#define DENOM_EPSILON 1.0e-7

/// GLSL floating point pseudorandom number generator, from
/// "Implementing a Photorealistic Rendering System using GLSL", Toshiya Hachisuka
/// http://arxiv.org/pdf/1505.06022.pdf
float rand(inout vec4 rnd)
{
    const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
    const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
    const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
    const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);
    vec4 beta = floor(rnd/q);
    vec4 p = a*(rnd - beta*q) - beta*r;
    beta = (1.0 - sign(p))*0.5*m;
    rnd = p + beta;
    return fract(dot(rnd/m, vec4(1.0, -1.0, 1.0, -1.0)));
}

bool boundsIntersect( in vec3 rayPos, in vec3 rayDir, in vec3 bbMin, in vec3 bbMax,
                      inout float t0, inout float t1 )
{
    vec3 dX = vec3(1.0f/rayDir.x, 1.0f/rayDir.y, 1.0f/rayDir.z);
    vec3 lo = (bbMin - rayPos) * dX;
    vec3 hi = (bbMax - rayPos) * dX;
    sort2(lo, hi);
    bool hit = !( lo.x>hi.y || lo.y>hi.x || lo.x>hi.z || lo.z>hi.x || lo.y>hi.z || lo.z>hi.y );
    t0 = max(max(lo.x, lo.y), lo.z);
    t1 = min(min(hi.x, hi.y), hi.z);
    return hit;
}

float shadowHit(in vec3 rayPos, in vec3 rayDir, in vec3 bbMin, in vec3 bbMax)
{
    // (get first hit assuming rayPos is interior to the AABB)
    vec3 dX = vec3(1.0f/rayDir.x, 1.0f/rayDir.y, 1.0f/rayDir.z);
    vec3 lo = (bbMin - rayPos) * dX;
    vec3 hi = (bbMax - rayPos) * dX;
    sort2(lo, hi);
    return min(min(hi.x, hi.y), hi.z);
}

void constructPrimaryRay(in vec2 frag,
                         inout vec3 rayPos, inout vec3 rayDir)
{
    // Compute world ray direction for given (possibly jittered) fragment
    vec2 ndc = -1.0 + 2.0*(frag/resolution.xy);
    float fh = tan(0.5*radians(camFovy)); // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    rayDir = normalize(camDir + s);
    rayPos = camPos;
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

ivec2 mapVsToFrag(in ivec3 vsP)
{
    // map integer voxel space coords to the corresponding fragment coords
    int i = vsP.x;
    int j = vsP.y;
    int k = vsP.z;
    int ui = Nx*j + i;
    int vi = k;
    return ivec2(ui, vi);
}

vec3 sun_transmittance(in vec3 pW, float stepSize, vec4 rnd)
{
    vec3 Tr = vec3(1.0); // transmittance
    float t = shadowHit(pW, sunDir, volMin, volMax);
    float xi = rand(rnd);
    int numSteps = int(ceil(t/stepSize+xi));
    for (int n=0; n<256; ++n)
    {
        if (n>=numSteps) break;
        float dt;
        float tmid;
        if (numSteps==1)
        {
            dt = t;
            tmid = dt*xi;
        }
        else
        {
            float tmin = (float(n) - xi)*stepSize;
            float tmax = min(tmin + stepSize, t);
            tmin = max(0.0, tmin);
            dt = tmax - tmin;
            tmid = 0.5*(tmin+tmax);
        }
        vec3 pMarch = pW + tmid*sunDir;
        vec3 wsP = pMarch - volMin; // transform pMarch into simulation domain:
        vec3 absorption = interp(absorption_sampler, clampToBounds(wsP)).rgb;
        vec3 scattering = interp(scattering_sampler, clampToBounds(wsP)).rgb;
        vec3 sigma_t = (absorption + scattering) * extinctionScale;
        Tr *= exp(-sigma_t*dt); // transmittance over step
    }
    return Tr;
}

float _sdBox(vec3 pW, vec3 bmin, vec3 bmax)
{
    vec3 d = abs(pW-0.5*(bmin+bmax)) - 0.5*(bmax-bmin);
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float _collisionSDF(vec3 wsP)
{
    float s = _sdBox(wsP, volMin, volMax);
    return max(s, collisionSDF(wsP, time, L, dL));
}

bool traceSDF(in vec3 start, in vec3 dir, float tend, float lengthScale, inout float t)
{
    float minMarch = 1.0e-3 * lengthScale;
    const float HUGE_VAL = 1.0e20;
    vec3 pW = start;
    vec3 wsP = pW - volMin;
    float sdf = _collisionSDF(wsP);
    float InitialSign = sign(sdf);
    t = InitialSign * sdf;
    if (t>=tend) return false;
    for (int n=0; n<256; n++)
    {
        vec3 pW = start + t*dir;
        wsP = pW - volMin;
        sdf = abs(_collisionSDF(wsP));
        if (sdf<minMarch)
            return true;
        t += InitialSign * sdf;
        if (t>=tend) return false;
    }
    return false;
}

vec3 normalSDF(in vec3 pW, float lengthScale)
{
    // Compute normal as gradient of SDF
    float normalEpsilon = 2.0e-3 * lengthScale;
    vec3 e = vec3(normalEpsilon, 0.0, 0.0);
    vec3 Xp = pW+e.xyy; vec3 Xn = pW-e.xyy;
    vec3 Yp = pW+e.yxy; vec3 Yn = pW-e.yxy;
    vec3 Zp = pW+e.yyx; vec3 Zn = pW-e.yyx;
    vec3 N;
    N = vec3(   _collisionSDF(Xp) - _collisionSDF(Xn),
                _collisionSDF(Yp) - _collisionSDF(Yn),
                _collisionSDF(Zp) - _collisionSDF(Zn));
    return normalize(N);
}


vec3 colliderRadiance(in vec3 pW, in vec3 rayDir, float lengthScale, float stepSize, inout vec4 rnd)
{
    vec3 N = normalSDF(pW, lengthScale); // local SDF normal
    vec3 L = sunDir;
    float LN = max(0.0, dot(L, N));
    vec3 Li = sunPower * sunColor * sun_transmittance(pW, 3.0*stepSize, rnd); // sun radiance
    vec3 R = -reflect(L, N);
    vec3 V = -rayDir;
    vec3 phong = colliderDiffuse*LN + colliderSpec*pow(max(0.0, dot(R, V)), 1.0/colliderRoughness);
    vec3 ambient = 0.05*colliderDiffuse;
    return (phong*Li + ambient);
}

void main()
{
    vec2 frag = gl_FragCoord.xy;
    vec3 rayPos, rayDir;
    constructPrimaryRay(frag, rayPos, rayDir);

    vec3 L = vec3(0.0);
    vec3 Lbackground = skyColor;
    vec3 Tr = vec3(1.0); // transmittance along camera ray

    // jitter raymarch randomly:
    vec4 rnd = texture(RngData, vTexCoord);
    float xi = rand(rnd);
    float lengthScale = length(volMax - volMin);
    float stepSize = lengthScale / float(Nraymarch);

    // intersect ray with volume bounds
    float t0, t1;
    if ( boundsIntersect(rayPos, rayDir, volMin, volMax, t0, t1) )
    {
        t0 = max(0.0, t0);
        float t01 = t1 - t0;

        // Check for hit with collision object
        vec3 trace_start = rayPos + t0*rayDir;
        float t_trace;
        if (traceSDF(trace_start, rayDir, t01, lengthScale, t_trace))
        {
            vec3 hit = trace_start + t_trace*rayDir;
            Lbackground = colliderRadiance(hit, rayDir, lengthScale, stepSize, rnd);
            t1 = t0 + t_trace;
            t01 = t1 - t0;
        }

        int numSteps = int(ceil(t01/stepSize+xi));
        numSteps = min(256, numSteps);
        for (int n=0; n<1024; ++n)
        {
            if (n>=numSteps) break;

            // Compute step bounds
            float dt;
            float tmid;
            if (numSteps==1)
            {
                dt = t01;
                tmid = t0 + dt*xi;
            }
            else
            {
                float tmin = t0 + (float(n) - xi)*stepSize;
                float tmax = min(tmin + stepSize, t1);
                tmin = max(t0, tmin);
                dt = tmax - tmin;
                tmid = 0.5*(tmin+tmax);
            }

            // transform pMarch into simulation domain:
            vec3 pMarch = rayPos + tmid*rayDir;
            vec3 wsP = pMarch - volMin;

            // Compute extinction and albedo at step midpoint
            vec3 absorption = interp(absorption_sampler, clampToBounds(wsP)).rgb;
            vec3 scattering = interp(scattering_sampler, clampToBounds(wsP)).rgb;
            vec3 mediumExtinction = (absorption + scattering);
            vec3 mediumAlbedo = scattering / max(mediumExtinction, vec3(DENOM_EPSILON));
            mediumRemap(mediumExtinction, mediumAlbedo);
            mediumExtinction *= extinctionScale;

            // Compute in-scattered sunlight
            vec3 Li = sunPower * sunColor * sun_transmittance(pMarch, 3.0*stepSize, rnd);
            vec3 dTr = exp(-mediumExtinction*dt); // transmittance over step
            vec3 J = Li * mediumAlbedo * Tr * (vec3(1.0) - dTr); // scattering term integrated over step, assuming constant Li
            L += J * phaseFunction(dot(sunDir, rayDir), anisotropy);

            // Emit radiation from hot air
            vec4 T = interp(Tair_sampler, clampToBounds(wsP));
            vec3 emission = emissionScale * temperatureToEmission(T);
            L += Tr * emission * dt;

            // Update transmittance for start of next step
            Tr *= dTr;
        }
    }

    // Add background "sky" radiance, attenuated by volume
    L += Tr * Lbackground;

    // Add latest radiance estimate to the the Monte Carlo average
    vec4 oldL = vec4(0.0);
    if (spp>0)
        oldL = texture(Radiance, vTexCoord);
    vec3 newL = (float(spp)*oldL.rgb + L) / (float(spp) + 1.0);

    // Output radiance
    gbuf_rad = vec4(newL, 1.0);
    gbuf_rng = rnd;
}
`,

'volume-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

out vec2 vTexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
    vTexCoord = TexCoord;
}
`,

'vorticity-fragment-shader': `#version 300 es
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
`,

'vorticity-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

}