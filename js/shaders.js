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
uniform float timestep;
uniform float buoyancy;
uniform float gravity;
uniform float radiationLoss;
uniform float T0;             // reference temperature for buoyancy
uniform float Tambient;       // temperature which generates buoyancy which balances gravity

in vec2 v_texcoord;

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field
uniform sampler2D Pair_sampler; // 1, float pressure field
uniform sampler2D Tair_sampler; // 2, float temperature field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Pair_output;
layout(location = 2) out vec4 Tair_output;


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

vec3 back_advect(vec3 wsX, vec3 vX, float h)
{
    // RK4 integration for position wsX advected backwards through time h:
    vec3 k1 = vX;
    vec3 k2 = interp(Vair_sampler, wsX - 0.5*h*k1).xyz;
    vec3 k3 = interp(Vair_sampler, wsX - 0.5*h*k2).xyz;
    vec3 k4 = interp(Vair_sampler, wsX -     h*k3).xyz;
    return clamp(wsX - h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0, vec3(dL), L-vec3(dL));
}

void main()
{
    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsX = mapFragToVs(frag);
    vec3 wsX = vsX*dL;

    // Apply semi-Lagrangian advection
    vec3 v0  = texelFetch(Vair_sampler, frag, 0).xyz;
    vec3 wsXp = back_advect(wsX, v0, timestep); // @todo: clamp to stay in grid bounds
    vec3  v = interp(Vair_sampler, wsXp).rgb;
    float P = interp(Pair_sampler, wsXp).x;
    float T = interp(Tair_sampler, wsXp).x;

    // Apply thermal relaxation to ambient temperature due to "radiation loss"
    float dT = T - Tambient;
    T = Tambient + dT*exp(-radiationLoss);

    // Apply gravity + buoyancy force
    float buoyancy_expansion = buoyancy*(T - T0);
    v.y += timestep * gravity * (1.0 - buoyancy_expansion);

    // Apply solid boundary condition at y=0
    ivec3 vsI = ivec3(vsX);
    if (vsI.y==0) v.y = 0.0;

    Vair_output = vec4(v, 0.0);
    Pair_output = vec4(P);
    Tair_output = vec4(T);
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

'debris-fragment-shader': `#version 300 es
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
uniform float timestep;

in vec2 v_texcoord;

/////// input buffers ///////
uniform sampler2D debris_sampler; // 0, float debris density field
uniform sampler2D Vair_sampler;   // 1, vec3 velocity field

/////// output buffers ///////
layout(location = 0) out vec4 debris_output;

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

vec3 back_advect(vec3 wsX, vec3 vX, float h)
{
    // RK4 integration for position wsX advected backwards through time h:
    vec3 k1 = vX;
    vec3 k2 = interp(Vair_sampler, wsX - 0.5*h*k1).xyz;
    vec3 k3 = interp(Vair_sampler, wsX - 0.5*h*k2).xyz;
    vec3 k4 = interp(Vair_sampler, wsX -     h*k3).xyz;
    return clamp(wsX - h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0, vec3(dL), L-vec3(dL));
}

void main()
{
    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsX = mapFragToVs(frag);
    vec3 wsX = vsX*dL;

    // Apply semi-Lagrangian advection
    vec3 v0 = texelFetch(Vair_sampler, frag, 0).rgb;
    vec3 wsXp = back_advect(wsX, v0, timestep);
    vec3 debris_p = interp(debris_sampler, wsXp).rgb;
    debris_output = vec4(debris_p, 1.0);
}
`,

'debris-vertex-shader': `#version 300 es
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

'div-fragment-shader': `#version 300 es
precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform float dL;

/////// input buffers ///////
uniform sampler2D Vair_sampler; // 0, vec3 velocity field

/////// output buffers ///////
layout(location = 0) out vec4 divVair_output;


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

    vec4 V_xp = texelFetch(Vair_sampler, X_ip, 0);
    vec4 V_xn = texelFetch(Vair_sampler, X_in, 0);
    vec4 V_yp = texelFetch(Vair_sampler, X_jp, 0);
    vec4 V_yn = texelFetch(Vair_sampler, X_jn, 0);
    vec4 V_zp = texelFetch(Vair_sampler, X_kp, 0);
    vec4 V_zn = texelFetch(Vair_sampler, X_kn, 0);

    float divVair = 0.5 * (V_xp.x - V_xn.x + V_yp.y - V_yn.y + V_zp.z - V_zn.z) / dL;
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

'line-fragment-shader': `#version 300 es
precision highp float;

out vec4 outputColor;
uniform vec3 color;

void main() 
{
	outputColor = vec4(color, 1.0);
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

    // Get air pressure, temperature and velocity divergence at voxel
    float P    = texelFetch(Pair_sampler, frag, 0).x;
    float T    = texelFetch(Tair_sampler, frag, 0).x;
    float divv = texelFetch(divVair_sampler, frag, 0).x;

    // Introduce local expansion due to heated fluid
    float phi = timestep * expansion * T;
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

    // @todo: what is the "correct" normalization here?

    float P_relax = avgp - dL*dL*rhs/6.0; // (See Numerical Recipes 19.5.5, extended to 3d)
    //float P_relax = avgp - rhs/6.0; // (See Numerical Recipes 19.5.5, extended to 3d)

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
uniform float contrast;
uniform float saturation;

out vec4 g_outputColor;

float toneMap(float L)
{
  return L / (1.0 + L);
}

void main()
{
    vec3 L = texture(Radiance, vTexCoord).rgb;
    float X = L.x;
    float Y = L.y;
    float Z = L.z;
    
    // convert XYZ tristimulus to linear RGB color space
    vec3 RGB;
    RGB.r =  3.2406*X - 1.5372*Y - 0.4986*Z;
    RGB.g = -0.9689*X + 1.8758*Y + 0.0415*Z;
    RGB.b =  0.0557*X - 0.2040*Y + 1.0570*Z;

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

    // apply contrast
    dR = R - 0.5;
    dG = G - 0.5;
    dB = B - 0.5;
    R = 0.5 + sign(dR)*pow(abs(dR), 1.0/contrast);
    G = 0.5 + sign(dG)*pow(abs(dG), 1.0/contrast);
    B = 0.5 + sign(dB)*pow(abs(dB), 1.0/contrast);

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
uniform float dL;
uniform int Nraymarch;

// Physics
uniform float debrisExtinction;
uniform float blackbodyEmission;
uniform float TtoKelvin;
uniform float exposure;
uniform float invGamma;

/////// input buffers ///////
uniform sampler2D debris_sampler; // 0, float debris density field
uniform sampler2D Tair_sampler;   // 1, float temperature field
uniform sampler2D Vair_sampler;   // 2, vec3 velocity field (for debug, for now)

/////// output buffers ///////
layout(location = 0) out vec4 gbuf_rad;

#define DENOM_EPS 1.0e-7
#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }

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

void planckianLocus(float T_kelvin, inout float xc, inout float yc)
{
    float thOvT = 1000.0/T_kelvin;
    float thOvT2 = thOvT*thOvT;
    float thOvT3 = thOvT2*thOvT;
    if      (T_kelvin<4000.0) xc = -0.2661239*thOvT3 - 0.2343580*thOvT2 + 0.8776956*thOvT + 0.179910;
    else                      xc = -3.0258469*thOvT3 + 2.1070379*thOvT2 + 0.2226347*thOvT + 0.240390;
    float xc2 = xc * xc;
    float xc3 = xc2 * xc;
    if      (T_kelvin<2222.0) yc = -1.1063814*xc3 - 1.34811020*xc2 + 2.18555832*xc - 0.20219683;
    else if (T_kelvin<4000.0) yc = -0.9549476*xc3 - 1.37418593*xc2 + 2.09137015*xc - 0.16748867;
    else                      yc =  3.0817580*xc3 - 5.87338670*xc2 + 3.75112997*xc - 0.37001483;
}

vec3 tempToRGB(float T_kelvin)
{
    if (T_kelvin <= 1000.0) return vec3(T_kelvin/1000.0, 0.0, 0.0);
    float x, y;
    planckianLocus(T_kelvin, x, y);
    float X = x/y;
    float Y = 1.0;
    float Z = (1.f - x - y)/y;
    float R = max(0.0,  3.2410*X - 1.5374*Y - 0.4986*Z);
    float G = max(0.0, -0.9682*X + 1.8760*Y + 0.0416*Z);
    float B = max(0.0,  0.0556*X - 0.2040*Y + 1.0570*Z);
    return vec3(R, G, B);
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

void main()
{
    vec2 frag = gl_FragCoord.xy;
    vec3 rayPos, rayDir;
    constructPrimaryRay(frag, rayPos, rayDir);

    vec3 L = vec3(0.0);
    vec3 Lbackground = vec3(0.2, 0.25, 0.5);
    vec3 Tr = vec3(1.0); // transmittance

    // intersect ray with volume bounds
    float t0, t1;
    if ( boundsIntersect(rayPos, rayDir, volMin, volMax, t0, t1) )
    {
        // Raymarch
        float dl = (t1 - t0) / float(Nraymarch);
        vec3 pMarch = rayPos + (t0+0.5*dl)*rayDir;

        for (int n=0; n<Nraymarch; ++n)
        {
            // transform pMarch into simulation domain:
            vec3 wsP = pMarch - volMin;

            // Absorption by dust
            vec3 debris = interp(debris_sampler, wsP).rgb;
            vec3 sigma = debrisExtinction * debris;
            Tr.r *= exp(-sigma.r*dl);
            Tr.g *= exp(-sigma.g*dl);
            Tr.b *= exp(-sigma.b*dl);

            // Emit blackbody radiation from hot air
            float T = interp(Tair_sampler, wsP).r;

            vec3 blackbody_color = tempToRGB(T * TtoKelvin);
            vec3 emission = pow(blackbodyEmission*T, 4.0) * blackbody_color * Tr;
            L += emission * dl;
            pMarch += rayDir*dl;
        }
    }

    L += Tr * Lbackground;

    // apply gamma correction to convert linear RGB to sRGB
    L = pow(L, vec3(invGamma));
    L *= pow(2.0, exposure);

    vec2 f = frag/resolution.xy;
    gbuf_rad = vec4(L, 1.0);
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

}