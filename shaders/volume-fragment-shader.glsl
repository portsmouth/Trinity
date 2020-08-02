
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
uniform float extinctionScale;
uniform float emissionScale;
uniform float anisotropy;

// Lighting
uniform vec3 skyColor;
uniform vec3 sunColor;
uniform float sunPower;
uniform vec3 sunDir;

// Tonemapping
uniform float exposure;
uniform float invGamma;

/////// input buffers ///////
uniform sampler2D absorption_sampler; // 0, vec3 absorption field
uniform sampler2D scattering_sampler; // 1, vec3 scattering field
uniform sampler2D Tair_sampler;       // 2, float temperature field

/////// output buffers ///////
layout(location = 0) out vec4 gbuf_rad;

/////////////////////// user-defined code ///////////////////////
_USER_CODE_
/////////////////////// user-defined code ///////////////////////


#define M_PI 3.141592653589793
#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }
#define DENOM_EPSILON 1.0e-7

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

vec3 clampToBounds(in vec3 wsX)
{
    vec3 halfVoxel = vec3(0.5*dL);
    return clamp(wsX, halfVoxel, L-halfVoxel);
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

vec3 sun_transmittance(in vec3 pos, float stepSize)
{
    vec3 Tr = vec3(1.0); // transmittance
    float t = shadowHit(pos, sunDir, volMin, volMax);
    int numSteps = int(floor(t/(4.0*stepSize)));
    for (int n=0; n<256; ++n)
    {
        if (n>=numSteps) break;
        float tmin = float(n)*stepSize;
        float tmax = min(tmin + stepSize, t);
        float dt = tmax - tmin;
        float tmid = 0.5*(tmin + tmax);
        vec3 pMarch = pos + tmid*sunDir;
        vec3 wsP = pMarch - volMin; // transform pMarch into simulation domain:
        vec3 absorption = interp(absorption_sampler, clampToBounds(wsP)).rgb;
        vec3 scattering = interp(scattering_sampler, clampToBounds(wsP)).rgb;
        vec3 sigma_t = (absorption + scattering) * extinctionScale;
        Tr *= exp(-sigma_t*dt); // transmittance over step
    }
    return Tr;
}

float phaseFunction(float mu)
{
    float g = anisotropy;
    float gSqr = g*g;
    return (1.0/(4.0*M_PI)) * (1.0 - gSqr) / pow(1.0 - 2.0*g*mu + gSqr, 1.5);
}

void main()
{
    vec2 frag = gl_FragCoord.xy;
    vec3 rayPos, rayDir;
    constructPrimaryRay(frag, rayPos, rayDir);

    float lengthScale = length(volMax - volMin);
    float stepSize = lengthScale / float(Nraymarch);

    vec3 L = vec3(0.0);
    vec3 Lbackground = skyColor;
    vec3 Tr = vec3(1.0); // transmittance along camera ray

    // intersect ray with volume bounds
    float t0, t1;
    if ( boundsIntersect(rayPos, rayDir, volMin, volMax, t0, t1) )
    {
        float t01 = t1 - t0;
        int numSteps = 1 + int(floor(t01/stepSize));
        numSteps = min(256, numSteps);

        for (int n=0; n<1024; ++n)
        {
            if (n>=numSteps)
                break;
            float tmin = t0 + float(n)*stepSize;
            float tmax = min(tmin + stepSize, t1);
            float dt = tmax - tmin;
            float tmid = 0.5*(tmin + tmax);
            vec3 pMarch = rayPos + tmid*rayDir;
            vec3 wsP = pMarch - volMin; // transform pMarch into simulation domain:

            // Compute extinction and albedo at step midpoint
            vec3 absorption = interp(absorption_sampler, clampToBounds(wsP)).rgb;
            vec3 scattering = interp(scattering_sampler, clampToBounds(wsP)).rgb;
            vec3 sigma_t = (absorption + scattering) * extinctionScale;
            vec3 albedo = scattering / max(absorption, vec3(DENOM_EPSILON));

            // Compute in-scattered sunlight
            vec3 Li = sunPower * sunColor * sun_transmittance(pMarch, stepSize);
            vec3 dTr = exp(-sigma_t*dt); // transmittance over step
            vec3 J = Li * albedo * Tr * (vec3(1.0) - dTr); // scattering term integrated over step, assuming constant Li
            L += J * phaseFunction(dot(sunDir, rayDir));

            // Emit blackbody radiation from hot air
            float T = interp(Tair_sampler, clampToBounds(wsP)).r;
            vec3 emission = emissionScale * temperatureToEmission(T);
            L += Tr * emission * dt;

            // Update transmittance for start of next step
            Tr *= dTr;
        }
    }

    L += Tr * Lbackground;

    // apply gamma correction to convert linear RGB to sRGB
    L = pow(L, vec3(invGamma));
    L *= pow(2.0, exposure);

    vec2 f = frag/resolution.xy;
    gbuf_rad = vec4(L, 1.0);
}

