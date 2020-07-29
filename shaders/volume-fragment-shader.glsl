
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
uniform float blackbodyEmission;
uniform float TtoKelvin;
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
uniform sampler2D debris_sampler; // 0, float debris density field
uniform sampler2D Tair_sampler;   // 1, float temperature field
uniform sampler2D Vair_sampler;   // 2, vec3 velocity field (for debug, for now)

/////// output buffers ///////
layout(location = 0) out vec4 gbuf_rad;

#define M_PI 3.141592653589793
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

float shadowHit(in vec3 rayPos, in vec3 rayDir, in vec3 bbMin, in vec3 bbMax)
{
    // (get first hit assuming rayPos is interior to the AABB)
    vec3 dX = vec3(1.0f/rayDir.x, 1.0f/rayDir.y, 1.0f/rayDir.z);
    vec3 lo = (bbMin - rayPos) * dX;
    vec3 hi = (bbMax - rayPos) * dX;
    sort2(lo, hi);
    return min(min(hi.x, hi.y), hi.z);
}

/*
bool inBounds(in vec3 pos, in vec3 bbMin, in vec3 bbMax)
{
    vec3 s = step(bbMin, pos) - step(bbMax, pos);
    return bool(s.x * s.y * s.z);
}
*/

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
        vec3 debrisExtinction = interp(debris_sampler, clampToBounds(wsP)).rgb;
        vec3 sigma_t = debrisExtinction * extinctionScale;
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

            // Dust extinction and albedo at step midpoints
            vec3 debrisExtinction = interp(debris_sampler, clampToBounds(wsP)).rgb;
            vec3 debrisAlbedo = vec3(0.5, 0.5, 0.5); // @todo: take from simulation
            vec3 sigma_t = debrisExtinction * extinctionScale;

            // Incident sunlight radiance at scattering point
            vec3 Li = sunPower * sunColor * sun_transmittance(pMarch, stepSize);

            vec3 dTr = exp(-sigma_t*dt); // transmittance over step
            vec3 J = Li * debrisAlbedo * Tr * (vec3(1.0) - dTr); // scattering term integrated over step, assuming constant Li
            L += J * phaseFunction(dot(sunDir, rayDir));

            // Emit blackbody radiation from hot air
            float T = interp(Tair_sampler, clampToBounds(wsP)).r;
            vec3 blackbody_color = tempToRGB(T * TtoKelvin);
            vec3 emission = pow(blackbodyEmission*T, 4.0) * blackbody_color;
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





