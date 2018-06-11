
precision highp float;

in vec2 vTexCoord;

layout(location = 0) out vec4 gbuf_rad;

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
uniform float volRadius;
uniform float volHeight;
uniform int Nr;
uniform int Ny;
uniform int Nraymarch;

uniform float extinctionMultiplier;
uniform float emissionMultiplier;
uniform float tempMultiplier;

uniform sampler2D Qair;    // air simulation
uniform sampler2D Qdebris; // debris simulation

#define DENOM_EPS 1.0e-7
#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }

bool boundsIntersect( in vec3 rayPos, in vec3 rayDir, in vec3 bbMin, in vec3 bbMax,
                      inout float t0, inout float t1 )
{
    vec3 dL = vec3(1.0f/rayDir.x, 1.0f/rayDir.y, 1.0f/rayDir.z);
    vec3 lo = (bbMin - rayPos) * dL;
    vec3 hi = (bbMax - rayPos) * dL;
    sort2(lo, hi);
    bool hit = !( lo.x>hi.y || lo.y>hi.x || lo.x>hi.z || lo.z>hi.x || lo.y>hi.z || lo.z>hi.y );
    t0 = max(max(lo.x, lo.y), lo.z);
    t1 = min(min(hi.x, hi.y), hi.z);
    return hit;
}

void constructPrimaryRay(in vec2 pixel,
                         inout vec3 rayPos, inout vec3 rayDir)
{
    // Compute world ray direction for given (possibly jittered) fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
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
    if (T_kelvin <= 1000.0) return vec3(0.0);
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

void main()
{
    vec2 pixel = gl_FragCoord.xy;
    vec3 rayPos, rayDir;
    constructPrimaryRay(pixel, rayPos, rayDir);

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
            float y = pMarch.y - volMin.y;
            float r = length((pMarch - volCenter).xz);
            if (r<=volRadius)
            {
                int ir = clamp(int(floor(r)), 0, Nr-1);
                int iy = clamp(int(floor(y)), 0, Ny-1);
                float u = r/volRadius;
                float v = y/volHeight;

                // Absorption by dust
                vec4 Qdebris_  = texture(Qdebris, vec2(u, v));
                vec3 sigma = extinctionMultiplier * Qdebris_.r * vec3(Qdebris_.g, Qdebris_.b, Qdebris_.a);
                Tr.r *= exp(-sigma.r*dl);
                Tr.g *= exp(-sigma.g*dl);
                Tr.b *= exp(-sigma.b*dl);

                // Emit blackbody radiation from hot air
                vec4 Qair_ = texture(Qair, vec2(u, v));
                float T = tempMultiplier * Qair_.b;
                vec3 blackbody_color = tempToRGB(T * 300.0 * tempMultiplier);
                vec3 emission = Tr * emissionMultiplier * 0.001 * T * blackbody_color;
            
                L += emission * dl;
            }
            pMarch += rayDir*dl;
        }  
    }

    L += Tr * Lbackground;
    
    gbuf_rad = vec4(L, 1.0);
}





