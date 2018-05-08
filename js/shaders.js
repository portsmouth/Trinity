var Shaders = {

'advect-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Qair;

uniform int Nr;
uniform int Ny;
uniform float Delta;
uniform float g;        // gravitational acceleration
uniform float beta;     // buoyancy

out vec2 v_texcoord;
out vec4 Qair_next;

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    vec4 Q = texelFetch(Qair, X,    0);
    float vr = Q.r;
    float vy = Q.g;
    float  T = Q.b;
    float  p = Q.a;

    // Semi-Lagrangian advection 
    float u_advect = clamp(v_texcoord.x - vr/float(Nr), 0.0, 1.0);
    float v_advect = clamp(v_texcoord.y - vy/float(Ny), 0.0, 1.0);
    vec4 Q_advect = texture(Qair, vec2(u_advect, v_advect),    0);

    // vr advect
    Qair_next.r = Q_advect.r;

    // vz advect + force
    Qair_next.g = Q_advect.g + g + beta * (Q_advect.b); // - T0

    // T advect
    Qair_next.b = Q_advect.b;

    // p copy
    Qair_next.a = Q.a;
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

uniform sampler2D Qred;
uniform sampler2D Qblack;

out vec4 Qcopy;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
    int black = int(mod(p.x + mod(p.y, 2), 2)); 
    if (black)
    {
        Qcopy = texelFetch(Qblack, X,    0);
    }
    else
    {
        Qcopy = texelFetch(Qred, X,    0);
    }
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

'project-black-fragment.shader': `#version 300 es
precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;
uniform float Delta;
uniform int writeBlack;

out vec4 Qout;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
    int black = int(mod(p.x + mod(p.y, 2), 2)); 

    int ir = X.x;
    int iy = X.y;
    float r = (0.5 + float(ir))*Delta;
    float y = (0.5 + float(iy))*Delta;

    // Neumann boundary conditions
    ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
    ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
    ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
    ivec2 X_yn = ivec2(ir, max(iy-1, 0));

    // Get pressure stencil
    vec4 Q   = texelFetch(Qair, X,    0);
    Qout = Q;

    if (black)
    {
        vec4 Q_rp = texelFetch(Qair, X_rp, 0);
        vec4 Q_rn = texelFetch(Qair, X_rn, 0);
        vec4 Q_yp = texelFetch(Qair, X_yp, 0);
        vec4 Q_yn = texelFetch(Qair, X_yn, 0);
        float divv = 0.5*(Q_rp.r - Q_rn.r + Q_yp.g - Q_yn.g)/Delta;
        float avgp = 0.5*(Q_rp.a + Q_rn.a + Q_yp.a + Q_yn.a);
        float p = avgp + 0.125*Delta*(Q_rp.a - Q_rn.a)/r - 0.25*Delta*Delta*divv;
        Qout = Q;
        Qout.a = p;
    }
}
`,

'project-black-vertex.shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'project-red-fragment.shader': `#version 300 es
precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;
uniform float Delta;
uniform int writeBlack;

out vec4 Qout;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
    int black = int(mod(p.x + mod(p.y, 2), 2)); 

    int ir = X.x;
    int iy = X.y;
    float r = (0.5 + float(ir))*Delta;
    float y = (0.5 + float(iy))*Delta;

    // Neumann boundary conditions
    ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
    ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
    ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
    ivec2 X_yn = ivec2(ir, max(iy-1, 0));

    // Get pressure stencil
    vec4 Q   = texelFetch(Qair, X,    0);
    Qout = Q;

    if (!black)
    {
        vec4 Q_rp = texelFetch(Qair, X_rp, 0);
        vec4 Q_rn = texelFetch(Qair, X_rn, 0);
        vec4 Q_yp = texelFetch(Qair, X_yp, 0);
        vec4 Q_yn = texelFetch(Qair, X_yn, 0);
        float divv = 0.5*(Q_rp.r - Q_rn.r + Q_yp.g - Q_yn.g)/Delta;
        float avgp = 0.5*(Q_rp.a + Q_rn.a + Q_yp.a + Q_yn.a);
        float p = avgp + 0.125*Delta*(Q_rp.a - Q_rn.a)/r - 0.25*Delta*Delta*divv;
        Qout = Q;
        Qout.a = p;
    }
}
`,

'project-red-vertex.shader': `#version 300 es
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

uniform sampler2D Qair;

uniform int Nr;
uniform int Ny;
uniform float Delta;

out vec4 Qout;

void main()
{
    // Setup local stencil:
    ivec2 X = ivec2(gl_FragCoord.xy);
    int ir = X.x;
    int iy = X.y;

    // Neumann boundary conditions
    ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
    ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
    ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
    ivec2 X_yn = ivec2(ir, max(iy-1, 0));

    // Compute gradient of pressure field
    vec4 Q   = texelFetch(Qair, X,    0);
    vec4 Q_rp = texelFetch(Qair, X_rp, 0);
    vec4 Q_rn = texelFetch(Qair, X_rn, 0);
    vec4 Q_yp = texelFetch(Qair, X_yp, 0);
    vec4 Q_yn = texelFetch(Qair, X_yn, 0);
    float dpdr = 0.5*(Q_rp.a - Q_rn.a)/Delta;
    float dpdy = 0.5*(Q_yp.a - Q_yn.a)/Delta;

    // Update velocity accordingly
    Qout.r = Qin.r - dpdr;
    Qout.g = Qin.g - dpdy;
    Qout.ba = Qin.ba;
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
uniform float dr;
uniform float dy;
uniform int Nr;
uniform int Ny;
uniform int Nraymarch;
uniform float cv; // specific heat capacity
uniform vec3 extinctionMultiplier;
uniform float emissionMultiplier;

uniform sampler2D Qair; // air simulation
uniform sampler2D Qpar; // particles simulation

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

float airTemp(vec4 Qair_)
{
    float rho_air = max(Qair_.x, DENOM_EPS);
    float E_air   = Qair_.y;
    vec2 vrz_air  = Qair_.zw / rho_air;
    float e_air = max(0.0, E_air/rho_air - 0.5*dot(vrz_air, vrz_air));  
    return e_air / cv;  
}

void main()
{
    vec2 pixel = gl_FragCoord.xy;
    vec3 rayPos, rayDir;
    constructPrimaryRay(pixel, rayPos, rayDir);

    vec3 L = vec3(0.0);

    // intersect ray with volume bounds
    float t0, t1;
    if ( boundsIntersect(rayPos, rayDir, volMin, volMax, t0, t1) )
    {
        // Raymarch
        float dl = (t1 - t0) / float(Nraymarch);
        vec3 pMarch = rayPos + (t0+0.5*dl)*rayDir;

        vec3 Transmission = vec3(1.0); // transmittance
        for (int n=0; n<Nraymarch; ++n)
        {   
            // transform pMarch into simulation domain:
            float y = pMarch.y - volMin.y;
            float r = length((pMarch - volCenter).xz);
            if (r<=volRadius)
            {
                int ir = clamp(int(floor(r/dr)), 0, Nr-1);
                int iy = clamp(int(floor(y/dy)), 0, Ny-1);
                vec4 Qair_  = texelFetch(Qair, ivec2(ir, iy), 0);
                vec4 Qpar_  = texelFetch(Qpar, ivec2(ir, iy), 0);
                
                // Emission to blackbody radiation from hot air
                float T_air = airTemp(Qair_);
                vec3 blackbody_color = tempToRGB(T_air);
                float rho_air = max(Qair_.x, DENOM_EPS);
                const float sigma = 5.67e-8; // Stefan–Boltzmann constant
                vec3 emission = Transmission * emissionMultiplier * sigma * pow(T_air, 4.0) * blackbody_color;

                // Absorption by particles ("dust")
                float rho_par = Qpar_.x;
                vec3 extinction = extinctionMultiplier * rho_par;

                Transmission *= exp(-extinction*dl);
                L += emission * dl;
            }
            pMarch += rayDir*dl;
        }        
    }
    
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