var Shaders = {

'advect-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;
uniform float Delta;
uniform float g;        // gravitational acceleration (in -y dir)
uniform float beta;     // buoyancy

in vec2 v_texcoord;
out vec4 Qout;

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    vec4 Q = texelFetch(Qin, X,    0);

    float vr = Q.r; // in voxels/timestep
    float vy = Q.g; // in voxels/timestep
    float  T = Q.b;
    float  p = Q.a;

    // Semi-Lagrangian advection 
    float u_advect = clamp(v_texcoord.x - vr/float(Nr), 0.0, 1.0);
    float v_advect = clamp(v_texcoord.y - vy/float(Ny), 0.0, 1.0);
    vec4 Q_advect = texture(Qin, vec2(u_advect, v_advect));

    float vr_advect = Q_advect.r;
    float vy_advect = Q_advect.g;
    float  T_advect = Q_advect.b;
    float  p_advect = Q_advect.a;

    // vr advect
    Qout.r = vr_advect;
    
    // Solid boundary condition at y=0
    int iy = X.y;
    if (iy==0)
    {
        Qout.g = 0.0;
    }
    else
    {
        // vy advect + force
        Qout.g = vy_advect - g + beta*T_advect; // - T0
    }

    // T advect
    Qout.b = T_advect;

    // p copy
    Qout.a = p_advect;
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

uniform sampler2D Qin;
out vec4 Qcopy;

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    Qcopy = texelFetch(Qin, X,    0);
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

uniform sampler2D Qdebris;
uniform sampler2D Qair;

uniform int Nr;
uniform int Ny;

in vec2 v_texcoord;
out vec4 Qout;

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    vec4 Qair_ = texelFetch(Qair, X, 0);
    float vr = Qair_.r; // in voxels/timestep
    float vy = Qair_.g; // in voxels/timestep

    // Semi-Lagrangian advection 
    float u_advect = clamp(v_texcoord.x - vr/float(Nr), 0.0, 1.0);
    float v_advect = clamp(v_texcoord.y - vy/float(Ny), 0.0, 1.0);
    Qout = texture(Qdebris, vec2(u_advect, v_advect));
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

'merge-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Qred;
uniform sampler2D Qblack;
out vec4 Qmerge;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
    int black = int(mod(p.x + mod(p.y, 2.0), 2.0)); 
    if (black==0)
    {
        Qmerge = texelFetch(Qred, X, 0);
    }
    else
    {
        Qmerge = texelFetch(Qblack, X, 0);
    }
}
`,

'merge-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'project-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;
uniform float Delta;

out vec4 Qout;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
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
    vec4 Q   = texelFetch(Qin, X,    0);
    Qout = Q;

    vec4 Q_rp = texelFetch(Qin, X_rp, 0);
    vec4 Q_rn = texelFetch(Qin, X_rn, 0);
    vec4 Q_yp = texelFetch(Qin, X_yp, 0);
    vec4 Q_yn = texelFetch(Qin, X_yn, 0);
    float divv = 0.5*(Q_rp.r - Q_rn.r + Q_yp.g - Q_yn.g)/Delta;
    float avgp = 0.25*(Q_rp.a + Q_rn.a + Q_yp.a + Q_yn.a);
    float pressure = avgp - 0.25*Delta*Delta*divv + 0.125*Delta*(Q_rp.a - Q_rn.a)/r;
    Qout = Q;
    Qout.a = pressure;    
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

uniform sampler2D Qin;

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
    vec4 Q    = texelFetch(Qin, X,    0);
    vec4 Q_rp = texelFetch(Qin, X_rp, 0);
    vec4 Q_rn = texelFetch(Qin, X_rn, 0);
    vec4 Q_yp = texelFetch(Qin, X_yp, 0);
    vec4 Q_yn = texelFetch(Qin, X_yn, 0);
    float dpdr = 0.5*(Q_rp.a - Q_rn.a);
    float dpdy = 0.5*(Q_yp.a - Q_yn.a);

    // Update velocity accordingly
    Qout.r = Q.r - dpdr;
    Qout.g = Q.g - dpdy;
    Qout.ba = Q.ba;
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
uniform float Delta; // (cubic) voxel size
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
                int ir = clamp(int(floor(r/Delta)), 0, Nr-1);
                int iy = clamp(int(floor(y/Delta)), 0, Ny-1);
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
                //vec3 blackbody_color = tempToRGB(T * 300.0 * tempMultiplier);
                vec3 emission = Tr * emissionMultiplier * 0.001 * T * vec3(1.0, 0.5, 0.1); //blackbody_color;

                L += emission * dl;
            }
            pMarch += rayDir*dl;
        }  
    }

    L += Tr * Lbackground;
    
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