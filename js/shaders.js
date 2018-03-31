var Shaders = {

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

'solve-energy-fragment-shader': `#version 300 es
precision highp float;

uniform sampler2D Q;

uniform int Nr;
uniform int Ny;
uniform float dt;
uniform float dr;
uniform float dy;
uniform float gammaMO; // (gamma-1)
uniform float g;        // gravitational acceleration

out vec4 Qnext;

#define DENOM_EPS 1.0e-7

// energy injection term
float G(in ivec2 X)
{
	// @todo: inject non-zero energy in the first timestep
	//        within a sphere (which encloses at least one voxel center)
	// @todo: for now, place this detonation point in the grid center voxel
	return 0.0;
}

// Equation of state
float pressure(in vec4 Q_)
{
	float rho = Q_.x; // rho
	float E   = Q_.y; // total energy density
	float vr  = Q_.z / max(DENOM_EPS, rho);
	float vz  = Q_.w / max(DENOM_EPS, rho);
	// p = (gamma-1) * rho * e, where e = E/rho - v^2/2 (specific internal energy)
	float p = gammaMO * max(0.0, E - 0.5*rho*(vr*vr + vz*vz));	
	return p;
}

// r-flux
vec4 R(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(Q_);
	vec4 R_;
	R_.x = q2;
	R_.y = (q1 + p) * q2/q0;
	R_.z = p + q2*q2/q0;
	R_.w = q2*q3 / q0;
	return R_;
}

// y-flux
vec4 Y(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(Q_);
	vec4 Y_;
	Y_.x = q3;
	Y_.y = (q1 + p) * q3/q0;
	Y_.z = q2*q3/q0;
	Y_.w = p + q3*q3 / q0;
	return Y_;
}


// source term
vec4 S(in vec4 Q_, in ivec2 X)
{
	float r = (0.5 + float(X.x))*dr;
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(Q_);
	vec4 S_;
	S_.x = -q2/r;
	S_.y = -(q1 + p)*q2/(q0*r) - g*q3/q0 + G(X)/r;
	S_.z = -q2*q2 / (r*q0);
	S_.w = -q2*q3/(r*q0) - q0*g;
	return S_;
}	


void main()
{
	// Setup local stencil:
	ivec2 X = ivec2(gl_FragCoord.xy);
	int ir = X.x;
	int iy = X.y;
	ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
	ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
	ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
	ivec2 X_yn = ivec2(ir, max(iy-1, 0));

	// Update Q via Lax-Friedrichs method:
	float lambda_r = dt/dr;
	float lambda_y = dt/dy;
	vec4 Q_   = texelFetch(Q, X,    0);
	vec4 Q_rp = texelFetch(Q, X_rp, 0);
	vec4 Q_rn = texelFetch(Q, X_rn, 0);
	vec4 Q_yp = texelFetch(Q, X_yp, 0);
	vec4 Q_yn = texelFetch(Q, X_yn, 0);
	vec4 Qavg = 0.25 * (Q_rp + Q_rn + Q_yp + Q_yn);

	Qnext = Qavg - 0.5 * lambda_r * ( R(Q_rp) - R(Q_rn) ) 
	             - 0.5 * lambda_y * ( Y(Q_yp) - Y(Q_yn) );
	             + dt * S(Q_, X);
}
`,

'solve-energy-vertex-shader': `#version 300 es
precision highp float;

in vec3 Position;
in vec2 TexCoord;

void main() 
{
    gl_Position = vec4(Position, 1.0);
}
`,

'solve-velocity-fragment-shader': `#version 300 es

`,

'solve-velocity-vertex-shader': `#version 300 es
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

'solver-fragment-shader': `#version 300 es

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

uniform sampler2D Q;

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

    // intersect ray with volume bounds
    float t0, t1;
    if ( boundsIntersect(rayPos, rayDir, volMin, volMax, t0, t1) )
    {
        // Raymarch
        float dl = (t1 - t0) / float(Nraymarch);
        vec3 x = rayPos + (t0+0.5*dl)*rayDir;
        for (int n=0; n<Nraymarch; ++n)
        {   
            // transform x into simulation domain:
            float y = x.y - volMin.y;
            float r = length((x - volCenter).xz);
            if (r<=volRadius)
            {
                int ir = clamp(int(floor(r/dr)), 0, Nr-1);
                int iy = clamp(int(floor(y/dy)), 0, Ny-1);
                vec4 Q_  = texelFetch(Q, ivec2(ir, iy), 0);

                float rho = max(Q_.x, DENOM_EPS);
                float E   = Q_.y;
                float vr = Q_.z / rho; // rho * vr
                float vy = Q_.w / rho; // rho * vy
                float e = max(0.0, E - 0.5*(vr*vr + vy*vy));  
                float T = e / cv;
                vec3 blackbody_color = tempToRGB(T);

                float emission = 1.0e-6 * pow(T, 4.0);
                L += emission * blackbody_color;

            }

            x += rayDir*dl;
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