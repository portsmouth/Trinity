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

uniform float dt;
uniform float dr;
uniform float dz;
uniform float gammaMO; // (gamma-1)
uniform float g;        // gravitational acceleration

out vec4 Qnext;


// energy injection term
float G(in vec2 X)
{
	// @todo: inject non-zero energy in the first timestep
	//        within a sphere (which encloses at least one voxel center)
	// @todo: for now, place this detonation point in the grid center voxel
	return 0.0;
}

// Equation of state
float pressure(float rho, float E)
{
	float e = max(0.0, E - 0.5*(vr*vr + vz*vz));	
	return gammaMO * rho * e;
}

// r-flux
vec4 R(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(q0, q1);
	vec4 R_;
	R_.x = q2;
	R_.y = (q1 + p) * q2/q0;
	R_.z = p + q2*q2/q0;
	R_.w = q2*q3 / q0;
	return R;
}

// z-flux
vec4 Z(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(q0, q1);
	vec4 Z_;
	Z_.x = q3;
	Z_.y = (q1 + p) * q3/q0;
	Z_.z = q2*q3/q0;
	Z_.w = p + q3*q3 / q0;
	return R;
}


// source term
vec4 S(in vec4 Q_, in vec2 X)
{
	float r = (0.5 + X.x)*dr;
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(q0, q1);
	vec4 S_;
	S_.x = - q2/r;
	S_.y = -(q1 + p)*q2/(q0*r) - g*q3/q0 + G(X)/r;
	S_.z = -q2*q2 / (r*q0);
	S_.w = -q2*q3/(r*q0) - q0*g;
	return S_;
}	


void main()
{
	// Setup local stencil:
	vec2 size = textureSize(U_v);
	int Nr = size.x;
	int Nz = size.y;
	ivec2 X = gl_FragCoord.xy;
	int ir = X.x;
	int iz = X.y;
	ivec2 X_rp = ivec2(min(ir+1, Nr-1), iz);
	ivec2 X_rn = ivec2(max(ir-1, 0),    iz);
	ivec2 X_zp = ivec2(ir, min(iz+1, Nz-1));
	ivec2 X_zn = ivec2(ir, max(iz-1, 0));

	// Update Q via Lax-Friedrichs method:
	float lambda_r = dt/dr;
	float lambda_z = dt/dz;
	vec4 Q_ = texelFetch(Q, X).rgb;
	vec4 Q_rp = texelFetch(Q, X_rp).rgb;
	vec4 Q_rn = texelFetch(Q, X_rn).rgb;
	vec4 Q_zp = texelFetch(Q, X_zp).rgb;
	vec4 Q_zn = texelFetch(Q, X_zn).rgb;
	vec4 Qavg = 0.25 * (Q_rp + Q_rn + Q_zp + Q_zn);
	Qnext = Qavg - 0.5*lambda_r*(R(Q_rp)-R(Q_rn)) - 0.5*lambda_z*(R(Q_zp)-R(Q_zn)) + dt*S(Q_, X);
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


void constructPrimaryRay(in vec2 pixel, inout vec4 rnd,
                         inout vec3 primaryStart, inout vec3 primaryDir)
{
    // Compute world ray direction for given (possibly jittered) fragment
    vec2 ndc = -1.0 + 2.0*(pixel/resolution.xy);
    float fh = tan(0.5*radians(camFovy)); // frustum height
    float fw = camAspect*fh;
    vec3 s = -fw*ndc.x*camX + fh*ndc.y*camY;
    primaryDir = normalize(camDir + s);
    primaryStart = camPos;
}


void main()
{
    vec2 pixel = gl_FragCoord.xy;

    vec3 RGB = vec3(1.0);
    
    gbuf_rad = vec4(RGB, 1.0);
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