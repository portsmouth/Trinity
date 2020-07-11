
precision highp float;

in vec2 vTexCoord;

layout(location = 0) out vec4 gbuf_rad;

uniform vec2 resolution;

uniform float debrisExtinction;
uniform float blackbodyEmission;
uniform float T0;
uniform float exposure;
uniform float invGamma;

uniform sampler2D Qair;    // air simulation
uniform sampler2D Qdebris; // debris simulation

#define DENOM_EPS 1.0e-7
#define sort2(a,b) { vec3 tmp=min(a,b); b=a+b-tmp; a=tmp; }

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
    if (T_kelvin <= 1000.0) T_kelvin = 1000.0;
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
    ivec2 X = ivec2(gl_FragCoord.xy);

    vec4 Qair_    = texelFetch(Qair, X,    0);
    vec4 Qdebris_ = texelFetch(Qdebris, X,    0);

    float T = Qair_.b;
    vec3 blackbody_color = tempToRGB(T/T0 * 300.0);

    vec3 emission = blackbodyEmission * pow((T-T0)/T0, 4.0) * blackbody_color;
    vec3 sigma = debrisExtinction * Qdebris_.r * vec3(Qdebris_.g, Qdebris_.b, Qdebris_.a);

    vec3 L0 = vec3(0.5, 0.6, 0.9);
    vec3 L = (L0 + emission) * exp(-sigma);

    // apply gamma correction to convert linear RGB to sRGB
    L = pow(L, vec3(invGamma));
    L *= pow(2.0, exposure);

    gbuf_rad = vec4(L, 1.0);
}





