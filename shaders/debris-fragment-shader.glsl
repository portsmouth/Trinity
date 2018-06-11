precision highp float;

uniform sampler2D Qdebris;
uniform sampler2D Qair;

uniform int Nr;
uniform int Ny;
uniform float volRadius;
uniform float volHeight;
uniform float timestep;

in vec2 v_texcoord;
out vec4 Qout;

vec2 RK4(vec2 p)
{
    float h = timestep;
    vec2 res = vec2(Nr, Ny);

    vec2 uv1 = p/res; 
    vec2 k1 = texture(Qair, uv1).xy;

    vec2 uv2 = (p - 0.5*h*k1)/res; uv2.x = abs(uv2.x);
    vec2 k2 = texture(Qair, uv2).xy;

    vec2 uv3 = (p - 0.5*h*k2)/res; uv3.x = abs(uv3.x);
    vec2 k3 = texture(Qair, uv3).xy;

    vec2 uv4 = (p - h*k3)/res; uv4.x = abs(uv3.x);
    vec2 k4 = texture(Qair, uv4).xy;

    return h/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    vec4 Qair_ = texelFetch(Qair, X, 0);
    float vr = Qair_.r; // in voxels/timestep
    float vy = Qair_.g; // in voxels/timestep

     // Semi-Lagrangian advection 
    vec2 C = gl_FragCoord.xy;
    vec2 res = vec2(Nr, Ny);
    Qout = texture(Qdebris, (C - RK4(C))/res);
}





