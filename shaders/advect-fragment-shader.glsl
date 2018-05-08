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





