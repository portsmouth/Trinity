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





