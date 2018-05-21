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





