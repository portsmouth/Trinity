precision highp float;

uniform sampler2D Qair;
uniform sampler2D Qdust;

uniform int Nr;
uniform int Ny;
uniform float volRadius;
uniform float volHeight;
uniform float timestep;
uniform float buoyancy;
uniform float dustWeight;
uniform float radiationLoss;
uniform float T0;

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
    int ir = X.x;
    int iy = X.y;

    // Axis boundary condition 
    // (on symmetry axis, the cells are updated to match the adjacent column,
    //  except for the radial velocity, which is zeroed)
    if (ir==0) X.x = 1;
        
    // Current frame air variables (vr, vy, T, p)
    vec4 Q = texelFetch(Qair, X, 0);
    
    // Semi-Lagrangian advection 
    vec2 C = gl_FragCoord.xy;
    vec2 res = vec2(Nr, Ny);
    vec4 Q_advect = texture(Qair, (C - RK4(C))/res);
    float vr = Q_advect.r;
    float vy = Q_advect.g;
    float  T = Q_advect.b;
    float  p = Q_advect.a;

    // Apply cooling due to "radiation loss"
    float DT = T - T0;
    if (DT > 0.0)
    {   
        DT *= exp(-radiationLoss);
        T = T0 + DT;
    }

    // Apply external force
    vec4 Qdust_ = texelFetch(Qdust, X, 0);
    float fy = timestep * (buoyancy*DT - dustWeight*Qdust_.r);
    vy += fy;

    // Apply boundary conditions
    if (ir==0) vr = 0.0;     // axis boundary condition at r=0
    if (iy==0) vy = abs(vy); // solid boundary condition at y=0

    // vy advected + y-force
    Qout.r = vr;
    Qout.g = vy;
    Qout.b = T;
    Qout.a = p;
}





