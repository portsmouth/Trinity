precision highp float;

uniform sampler2D Qin;
uniform sampler2D div;

uniform int Nr;
uniform int Ny;
uniform float volRadius;
uniform float volHeight;
uniform float timestep;
uniform float expansion;

out vec4 Qout;

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    int ir = X.x;
    int iy = X.y;
    
    // Neumann boundary conditions
    ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
    ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
    ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
    ivec2 X_yn = ivec2(ir, max(iy-1, 0));

    // Get pressure stencil
    vec4 Q = texelFetch(Qin, X, 0);
    Qout = Q;

    // Read divergence
    float divv = texelFetch(div, X, 0).x;

    // Introduce local expansion due to heated fluid
    float  T = Q.b;
    float phi = timestep/float(Nr) * expansion * T;
    divv -= phi;

    // Update local pressure according to Poisson equation (in cylindrical polars)
    vec4 Q_rp = texelFetch(Qin, X_rp, 0);
    vec4 Q_rn = texelFetch(Qin, X_rn, 0);
    vec4 Q_yp = texelFetch(Qin, X_yp, 0);
    vec4 Q_yn = texelFetch(Qin, X_yn, 0);
    float avgp = 0.25*(Q_rp.a + Q_rn.a + Q_yp.a + Q_yn.a);
    float r = (0.5 + float(ir));
    float pressure = avgp - 0.25*divv + 0.125*(Q_rp.a - Q_rn.a)/r;
    Qout.a = pressure; 
}





