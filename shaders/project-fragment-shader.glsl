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





