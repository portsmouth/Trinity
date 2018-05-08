precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;
uniform float Delta;
uniform int writeBlack;

out vec4 Qout;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
    int black = int(mod(p.x + mod(p.y, 2), 2)); 

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
    vec4 Q   = texelFetch(Qair, X,    0);
    Qout = Q;

    if (black)
    {
        vec4 Q_rp = texelFetch(Qair, X_rp, 0);
        vec4 Q_rn = texelFetch(Qair, X_rn, 0);
        vec4 Q_yp = texelFetch(Qair, X_yp, 0);
        vec4 Q_yn = texelFetch(Qair, X_yn, 0);
        float divv = 0.5*(Q_rp.r - Q_rn.r + Q_yp.g - Q_yn.g)/Delta;
        float avgp = 0.5*(Q_rp.a + Q_rn.a + Q_yp.a + Q_yn.a);
        float p = avgp + 0.125*Delta*(Q_rp.a - Q_rn.a)/r - 0.25*Delta*Delta*divv;
        Qout = Q;
        Qout.a = p;
    }
}





