precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;
uniform float Delta;

out vec4 Qout;

void main()
{
    // Setup local stencil:
    ivec2 X = ivec2(gl_FragCoord.xy);
    int ir = X.x;
    int iy = X.y;

    // Neumann boundary conditions
    ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
    ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
    ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
    ivec2 X_yn = ivec2(ir, max(iy-1, 0));

    // Compute gradient of pressure field
    vec4 Q    = texelFetch(Qin, X,    0);
    vec4 Q_rp = texelFetch(Qin, X_rp, 0);
    vec4 Q_rn = texelFetch(Qin, X_rn, 0);
    vec4 Q_yp = texelFetch(Qin, X_yp, 0);
    vec4 Q_yn = texelFetch(Qin, X_yn, 0);
    float dpdr = 0.5*(Q_rp.a - Q_rn.a);
    float dpdy = 0.5*(Q_yp.a - Q_yn.a);

    // Update velocity accordingly
    Qout.r = Q.r - dpdr;
    Qout.g = Q.g - dpdy;
    Qout.ba = Q.ba;
}





