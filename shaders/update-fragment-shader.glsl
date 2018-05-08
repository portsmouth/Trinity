precision highp float;

uniform sampler2D Qair;

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
    vec4 Q   = texelFetch(Qair, X,    0);
    vec4 Q_rp = texelFetch(Qair, X_rp, 0);
    vec4 Q_rn = texelFetch(Qair, X_rn, 0);
    vec4 Q_yp = texelFetch(Qair, X_yp, 0);
    vec4 Q_yn = texelFetch(Qair, X_yn, 0);
    float dpdr = 0.5*(Q_rp.a - Q_rn.a)/Delta;
    float dpdy = 0.5*(Q_yp.a - Q_yn.a)/Delta;

    // Update velocity accordingly
    Qout.r = Qin.r - dpdr;
    Qout.g = Qin.g - dpdy;
    Qout.ba = Qin.ba;
}





