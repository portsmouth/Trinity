precision highp float;

uniform sampler2D Qin;

uniform int Nr;
uniform int Ny;

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

    vec4 Q_rp = texelFetch(Qin, X_rp, 0);
    vec4 Q_rn = texelFetch(Qin, X_rn, 0);
    vec4 Q_yp = texelFetch(Qin, X_yp, 0);
    vec4 Q_yn = texelFetch(Qin, X_yn, 0);

    float divv = 0.5 * (Q_rp.r - Q_rn.r + Q_yp.g - Q_yn.g);
    Qout.r = divv;    
}





