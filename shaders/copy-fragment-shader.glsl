precision highp float;

uniform sampler2D Qin;
out vec4 Qcopy;

void main()
{
    ivec2 X = ivec2(gl_FragCoord.xy);
    Qcopy = texelFetch(Qin, X,    0);
}



