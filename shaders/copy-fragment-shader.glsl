precision highp float;

uniform sampler2D Qred;
uniform sampler2D Qblack;

out vec4 Qcopy;

void main()
{
    vec2 p = vec2(floor(gl_FragCoord.x), floor(gl_FragCoord.y));
    ivec2 X = ivec2(gl_FragCoord.xy);
    int black = int(mod(p.x + mod(p.y, 2), 2)); 
    if (black)
    {
        Qcopy = texelFetch(Qblack, X,    0);
    }
    else
    {
        Qcopy = texelFetch(Qred, X,    0);
    }
}



