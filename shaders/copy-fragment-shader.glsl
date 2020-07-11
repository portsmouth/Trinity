precision highp float;

/////// input buffers ///////
uniform sampler2D Qin;

/////// output buffers ///////
layout(location = 0) out vec4 Qcopy;


void main()
{
    ivec2 frag = ivec2(gl_FragCoord.xy);
    Qcopy = texelFetch(Qin, frag, 0);
}



