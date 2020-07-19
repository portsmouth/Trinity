precision highp float;

// Geometry
uniform int Nx;
uniform int Ny;
uniform int Nz;
uniform int Ncol;
uniform int W;
uniform int H;
uniform vec3 L; // world-space extents of grid (also the upper right corner in world space)
uniform float dL;

// Physics
uniform float timestep;
uniform float buoyancy;
uniform float gravity;
uniform float radiationLoss;
uniform float T0;             // reference temperature for buoyancy
uniform float Tambient;       // temperature which generates buoyancy which balances gravity

in vec2 v_texcoord;

/////// input buffers ///////
uniform sampler2D Vair_sampler;   // 0, vec3 velocity field
uniform sampler2D Pair_sampler;   // 1, float pressure field
uniform sampler2D Tair_sampler;   // 2, float temperature field
uniform sampler2D debris_sampler; // 3, float debris density field

/////// output buffers ///////
layout(location = 0) out vec4 Vair_output;
layout(location = 1) out vec4 Pair_output;
layout(location = 2) out vec4 Tair_output;
layout(location = 3) out vec4 debris_output;

vec3 mapFragToVs(in ivec2 frag)
{
    // map fragment coord in [W, H] to continuous position of corresponding voxel center in voxel space
    int iu = frag.x;
    int iv = frag.y;
    int row = int(floor(float(iv)/float(Nz)));
    int col = int(floor(float(iu)/float(Nx)));
    int i = iu - col*Nx;
    int j = col + row*Ncol;
    int k = iv - row*Nz;
    return vec3(ivec3(i, j, k)) + vec3(0.5);
}

vec2 slicetoUV(int j, vec3 vsP)
{
    // Given y-slice index j, and continuous voxel space xz-location,
    // return corresponding continuous frag UV for interpolation within this slice
    int row = int(floor(float(j)/float(Ncol)));
    int col = j - row*Ncol;
    vec2 uv_ll = vec2(float(col*Nx)/float(W), float(row*Nz)/float(H));
    float du = vsP.x/float(W);
    float dv = vsP.z/float(H);
    return uv_ll + vec2(du, dv);
}

vec4 interp(in sampler2D S, in vec3 wsP)
{
    vec3 vsP = wsP / dL;
    float pY = vsP.y - 0.5;
    int jlo = clamp(int(floor(pY)), 0, Ny-1); // lower j-slice
    int jhi = clamp(         jlo+1, 0, Ny-1); // upper j-slice
    float flo = float(jhi) - pY;              // lower j fraction
    float fhi = 1.0 - flo;                    // upper j fraction
    vec2 uv_lo = slicetoUV(jlo, vsP);
    vec2 uv_hi = slicetoUV(jhi, vsP);
    vec4 Slo = texture(S, uv_lo);
    vec4 Shi = texture(S, uv_hi);
    return flo*Slo + fhi*Shi;
}

vec3 clampToBounds(in vec3 wsX)
{
    vec3 halfVoxel = vec3(dL);
    return clamp(wsX, halfVoxel, L-halfVoxel);
}

vec3 back_advect(in vec3 wsX, in vec3 vX, float h)
{
    // RK4 integration for position wsX advected backwards through time h:
    vec3 k1 = vX;
    vec3 k2 = interp(Vair_sampler, clampToBounds(wsX - 0.5*h*k1)).xyz;
    vec3 k3 = interp(Vair_sampler, clampToBounds(wsX - 0.5*h*k2)).xyz;
    vec3 k4 = interp(Vair_sampler, clampToBounds(wsX -     h*k3)).xyz;
    return clampToBounds(wsX - h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0);
}

bool isSolidCell(in ivec3 vsPi)
{
    // @todo: expose for generalization
    vec3 vsP = vec3(float(vsPi.x)+0.5, float(vsPi.y)+0.5,float(vsPi.z)+0.5);
    vec3 wsP = vsP*dL;
    vec3 C = L/2.0;
    float r = length(wsP - C);
    return r <= (L.x/2.5);
    //return vsPi.y<=1;
}

void main()
{
    // fragment range over [N*N, N] space
    ivec2 frag = ivec2(gl_FragCoord.xy);
    vec3 vsX = mapFragToVs(frag);
    ivec3 vsXi = ivec3(floor(vsX));

    /*
    if (isSolidCell(vsXi))
    {
        Vair_output = vec4(0.0);
        Pair_output = vec4(0.0);
        Tair_output = vec4(0.0);
        debris_output = vec4(0.0);
    }
    else
    */
    {
        // Apply semi-Lagrangian advection
        vec3 v0 = texelFetch(Vair_sampler, frag, 0).xyz;
        vec3 wsX = vsX*dL;
        vec3 wsXp = back_advect(wsX, v0, timestep);

       //ivec3 vsXpi = ivec3(floor(wsXp/dL));
        //if (isSolidCell(vsXpi))


        vec3  v = interp(Vair_sampler, wsXp).rgb;
        float P = interp(Pair_sampler, wsXp).x;
        float T = interp(Tair_sampler, wsXp).x;
        float debris = interp(debris_sampler, wsXp).x;

        // Apply thermal relaxation to ambient temperature due to "radiation loss"
        float dT = T - Tambient;
        T = Tambient + dT*exp(-radiationLoss);

        // Apply gravity + buoyancy force
        float buoyancy_force = -gravity + gravity*buoyancy*(T - T0);
        v.y += timestep * buoyancy_force;

        Vair_output = vec4(v, 0.0);
        Pair_output = vec4(P);
        Tair_output = vec4(T);
        debris_output = vec4(debris);
    }
}





