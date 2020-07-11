

/** @constructor
*/
var Solver = function()
{
    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource({
        'advect'         : {'v': 'advect-vertex-shader',       'f': 'advect-fragment-shader'       },
        'div'            : {'v': 'div-vertex-shader',          'f': 'div-fragment-shader'          },
        'project'        : {'v': 'project-vertex-shader',      'f': 'project-fragment-shader'      },
        'copy'           : {'v': 'copy-vertex-shader',         'f': 'copy-fragment-shader'         },
        'update'         : {'v': 'update-vertex-shader',       'f': 'update-fragment-shader'       },
        'debris'         : {'v': 'debris-vertex-shader',       'f': 'debris-fragment-shader'       }
    });

    this.fbo = new GLU.RenderTarget();
    this.quadVbo = this.createQuadVbo();
    this.compileShaders();

    // Defaults
    this.timestep = 0.3;
    this.blastHeight = 0.1;
    this.blastRadius = 0.05;
    this.blastTemperature = 100.0;  // initial temperature of fireball relative to ambient
    this.blastVelocity = 50.0;    // outward blast speed in voxels/timestep
    this.debrisHeight = 0.2;       // maximum height of dust layer, as a fraction of domain height
    this.debrisFalloff = 0.01;      // fall-off exponent within dust layer
    this.T0 = 1.0;                 // nominal reference temperature for buoyancy
    this.buoyancy = 0.003;         // initial buoyancy (thermal expansion coeff. of air)
    this.expansion = 0.001;
    this.gravity = -1.0;
    this.radiationLoss = 0.01;    // radiation loss rate (per timestep fractional absorption)

    // Initialize solver
    this.resize(128);
}

Solver.prototype.compileShaders = function()
{
	this.advect_program       = new GLU.Shader('advect',        this.shaderSources, null);
    this.project_program      = new GLU.Shader('project',       this.shaderSources, null);
    this.div_program          = new GLU.Shader('div',           this.shaderSources, null);
    this.copy_program         = new GLU.Shader('copy',          this.shaderSources, null);
    this.update_program       = new GLU.Shader('update',        this.shaderSources, null);
    this.debris_program       = new GLU.Shader('debris',        this.shaderSources, null);
}

Solver.prototype.resize = function(N)
{
    this.N = N;
    this.reset();
}

Solver.prototype.reset = function()
{
	this.initialize();
}


Solver.prototype.getDomain = function()
{
	return this.domain;
}

Solver.prototype.smoothSphere = function(x, c, R)
{
	// x = [r, y]
	let rho = Math.sqrt(Math.pow(x[0]-c[0], 2.0) + Math.pow(x[1]-c[1], 2.0) + Math.pow(x[2]-c[2], 2.0));
	let t = rho/R;
	if (t<=1.0) return Math.max(0.0, 1.0 - t*t*(3.0 - 2.0*t));
	return 0.0;
}

Solver.prototype.radialFlow = function(x, c, R, V)
{
    // x = [r, y]
	let rho = Math.sqrt(Math.pow(x[0]-c[0], 2.0) + Math.pow(x[1]-c[1], 2.0) + Math.pow(x[2]-c[2], 2.0));
	let t = rho/R;
    let v = (t<=1.0) ? V * Math.max(0.0, 1.0 - t*t*(3.0 - 2.0*t)) : 0.0;
    let vx = v * (x[0]-c[0])/rho
    let vy = v * (x[1]-c[1])/rho
    let vz = v * (x[2]-c[2])/rho
    return [vx, vy, vz];
}

Solver.prototype.mapVsToFrag = function(vsP, N)
{
    let i = vsP[0];
    let j = vsP[1];
    let k = vsP[2];
    let ui = N*j + i;
    let vi = k;
    return [ui, vi];
}

Solver.prototype.mapFragtoBufferIndex = function(frag, N, numChannels)
{
    let iu = frag[0];
    let iv = frag[1];
    let width = N*N;
    let index = iv*width + iu;
    return numChannels*index;
}

Solver.prototype.initialize = function()
{
    let N = this.N;

    let Vair0    = new Float32Array(N*N * N * 4);  // air velocity field
    let divVair0 = new Float32Array(N*N * N * 4);  // air velocity divergence field
    let Pair0    = new Float32Array(N*N * N * 4);  // air pressure field
    let Tair0    = new Float32Array(N*N * N * 4);  // air temperature field
    let debris0  = new Float32Array(N*N * N * 4);  // debris extinction field

    let L = 100.0; // domain size in world units
    this.dL = L/N;
    let dL = this.dL;  // voxel size in world units

    noise.seed(Math.random());

    let blast_center = [0.5*L, L*this.blastHeight, 0.5*L];

    this.Tambient = this.T0 + 1.0/this.buoyancy; // ambient temp which balances buoyancy and gravity

    for (let iy=0; iy<N; ++iy)
    {
        for (var iz=0; iz<N; ++iz)
         {
            for (var ix=0; ix<N; ++ix)
            {
                let x = (0.5 + ix)*dL;
                let y = (0.5 + iy)*dL;
                let z = (0.5 + iz)*dL;
                let wsP = [x,y,z];      // world space location of voxel center
                let vsP = [ix, iy, iz]  // voxel space location of lower-left corner (== voxel index)

                let r = Math.sqrt(Math.pow(x-blast_center[0], 2.0) + Math.pow(y-blast_center[1], 2.0)+ Math.pow(z-blast_center[2], 2.0));

                // initialize temperature profile
                let T = this.Tambient * (1.0 + this.blastTemperature*this.smoothSphere(wsP, blast_center, this.blastRadius*L));

                // initialize velocity to a spherical blast with a finite region:
                V = this.radialFlow(wsP, blast_center, this.blastRadius*L, this.blastVelocity);

                let frag = this.mapVsToFrag(vsP, N);
                let b = this.mapFragtoBufferIndex(frag, N, 4);

                // Initialize air variables
                Vair0[b+0]  = V[0];
                Vair0[b+1]  = V[1];
                Vair0[b+2]  = V[2];
                divVair0[b] = 0.0;
                Pair0[b]    = 0.0;
                Tair0[b]    = T;

                // Initialize debris density field
                let density = y < this.debrisHeight*L ? Math.exp(-3.0*this.debrisFalloff*y/(this.debrisHeight*L)) : 0.0;
                debris0[b+0] = density;
                debris0[b+1] = density;
                debris0[b+2] = density;
            }
        }
    }

    // Initial velocity texture
    this.Vair = [ new GLU.Texture(N*N, N, 4, true, true, true, Vair0),
                  new GLU.Texture(N*N, N, 4, true, true, true, Vair0) ];

    // Initial velocity divergence texture
    this.divVair = new GLU.Texture(N*N, N, 4, true, true, true, divVair0);

    // Initial pressure texture
    this.Pair = [ new GLU.Texture(N*N, N, 4, true, true, true, Pair0),
                  new GLU.Texture(N*N, N, 4, true, true, true, Pair0) ];

    // Initial temperature texture
    this.Tair = [ new GLU.Texture(N*N, N, 4, true, true, true, Tair0),
                  new GLU.Texture(N*N, N, 4, true, true, true, Tair0) ];

    // Initial debris texture
    this.debris = [ new GLU.Texture(N*N, N, 4, true, true, true, debris0),
                    new GLU.Texture(N*N, N, 4, true, true, true, debris0) ];

    this.frame = 0;

    this.domain = { N: N,
                    L: this.L,
                    dL: this.dL,
                    boundsMin: [0.0, 0.0, 0.0 ],
                    boundsMax: [L, L, L],
                    center:    [L/2.0, L/2.0, L/2.0],
                    T0: this.T0
                  };
    }


Solver.prototype.createQuadVbo = function()
{
	let gl = GLU.gl;

    var vbo = new GLU.VertexBuffer();
    vbo.addAttribute("Position", 3, gl.FLOAT, false);
    vbo.addAttribute("TexCoord", 2, gl.FLOAT, false);
    vbo.init(4);
    vbo.copy(new Float32Array([
         1.0,  1.0, 0.0, 1.0, 1.0,
        -1.0,  1.0, 0.0, 0.0, 1.0,
        -1.0, -1.0, 0.0, 0.0, 0.0,
         1.0, -1.0, 0.0, 1.0, 0.0
    ]));
    return vbo;
}

// air velocity texture
Solver.prototype.getVair = function()
{
	return this.Vair[0];
}

// air temperature texture
Solver.prototype.getTair = function()
{
	return this.Tair[0];
}

// debris density texture
Solver.prototype.getDebris = function()
{
	return this.debris[0];
}

Solver.prototype.step = function()
{
    console.log("Running frame: ", this.frame);

    let gl = GLU.gl;

    let N = this.N;
    gl.viewport(0, 0, N*N, N);
    this.quadVbo.bind();

    // Run air force/advection step, writing 0 -> 1
    {
        this.advect_program.bind();
        this.advect_program.uniformI("N",             this.domain.N);
        this.advect_program.uniformF("dL",            this.domain.dL);
        this.advect_program.uniformF("timestep",      this.timestep);
        this.advect_program.uniformF("buoyancy",      this.buoyancy);
        this.advect_program.uniformF("gravity",      this.gravity);
        this.advect_program.uniformF("radiationLoss", this.radiationLoss * 1.0e-2);
        this.advect_program.uniformF("Tambient",      this.Tambient);

        this.fbo.bind();
        this.fbo.drawBuffers(3);
        this.fbo.attachTexture(this.Vair[1], 0); // write to  Vair[1]
        this.fbo.attachTexture(this.Pair[1], 1); // write to  Pair[1]
        this.fbo.attachTexture(this.Tair[1], 2); // write to  Tair[1]
        this.Vair[0].bind(0);                    // read from Vair[0]
        this.Pair[0].bind(1);                    // read from Pair[0]
        this.Tair[0].bind(2);                    // read from Tair[0]
        this.advect_program.uniformTexture("Vair_sampler", this.Vair[0]);
        this.advect_program.uniformTexture("Pair_sampler", this.Pair[0]);
        this.advect_program.uniformTexture("Tair_sampler", this.Tair[0]);
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.unbind();

        // copy temperature 1 -> 0
        this.copy_program.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Tair[0], 0); // write to  Tair[0]
        this.Tair[1].bind(0);                    // read from Tair[1]
        this.copy_program.uniformTexture("Qin", this.Tair[1]);
        this.quadVbo.draw(this.copy_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();

    }

    gl.bindTexture(gl.TEXTURE_2D, null);

    // Compute velocity divergence
    {
        this.div_program.bind();
        this.div_program.uniformI("N",             this.domain.N);
        this.div_program.uniformF("dL",            this.domain.dL);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.divVair, 0); // write to divVair
        this.Vair[1].bind(0);                    // read from Vair[1]
        this.div_program.uniformTexture("Vair_sampler", this.Vair[1]);
        this.quadVbo.draw(this.div_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

    gl.bindTexture(gl.TEXTURE_2D, null);

    // Run air pressure projection step
    {
        this.project_program.bind();
        this.project_program.uniformI("N",             this.domain.N);
        this.project_program.uniformF("dL",            this.domain.dL);
        this.project_program.uniformF("timestep", this.timestep);
        this.project_program.uniformF("expansion", this.expansion);

        // Update pressure field by Jacobi iteration
        // (using last frame pressure as a warm-start)
        // Sequence:
        //   1  ->  0  (read black from 1, write red   to 0)
        //   0  ->  1  (read red   from 0, write black to 1)
        const numIterations = 8;
        for (let n=0; n<numIterations; ++n)
        {
            this.fbo.bind();
            this.fbo.drawBuffers(1);
            this.fbo.attachTexture(this.Pair[0], 0);          // write to Pair[0]
            this.Pair[1].bind(0);                             // read from Pair[1]
            this.Tair[0].bind(1);                             // read from Tair[0]
            this.divVair.bind(2);                             // read divVair
            this.project_program.uniformTexture("Pair_sampler", this.Pair[1]);
            this.project_program.uniformTexture("Tair_sampler", this.Tair[0]);
            this.project_program.uniformTexture("divVair_sampler", this.divVair);
            this.quadVbo.draw(this.project_program, gl.TRIANGLE_FAN);
            this.fbo.detachTexture(0);
            this.fbo.unbind();

            gl.bindTexture(gl.TEXTURE_2D, null);

            this.fbo.bind();
            this.fbo.drawBuffers(1);
            this.fbo.attachTexture(this.Pair[1], 0);          // write to Pair[1]
            this.Pair[0].bind(0);                             // read from Pair[0]
            this.Tair[0].bind(1);                             // read from Tair[0]
            this.divVair.bind(2);                             // read divVair
            this.project_program.uniformTexture("Pair_sampler", this.Pair[0]);
            this.project_program.uniformTexture("Tair_sampler", this.Tair[0]);
            this.project_program.uniformTexture("divVair_sampler", this.divVair);
            this.quadVbo.draw(this.project_program, gl.TRIANGLE_FAN);
            this.fbo.detachTexture(0);
            this.fbo.unbind();

            gl.bindTexture(gl.TEXTURE_2D, null);
        }
    }

    // Update air velocity from 1 -> 0
    {
        this.update_program.bind();
        this.update_program.uniformI("N",             this.domain.N);
        this.update_program.uniformF("dL",            this.domain.dL);
        this.update_program.uniformF("timestep", this.timestep);

        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Vair[0], 0);          // write to  Vair[0]
        this.Vair[1].bind(0);                             // read from Vair[1]
        this.Pair[1].bind(1);                             // read from Pair[1]
        this.update_program.uniformTexture("Vair_sampler", this.Vair[1]);
        this.update_program.uniformTexture("Pair_sampler", this.Pair[1]);
        this.quadVbo.draw(this.update_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

    gl.bindTexture(gl.TEXTURE_2D, null);
    this.fbo.detachTexture(0);
    this.fbo.unbind();

    // Advect debris in air flow
    {
        // write debris 0 -> 1
        this.debris_program.bind();
        this.debris_program.uniformI("N",             this.domain.N);
        this.debris_program.uniformF("dL",            this.domain.dL);
        this.debris_program.uniformF("timestep", this.timestep);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.debris[1], 0); // write to  debris[1]
        this.debris[0].bind(0);                    // read from debris[0]
        this.Vair[1].bind(1);                      // read from Vair[1]
        this.debris_program.uniformTexture("debris_sampler", this.debris[0]);
        this.debris_program.uniformTexture("Vair_sampler",   this.Vair[1]);
        this.quadVbo.draw(this.debris_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();

        gl.bindTexture(gl.TEXTURE_2D, null);

        // copy debris 1 -> 0
        this.copy_program.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.debris[0], 0); // write to  debris[0]
        this.debris[1].bind(0);                    // read from debris[1]
        this.copy_program.uniformTexture("Qin", this.debris[1]);
        this.quadVbo.draw(this.copy_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

    gl.bindTexture(gl.TEXTURE_2D, null);

    /*
    this.fbo.bind();
    let status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    switch (status)
    {
    	case gl.FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         console.log("FRAMEBUFFER_INCOMPLETE_ATTACHMENT");         break;
    	case gl.FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: console.log("FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT"); break;
    	case gl.FRAMEBUFFER_INCOMPLETE_DIMENSIONS:         console.log("FRAMEBUFFER_INCOMPLETE_DIMENSIONS");         break;
    	case gl.FRAMEBUFFER_UNSUPPORTED:                   console.log("FRAMEBUFFER_UNSUPPORTED");                   break;
    	case gl.FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:        console.log("FRAMEBUFFER_INCOMPLETE_MULTISAMPLE");        break;
    	case gl.RENDERBUFFER_SAMPLES:                      console.log("RENDERBUFFER_SAMPLES");                      break;
    	default: break;
    }

	this.fbo.unbind();
    gl.bindTexture(gl.TEXTURE_2D, null);
    */

    this.frame++;
}
