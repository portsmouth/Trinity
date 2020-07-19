

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
    this.timestep = 2.0;
    this.NprojSteps = 16;
    this.blastHeight = 0.1;
    this.blastRadius = 0.15;
    this.blastTemperature = 500.0;   // initial temperature of fireball relative to ambient
    this.blastVelocity = 100.0;        // outward blast speed in voxels/timestep
    this.debrisHeight = 0.1;        // maximum height of dust layer, as a fraction of domain height
    this.debrisFalloff = 0.5;       // fall-off exponent within dust layer
    this.T0       = 266.0;                  // nominal reference temperature for buoyancy
    this.Tambient = 270.0;                  // nominal reference temperature for buoyancy
    this.buoyancy = 0.2;           // initial buoyancy (thermal expansion coeff. of air)
    this.expansion = 0.0;
    this.gravity = 0.003; //0.5;
    this.radiationLoss = 0.0;        // radiation loss rate (per timestep fractional absorption)

    // Initialize solver
    this.resize(128, 512, 128);
}

Solver.prototype.compileShaders = function()
{
	this.advect_program       = new GLU.Shader('advect',        this.shaderSources, null);
    this.project_program      = new GLU.Shader('project',       this.shaderSources, null);
    this.div_program          = new GLU.Shader('div',           this.shaderSources, null);
    this.copy_program         = new GLU.Shader('copy',          this.shaderSources, null);
    this.update_program       = new GLU.Shader('update',        this.shaderSources, null);
}

Solver.prototype.resize = function(Nx, Ny, Nz)
{
    this.Nx = Nx;
    this.Ny = Ny;
    this.Nz = Nz;

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

Solver.prototype.mapVsToFrag = function(vsP, Ncol, Nx, Ny, Nz)
{
    let i = vsP[0];
    let j = vsP[1];
    let k = vsP[2];
    let row = Math.floor(j/Ncol);
    let col = j - row*Ncol;
    let iu = col*Nx + i;
    let iv = row*Nz + k;
    return [iu, iv];
}

Solver.prototype.mapFragtoBufferIndex = function(frag, W, numChannels)
{
    let iu = frag[0];
    let iv = frag[1];
    let index = iv*W + iu;
    return numChannels*index;
}

Solver.prototype.initialize = function()
{
    let Nx = this.Nx;
    let Ny = this.Ny;
    let Nz = this.Nz;

    let W = 16384; // maximum texture width available in WebGL2
    let Ncol = Math.floor(W/Nx);
    if (Ncol<1)
        GLU.fail("Resolution too high: Nx=", Nx);
    let Nrow = Math.ceil(Ny/Ncol);
    let H = Nrow*Nz;
    if (H > W)
        GLU.fail("Resolution too high: Ny=", Ny);
    this.W = W;
    this.H = H;

    let Vair0    = new Float32Array(W*H*4);  // air velocity field
    let divVair0 = new Float32Array(W*H*4);  // air velocity divergence field
    let Pair0    = new Float32Array(W*H*4);  // air pressure field
    let Tair0    = new Float32Array(W*H*4);  // air temperature field
    let debris0  = new Float32Array(W*H*4);  // debris extinction field

    let dL = 1.0; // voxel size in world units
    this.dL = dL;

    let L = [dL*Nx, dL*Ny, dL*Nz];
    this.L = L;
    let lengthScale = Math.min(L[0], L[1], L[2]);

    let blast_center = [0.5*L[0], L[1]*this.blastHeight, 0.5*L[2]];

    //this.Tambient = this.T0 + 1.0/this.buoyancy; // ambient temp which balances buoyancy and gravity

    for (let iy=0; iy<Ny; ++iy)
    {
        for (var iz=0; iz<Nz; ++iz)
         {
            for (var ix=0; ix<Nx; ++ix)
            {
                let x = (0.5 + ix)*dL;
                let y = (0.5 + iy)*dL;
                let z = (0.5 + iz)*dL;
                let wsP = [x, y, z];      // world space location of voxel center
                let vsP = [ix, iy, iz]  // voxel space location of lower-left corner (== voxel index)

                let r = Math.sqrt(Math.pow(x-blast_center[0], 2.0) + Math.pow(y-blast_center[1], 2.0)+ Math.pow(z-blast_center[2], 2.0));

                // initialize temperature profile
                let T = this.Tambient * (1.0 + this.blastTemperature*this.smoothSphere(wsP, blast_center, this.blastRadius*lengthScale));

                // initialize velocity to a spherical blast with a finite region:
                V = this.radialFlow(wsP, blast_center, this.blastRadius*lengthScale, this.blastVelocity);

                let frag = this.mapVsToFrag(vsP, Ncol, Nx, Ny, Nz);
                let    b = this.mapFragtoBufferIndex(frag, W, 4);

                // Initialize air variables
                Vair0[b+0]  = V[0];
                Vair0[b+1]  = V[1];
                Vair0[b+2]  = V[2];
                divVair0[b] = 0.0;
                Pair0[b]    = 0.0;
                Tair0[b]    = T;

                // Initialize debris density field
                let density = y < this.debrisHeight*L[1] ? Math.exp(-3.0*this.debrisFalloff*y/(this.debrisHeight*L[1])) : 0.0;
                debris0[b+0] = density;
                debris0[b+1] = density;
                debris0[b+2] = density;
            }
        }
    }

    // Initial velocity texture
    this.Vair = [ new GLU.Texture(W, H, 4, true, true, true, Vair0),
                  new GLU.Texture(W, H, 4, true, true, true, Vair0) ];

    // Initial velocity divergence texture
    this.divVair = new GLU.Texture(W, H, 4, true, true, true, divVair0);

    // Initial pressure texture
    this.Pair = [ new GLU.Texture(W, H, 4, true, true, true, Pair0),
                  new GLU.Texture(W, H, 4, true, true, true, Pair0) ];

    // Initial temperature texture
    this.Tair = [ new GLU.Texture(W, H, 4, true, true, true, Tair0),
                  new GLU.Texture(W, H, 4, true, true, true, Tair0) ];

    // Initial debris texture
    this.debris = [ new GLU.Texture(W, H, 4, true, true, true, debris0),
                    new GLU.Texture(W, H, 4, true, true, true, debris0) ];

    this.frame = 0;

    this.domain = { Nx: Nx,
                    Ny: Ny,
                    Nz: Nz,
                    Ncol: Ncol,
                    Nrow: Nrow,
                    dL: this.dL,
                    L: this.L,
                    W: this.W,
                    H: this.H,
                    boundsMin:    [0.0,       0.0,       0.0      ],
                    boundsMax:    [dL*Nx,     dL*Ny,     dL*Nz    ],
                    boundsCenter: [dL*Nx/2.0, dL*Ny/2.0, dL*Nz/2.0],
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

    gl.viewport(0, 0, this.W, this.H);
    this.quadVbo.bind();

    // Run air/debris force-advection step, writing 0 -> 1
    {
        this.advect_program.bind();
        this.advect_program.uniformI("Nx",             this.domain.Nx);
        this.advect_program.uniformI("Ny",             this.domain.Ny);
        this.advect_program.uniformI("Nz",             this.domain.Nz);
        this.advect_program.uniformI("Ncol",           this.domain.Ncol);
        this.advect_program.uniformI("W",              this.domain.W);
        this.advect_program.uniformI("H",              this.domain.H);
        this.advect_program.uniform3Fv("L",            this.domain.L);
        this.advect_program.uniformF("dL",             this.domain.dL);
        this.advect_program.uniformF("timestep",      this.timestep);
        this.advect_program.uniformF("buoyancy",      this.buoyancy);
        this.advect_program.uniformF("gravity",       this.gravity);
        this.advect_program.uniformF("radiationLoss", this.radiationLoss * 1.0e-2);
        this.advect_program.uniformF("T0",            this.domain.T0);
        this.advect_program.uniformF("Tambient",      this.Tambient);

        this.fbo.bind();
        this.fbo.drawBuffers(4);
        this.fbo.attachTexture(this.Vair[1], 0);   // write to  Vair[1]
        this.fbo.attachTexture(this.Pair[1], 1);   // write to  Pair[1]
        this.fbo.attachTexture(this.Tair[1], 2);   // write to  Tair[1]
        this.fbo.attachTexture(this.debris[1], 3); // write to  debris[1]

        this.Vair[0].bind(0);                    // read from Vair[0]
        this.Pair[0].bind(1);                    // read from Pair[0]
        this.Tair[0].bind(2);                    // read from Tair[0]
        this.debris[0].bind(3);                  // read from debris[0]
        this.advect_program.uniformTexture("Vair_sampler", this.Vair[0]);
        this.advect_program.uniformTexture("Pair_sampler", this.Pair[0]);
        this.advect_program.uniformTexture("Tair_sampler", this.Tair[0]);
        this.advect_program.uniformTexture("debris_sampler", this.debris[0]);
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.detachTexture(3);
        this.fbo.unbind();

        // copy temperature 1 -> 0
        // @todo: get rid of copy and do ping-pong
        this.copy_program.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Tair[0], 0); // write to  Tair[0]
        this.Tair[1].bind(0);                    // read from Tair[1]
        this.copy_program.uniformTexture("Qin", this.Tair[1]);
        this.quadVbo.draw(this.copy_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();

        // copy density 1 -> 0
        // @todo: get rid of copy and do ping-pong
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

    // Compute velocity divergence
    {
        this.div_program.bind();
        this.div_program.uniformI("Nx",            this.domain.Nx);
        this.div_program.uniformI("Ny",            this.domain.Ny);
        this.div_program.uniformI("Nz",            this.domain.Nz);
        this.div_program.uniformI("Ncol",          this.domain.Ncol);
        this.div_program.uniform3Fv("L",           this.domain.L);
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
        this.project_program.uniformI("Nx",            this.domain.Nx);
        this.project_program.uniformI("Ny",            this.domain.Ny);
        this.project_program.uniformI("Nz",            this.domain.Nz);
        this.project_program.uniformI("Ncol",          this.domain.Ncol);
        this.project_program.uniform3Fv("L",           this.domain.L);
        this.project_program.uniformF("dL",            this.domain.dL);
        this.project_program.uniformF("timestep",  this.timestep);
        this.project_program.uniformF("expansion", this.expansion);

        // Update pressure field by Jacobi iteration
        // (using last frame pressure as a warm-start)
        // Sequence:
        //   1  ->  0  (read black from 1, write red   to 0)
        //   0  ->  1  (read red   from 0, write black to 1)
        for (let n=0; n<Math.floor(this.NprojSteps); ++n)
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
        this.update_program.uniformI("Nx",            this.domain.Nx);
        this.update_program.uniformI("Ny",            this.domain.Ny);
        this.update_program.uniformI("Nz",            this.domain.Nz);
        this.update_program.uniformI("Ncol",          this.domain.Ncol);
        this.update_program.uniform3Fv("L",           this.domain.L);
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
