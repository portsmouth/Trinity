

/** @constructor
*/
var Solver = function()
{
    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource({
        'advect'         : {'v': 'advect-vertex-shader',       'f': 'advect-fragment-shader'       },
        'project'        : {'v': 'project-vertex-shader',      'f': 'project-fragment-shader'      },
        'copy'           : {'v': 'copy-vertex-shader',         'f': 'copy-fragment-shader'         },
        'merge'          : {'v': 'merge-vertex-shader',        'f': 'merge-fragment-shader'        },
        'update'         : {'v': 'update-vertex-shader',       'f': 'update-fragment-shader'       },
        'debris'         : {'v': 'debris-vertex-shader',       'f': 'debris-fragment-shader'       }
    });

    this.fbo = new GLU.RenderTarget();
    this.quadVbo = this.createQuadVbo();
    this.compileShaders();
}

Solver.prototype.compileShaders = function()
{
	this.advect       = new GLU.Shader('advect',        this.shaderSources, null);
    this.project      = new GLU.Shader('project',       this.shaderSources, null);
    this.copy         = new GLU.Shader('copy',          this.shaderSources, null);
    this.merge        = new GLU.Shader('merge',         this.shaderSources, null);
    this.update       = new GLU.Shader('update',        this.shaderSources, null);
    this.debris       = new GLU.Shader('debris',        this.shaderSources, null);
}


Solver.prototype.resize = function(Nr, Ny)
{
    this.Nr = Nr;
    this.Ny = Ny;
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
	let rho= Math.sqrt(Math.pow(x[0]-c[0], 2.0) + Math.pow(x[1]-c[1], 2.0));
	let t = rho/R;
	if (t<=1.0) return Math.max(0.0, 1.0 - t*t*(3.0 - 2.0*t));
	return 0.0;
}

Solver.prototype.radialFlow = function(x, c, R, V)
{
    // x = [r, y]
	let rho = Math.sqrt(Math.pow(x[0]-c[0], 2.0) + Math.pow(x[1]-c[1], 2.0));
	let t = rho/R;
    let v = (t<=1.0) ? V * Math.max(0.0, 1.0 - t*t*(3.0 - 2.0*t)) : 0.0;
    let vr = v * (x[0]-c[0])/rho
    let vy = v * (x[1]-c[1])/rho
    return [vr, vy]
}

Solver.prototype.initialize = function()
{
	let Nr = this.Nr;
	let Ny = this.Ny;

	let Qair0    = new Float32Array(Nr*Ny*4); // air
	let Qdebris0 = new Float32Array(Nr*Ny*4); // particulate matter

    this.Delta = 1.0;     // voxel size is 1
    let Delta = this.Delta
	let R = Delta * Nr;  // full domain radius in voxels
	let H = Delta * Ny;  // full domain height in voxels

    let blastHeight = 0.1*H;
    let blastRadius = 10.0*Delta;
    let debrisHeight = H;

    let T0 = 1.0;        // nominal ambient temperature
    let p = 1.0;         // initialize pressure to zero (@todo: warm start?)
    this.g = 0.01;        // nominal gravitational acceleration
    this.beta = 0.02;     // nominal buoyancy

    for (let iy=0; iy<Ny; ++iy)
    {
        // ground height (domain bottom) is y = 0.0
        let y = (0.5 + iy)*Delta;
        for (var ir=0; ir<Nr; ++ir)
        {
            let r = (0.5 + ir)*Delta;

            // initial temperature profile
            let T = T0 * (1.0 + 500.0*this.smoothSphere([r, y], [0.0, blastHeight], 2.0*blastRadius)) * Math.exp(-y/H);

            // Initialize thermal medium:
            let index = iy*Nr + ir;

            // Initialize velocity to a spherical blast with a finite region:
            let speed = 0.5*Delta; // blast velocity in voxels / timestep
            V = this.radialFlow([r, y], [0.0, blastHeight], 2.0*blastRadius, speed);
            let vr = V[0]
            let vy = V[1]

            Qair0[4*index+0] = vr; // vr, air
            Qair0[4*index+1] = vy; // vy, air
            Qair0[4*index+2] = T;  // T air
            Qair0[4*index+3] = p;  // p air

            // Initialize debris
            Qdebris0[4*index+0] = 10.0 * Math.exp(-30.0*y/debrisHeight) / (2.0 * R); //(y < debrisHeight) ? 1.0 : 0.0; // density
            Qdebris0[4*index+1] = 0.20;  // relative extinction in R channel
            Qdebris0[4*index+2] = 0.75; // relative extinction in G channel
            Qdebris0[4*index+3] = 1.00;  // relative extinction in B channel
        }
    }

	this.Qair = [ new GLU.Texture(Nr, Ny, 4, true, false, true, Qair0),  
                  new GLU.Texture(Nr, Ny, 4, true, false, true, Qair0),
                  new GLU.Texture(Nr, Ny, 4, true, false, true, Qair0) ];

	this.Qdebris = [ new GLU.Texture(Nr, Ny, 4, true, false, true, Qdebris0), 
                     new GLU.Texture(Nr, Ny, 4, true, false, true, Qdebris0) ];

	this.timestep = 0;

	this.domain = { boundsMin: [-R, 0.0, -R],
	                boundsMax: [ R,  H,   R],
	                center:    [0.0, H/2.0, 0.0],
                    Delta:  Delta,
	                radius: R,
	                height: H,
	                Nr: Nr,
	                Ny: Ny }
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

// air simulation
Solver.prototype.getQair = function()
{
	return this.Qair[0];
}	

// particulate matter simulation
Solver.prototype.getQdebris = function()
{
	return this.Qdebris[0];
}	

Solver.prototype.step = function()
{
    console.log("Running timestep: ", this.timestep);

    let gl = GLU.gl;

    gl.viewport(0, 0, this.Nr, this.Ny);
    this.quadVbo.bind();

    // Run air force/advection step, writing 0 -> 1
    {
        this.advect.bind();
        this.advect.uniformI("Nr",       this.Nr);
        this.advect.uniformI("Ny",       this.Ny);
        this.advect.uniformF("Delta",    this.Delta);
        this.advect.uniformF("g",        this.g);
        this.advect.uniformF("beta",     this.beta);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Qair[1], 0); // write to  Qair[1]
        this.Qair[0].bind(0);                    // read from Qair[0]
        this.advect.uniformTexture("Qin", this.Qair[0]);
        this.quadVbo.draw(this.advect, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Run air pressure projection step
    {
        // Update pressure field by red-black relaxation
        // (using last frame pressure as a warm-start)
        // Sequence:
        //   1  ->  0  (read black from 1, write red   to 0)
        //   0  ->  1  (read red   from 0, write black to 1)
        const numIterations = 50;
        for (let n=0; n<numIterations; ++n)
        {
            this.project.bind();
            this.project.uniformI("Nr",       this.Nr);
            this.project.uniformI("Ny",       this.Ny);
            this.project.uniformF("Delta",    this.Delta);
            this.fbo.bind();
            this.fbo.attachTexture(this.Qair[0], 0); // write red to 0
            this.Qair[1].bind(0); // read black from 1
            this.project.uniformTexture("Qin", this.Qair[1]);
            this.quadVbo.draw(this.project, gl.TRIANGLE_FAN);
            this.fbo.unbind();

            this.project.bind();
            this.project.uniformI("Nr",       this.Nr);
            this.project.uniformI("Ny",       this.Ny);
            this.project.uniformF("Delta",    this.Delta);
            this.fbo.bind();
            this.fbo.attachTexture(this.Qair[1], 0); // write black to 1
            this.Qair[0].bind(0); // read red from 0
            this.project.uniformTexture("Qin", this.Qair[0]);
            this.quadVbo.draw(this.project, gl.TRIANGLE_FAN);
            this.fbo.unbind();
        }
    }

    // Copy black air cells from 1 -> 2 (scratch)
    {
        this.copy.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Qair[2], 0); // write to  Qair[2]
        this.Qair[1].bind(0);                    // read from Qair[1]
        this.copy.uniformTexture("Qin",   this.Qair[1]);
        this.quadVbo.draw(this.copy, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Merge red (0) and black (2) air cells into 1
    {
        this.merge.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Qair[1], 0); // write to  Qair[1]
        this.Qair[0].bind(0);                    // read from Qair[0]
        this.Qair[2].bind(1);                    // read from Qair[2]
        this.merge.uniformTexture("Qred",   this.Qair[0]);
        this.merge.uniformTexture("Qblack", this.Qair[2]);
        this.quadVbo.draw(this.merge, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Update air velocity from 1 -> 0
    {
        this.update.bind();
        this.update.uniformI("Nr",       this.Nr);
        this.update.uniformI("Ny",       this.Ny);
        this.update.uniformF("Delta",    this.Delta);
        this.fbo.bind();
        this.fbo.attachTexture(this.Qair[0], 0);         // write to  Qair[0]
        this.Qair[1].bind(0);                            // read from Qair[1]
        this.update.uniformTexture("Qin", this.Qair[1]);
        this.quadVbo.draw(this.update, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Advect debris in air flow
    {
        // write debris 0 -> 1
        this.debris.bind();
        this.debris.uniformI("Nr",       this.Nr);
        this.debris.uniformI("Ny",       this.Ny);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Qdebris[1], 0); // write to  Qdebris[1]
        this.Qdebris[0].bind(0);                    // read from Qdebris[0]
        this.Qair[1].bind(1);                       // read from Qair[1]
        this.debris.uniformTexture("Qdebris", this.Qdebris[0]);
        this.debris.uniformTexture("Qair", this.Qair[1]);
        this.quadVbo.draw(this.debris, gl.TRIANGLE_FAN);
        this.fbo.unbind();

        // copy debris 1 -> 0
        this.copy.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Qdebris[0], 0); // write to  Qdebris[0]
        this.Qdebris[1].bind(0);                    // read from Qdebris[1]
        this.copy.uniformTexture("Qin",   this.Qdebris[1]);
        this.quadVbo.draw(this.copy, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

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

    // @todo: advect (and inject) particulate matter passively in the flow, 
    //        to model dust (and fission products)

    this.timestep++;
}
