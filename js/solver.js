

/** @constructor
*/
var Solver = function()
{
    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource({
    	'advect'         : {'v': 'advect-vertex-shader',       'f': 'advect-fragment-shader'},
        'project-black'  : {'v': 'project-black-vertex-shader','f': 'project-black-fragment-shader'},
        'project-red'    : {'v': 'project-red-vertex-shader',  'f': 'project-red-fragment-shader'},
        'copy'           : {'v': 'copy-vertex-shader',         'f': 'copy-fragment-shader'},
        'update'         : {'v': 'update-vertex-shader',       'f': 'update-fragment-shader'}
    });

    this.fbo = new GLU.RenderTarget();
    this.quadVbo = this.createQuadVbo();
    this.compileShaders();
}

Solver.prototype.compileShaders = function()
{
	this.advect       = new GLU.Shader('advect',  this.shaderSources, null);
    this.projectblack = new GLU.Shader('project-black', this.shaderSources, null);
    this.projectred   = new GLU.Shader('project-red', this.shaderSources, null);
    this.copy         = new GLU.Shader('copy',    this.shaderSources, null);
    this.update       = new GLU.Shader('update',  this.shaderSources, null);
    this.temperature  = new GLU.Shader('temperature',  this.shaderSources, null);
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

	let Qair0 = new Float32Array(Nr*Ny*4); // air
	let Qpar0 = new Float32Array(Nr*Ny*4); // particulate matter

    this.Delta = 1.0;     // voxel size is 1
    let Delta = this.Delta
	let R = Delta * Nr;  // full domain radius in voxels
	let H = Delta * Ny;  // full domain height in voxels

    let T0 = 1.0;        // nominal ambient temperature
    let p = 0;           // initialize pressure to zero (@todo: warm start?)
    this.g = [0, 0, -1]; // nominal gravitational acceleration
    this.beta = 1.0;     // nominal buoyancy

	for (let iy=0; iy<Ny; ++iy)
	{
		// ground height (domain bottom) is y = 0.0
		let y = (0.5 + iy)*dy;
		for (var ir=0; ir<Nr; ++ir)
    	{
    		let r = (0.5 + ir)*dr;

    		// initial temperature profile
    		let T = T0;

    		// Initialize thermal medium:
    		let index = iy*Nr + ir;

            // Initialize velocity to a spherical blast with a finite region:
            let V = Delta; // blast velocity in voxels / timestep
            let [vr, vy] = radialFlow([r, y], [0.0, this.blastHeight], 2.0*this.blastRadius, V);

 			// @todo: deposit blast energy as an initial condition here, rather than in solver
    		Qair0[4*index+0] = vr; // vr, air, initially zero
    		Qair0[4*index+1] = vy; // vy, air, initially zero
    		Qair0[4*index+2] = T;   // T air
    		Qair0[4*index+3] = p;   // p air

    		// Initialize particulate matter density:
    		Qpar0[4*index+0] = 0.0; // unused
    		Qpar0[4*index+1] = 0.0; // unused
    		Qpar0[4*index+2] = 0.0; // unused
    		Qpar0[4*index+3] = 0.0; // unused
    	}
	}

	this.Qair = [ new GLU.Texture(Nr, Ny, 4, true, false, true, Qair0), new GLU.Texture(Nr, Ny, 4, true, false, true, Qair0) ];
	//this.Qpar = [ new GLU.Texture(Nr, Ny, 4, true, false, true, Qpar0), new GLU.Texture(Nr, Ny, 4, true, false, true, Qpar0) ];

	this.timestep = 0;

	this.domain = { boundsMin: [-R, 0.0, -R],
	                boundsMax: [ R,  H,   R],
	                center:    [0.0, H/2.0, 0.0],
	                radius: R,
	                height: H,
	                Nr: Nr,
	                Ny: Ny }

	// constants
	this.dr = dr;
	this.dy = dy;
	this.dt = dt; 
	this.g = g;
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
Solver.prototype.getQpar = function()
{
	return this.Qpar[0];
}	

Solver.prototype.step = function()
{
	console.log("Running timestep: ", this.timestep);

	let gl = GLU.gl;

    gl.viewport(0, 0, this.Nr, this.Ny);
    this.quadVbo.bind();

	// Run force/advection step, writing 0 -> 1
	{
		this.advect.bind();
		this.advect.uniformI("Nr",       this.Nr);
		this.advect.uniformI("Ny",       this.Ny);
		this.advect.uniformF("Delta",    this.Delta);
		this.advect.uniformF("g",        this.g);
        this.advect.uniformF("beta",     this.beta);

		this.fbo.bind();
		this.fbo.drawBuffers(1);
		this.fbo.attachTexture(this.Qair[1], 0);         // write to  Qair[1]
	    this.Qair[currIndex].bind(0);                    // read from Qair[0]
	    this.advect.uniformTexture("Qair", this.Qair[0]);
	    this.quadVbo.draw(this.advect, gl.TRIANGLE_FAN);
        this.fbo.unbind();
	}

    // Run pressure projection step
	{
        this.project.bind();
		this.project.uniformI("Nr",       this.Nr);
		this.project.uniformI("Ny",       this.Ny);
		this.project.uniformF("Delta",    this.Delta);
        this.fbo.bind();

        // Update pressure field by red-black relaxation
        // (using last frame pressure as a warm-start)
        // Sequence:
        //   1  ->  0  (read black from 1, write red   to 0)
        //   0  ->  1  (read red   from 0, write black to 1)
        const numIterations = 10;
        for (let n=0; n<numIterations; ++n)
        {
            this.fbo.attachTexture(this.Qair[0], 0); // write red to 0
	        this.Qair[1].bind(0); // read black from 1
            this.projectred.uniformTexture("Qair", this.Qair[1]);
            this.quadVbo.draw(this.project, gl.TRIANGLE_FAN);

            this.fbo.attachTexture(this.Qair[1], 0); // write black to 1
	        this.Qair[0].bind(0); // read red from 0
            this.projectblack.uniformTexture("Qair", this.Qair[0]);
            this.quadVbo.draw(this.project, gl.TRIANGLE_FAN);
        }
        this.fbo.unbind();
    }

    // Copy red from 0 -> 1
    {
        this.copy.bind();
        this.fbo.bind();
		this.fbo.drawBuffers(1);
		this.fbo.attachTexture(this.Qair[1], 0); // write to  Qair[1]
	    this.Qair[0].bind(0);                    // read from Qair[0]
        this.Qair[1].bind(0);                    // read from Qair[1]
	    this.copy.uniformTexture("Qred",   this.Qair[0]);
        this.copy.uniformTexture("Qblack", this.Qair[1]);
	    this.quadVbo.draw(this.copy, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Update velocity from 1 -> 0
    {
        this.update.bind();
        this.update.uniformI("Nr",       this.Nr);
		this.update.uniformI("Ny",       this.Ny);
        this.update.uniformF("Delta",    this.Delta);
        this.fbo.bind();
        this.fbo.attachTexture(this.Qair[0], 0);         // write to  Qair[0]
        this.Qair[1].bind(0);                            // read from Qair[1]
        this.update.uniformTexture("Qair", this.Qair[1]);
        this.quadVbo.draw(this.update, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    /*
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
    */

	this.fbo.unbind();
    gl.bindTexture(gl.TEXTURE_2D, null);

    // @todo: advect (and inject) particulate matter passively in the flow, 
    //        to model dust (and fission products)

    this.timestep++;
}
