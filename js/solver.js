

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
    this.timestep = 1.0;
    this.blastHeight = 0.1;
    this.blastRadius = 0.01;
    this.blastTemperature = 2.0; // initial temperature of fireball relative to ambient
    this.blastVelocity = 5.0;     // outward blast speed in voxels/timestep
    this.debrisHeight = 0.3;      // maximum height of dust layer, as a fraction of domain height 
    this.debrisFalloff = 0.5;     // fall-off exponent within dust layer
    this.T0 = 1.0;                // nominal ambient temperature
    this.buoyancy = 0.33;         // nominal buoyancy
    this.expansion = 7.0;
    this.dustWeight = 0.3;        // controls force of gravity on dust
    this.radiationLoss = 0.08;    // radiation loss rate (per timestep fractional absorption)

    // Initialize solver
    this.resize(400, 1000);
}

Solver.prototype.compileShaders = function()
{
	this.advect       = new GLU.Shader('advect',        this.shaderSources, null);
    this.project      = new GLU.Shader('project',       this.shaderSources, null);
    this.div          = new GLU.Shader('div',           this.shaderSources, null);
    this.copy         = new GLU.Shader('copy',          this.shaderSources, null);
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
    let div0     = new Float32Array(Nr*Ny*4); // divergence

	let R = Nr;  // full domain radius in voxels
	let H = Ny;  // full domain height in voxels

    noise.seed(Math.random());

    for (let iy=0; iy<Ny; ++iy)
    {
        // ground height (domain bottom) is y = 0.0
        let y = (0.5 + iy);
        for (var ir=0; ir<Nr; ++ir)
        {
            let r = (0.5 + ir);

            // initial temperature profile
            let T = this.T0 * (1.0 + this.blastTemperature*this.smoothSphere([r, y], [0.0, this.blastHeight*H], 2.0*this.blastRadius*R));

            // Initialize thermal medium:
            let index = iy*Nr + ir;

            // Initialize velocity to a spherical blast with a finite region:
            V = this.radialFlow([r, y], [0.0, this.blastHeight*H], 2.0*this.blastRadius*R, this.blastVelocity);
            let vr = V[0]
            let vy = V[1]

            Qair0[4*index+0] = vr; // vr, air
            Qair0[4*index+1] = vy; // vy, air
            Qair0[4*index+2] = T;  // T air
            Qair0[4*index+3] = 0.0;

            div0[4*index+0] = 0.0;
            div0[4*index+1] = 0.0;
            div0[4*index+2] = 0.0;
            div0[4*index+3] = 0.0;

            // Initialize debris
            Qdebris0[4*index+0] = y < this.debrisHeight*H ? Math.exp(-10.0*this.debrisFalloff*y/(this.debrisHeight*H)) : 0.0;
            Qdebris0[4*index+1] = 0.6 * (1.0 + 0.1*noise.simplex3(32.0*y/H, 32.0*r/R, 0.0));
            Qdebris0[4*index+2] = 0.7 * (1.0 + 0.1*noise.simplex3(32.0*y/H, 32.0*r/R, 1.0));
            Qdebris0[4*index+3] = 0.8 * (1.0 + 0.1*noise.simplex3(32.0*y/H, 32.0*r/R, 2.0));
        }
    }

	this.Qair = [ new GLU.Texture(Nr, Ny, 4, true, true, true, Qair0),  
                  new GLU.Texture(Nr, Ny, 4, true, true, true, Qair0),
                  new GLU.Texture(Nr, Ny, 4, true, true, true, Qair0) ];

    this.divv = new GLU.Texture(Nr, Ny, 4, true, true, true, div0);

	this.Qdebris = [ new GLU.Texture(Nr, Ny, 4, true, true, true, Qdebris0), 
                     new GLU.Texture(Nr, Ny, 4, true, true, true, Qdebris0) ];

	this.frame = 0;

	this.domain = { boundsMin: [-R, 0.0, -R],
	                boundsMax: [ R,  H,   R],
	                center:    [0.0, H/2.0, 0.0],
	                radius: R,
	                height: H,
	                Nr: Nr,
	                Ny: Ny,
                    T0: this.T0 }
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
    console.log("Running frame: ", this.frame);

    let gl = GLU.gl;

    gl.viewport(0, 0, this.Nr, this.Ny);
    this.quadVbo.bind();

    // Run air force/advection step, writing 0 -> 1
    {
        this.advect.bind();
        this.advect.uniformI("Nr",            this.Nr);
        this.advect.uniformI("Ny",            this.Ny);
        this.advect.uniformF("volRadius",     this.domain.radius);
        this.advect.uniformF("volHeight",     this.domain.height);
        this.advect.uniformF("buoyancy",      this.buoyancy);
        this.advect.uniformF("dustWeight",    this.dustWeight);
        this.advect.uniformF("radiationLoss", this.radiationLoss * 1.0e-2);
        this.advect.uniformF("T0",            this.T0);
        this.advect.uniformF("timestep",      this.timestep);
        
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Qair[1], 0); // write to  Qair[1]
        this.Qair[0].bind(0);                    // read from Qair[0]
        this.advect.uniformTexture("Qin", this.Qair[0]);

        this.Qdebris[0].bind(1);                 // read from Qdebris[0]
        this.advect.uniformTexture("Qdust", this.Qdebris[0]); 

        this.quadVbo.draw(this.advect, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Compute velocity divergence
    {
        this.div.bind();
        this.div.uniformI("Nr",       this.Nr);
        this.div.uniformI("Ny",       this.Ny);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.divv, 0); // write to this.divv
        this.Qair[1].bind(0);                 // read from Qair[1]
        this.div.uniformTexture("Qin", this.Qair[1]);
        this.quadVbo.draw(this.div, gl.TRIANGLE_FAN);
        this.fbo.unbind();
    }

    // Run air pressure projection step
    {
        this.project.bind();
        this.project.uniformI("Nr",       this.Nr);
        this.project.uniformI("Ny",       this.Ny);
        this.project.uniformF("volRadius", this.domain.radius);
        this.project.uniformF("volHeight", this.domain.height);
        this.project.uniformF("expansion", this.expansion);
        this.project.uniformF("timestep", this.timestep);
        
        // Update pressure field by Jacobi iteration
        // (using last frame pressure as a warm-start)
        // Sequence:
        //   1  ->  0  (read black from 1, write red   to 0)
        //   0  ->  1  (read red   from 0, write black to 1)
        const numIterations = 16;
        for (let n=0; n<numIterations; ++n)
        {
            this.fbo.bind();
            this.fbo.attachTexture(this.Qair[0], 0);          // write to 0
            this.Qair[1].bind(0);                             // read from 1
            this.project.uniformTexture("Qin", this.Qair[1]);
            this.divv.bind(1);                               // read divv
            this.project.uniformTexture("div", this.divv);
            this.quadVbo.draw(this.project, gl.TRIANGLE_FAN);
            this.fbo.unbind();

            this.fbo.bind();
            this.fbo.attachTexture(this.Qair[1], 0);          // write to 1
            this.Qair[0].bind(0);                             // read from 0
            this.project.uniformTexture("Qin", this.Qair[0]);
            this.divv.bind(1);                                // read divv
            this.project.uniformTexture("div", this.divv);
            this.quadVbo.draw(this.project, gl.TRIANGLE_FAN);
            this.fbo.unbind();
        }
    }

    // Update air velocity from 1 -> 0
    {
        this.update.bind();
        this.update.uniformI("Nr",       this.Nr);
        this.update.uniformI("Ny",       this.Ny);
        this.update.uniformF("volRadius", this.domain.radius);
        this.update.uniformF("volHeight", this.domain.height);
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
        this.debris.uniformF("timestep",      this.timestep);
        this.debris.uniformF("volRadius", this.domain.radius);
        this.debris.uniformF("volHeight", this.domain.height);
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

    this.frame++;
}
