

/** @constructor
*/
var Solver = function()
{
    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource(["initial",
                                                  "inject",
                                                  "advect",
                                                  "div",
                                                  "project",
                                                  "copy",
                                                  "update",
                                                  "debris",
                                                  "vorticity"
                                                ]);

    this.fbo = new GLU.RenderTarget();
    this.quadVbo = this.createQuadVbo();

    // Default user-adjustable properties
    this.settings = {};

    this.settings.timestep = 1.0;
    this.settings.NprojSteps = 16;
    this.settings.vorticity_scale = 0.2;
    this.settings.Nx = 0;
    this.settings.Ny = 0;
    this.settings.Nz = 0;

    // @todo: all these will be defined in user code:
    this.blastHeight = 0.1;
    this.blastRadius = 0.15;
    this.blastTemperature = 500.0;   // initial temperature of fireball relative to ambient
    this.blastVelocity = 100.0;      // outward blast speed in voxels/timestep
    this.debrisHeight = 0.1;         // maximum height of dust layer, as a fraction of domain height
    this.debrisFalloff = 0.5;        // fall-off exponent within dust layer
    this.T0       = 266.0;           // nominal reference temperature for buoyancy
    this.Tambient = 270.0;           // nominal reference temperature for buoyancy
    this.buoyancy = 0.2;             // initial buoyancy (thermal expansion coeff. of air)
    this.expansion = 0.0;
    this.gravity = 0.003;
    this.radiationLoss = 0.0;        // radiation loss rate (per timestep fractional absorption)

    this.uniforms_float = {};
    this.uniforms_vec3 = {};

    this.compiled_successfully = false;
    this.compileShaders();

    this.resize(128, 512, 128);
    this.restart();
    this.paused = false;
}

Solver.prototype.syncFloatToShader = function(name, value)
{
    this.uniforms_float[name] = value;
}

Solver.prototype.syncColorToShader = function(name, color)
{
    this.uniforms_vec3[name] = color;
}

Solver.prototype.syncUserUniforms = function(program)
{
    for (const key of Object.keys(this.uniforms_float))
    {
        let uniform_name = key;
        let float_value = this.uniforms_float[key];
        program.uniformF(uniform_name, float_value);
    }
    for (const key of Object.keys(this.uniforms_vec3))
    {
        let uniform_name = key;
        let vec3_value = this.uniforms_vec3[key];
        program.uniform3Fv(uniform_name, vec3_value);
    }
}

Solver.prototype.compileShaders = function()
{
    console.warn("[Trinity] Solver.prototype.compileShaders");

    let glsl = trinity.getGlsl();
    if (glsl === undefined)
        return;

    this.uniforms_float = {};
    this.uniforms_vec3 = {};

    trinity.getGUI().refresh();

    replacements = {};
    replacements._USER_CODE_ = glsl;

    this.compiled_successfully = false;

    this.initial_program   = new GLU.Shader('initial',    this.shaderSources, replacements);
    if (!this.initial_program.program) { this.initial_program = null; return; }

    this.inject_program       = new GLU.Shader('inject',        this.shaderSources, replacements);
    if (!this.inject_program.program) { this.inject_program = null; return; }

    this.advect_program       = new GLU.Shader('advect',        this.shaderSources, replacements);
    if (!this.advect_program.program) { this.advect_program = null; return; }

    this.project_program      = new GLU.Shader('project',       this.shaderSources, replacements);
    if (!this.project_program.program) { this.project_program = null; return; }

    this.div_program          = new GLU.Shader('div',           this.shaderSources, replacements);
    if (!this.div_program.program) { this.div_program = null; return; }

    this.update_program       = new GLU.Shader('update',        this.shaderSources, replacements);
    if (!this.update_program.program) { this.update_program = null; return; }

    this.copy_program         = new GLU.Shader('copy',          this.shaderSources, null);
    if (!this.copy_program.program) { this.copy_program = null; return; }

    this.vorticity_program    = new GLU.Shader('vorticity',     this.shaderSources, null);
    if (!this.vorticity_program.program) { this.vorticity_program = null; return; }

    this.compiled_successfully = true;
    trinity.hide_errors();
}

Solver.prototype.updateShaders = function()
{
    console.warn("[Trinity] Solver.prototype.reset");

    this.compileShaders();
    if (!this.compiled_successfully)
        return;

    this.frame = 0;
    this.time = 0.0;
}

Solver.prototype.getDomain = function()
{
	return this.domain;
}

Solver.prototype.pauseToggle = function()
{
    this.paused = !this.paused;
}

Solver.prototype.restart = function()
{
    this.frame = 0;
    this.time = 0.0;
}

Solver.prototype.resize = function(Nx, Ny, Nz)
{
    if ((Nx == this.settings.Nx) &&
        (Ny == this.settings.Ny) &&
        (Nz == this.settings.Nz))
        return;

    this.settings.Nx = Nx;
    this.settings.Ny = Ny;
    this.settings.Nz = Nz;

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
    let vorticity = new Float32Array(W*H*4);  // vorticity field

    let dL = 1.0; // voxel size in world units
    this.dL = dL;

    let L = [dL*Nx, dL*Ny, dL*Nz];
    this.L = L;
    let lengthScale = Math.min(L[0], L[1], L[2]);
    let blast_center = [0.5*L[0], L[1]*this.blastHeight, 0.5*L[2]];

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

    this.vorticity = new GLU.Texture(W, H, 4, true, true, true, vorticity);

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

    this.frame = 0;
    this.time = 0.0;

    trinity.getRenderer().setBounds(this.domain);
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
    //console.warn("[Trinity] Solver.prototype.step");

    if (!this.compiled_successfully)
        return;

    if (this.paused)
        return;

    let gl = GLU.gl;

    gl.viewport(0, 0, this.W, this.H);
    this.quadVbo.bind();

    // Run initial conditions shader on first frame
    if (this.frame==0)
    {
        this.initial_program.bind();
        this.initial_program.uniformI("Nx",             this.domain.Nx);
        this.initial_program.uniformI("Ny",             this.domain.Ny);
        this.initial_program.uniformI("Nz",             this.domain.Nz);
        this.initial_program.uniformI("Ncol",           this.domain.Ncol);
        this.initial_program.uniformI("W",              this.domain.W);
        this.initial_program.uniformI("H",              this.domain.H);
        this.initial_program.uniform3Fv("L",            this.domain.L);
        this.initial_program.uniformF("dL",             this.domain.dL);
        this.syncUserUniforms(this.initial_program);

        this.fbo.bind();
        this.fbo.drawBuffers(4);
        this.fbo.attachTexture(this.Vair[0], 0);   // write to  Vair[0]
        this.fbo.attachTexture(this.Pair[0], 1);   // write to  Pair[0]
        this.fbo.attachTexture(this.Tair[0], 2);   // write to  Tair[0]
        this.fbo.attachTexture(this.debris[0], 3); // write to  debris[0]
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.detachTexture(3);
        this.fbo.unbind();
    }

    // Inject fresh mass and heat:
    {
        this.inject_program.bind();
        this.inject_program.uniformI("Nx",             this.domain.Nx);
        this.inject_program.uniformI("Ny",             this.domain.Ny);
        this.inject_program.uniformI("Nz",             this.domain.Nz);
        this.inject_program.uniformI("Ncol",           this.domain.Ncol);
        this.inject_program.uniformI("W",              this.domain.W);
        this.inject_program.uniformI("H",              this.domain.H);
        this.inject_program.uniform3Fv("L",            this.domain.L);
        this.inject_program.uniformF("dL",             this.domain.dL);
        this.inject_program.uniformF("time",           this.time);
        this.syncUserUniforms(this.inject_program);

        this.fbo.bind();
        this.fbo.drawBuffers(3);
        this.fbo.attachTexture(this.Vair[1], 0);   // write to  Vair[1]
        this.fbo.attachTexture(this.Tair[1], 1);   // write to  Tair[1]
        this.fbo.attachTexture(this.debris[1], 2); // write to  debris[1]
        this.Vair[0].bind(0);                    // read from Vair[0]
        this.Tair[0].bind(1);                    // read from Tair[0]
        this.debris[0].bind(2);                  // read from debris[0]
        this.inject_program.uniformTexture("Vair_sampler", this.Vair[0]);
        this.inject_program.uniformTexture("Tair_sampler", this.Tair[0]);
        this.inject_program.uniformTexture("debris_sampler", this.debris[0]);
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.unbind();

        // copy velocity 1 -> 0
        // @todo: get rid of copy and do ping-pong
        this.copy_program.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Vair[0], 0); // write to  Vair[0]
        this.Vair[1].bind(0);                    // read from Vair[1]
        this.copy_program.uniformTexture("Qin", this.Vair[1]);
        this.quadVbo.draw(this.copy_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
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

        // copy debris 1 -> 0
        // @todo: get rid of copy and do ping-pong
        this.copy_program.bind();
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.debris[0], 0); // write to  debris[0]
        this.debris[1].bind(0);                    // read from debris[1]
        this.copy_program.uniformTexture("Qin", this.Vair[1]);
        this.quadVbo.draw(this.copy_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

    // Compute vorticity
    // @todo (only if vorticity confinement enabled)
    {
        this.vorticity_program.bind();
        this.vorticity_program.uniformI("Nx",             this.domain.Nx);
        this.vorticity_program.uniformI("Ny",             this.domain.Ny);
        this.vorticity_program.uniformI("Nz",             this.domain.Nz);
        this.vorticity_program.uniformI("Ncol",           this.domain.Ncol);
        this.vorticity_program.uniformI("W",              this.domain.W);
        this.vorticity_program.uniformI("H",              this.domain.H);
        this.vorticity_program.uniform3Fv("L",            this.domain.L);
        this.vorticity_program.uniformF("dL",             this.domain.dL);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.vorticity, 0);        // write to vorticity
        this.Vair[0].bind(0);                             // read from Vair[0]
        this.vorticity_program.uniformTexture("Vair_sampler", this.Vair[0]);
        this.quadVbo.draw(this.update_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

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
        this.advect_program.uniformF("time",           this.time);
        this.advect_program.uniformF("timestep",        this.settings.timestep);
        this.advect_program.uniformF("vorticity_scale", this.settings.vorticity_scale);
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
        this.vorticity.bind(4);                  // read from vorticity
        this.advect_program.uniformTexture("Vair_sampler", this.Vair[0]);
        this.advect_program.uniformTexture("Pair_sampler", this.Pair[0]);
        this.advect_program.uniformTexture("Tair_sampler", this.Tair[0]);
        this.advect_program.uniformTexture("debris_sampler", this.debris[0]);
        this.advect_program.uniformTexture("vorticity_sampler", this.vorticity);
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
        this.project_program.uniformF("timestep",  this.settings.timestep);
        this.project_program.uniformF("expansion", this.expansion);

        // Update pressure field by Jacobi iteration
        // (using last frame pressure as a warm-start)
        // Sequence:
        //   1  ->  0  (read black from 1, write red   to 0)
        //   0  ->  1  (read red   from 0, write black to 1)
        for (let n=0; n<Math.floor(this.settings.NprojSteps); ++n)
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
    this.time += this.settings.timestep;
}
