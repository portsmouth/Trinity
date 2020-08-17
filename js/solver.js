

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
    this.settings.max_timesteps = 100;
    this.settings.expansion = 0.001;

    this.uniforms_float = {};
    this.uniforms_vec3 = {};

    this.compiled_successfully = false;
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

    let common_glsl  = trinity.getGlsl('common') + '\n';
    let collide_glsl    = common_glsl  + trinity.getGlsl('collide') + '\n';
    let initial_glsl    = common_glsl  + trinity.getGlsl('initial');
    let inject_glsl     = collide_glsl + trinity.getGlsl('inject');
    let influence_glsl  = collide_glsl + trinity.getGlsl('influence');

    this.uniforms_float = {};
    this.uniforms_vec3 = {};

    trinity.getGUI().refresh();
    trinity.show_errors();
    this.compiled_successfully = false;

    this.initial_program   = new GLU.Shader('initial',   this.shaderSources,  { _USER_CODE_: initial_glsl });
    if (!this.initial_program.program) { this.initial_program = null; return; }

    this.inject_program    = new GLU.Shader('inject',    this.shaderSources,  { _USER_CODE_: inject_glsl });
    if (!this.inject_program.program) { this.inject_program = null; return; }

    this.advect_program    = new GLU.Shader('advect',    this.shaderSources,  { _USER_CODE_: influence_glsl });
    if (!this.advect_program.program) { this.advect_program = null; return; }

    this.project_program   = new GLU.Shader('project',   this.shaderSources,  { _USER_CODE_: collide_glsl });
    if (!this.project_program.program) { this.project_program = null; return; }

    this.div_program       = new GLU.Shader('div',       this.shaderSources,  { _USER_CODE_: collide_glsl });
    if (!this.div_program.program) { this.div_program = null; return; }

    this.update_program    = new GLU.Shader('update',    this.shaderSources,  { _USER_CODE_: collide_glsl });
    if (!this.update_program.program) { this.update_program = null; return; }

    this.copy_program      = new GLU.Shader('copy',      this.shaderSources, null);
    if (!this.copy_program.program) { this.copy_program = null; return; }

    this.vorticity_program = new GLU.Shader('vorticity', this.shaderSources, null);
    if (!this.vorticity_program.program) { this.vorticity_program = null; return; }

    this.compiled_successfully = true;
    trinity.hide_errors();
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

Solver.prototype.alloc_texture_pair = function(W, H)
{
    return  [ new GLU.Texture(W, H, 4, true, true, true, null),
              new GLU.Texture(W, H, 4, true, true, true, null) ];
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
    let dL = 1.0; // voxel size in world units
    this.dL = dL;
    let L = [dL*Nx, dL*Ny, dL*Nz];
    this.L = L;

    this.Vair       = this.alloc_texture_pair(W, H); // air velocity field
    this.Pair       = this.alloc_texture_pair(W, H); // air pressure field
    this.Tair       = this.alloc_texture_pair(W, H); // air temperature field
    this.absorption = this.alloc_texture_pair(W, H); // medium absorption field
    this.scattering = this.alloc_texture_pair(W, H); // medium scattering field
    this.divVair    = new GLU.Texture(W, H, 4, true, true, true, null); // divergence field
    this.vorticity  = new GLU.Texture(W, H, 4, true, true, true, null); // vorticity field

    // object which describes the grid geometry
    this.domain = {
        Nx: Nx,
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
    };

    this.frame = 0;
    this.time = 0.0;

    trinity.getRenderer().setBounds(this.domain);
    trinity.render_dirty();
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

Solver.prototype.getAbsorption = function()
{
	return this.absorption[0];
}

Solver.prototype.getScattering = function()
{
	return this.scattering[0];
}

Solver.prototype.step = function()
{
    if (!this.compiled_successfully)
        return;
    if (this.paused)
        return;
    if (this.frame >= this.settings.max_timesteps)
        this.restart();

    let gl = GLU.gl;
    gl.viewport(0, 0, this.W, this.H);
    this.quadVbo.bind();

    let BUFFER_0 = this.frame % 2; // ping-pong the index 0 buffer for v, P
    let BUFFER_1 = 1 - BUFFER_0;   // ping-pong the index 1 buffer for v, P

    // Run initial conditions shader on first frame
    // -> output(v0, P0, T0, M0) (where "M" means the medium fields absorption/scattering)
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
        this.fbo.drawBuffers(5);
        this.fbo.attachTexture(this.Vair[BUFFER_0], 0); // write to  Vair[BUFFER_0]
        this.fbo.attachTexture(this.Pair[BUFFER_0], 1); // write to  Pair[BUFFER_0]
        this.fbo.attachTexture(this.Tair[0], 2);        // write to  Tair[0]
        this.fbo.attachTexture(this.absorption[0], 3);  // write to absorption[0]
        this.fbo.attachTexture(this.scattering[0], 4);  // write to scattering[0]
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.detachTexture(3);
        this.fbo.detachTexture(4);
        this.fbo.unbind();
    }

    // Inject fresh mass, heat and medium:
    // input(v0, T0, M0) -> output(v1, T1, M1)
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
        this.inject_program.uniformF("timestep",       this.settings.timestep);
        this.syncUserUniforms(this.inject_program);

        this.fbo.bind();
        this.fbo.drawBuffers(4);
        this.fbo.attachTexture(this.Vair[BUFFER_1], 0); // write to  Vair[BUFFER_1]
        this.fbo.attachTexture(this.Tair[1], 1);        // write to  Tair[1]
        this.fbo.attachTexture(this.absorption[1], 2);  // write to  absorption[1]
        this.fbo.attachTexture(this.scattering[1], 3);  // write to  scattering[1]
        this.Vair[BUFFER_0].bind(0);                    // read from Vair[BUFFER_0]
        this.Tair[0].bind(1);                           // read from Tair[0]
        this.absorption[0].bind(2);                     // read from absorption[0]
        this.scattering[0].bind(3);                     // read from scattering[0]
        this.inject_program.uniformTexture("Vair_sampler", this.Vair[BUFFER_0]);
        this.inject_program.uniformTexture("Tair_sampler", this.Tair[0]);
        this.inject_program.uniformTexture("absorption_sampler", this.absorption[0]);
        this.inject_program.uniformTexture("scattering_sampler", this.scattering[0]);
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.detachTexture(3);
        this.fbo.unbind();
    }

    // Compute vorticity
    // input(v1) -> output(vorticity)
    if (this.settings.vorticity_scale > 0.0)
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
        this.fbo.attachTexture(this.vorticity, 0);    // write to vorticity
        this.Vair[BUFFER_1].bind(0);                  // read from Vair[BUFFER_1]
        this.vorticity_program.uniformTexture("Vair_sampler", this.Vair[BUFFER_1]);
        this.quadVbo.draw(this.update_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

    // Run air/debris force-advection step, writing 1 -> 0
    // input(v1, P0, T1, M1, vorticity) -> output(v0, P1, T0, M0)
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
        this.syncUserUniforms(this.advect_program);

        this.fbo.bind();
        this.fbo.drawBuffers(5);
        this.fbo.attachTexture(this.Vair[BUFFER_0], 0); // write to  Vair[BUFFER_0]
        this.fbo.attachTexture(this.Pair[BUFFER_1], 1); // write to  Pair[BUFFER_1]
        this.fbo.attachTexture(this.Tair[0], 2);        // write to  Tair[0]
        this.fbo.attachTexture(this.absorption[0], 3);  // write to  absorption[0]
        this.fbo.attachTexture(this.scattering[0], 4);  // write to  scattering[0]
        this.Vair[BUFFER_1].bind(0);                    // read from Vair[BUFFER_1]
        this.Pair[BUFFER_0].bind(1);                    // read from Pair[BUFFER_0]
        this.Tair[1].bind(2);                           // read from Tair[1]
        this.absorption[1].bind(3);                     // read from absorption[1]
        this.scattering[1].bind(4);                     // read from scattering[1]
        this.vorticity.bind(5);                         // read from vorticity
        this.advect_program.uniformTexture("Vair_sampler", this.Vair[BUFFER_1]);
        this.advect_program.uniformTexture("Pair_sampler", this.Pair[BUFFER_0]);
        this.advect_program.uniformTexture("Tair_sampler", this.Tair[1]);
        this.advect_program.uniformTexture("absorption_sampler", this.absorption[1]);
        this.advect_program.uniformTexture("scattering_sampler", this.scattering[1]);
        this.advect_program.uniformTexture("vorticity_sampler", this.vorticity);
        this.quadVbo.draw(this.advect_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.detachTexture(1);
        this.fbo.detachTexture(2);
        this.fbo.detachTexture(3);
        this.fbo.detachTexture(4);
        this.fbo.unbind();
    }

    gl.bindTexture(gl.TEXTURE_2D, null);

    // Compute velocity divergence
    // input(v0) -> output(divergence)
    {
        this.div_program.bind();
        this.div_program.uniformI("Nx",            this.domain.Nx);
        this.div_program.uniformI("Ny",            this.domain.Ny);
        this.div_program.uniformI("Nz",            this.domain.Nz);
        this.div_program.uniformI("Ncol",          this.domain.Ncol);
        this.div_program.uniform3Fv("L",           this.domain.L);
        this.div_program.uniformF("dL",            this.domain.dL);
        this.div_program.uniformF("time",          this.time);
        this.syncUserUniforms(this.div_program);
        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.divVair, 0); // write to divVair
        this.Vair[BUFFER_0].bind(0);             // read from Vair[BUFFER_0]
        this.div_program.uniformTexture("Vair_sampler", this.Vair[BUFFER_0]);
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
        this.project_program.uniformF("time",          this.time);
        this.project_program.uniformF("timestep",  this.settings.timestep);
        this.project_program.uniformF("expansion", 0.01*this.settings.expansion);
        this.syncUserUniforms(this.project_program);

        // Update pressure field by Jacobi iteration
        // (using last frame pressure as a warm-start)
        for (let n=0; n<Math.floor(this.settings.NprojSteps); ++n)
        {
            // input(P1, T0, divergence) -> output(P0)
            this.fbo.bind();
            this.fbo.drawBuffers(1);
            this.fbo.attachTexture(this.Pair[BUFFER_0], 0);   // write to Pair[BUFFER_0]
            this.Pair[BUFFER_1].bind(0);                      // read from Pair[BUFFER_1]
            this.Tair[0].bind(1);                             // read from Tair[0]
            this.divVair.bind(2);                             // read divVair
            this.project_program.uniformTexture("Pair_sampler", this.Pair[BUFFER_1]);
            this.project_program.uniformTexture("Tair_sampler", this.Tair[0]);
            this.project_program.uniformTexture("divVair_sampler", this.divVair);
            this.quadVbo.draw(this.project_program, gl.TRIANGLE_FAN);
            this.fbo.detachTexture(0);
            this.fbo.unbind();
            gl.bindTexture(gl.TEXTURE_2D, null);

            // input(P0, T0, divergence) -> output(P1) [end of frame P output -> ping-pong P]
            this.fbo.bind();
            this.fbo.drawBuffers(1);
            this.fbo.attachTexture(this.Pair[BUFFER_1], 0);   // write to Pair[BUFFER_1]
            this.Pair[BUFFER_0].bind(0);                      // read from Pair[BUFFER_0]
            this.Tair[0].bind(1);                             // read from Tair[0]
            this.divVair.bind(2);                             // read divVair
            this.project_program.uniformTexture("Pair_sampler", this.Pair[BUFFER_0]);
            this.project_program.uniformTexture("Tair_sampler", this.Tair[0]);
            this.project_program.uniformTexture("divVair_sampler", this.divVair);
            this.quadVbo.draw(this.project_program, gl.TRIANGLE_FAN);
            this.fbo.detachTexture(0);
            this.fbo.unbind();
            gl.bindTexture(gl.TEXTURE_2D, null);
        }
    }

    // Update air velocity from 1 -> 0
    // input(v0, P1) -> output(v1) [end of frame v output -> ping-pong v]
    {
        this.update_program.bind();
        this.update_program.uniformI("Nx",            this.domain.Nx);
        this.update_program.uniformI("Ny",            this.domain.Ny);
        this.update_program.uniformI("Nz",            this.domain.Nz);
        this.update_program.uniformI("Ncol",          this.domain.Ncol);
        this.update_program.uniform3Fv("L",           this.domain.L);
        this.update_program.uniformF("dL",            this.domain.dL);
        this.update_program.uniformF("timestep", this.timestep);
        this.syncUserUniforms(this.update_program);

        this.fbo.bind();
        this.fbo.drawBuffers(1);
        this.fbo.attachTexture(this.Vair[BUFFER_1], 0); // write to  Vair[BUFFER_1]
        this.Vair[BUFFER_0].bind(0);                    // read from Vair[BUFFER_0]
        this.Pair[BUFFER_1].bind(1);                    // read from Pair[BUFFER_1]
        this.update_program.uniformTexture("Vair_sampler", this.Vair[BUFFER_0]);
        this.update_program.uniformTexture("Pair_sampler", this.Pair[BUFFER_1]);
        this.quadVbo.draw(this.update_program, gl.TRIANGLE_FAN);
        this.fbo.detachTexture(0);
        this.fbo.unbind();
    }

    gl.bindTexture(gl.TEXTURE_2D, null);

    this.frame++;
    this.time += this.settings.timestep;
}
