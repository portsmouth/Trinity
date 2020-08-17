
/** @constructor
*/
var Renderer = function()
{
	let gl = GLU.gl;

    // Default user-adjustable properties
    this.settings = {};

    // Raymarch resolution
    this.settings.Nraymarch = 100;
    this.settings.spp_per_frame = 1;
    this.settings.max_spp = 32;
    this.spp = 0;

    // Bounds
    this.settings.show_bounds = true;

    // Lighting
    this.settings.skyColor = [0.5, 0.5, 1.0];
    this.settings.sunColor = [1.0, 1.0, 0.5];
    this.settings.sunPower = 1.0;
    this.settings.sunLatitude   = 60.0;
    this.settings.sunLongitude  = 0.0;
    this.settings.colliderSpec = [0.5, 0.5, 0.5];
    this.settings.colliderDiffuse = [0.8, 0.2, 0.2];
    this.settings.colliderRoughness = 0.2;

    // Tonemapping
    this.settings.exposure = -1.0;
    this.settings.gamma = 2.2;
    this.settings.saturation = 1.0;

    // Phase function
    this.settings.anisotropy = 0.0;

    // Tweak factors
    this.settings.extinctionScale = 0.0;
    this.settings.emissionScale = 3.0;

    // Internal buffers and programs
    this.boxVbo = null;
    this.quadVbo = this.createQuadVbo();

    this.volumeProgram    = null;
    this.lineProgram      = null;
    this.tonemapProgram   = null;

    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource(["volume", "tonemap", "line"]);

	var render_canvas = trinity.render_canvas;
    render_canvas.width  = window.innerWidth;
    render_canvas.height = window.innerHeight;
    this._width = render_canvas.width;
    this._height = render_canvas.height;

    this.compiled_successfully = false;

     // Trigger initial buffer generation
    this.resize(this._width, this._height);
}


Renderer.prototype.createQuadVbo = function()
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

Renderer.prototype.createBoxVbo = function(origin, extents)
{
	let gl = GLU.gl;
	let o = origin;
	let e = extents;

	var corners = [
	  o[0],        o[1],        o[2],
	  o[0] + e[0], o[1],        o[2],
	  o[0] + e[0], o[1],        o[2],
	  o[0] + e[0], o[1] + e[1], o[2],
	  o[0] + e[0], o[1] + e[1], o[2],
	  o[0],        o[1] + e[1], o[2],
	  o[0],        o[1] + e[1], o[2],
	  o[0],        o[1],        o[2],
	  o[0],        o[1],        o[2] + e[2],
	  o[0] + e[0], o[1],        o[2] + e[2],
	  o[0] + e[0], o[1],        o[2] + e[2],
	  o[0] + e[0], o[1] + e[1], o[2] + e[2],
	  o[0] + e[0], o[1] + e[1], o[2] + e[2],
	  o[0],        o[1] + e[1], o[2] + e[2],
	  o[0],        o[1] + e[1], o[2] + e[2],
	  o[0],        o[1],        o[2] + e[2],
	  o[0],        o[1],        o[2],
	  o[0],        o[1],        o[2] + e[2],
	  o[0] + e[0], o[1],        o[2],
	  o[0] + e[0], o[1],        o[2] + e[2],
	  o[0],        o[1] + e[1], o[2],
	  o[0],        o[1] + e[1], o[2] + e[2],
	  o[0] + e[0], o[1] + e[1], o[2],
	  o[0] + e[0], o[1] + e[1], o[2] + e[2],
	];

	var vbo = new GLU.VertexBuffer();
    vbo.addAttribute("Position", 3, gl.FLOAT, false);
    vbo.init(24);
    vbo.copy(new Float32Array(corners));

    return vbo;
}



Renderer.prototype.setBounds = function(domain)
{
    let o = domain.boundsMin;
    let c = domain.boundsMax;
    let extents = [c[0]-o[0], c[1]-o[1], c[2]-o[2]];
    this.boxVbo = this.createBoxVbo(o, extents)
}

Renderer.prototype.compileShaders = function()
{
    console.warn("[Trinity] Renderer.prototype.compileShaders");

    let common_glsl  = trinity.getGlsl('common') + '\n';
    let collide_glsl = common_glsl  + trinity.getGlsl('collide')  + '\n';
    let render_glsl = collide_glsl + trinity.getGlsl('render');

    trinity.show_errors();
    this.compiled_successfully = false;

    this.volumeProgram  = new GLU.Shader('volume',  this.shaderSources, { _USER_CODE_: render_glsl });
    if (!this.volumeProgram.program) { this.volumeProgram = null; return; }

    this.lineProgram    = new GLU.Shader('line',    this.shaderSources, null);
    if (!this.lineProgram.program) { this.lineProgram = null; return; }

    this.tonemapProgram = new GLU.Shader('tonemap', this.shaderSources, null);
    if (!this.tonemapProgram.program) { this.tonemapProgram = null; return; }

    this.compiled_successfully = true;
    trinity.hide_errors();
}

Renderer.prototype.render = function(solver)
{
    if (!this.compiled_successfully)
        return;

    if (this.spp >= this.settings.max_spp)
        return;

    let gl = GLU.gl;

    gl.disable(gl.DEPTH_TEST);
    gl.disable(gl.BLEND);
    gl.viewport(0.0, 0.0, this._width, this._height);

	var camera = trinity.getCamera();
    var camPos = camera.position.clone();
    var camDir = camera.getWorldDirection();
    var camUp = camera.up.clone();
    camUp.transformDirection( camera.matrixWorld );
    var camX = new THREE.Vector3();
    camX.crossVectors(camUp, camDir);

    let domain = solver.getDomain();
    let Tair   = solver.getTair();           // (the air temperature texture)
    let absorption = solver.getAbsorption(); // (the absorption texture)
    let scattering = solver.getScattering(); // (the scattering texture)

    let LAST_WRITE_BUFF = 0;

    // Volume render
    {
        this.volumeProgram.bind();

        // Bind textures for air temperature and medium absorption/scattering
        absorption.bind(2);
        scattering.bind(3);
        Tair.bind(4);
        this.volumeProgram.uniformTexture("absorption_sampler", absorption);
        this.volumeProgram.uniformTexture("scattering_sampler", scattering);
        this.volumeProgram.uniformTexture("Tair_sampler", Tair);

        // sync UI-bound user-code uniforms to the volume program
        trinity.getSolver().syncUserUniforms(this.volumeProgram);

        // Camera
        this.volumeProgram.uniform2Fv("resolution", [this._width, this._height]);
        this.volumeProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
        this.volumeProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
        this.volumeProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
        this.volumeProgram.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
        this.volumeProgram.uniformF("camFovy", camera.fov);
        this.volumeProgram.uniformF("camAspect", camera.aspect);
        this.volumeProgram.uniform3Fv("volMin", domain.boundsMin);
        this.volumeProgram.uniform3Fv("volMax", domain.boundsMax);
        this.volumeProgram.uniform3Fv("volCenter", domain.boundsCenter);

        // Geometry
        this.volumeProgram.uniformI("Nx",             domain.Nx);
        this.volumeProgram.uniformI("Ny",             domain.Ny);
        this.volumeProgram.uniformI("Nz",             domain.Nz);
        this.volumeProgram.uniformI("Ncol",           domain.Ncol);
        this.volumeProgram.uniformI("W",              domain.W);
        this.volumeProgram.uniformI("H",              domain.H);
        this.volumeProgram.uniform3Fv("L",            domain.L);
        this.volumeProgram.uniformF("dL",             domain.dL);
        this.volumeProgram.uniformI("Nraymarch", this.settings.Nraymarch);

        // Physics
        this.volumeProgram.uniformF("time",           solver.time);
        this.volumeProgram.uniformF("anisotropy", this.settings.anisotropy);
        this.volumeProgram.uniformF("extinctionScale", Math.pow(10.0, this.settings.extinctionScale));
        this.volumeProgram.uniformF("emissionScale", Math.pow(10.0, this.settings.emissionScale));

        let latTheta = (90.0-this.settings.sunLatitude) * Math.PI/180.0;
        let lonPhi = this.settings.sunLongitude * Math.PI/180.0;
        let costheta = Math.cos(latTheta);
        let sintheta = Math.sin(latTheta);
        let cosphi = Math.cos(lonPhi);
        let sinphi = Math.sin(lonPhi);
        let x = sintheta * cosphi;
        let z = sintheta * sinphi;
        let y = costheta;
        let sunDir = [x, y, z];

        this.volumeProgram.uniform3Fv("skyColor", this.settings.skyColor);
        this.volumeProgram.uniform3Fv("sunColor", this.settings.sunColor);
        this.volumeProgram.uniformF("sunPower", Math.pow(10.0, this.settings.sunPower));
        this.volumeProgram.uniform3Fv("sunDir", sunDir);

        this.volumeProgram.uniform3Fv("colliderDiffuse", this.settings.colliderDiffuse);
        this.volumeProgram.uniform3Fv("colliderSpec", this.settings.colliderSpec);
        this.volumeProgram.uniformF("colliderRoughness", this.settings.colliderRoughness);

        // Render into progressive radiance buffer spp_per_frames times:
        for (let n=0; n<this.settings.spp_per_frame; ++n)
        {
            this.fbo.bind();
            this.fbo.drawBuffers(2);
            let READ_BUFF  = this.spp % 2;
            let WRITE_BUFF = 1 - READ_BUFF;
            LAST_WRITE_BUFF = WRITE_BUFF;
            this.progressiveStates[READ_BUFF].bind(this.volumeProgram);
            this.progressiveStates[WRITE_BUFF].attach(this.fbo);
            this.volumeProgram.uniformI("spp", this.spp);
            this.quadVbo.bind();
            this.quadVbo.draw(this.volumeProgram, gl.TRIANGLE_FAN);
            this.spp++;
        }

        if (this.spp > this.settings.max_spp)
            this.spp = this.settings.max_spp;

        this.fbo.unbind();
        gl.bindTexture(gl.TEXTURE_2D, null);
    }

    // Tonemapping / compositing
    {
        this.tonemapProgram.bind();
        let radianceTexCurrent = this.progressiveStates[LAST_WRITE_BUFF].radianceTex;
        radianceTexCurrent.bind(0);
        this.tonemapProgram.uniformTexture("Radiance", radianceTexCurrent);
        this.tonemapProgram.uniformF("exposure", this.settings.exposure);
        this.tonemapProgram.uniformF("invGamma", 1.0/this.settings.gamma);
        this.tonemapProgram.uniformF("saturation", this.settings.saturation);

        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
        gl.blendEquation(gl.FUNC_ADD);
        this.quadVbo.bind();
        this.quadVbo.draw(this.tonemapProgram, gl.TRIANGLE_FAN);
        gl.disable(gl.BLEND);
    }

    // draw simulation bounds
    if (this.settings.show_bounds)
    {
        gl.enable(gl.BLEND);
        gl.disable(gl.DEPTH_TEST);
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

        this.lineProgram.bind();

        // Setup projection matrix
        var projectionMatrix = camera.projectionMatrix.toArray();
        var projectionMatrixLocation = this.lineProgram.getUniformLocation("u_projectionMatrix");
        gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

        this.lineProgram.uniform4Fv("color", [0.5, 0.5, 0.6, 0.15]);

        // Setup modelview matrix (to match camera)
        camera.updateMatrixWorld();
        var matrixWorldInverse = new THREE.Matrix4();
        matrixWorldInverse.getInverse( camera.matrixWorld );
        var modelViewMatrix = matrixWorldInverse.toArray();
        var modelViewMatrixLocation = this.lineProgram.getUniformLocation("u_modelViewMatrix");
        gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

        this.boxVbo.bind();
        this.boxVbo.draw(this.lineProgram, gl.LINES);
        gl.disable(gl.BLEND);
    }

    gl.finish();
}



Renderer.prototype.resize = function(width, height)
{
    this._width = width;
    this._height = height;

    if (this.fbo != undefined)
    this.fbo.unbind();
    this.fbo = new GLU.RenderTarget();

    // Progressive rendering radiance/rng states, for read/write ping-pong
    this.progressiveStates = [new ProgressiveState(this._width, this._height),
                              new ProgressiveState(this._width, this._height)];

    this.quadVbo = this.createQuadVbo();
    this.reset();
}


Renderer.prototype.reset = function()
{
    this.spp = 0;
}


//////////////////////////////////////////////////
// ProgressiveState
//////////////////////////////////////////////////

var ProgressiveState = function(width, height)
{
    var radianceData = new Float32Array(width*height*4); // Path radiance, and sample count
    var rngData      = new Float32Array(width*height*4); // Random number seed
    for (var i=0; i<width*height; ++i)
    {
        rngData[i*4 + 0] = Math.random()*4194167.0;
        rngData[i*4 + 1] = Math.random()*4194167.0;
        rngData[i*4 + 2] = Math.random()*4194167.0;
        rngData[i*4 + 3] = Math.random()*4194167.0;
    }
    this.radianceTex = new GLU.Texture(width, height, 4, true, false, true, radianceData);
    this.rngTex      = new GLU.Texture(width, height, 4, true, false, true, rngData);
}

ProgressiveState.prototype.bind = function(shader)
{
    this.radianceTex.bind(0);
    this.rngTex.bind(1);
    shader.uniformTexture("Radiance", this.radianceTex);
    shader.uniformTexture("RngData", this.rngTex);
}

ProgressiveState.prototype.attach = function(fbo)
{
    var gl = GLU.gl;
    fbo.attachTexture(this.radianceTex, 0);
    fbo.attachTexture(this.rngTex, 1);
    /*
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) != gl.FRAMEBUFFER_COMPLETE)
    {
        GLU.fail("Invalid framebuffer");
    }
    */
}

ProgressiveState.prototype.detach = function(fbo)
{
    var gl = GLU.gl;
    fbo.detachTexture(0);
    fbo.detachTexture(1);
}

ProgressiveState.prototype.clear = function(fbo)
{
    // clear radiance buffer
    var gl = GLU.gl;
    fbo.bind();
    fbo.drawBuffers(1);
    fbo.attachTexture(this.radianceTex, 0);
    gl.clearColor(0.0, 0.0, 0.0, 0.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    fbo.unbind();
}







