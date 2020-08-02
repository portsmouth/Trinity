
/** @constructor
*/
var Renderer = function()
{
	let gl = GLU.gl;

    // Default user-adjustable properties
    this.settings = {};

    // Raymarch resolution
    this.settings.Nraymarch = 256;

    // Lighting
    this.settings.skyColor = [0.5, 0.5, 1.0];
    this.settings.sunColor = [1.0, 1.0, 0.5];
    this.settings.sunPower = 1.0;
    this.settings.sunLatitude   = 60.0;
    this.settings.sunLongitude  = 0.0;

    // Tonemapping
    this.settings.exposure = -1.0;
    this.settings.gamma = 2.2;

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

Renderer.prototype.resize = function(width, height)
{
    this._width = width;
    this._height = height;
    this.quadVbo = this.createQuadVbo();
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
    let render_glsl = common_glsl + trinity.getGlsl('render');

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
    // @todo:  render only on dirty
    if (!this.compiled_successfully)
        return;

	let gl = GLU.gl;

    gl.disable(gl.DEPTH_TEST);
    gl.viewport(0.0, 0.0, this._width, this._height);
	gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

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

    // Volume render
    {
        this.volumeProgram.bind();

        // Bind textures for air temperature and medium absorption/scattering
        absorption.bind(0);
        scattering.bind(1);
        Tair.bind(2);
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
        this.volumeProgram.uniformF("anisotropy", this.settings.anisotropy);
        this.volumeProgram.uniformF("extinctionScale", Math.pow(10.0, this.settings.extinctionScale));
        this.volumeProgram.uniformF("emissionScale", Math.pow(10.0, this.settings.emissionScale));

        // Tonemapping
        this.volumeProgram.uniformF("exposure", this.settings.exposure);
        this.volumeProgram.uniformF("invGamma", 1.0/this.settings.gamma);

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

        this.quadVbo.bind();
        this.quadVbo.draw(this.volumeProgram, gl.TRIANGLE_FAN);
    }

    // draw simulation bounds
    {
        this.lineProgram.bind();

        // Setup projection matrix
        var projectionMatrix = camera.projectionMatrix.toArray();
        var projectionMatrixLocation = this.lineProgram.getUniformLocation("u_projectionMatrix");
        gl.uniformMatrix4fv(projectionMatrixLocation, false, projectionMatrix);

        this.lineProgram.uniform3Fv("color", [0.5, 0.0, 0.0]);

        // Setup modelview matrix (to match camera)
        camera.updateMatrixWorld();
        var matrixWorldInverse = new THREE.Matrix4();
        matrixWorldInverse.getInverse( camera.matrixWorld );
        var modelViewMatrix = matrixWorldInverse.toArray();
        var modelViewMatrixLocation = this.lineProgram.getUniformLocation("u_modelViewMatrix");
        gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

        this.boxVbo.bind();
        this.boxVbo.draw(this.lineProgram, gl.LINES);
    }

    gl.finish();
}









