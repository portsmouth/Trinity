
/** @constructor
*/
var Renderer = function()
{
	let gl = GLU.gl;

    // Default user-adjustable properties
    this.exposure = -2.0;
    this.gamma = 2.2;
    this.colorA = [1.0,0.8,0.5];
    this.colorB = [1.0,0.8,0.5];
    this.blackbodyEmission = -1.0;
    this.debrisExtinction = 10.0;
    this.tempMultiplier = 1.0;
    this.TtoKelvin = 4.0;

    // Internal buffers and programs
    this.boxVbo = null;
    this.quadVbo = this.createQuadVbo();
    this.fbo = new GLU.RenderTarget();

    this.volumeProgram    = null;
    this.sliceProgram     = null;
    this.lineProgram      = null;
    this.tonemapProgram   = null;

    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource({
        'volume'         : {'v': 'volume-vertex-shader',  'f': 'volume-fragment-shader'},
        'tonemap'        : {'v': 'tonemap-vertex-shader', 'f': 'tonemap-fragment-shader'},
        'line'           : {'v': 'line-vertex-shader',    'f': 'line-fragment-shader'}
    });

	var render_canvas = trinity.render_canvas;
    render_canvas.width  = window.innerWidth;
    render_canvas.height = window.innerHeight;
    this._width = render_canvas.width;
    this._height = render_canvas.height;

    this.compileShaders();

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

Renderer.prototype.reset = function(no_recompile = false)
{
    this.numFramesSinceReset = 0;
    if (!no_recompile) this.compileShaders();
}


Renderer.prototype.resize = function(width, height)
{
    this._width = width;
    this._height = height;

    if (this.fbo)
	    this.fbo.unbind();
    this.quadVbo = this.createQuadVbo();
    //this.fbo = new GLU.RenderTarget();

    this.reset(true);
}


/**
* Restart accumulating samples.
* @param {Boolean} [no_recompile=false] - set to true if shaders need recompilation too
*/
Renderer.prototype.reset = function(no_recompile = false)
{
    if (!no_recompile) this.compileShaders();
}

Renderer.prototype.compileShaders = function()
{
	this.volumeProgram  = new GLU.Shader('volume',  this.shaderSources, null);
	this.lineProgram    = new GLU.Shader('line',    this.shaderSources, null);
	this.tonemapProgram = new GLU.Shader('tonemap', this.shaderSources, null);
}

Renderer.prototype.render = function(solver)
{
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
    let Vair   = solver.getVair();    // (the air velocity texture)
    let Tair   = solver.getTair();    // (the air temperature texture)
    let debris = solver.getDebris(); // (the debris density texture)

    // Volume render
    {
        this.volumeProgram.bind();

        debris.bind(0);
        this.volumeProgram.uniformTexture("debris_sampler", debris);

        Tair.bind(1);
        this.volumeProgram.uniformTexture("Tair_sampler", Tair);

        Vair.bind(2);
        this.volumeProgram.uniformTexture("Vair_sampler", Vair);

        this.volumeProgram.uniform2Fv("resolution", [this._width, this._height]);
        this.volumeProgram.uniform3Fv("camPos", [camPos.x, camPos.y, camPos.z]);
        this.volumeProgram.uniform3Fv("camDir", [camDir.x, camDir.y, camDir.z]);
        this.volumeProgram.uniform3Fv("camX", [camX.x, camX.y, camX.z]);
        this.volumeProgram.uniform3Fv("camY", [camUp.x, camUp.y, camUp.z]);
        this.volumeProgram.uniformF("camFovy", camera.fov);
        this.volumeProgram.uniformF("camAspect", camera.aspect);
        this.volumeProgram.uniform3Fv("volMin", domain.boundsMin);
        this.volumeProgram.uniform3Fv("volMax", domain.boundsMax);
        this.volumeProgram.uniform3Fv("volCenter", domain.center);

        this.volumeProgram.uniformI("N", domain.N);
        this.volumeProgram.uniformF("dL", domain.dL);
        this.volumeProgram.uniformI("Nraymarch", 256);

        this.volumeProgram.uniformF("debrisExtinction", this.debrisExtinction);
        this.volumeProgram.uniformF("blackbodyEmission", Math.pow(2.0, this.blackbodyEmission));
        this.volumeProgram.uniformF("TtoKelvin", this.TtoKelvin);
        this.volumeProgram.uniformF("exposure", this.exposure);
        this.volumeProgram.uniformF("invGamma", 1.0/this.gamma);

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

        this.lineProgram.uniform3Fv("color", [1.0, 0.0, 0.0]);

        // Setup modelview matrix (to match camera)
        camera.updateMatrixWorld();
        var matrixWorldInverse = new THREE.Matrix4();
        matrixWorldInverse.getInverse( camera.matrixWorld );
        var modelViewMatrix = matrixWorldInverse.toArray();
        var modelViewMatrixLocation = this.lineProgram.getUniformLocation("u_modelViewMatrix");
        gl.uniformMatrix4fv(modelViewMatrixLocation, false, modelViewMatrix);

        if (!this.boxVbo)
        {
            let o = domain.boundsMin;
            let c = domain.boundsMax;
            let extents = [c[0]-o[0], c[1]-o[1], c[2]-o[2]];
            this.boxVbo = this.createBoxVbo(o, extents)
        }

        this.boxVbo.bind();
        this.boxVbo.draw(this.lineProgram, gl.LINES);
    }

    gl.finish();
}









