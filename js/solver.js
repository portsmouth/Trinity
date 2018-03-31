
    /*

 	// Assume the blast consists of an injection of energy at a point in ambient air

	// physical fields:
	//   T = temperature
	//   V = velocity,  (vr, vphi, vz)
	//   rho = density

	// Then:
	//     e = cv * T,   e = specific energy (energy per unit mass)
	//     E = e + |v|^2/2
	//     p = (gamma - 1) rho * e

	// Evolve  rho, V, E
	// (enforce dirichlet boundary conditions (V=0) on ground layer)
	// -> solve for T each time-step
	 
    // rho, T are used to volume render the incandescent air blackbody radiation

    // Other processes which produce the "smoke";

    //  - Ground debris:
    //  - 	contact of super heated blast with the ground vaporizes earth, which is sucked into the mushroom stem
    //  - 	can model this as a ground layer of dense particulate advected with the flow according to the continuity equation
    
    //  - nuclear isotopes from the weapon itself produce radioactive particulate matter
    //  - heating the air generates chemical reactions which produce various molecules, e.g. nitric oxide.
    //  - this generates "smoke" (i.e. opaque particulate matter). This can be modelled by injecting particulate matter
    //    into the flow according to the local temperature (requires some reasonable model of this, can be developed iteratively)


    // More ambitiously, later could attempt to incorporate nuclear physics/chemistry to determine the actual
    // isotopic constituents at each point. A simpler approach might be to assume the particulate matter is generated at the
    // blast point and subsequently advected/mixed with the surrounding air, which constitutes the radioactive component of the medium.

	// To extend to 3d solver,  since webGL does not allow rendering into 3d textures,
	// we can keep all slices in a 2d texture with appropriate indexing,
	// e.g. a 256^3 3d texture can be packed into a 4096^2 2d texture.
	// The caveat is, no hardware filtering, amd possibly slow texture fetches. 

	*/




/** @constructor
*/
var Solver = function()
{
    // Specify shaders
    this.shaderSources = GLU.resolveShaderSource({
    	'solve-energy'         : {'v': 'solve-energy-vertex-shader',   'f': 'solve-energy-fragment-shader'}
        'solve-velocity'       : {'v': 'solve-velocity-vertex-shader', 'f': 'solve-velocity-fragment-shader'}
    });

    this.fbo = new GLU.RenderTarget();

    this.compileShaders();

    this.resize(64, 64);
}

Solver.prototype.compileShaders = function()
{
	this.solveEnergy   = new GLU.Shader('solve-energy',    this.shaderSources, null);
	this.solveVelocity = new GLU.Shader('solve-velocity',  this.shaderSources, null);
}


Solver.prototype.resize = function(Nr, Nz)
{
    this.Nr = Nr;
    this.Nz = Nz;

    this.initialize(this.Nr, this.Nz);

    this.quadVbo = this.createQuadVbo();
}


Solver.prototype.initialize = function(Nr, Nz)
{
	this.U_rho_tex   = null;
	this.U_rho_p_tex = null;
	this.U_E_tex     = null;
	this.U_v_tex     = null;

	// thermal medium conservation vector:
	let Q0 = new Float32Array(Nr*Nz*4);

	// passive particulate matter density (A = advected particles)
	let A0 = new Float32Array(Nr*Nz*4);

	let R = 1000; // full domain radius in metres
	let H = 1000; // full domain height in metres
	let dr = R/Nr;
	let dz = H/Nz;

	let p0 = 1.0e5;   // ground-level atmospheric pressure, in N / m^2
	let T0 = 300.0;   // ground-level ambient temperature, in K
	let Rs = 287.0;   // specific gas constant, in J / (kg K)
	let rho0 = p0 / (Rs * T0);  // ground-level atmosphere density, in kg / m^3
	let g = 9.81;     // gravitational acceleration, in m / s^2
	let cv = 0.718e3; // specific heat capacity (constant volume), in J / (K gg)
	let gamma = 1.4;  // ratio of specific heats, cp/cv
	let ground_depth = 10.0; // depth of ground particulate matter
	let ground_density = 1.0; // number density of ground particulate matter, arbitrary units
	                          // (extinction coefficient will be given by some adjustable proportionality factor)
	let lapse = Rs * T0 / g; 

	// atmospheric sound speed
	let c = Math.sqrt(gamma*p0/rho0);
	let dt = 0.5 * Math.min(dr, dz) / c; // approximate timestep satisfying Courant condition

	for (let iz=0; iz<Nz; ++iz)
	{
		// ground height (domain bottom) is z = 0
		let z = iz*dH;
		for (var ir=0; ir<Nr; ++ir)
    	{
    		let index = iz*Nr + ir;
    		let r = ir*dR;

    		// initial atmosphere profile
    		let adiabatic_falloff = 1.0 -  ((gamma-1.0)/gamma)*(z/lapse);
    		let T = T0 * adiabatic_falloff;
    		let p = p0 * Math.pow(adiabatic_falloff, gamma/(gamma-1.0));

    		// atmosphere density follows from equation of state:
    		let rho = p / ((gamma-1.0)*e);

    		// specific internal energy, e
    		let e = cv * T;

    		// assume initially static, so total energy density E is:
    		let E = rho * e;

    		// Initialize thermal medium:
    		Q0[3*index+0] = rho;
    		Q0[3*index+1] = E;
    		Q0[3*index+2] = 0.0; // rho * v_r
    		Q0[3*index+3] = 0.0; // rho * v_z
   
    		// Initialize particulate matter density:
    		A0[3*index+0] = (z < ground_depth) ? ground_density : 0.0;
    		A0[3*index+0] = 0.0; // (unused for now)
    		A0[3*index+1] = 0.0; // (unused for now)
    		A0[3*index+2] = 0.0; // (unused for now)
    	}
	}

	this.Q = [ new GLU.Texture(Nr, Nz, 4, true, false, true, Q0), new GLU.Texture(Nr, Nz, 4, true, false, true, Q0) ];
	this.V = [ new GLU.Texture(Nr, Nz, 4, true, false, true, A0), new GLU.Texture(Nr, Nz, 4, true, false, true, V0) ];

	this.timestep = 0;

	// constants
	this.dr = dr;
	this.dz = dz;
	this.dt = dt;
	this.gammaMO = gamma-1.0; 
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

Solver.prototype.step = function()
{
	console.log("Running timestep: ", this.timestep);

	let gl = GLU.gl;

	let currIndex = (this.timestep+0) % 2;
	let nextIndex = (this.timestep+1) % 2;

	// Set shader uniforms
	this.solveEnergy.bind();
	this.solveEnergy.uniformF("dt",      this.dt);
	this.solveEnergy.uniformF("dr",      this.dr);
	this.solveEnergy.uniformF("dz",      this.dz);
	this.solveEnergy.uniformF("gammaMO", this.gammaMO);
	this.solveEnergy.uniformF("g",       this.g);

	// Run timestep update
	this.fbo.bind();
	this.fbo.drawBuffers(1);
    this.Q[currIndex].bind(0);                                   // read current frame V
    this.solveEnergy.uniformTexture("Q", this.Q[currIndex]);
    this.fbo.attachTexture(this.Q[nextIndex], 0);                 // write next frame V
    this.quadVbo.bind();
    this.quadVbo.draw(this.solveDensity, gl.TRIANGLE_FAN);
	this.fbo.unbind();
    gl.bindTexture(gl.TEXTURE_2D, null);

    // @todo: advect (and inject) particulate matter passively in the flow, 
    //        to model dust (and fission products)

    this.timestep++;
}
