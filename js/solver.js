
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
    });

    this.fbo = new GLU.RenderTarget();
    this.quadVbo = this.createQuadVbo();
    this.compileShaders();
}

Solver.prototype.compileShaders = function()
{
	this.solveEnergy   = new GLU.Shader('solve-energy', this.shaderSources, null);
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

Solver.prototype.initialize = function()
{
	let Nr = this.Nr;
	let Ny = this.Ny;

	this.U_rho_tex   = null;
	this.U_rho_p_tex = null;
	this.U_E_tex     = null;
	this.U_v_tex     = null;

	// thermal medium conservation vector:
	let Q0 = new Float32Array(Nr*Ny*4);

	// passive particulate matter density (A = advected particles)
	let A0 = new Float32Array(Nr*Ny*4);

	let R = 1000; // full domain radius in metres
	let H = 3000; // full domain height in metres
	let dr = R/Nr;
	let dy = H/Ny;

	let p0 = 1.0e5;                      // ground-level atmospheric pressure, in N / m^2
	let T0 = 300.0;                      // ground-level ambient temperature, in K
	let Rs = 287.0;                      // specific gas constant, in J / (kg K)
	let g = 9.81;                        // gravitational acceleration, in m / s^2
	let gamma = 1.4;                     // ratio of specific heats, cp/cv
	let ground_depth = 10.0;             // depth of ground particulate matter
	let ground_density = 1.0;            // number density of ground particulate matter, arbitrary units
	                                     // (extinction coefficient will be given by some adjustable proportionality factor)
	let lapse = Rs * T0 / g;             // atmospheric temperature scale height, in m
	let rho0 = p0 / (Rs * T0);           // ground-level atmosphere density, in kg / m^3
	let cv = Rs / (gamma-1.0);           // specific heat capacity (constant volume), in J / (K gg)
	let c = Math.sqrt(gamma*p0/rho0);    // atmospheric sound speed in m/s
	let dt = 0.5 * Math.min(dr, dy) / c; // approximate timestep satisfying Courant condition

	for (let iy=0; iy<Ny; ++iy)
	{
		// ground height (domain bottom) is y = 0
		let y = (0.5 + iy)*dy;
		for (var ir=0; ir<Nr; ++ir)
    	{
    		let r = (0.5 + ir)*dr;

    		// initial atmosphere profile
    		let adiabatic_falloff = 1.0 -  ((gamma-1.0)/gamma)*(y/lapse);
    		let T = T0 * adiabatic_falloff;
    		let p = p0 * Math.pow(adiabatic_falloff, gamma/(gamma-1.0));
    		let e = cv * T;                // specific internal energy, e
    		let rho = p / ((gamma-1.0)*e); // atmosphere density follows from equation of state:
    		let E = rho * e;               // assume initially static, so total energy density E is:

    		// Initialize thermal medium:
    		let index = iy*Nr + ir;
 
    		Q0[4*index+0] = rho;
    		Q0[4*index+1] = E;
    		Q0[4*index+2] = 0.0; // rho * v_r
    		Q0[4*index+3] = 0.0; // rho * v_z

    		// Initialize particulate matter density:
    		A0[4*index+0] = (y < ground_depth) ? ground_density : 0.0;
    		A0[4*index+0] = 0.0; // (unused for now)
    		A0[4*index+1] = 0.0; // (unused for now)
    		A0[4*index+2] = 0.0; // (unused for now)
    	}
	}

	this.Q = [ new GLU.Texture(Nr, Ny, 4, true, false, true, Q0), new GLU.Texture(Nr, Ny, 4, true, false, true, Q0) ];
	this.A = [ new GLU.Texture(Nr, Ny, 4, true, false, true, A0), new GLU.Texture(Nr, Ny, 4, true, false, true, A0) ];

	this.timestep = 0;

	this.domain = { boundsMin: [-R, 0.0, -R],
	                boundsMax: [ R,  H,   R],
	                center:    [0.0, H/2.0, 0.0],
	                radius: R,
	                height: H,
	                Nr: Nr,
	                Ny: Ny,
	                cv: cv }

	// constants
	this.dr = dr;
	this.dy = dy;
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

Solver.prototype.getQ = function()
{
	let previousIndex = this.timestep % 2;
	return this.Q[previousIndex];
}	


Solver.prototype.step = function()
{
	console.log("Running timestep: ", this.timestep);

	let gl = GLU.gl;

	let currIndex = (this.timestep+0) % 2;
	let nextIndex = (this.timestep+1) % 2;

	// Set shader uniforms
	this.solveEnergy.bind();

	this.solveEnergy.uniformI("Nr",      this.Nr);
	this.solveEnergy.uniformI("Ny",      this.Ny);
	this.solveEnergy.uniformF("dt",      this.dt);
	this.solveEnergy.uniformF("dr",      this.dr);
	this.solveEnergy.uniformF("dy",      this.dy);
	this.solveEnergy.uniformF("gammaMO", this.gammaMO);
	this.solveEnergy.uniformF("g",       this.g);

	// Run timestep update
	gl.viewport(0, 0, this.Nr, this.Ny);
	this.fbo.bind();
	this.fbo.drawBuffers(1);
	this.fbo.attachTexture(this.Q[nextIndex], 0);            // write to Q[nextIndex]
    this.Q[currIndex].bind(0);                               // read from Q[currIndex]
    this.solveEnergy.uniformTexture("Q", this.Q[currIndex]);
    
    this.quadVbo.bind();
    this.quadVbo.draw(this.solveEnergy, gl.TRIANGLE_FAN);

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
