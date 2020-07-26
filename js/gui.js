
/** @constructor 
* Interface to the dat.GUI UI.
*/
var GUI = function(visible = true) 
{
	// Create dat gui
	this.gui = new dat.GUI();
	this.gui.domElement.id = 'gui';
	var gui = this.gui;
	this.visible = visible;

	this.createguiSettings();
	if (!visible)
		this.gui.__proto__.constructor.toggleHide();
}

function updateDisplay(gui) 
{
    for (var i in gui.__controllers) {
        gui.__controllers[i].updateDisplay();
    }
    for (var f in gui.__folders) {
        updateDisplay(gui.__folders[f]);
    }
}

/**
* Call to explicitly force the GUI to synchronize with the
* current parameter values, if they have been changed programmatically.
*/
GUI.prototype.sync = function()
{
	updateDisplay(this.gui);
}

dat.GUI.prototype.removeFolder = function(name) {
	var folder = this.__folders[name];
	if (!folder) {
	  return;
	}
	folder.close();
	this.__ul.removeChild(folder.domElement.parentNode);
	delete this.__folders[name];
	this.onResize();
  }

  
GUI.prototype.refresh = function()
{
    if (this.presetsFolder    != undefined) this.gui.removeFolder(this.presetsFolder.name);
    if (this.simulationFolder != undefined) this.gui.removeFolder(this.simulationFolder.name);
    if (this.rendererFolder   != undefined) this.gui.removeFolder(this.rendererFolder.name);

    this.createguiSettings();
}

GUI.prototype.toggleHide = function()
{
	this.visible = !this.visible;
}

function hexToRgb(hex) 
{
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}


GUI.prototype.createguiSettings = function()
{
    this.guiSettings = {};

    // Presets folder
    this.presetsFolder = this.gui.addFolder('Presets');
    let preset_names = trinity.presets.get_preset_names();
    this.presetSettings = {};
    this.presetSettings["preset"] = trinity.preset_selection;
    var presetItem = this.presetsFolder.add(this.presetSettings, 'preset', preset_names);
    presetItem.onChange(function(preset_name) { trinity.presets.load_preset(preset_name); });

	this.createSimulationSettings();
	this.createRendererSettings();
}

GUI.prototype.createRendererSettings = function()
{
    this.rendererFolder = this.gui.addFolder('Raymarcher');
    var renderer = trinity.getRenderer();
    var camera = trinity.getCamera();

    this.rendererFolder.add(renderer.settings, 'exposure', -10.0, 10.0);
    this.rendererFolder.add(renderer.settings, 'gamma', 1.0, 3.0);
    this.rendererFolder.add(renderer.settings, 'debrisExtinction', 0.0, 100.0);
    this.rendererFolder.add(renderer.settings, 'blackbodyEmission', -20.0, 20.0);
    this.rendererFolder.add(renderer.settings, 'TtoKelvin', 0.0, 10.0);

	this.rendererFolder.close();
}

GUI.prototype.createSimulationSettings = function()
{
    var solver = trinity.getSolver();

    this.simulationFolder = this.gui.addFolder('Simulation');

    this.simulationFolder.Nx = solver.settings.Nx;
    this.simulationFolder.Ny = solver.settings.Ny;
    this.simulationFolder.Nz = solver.settings.Nz;

    this.simulationFolder.add(this.simulationFolder, 'Nx', 32, 512).onChange(                function(Nx) { solver.resize(Math.floor(Nx), solver.settings.Ny, solver.settings.Nz); } );
    this.simulationFolder.add(this.simulationFolder, 'Ny', 32, 512).onChange(                function(Ny) { solver.resize(solver.settings.Nx, Math.floor(Ny), solver.settings.Nz); } );
    this.simulationFolder.add(this.simulationFolder, 'Nz', 32, 512).onChange(                function(Nz) { solver.resize(solver.settings.Nx, solver.settings.Ny, Math.floor(Nz)); } );
    this.simulationFolder.add(solver.settings, 'NprojSteps', 1, 256).onChange(         function(NprojSteps) { solver.NprojSteps = NprojSteps; solver.restart(); } );
    this.simulationFolder.add(solver.settings, 'timestep', 0.0, 10.0).onChange(        function() { solver.restart(); } );
    this.simulationFolder.add(solver.settings, 'vorticity_scale', 0.0, 0.99).onChange( function() { solver.restart(); } );

    // @todo: move to user code
    this.simulationFolder.add(solver, 'blastHeight', 0.0, 1.0).onChange(           function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'blastRadius', 0.0, 1.0).onChange(           function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'blastTemperature', 0.0, 100000.0).onChange( function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'blastVelocity', 0.0, 1000.0).onChange(      function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'debrisHeight', 0.0, 1.0).onChange(          function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'debrisFalloff', 0.0, 1.0).onChange(         function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'gravity', 0.0, 1.0).onChange(               function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'T0', 0.0, 1000.0).onChange(                 function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'Tambient', 0.0, 1000.0).onChange(           function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'buoyancy', 0.0, 1.0).onChange(              function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'expansion', 0.0, 1.0).onChange(             function() { solver.restart(); } );
    this.simulationFolder.add(solver, 'radiationLoss', 0.0, 1.0).onChange(         function() { solver.restart(); } );

    let button = { nuke:function(){ solver.reset(); }};
    let button_ui = this.simulationFolder.add(button, 'nuke').name('NUKE');
    button_ui.__li.className = 'cr link footer';

    this.simulationFolder.open();
}

/**
 * Add a dat.GUI UI slider to control a float parameter.
 * The scene parameters need to be organized into an Object as
 * key-value pairs, for supply to this function.
 * @param {Object} parameters - the parameters object for the scene, with a key-value pair (where value is number) for the float parameter name
 * @param {Object} param - the slider range for this parameter, in the form `{name: 'foo', min: 0.0, max: 100.0, step: 1.0, recompile: true}` (step is optional, recompile is optional [default is false])
 * @param {Object} folder - optionally, pass the dat.GUI folder to add the parameter to (defaults to the main scene folder)
 * @returns {Object} the created dat.GUI slider item
 * @example
 *		Scene.prototype.initGui = function(gui)
 *		{
 *			gui.addSlider(this.parameters, c);
 *			gui.addSlider(this.parameters, {name: 'foo2', min: 0.0, max: 1.0});
 *			gui.addSlider(this.parameters, {name: 'bar', min: 0.0, max: 3.0, recompile: true});
 *		}
 */
GUI.prototype.addSlider = function(parameters, param, folder=undefined)
{
	let _f = this.userFolder;
	if (typeof folder !== 'undefined') _f = folder;
	var name = param.name;
	var min  = param.min;
	var max  = param.max;
	var step = param.step;
	var recompile = param.recompile;
	var no_recompile = true;
	if (!(recompile==null || recompile==undefined)) no_recompile = !recompile;
	var item;
	if (step==null || step==undefined) { item = _f.add(parameters, name, min, max, step); }
	else                               { item = _f.add(parameters, name, min, max);       }
	item.onChange( function(value) { gravy.reset(no_recompile); gravy.camera.enabled = false; } );
	item.onFinishChange( function(value) { gravy.camera.enabled = true; } );
	return item;
}

/** 
 * Add a dat.GUI UI color picker to control a 3-element array parameter (where the RGB color channels are mapped into [0,1] float range)
 * @param {Object} parameters - the parameters object for the scene, with a key-value pair (where value is a 3-element array) for the color parameter name
 * @param {Object} name - the color parameter name
 * @param {Object} folder - optionally, pass a scale factor to apply to the RGB color components to calculate the result (defaults to 1.0)
 * @param {Object} folder - optionally, pass the dat.GUI folder to add the parameter to (defaults to the main scene folder)
 * @returns {Object} the created dat.GUI color picker item
*/
GUI.prototype.addColor = function(parameters, name, scale=1.0, folder=undefined)
{
	let _f = this.userFolder;
	if (typeof folder !== 'undefined') _f = folder;
	_f[name] = [parameters[name][0]*255.0, parameters[name][1]*255.0, parameters[name][2]*255.0];
	var item = _f.addColor(_f, name);
	item.onChange( function(color) {
								if (typeof color==='string' || color instanceof String)
								{
									var C = hexToRgb(color);
									parameters[name][0] = scale * C.r / 255.0;
									parameters[name][1] = scale * C.g / 255.0;
									parameters[name][2] = scale * C.b / 255.0;
								}
								else
								{
									parameters[name][0] = scale * color[0] / 255.0;
									parameters[name][1] = scale * color[1] / 255.0;
									parameters[name][2] = scale * color[2] / 255.0;
								}
								gravy.reset(true);
							} );
	return item;
}

// (deprecated)
GUI.prototype.addParameter = function(parameters, param)
{
	this.addSlider(parameters, param);
}

/**
* Access to internal dat.GUI object
* @returns {dat.GUI}
*/
GUI.prototype.getGUI = function()
{
	return this.gui;
}

/**
* Access to dat.GUI object folder object containing user UI parameters
* @returns {dat.GUI}
*/
GUI.prototype.getUserFolder = function()
{
	return this.userFolder;
}



