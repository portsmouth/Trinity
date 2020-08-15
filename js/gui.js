
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


GUI.prototype.wipe = function()
{
    if (this.presetsFolder    != undefined) this.gui.removeFolder(this.presetsFolder.name);
    this.gui.removeFolder('Simulation');
    if (this.solverFolder     != undefined) this.gui.removeFolder(this.solverFolder.name);
    if (this.rendererFolder   != undefined) this.gui.removeFolder(this.rendererFolder.name);
}

GUI.prototype.refresh = function()
{
    if (this.simulationFolder != undefined)
        this.simulationFolderParametersCached = jQuery.extend(true, {}, this.simulationFolder.parameters);

    this.wipe();
    this.addPresetsFolder();
    this.generateSimulationSettings();
    this.createSolverSettings();
    this.createRendererSettings();
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

GUI.prototype.generateSimulationSettings = function()
{
    let simulationFolder = this.gui.addFolder('Simulation');
    this.simulationFolder = simulationFolder;
    simulationFolder.parameters = {};

    let glsl = '';
    glsl += trinity.getGlsl('common') + '\n';
    glsl += trinity.getGlsl('initial') + '\n';
    glsl += trinity.getGlsl('inject') + '\n';
    glsl += trinity.getGlsl('influence') + '\n';
    glsl += trinity.getGlsl('collide') + '\n';

    // Parse code for uniform bindings:
    let SOLVER = trinity.getSolver();
    let GUI = this;
    let CACHED_PARAMS = this.simulationFolderParametersCached;

    let lines = glsl.split("\n");
    lines.forEach(function(line) {
        line = line.trim();
        if (line.startsWith('uniform'))
        {
            // find beginning of comment
            let bind_statement = line.substr(line.indexOf('//')+2);
            try {
                let bind_object   = JSON.parse(bind_statement);
                let param_name    = bind_object.name;
                let param_default = bind_object.default;
                if ((param_name    !== undefined) &&
                    (param_default !== undefined))
                {
                    if (typeof param_default == 'number')
                    {
                        // add numeric slider
                        let min_val = 0.0;
                        let max_val = param_default * 5.0;
                        let step_val = 0.001;
                        if (typeof bind_object.min == 'number') min_val = bind_object.min;
                        if (typeof bind_object.max == 'number') max_val = bind_object.max;
                        if (typeof bind_object.step == 'number') step_val = bind_object.step;

                        let param_value = param_default;
                        if (CACHED_PARAMS != undefined)
                        {
                            if (CACHED_PARAMS.hasOwnProperty(param_name))
                                param_value = CACHED_PARAMS[param_name];
                        }
                        simulationFolder.parameters[param_name] = param_value;

                        let syncToShader = true;
                        GUI.addSlider(simulationFolder.parameters,
                                                   {name: param_name, min: min_val, max: max_val, step: step_val},
                                                   simulationFolder,
                                                   syncToShader);
                        SOLVER.syncFloatToShader(param_name, param_value);
                    }
                    else if (Array.isArray(param_default) && param_default.length==3)
                    {
                        // add color picker
                        let scale = 1.0;
                        if (typeof bind_object.scale == 'number')
                            scale = bind_object.scale;

                        let param_value = param_default;
                        if (CACHED_PARAMS != undefined)
                        {
                            if (CACHED_PARAMS.hasOwnProperty(param_name))
                                param_value = CACHED_PARAMS[param_name];
                        }
                        simulationFolder.parameters[param_name] = param_value;

                        let syncToShader = true;
                        GUI.addColor(simulationFolder.parameters,
                                                  param_name,
                                                  scale,
                                                  simulationFolder,
                                                  syncToShader);
                        SOLVER.syncColorToShader(param_name, param_value);
                    }
                }
            } catch(e) {
                console.log(e);
            }
        }
    });

    // Retain UI simulation parameters after a recompile
    if (this.simulationFolderParametersCached != undefined)
    {
        for (const param_name of Object.keys(this.simulationFolderParametersCached))
        {
            if (!simulationFolder.parameters.hasOwnProperty(param_name))
                continue;
            let param_value = this.simulationFolderParametersCached[param_name];
            simulationFolder.parameters[param_name] = param_value;
            if (typeof param_value == 'number')
            {
                SOLVER.syncFloatToShader(param_name, param_value);
            }
            if (Array.isArray(param_value) && param_value.length==3)
            {
                SOLVER.syncColorToShader(param_name, param_value);
            }
        }
        this.sync();
    }

    simulationFolder.open();
}

GUI.prototype.addPresetsFolder = function()
{
    // Presets folder
    this.presetsFolder = this.gui.addFolder('Presets');
    let preset_names = trinity.presets.get_preset_names();
    this.presetSettings = {};
    this.presetSettings["preset"] = trinity.preset_selection;
    var presetItem = this.presetsFolder.add(this.presetSettings, 'preset', preset_names);
    presetItem.onChange(function(preset_name) { trinity.presets.load_preset(preset_name); });
    this.presetsFolder.open();
}

GUI.prototype.createSolverSettings = function()
{
    let solver = trinity.getSolver();
    let renderer = trinity.getRenderer();

    this.solverFolder = this.gui.addFolder('Solver');

    this.solverFolder.Nx = solver.settings.Nx;
    this.solverFolder.Ny = solver.settings.Ny;
    this.solverFolder.Nz = solver.settings.Nz;

    this.solverFolder.add(this.solverFolder, 'Nx', 2, 2048).onChange(             function(Nx) { solver.resize(Math.floor(Nx), solver.settings.Ny, solver.settings.Nz); } );
    this.solverFolder.add(this.solverFolder, 'Ny', 2, 2048).onChange(             function(Ny) { solver.resize(solver.settings.Nx, Math.floor(Ny), solver.settings.Nz); } );
    this.solverFolder.add(this.solverFolder, 'Nz', 2, 2048).onChange(             function(Nz) { solver.resize(solver.settings.Nx, solver.settings.Ny, Math.floor(Nz)); } );
    this.solverFolder.add(solver.settings, 'NprojSteps', 1, 256).onChange(         function(NprojSteps) { solver.settings.NprojSteps = Math.floor(NprojSteps); trinity.render_dirty(); } );
    this.solverFolder.add(solver.settings, 'timestep', 0.0, 10.0).onChange(        function() { trinity.render_dirty(); } );
    this.solverFolder.add(solver.settings, 'max_timesteps', 0.0, 10000.0).onChange( function(max_timesteps) { solver.settings.max_timesteps = Math.floor(max_timesteps); trinity.render_dirty(); } );
    this.solverFolder.add(solver.settings, 'vorticity_scale', 0.0, 0.99).onChange( function() { trinity.render_dirty(); } );
    this.solverFolder.add(solver.settings, 'expansion', 0.0, 1.0).onChange(        function() { trinity.render_dirty(); } );

    this.solverFolder.open();
}

GUI.prototype.createRendererSettings = function()
{
    this.rendererFolder = this.gui.addFolder('Raymarcher');
    var renderer = trinity.getRenderer();
    var camera = trinity.getCamera();

    this.rendererFolder.add(renderer.settings, 'extinctionScale', -6.0, 6.0).onChange(function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'emissionScale', -10.0, 10.0).onChange(function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'anisotropy', -1.0, 1.0).onChange(     function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'exposure', -10.0, 10.0).onChange(     function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'gamma', 0.0, 3.0).onChange(           function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'saturation', 0.0, 3.0).onChange(      function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'sunLatitude', -90.0, 90.0).onChange(  function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'sunLongitude', 0.0, 360.0).onChange(  function() { trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'sunPower', -3.0, 3.0).onChange(       function() { trinity.render_dirty(); });

    this.rendererFolder.sunColor = [renderer.settings.sunColor[0]*255.0,
                                    renderer.settings.sunColor[1]*255.0,
                                    renderer.settings.sunColor[2]*255.0];
    let sunColorItem = this.rendererFolder.addColor(this.rendererFolder, 'sunColor');
    sunColorItem.onChange( function(C) {
                            if (typeof C==='string' || C instanceof String)
                            {
                                var color = hexToRgb(C);
                                renderer.settings.sunColor[0] = color.r / 255.0;
                                renderer.settings.sunColor[1] = color.g / 255.0;
                                renderer.settings.sunColor[2] = color.b / 255.0;
                            }
                            else
                            {
                                renderer.settings.sunColor[0] = C[0] / 255.0;
                                renderer.settings.sunColor[1] = C[1] / 255.0;
                                renderer.settings.sunColor[2] = C[2] / 255.0;
                            }
                            trinity.render_dirty();
                        } );

    this.rendererFolder.skyColor = [renderer.settings.skyColor[0]*255.0,
                                    renderer.settings.skyColor[1]*255.0,
                                    renderer.settings.skyColor[2]*255.0];
    let skyColorItem = this.rendererFolder.addColor(this.rendererFolder, 'skyColor');
    skyColorItem.onChange( function(C) {
                            if (typeof C==='string' || C instanceof String)
                            {
                                var color = hexToRgb(C);
                                renderer.settings.skyColor[0] = color.r / 255.0;
                                renderer.settings.skyColor[1] = color.g / 255.0;
                                renderer.settings.skyColor[2] = color.b / 255.0;
                            }
                            else
                            {
                                renderer.settings.skyColor[0] = C[0] / 255.0;
                                renderer.settings.skyColor[1] = C[1] / 255.0;
                                renderer.settings.skyColor[2] = C[2] / 255.0;
                            }
                            trinity.render_dirty();
                        } );

    this.rendererFolder.add(renderer.settings, 'Nraymarch', 8, 1024).onChange(       function(Nraymarch)     { renderer.settings.Nraymarch     = Math.floor(Nraymarch);     trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'max_spp', 1, 256).onChange(           function(max_spp)       { renderer.settings.max_spp       = Math.floor(max_spp);       trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'spp_per_frame', 1, 8).onChange(       function(spp_per_frame) { renderer.settings.spp_per_frame = Math.floor(spp_per_frame); trinity.render_dirty(); });
    this.rendererFolder.add(renderer.settings, 'show_bounds').onChange(               function() { trinity.render_dirty(); });

    this.rendererFolder.close();
}

GUI.prototype.addSlider = function(parameters, param, folder=undefined, syncToShader=false)
{
    let _f = this.userFolder;
    if (typeof folder !== 'undefined') _f = folder;
    var name = param.name;
    var min  = param.min;
    var max  = param.max;
    var step = param.step;
    var item;
    if (step==null || step==undefined) { item = _f.add(parameters, name, min, max, step); }
    else                               { item = _f.add(parameters, name, min, max);       }
    item.onChange( function(value) {
        trinity.camera.enabled = false;
        if (syncToShader) trinity.getSolver().syncFloatToShader(name, value);
        trinity.render_dirty();
    }
    );
    item.onFinishChange( function(value) {
        if (syncToShader) trinity.getSolver().syncFloatToShader(name, value);
        trinity.render_dirty();
        trinity.camera.enabled = true; } );
    return item;
}

GUI.prototype.addColor = function(parameters, name, scale=1.0, folder=undefined, syncToShader=false)
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
                                if (syncToShader) trinity.getSolver().syncColorToShader(name, parameters[name]);
                                trinity.render_dirty();
                            } );
    item.onFinishChange( function(value) {
        if (syncToShader) trinity.getSolver().syncColorToShader(name, parameters[name]);
        trinity.render_dirty();
        trinity.render_dirty();
    });
    return item;
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



