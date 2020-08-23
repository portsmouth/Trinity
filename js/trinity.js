
/**
* Trinity is the global object providing access to all functionality in the system.
* @constructor
* @param {Simulation} sceneObj - The object defining the simulation setup
*/

var trinity;

var Trinity = function(editor, error_editor)
{
    trinity = this;

    this.editor = editor;
    this.error_editor = error_editor;
    $(this.error_editor.getWrapperElement()).hide();

    this.editor_docs = {};
    this.editor_docs['common']    = CodeMirror.Doc('', "x-shader/x-fragment");
    this.editor_docs['initial']   = CodeMirror.Doc('', "x-shader/x-fragment");
    this.editor_docs['inject']    = CodeMirror.Doc('', "x-shader/x-fragment");
    this.editor_docs['influence'] = CodeMirror.Doc('', "x-shader/x-fragment");
    this.editor_docs['collide']   = CodeMirror.Doc('', "x-shader/x-fragment");
    this.editor_docs['render']    = CodeMirror.Doc('', "x-shader/x-fragment");

    var doc = this.editor_docs['initial'];
    this.editor.swapDoc(doc);

    let container = document.getElementById("container");
    this.container = container;

    var render_canvas = document.getElementById('render-canvas');
    this.render_canvas = render_canvas;
    this.width = render_canvas.width;
    this.height = render_canvas.height;
    render_canvas.style.width = render_canvas.width;
    render_canvas.style.height = render_canvas.height;

    var text_canvas = document.getElementById('text-canvas');
    this.text_canvas = text_canvas;
    this.text_canvas.style.width = render_canvas.width;
    this.text_canvas.style.height = render_canvas.height;
    this.textCtx = text_canvas.getContext("2d");
    this.onTrinityLink = false;
    this.onUserLink = false;

    window.addEventListener( 'resize', this, false );

    // Setup THREE.js orbit camera
    var VIEW_ANGLE = 45;
    var ASPECT = this.width / this.height;
    var NEAR = 1.0;
    var FAR = 20000;
    this.camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT, NEAR, FAR);
    this.camera.up.set(0, 1, 0);
    this.camera.position.set(1, 1, 1);
    this.camera.lookAt(new THREE.Vector3(0, 0, 0));

    this.camControls = new THREE.OrbitControls(this.camera, this.container);
    this.camControls.zoomSpeed = 2.0;
    this.camControls.flySpeed = 0.01;
    this.disable_reset = false;
    this.camControls.addEventListener('change',
        function() {
            if (trinity.disable_reset) return;
            trinity.render_dirty();
        });
    this.camControls.keyPanSpeed = 100.0;

    this.camControls.saveState = function () {
        this.target0.copy( this.target );
        this.position0.copy( this.object.position );
        this.zoom0 = this.object.zoom;
    };

    this.gui = null;
    this.guiVisible = true;

    // Instantiate renderer
    this.renderer = new Renderer();
    let RENDERER = this.renderer;

    // Instantiate solver
    this.solver = new Solver();

    // Field presets
    this.presets = new Presets();
    this.preset_selection = 'None';

    // Create dat gui
    this.gui = new GUI(this.guiVisible);

    // Setup codemirror events:
    let TRINITY = this;
    let SOLVER = this.solver;
    this.editor.on("change", function(cm, n) {
        if (!TRINITY.loading)
        {
            SOLVER.compileShaders();
            RENDERER.compileShaders();
        }
    });
    this.editing = false;
    this.loading = false;

    // Setup keypress and mouse events
    window.addEventListener( 'mousemove', this, false );
    window.addEventListener( 'mousedown', this, false );
    window.addEventListener( 'mouseup',   this, false );
    window.addEventListener( 'contextmenu',   this, false );
    window.addEventListener( 'click', this, false );
    window.addEventListener( 'keydown', this, false );

    // Attempt to take preset from current URL
    if (!this.load_url(window.location.href))
    {
        // Otherwise default to one of the presets
        this.presets.load_preset('Plume & sphere collider I');
    }
}

/**
* Returns the current version number of the Trinity system, in the format [1, 2, 3] (i.e. major, minor, patch version)
*  @returns {Array}
*/
Trinity.prototype.getVersion = function()
{
    return [1, 0, 0];
}

/**
* Access to the GUI object
*  @returns {GUI}
*/
Trinity.prototype.getGUI = function()
{
    return this.gui;
}

Trinity.prototype.getRenderer = function()
{
    return this.renderer;
}

Trinity.prototype.getSolver = function()
{
    return this.solver;
}

/**
* Access to the camera object
* @returns {THREE.PerspectiveCamera}.
*/
Trinity.prototype.getCamera = function()
{
    return this.camera;
}

/**
* Access to the camera controller object
* @returns {THREE.OrbitControls}
*/
Trinity.prototype.getCameraControls = function()
{
    return this.camControls;
}

/**
 * @returns {WebGLRenderingContext} The webGL context
 */
Trinity.prototype.getGLContext = function()
{
    return GLU.gl;
}

/**
* Programmatically show or hide the dat.GUI UI
* @param {Boolean} show - toggle
*/
Trinity.prototype.showGUI = function(show)
{
    this.guiVisible = show;
}


///////////////////////////////////////////////////////////////
// Camera logic
///////////////////////////////////////////////////////////////

Trinity.prototype.resetCam = function()
{
    let domain = this.getSolver().getDomain();
    let min     = domain.boundsMin;
    let max     = domain.boundsMax;
    let center  = domain.boundsCenter;
    let extents = [    (max[0]-min[0]),     (max[1]-min[1]),     (max[2]-min[2])];
    let relDist = [1.0, 0.666, 5.0];
    this.camera.position.set(center[0]+relDist[0]*extents[0], 
                             center[1]+relDist[1]*extents[1], 
                             center[2]+relDist[2]*extents[2]);
    this.camControls.target.set(center[0], center[1], center[2]);
    this.camControls.update();
}


///////////////////////////////////////////////////////////////
// App state management
///////////////////////////////////////////////////////////////

Trinity.prototype.get_escaped_stringified_state = function(state)
{
    let json_str = JSON.stringify(state);
    var json_str_escaped = json_str.replace(/[\\]/g, '\\\\')
                                    .replace(/[\b]/g, '\\b')
                                    .replace(/[\f]/g, '\\f')
                                    .replace(/[\n]/g, '\\n')
                                    .replace(/[\r]/g, '\\r')
                                    .replace(/[\t]/g, '\\t');
    return json_str_escaped;
}

Trinity.prototype.get_state = function()
{
    let camPos = this.camera.position;
    let camTar = this.camControls.target;
    let camera_settings = { pos: [camPos.x, camPos.y, camPos.z],
                            tar: [camTar.x, camTar.y, camTar.z],
                            near: this.camera.near,
                            far:  this.camera.far
    };

    let editor_settings = {
        common_glsl:     trinity.getGlsl('common'),
        initial_glsl:    trinity.getGlsl('initial'),
        inject_glsl:     trinity.getGlsl('inject'),
        influence_glsl:  trinity.getGlsl('influence'),
        collide_glsl:    trinity.getGlsl('collide'),
        render_glsl:     trinity.getGlsl('render')
    };

    // Extract user simulation settings from GUI parameters
    let simulationParameters = {};
    let simulation_folder = this.getGUI().simulationFolder;
    if (simulation_folder != undefined)
    {
        simulationParameters = simulation_folder.parameters;
    }

    let gui_settings = { visible: this.getGUI().visible };
    let state = { RENDERER_STATE:   this.renderer.settings,
                  SOLVER_STATE:     this.solver.settings,
                  SIMULATION_STATE: simulationParameters,
                  CAMERA_STATE:     camera_settings,
                  GUI_STATE:        gui_settings,
                  EDITOR_STATE:     editor_settings };

    return state;
}

Trinity.prototype.load_scene = function()
{
    var input = document.createElement('input');
    input.type = 'file';
    input.accept = '.json';
    input.onchange = e =>
    {
        let file = e.target.files[0];
        var reader = new FileReader();
        reader.readAsText(file,'UTF-8');
        reader.onload = readerEvent => {
            var file_contents = readerEvent.target.result;
            let state = null;
            try {
                state = JSON.parse(file_contents);
            } catch(e) {
                console.log("Failed to parse JSON from file ", file.name, ": ", e);
                return;
            }
            if (!this.load_state(state))
                console.log("Failed to load scene from file ", file.name);
            console.log("Loaded scene from file ", file.name);
         }
    }
    input.click();
}

Trinity.prototype.getQueryParam = function(url, key)
{
    var queryStartPos = url.indexOf('?');
    if (queryStartPos === -1) {
        return null;
    }
    var params = url.substring(queryStartPos + 1).split('&');
    for (var i = 0; i < params.length; i++) {
        var pairs = params[i].split('=');
        if (decodeURIComponent(pairs.shift()) == key) {
            return decodeURIComponent(pairs.join('='));
        }
    }
    return null;
}

Trinity.prototype.load_url = function(url)
{
    let URL = url;
    let preset_name = this.getQueryParam(URL, 'preset');
    preset_name = preset_name.replace(/^"|"$/g, '');
    if (!preset_name) return false;
    return this.presets.load_preset(preset_name);
}

Trinity.prototype.save_scene = function()
{
    let state = this.get_state();
    let json_str = JSON.stringify(state);
    //let objJsonStr = this.get_escaped_stringified_state(state);
    var currentdate = new Date(); 
    var datetime = currentdate.getDate() + "-" + (currentdate.getMonth()+1)  + "-" + currentdate.getFullYear() + "_"  
                 + currentdate.getHours() + "-" + currentdate.getMinutes() + "-" + currentdate.getSeconds();
    let filename = `trinity-scene-${datetime}.json`;
    let link = document.createElement('a');
    link.download = filename;
    let blob = new Blob([json_str], {type: 'text/csv'});
    link.href = URL.createObjectURL(blob);
    let event = new MouseEvent('click');
    link.dispatchEvent(event);link.dispatchEvent(event);
    console.log("Saved scene to file ", file);
}

Trinity.prototype.load_state = function(state)
{
    this.loading = true;

    if (!('CAMERA_STATE'     in state)) return false;
    if (!('GUI_STATE'        in state)) return false;
    if (!('EDITOR_STATE'     in state)) return false;
    if (!('RENDERER_STATE'   in state)) return false;
    if (!('SOLVER_STATE'     in state)) return false;
    if (!('SIMULATION_STATE' in state)) return false;

    // Load camera state
    let camera_settings = state.CAMERA_STATE;
    let P = camera_settings.pos;
    let T = camera_settings.tar;
    let near = camera_settings.near;
    let far = camera_settings.far;
    this.camera.position.copy(new THREE.Vector3(P[0], P[1], P[2]));
    this.camera.near = near;
    this.camera.far = far;
    this.camControls.target.copy(new THREE.Vector3(T[0], T[1], T[2]));
    this.camControls.saveState();
    this.initial_camera_position = this.camera.position.clone();
    this.initial_camera_target = this.camControls.target.clone();

    // Load GUI settings
    let gui_settings = state.GUI_STATE;
    if (gui_settings)
    {
        this.showGUI(gui_settings.visible);
    }

    // Load editor state
    if (state.EDITOR_STATE.common_glsl     === undefined) state.EDITOR_STATE.common_glsl  = '';
    if (state.EDITOR_STATE.initial_glsl    === undefined) state.EDITOR_STATE.initial_glsl = '';
    if (state.EDITOR_STATE.inject_glsl     === undefined) state.EDITOR_STATE.inject_glsl  = '';
    if (state.EDITOR_STATE.influence_glsl  === undefined) state.EDITOR_STATE.influence_glsl  = '';
    if (state.EDITOR_STATE.collide_glsl    === undefined) state.EDITOR_STATE.collide_glsl = '';
    if (state.EDITOR_STATE.render_glsl     === undefined) state.EDITOR_STATE.render_glsl = '';
    trinity.getDoc('common').setValue( state.EDITOR_STATE.common_glsl);
    trinity.getDoc('initial').setValue(state.EDITOR_STATE.initial_glsl);
    trinity.getDoc('inject').setValue( state.EDITOR_STATE.inject_glsl);
    trinity.getDoc('influence').setValue( state.EDITOR_STATE.influence_glsl);
    trinity.getDoc('collide').setValue(state.EDITOR_STATE.collide_glsl);
    trinity.getDoc('render').setValue(state.EDITOR_STATE.render_glsl);

    // Load renderer state
    this.renderer.settings = Object.assign(this.renderer.settings, state.RENDERER_STATE);

    // Load solver state
    let Nx = state.SOLVER_STATE.Nx;
    let Ny = state.SOLVER_STATE.Ny;
    let Nz = state.SOLVER_STATE.Nz;
    this.solver.resize(Nx, Ny, Nz);
    this.solver.settings = Object.assign(this.solver.settings, state.SOLVER_STATE);

    // Load GUI simulation settings
    if (!jQuery.isEmptyObject(state.SIMULATION_STATE))
    {
        let simulationParameters = {};
        simulationParameters = Object.assign(simulationParameters, state.SIMULATION_STATE);
        this.getGUI().simulationFolder = { parameters: simulationParameters };
    }

    // Sync logic
    this.camControls.update();
    this.gui.refresh();

    this.loading = false;

    // Do initial shader compilation (and start solver)
    this.solver.compileShaders();
    this.renderer.compileShaders();
    this.solver.restart();

    this.resize();
    return true;
}

///////////////////////////////////////////////////////////////
// user code management
///////////////////////////////////////////////////////////////


Trinity.prototype.getDoc = function(program_name)
{
    if (program_name in this.editor_docs)
    {
        return this.editor_docs[program_name];
    }
    return null;
}

Trinity.prototype.getGlsl = function(program_name)
{
    if (program_name in this.editor_docs)
    {
        let doc = this.editor_docs[program_name];
        return doc.getValue();
    }
    return '';
}

Trinity.prototype.show_errors = function()
{
    $(this.error_editor.getWrapperElement()).show();
}

Trinity.prototype.hide_errors = function()
{
    if (this.getRenderer().compiled_successfully &&
        this.getSolver().compiled_successfully)
    {
        $(this.error_editor.getWrapperElement()).hide();
    }
}

Trinity.prototype.num_lines  = function(str)
{
    return str.split(/\r\n|\r|\n/).length;
}

Trinity.prototype.compile_error = function(shaderName, shaderTypeStr, error_log)
{
    this.show_errors();

    this.error_editor.setValue('');
    this.error_editor.clearHistory();

    // copied from Solver.prototype.compileShader
    let common_glsl  = trinity.getGlsl('common') + '\n';
    let collide_glsl = common_glsl  + trinity.getGlsl('collide')  + '\n';

    let program_line_starts =
    {
        // (these numbers are the line numbers on which _USER_CODE_ appears)
        'initial' : this.num_lines(common_glsl)  + 21,
        'inject'  : this.num_lines(collide_glsl) + 30,
        'advect'  : this.num_lines(collide_glsl) + 36,
        'volume'  : this.num_lines(collide_glsl) + 59
    }

    program_line_start = 0;
    if (program_line_starts.hasOwnProperty(shaderName))
        program_line_start = program_line_starts[shaderName];

    // Map internal shader names to user-facing codeMirror doc names
    if (shaderName == 'advect') shaderName = 'influence';
    if (shaderName == 'volume') shaderName = 'render';

    let prefix = "[" + shaderName + " " + shaderTypeStr + " shader]";
    let errStr = '';

    const errorRE = /\d+:(\d+):/;
    error_log.split('\n').forEach( function(error, ndx) {
        const m = errorRE.exec(error);
        if (m)
        {
            let lineNum = Math.max(m ? parseInt(m[1]) : 0, 0) - program_line_start;
            lineNum += 1; // (since lines are 1-based in codeMirror)
            error = error.replace(errorRE, "");
            error = error.replace('ERROR:', "");
            error = error.trim();
            if (error)
                errStr += prefix + ' Error on line ' + lineNum + ': ' + error + '\n';
        }
    });

    console.log(errStr);
    this.error_editor.setValue(errStr);
}

Trinity.prototype.link_error = function(program_info, error_log)
{
    console.log("Link error: ");
    console.log("\t\tprogram_info: ", program_info);
    console.log("\t\terror_log: ", error_log);
}


///////////////////////////////////////////////////////////////
// Simulation controls
///////////////////////////////////////////////////////////////

Trinity.prototype.pauseSimToggle = function()
{
    this.solver.pauseToggle();
}

Trinity.prototype.restartSim = function()
{
    this.solver.restart();
}

// Timestep the simulation
Trinity.prototype.runFrame = function()
{
    // Run a simulation timestep
    if (!this.solver.paused)
    {
        this.solver.step();
        this.render_dirty();
    }

    // Volume render
    this.renderer.render(this.solver);

    this.update_hud();
}

Trinity.prototype.render_dirty = function()
{
    this.renderer.reset();
}

Trinity.prototype.update_hud = function()
{
    // Update HUD text canvas
    if (this.textCtx)
    {
        this.textCtx.textAlign = "left";   	// This determines the alignment of text, e.g. left, center, right
        this.textCtx.textBaseline = "middle";	// This determines the baseline of the text, e.g. top, middle, bottom
        this.textCtx.font = '12px monospace';	// This determines the size of the text and the font family used
        this.textCtx.clearRect(0, 0, this.textCtx.canvas.width, this.textCtx.canvas.height);
        this.textCtx.globalAlpha = 0.95;
        this.textCtx.strokeStyle = 'black';
        this.textCtx.lineWidth  = 2;
        if (this.getGUI().visible)
        {
            if (this.onTrinityLink) this.textCtx.fillStyle = "#ff5500";
            else                    this.textCtx.fillStyle = "#ffff00";
            let ver = this.getVersion();
            let linkWidth = this.textCtx.measureText('Trinity vX.X.X').width;
            let trinity_version_str = 'Trinity v'+ver[0]+'.'+ver[1]+'.'+ver[2];
            this.textCtx.strokeText(trinity_version_str, this.textCtx.canvas.width - linkWidth - 14, this.textCtx.canvas.height-20);
            this.textCtx.fillText(  trinity_version_str, this.textCtx.canvas.width - linkWidth - 14, this.textCtx.canvas.height-20);
            {
                // status text
                this.textCtx.fillStyle = "#ffaa22";
                {
                    let spp           = this.renderer.spp;
                    let max_spp       = this.renderer.settings.max_spp;
                    let timestep      = this.solver.frame;
                    let max_timesteps = this.solver.settings.max_timesteps;
                    let time          = Number(this.solver.time.toPrecision(3));
                    let status_txt = '';
                    if (!this.solver.paused)
                        status_txt = 'timestep ' + timestep + '/' + max_timesteps + ", time " + time;
                    else
                        status_txt = 'timestep ' + timestep + '/' + max_timesteps + ", time " + time + ' [' + spp + '/' + max_spp + ' spp]';
                    this.textCtx.strokeText(status_txt, 14, this.textCtx.canvas.height-25);
                    this.textCtx.fillText(  status_txt, 14, this.textCtx.canvas.height-25);
                }

                // status text
                this.textCtx.fillStyle = "#ccaaff";
                {
                    let help_txt = '[SPC to pause/play, ESC to restart]';
                    this.textCtx.strokeText(help_txt, 14, this.textCtx.canvas.height-42);
                    this.textCtx.fillText(  help_txt, 14, this.textCtx.canvas.height-42);
                }
            }
        }
    }
}



///////////////////////////////////////////////////////////////
// Window resize logic
///////////////////////////////////////////////////////////////

Trinity.prototype._resize = function(width, height)
{
    this.width = width;
    this.height = height;

    let render_canvas = this.render_canvas;
    render_canvas.width  = width;
    render_canvas.height = height;
    render_canvas.style.width = width;
    render_canvas.style.height = height;

    var text_canvas = this.text_canvas;
    text_canvas.width  = width;
    text_canvas.height = height;
    text_canvas.style.width = width;
    text_canvas.style.height = height;

    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.camControls.update();

    this.renderer.resize(width, height);
}

Trinity.prototype.resize = function()
{
    let width = window.innerWidth;
    let height = window.innerHeight;
    this._resize(width, height);
}


Trinity.prototype.handleEvent = function(event)
{
    switch (event.type)
    {
        case 'resize':      this.resize();  break;
        case 'mousemove':   this.onDocumentMouseMove(event);  break;
        case 'mousedown':   this.onDocumentMouseDown(event);  break;
        case 'mouseup':     this.onDocumentMouseUp(event);    break;
        case 'contextmenu': this.onDocumentRightClick(event); break;
        case 'click':       this.onClick(event);  break;
        case 'keydown':     this.onkeydown(event);  break;
    }
}

Trinity.prototype.onClick = function(event)
{
    event.preventDefault();
}

Trinity.prototype.onDocumentMouseMove = function(event)
{
    let u = event.clientX/window.innerWidth;
    let v = event.clientY/window.innerHeight;

    // check not within editor region
    var cmEl = document.querySelector('.CodeMirror');
    var edRect = cmEl.getBoundingClientRect();
    if (event.clientX >= edRect.left && event.clientX <= edRect.right &&
        event.clientY >= edRect.top  && event.clientY <= edRect.bottom)
    {
        this.disable_reset = true;
        this.camControls.reset();
        this.camControls.update();
        this.disable_reset = false;
        this.camControls.enabled = false;
        return;
    }

    // Check whether user is trying to click the Trinity home link
    var textCtx = this.textCtx;
    if (textCtx)
    {
        var x = event.pageX;
        var y = event.pageY;
        let linkWidth = this.textCtx.measureText('Trinity vX.X.X').width;

        let xmin = this.textCtx.canvas.width - linkWidth - 14;
        let xmax = xmin + linkWidth;
        let ymin = this.textCtx.canvas.height-25;
        let ymax = this.textCtx.canvas.height-10;
        if (x>=xmin && x<=xmax && y>=ymin && y<=ymax) this.onTrinityLink = true;
        else this.onTrinityLink = false;
    }

    this.camControls.update();
    this.camControls.saveState();
}

Trinity.prototype.onDocumentMouseDown = function(event)
{
    this.camControls.update();
}

Trinity.prototype.onDocumentMouseUp = function(event)
{
    this.camControls.update();
    event.preventDefault();
}

Trinity.prototype.onDocumentRightClick = function(event) {}

Trinity.prototype.onClick = function(event)
{
    if (this.onTrinityLink)
    {
        window.open("https://github.com/portsmouth/Trinity");
    }
    event.preventDefault();
}



Trinity.prototype.onkeydown = function(event)
{
    var charCode = (event.which) ? event.which : event.keyCode;
    switch (charCode)
    {
        case 122: // F11 key: go fullscreen
            var element	= document.body;
            if      ( 'webkitCancelFullScreen' in document ) element.webkitRequestFullScreen();
            else if ( 'mozCancelFullScreen'    in document ) element.mozRequestFullScreen();
            else console.assert(false);
            break;

        case 72: // H key: toggle hide/show dat gui
        if (!this.camControls.enabled || trinity.editing) break;
            this.guiVisible = !this.guiVisible;
            trinity.getGUI().toggleHide();
            break;

        case 79: // O key: dump JSON application state to console
            if (!this.camControls.enabled || trinity.editing) break;
            console.log(this.get_escaped_stringified_state(this.get_state()));
            break;

        case 83: // S key: save application state to file
            if (!this.camControls.enabled || trinity.editing) break;
            this.save_scene();
            break;

        case 70: // F key: move cam to standard orientation
            if (!this.camControls.enabled || trinity.editing) break;
            this.resetCam();
            break;

        case 32: // space bar: pause toggle
            if (!this.camControls.enabled || trinity.editing) break;
            this.pauseSimToggle();
            break;

        case 27: // escape key: restart sim from t=0
            this.restartSim();
            break;

    }
}




