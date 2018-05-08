
/**
* Trinity is the global object providing access to all functionality in the system.
* @constructor
* @param {Simulation} sceneObj - The object defining the simulation setup
*/

var trinity;

var Trinity = function()
{
    trinity = this;

    this.initialized = false;
    this.terminated = false;
    this.rendering = false;
    
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
    this.textCtx = text_canvas.getContext("2d");
    this.onGravyLink = false;
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
    this.camControls.addEventListener('change', camChanged);
    this.camControls.keyPanSpeed = 100.0;

    this.gui = null;
    this.guiVisible = true;

    // Instantiate solver
    this.solver = new Solver();

    // Instantiate renderer
    this.renderer = new Renderer();
    
    // initialize
    this.init();
  
    // Do initial resize:
    this.resize();

    // Create dat gui
    this.gui = new GUI(this.guiVisible);

    // Setup keypress and mouse events
    window.addEventListener( 'mousemove', this, false );
    window.addEventListener( 'mousedown', this, false );
    window.addEventListener( 'mouseup',   this, false );
    window.addEventListener( 'contextmenu',   this, false );
    window.addEventListener( 'click', this, false );
    window.addEventListener( 'keydown', this, false );

    this.initialized = true;
}

/**
* Returns the current version number of the Trinity system, in the format [1, 2, 3] (i.e. major, minor, patch version)
*  @returns {Array}
*/
Trinity.prototype.getVersion = function()
{
    return [1, 0, 0];
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

Trinity.prototype.getSimulation = function()
{
    return this.simulation;
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
Trinity.prototype.getControls = function()
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



Trinity.prototype.init = function()
{
    // Initialize solver
    this.solver.resize(256, 256);
    
    // Initialize camera
    let domain = this.solver.getDomain();   
    let min     = domain.boundsMin;
    let max     = domain.boundsMax;
    let center  = domain.center;
    let extents = [    (max[0]-min[0]),     (max[1]-min[1]),     (max[2]-min[2])];
    let relDist = [1.5, 0.75, 1.3];
    this.camera.position.set(center[0]+relDist[0]*extents[0], 
                             center[1]+relDist[1]*extents[1], 
                             center[2]+relDist[2]*extents[2]);
    this.camControls.target.set(center[0], center[1], center[2]);
}


// Renderer reset on camera or other parameters update
Trinity.prototype.reset = function()
{
    if (!this.initialized || this.terminated) return;
    this.solver.reset();
    this.renderer.reset();
}
   

// Timestep the simulation
Trinity.prototype.step = function()
{
    // Run a simulation timestep
    this.solver.step();
    
    // Volume render simulation 
    this.renderer.render(this.solver);
}


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
    text_canvas.height = height

    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.camControls.update();


}

Trinity.prototype.resize = function()
{

}


Trinity.prototype.onClick = function(event)
{
    
}

Trinity.prototype.onDocumentMouseMove = function(event)
{


    this.camControls.update();
}

Trinity.prototype.onDocumentMouseDown = function(event)
{
    this.camControls.update();
}

Trinity.prototype.onDocumentMouseUp = function(event)
{
    this.camControls.update();
}

Trinity.prototype.onDocumentRightClick = function(event)
{

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
            this.guiVisible = !this.guiVisible;
            gravy.getGUI().toggleHide();
            break;
        
        case 79: // O key: output scene settings code to console
            let code = this.dumpScene();
            console.log(code);
            break;
    }
}



// Renderer reset on camera or other parameters update
Trinity.prototype.reset = function(no_recompile = false)
{
    if (!this.initialized || this.terminated) return;

    this.renderer.reset(no_recompile);
}

function camChanged()
{
    if (!trinity.rendering)
    {
        var no_recompile = true;
        trinity.reset(no_recompile);
    }
}
