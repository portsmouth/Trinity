
var _presets_table = {

    'Basic plume': `
    {
        "R":{"exposure":-1,"gamma":2.2,"blackbodyEmission":-9.5,"debrisExtinction":20,"TtoKelvin":1},
        "S":{"timestep":1,"NprojSteps":16,"vorticity_scale":0.2,"Nx":128,"Ny":512,"Nz":128},
        "C":{"pos":[255.99999999999994,640,230.4],"tar":[64,256,64],"near":1,"far":20000},
        "G":{"visible":true},
        "E":{"code":"\\nbool isSolid(in vec3 wsP, // world space point of current voxel\\n             in vec3 L)   // world-space extents of grid\\n{\\n    // define regions which are solid (static) obstacles\\n    return false;\\n}"}
    }`

};


var Presets = function()
{
    this.preset_names = [];
    for (var preset_name in _presets_table) {
        if (_presets_table.hasOwnProperty(preset_name)) {
            this.preset_names.push(preset_name);
        }
    }
}

Presets.prototype.get_preset_names = function()
{
    return this.preset_names;
}

Presets.prototype.get_preset = function(preset_name)
{
    return this.preset_names[preset_name];
}

Presets.prototype.load_preset = function(preset_name)
{
    if (preset_name in _presets_table)
    {
        let preset = _presets_table[preset_name];
        let state = JSON.parse(preset);
        trinity.preset_selection = preset_name;
        trinity.load_state(state);
    }
}

