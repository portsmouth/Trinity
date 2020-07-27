
var _presets_table = {

    'Basic plume': `
{"R":{"Nraymarch":128,"exposure":-1,"gamma":2.2,"blackbodyEmission":-9.5,"debrisExtinction":36.39054331081162,"TtoKelvin":1},"S":{"timestep":1,"NprojSteps":16,"vorticity_scale":0.2,"Nx":128,"Ny":512,"Nz":128},"C":{"pos":[-619.0635843562742,426.2684503298603,647.5154673154354],"tar":[62.9244116822174,225.30502033792672,-58.19073281589177],"near":1,"far":20000},"G":{"visible":true},"E":{"common_glsl":"// \\"Physics\\"\\nuniform float gravity;                    // {\\"label\\":\\"gravity\\",                    \\"min\\":0.0, \\"max\\":0.1,    \\"step\\":0.001, \\"default\\":0.05}\\nuniform float buoyancy;                   // {\\"label\\":\\"buoyancy\\",                   \\"min\\":0.0, \\"max\\":1.0,    \\"step\\":0.001, \\"default\\":0.5}\\nuniform float Tambient;                   // {\\"label\\":\\"Tambient\\",                   \\"min\\":0.0, \\"max\\":1000.0, \\"step\\":0.001, \\"default\\":300.0}\\nuniform float Tbuoyancy;                  // {\\"label\\":\\"Tbuoyancy\\",                  \\"min\\":0.0, \\"max\\":1000.0, \\"step\\":0.001, \\"default\\":300.0}\\nuniform float radiationLoss;              // {\\"label\\":\\"radiationLoss\\",              \\"min\\":0.0, \\"max\\":1.0,    \\"step\\":0.001, \\"default\\":0.999}\\n\\n// Blast geometry\\nuniform float blast_height;               // {\\"label\\":\\"blast_height\\",               \\"min\\":0.0, \\"max\\":1.0,    \\"step\\":0.001, \\"default\\":0.25}\\nuniform float blast_radius;               // {\\"label\\":\\"blast_radius\\",               \\"min\\":0.0, \\"max\\":1.0,    \\"step\\":0.001, \\"default\\":0.1}\\nuniform float blast_velocity;             // {\\"label\\":\\"blast_velocity\\",             \\"min\\":0.0, \\"max\\":100.0,  \\"step\\":0.001, \\"default\\":50.0}\\nuniform float blast_temperature_contrast; // {\\"label\\":\\"blast_temperature_contrast\\", \\"min\\":0.0, \\"max\\":1000.0, \\"step\\":0.001, \\"default\\":100.0}\\n\\nuniform float debris_inflow_rate;         // {\\"label\\":\\"debris_inflow_rate\\",         \\"min\\":0.0, \\"max\\":10.0,   \\"step\\":0.01, \\"default\\":1.0}\\nuniform vec3 debris_albedo;               // {\\"label\\":\\"debris_albedo\\",              \\"default\\":[0.5, 0.5, 0.5], \\"scale\\":1.0}\\n\\nvoid init()\\n{\\n    // Any global constants defined here are available in all functions\\n}","initial_glsl":"// Specify velocity, temperature, and debris density/albedo at time=0\\nvoid initial_conditions(in vec3 wsP,                // world space point of current voxel\\n                        in vec3 L,                  // world-space extents of grid\\n                        inout vec3 velocity,        // initial velocity\\n                        inout float temperature,    // initial temperature\\n                        inout vec3 density,         // initial per-channel debris extinction\\n                        inout vec3 albedo)          // initial per-channel debris albedo\\n{\\n    velocity = vec3(0.0);\\n    temperature = Tambient;\\n    density = vec3(0.0);\\n    albedo = vec3(0.0);      \\n}\\n\\n","inject_glsl":"void inject(in vec3 wsP,        // world space point of current voxel\\n            in float time,      // time in units of frames\\n            in vec3 L,          // world-space extents of grid\\n            inout vec3 dv,      // injected velocity in voxels/frame (defaults to 0)\\n            inout float T,      // temperature updated in-place for heating/cooling (no update by default)\\n            inout vec3 drho,    // injected per-channel debris extinction (defaults to 0)\\n            inout vec3 albedo)  // albedo of the injected debris (defaults to grey)\\n{\\n    vec3 blast_center = vec3(0.5*L.x, blast_height*L.y, 0.5*L.z);\\n    vec3 dir = wsP - blast_center;\\n    float r = length(dir);\\n    dir /= r;\\n    float rt = r/(blast_radius*L.y);\\n    if (rt <= 1.0 && time==0.0)\\n    {\\n        // Within blast radius: add velocity and set temperature\\n        float radial_falloff = max(0.0, 1.0 - rt*rt*(3.0 - 2.0*rt));\\n        dv = dir * blast_velocity * radial_falloff;\\n        T = Tambient * (1.0 + blast_temperature_contrast*radial_falloff);\\n        drho = vec3(1.0) * debris_inflow_rate * radial_falloff;\\n        albedo = debris_albedo;\\n    }\\n  \\telse\\n  \\t{\\n        // Apply thermal relaxation to ambient temperature Tambient due to \\"radiation loss\\"\\n        float dT = T - Tambient;\\n        T = Tambient + dT*radiationLoss;\\n    }\\n}\\n","advect_glsl":"\\nvec3 externalForces(in vec3 wsP,                       // world space point of current voxel\\n                    in vec3 v, in float P, in float T, // velocity, pressure, temperature of current voxel\\n                    in vec3 L)                         // world-space extents of grid\\n{\\n    // Tbuoyancy is the reference temperature for buoyancy\\n    float buoyancy_force = -gravity + gravity*buoyancy*(T - Tbuoyancy);\\n    return vec3(0.0, buoyancy_force, 0.0);\\n}","collide_glsl":"bool isSolid(in vec3 wsP, // world space point of current voxel\\n             in vec3 L)   // world-space extents of grid\\n{\\n    // define regions which are solid (static) obstacles\\n    return false;\\n}\\n"}}`

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

