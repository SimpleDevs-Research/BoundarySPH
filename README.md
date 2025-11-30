# BoundarySPH

_Expected to run with Unity Ver. `2021.3.11f1`_

## Instructions for Installation

1. Clone this repository!
2. Before opening the repository as a Unity project, make sure to either clone or copy [UnityUtils V1.2](https://github.com/SimpleDevs-Tools/UnityUtils) into the `Assets/Scripts/` directory.
3. ... And that's it! You are done!

## How to Run

### Running Demo Scenes

Inside of `Assets/BSPH/Scenes/`, you have several demo scenes. These scenes should be able to run without any modification, so play away! 

### How Demo Scenes Work

If you open any of these scenes, you will notice the following:

#### `BufferManager.cs`

This script acts as an data management system for all data buffers that are passed to the GPU. The only things you may need to touch are:

- `PARTICLE_CONTROLLER`
- `BOIDS_CONTROLLER` - _Optional_
- `OBSTACLES_CONTROLLER`

Let's look at these references.

#### `PARTICLE_CONTROLLER`: `Simulation3D.cs`

This script actually controls the particles and simulation. The hierarchy gives you several references you need to manage as well as controls that affect the simulation proper.

This README only really covers the core params you will likely touch. Feel free to experiment with the other params that aren't mentioned here explicitly.

- **References**:
    - `BM`: `BufferManager.cs` - you have to link it back to the same `BufferManager` as above.
    - `SPAWNER`: References a `ParticleSpawner.cs` object in the scene.
    - `Compute`: This actually references `Assets/BSPH/Scripts/SPHCalculations.compute`, which is a compute shader. This shader handles all SPH-based calculations.
- **SPH Configurations**:
    - `Dt`: The delta time for the simulation, in seconds? Default = `0.00825`
    - `Particle Mass`: How much mass should each particle represent? Default = `1`
    - `external_force`: Basically any kind of external influence on the particles. Typically gravity, at 9.8 m/s^2 downward.
    - `kernel_radius`: The radius (m) of each particle
    - `rest_density`: The density of the fluid at rest.
    - `viscosity_influence`: the viscosity coefficient.
    - `pressure_influence`: a numeric weight to control how much pressure a particle exerts on particles that are on the outer parts of each particle's kernel of influence.
    - `near_pressure_influence`: a numeric weight to control how much pressure a particle exerts on particles that are close to each particle within the kernel sphere of influence.

#### `ParticleSpawner.cs`

This script is referenced by `Simulation3D.cs`. This script lets you control the spawn area of particles. Note that you don't designate how many particles are spawned; rather, this number is auto-calculated based on the dimensions you feed it inside of `num_particles_per_axis`.

To adjust the spawn behavior, you generally want to do the following:

1. Position the GameObject that this script is attached to within the 3D world space. All particles will spawn relative to the position of this GameObject's position.
2. Adjust `num_particles_per_axis` to determine the initial grid of particles to be spawned when the simulation starts.
3. As you adjust `num_particles_per_axis`, you can visually see the spawn area of the particles as a wire box Gizmo. If you don't see any such grid, make sure you have turned Gizmos on. You can also adjust the appearance of the wire box via `show_spawn_bounds` and `spawn_bounds_color`.

#### `OBSTACLES_CONTROLLER` = `MeshObsGPU.cs`

This script is where you would define all the meshed obstacles that you would want your particles to interact with. Similar to `Simulation3D.cs`, you have some references to be aware of and some values you can control via the Inspector:

- **References**:
    - `BM`: The Buffer Manager
    - `GRID`: A `Grid.cs` component. Can be attached to the same GameObject as your `MeshObsGPU.cs`.
    - `SHADER`: A compute shader to handle mesh interactions. By default, you want to set this to `Assets/BSPH/Scripts/MeshObsShader.compute`.
    - `PARTICLE_CONTROLELR`: The `Simulation3D.cs` script.
    - `BOIDS_CONTROLLER`: Optional, can be left empty.
- **Obstacles and Boids**:
    - `Obstacles`: An array of all meshed obstacles in your scene. Adjust this array to add/remove objects that will affect the simulation.
    - (You can ignore all the boid-related stuff)

... It's actually pretty simple all things considered.


