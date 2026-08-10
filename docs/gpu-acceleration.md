# GPU fluid acceleration

ChemCAD has an optional OpenGL compute backend in `gfx::FluidGpuSolver`. It is not a replacement for `fluid::Solver`: the CPU solver remains the reference implementation and the runtime keeps using it unless the compute probe, shader compilation, and buffer setup all succeed.

## Work moved to the GPU

Particle position, velocity, acceleration, number density, pressure, colour field, interface normal, and phase index are stored in shader-storage buffer objects (SSBOs). They remain there across pressure iterations and substeps. The CPU downloads the particle buffer only when the caller publishes a rendering or diagnostic snapshot.

The compute passes perform:

- counting-sort uniform-grid construction: parallel cell counts, a prefix pass, and parallel scatter;
- Wendland C2 number-density gathering over the 27 neighbouring cells;
- colour-field gradient with the Bonet-Lok kernel-gradient correction, and the curvature that follows from it;
- harmonic-pair physical viscosity;
- Continuum Surface Force interfacial tension, `sigma * kappa * grad c`, with the interface thresholds compiled into the shader prelude from `fluid::kInterfaceGradientFloor` and `fluid::kInterfaceCorrectionDeterminant` so the two paths cannot disagree about where an interface is;
- PCISPH prediction, wall-density correction, non-negative pressure update, and corrected symmetric pressure force;
- frame, Coriolis, Euler, and centrifugal acceleration;
- speed/displacement clamping, four-pass analytic-SDF contact projection, and integration.

The vessel profile and boundary-density lookup table are sampled from `fluid::VesselBoundary` at the same 1,025 and 4,097 points used by the CPU implementation, then uploaded once per boundary change. Pressure prediction uses the same frozen wall tangent and frozen neighbour set as the CPU path. The final contact projection scans the complete profile rather than substituting a cylinder or box.

## Work retained on the CPU

Substep selection, CFL/acceleration/transport bounds, PCISPH convergence and stall decisions, frame-motion evaluation, and pressure-stiffness calibration remain on the CPU. They are low-volume control work, contain useful validated reference code, and do not justify a second GPU implementation. A four-word reduction buffer carries only maximum speed, acceleration, compression, deficit, rejection, and clamp counters back to the CPU; particle state does not make a per-substep round trip.

The optional XSPH display-only velocity smoothing mode is not ported. It defaults to zero. Requesting a non-zero `SolverConfig::xsphSmoothing` makes GPU setup unavailable so the caller uses the CPU solver rather than silently changing the requested model.

## Numerical equivalence and determinism

Both paths use the same kernel equations, support radius, density-contrast pressure operator, phase masses, harmonic viscosity, surface-tension model, boundary-density correction, substep ladder, iteration limits, pressure relaxation, and transport clamps. Interfacial tension needs no calibration on either path: `sigma` enters in N/m and `tests/test_fluid_solver.cpp` holds the CPU implementation to `dp = 2 sigma / R` directly.

The GPU path is not bit-identical to the CPU path. GLSL accumulates in 32-bit float, while the CPU reference performs reductions in double, and atomic counting-sort scatter does not promise the CPU solver's ascending-particle neighbour order. These differences can move the last few bits and the exact iteration at which a tolerance is crossed. They do not deliberately change the physical model. Reproducibility tests continue to use `fluid::Solver`.

## Requirements and fallback

The backend requires:

- OpenGL 4.3 core or newer (compute shaders and SSBOs are core; no vendor extension or CUDA is required);
- `GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS >= 128` and an x work-group size of at least 128;
- at least 12 compute-stage SSBO blocks and 12 SSBO binding points;
- `GL_MAX_SHADER_STORAGE_BLOCK_SIZE >= 16 MiB`;
- the compute, memory-barrier, SSBO-binding, buffer-clear, and buffer-mapping entry points loaded by `gfx::loadGl`.

`gfx::probeCompute()` is safe before GL initialisation. No loader, no current context, an older driver, insufficient limits, allocation failure, or any shader compile/link failure leaves `FluidGpuSolver::available()` false. The caller must then continue with the existing CPU solver. This is also the normal behaviour of headless tests.

## Measured

RTX 2070 Max-Q, 610.74, separatory funnel, 120 mL of water over dichloromethane, Quality budget (0.5% density tolerance, up to 40 pressure iterations). `tests/test_fluid_gpu.cpp` produces these numbers on every run and asserts the contracts behind them.

| spacing | particles | speedup | GPU compression | CPU compression |
| --- | --- | --- | --- | --- |
| 6 mm | 556 | 1.3x | 0.48% | 10.5% (3 stalled solves) |
| 3 mm | 4444 | 2.6x | 0.50% | 0.92% (5 stalled solves) |

Two things that measurement settled. Speedup below roughly 600 particles is dispatch and convergence-readback overhead rather than throughput, so `makeFluidAccelerator` leaves smaller charges on the CPU and the Interactive preset never reaches the device. And the device is not merely as accurate as the reference here but more so: it takes more pressure iterations per substep and reaches the requested tolerance, where the CPU solver's plateau detector gives up short of it.

An earlier measurement that showed the GPU at 2-5x WORSE compression was an artefact of capping both solvers at 12 pressure iterations; that cap bound the GPU and not the CPU, so it compared one converged solve against one truncated one.

## Runtime integration

`fluid::Simulation` owns a dedicated physics worker, and an OpenGL context belongs to one thread at a time -- the renderer keeps the main one. The application therefore creates a second, invisible, unshared context and hands it to `gfx::makeFluidAccelerator`, which binds it on the worker thread the first time a charge is large enough to be worth the trip. Nothing is shared between the two contexts: the compute buffers are private and results come back as ordinary memory.

The simulation re-uploads whenever the vessel, resolution, materials or charge change, falls back to `fluid::Solver` for any step the device declines, and names the device in its status line. `CHEMCAD_NO_FLUID_GPU=1` keeps the solve on the CPU. A machine without a 4.3 context gets a 3.3 one and the CPU solver, not an error.

