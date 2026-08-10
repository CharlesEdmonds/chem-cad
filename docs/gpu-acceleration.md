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

Set `CHEMCAD_FLUID_GPU=0` (also accepts `false`, `off`, or `cpu`) before starting ChemCAD to force the CPU path. The environment check occurs before shader compilation.

## Expected speedup

The estimate is analytical, not a measured benchmark of this unintegrated backend. The CPU measurements recorded beside the solver are approximately 5.3 ms per substep for 379 particles and 14.6 ms for 926 particles, with neighbour gathers and repeated pressure iterations dominating. Those gathers expose one independent invocation per particle and roughly 30-45 neighbours, a good fit for a Turing GPU. Dispatch and scalar readback overhead dominate small charges, while arithmetic dominates larger charges.

On an RTX 2070/2080-class Turing GPU, a reasonable initial expectation is 1.5-4x end-to-end solver speedup around the 900-particle balanced workload and 4-8x once several thousand particles keep the GPU occupied. Charges below roughly 256 particles may see little or no gain. These ranges include the many compute dispatches and one small convergence readback per pressure iteration; they intentionally do not quote the much larger raw ALU-throughput ratio. Actual numbers must be measured after the orchestrator integrates the runtime switch and snapshot cadence.
