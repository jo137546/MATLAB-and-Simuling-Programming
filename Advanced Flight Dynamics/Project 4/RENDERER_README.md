# Project 4 arcade orbital renderer

1. Extract this entire folder.
2. For an immediate preview, double-click **Orbit_Renderer.html**. The supplied preview is clearly labelled and uses an independent integration of the same deterministic CW/LQ model; it is not a recorded Simulink run.
3. In MATLAB, make this folder the Current Folder and run **Project4_Main**. After simulation, it replaces the preview with your actual `time` and `x` histories and opens it in your system browser.

If you already have a completed simulation in the MATLAB workspace, use:

```matlab
Render_Orbits(time, x, rho_f)
```

No web server, internet connection, JavaScript package install, or additional MATLAB graphics toolbox is needed. Chrome or Edge with WebGL enabled is suitable. The simulation still requires your existing MATLAB/Simulink and control tools.

The scene opens on a visible low-poly Earth, with the leader reference orbit, equatorial reference and ascending/descending nodes. Switch to Formation for a close-up in the leader's rotating LVLH frame. Use the playback controls to play/pause, change speed, restart or scrub to any recorded time. Left-drag to orbit the camera, right-drag or Shift-drag to pan, and scroll to zoom. Reset camera restores the original framing.

The arcade renderer uses a faceted Earth, angular spacecraft and flat colours. Pixel size 2 is the default; choose 3 or 4 for larger pixels and less rendering work, or 1 for finer edges. It caps the drawing buffer at 960 by 720 pixels, disables antialiasing, uses simple vertex lighting, and redraws a paused scene only when something changes. Telemetry updates at 10 Hz while playing. These changes reduce graphics work; actual frame rate depends on the browser/device and was not measured here.

Satellite shapes and orientation are illustrative: the simulation supplies relative translation, not spacecraft dimensions or attitude. At Earth scale, the pair is represented by one enlarged spacecraft glyph because the true separation is smaller than a pixel. The two recorded trails retain their real positions. Use Formation to see both spacecraft and their separation in metres.

The circular leader orbit uses `a=6898 km`, `i=23 deg`, `RAAN=10 deg`, `argument of periapsis=90 deg`, and `initial true anomaly=0`. The LVLH axes are radial outward, along track and orbit normal. With `Q=Rz(RAAN)*Rx(i)*Rz(argument of periapsis)` and `theta=nu0+sqrt(mu/a^3)*(t-t0)`:

```text
eR = Q*[cos(theta); sin(theta); 0]
eT = Q*[-sin(theta); cos(theta); 0]
eN = Q*[0; 0; 1]
rLeader   = a*eR
rFollower = rLeader + x*eR + y*eT + z*eN
```

This is an inertial reconstruction of the existing CW relative-motion model. It does not add a separate nonlinear orbital propagation model. The Earth surface is illustrative, not a dated geographical/ground-track reference. The given radius is already a low Earth orbit: 6898 - 6378.1 = 519.9 km altitude. The leader period is about 5701.59 seconds. Only the `tf` line was changed in the main script for this update, making a run last 5702 seconds (approximately one complete orbit); set `tf=1000` for a shorter run. Rendering cannot supply controlled follower motion beyond the simulated interval.

Ascending and descending nodes mark the **leader's reference orbit** crossing the ECI equatorial plane, z=0. At the ascending node the leader moves from negative to positive z; at the descending node it moves from positive to negative z. Their positions are `a*[cos(RAAN);sin(RAAN);0]` and its negative. With RAAN=10 degrees, the node directions are 10 and 190 degrees from the ECI x axis. These are not separately estimated nodes of a controlled follower osculating orbit.

The controller, noise and Kalman-filter design are unchanged by this rendering addition. Preview trajectories do not certify the homework performance requirements.

Set `render_orbits = false` near the top of `Project4_Main.m` to skip the browser export. The existing MATLAB plots remain available. For very long runs, the export keeps up to 12000 evenly spaced recorded samples; playback interpolates between them for display only.

The renderer mathematics, interpolation and control-handler logic were checked independently. A JavaScript test harness used stubbed graphics calls; it did not verify WebGL appearance. A live browser preview could not be completed because the computer controls could not reliably identify the browser URL. MATLAB could not be launched in the development environment because its MathWorks services connection returned error 5202, so the new MATLAB export call has not been executed there.
