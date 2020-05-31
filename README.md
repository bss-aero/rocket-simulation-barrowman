# Rocket Simulation - Barrowman

The implementation was made in Matlab 2018. To run the files is necessary to install the software. Different versions may not be able to run the simulation.slx file correctly.

In order to simulate a new rocket flight follow the steps below:

  1. Inside inputsSimulink.m, change the rocket geometry, thrust curve, and environment settings. Run the file to load the inputs into the workspace.
  2. Open simulation.slx and run the file.
  3. The results can be seen inside DataView subsystem in the main view.

To change the wind settings is necessary to enter simulation.slx and go to "Translational motion > Wind setting > Wind". Inside the wind.m file you can set the wind as a function of time and position of the rocket.
