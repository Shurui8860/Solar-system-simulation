# Celestial Simulation Project

## Overview

This project simulates the motion of celestial bodies in a solar system-like environment using numerical integration techniques. It provides insights into orbital dynamics, including the comparison between simulated and actual orbital periods. The simulation also includes an animation of the celestial motions.

## Files and Directories

### Python Files

- ``

  - Serves as the central script for the project, orchestrating the execution of the simulation. It integrates the functionality from all other modules, initializes parameters, and displays the simulation results and animations.

- ``

  - Contains the `Body` class, which models celestial objects such as the Sun, Earth, and other planets. The class includes attributes like mass, position, velocity, radius, and color. It also provides methods for updating the state of the celestial body during the simulation.

- ``

  - Manages the loading and parsing of input data from external files, such as `demo.txt`. This module ensures that simulation parameters and celestial body properties are correctly initialized.

- ``

  - Provides a framework for testing and analyzing different simulation scenarios. This file allows users to set up experimental conditions to evaluate how variations in parameters affect orbital dynamics and system behavior.

- ``

  - Implements the core simulation logic, including the numerical integration methods used to calculate the trajectories of celestial bodies. This module handles the iterative updates of positions and velocities based on physical laws and initial conditions.

### Data Files

- ``

  - Contains a detailed comparison of actual and simulated orbital periods for each celestial body (Mercury, Venus, Earth, Mars). This file serves as a reference to assess the accuracy of the simulation.
  - Example content:
    ```
    The actual orbital period of Mercury is 87.97 days.
    The average orbital period of Mercury in this simulation is 87.95405982905983 days.
    ```

- ``

  - Defines the initial setup for the simulation, including global parameters (e.g., timestep, number of iterations) and the properties of each celestial body (e.g., mass, initial position, velocity, radius, and color).
  - Example content:
    ```
    timestep: 5000
    number of iterations: 20000

    name: Mercury
    mass: 3.3010e23
    initial position: 5.7909e10,0
    initial velocity: 0,0
    radius: 2.4397e6
    colour: b
    ```

## Features

- **Numerical Integration:** Simulates orbital dynamics using techniques like the Beeman scheme.
- **Celestial Animations:** Visualizes the motion of celestial bodies based on their calculated trajectories.
- **Orbital Period Analysis:** Compares simulated periods with real-world values for validation.
- **Scalable Parameters:** Easily adjustable parameters (e.g., timestep, number of iterations) for experimentation.

## Requirements

- Python 3.8+
- Required libraries:
  - `numpy`
  - `matplotlib`
  - Any additional dependencies listed in the `requirements.txt` file (if available).

## How to Run

1. Install the required libraries:
   ```bash
   pip install -r requirements.txt
   ```
2. Run the main script:
   ```bash
   python main.py
   ```
3. View results and animations as they are generated.

## File Descriptions

- **Input Data:** Located in `demo.txt`. Modify this file to change simulation parameters or celestial properties.
- **Output Data:** Orbital period comparisons are stored in `average_periods.txt` for analysis.

## Customization

- Update the `demo.txt` file to add or modify celestial bodies.
- Change the timestep or iteration count for more detailed or faster simulations.

## Future Enhancements

- Extend support for additional celestial bodies.
- Implement more advanced numerical integration techniques.
- Include real-time 3D visualization using libraries like `pythreejs` or `vpython`.

## Credits

Developed by **Sian Shu Wang** and contributors.

## License

This project is licensed under the MIT License. See `LICENSE` for details.

