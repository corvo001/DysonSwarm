<h1 align="center">Dyson Swarm</h1>

<p align="center">
  <img src="https://github.com/user-attachments/assets/f52e7d8d-4ea9-413c-b7c7-97f14cc073ad" alt="Dyson Swarm orbital simulation in motion" width="800" />
</p>

## Description

A theoretical **Dyson Swarm simulator** based on a modular network of collector satellites orbiting a star.

Developed by [Daniel Cuervo](https://corvo01.com/es/sobre-mi/) under the [CORVO](https://corvo01.com/) brand.

[Proyecto en español](https://corvo01.com/es/proyectos/dyson-swarm/) ·
[Project in English](https://corvo01.com/en/work/dyson-swarm/)

## What the simulator models

Each satellite follows an independent elliptical orbit around a star. The interface can generate a new distribution and show a simplified collection estimate, both as a total and as an average per satellite. These values are in model units. Mutual gravity and collisions are outside the current model.

See [the orbital and collection calculations](https://github.com/corvo001/DysonSwarm/blob/7d64adfdb06aca99f2e0bacd2e879b9deaf8e3ab/main.cpp#L176-L258) for the implementation.

## Roadmap

1. Basic orbital simulation
2. Collector deployment
3. Swarm coordination and collision avoidance
4. Energy network
5. Large-scale expansion

## Future

- Autonomous replication
- Interstellar swarms
- Shkadov thruster concepts

## License

MIT License.
