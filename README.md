# Finite Elements Project

## Authors
- *Agla Descamp*
- *Henri Cornette*

## Course Information
- *Course*: LEPL1110 - Finite Elements
- *Institution*: Université catholique de Louvain (UCLouvain)

## Description

This folder contains the processor for the LEPL1110 course at UCLouvain on Finite Elements.
The processor is used for solving problems related to linear elasticity.

## Build Instructions

To build the processor, use the following command:

```bash
mkdir build && cd build && cmake .. && make && ./myFem
```

## The solver

This program is run with a band solver and a RCM renumeration

## Mesh Types

You can decide two types of mesh, triangle and quadrilateral, you can change them by modifying the elementType variable in the main.c file:

- FEM_TRIANGLE: Triangle mesh.
- FEM_QUAD: Quadrilateral mesh.