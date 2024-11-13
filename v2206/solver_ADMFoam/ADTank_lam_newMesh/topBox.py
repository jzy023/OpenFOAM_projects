import numpy as np
from stl import mesh

def create_box_stl(x0, y0, z0, x1, y1, z1, filename):
    # Define the 8 vertices of the box
    vertices = np.array([
        [x0, y0, z0],  # Vertex 0
        [x1, y0, z0],  # Vertex 1
        [x1, y1, z0],  # Vertex 2
        [x0, y1, z0],  # Vertex 3
        [x0, y0, z1],  # Vertex 4
        [x1, y0, z1],  # Vertex 5
        [x1, y1, z1],  # Vertex 6
        [x0, y1, z1],  # Vertex 7
    ])
    
    # Define the 12 triangles composing the box (two triangles per face)
    faces = np.array([
        [0, 3, 1], [1, 3, 2],  # Bottom face
        [0, 1, 4], [4, 1, 5],  # Front face
        [1, 2, 5], [5, 2, 6],  # Right face
        [2, 3, 6], [6, 3, 7],  # Back face
        [0, 4, 3], [3, 4, 7],  # Left face
        [4, 5, 7], [7, 5, 6]   # Top face
    ])
    
    # Create the mesh
    box = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, face in enumerate(faces):
        for j in range(3):
            box.vectors[i][j] = vertices[face[j], :]

    # Write the mesh to an STL file
    box.save(filename)
    print(f"STL file '{filename}' created successfully.")

# Example usage
x0, y0, z0 = -150, -150, 118
x1, y1, z1 =  150, 150, 180
create_box_stl(x0, y0, z0, x1, y1, z1, "constant/triSurface/GeoTop.stl")
